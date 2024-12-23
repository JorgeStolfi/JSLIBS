/* See {pst_graph_integrate.h}. */
/* Last edited on 2024-12-22 21:40:34 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include <assert.h>
#include <affirm.h>
#include <jsfile.h>
#include <rn.h>

#include <pst_graph.h>
#include <pst_imgsys.h>
#include <pst_slope_map.h>

#include <pst_graph_integrate.h>

pst_imgsys_t *pst_graph_build_integration_system(pst_graph_t *g, uint32_t NX_Z, uint32_t NY_Z)
  {
    uint32_t NXY_Z = NX_Z*NY_Z;
    int32_t *ind_ix = talloc(g->num_vertex, int32_t); /* Maps vertex index to unknown index, or {-1}. */
    pst_imgsys_equation_t *eq = talloc(g->num_vertex, pst_imgsys_equation_t);
    uint32_t N = 0;
    for (uint32_t i = 0; i < g->num_vertex; i++)
      { pst_imgsys_equation_t *eqk = &(eq[N]);

        uint32_t nt = 0; /* Number of terms in equation. */
        eqk->rhs = 0.0;
        eqk->ix[nt] = i; eqk->cf[nt] = 1.00; nt++;
        ind_ix[i] = (int32_t)N;
        pst_vertex_t *v = &(g->vertices[i]);
        uint32_t has_xp = (uint32_t)(v->vtxp != -1); 
        uint32_t has_xm = (uint32_t)(v->vtxm != -1);
        uint32_t has_yp = (uint32_t)(v->vtyp != -1);
        uint32_t has_ym = (uint32_t)(v->vtym != -1);

        uint32_t num_neighbours = has_xp + has_xm + has_yp + has_ym;
        if (num_neighbours > 0)
          { if (has_xm)
              { int32_t vxm;
                double dxm, wxm;
                pst_graph_vertex_get_neighbour(g, i, -1, 0, &vxm, &dxm, &wxm);
                assert(vxm != -1);
                eqk->ix[nt] = (uint32_t)vxm;
                eqk->cf[nt] = -wxm;
                eqk->rhs += +wxm*dxm;
                nt++;
              }
            if (has_xp)
              { int32_t vxp;
                double dxp, wxp;
                pst_graph_vertex_get_neighbour(g, i, 1, 0, &vxp, &dxp, &wxp);
                assert(vxp != -1);
                eqk->ix[nt] = (uint32_t)vxp;
                eqk->cf[nt] = -wxp;
                eqk->rhs += +wxp*dxp;
                nt++;
              }
            if (has_ym)
              { int32_t vym;
                double dym, wym;
                pst_graph_vertex_get_neighbour(g, i, 0, -1, &vym, &dym, &wym);
                assert(vym != -1);
                eqk->ix[nt] = (uint32_t)vym;
                eqk->cf[nt] = -wym;
                eqk->rhs += +wym*dym;
                nt++;
              }
            if (has_yp)
              { int32_t vyp;
                double dyp, wyp;
                pst_graph_vertex_get_neighbour(g, i, 0, 1, &vyp, &dyp, &wyp);
                assert(vyp != -1);
                eqk->ix[nt] = (uint32_t)vyp;
                eqk->cf[nt] = -wyp;
                eqk->rhs += +wyp*dyp;
                nt++;
              }
            double wtot = 0;
            for (uint32_t j = 1; j < nt; j++) { assert(eqk->cf[j] < 0); wtot+= -eqk->cf[j]; } 
            assert(wtot > 0);
            for (uint32_t j = 1; j < nt; j++) { eqk->cf[j]/=wtot;  }
            eqk->rhs /= wtot;
            eqk->nt = nt;
            N++;
          }
        else
          { ind_ix[i] = -1; }
      }

    for (uint32_t k = 0; k < N; k++)
      { pst_imgsys_equation_t *eqk = &(eq[k]);
        uint32_t nt = eqk->nt;
        uint32_t mt = 0;
        for (uint32_t i = 0; i < nt; i++)
          { /* Get the temporay index {xyi}: */
            uint32_t xyi = eqk->ix[i];
            /* Get the definitive index {ki}: */
            int32_t ki = ind_ix[xyi];
            if (ki >= 0)
              { /* Append the term to the equation: */
                int32_t j = (int32_t)mt;
                eqk->ix[j] = (uint32_t)ki;
                eqk->cf[j] = eqk->cf[i];
                assert(!isnan(eqk->cf[j]));
                mt++;
              }
          }
        eqk->nt = mt;
      }

    /* Build the inverse tables {col[0..N-1], row[0..N-1]}: */
    uint32_t *col = talloc(N, uint32_t);
    uint32_t *row = talloc(N, uint32_t);
    int32_t *ix = talloc(NXY_Z, int32_t);
      { for (uint32_t xy = 0; xy < NXY_Z; xy++) { ix[xy] = -1; }
        uint32_t count_idx = 0;
        for (uint32_t xy = 0; xy <g->num_vertex; xy++) 
          { if (ind_ix[xy] >= 0)
              { int32_t k = ind_ix[xy];
                int32_t x, y;
                pst_graph_restore_vertex_index((int32_t)g->vertices[xy].id, NX_Z, NY_Z, &x, &y);
                assert((x >= 0) && (x < NX_Z));
                assert((y >= 0) && (y < NY_Z));
                col[count_idx] = (uint32_t)x;
                row[count_idx] = (uint32_t)y;
                assert((k >> 30) == 0); /* Should never happen. */
                ix[g->vertices[xy].id] = (int32_t)k;
                count_idx++;
              }
          }
        assert(count_idx == N);
      }

    free(ind_ix);
    /* Now package the equations as a system: */
    pst_imgsys_t *S = pst_imgsys_from_eqs(NX_Z, NY_Z, N, eq, ix, col, row);  
    return S;
  }

void pst_graph_solve_system
  ( pst_graph_t *g,
    pst_imgsys_t *S,
    float_image_t *OZ, 
    uint32_t maxIter, 
    double convTol, 
    bool_t para, 
    bool_t szero, 
    bool_t verbose 
  )
  { double *Z = rn_alloc(S->N);
    pst_slope_map_copy_height_map_to_sol_vec(S, OZ, Z);
    pst_imgsys_solve(S, Z, NULL, maxIter, convTol, para, szero, verbose, 0, NULL);
    pst_slope_map_copy_sol_vec_to_height_map(S, Z, OZ);
    free(Z);
  }

