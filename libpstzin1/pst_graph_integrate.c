/* See {pst_graph_integrate.h}. */
/* Last edited on 2025-01-13 15:22:39 by stolfi */
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
#include <pst_graph_from_maps.h>

#include <pst_graph_integrate.h>

pst_imgsys_t *pst_graph_build_system(pst_graph_t *g, int32_t NX_Z, int32_t NY_Z)
  {
    int32_t NXY_Z = NX_Z*NY_Z;
    int32_t *ind_ix = talloc(NXY_Z, int32_t);
    /* Build a table {uid} that maps vertex index to equation/height index, or {-1}: */
    uint32_t NZ = 0; /* Number of equations and variables. */
    int32_t *uid = talloc(g->NV, int32_t);
    for (uint32_t kv = 0; kv < g->NV ; kv++)
      { pst_vertex_t *vdk = &(g->vertices[vk]);
        if (vdk->vmark == DELETED) 
          { uid[vk] = -1; }
        else
          { uid[vk] = (int32_t)NZ; NZ++; }
      }
    /* Now build the equations: */
    pst_imgsys_equation_t *eq = talloc(NZ, pst_imgsys_equation_t);
    for (uint32_t kv = 0; kv < g->NV; kv++)
      { int32_y uidk = uid[kv];
        if (uidk != -1)
          { assert((uidk >= 0) && (uidk < NZ));
            pst_imgsys_equation_t *eqk = &(eq[N]);
            eqk->nt = 0; /* Number of terms in equation. */
            eqk->rhs = 0.0;
            eqk->uid[nt] = kv; eqk->cf[nt] = 0.00; nt++;
            pst_vertex_t *vk = &(g->vertices[kv]);
            uint32_t has_xp = (uint32_t)(vk->vtxp != -1); 
            uint32_t has_xm = (uint32_t)(vk->vtxm != -1);
            uint32_t has_yp = (uint32_t)(vk->vtyp != -1);
            uint32_t has_ym = (uint32_t)(vk->vtym != -1);

        uint32_t num_neighbours = has_xp + has_xm + has_yp + has_ym;
        if (num_neighbours > 0)
          { if (has_xm)
              { int32_t vxm;
                double dxm, wxm;
                pst_graph_vertex_get_neighbour(g, kv, -1, 0, &vxm, &dxm, &wxm);
                assert(vxm != -1);
                eqk->uid[nt] = (uint32_t)vxm;
                eqk->cf[nt] = -wxm;
                eqk->rhs += +wxm*dxm;
                nt++;
              }
            if (has_xp)
              { int32_t vxp;
                double dxp, wxp;
                pst_graph_vertex_get_neighbour(g, kv, 1, 0, &vxp, &dxp, &wxp);
                assert(vxp != -1);
                eqk->uid[nt] = (uint32_t)vxp;
                eqk->cf[nt] = -wxp;
                eqk->rhs += +wxp*dxp;
                nt++;
              }
            if (has_ym)
              { int32_t vym;
                double dym, wym;
                pst_graph_vertex_get_neighbour(g, kv, 0, -1, &vym, &dym, &wym);
                assert(vym != -1);
                eqk->uid[nt] = (uint32_t)vym;
                eqk->cf[nt] = -wym;
                eqk->rhs += +wym*dym;
                nt++;
              }
            if (has_yp)
              { int32_t vyp;
                double dyp, wyp;
                pst_graph_vertex_get_neighbour(g, kv, 0, 1, &vyp, &dyp, &wyp);
                assert(vyp != -1);
                eqk->uid[nt] = (uint32_t)vyp;
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
          { ind_uid[kv] = -1; }
      }

    for (uint32_t k = 0; k < N; k++)
      { pst_imgsys_equation_t *eqk = &(eq[k]);
        uint32_t nt = eqk->nt;
        uint32_t mt = 0;
        for (uint32_t kv = 0; kv < nt; kv++)
          { /* Get the temporay index {xyi}: */
            uint32_t xyi = eqk->uid[kv];
            /* Get the definitive index {ki}: */
            int32_t ki = ind_uid[xyi];
            if (ki >= 0)
              { /* Append the term to the equation: */
                int32_t j = (int32_t)mt;
                eqk->uid[j] = (uint32_t)ki;
                eqk->cf[j] = eqk->cf[kv];
                assert(!isnan(eqk->cf[j]));
                mt++;
              }
          }
        eqk->nt = mt;
      }

    /* Build the inverse tables {col[0..N-1], row[0..N-1]}: */
    uint32_t *col = talloc(N, uint32_t);
    uint32_t *row = talloc(N, uint32_t);
    int32_t *uid = talloc(NXY_Z, int32_t);
      { for (uint32_t xy = 0; xy < NXY_Z; xy++) { uid[xy] = -1; }
        uint32_t count_idx = 0;
        for (uint32_t xy = 0; xy <g->NV; xy++) 
          { if (ind_uid[xy] >= 0)
              { int32_t k = ind_uid[xy];
                int32_t x, y;
                pst_graph_restore_vertex_index((int32_t)g->vertices[xy].id, NX_Z, NY_Z, &x, &y);
                assert((x >= 0) && (x < NX_Z));
                assert((y >= 0) && (y < NY_Z));
                col[count_idx] = (uint32_t)x;
                row[count_idx] = (uint32_t)y;
                assert((k >> 30) == 0); /* Should never happen. */
                uid[g->vertices[xy].id] = (int32_t)k;
                count_idx++;
              }
          }
        assert(count_idx == N);
      }

    free(ind_uid);
    /* Now package the equations as a system: */
    pst_imgsys_t *S = pst_imgsys_from_eqs(NX_Z, NY_Z, N, eq, uid, col, row);  
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
    pst_imgsys_copy_image_to_sol_vec(S, OZ, Z);
    pst_imgsys_solve_iterative(S, Z, NULL, maxIter, convTol, para, szero, verbose, 0, NULL);
    pst_imgsys_copy_sol_vec_to_image(S, Z, OZ);
    free(Z);
  }

