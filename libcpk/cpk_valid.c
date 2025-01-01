/* See cpk_valid.h. */
/* Last edited on 2024-12-31 16:36:21 by stolfi */ 

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <r2.h>

#include <cpk_valid.h>
#include <cpk_main.h>
#include <cpk_io.h>
#include <cpk_debug.h>
#include <cpk_basic.h>

bool_t cpk_check_edges (r2_vec_t *V, ui2_vec_t *E, double dMin)
  {
    uint32_t nV = V->ne;
    uint32_t nE = E->ne;
    /* Get the incompatibility edges the slow way: */
    ui2_vec_t ES = cpk_slow_get_proximity_edges(V, dMin);

    /* Print and check the edges */
    if (nE != ES.ne) 
      { fprintf(stderr, "** got %d incompatible pairs", nE);
        fprintf(stderr, "; should have got %d.\n", ES.ne);
        return FALSE;
      }
    for (int32_t ie = 0; ie < nE; ie++)
      { bool_t edge_ok = TRUE;
        ui2_t *Ei = &(E->e[ie]);
        uint32_t u = ORG(*Ei), v = DST(*Ei); 
        if ((u < 0) || (u >= nV)) 
          { fprintf(stderr, "** bad edge origin"); edge_ok = FALSE; }
        if ((v < 0) || (v >= nV)) 
          { fprintf(stderr, "** bad edge destination"); edge_ok = FALSE; }
        r2_t *up = &(V->e[u]), *vp = &(V->e[v]);
        ui2_t *ESi = &(ES.e[ie]);
        uint32_t uS = ORG(*ESi), vS = DST(*ESi);
        if ((u != uS) || (v != vS))
          { fprintf(stderr, "** should be ");
            cpk_ui2_print(stderr, ESi, VTX_FFMT);
            edge_ok = FALSE;
          }
        double duv = r2_dist(up, vp);
        if (duv >= dMin) 
          { fprintf(stderr, "** distance is " XY_FFMT ", min was " XY_FFMT "\n", duv, dMin);
            edge_ok = FALSE;
          }
        if (! edge_ok) 
          { cpk_ui2_print(stderr, Ei, VTX_FFMT);
            fprintf(stderr, "  at  ");
            r2_gen_print(stderr, up, XY_FFMT, "(", ",", ")");
            fprintf(stderr, " -- ");
            r2_gen_print(stderr, vp, XY_FFMT, "(", ",", ")");
            fprintf(stderr, "\n");
            return(FALSE);
          }
      }
    return TRUE;
  }

bool_t cpk_check_indep_set(cpk_graph_t *G, uint32_t nJ, uint32_t *J, bool_t sel[])
  {
    uint32_t nV = G->nV; 

    bool_t ok = TRUE;

    int32_t where[nV]; /* Inverse map {where[v] = -1} except {where[J[i]] = i}. */
    for (int32_t v = 0; v < nV; v++) { where[v] = -1; }
    for (int32_t i = 0; i < nJ; i++)
      { uint32_t v = J[i];
        if ((v < 0) || (v >= nV))
          { fprintf(stderr, "** vertex J[%d] = %d out of bounds\n", i, v);
            ok = FALSE;
          }
        else if (where[v] >= 0)
          { fprintf(stderr, "** vertex J[%d] = %d also in J[%d]\n", i, v, where[v]);
            ok = FALSE;
          }
        else if ((sel != NULL) && (! sel[v]))
          { fprintf(stderr, "** vertex J[%d] = %d not marked in bool_t vector\n", i, v);
            ok = FALSE;
          }
        else
          { /* Check for independence: */
            uint32_t *nv = &(G->nbr[G->fnb[v]]), dv = G->deg[v];
            for (int32_t k = 0; k < dv; k++)
              { uint32_t u = nv[k];
                if ((u < v) && where[u] >= 0)
                  { fprintf(
                      stderr, 
                      "** incompat vertices J[%d] = %d and J[%d] = %d\n", 
                      i, v, where[u], u
                    );
                    ok = FALSE;
                  }
              }
          }
        where[v] = i;
      }
    /* Check for spurious marked vertices: */
    if (sel != NULL)
      { for (int32_t v = 0; v < nV; v++)
          { if ((where[v] < 0) && sel[v])
              { fprintf(stderr, "** vertex %d incorrectly marked in bool_t vector\n", v);
                ok = FALSE;
              }
          }
      }

    /* Print {J} if check failed: */
    if (! ok) 
      { cpk_print_vertex_set(stderr, "supposedly independent set", nJ, J, NULL, UINT32_MAX); }
    return ok;
  }

bool_t cpk_check_solution
  ( cpk_domain_t *C,
    cpk_policy_t *P,
    r2_vec_t *V,
    uint32_vec_t *J
  )
  { 
    uint32_t nJ = J->ne;
    
    /* Print the indep vertices and their nearest neighbors */
    fprintf(stderr, "Got %d compatible auction points\n", nJ);

    double dAAMin = P->dMin + 2*P->rAuc;

    bool_t ok = TRUE;
    for (int32_t uj = 0; uj < nJ; uj++)
      { uint32_t u = J->e[uj]; 
        r2_t *up = &(V->e[u]);

        /* !!! Should check whether {up} satisfies the constraints {C}. */

        /* Find nearest neighbor {v = J[vj]} of {u} in {J}: */
        uint32_t v = UINT32_MAX, vj = UINT32_MAX;
        r2_t *vp = NULL;
        double dv = INF;
        { for (uint32_t wj = 0; wj < nJ; wj++)
            { if (wj != uj) 
                { uint32_t w = J->e[wj];
                  r2_t *wp = &(V->e[w]);
                  double dw = r2_dist(up, wp); 
                  if (dw < dv) { v = w; vj = wj; vp = wp; dv = dw; }
                }
            }
        }

        /* Check min distance constraint: */
        if (dv < dAAMin)
          { fprintf(stderr, "** vertices J[%d]=%d ", uj, u);
            r2_gen_print(stderr, up, XY_FFMT, " = (", ",", ")");
            fprintf(stderr, " and J[%d] = %d ", vj, v);
            r2_gen_print(stderr, vp, XY_FFMT, " = (", ",", ")");
            fprintf(stderr, " (dist = " XY_FFMT ") are too close\n", dv);
            ok = FALSE;
          }
      }
    return ok;
  }

ui2_vec_t cpk_slow_get_proximity_edges (r2_vec_t *V, double dLim)
  { 
    /* Slow and dumb: */
    ui2_vec_t E = ui2_vec_new(V->ne);
    uint32_t nE = 0;
    for (uint32_t i = 0; i < V->ne; i++)
      { r2_t *Vi = &(V->e[i]);
        for (uint32_t j = i+1; j < V->ne; j++)
          { r2_t *Vj = &(V->e[j]);
            if (r2_dist(Vi, Vj) < dLim)
              { ui2_vec_expand(&E, (vec_index_t)nE);
                E.e[nE] = (ui2_t){{i,j}};
                nE++;
              }
          }
      }
    ui2_vec_trim(&E, nE);
    return E;
  }

