/* see cpk_debug.h */
/* Last edited on 2024-12-31 16:34:11 by stolfi */ 

#include <stdint.h>
#include <stdio.h>
#include <assert.h>

#include <r2.h>
#include <pqueue.h>

#include <cpk_debug.h>
#include <cpk_io.h>
#include <cpk_graph.h>
#include <cpk_basic.h>

void cpk_print_vertices(FILE *wr, char *title, r2_vec_t *V, double_vec_t *W, uint32_t nMax)
  { uint32_t nV = V->ne;
    fprintf(wr, "\n--------------------\n");
    fprintf(wr, "%s (count = %d):\n", title, nV);
    for (uint32_t u = 0; u < nV; u++)
      { if ((nMax < 0) || (u < nMax) || (u >= nV - nMax))
          { fprintf(wr, "  " VTX_FMT, u); 
            r2_t *up = &(V->e[u]);
            r2_gen_print(wr, up, XY_FMT, " = (", ",", ")");
            if (W != NULL) 
              { fprintf(wr, " W = " WT_FMT, W->e[u]); }
            fprintf(wr, "\n");
          }
        else if (u == nMax)
          { fprintf(wr, "  ...\n"); }
      }
    fprintf(wr, "\n--------------------\n");
    fprintf(wr, "\n");
  }

void cpk_print_edges(FILE *wr, char *title, ui2_vec_t *E, r2_vec_t *V, uint32_t nMax)
  { uint32_t nE = E->ne;
    fprintf(wr, "\n--------------------\n");
    fprintf(wr, "%s (count = %d):\n", title, E->ne);
    for (uint32_t ie = 0; ie < nE; ie++)
      { ui2_t *Ei = &(E->e[ie]);
        uint32_t u = ORG(*Ei);
        uint32_t v = DST(*Ei);
        if ((nMax < 0) || (ie < nMax) || (ie >= nE - nMax))
          { cpk_ui2_print(wr, Ei, VTX_FMT);
            if (V != NULL)
              { r2_t *up = &(V->e[u]);
                r2_t *vp = &(V->e[v]);
                fprintf(wr, "  at  ");
                r2_gen_print(wr, up, XY_FMT, "(", ",", ")");
                fprintf(wr, " -- ");
                r2_gen_print(wr, vp, XY_FMT, "(", ",", ")");
                fprintf(wr, "  d = %8.5f", r2_dist(up, vp));
              }
            fprintf(wr, "\n");
          }
        else if (ie == nMax)
          { fprintf(wr, "  ...\n"); }
      }
    fprintf(wr, "\n--------------------\n");
    fprintf(wr, "\n");
  }

void cpk_print_vertex_set(FILE *wr, char *title, uint32_t nJ, uint32_t *J, r2_vec_t *V, uint32_t nMax)
  { 
    fprintf(wr, "\n--------------------\n");
    fprintf(wr, "%s (size = %d):%s", title, nJ, (V == NULL ? "" : "\n"));
    for(uint32_t i = 0; i < nJ; i++) 
      { if ((i < nMax) || (nJ - i <= nMax))
          { uint32_t u = J[i];
            if (V == NULL)
              { fprintf(wr, " " VTX_FFMT, u); }
            else
              { fprintf(wr, "  " VTX_FMT, u);
                r2_t *up = &(V->e[u]);
                fprintf(wr, "  at  ");
                r2_gen_print(wr, up, XY_FMT, "(", ",", ")");
                fprintf(wr, "\n");
              }
          }
        else if (i == nMax)
          { if (V == NULL)
              { fprintf(wr, " ... "); }
            else
              { fprintf(wr, "  ...\n"); }
          }
      }
    if (V == NULL) { fprintf(wr, "\n"); }
    fprintf(wr, "\n--------------------\n");
    fprintf(wr, "\n");
  }


void TRACE_SOL(char *title, uint32_t nS, double WS)
  { fprintf(stderr, "%s: tam = %d peso = " WT_FFMT "\n", title, nS, WS); }

void TRACE_Q(char *title, pqueue_t *Q)
  { fprintf(stderr, " %s:", title);
    uint32_t nQ = pqueue_count(Q);
    for (int32_t i = 0; i < nQ; i++)
      { pqueue_item_t u = pqueue_item(Q, (uint32_t)i); 
        double eu = pqueue_value(Q,u); 
        fprintf(stderr, " %d:" WT_FFMT, u, eu); 
      }
    fprintf(stderr, "\n");
  }
