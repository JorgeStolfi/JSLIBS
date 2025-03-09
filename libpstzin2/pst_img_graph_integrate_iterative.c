/* See {pst_img_graph.h} */
/* Last edited on 2025-02-24 06:40:28 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <haf.h>
#include <float_image.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <rn.h>

#include <pst_imgsys.h>
#include <pst_img_graph.h>
#include <pst_img_graph_integrate.h>

#include <pst_img_graph_integrate_iterative.h>

void pst_img_graph_integrate_iterative
  ( pst_img_graph_t *g,
    double Z[],
    bool_t sortSys,
    uint32_t maxIter,
    double convTol, 
    bool_t para, 
    bool_t verbose,
    int32_t level,
    pst_imgsys_report_sys_proc_t *reportSys,
    uint32_t reportStep,
    pst_integrate_report_heights_proc_t *reportHeights
  )
  { int32_t indent = (level < -1 ? 0 : 2*level+2);

    if (verbose) { fprintf(stderr, "%*sbuilding the system ...\n", indent, ""); }
    int32_t *zid_from_vid = talloc(g->NV, int32_t);
    pst_imgsys_t *S = pst_img_graph_integrate_build_system(g, zid_from_vid, verbose);
    uint32_t NZ = S->N;
    if (verbose) { fprintf(stderr, "%*ssystem has %d equations and %d variables\n", indent, "", NZ, NZ); }

    if (verbose) { fprintf(stderr, "%*scopying initial guess ...\n", indent, ""); }
    double *h = talloc(NZ, double);
    for (uint32_t kv = 0; kv < g->NV; kv++)
      { int32_t kz = zid_from_vid[kv];
        if (kz != -1)
          { /* Grab initial guess: */
            assert((kz >= 0) && (kz < NZ));
            h[kz] = Z[kv];
          }
      }
    
    if (verbose) { fprintf(stderr, "%*ssolving the system ...\n", indent, ""); }
    
    float_image_t *Zimg = float_image_new(1, NX, NY);

    auto void reportSol(int32_t level, int32_t iter, double change, bool_t final, uint32_t N, double Z[]);

    uint32_t *ord = NULL;
    if (sortSys) { ord = pst_imgsys_sort_equations(S); }
    bool_t szero = TRUE;
    pst_imgsys_solve_iterative
      ( S, h, ord, maxIter, convTol, 
        para, szero, verbose, level, 
        reportStep,
        (reportHeights == NULL ? NULL : &reportSol)
      );

    if (verbose) { fprintf(stderr, "%*sreturning solution ...\n", indent, ""); }
    for (uint32_t kv = 0; kv < g->NV; kv++)
      { int32_t kz = zid_from_vid[kv];
        if (kz == -1)
          { Z[kv] = 0.0; }
        else
          { assert((kz >= 0) && (kz < NZ));
            Z[kv] = h[kz];
          }
      }

    if(ord != NULL) { free(ord); }
    free(h);
    pst_imgsys_free(S);
    float_image_free(Zimg);
    return;
    
    void reportSol(int32_t level, int32_t iter, double change, bool_t final, uint32_t N, double Z[])
      { assert(reportHeights != NULL);
        pst_imgsys_copy_sol_vec_to_image(S, Z, Zimg, 0.0);
        pst_imgsys_extract_system_eq_tot_weight_image(S, Uimg, 0.0)
        reportHeights(level, iter, change, final, Zimg, Uimg);
      }
  }
