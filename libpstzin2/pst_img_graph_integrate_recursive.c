/* See {pst_img_graph_integration_recursive.h} */
/* Last edited on 2025-01-14 17:34:22 by stolfi */
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
#include <pst_img_graph_shrink.h>
#include <pst_img_graph_integrate.h>

#include <pst_img_graph_integrate_recursive.h>

#define DELETED pst_img_graph_mark_DELETED

void pst_img_graph_copy_heights_from_shrunk(pst_img_graph_t *jg, double jZ[], pst_img_graph_t *ig, double iZ[], uint32_t iv_from_jv[]);
  /* Copies the heights {jZ[0..jg.NV-1]} computed on the shrunk graph
    {jg} to the heights vector {iZ[0..ig.NV-1]} for the original graph
    {iG}. Assumes that {iv_from_jv[jv]} is the index in {ig} of the
    vertex with index {jv} in {ig}. If {iv} is the index of a vertex of
    {ig} that is not a vertex of {jg}, the entry {iZ[iv]} is set to
    {NAN}. */

void pst_img_graph_get_height_estimates_from_shrunk(pst_img_graph_t *jg, double jZ[], pst_img_graph_t *g, double iZ[], uint32_t iv_from_jv[]);
  /* Copies the heights {jZ[0..jg.NV-1]} to the corresponding entries of {iZ[0..ig.NV-1]}, as per
    {pst_img_graph_copy_solution_from_shrunk}.  Then, for every vertex {iv} of {ig}
    that is not in {jg}, replaces the entry {iZ[iv]} (which will be {NAN}) by the 
    weighted average of the heights of it sneighbors.  If the vertex {iv} has no neighbors,
    sets {iZ[iv]} to zero. */
  
/* IMPLEMENTATIONS */
  
void pst_img_graph_copy_heights_from_shrunk(pst_img_graph_t *jg, double jZ[], pst_img_graph_t *ig, double iZ[], uint32_t iv_from_jv[])
  {
    assert(ig->NV >= jg->NV);
    for (uint32_t iv = 0; iv < ig->NV; iv++) { iZ[iv] = NAN; }
    for (uint32_t jv = 0; jv < jg->NV; jv++)
      { uint32_t iv = iv_from_jv[jv];
        iZ[iv] = jZ[jv];
      }
  }

void pst_img_graph_get_height_estimates_from_shrunk(pst_img_graph_t *jg, double jZ[], pst_img_graph_t *ig, double iZ[], uint32_t iv_from_jv[])
  {
    assert(ig->NV >= jg->NV);
    pst_img_graph_copy_heights_from_shrunk(jg, jZ, ig, iZ, iv_from_jv);
    /* Replace {NAN}s by averages of neighbors: */
    for (uint32_t iv = 0; iv < ig->NV; iv++)
      { pst_img_graph_vertex_data_t *vdi = &(ig->vdata[iv]);
        if (isnan(iZ[iv]))
          { haf_arc_t a0 = vdi->aout;
            if (a0 == NULL) 
              { /* Isolated vertex: */
                iZ[iv] = 0;
              }
            else
              { haf_arc_t a =  a0;
                double sW = 0;
                double sWZ = 0;
                do
                  { uint32_t iv_dst = pst_img_graph_get_arc_origin(ig, haf_sym(a));
                    assert(! isnan(iZ[iv_dst]));
                    double d = pst_img_graph_get_arc_delta(ig, a);
                    double w = pst_img_graph_get_edge_weight(ig, a);
                    sW += w;
                    sWZ += w*(iZ[iv_dst] - d);
                    a = haf_onext(a);
                  } while (a0 != a);
                assert(sW > 0);
                iZ[iv] = sWZ/sW;
              }
          }
      }
  }

void pst_img_graph_integration_recursive
  ( pst_img_graph_t* g,
    double Z[],
    int32_t level,
    uint32_t maxIter,
    double convTol, 
    bool_t para, 
    bool_t verbose,
    pst_img_graph_integrate_report_data_proc_t *reportData,
    pst_imgsys_report_sys_proc_t *reportSys,
    uint32_t reportStep,
    pst_integrate_report_heights_proc_t *reportHeights
  )
  {
    int32_t indent = (level < -1 ? 0 : 2*level + 2);
    if (verbose) { fprintf(stderr,"%*sstarting level %d with %d vertices\n", indent,"", level, g->NV); }
    if (reportData != NULL) { reportData(level, g); }

    if ( g->NV >= 2)
      {
        if (verbose) { fprintf(stderr,"%*sreducing graph ...\n", indent,""); }
        int32_t *jv_from_iv = talloc(ig->NV, int32_t);
        pst_img_graph_t *jg = pst_img_graph_shrink(g, jv_from_iv);
        double *jZ = rn_alloc(jg->NV);
        double ratio = (((double)g->NV)/((double)jg->NV));
        double newConvTol = convTol/sqrt(ratio);
        uint32_t newmaxIter = (uint32_t)ceil(sqrt(ratio)*(double)maxIter);
        if (verbose) 
          { fprintf(stdout,"%*sg.N = %d jg.N = %d ratio = %9.6lf\n", indent,"", ratio g->NV, jg->NV);
            fprintf(stderr,"%*snewConvTol = %.4e newMaxIter = %d\n", indent,"", newConvTol, newmaxIter);
          }
        pst_img_graph_integration_recursive
          ( jg, jZ, level+1, newMaxIter, newConvTol, para,  
            verbose, reportData, reportSys, reportStep, reportHeights
          );
        pst_img_graph_get_height_estimates_from_shrunk(jg, jZ, g, iZ, jv_from_iv);
        /* The vertices must be painted now*/
        free(jZ);
        free(jv_from_iv);
      }
    else
      { fprintf(stderr,"%*send of recursion - returning\n", indent,"");
        for (uint32_t iv = 0; iv < g->NV; iv++) { iZ[iv] = 0; }
      }

    fprintf(stderr,"%*ssolving level %d with %dvertices\n", indent,"", level, g->NV);

    pst_img_graph_integrate_iterative
      ( g, NX, NY, iZ, topoSort, convTol, para, szero, verbose, 
        level, reportSys, reportStep, reportHeights 
      );
  }

