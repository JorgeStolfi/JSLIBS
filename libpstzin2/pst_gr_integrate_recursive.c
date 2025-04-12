/* See {pst_gr_integration_recursive.h} */
/* Last edited on 2025-03-11 14:55:17 by stolfi */
/* Created by Rafael F. V. Saracchini */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <rn.h>

#include <pst_imgsys.h>
#include <pst_gr.h>
#include <pst_gr_shrink.h>
#include <pst_gr_integrate.h>

#include <pst_gr_integrate_recursive.h>

#define DELETED pst_gr_mark_DELETED

void pst_gr_copy_heights_from_shrunk(pst_gr_t *grj, double jZ[], pst_gr_t *gri, double iZ[], uint32_t vi_from_jv[]);
  /* Copies the heights {jZ[0..grj.NV-1]} computed on the shrunk graph
    {grj} to the heights vector {iZ[0..gri.NV-1]} for the original graph
    {gri}. Assumes that {vi_from_vj[vj]} is the index in {gri} of the
    vertex with index {vj} in {gri}. If {vi} is the index of a vertex of
    {gri} that is not a vertex of {grj}, the entry {iZ[vi]} is set to
    {NAN}. */

void pst_gr_get_height_estimates_from_shrunk(pst_gr_t *grj, double jZ[], pst_gr_t *gr, double iZ[], uint32_t vi_from_vj[]);
  /* Copies the heights {jZ[0..grj.NV-1]} to the corresponding entries of {iZ[0..gri.NV-1]}, as per
    {pst_gr_copy_solution_from_shrunk}.  Then, for every vertex {vi} of {gri}
    that is not in {grj}, replaces the entry {iZ[vi]} (which will be {NAN}) by the 
    weighted average of the heights of it sneighbors.  If the vertex {vi} has no neighbors,
    sets {iZ[vi]} to zero. */
  
/* IMPLEMENTATIONS */
  
void pst_gr_copy_heights_from_shrunk(pst_gr_t *grj, double jZ[], pst_gr_t *gri, double iZ[], uint32_t vi_from_vj[])
  {
    assert(gri->NV >= grj->NV);
    for (uint32_t vi = 0; vi < gri->NV; vi++) { iZ[vi] = NAN; }
    for (uint32_t vj = 0; vj < grj->NV; vj++)
      { uint32_t vi = vi_from_vj[vj];
        iZ[vi] = jZ[vj];
      }
  }

void pst_gr_get_height_estimates_from_shrunk(pst_gr_t *grj, double jZ[], pst_gr_t *gri, double iZ[], uint32_t vi_from_vj[])
  {
    assert(gri->NV >= grj->NV);
    pst_gr_copy_heights_from_shrunk(grj, jZ, gri, iZ, vi_from_vj);
    /* Replace {NAN}s by averages of neighbors: */
    for (uint32_t vi = 0; vi < gri->NV; vi++)
      { pst_gr_vertex_data_t *vdi = &(gri->vdata[vi]);
        if (isnan(iZ[vi]))
          { pst_gr_arc_t a0 = vdi->aout;
            if (a0 == NULL) 
              { /* Isolated vertex: */
                iZ[vi] = 0;
              }
            else
              { pst_gr_arc_t a =  a0;
                double sW = 0;
                double sWZ = 0;
                do
                  { uint32_t vi_dst = pst_gr_arc_org(gri, pst_gr_arc_sym(a));
                    assert(! isnan(iZ[vi_dst]));
                    double d = pst_gr_arc_delta(gri, a);
                    double w = pst_gr_arc_weight(gri, a);
                    sW += w;
                    sWZ += w*(iZ[vi_dst] - d);
                    a = pst_gr_arc_onext(a);
                  } while (a0 != a);
                assert(sW > 0);
                iZ[vi] = sWZ/sW;
              }
          }
      }
  }

void pst_gr_integration_recursive
  ( pst_gr_t* gr,
    double Z[],
    int32_t level,
    uint32_t maxIter,
    double convTol, 
    bool_t para, 
    bool_t verbose,
    pst_gr_integrate_report_data_proc_t *reportData,
    pst_imgsys_report_sys_proc_t *reportSys,
    uint32_t reportStep,
    pst_integrate_report_heights_proc_t *reportHeights
  )
  {
    int32_t indent = (level < -1 ? 0 : 2*level+2);
    if (verbose) { fprintf(stderr,"%*sstarting level %d with %d vertices\n", indent,"", level, gr->NV); }
    if (reportData != NULL) { reportData(level, gr); }

    if ( gr->NV >= 2)
      {
        if (verbose) { fprintf(stderr,"%*sreducing graph ...\n", indent,""); }
        int32_t *vj_from_vi = talloc(gri->NV, int32_t);
        pst_gr_t *grj = pst_gr_shrink(gr, vj_from_vi);
        double *jZ = rn_alloc(grj->NV);
        double ratio = (((double)gr->NV)/((double)grj->NV));
        double newConvTol = convTol/sqrt(ratio);
        uint32_t newmaxIter = (uint32_t)ceil(sqrt(ratio)*(double)maxIter);
        if (verbose) 
          { fprintf(stdout,"%*sg.N = %d grj.N = %d ratio = %9.6lf\n", indent,"", ratio gr->NV, grj->NV);
            fprintf(stderr,"%*snewConvTol = %.4e newMaxIter = %d\n", indent,"", newConvTol, newmaxIter);
          }
        pst_gr_integration_recursive
          ( grj, jZ, level+1, newMaxIter, newConvTol, para,  
            verbose, reportData, reportSys, reportStep, reportHeights
          );
        pst_gr_get_height_estimates_from_shrunk(grj, jZ, gr, iZ, vj_from_vi);
        /* The vertices must be painted now*/
        free(jZ);
        free(vj_from_vi);
      }
    else
      { fprintf(stderr,"%*send of recursion - returning\n", indent,"");
        for (uint32_t vi = 0; vi < gr->NV; vi++) { iZ[vi] = 0; }
      }

    fprintf(stderr,"%*ssolving level %d with %dvertices\n", indent,"", level, gr->NV);

    pst_gr_integrate_iterative
      ( gr, NX, NY, iZ, sortSys, convTol, para, szero, verbose, 
        level, reportSys, reportStep, reportHeights 
      );
  }

