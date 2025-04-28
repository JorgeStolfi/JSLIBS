/* See {psgr_integration_recursive.h} */
/* Last edited on 2025-04-27 02:55:09 by stolfi */
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

#include <psgr_types.h>
#include <psgr.h>
#include <psgr_shrink.h>
#include <psgr_integrate.h>

#include <psgr_integrate_recursive.h>

void psgr_copy_heights_from_shrunk(psgr_t *grj, double jZ[], psgr_t *gri, double iZ[], uint32_t vi_from_jv[]);
  /* Copies the heights {jZ[0..grj.NV-1]} computed on the shrunk graph
    {grj} to the heights vector {iZ[0..gri.NV-1]} for the original graph
    {gri}. Assumes that {vi_from_vj[vj]} is the index in {gri} of the
    vertex with index {vj} in {gri}. If {vi} is the index of a vertex of
    {gri} that is not a vertex of {grj}, the entry {iZ[vi]} is set to
    {NAN}. */

void psgr_get_height_estimates_from_shrunk(psgr_t *grj, double jZ[], psgr_t *gr, double iZ[], uint32_t vi_from_vj[]);
  /* Copies the heights {jZ[0..grj.NV-1]} to the corresponding entries of {iZ[0..gri.NV-1]}, as per
    {psgr_copy_solution_from_shrunk}.  Then, for every vertex {vi} of {gri}
    that is not in {grj}, replaces the entry {iZ[vi]} (which will be {NAN}) by the 
    weighted average of the heights of it sneighbors.  If the vertex {vi} has no neighbors,
    sets {iZ[vi]} to zero. */
  
/* IMPLEMENTATIONS */
  
void psgr_copy_heights_from_shrunk(psgr_t *grj, double jZ[], psgr_t *gri, double iZ[], uint32_t vi_from_vj[])
  {
    assert(gri->NV >= grj->NV);
    for (uint32_t vi = 0; vi < gri->NV; vi++) { iZ[vi] = NAN; }
    for (uint32_t vj = 0; vj < grj->NV; vj++)
      { uint32_t vi = vi_from_vj[vj];
        iZ[vi] = jZ[vj];
      }
  }

void psgr_get_height_estimates_from_shrunk(psgr_t *grj, double jZ[], psgr_t *gri, double iZ[], uint32_t vi_from_vj[])
  {
    assert(gri->NV >= grj->NV);
    psgr_copy_heights_from_shrunk(grj, jZ, gri, iZ, vi_from_vj);
    /* Replace {NAN}s by averages of neighbors: */
    for (uint32_t vi = 0; vi < gri->NV; vi++)
      { psgr_vertex_data_t *vdi = &(gri->vdata[vi]);
        if (isnan(iZ[vi]))
          { psgr_arc_t a0 = vdi->aout;
            if (a0 == NULL) 
              { /* Isolated vertex: */
                iZ[vi] = 0;
              }
            else
              { psgr_arc_t a =  a0;
                double sW = 0;
                double sWZ = 0;
                do
                  { uint32_t vi_dst = psgr_arc_org(gri, psgr_arc_sym(a));
                    assert(! isnan(iZ[vi_dst]));
                    double d = psgr_arc_delta(gri, a);
                    double w = psgr_arc_weight(gri, a);
                    sW += w;
                    sWZ += w*(iZ[vi_dst] - d);
                    a = psgr_arc_onext(a);
                  } while (a0 != a);
                assert(sW > 0);
                iZ[vi] = sWZ/sW;
              }
          }
      }
  }

void psgr_integration_recursive
  ( psgr_t* gr,
    double Z[],
    int32_t level,
    uint32_t maxIter,
    double convTol, 
    bool_t para, 
    bool_t verbose,
    psgr_integrate_report_data_proc_t *reportData,
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
        psgr_t *grj = psgr_shrink(gr, vj_from_vi);
        double *jZ = rn_alloc(grj->NV);
        double ratio = (((double)gr->NV)/((double)grj->NV));
        double newConvTol = convTol/sqrt(ratio);
        uint32_t newmaxIter = (uint32_t)ceil(sqrt(ratio)*(double)maxIter);
        if (verbose) 
          { fprintf(stdout,"%*sg.N = %d grj.N = %d ratio = %9.6lf\n", indent,"", ratio gr->NV, grj->NV);
            fprintf(stderr,"%*snewConvTol = %.4e newMaxIter = %d\n", indent,"", newConvTol, newmaxIter);
          }
        psgr_integration_recursive
          ( grj, jZ, level+1, newMaxIter, newConvTol, para,  
            verbose, reportData, reportSys, reportStep, reportHeights
          );
        psgr_get_height_estimates_from_shrunk(grj, jZ, gr, iZ, vj_from_vi);
        /* The vertices must be painted now*/
        free(jZ);
        free(vj_from_vi);
      }
    else
      { fprintf(stderr,"%*send of recursion - returning\n", indent,"");
        for (uint32_t vi = 0; vi < gr->NV; vi++) { iZ[vi] = 0; }
      }

    fprintf(stderr,"%*ssolving level %d with %dvertices\n", indent,"", level, gr->NV);

    psgr_integrate_iterative
      ( gr, NX, NY, iZ, sortSys, convTol, para, szero, verbose, 
        level, reportSys, reportStep, reportHeights 
      );
  }

??? ----------------------------------------------------------------------
  
  
  fprintf(stderr,"Solving Graph level[%d] with [%d] vertices [%d] valid\n",level, gr->NV,gr->NV_valid);

  int32_t *ref_tab;
  pst_imgsys_t *S = psgr_integrate_build_system(gr, iW, (int32_t)OZ->sz[1], (int32_t)OZ->sz[2], &ref_tab);
  
  if (debug)
    {
      char *filename = jsprintf("%s-%02d-S.txt",out_prefix,level);
      FILE *wr_dump = open_write(filename,FALSE);
      pst_imgsys_write(wr_dump,S,"%+10.6f");
      fclose(wr_dump);
      free(filename);

      psgr_put_solution_into_image(gr,iZ,OZ);
      filename = jsprintf("%s-%02d-iZ.fni",out_prefix,level);
      wr_dump = open_write(filename,FALSE);
      float_image_write(wr_dump,OZ);
      fclose(wr_dump);
      free(filename);

      psgr_put_solution_into_image(gr,iW,OZ);
      filename = jsprintf("%s-%02d-iW.fni",out_prefix,level);
      wr_dump = open_write(filename,FALSE);
      float_image_write(wr_dump,OZ);
      fclose(wr_dump);
      free(filename);
    }
    
    if (RZ != NULL)
{
      psgr_put_error_into_image(gr,iZ,iW,RZ,OZ);
      filename = NULL;
      char *filename = jsprintf("%s-%02d-iE.fni",out_prefix,level);
      wr_dump = open_write(filename,FALSE);
      float_image_write(wr_dump,OZ);
      fclose(wr_dump);
      free(filename);
    }
  }
  psgr_solve_system(gr,S,iZ,iW,ref_tab,maxIter,convTol,para,szero,verbose);
  

  if (debug)
{
    
    if (RZ != NULL)
{
      psgr_put_error_into_image(gr,iZ,iW,RZ,OZ);
      char *filename = jsprintf("%s-%02d-oE.fni",out_prefix,level);
      FILE *wr_dump = open_write(filename,FALSE);
      float_image_write(wr_dump,OZ);
      fclose(wr_dump);
      free(filename);
    }

  }
  psgr_put_solution_into_image(gr,iZ,OZ);
  if (debug)
  {   char *filename = jsprintf("%s-%02d-oZ.fni",out_prefix,level);
      FILE *wr_dump = open_write(filename,FALSE);
      float_image_write(wr_dump,OZ);
      fclose(wr_dump);
      free(filename);
    }
    free(ref_tab);
    pst_imgsys_free(S);
  }
