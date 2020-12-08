/* See {mkgr_mark_grid.h} */
/* Last edited on 2020-12-08 00:00:11 by jstolfi */

#define _GNU_SOURCE
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>

#include <bool.h>
#include <vec.h>
#include <r2.h>
#include <i2.h>
#include <affirm.h>
#include <frgb.h>
#include <epswr.h>

#include <mkgr_mark.h>
#include <mkgr_mark_grid.h>

mkgr_mark_grid_t *mkgr_mark_grid_new(void)
  {
    mkgr_mark_grid_t *gr = (mkgr_mark_grid_t *)notnull(malloc(sizeof(mkgr_mark_grid_t)), "no mem");
    gr->NM = 0;
    gr->mark = mkgr_mark_vec_new(0);
    /* Empty bounding box: */
    gr->pMin = (r2_t){{ +INF, +INF }};
    gr->pMax = (r2_t){{ -INF, -INF }};
    return gr;
  }
  
void mkgr_mark_grid_free(mkgr_mark_grid_t *gr)
  { 
    free(gr->mark.e);
    free(gr);
  }

void mkgr_mark_grid_append_mark(mkgr_mark_grid_t *gr, mkgr_mark_t *mk)
  {
    double trad = mk->rad + mk->lwd/2; /* Total radius including outline. */
    if (trad > 0.0)
      { int32_t km = gr->NM;
        mkgr_mark_vec_expand(&(gr->mark), km);
        gr->mark.e[km] = (*mk);
        /* Update the bounding box: */
        for (int32_t j = 0; j < 2; j++)
          { /* Interval of mark's nominal bounding box projected on axis {j}: */
            double ctrj = mk->ctr.c[j];
            double ploj = ctrj - trad;
            double phij = ctrj + trad;
            /* Update the box along axis {j}: */
            if (ploj < gr->pMin.c[j]) { gr->pMin.c[j] = ploj; }
            if (phij > gr->pMax.c[j]) { gr->pMax.c[j] = phij; }
          }
        gr->NM++;
      }
  }

void mkgr_mark_grid_append_marks
  ( mkgr_mark_grid_t *gr,
    i2_t iMin, 
    i2_t iMax, 
    mkgr_mark_grid_def_proc_t *defmark
  )
  {
    
    /* Generate the marks: */
    for (int32_t iy = iMin.c[1]; iy <= iMax.c[1]; iy++)
      { for (int32_t ix = iMin.c[0]; ix <= iMax.c[0]; ix++) 
          { mkgr_mark_t mk = defmark(ix, iy);
            mkgr_mark_grid_append_mark(gr, &mk);
          }
      }
  }
  
void mkgr_mark_grid_write(FILE *wr, mkgr_mark_grid_t *gr)
  {
    int32_t NM = gr->NM;
    fprintf(wr,"NM = %d\n", gr->NM);
    fprintf(wr,"pMin = %9.6f %9.6f\n", gr->pMin.c[0], gr->pMin.c[1]);
    fprintf(wr,"pMax = %9.6f %9.6f\n", gr->pMax.c[0], gr->pMax.c[1]);
    for (int32_t km = 0; km < NM; km++)
      { mkgr_mark_t *mk = &(gr->mark.e[km]);
        bool_t cross = mk->cross;
        double xc = mk->ctr.c[0];
        double yc = mk->ctr.c[1];
        double rad = mk->rad;
        double lwd = mk->lwd;
        double ang = mk->ang;
        fprintf(wr, "%2d %+11.6f %+11.6f  %d  %9.6f  %9.6f  %9.6f\n", km, xc, yc, cross, rad, lwd, ang);
      }
    fflush(wr);
  }
