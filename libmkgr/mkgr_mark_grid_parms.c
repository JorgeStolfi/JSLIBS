/* See {mkgr_mark_grid_parms.h} */
/* Last edited on 2020-11-30 13:56:56 by jstolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <affirm.h>
#include <epswr.h>

#include <mkgr_mark_grid_parms.h>
  
void mkgr_mark_grid_parms_show(FILE *wr, char *kind, mkgr_mark_grid_parms_t *parms)
  {
    fprintf(wr, "%s dimensions (mm):\n", kind);
    fprintf(wr, "  cross radius = %.4f\n", parms->rcross);
    fprintf(wr, "  circle radius = %.4f\n", parms->rcircle);
    fprintf(wr, "  dot radius = %.4f\n", parms->rdot);
    fprintf(wr, "  bg circle radius = %.4f\n", parms->rbgcirc);
    fprintf(wr, "  line width = %.4f\n", parms->lwd);
    fprintf(wr, "  spacing = %.4f\n", parms->spacing);
  }

void mkgr_mark_grid_parms_round(mkgr_mark_grid_parms_t *parms, double dots, double unit)
  { 
    double dpi = (25.4*dots)/unit; /* resolution in pixels per inch. */
    mkgr_mark_grid_parms_show(stderr, "nominal", parms);
    parms->rcross = epswr_round_dim_mm(parms->rcross, dpi);   
    parms->rcircle = epswr_round_dim_mm(parms->rcircle, dpi); 
    parms->rdot = epswr_round_dim_mm(parms->rdot, dpi);    
    parms->rbgcirc = epswr_round_dim_mm(parms->rbgcirc, dpi);    
    parms->lwd = 2*epswr_round_dim_mm(parms->lwd/2, dpi); /* Even number of pixels. */
    parms->spacing = 2*epswr_round_dim_mm(parms->spacing/2, dpi); /* Even number of pixels. */
    mkgr_mark_grid_parms_show(stderr, "rounded", parms);
  }
