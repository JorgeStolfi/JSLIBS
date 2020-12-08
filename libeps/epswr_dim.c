/* See epswr_dim.h */
/* Last edited on 2020-10-27 17:02:45 by jstolfi */

#define _GNU_SOURCE
#include <math.h>

#include <affirm.h>

#include <epswr.h>
#include <epswr_dim.h>
#include <epswr_dev_dim.h>

void epswr_dim_linear
  ( epswr_figure_t *eps, 
    double xa, double ya,
    double xb, double yb,
    double *dabP,
    double agap,
    double bgap, 
    double elen,
    bool_t inner,
    double dpos, double dlen, 
    double hoff, double voff,
    double *xrP, double *yrP,
    double *rotP
  )
  {
    /* Compute distance in Client coordinates: */
    double dab = hypot(xb-xa, yb-ya);
    demand(dab > 0.0, "coincident points");
    if (dabP != NULL) { (*dabP) = dab; }
    
    /* Convert oDevice coordinates: */
    double psxa; epswr_x_to_h_coord(eps, xa, &(psxa));
    double psya; epswr_y_to_v_coord(eps, ya, &(psya));
    double psxb; epswr_x_to_h_coord(eps, xb, &(psxb));
    double psyb; epswr_y_to_v_coord(eps, yb, &(psyb));
    
    double psagap; epswr_xy_to_hv_dist(eps, agap, &(psagap));
    double psbgap; epswr_xy_to_hv_dist(eps, bgap, &(psbgap));

    double pselen = elen * epswr_mm;
    
    double psdpos = dpos * epswr_mm;
    double psdlen = dlen * epswr_mm;
    double pshoff = hoff * epswr_mm;
    double psvoff = voff * epswr_mm;
    if (dabP != NULL) { (*dabP) = dab; }
    double psxr, psyr;
    epswr_dev_dim_linear
      ( eps, psxa, psya, psxb, psyb,
        psagap, psbgap, pselen, 
        inner,  psdpos, psdlen,
        pshoff, psvoff, &psxr, &psyr, rotP
      );
    if (xrP != NULL) { epswr_h_to_x_coord(eps, psxr, xrP); }
    if (yrP != NULL) { epswr_v_to_y_coord(eps, psyr, yrP); }
  }
  
