/* See {neuromat_eeg_stats.h}. */
/* Last edited on 2013-11-30 04:39:30 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <jsmath.h>
#include <affirm.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_stats.h>

void neuromat_eeg_stats_per_channel
  ( int nt, 
    int nc, 
    double **val, 
    bool_t zeroMean,
    double vavg[], 
    double vvar[],
    double vmin[], 
    double vmax[], 
    int ne, 
    double *minMinP, 
    double *maxMaxP, 
    double *maxDevP
  )
  {
    double vmin_min = +INF;     /* Minimum value among electrode channels. */
    double vmax_max = -INF;     /* Maximum value among electrode channels. */
    double vvar_max = 1.0e-30;  /* Maximum variance among electrode channels. */
    int it, ic;
    for (ic = 0; ic < nc; ic++)
      { /* Compute range {minv,maxv} of channel: */
        double minv = +INF;
        double maxv = -INF;
        for (it = -0; it < nt; it++) 
          { double v = val[it][ic];
            if (v < minv) { minv = v; }
            if (v > maxv) { maxv = v; }
          }
        vmin[ic] = minv;
        vmax[ic] = maxv;
        /* Define the mean {avg}: */
        double avg = NAN;
        if (zeroMean)
          { avg = 0.0; }
        else
          { /* Compute the signal {avg}: */
            double sumv = 0;
            for (it = -0; it < nt; it++) 
              { double v = val[it][ic];
                sumv += v;
              }
            avg = sumv/nt;
          }
        vavg[ic] = avg;
        /* Compute the variance {var} relative to {avg}: */
        double sumd2 = 0;
        for (it = -0; it < nt; it++) 
          { double d = val[it][ic] - avg; 
            sumd2 += d*d;
          }
        double var = sumd2/nt;
        vvar[ic] = var;
        /* Update the global extrema: */
        if (ic < ne)
          { if (vmin[ic] < vmin_min) { vmin_min = vmin[ic]; }
            if (vmax[ic] > vmax_max) { vmax_max = vmax[ic]; }
            if (vvar[ic] > vvar_max) { vvar_max = vvar[ic]; }
          }
      }
    (*minMinP) = vmin_min;
    (*maxMaxP) = vmax_max;
    double vdev_max = sqrt(vvar_max); /* Max standard deviation of all electrode signals. */
    (*maxDevP) = vdev_max;
  }

void neuromat_eeg_stats_per_channel_print
  ( FILE *wr, 
    int nc, 
    char *name[], 
    double vavg[], 
    double vvar[], 
    double vmin[], 
    double vmax[] 
  )
  {
    int ic;
    for (ic = 0; ic < nc; ic++) 
      { fprintf(wr, "channel %3d = %-4s", ic, name[ic]); 
        fprintf(wr, "  avg = %+14.5f  dev = %14.5f", vavg[ic], sqrt(vvar[ic])); 
        fprintf(wr, "  min = %+14.5f  max = %+14.5f\n", vmin[ic], vmax[ic]); 
      }
  }
  
double *neuromat_eeg_stats_covariance_matrix(int nt, int ne, double **val, double vavg[])
  {
    double *Cv = notnull(malloc(ne*ne*sizeof(double)), "no mem");
    int it, ie, je;
    for (ie = 0; ie < ne; ie++)
      { for (je = 0; je <= ie; je++)
          { double sum2 = 0;
            for (it = 0; it < nt; it++)
              { double vi = val[it][ie] - vavg[ie];
                double vj = val[it][je] - vavg[je];
                sum2 += vi*vj;
              }
            double Cij = sum2/nt;
            Cv[ie*ne + je] = Cij;
            Cv[je*ne + ie] = Cij;
          }
      }
    return Cv;
  }
