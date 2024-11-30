/* See {neuromat_eeg_channel_stats.h}. */
/* Last edited on 2021-09-01 20:25:04 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <jsmath.h>
#include <affirm.h>

#include <neuromat_eeg.h>
#include <neuromat_eeg_channel_stats.h>

void neuromat_eeg_channel_stats_extreme
  ( int32_t ne, 
    neuromat_eeg_channel_stats_t st[],
    neuromat_eeg_channel_stats_t *stg
  );
  /* Merges the data {st[0..ne-1]} into the overall 
    electrode channel summary
    record {*stg}, as described for
    {neuromat_eeg_channel_stats_gather_all}. */

void neuromat_eeg_channel_stats_extreme_print
  ( FILE *wr, 
    int32_t indent,
    int32_t ne,
    char *name[],
    neuromat_eeg_channel_stats_t *stg
  );
  /* Prints the overall 
    electrode channel summary
    record {*stg}, as described for
    {neuromat_eeg_channel_stats_print_all}. */
  
neuromat_eeg_channel_stats_t *neuromat_eeg_channel_stats_new(int32_t nc)
  { neuromat_eeg_channel_stats_t *st;
    st = notnull(malloc(nc*sizeof(neuromat_eeg_channel_stats_t)), "no mem");
    for (uint32_t ic = 0;  ic < nc; ic++) 
      { neuromat_eeg_channel_stats_clear(&(st[ic])); }
    return st;
  }
  
void neuromat_eeg_channel_stats_clear(neuromat_eeg_channel_stats_t *st)
  { st->num = 0;
    st->twt = 0.0;
    st->min = +INF; st->icmin = -1; st->itmin = -1;
    st->max = -INF; st->icmax = -1; st->itmax = -1;
    st->avg = NAN;
    st->msq = NAN;
    st->rms = NAN;
    st->var = NAN;
    st->dev = NAN;
  }

double neuromat_eeg_channel_stats_avg
  ( int32_t nt, 
    int32_t nc, 
    double **val,
    double wt[],
    int32_t ic
  )
  { demand ((ic >= 0) && (ic < nc), "invalid channel index {ic}");
    double sumw = 0.0;
    double sumwv = 0;
    for (uint32_t it = 0;  it < nt; it++) 
      { double v = val[it][ic];
        double w = (wt == NULL ? 1.0 : wt[it]);
        demand((! isnan(w)) && (w >= 0.0) && (w < +INF), "invalid weight");
        if (w > 0.0)
          { sumw += w;
            sumwv += w*v;
          }
      }
    if (sumw == 0.0)
      { return 0.0; }
    else
      { return sumwv/sumw; }
  }

void neuromat_eeg_channel_stats_gather
  ( int32_t nt, 
    int32_t nc, 
    double **val,
    double wt[],
    double eps,
    int32_t ic,
    neuromat_eeg_channel_stats_t *st
  )
  { demand ((! isnan(eps)) && (eps >= 0.0) && (eps < +INF), "invalid {eps}");
    demand ((ic >= 0) && (ic < nc), "invalid channel index {ic}");
    /* Compute {num,twt,min,icmin,itmin,max,icmax,itmax} of channel, and sum and sum of squares: */
    st->num = 0;
    st->min = +INF; st->icmin = -1; st->itmin = -1;
    st->max = -INF; st->icmax = -1; st->itmax = -1;
    double sumw = 0.0;
    double sumwv = 0;
    double sumwv2 = 0;
    for (uint32_t it = 0;  it < nt; it++) 
      { double v = val[it][ic];
        double w = (wt == NULL ? 1.0 : wt[it]);
        demand((! isnan(w)) && (w >= 0.0) && (w < +INF), "invalid weight");
        if (w > 0.0)
          { st->num++;
            double vlo = v - 3*eps;
            if (vlo < st->min) { st->min = vlo; st->icmin = ic; st->itmin = it; }
            double vhi = v + 3*eps;
            if (vhi > st->max) { st->max = vhi; st->icmax = ic; st->itmax = it; }
            sumw += w;
            sumwv += w*v;
            sumwv2 += w*(v*v + eps*eps);
          }
      }
    st->twt = sumw;
    /* Compute {avg,msq,rms,var,dev} of channel: */
    if (sumw == 0.0)
      { neuromat_eeg_channel_stats_clear(st); }
    else
      { st->avg = sumwv/sumw;
        st->msq = sumwv2/sumw;
        st->rms = sqrt(st->msq);
        /* Compute {var,dev} of channel: */
        double sumwd2 = 0.0;
        for (uint32_t it = 0;  it < nt; it++) 
          { double v = val[it][ic];
            double w = (wt == NULL ? 1.0 : wt[it]);
            if (w > 0.0) 
              { double d = v - st->avg; 
                sumwd2 += w*(d*d + eps*eps);
              }
          }   
        st->var = sumwd2/sumw;
        st->dev = sqrt(st->var);
      }
  }

void neuromat_eeg_channel_stats_gather_all
  ( int32_t nt, 
    int32_t nc, 
    double **val,
    double wt[],
    double eps,
    neuromat_eeg_channel_stats_t st[],
    int32_t ne, 
    neuromat_eeg_channel_stats_t *stg
  )
  {
    for (uint32_t ic = 0;  ic < nc; ic++)
      { neuromat_eeg_channel_stats_gather(nt, nc, val, wt, eps, ic, &(st[ic])); }
    if (stg != NULL)
      { neuromat_eeg_channel_stats_extreme(ne, st, stg); }
  }
  
void neuromat_eeg_channel_stats_extreme
  ( int32_t ne, 
    neuromat_eeg_channel_stats_t st[],
    neuromat_eeg_channel_stats_t *stg
  )
  {
    int32_t min_num = INT32_MAX;
    double min_twt = +INF;
    double min_min = +INF; int32_t min_icmin = -1; int32_t min_itmin = -1;
    double max_max = -INF; int32_t max_icmax = -1; int32_t max_itmax = -1;
    double max_var = 0.0;
    double max_dev = 0.0;
    double max_msq = 0.0;
    double max_rms = 0.0;
    for (uint32_t ie = 0;  ie < ne; ie++)
      { neuromat_eeg_channel_stats_t *sti = &(st[ie]);
        if (sti->num < min_num) { min_num = sti->num; }
        if (sti->twt < min_twt) { min_twt = sti->twt; }
        if (sti->min < min_min) { min_min = sti->min; min_icmin = sti->icmin; min_itmin = sti->itmin; }
        if (sti->max > max_max) { max_max = sti->max; max_icmax = sti->icmax; max_itmax = sti->itmax; }
        if (sti->var > max_var) { max_var = sti->var; }
        if (sti->dev > max_dev) { max_dev = sti->dev; }
        if (sti->msq > max_msq) { max_msq = sti->msq; }
        if (sti->rms > max_rms) { max_rms = sti->rms; }
      }
    if ((min_twt == 0.0) || (min_num == 0))
      { assert((min_twt == 0.0) && (min_num == 0));
        neuromat_eeg_channel_stats_clear(stg);
      }
    else
      { stg->num = min_num;
        stg->twt = min_twt;
        stg->min = min_min; stg->icmin = min_icmin; stg->itmin = min_itmin;
        stg->max = max_max; stg->icmax = max_icmax; stg->itmax = max_itmax;
        stg->avg = NAN;
        stg->var = max_var;
        stg->dev = max_dev;
        stg->msq = max_msq;
        stg->rms = max_rms;
      }
  }
  
void neuromat_eeg_channel_stats_print
  ( FILE *wr, 
    int32_t indent,
    int32_t ic,
    char *name,
    bool_t pnum,
    neuromat_eeg_channel_stats_t *st
  )
  {
    if (indent > 0) { fprintf(wr, "%*s", indent, ""); }
    fprintf(wr, "channel %3d = %-4s", ic, name); 
    if (pnum) { fprintf(wr, "  num = %6d  twt = %14.8f", st->num, st->twt); } 
    fprintf(wr, "  avg = %+14.5f  dev = %14.5f  rms = %14.5f", st->avg, st->dev, st->rms); 
    fprintf(wr, "  min = %+14.5f at frame %5d", st->min, st->itmin); 
    fprintf(wr, "  max = %+14.5f at frame %5d", st->max, st->itmax); 
    fprintf(wr, "\n"); 
  }

void neuromat_eeg_channel_stats_extreme_print
  ( FILE *wr, 
    int32_t indent,
    int32_t ne,
    char *name[],
    neuromat_eeg_channel_stats_t *stg
  )
  {
    char *namemin = (stg->icmin < 0 ? "????" : name[stg->icmin]);
    fprintf(wr, "%*smin value %+14.5f at channel %3d = %s frame %5d\n", indent, "", stg->min, stg->icmin, namemin, stg->itmin);
    char *namemax = (stg->icmax < 0 ? "????" : name[stg->icmax]);
    fprintf(wr, "%*smax value %+14.5f at channel %3d = %s frame %5d\n", indent, "", stg->max, stg->icmax, namemax, stg->itmax);
    fprintf(wr, "%*smax dev = %14.5f  max rms = %14.5f\n", indent, "", stg->dev, stg->rms); 
  }

void neuromat_eeg_channel_stats_print_all
  ( FILE *wr, 
    int32_t indent,
    int32_t nc, 
    char *name[],
    bool_t pnum,
    neuromat_eeg_channel_stats_t st[],
    int32_t ne,
    neuromat_eeg_channel_stats_t *stg
  )
  { for (uint32_t ic = 0;  ic < nc; ic++) 
      { neuromat_eeg_channel_stats_print(wr, indent, ic, name[ic], pnum, &(st[ic]));  }
    if ((stg != NULL) && (ne > 0))
      { fprintf(wr, "\n");
        neuromat_eeg_channel_stats_extreme_print(wr, indent, ne, name, stg);
      }
  }

double *neuromat_eeg_channel_stats_covariance_matrix
  ( int32_t nt, 
    int32_t ne, 
    double **val, 
    double vshift[],
    double wt[]
  )
  {
    double *Cv = notnull(malloc(ne*ne*sizeof(double)), "no mem");
    for (uint32_t ije = 0;  ije < ne*ne; ije++) { Cv[ije] = 0.0; }
    neuromat_eeg_channel_stats_accum_covariance_matrix(nt, ne, val, vshift, wt, 1.0, 0.0, Cv);
    return Cv;
  }
    
void neuromat_eeg_channel_stats_accum_covariance_matrix
  ( int32_t nt,      /* Number of frames. */
    int32_t ne,      /* Number of electrodes. */
    double **val,    /* The samples per frame and electrode. */
    double vshift[], /* Values to subtract from each channel. */
    double wt[],     /* Weight of each frame. */
    double cnew,     /* Coefficient for newly computed matrix. */
    double cold,     /* Coefficient for old matrix. */
    double *Cv       /* (IN/OUT) Covariance matrix. */
  )
  {
    for (uint32_t ie = 0;  ie < ne; ie++)
      { for (uint32_t je = 0;  je <= ie; je++)
          { double sumwvivj = 0;
            double sumw = 0.0;
            for (uint32_t it = 0;  it < nt; it++)
              { double vi = val[it][ie] - vshift[ie];
                double vj = val[it][je] - vshift[je];
                double w = (wt == NULL ? 1.0 : wt[it]);
                demand((! isnan(w)) && (w >= 0.0) && (w < +INF), "invalid weight");
                if (w > 0.0)
                  { sumw += w;
                    sumwvivj += w*vi*vj;
                  }
              }
            double Cij = (sumw > 0.0 ? sumwvivj/sumw : 0.0);
            int32_t ije = ie*ne + je;
            Cv[ije] = cold*Cv[ije] + cnew*Cij;
            if (ie != je)
              { int32_t jie = je*ne + ie;
                Cv[jie] = cold*Cv[jie] + cnew*Cij;
              }
          }
      }
  }
