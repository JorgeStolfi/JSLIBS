/* See {btc_bubble_eval_rms_log_error.h} */
/* Last edited on 2015-04-22 20:46:49 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>

#include <btc_bubble_t.h>

#include <btc_bubble_eval_rms_log_error.h>

double btc_bubble_eval_rms_log_error
  ( int nd, 
    double ap[],
    int id_ini,
    int id_fin,
    int nb,
    btc_bubble_t bp[],
    double bval[]
  )
  {
    bool_t verbose = FALSE;
    
    demand((0 <= id_ini) && (id_ini<=id_fin) && (id_fin < nd), "invalid eval date range");
    
    double log10 = log(10.0);
    double sum_d2 = 0.0;
    int ngood = 0; /* Number of non-zero prices in given series. */
    int id;
    for (id = id_ini; id <= id_fin; id++)
      { if (ap[id] > 0)
          { double mp = 0.0;
            int jb;
            for (jb = 0; jb < nb; jb++)  { mp += bp[jb].coef * bval[nb*id + jb]; }
            if (verbose) { fprintf(stderr, "      mp[%02d] = %25.16e  ap[%02d] = %25.16e\n", id, mp, id, ap[id]); } 
            double d = (log(mp) - log(ap[id]))/log10;
            sum_d2 += d*d;
            ngood ++;
          }
      }
    assert(ngood > 0);
    return sum_d2/ngood;
  }
