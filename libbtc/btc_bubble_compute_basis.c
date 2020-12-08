/* See {btc_bubble_compute_basis.h} */
/* Last edited on 2015-04-30 00:00:38 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>

#include <btc_bubble_t.h>
#include <btc_bubble_compute_basis.h>
#include <btc_price_series_smooth.h>

void btc_bubble_compute_basis(int nd, int nb, btc_bubble_t bp[], int hrad, double bval[])
  {
    bool_t debug = TRUE;
    
    if (debug) { fprintf(stderr, "computing the bubble basis\n"); }
    
    double* bvali = notnull(malloc(nd*sizeof(double)), "no mem"); /* Values of a single bubble, raw. */
    double* bvalo = notnull(malloc(nd*sizeof(double)), "no mem"); /* Values of a single bubble, smoothed. */
    int id, jb;
    for (jb = 0; jb < nb; jb++)
      { /* Compute the values of bubble {jb}: */
        btc_bubble_t* bpj = &(bp[jb]);
        double log_rtup = log(bpj->rt_up);
        double log_rtdn = log(bpj->rt_dn);
        int id_fin_up = bpj->id_fin_up;
        int id_ini_dn = id_fin_up + bpj->wd_plat;
        for (id = 0; id < nd; id++)
          { double val;
            if (id < id_fin_up)
              { val = exp(log_rtup*(id - id_fin_up)); }
            else if (id > id_ini_dn)
              { val = exp(log_rtdn*(id - id_ini_dn)); }
            else
              { val = 1.0; }
            bvali[id] = val;
          }
          
        /* Smooth the values with Hann window: */
        btc_price_series_smooth(nd, bvali, hrad, bvalo);
        
        /* Store in big array: */
        for (id = 0; id < nd; id++) { bval[nb*id + jb] = bvalo[id]; }
      }
    free(bvali);
    free(bvalo);
  }
    

