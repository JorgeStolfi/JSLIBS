/* See {float_image_hist.h} */
/* Last edited on 2024-12-24 15:57:09 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#include <vec.h>
#include <bool.h>
#include <jsfile.h>
#include <affirm.h>
#include <float_image.h>

#include <float_image_hist.h>
    
void float_image_hist_build
  ( float_image_t *A,
    int32_t c,
    uint32_t N,
    double *hmin_P, 
    double *hmax_P, 
    double_vec_t *hist_P,
    double_vec_t *cumh_P,
    uint32_t *ngud_P,
    uint32_t *nbad_P
  )
  { 
    demand(N >= 2, "invalid histogram size {N}.");
    int32_t NC, NX, NY;
    float_image_get_size(A, &NC, &NX, &NY);
    demand(c < NC, "invalid channel index {c}");

    /* Get the actual valid sample range {[smin _ smax]}: */
    float smin = +INF, smax = -INF;
    float_image_update_sample_range(A, c, &smin, &smax);
    
    /* Choose the histogram range {[hmin _ hmax]}: */
    float hmin, hmax;
    if (smin > smax)
      { /* No samples: */
        hmin = -1; hmax = +1;
      }
    else
      { double dv; /* Bin width. */
        if (smin == smax)
          { dv = hypot(smin, 1.0)/N;
            uint32_t kmid = N/2;
            hmin = (float)(smin - ((double)kmid + 0.500001)*dv);
            hmax = (float)(hmin + ((double)N + 0.000001)*dv);
          }
        else
          { dv = 1.000001*(smax - smin)/(N-1);
            hmin = (float)(smin - 0.500001*dv);
            hmax = (float)(smax + 0.500001*dv);
          }
      }
    (*hmin_P) = hmin;
    (*hmax_P) = hmax;

    /* Compute the actual bin width: */
    double dv = (hmax - hmin)/((double)N);
    
    /* Allocate the histogram if needed: */
    uint32_t NH = ((hist_P != NULL) || (cumh_P != NULL) ? N : 0);
    double_vec_t hist = double_vec_new(NH);
    for (uint32_t ih = 0; ih < N; ih++) { hist.e[ih] = 0.0; }
    
    /* Scan samples computing {ngud,nbad} and filling {hist}: */
    uint32_t ngud = 0, nbad = 0;
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float s = float_image_get_sample(A, c, x, y);
            if (isfinite(s))
              { ngud++;
                assert((s >= hmin) && (s <= hmax));
                /* Find the bin {k0} that has {s}: */
                int32_t k0 = (int32_t)floor((s - hmin)/dv);
                assert((k0 >= 0) && (k0 < N));
                /* Find the nearest bin {k1} adjacent to {k0}: */
                int32_t k1;
                double ctr = (double)hmin + (k0 + 0.5)*dv;
                double t; /* Fraction to add to bin {k1} instead of {k0}. */
                if (s >= ctr)
                  { k1 = k0 + 1;
                    t = (k1 >= N ? 0.0 : s - ctr);
                  }
                else
                  { k1 = k0 - 1;
                    t = (k1 < 0 ? 0.0 : ctr - s);
                  }
                if (NH > 0)
                  { hist.e[k0] += 1 - t;
                    if (t > 0) { hist.e[k1] += t; }
                  }
              }
            else
              { nbad++; }
          }
      }
    if (ngud_P != NULL) { (*ngud_P) = ngud; }
    if (nbad_P != NULL) { (*nbad_P) = nbad; }

    if (cumh_P != NULL)
      { /* Create the cumulative histogram: */
        assert(NH == N);
        double_vec_t cumh = double_vec_new(NH);
        double sum = 0;
        for (uint32_t ih = 0; ih < N; ih++)
          { sum += hist.e[ih];
            cumh.e[ih] = sum;
          }
        (*cumh_P) = cumh;
      }

    if (hist_P != NULL)
      { assert(NH == N);
        (*hist_P) = hist;
      }
    else if (hist.ne != 0)
      { free(hist.e); }
  }

void float_image_hist_write_file(FILE *wr, double hmin, double hmax, double_vec_t *hist)
  { uint32_t N = hist->ne;
    demand(N >= 2, "invalid histogram length");
    demand(isfinite(hmin) && isfinite(hmax) && (hmin< hmax), "invalid {hmin,hmax}");

    /* Print the histogram: */
    double sum = 0;
    double slo = hmin;
    for (int32_t k = -1; k <= (int32_t)N; k++)
      { double hk, shi;
        if ((k < 0) || (k == (int32_t)N))
          { shi = slo; hk = 0.0; }
        else
          { double r = ((double)k+1)/N;
            shi = (1-r)*hmin + r*hmax;
            hk = hist->e[k]; 
            sum += hk;
          }
        fprintf(wr, "%6d %+23.16e %+23.16e %14.3f %14.2f\n", k, slo, shi, hk, sum);
        slo = shi;
      }
    assert(slo == hmax);
    fflush(wr);
  }
  
void float_image_hist_write_named(char *fname, double hmin, double hmax, double_vec_t *hist)
  { FILE *wr = open_write(fname, TRUE);
    float_image_hist_write_file(wr, hmin, hmax, hist);
    fclose(wr);
  }
