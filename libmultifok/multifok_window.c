/* See {multifok_window.h}. */
/* Last edited on 2023-01-30 06:56:14 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <vec.h>
#include <wt_table.h>
#include <affirm.h>
#include <bool.h>

#include <multifok_window.h>
  
int32_t multifok_window_num_samples(int32_t NW)
  { demand(((NW % 2) == 1) && (NW >= 3), "window size must be odd");
    return NW*NW;
  }
 
double *multifok_window_sample_weights(int32_t NW)
  { demand(NW % 2 == 1, "{NW} must be odd");
    int32_t HW = (NW-1)/2;
    /* Create a unidimensional weight table: */
    double u[NW]; /* Unidimensional weights. */
    bool_t normalize = FALSE;
    wt_table_fill_binomial(NW, u, normalize);

    /* Normalize so that the central element is 1: */
    double umax = u[HW];
    for (int32_t ku = 0; ku < NW; ku++) { u[ku] /= umax; }
    
    /* Now fill the bidimensional table: */
    int32_t NS = multifok_window_num_samples(NW);
    double *ws = notnull(malloc(NS*sizeof(double)), "no mem");
    for (int32_t iy = 0; iy < NW; iy++)
      { for (int32_t ix = 0; ix < NW; ix++)
          { double wxy = u[iy]*u[ix];
            int32_t ks = iy*NW + ix;
            ws[ks] = wxy;
          }
      }

    return ws;
  }

void multifok_window_normalize_samples
  ( int32_t NW, 
    double s[], 
    double ws[], 
    double noise, 
    double *avg_P,
    double *dev_P
  )
  {
    int32_t NS = NW*NW;
    /* Compute weighted sample average: */
    double sumw_s = 0;
    double sum_w = 0;
    for (int32_t ks = 0; ks < NS; ks++) 
      { sumw_s += ws[ks]*s[ks]; 
        sum_w += ws[ks];
      }
    assert(sum_w > 0);
    double avg = sumw_s/sum_w;
    /* Compute the weighted deviation: */
    double sum_w_d2 = 0;
    for (int32_t ks = 0; ks < NS; ks++) 
      { double d = s[ks] - avg;
        sum_w_d2 += ws[ks]*d*d;
      }
    double dev = sqrt(sum_w_d2/sum_w);
    /* Adjust {dev} for the noise level: */
    noise = fmax(1.0e-200, noise); /* To avoid division of zero by zero. */
    double mag = hypot(dev, noise);
    /* Normalize samples: */
    for (int32_t ks = 0; ks < NS; ks++) 
      { s[ks] = (s[ks] - avg) / mag; }
    (*avg_P) = avg;
    (*dev_P) = dev;
  }

double multifok_window_prod(int32_t NW, double a[], double b[])
  { int32_t NS = NW*NW;
    double prod = 0.0;
    for (int32_t ks = 0; ks < NS; ks++) 
      { double ak = a[ks];
        double bk = b[ks];
        prod += ak*bk;
      }
    return prod;
  }

double multifok_window_dist_sqr(int32_t NW, double a[], double b[])
  { int32_t NS = NW*NW;
    double d2 = 0.0;
    for (int32_t ks = 0; ks < NS; ks++) 
      { double dk = a[ks] - b[ks];
        d2 += dk*dk;
      }
    return d2;
  }

char *multifok_window_mop_code(int32_t i)
  { char *u = NULL;
    char s = (i < 0 ? 'm' : (i > 0 ? 'p' : 'o'));
    if (abs(i) <= 1) { asprintf(&u, "%c", s); } else { asprintf(&u, "%c%d", s, abs(i)); }
    return u;
  }
  
char *multifok_window_sample_name(char *tag, int32_t ix, int32_t iy)
  {
    char *ux = multifok_window_mop_code(ix);
    char *uy = multifok_window_mop_code(iy);
    char *uxy = NULL;
    asprintf(&uxy, "%s%s%s", tag, ux, uy);
    free(ux); 
    free(uy);
    return uxy;
  }    

void multifok_window_sample_names(int32_t NW, char *tag, char *sname[])
  {
    assert((NW % 2 ) == 1);
    int32_t HW = NW/2;
    int32_t NS = NW*NW;
    
    /* Set {u[i+HW]} to text versions of coordinate {i}, for {i} in {-HW..+HW}: */
    char *u[NW];
    for (int32_t i = -HW; i <= +HW; i++)
      { u[i + HW] = multifok_window_mop_code(i); }
    
    /* Assemble sample names: */
    int32_t ks = 0; /* Windos sample index. */
    for (int32_t iy = -HW; iy <= +HW; iy++)
      { for (int32_t ix = -HW; ix <= +HW; ix++)
          { char *uxy = NULL;
            asprintf(&uxy, "%s%s%s", tag, u[ix+HW], u[iy+HW]);
            sname[ks] = uxy;
            ks++;
         }
      }
    assert(ks == NS);
    
    /* Release the auxiliary coordinate names: */
    for (int32_t i = -HW; i <= +HW; i++){ free(u[i + HW]); }
  }

void multifok_window_set_samples_3x3
  ( double s[], 
    double scale,
    double s00, double s01, double s02, 
    double s10, double s11, double s12, 
    double s20, double s21, double s22 
  )
  { s[0] = scale * s00;
    s[1] = scale * s01;
    s[2] = scale * s02;
    s[3] = scale * s10;
    s[4] = scale * s11;
    s[5] = scale * s12;
    s[6] = scale * s20;
    s[7] = scale * s21;
    s[8] = scale * s22;
  }

#define multifok_window_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

