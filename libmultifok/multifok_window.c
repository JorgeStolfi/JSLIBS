/* See {multifok_window.h}. */
/* Last edited on 2024-12-05 15:15:11 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <vec.h>
#include <jsprintf.h>
#include <wt_table.h>
#include <wt_table_binomial.h>
#include <affirm.h>
#include <bool.h>

#include <multifok_window.h>
  
uint32_t multifok_window_num_samples(uint32_t NW)
  { demand(((NW % 2) == 1) && (NW >= 3), "window size must be odd");
    return NW*NW;
  }
 
double *multifok_window_weights(uint32_t NW, multifok_window_type_t type)
  { double *ws;
    switch (type)
      { case multifok_window_type_BIN:
          ws = multifok_window_weights_binomial(NW);
          break;
        
        case multifok_window_type_GLD:
          ws = multifok_window_weights_golden(NW);
          break;
       
        default:
          affirm(FALSE, "invalid weights type");
      }
    return ws;
  }

double *multifok_window_weights_binomial(uint32_t NW)
  { demand(NW % 2 == 1, "{NW} must be odd");
    uint32_t HW = (uint32_t)(NW-1)/2;
    /* Create a unidimensional weight table: */
    double u[NW]; /* Unidimensional weights. */
    wt_table_binomial_fill(NW, u, NULL);

    /* Normalize so that the central element is 1: */
    double umax = u[HW];
    for (uint32_t ku = 0;  ku < NW; ku++) { u[ku] /= umax; }
    
    /* Now fill the bidimensional table: */
    uint32_t NS = multifok_window_num_samples(NW);
    double *ws = notnull(malloc(NS*sizeof(double)), "no mem");
    for (uint32_t iy = 0;  iy < NW; iy++)
      { for (uint32_t ix = 0;  ix < NW; ix++)
          { double wxy = u[iy]*u[ix];
            uint32_t ks = iy*NW + ix;
            ws[ks] = wxy;
          }
      }

    return ws;
  }
 
double *multifok_window_weights_golden(uint32_t NW)
  { demand(NW % 2 == 1, "{NW} must be odd");
    uint32_t HW = (uint32_t)(NW-1)/2;

    /* Normalize so that the central element is 1: */
    double C = (sqrt(5) - 1)/2;
    
    /* Now fill the bidimensional table: */
    uint32_t NS = multifok_window_num_samples(NW);
    double *ws = notnull(malloc(NS*sizeof(double)), "no mem");
    uint32_t ks = 0;
    for (int32_t iy = -(int32_t)HW; iy <= +(int32_t)HW; iy++)
      { for (int32_t ix = -(int32_t)HW; ix <= +(int32_t)HW; ix++)
          { ws[ks] = C/(C + ix*ix + iy*iy);
            ks++;
          }
      }

    return ws;
  }
char *multifok_window_type_to_text(multifok_window_type_t wType)
  { switch(wType)
      { case multifok_window_type_BIN: return "BIN";
        case multifok_window_type_GLD: return "GLD";
        default: assert(FALSE);
      }
  }

multifok_window_type_t multifok_window_type_from_text(char *name, bool_t fail)
  { if (strcmp(name, "BIN") == 0) { return multifok_window_type_BIN; }
    if (strcmp(name, "GLD") == 0) { return multifok_window_type_GLD; }
    if (fail) { demand(FALSE, "invalid window weighst type name"); } else { return -1; }
  }

void multifok_window_compute_average_gradient_and_deviation
  ( uint32_t NW, 
    double s[], 
    double ws[], 
    double *sAvg_P,
    double *sGrx_P,
    double *sGry_P,
    double *sDev_P
  )
  {
    uint32_t NS = NW*NW;

    /* Compute weighted sample average {sAvg}: */
    double sum_ws = 0;
    double sum_w = 0;
    for (uint32_t ks = 0;  ks < NS; ks++) 
      { sum_ws += ws[ks]*s[ks]; 
        sum_w += ws[ks];
      }
    assert(sum_w > 0);
    double sAvg = sum_ws/sum_w;
    
    /* Compute gradient {(sGrx,sGry)} by weighted dot product with basis {X,Y}: */
    uint32_t HW = (uint32_t)(NW-1)/2;
    double sum_wsx = 0;
    double sum_wsy = 0;
    double sum_wxx = 0;
    double sum_wyy = 0;
    for (uint32_t ks = 0;  ks < NS; ks++) 
      { double xk = (ks % NW) - HW;
        double yk = (ks / NW) - HW;
        double wk = ws[ks];
        double dk = s[ks] - sAvg;
        sum_wsx += wk*dk*xk; 
        sum_wsy += wk*dk*yk; 
        sum_wxx += wk*xk*xk;
        sum_wyy += wk*yk*yk;
      }
    assert(sum_wxx > 0);
    assert(sum_wyy > 0);
    double sGrx = sum_wsx/sum_wxx;
    double sGry = sum_wsy/sum_wyy;
    
    /* Compute deviation {sDev} of residual: */
    double sum_wd2 = 0;
    for (uint32_t ks = 0;  ks < NS; ks++) 
      { double xk = (ks % NW) - HW;
        double yk = (ks / NW) - HW;
        double wk = ws[ks];
        double dk = s[ks] - sAvg - sGrx*xk - sGry*yk;
        sum_wd2 += wk*dk*dk;
      }
    double sDev = sqrt(sum_wd2/sum_w);
    
    (*sAvg_P) = sAvg;
    (*sGrx_P) = sGrx;
    (*sGry_P) = sGry;
    (*sDev_P) = sDev;
  }

void multifok_window_remove_average_and_gradient
  ( uint32_t NW, 
    double s[], 
    double ws[], 
    double sAvg,
    double sGrx,
    double sGry
  )
  {
    uint32_t NS = NW*NW;
    uint32_t HW = (uint32_t)(NW-1)/2;
    for (uint32_t ks = 0;  ks < NS; ks++) 
      { double xk = (ks % NW) - HW;
        double yk = (ks / NW) - HW;
        s[ks] = s[ks] - sAvg - sGrx*xk - sGry*yk;
      }
  }
  
double multifok_window_deviation(uint32_t NW, double s[], double ws[])
  {
    uint32_t NS = NW*NW;
    double sum_w_d2 = 0;
    double sum_w = 0;
    for (uint32_t ks = 0;  ks < NS; ks++) 
      { double d = s[ks];
        sum_w_d2 += ws[ks]*d*d;
      }
    assert(sum_w > 0);
    double sDev = sqrt(sum_w_d2/sum_w);
    return sDev;
  }
 
void multifok_window_normalize_samples
  ( uint32_t NW, 
    double s[], 
    double ws[], 
    double noise, 
    double *sAvg_P,
    double *sGrx_P,
    double *sGry_P,
    double *sDev_P
  )
  {
    uint32_t NS = NW*NW;
    multifok_window_compute_average_gradient_and_deviation(NW, s, ws, sAvg_P, sGrx_P, sGry_P, sDev_P);
    multifok_window_remove_average_and_gradient(NW, s, ws, *sAvg_P, *sGrx_P, *sGry_P);
    noise = fmax(1.0e-200, noise); /* To avoid division of zero by zero. */
    double mag = hypot((*sDev_P), noise);
    for (uint32_t ks = 0;  ks < NS; ks++) { s[ks] /= mag; }
  }

double multifok_window_prod(uint32_t NW, double a[], double b[])
  { uint32_t NS = NW*NW;
    double prod = 0.0;
    for (uint32_t ks = 0;  ks < NS; ks++) 
      { double ak = a[ks];
        double bk = b[ks];
        prod += ak*bk;
      }
    return prod;
  }

double multifok_window_dist_sqr(uint32_t NW, double a[], double b[])
  { uint32_t NS = NW*NW;
    double d2 = 0.0;
    for (uint32_t ks = 0;  ks < NS; ks++) 
      { double dk = a[ks] - b[ks];
        d2 += dk*dk;
      }
    return d2;
  }

char *multifok_window_mop_code(int32_t i)
  { char *u = NULL;
    char s = (i < 0 ? 'm' : (i > 0 ? 'p' : 'o'));
    if (abs(i) <= 1) 
      { u = jsprintf("%c", s); } 
    else 
      { u = jsprintf("%c%d", s, abs(i)); }
    return u;
  }
  
char *multifok_window_sample_name(char *tag, int32_t ix, int32_t iy)
  {
    char *ux = multifok_window_mop_code(ix);
    char *uy = multifok_window_mop_code(iy);
    char *uxy = jsprintf("%s%s%s", tag, ux, uy);
    free(ux); 
    free(uy);
    return uxy;
  }    

void multifok_window_sample_names(uint32_t NW, char *tag, char *sname[])
  {
    assert((NW % 2 ) == 1);
    uint32_t HW = NW/2;
    uint32_t NS = NW*NW;
    
    /* Set {u[i+HW]} to text versions of coordinate {i}, for {i} in {-HW..+HW}: */
    char *u[NW];
    for (int32_t i = -(int32_t)HW; i <= +(int32_t)HW; i++)
      { u[i + (int32_t)HW] = multifok_window_mop_code(i); }
    
    /* Assemble sample names: */
    uint32_t ks = 0; /* Windos sample index. */
    for (int32_t iy = -(int32_t)HW; iy <= +(int32_t)HW; iy++)
      { for (int32_t ix = -(int32_t)HW; ix <= +(int32_t)HW; ix++)
          { char *uxy = jsprintf("%s%s%s", tag, u[ix+(int32_t)HW], u[iy+(int32_t)HW]);
            sname[ks] = uxy;
            ks++;
         }
      }
    assert(ks == NS);
    
    /* Release the auxiliary coordinate names: */
    for (int32_t i = -(int32_t)HW; i <= +(int32_t)HW; i++)
      { free(u[i + (int32_t)HW]); }
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

