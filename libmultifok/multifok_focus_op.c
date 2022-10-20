/* See {multifok_focus_op.h}. */
/* Last edited on 2022-10-20 06:27:05 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <vec.h>
#include <affirm.h>

#include <multifok_focus_op.h>

double multifok_focus_op_score_simple(int32_t NW, double v[], double *phi[], double w[], double noise)
  { 
    int32_t NS = multifok_focus_op_num_samples(NW);
    double c[NS];
    for (int32_t i = 0; i < NS; i++) { c[i] = multifok_focus_op_prod(NW, v, phi[i], w); }
    double m2 = c[0]*c[0]; /* Mean value squared. */
    double g2 = c[1]*c[1] + c[2]*c[2]; /* Gradient modulus squared. */
    /* double foc = (g2 + 1.0e-100)/(g2 + noise + 1.0e-100); */
    double foc = (g2 + 1.0e-100)/(g2 + m2 + noise + 1.0e-100);
    /* double p2 = c[3]*c[3] + c[4]*c[4] + c[5]*c[5]; */ /* Second-deg moduli squared. */
    return foc;
  }
  
int32_t multifok_focus_op_num_samples(int32_t NW)
  { demand((NW % 2) == 1, "window size must be odd");
    return NW*NW;
  }
 
double *multifok_focus_op_prod_weights(int32_t NW)
  { assert(NW == 3); /* For now. */
    int32_t NS = multifok_focus_op_num_samples(NW);
    double *w = notnull(malloc(NS*sizeof(double)), "no mem");
    for (int32_t dy = -1; dy <= +1; dy++)
      { for (int32_t dx = -1; dx <= +1; dx++)
          { double wy = (dy == 0 ? 0.50 : 0.25);
            double wx = (dx == 0 ? 0.50 : 0.25);
            w[(dy+1)*3 + (dx+1)] = wx*wy;
          }
      }
    return w;
  }

double multifok_focus_op_prod(int32_t NW, double x[], double y[], double w[])
  { int32_t NS = multifok_focus_op_num_samples(NW);
    double prod = 0.0;
    for (int32_t k = 0; k < NS; k++) 
      { double xk = x[k];
        double yk = y[k];
        double wk = w[k];
        prod += wk*xk*yk;
      }
    return prod;
  }

double multifok_focus_op_dist_sqr(int32_t NW, double x[], double y[], double w[])
  { int32_t NS = multifok_focus_op_num_samples(NW);
    double d2 = 0.0;
    for (int32_t k = 0; k < NS; k++) 
      { double dk = x[k] - y[k];
        double wk = w[k];
        d2 += wk*dk*dk;
      }
    return d2;
  }

double **multifok_focus_op_basis(int32_t NW)
  { assert(NW == 3); /* For now. */
    int32_t NS = multifok_focus_op_num_samples(NW);
    
    double **phi = notnull(malloc(NS*sizeof(double*)), "no mem"); 
    for (int32_t i = 0; i < NS; i++) { phi[i] = notnull(malloc(NS*sizeof(double)), "no mem"); }
    
    double sqrt2 = M_SQRT2;
    
    multifok_focus_op_set_samples_3x3(phi[0], 1.000,  +1.0f, +1.0f, +1.0f,  +1.0f, +1.0f, +1.0f,  +1.0f, +1.0f, +1.0f ); /* Mean. */
    multifok_focus_op_set_samples_3x3(phi[1], sqrt2,  -1.0f, 00.0f, +1.0f,  -1.0f, 00.0f, +1.0f,  -1.0f, 00.0f, +1.0f ); /* DX. */
    multifok_focus_op_set_samples_3x3(phi[2], sqrt2,  -1.0f, -1.0f, -1.0f,  00.0f, 00.0f, 00.0f,  +1.0f, +1.0f, +1.0f ); /* DY. */
    multifok_focus_op_set_samples_3x3(phi[3], 1.000,  -1.0f, +1.0f, -1.0f,  -1.0f, +1.0f, -1.0f,  -1.0f, +1.0f, -1.0f ); /* DXX. */
    multifok_focus_op_set_samples_3x3(phi[4], 1.000,  -1.0f, -1.0f, -1.0f,  +1.0f, +1.0f, +1.0f,  -1.0f, -1.0f, -1.0f ); /* DYY. */
    multifok_focus_op_set_samples_3x3(phi[5], 2.000,  -1.0f, 00.0f, +1.0f,  00.0f, 00.0f, 00.0f,  +1.0f, 00.0f, -1.0f ); /* DXY. */
    multifok_focus_op_set_samples_3x3(phi[6], sqrt2,  -1.0f, +1.0f, -1.0f,  00.0f, 00.0f, 00.0f,  +1.0f, -1.0f, +1.0f ); /* Saddle3. */
    multifok_focus_op_set_samples_3x3(phi[7], sqrt2,  +1.0f, 00.0f, -1.0f,  -1.0f, 00.0f, +1.0f,  +1.0f, 00.0f, -1.0f ); /* Cosaddle3. */
    multifok_focus_op_set_samples_3x3(phi[8], 1.000,  +1.0f, -1.0f, +1.0f,  -1.0f, +1.0f, -1.0f,  +1.0f, -1.0f, +1.0f ); /* Checker. */
    
    return phi;
  }

void multifok_focus_op_orthize(int32_t NW, int32_t k, double *phi[], double w[])
  { int32_t NS = multifok_focus_op_num_samples(NW);
    /* Make orthogonal:*/
    for (int32_t i = 0; i < k; i++)
      { double ci = multifok_focus_op_prod(NW, phi[k], phi[i], w);
        for (int32_t j = 0; j < NS; j++)
          { phi[k][j] = phi[k][j] - ci*phi[i][j]; }
      }
    /* Normalize: */
    double d2 = multifok_focus_op_prod(NW, phi[k], phi[k], w);
    double d = sqrt(d2);
    for (int32_t j = 0; j < NS; j++)
      { phi[k][j] = phi[k][j]/d; }
  }

void multifok_focus_op_set_samples_3x3
  ( double x[], 
    double scale,
    float x00, float x01, float x02, 
    float x10, float x11, float x12, 
    float x20, float x21, float x22 
  )
  { x[0] = scale * ((double)x00);
    x[1] = scale * ((double)x01);
    x[2] = scale * ((double)x02);
    x[3] = scale * ((double)x10);
    x[4] = scale * ((double)x11);
    x[5] = scale * ((double)x12);
    x[6] = scale * ((double)x20);
    x[7] = scale * ((double)x21);
    x[8] = scale * ((double)x22);
  }

void multifok_focus_op_basis_check(int32_t NW, double *phi[], double w[])
  { int32_t NS = multifok_focus_op_num_samples(NW);

    /* Adjust {phi[8]}: */
    /* multifok_focus_op_orthize(NW, 8, phi, w); */
    
    /* Print basis: */
    fprintf(stderr, "--- basis ---------------------------------------------\n");
    fprintf(stderr, "%3s", "");
    for (int32_t j = 0; j < NS; j++) { fprintf(stderr, " %10d", j); }
    fprintf(stderr, "\n");
    for (int32_t i = 0; i < NS; i++)
      { fprintf(stderr, "%3d", i);
        for (int32_t j = 0; j < NS; j++)
          { double pij = phi[i][j];
            fprintf(stderr, " %+10.6f", pij);
          }
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "-------------------------------------------------------\n");
    
    /* Print basis products: */
    fprintf(stderr, "--- momentum matrix -----------------------------------\n");
    fprintf(stderr, "%3s", "");
    for (int32_t j = 0; j < NS; j++) { fprintf(stderr, " %10d", j); }
    fprintf(stderr, "\n");
    for (int32_t i = 0; i < NS; i++)
      { fprintf(stderr, "%3d", i);
        for (int32_t j = 0; j <= i; j++)
          { double pij = multifok_focus_op_prod(NW, phi[i], phi[j], w);
            fprintf(stderr, " %+10.6f", pij);
          }
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "-------------------------------------------------------\n");
  }

void multifok_focus_op_remap_samples(int32_t NW, double v[], double *phi[], double w[], double c[])
  { int32_t NS = multifok_focus_op_num_samples(NW);
    for (int32_t i = 0; i < NS; i++) 
      { c[i] = multifok_focus_op_prod(NW, v, phi[i], w); }
  }
    
void multifok_focus_op_check(int32_t NW)
  { int32_t NS = multifok_focus_op_num_samples(NW);

    double *w = multifok_focus_op_prod_weights(NW);
    double **phi = multifok_focus_op_basis(NW);
    multifok_focus_op_basis_check(NW, phi, w);
    free(w);
    for (int32_t i = 0; i < NS; i++) { free(phi[i]); }
    free(phi);
  } 

#define multifok_focus_op_C_COPYRIGHT \
    "© 2017 by the State University of Campinas (UNICAMP)"

