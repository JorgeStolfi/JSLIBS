/* See {btc_price_series_local_fit_and_eval.h} */
/* Last edited on 2024-12-05 10:23:43 by stolfi */

#include <stdio.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>

#include <btc_price_series_local_fit_and_eval.h>

double btc_price_series_local_fit_and_eval(int deg, int hrad, double val[], double wht[])
  { 
    bool_t debug = FALSE;

    demand(deg == 1, "degree not implemented");
    
    /* Least squares matrix: */
    double A00 = 0; /* {<f0|f0> = SUM{wht[hrad+j]} */
    double A01 = 0; /* {<f0|f1> = SUM{wht[hrad+j]*j} */
    double A11 = 0; /* {<f1|f1> = SUM{wht[hrad+j]*j*j} */
    
    /* Independent term: */
    double b0 = 0; /* {<f0|val> = SUM{wht[hrad+j]*val[hrad+j]} */
    double b1 = 0; /* {<f1|val> = SUM{wht[hrad+j]*val[hrad+j]*j} */
    
    int j;
    for (j = -hrad; j <= hrad; j++)
      { double vj = val[hrad + j];
        double wj = wht[hrad + j];
        if (wj > 0.0)
          { demand(vj > 0.0, "non-positive price");
            vj = log(vj);
            A00 += wj;
            A01 += wj*j;
            A11 += wj*j*j;
            b0 += wj*vj;
            b1 += wj*vj*j;
          }
      }
    
    if (A00 < 0.000001) { return 0.0; }

    /* Solve the system {A*c = b}: */
    double det = A00*A11 - A01*A01;
    if (fabs(det) < 1.0e-6) { return 0.0; }
    /* Compute adjont of {A}: */
    double M00 = A11; 
    double M01 = -A01; 
    double M11 = A00;
    /* Compute the solution {c0,c1}: */
    double c0 = (M00*b0 + M01*b1)/det;
    double c1 = (M01*b0 + M11*b1)/det;
    /* The fitted polynomial is {P(j) = c0 + c1*j}. */

    if (debug)
      { fprintf(stderr, "(%8.5f) + (%8.5f)*j\n", c0, c1);
        for (j = -hrad; j <= +hrad; j++) 
          { double wj = wht[hrad + j];
            double vj = val[hrad + j];
            if (wj > 0.0)
              { vj = log(vj);
                fprintf(stderr, "  %+4d %8.5f %8.5f %8.5f\n", j, wj, vj, c0+c1*j);
              }
          }
      }
    
    /* Evaluate the polynomial {P} at {j = 0}: */
    return exp(c0);
  }

