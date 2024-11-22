#define PROG_NAME "test_hermite3"
#define PROG_DESC "test of {hermite3.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-20 06:55:05 by stolfi */ 
/* Created on 2012-03-04 by J. Stolfi, UNICAMP */

#define test_hermite3_COPYRIGHT \
  "Copyright © 2014  by the State University of Campinas (UNICAMP)"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <hermite3.h>

#include <bool.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>

int32_t main(int32_t argn, char **argv);

void do_deriv_test(int32_t g);
  /* Tests the derivative estimators with a polynomial of degree {g}. */

void do_plot_test(char *prefix);
  /* Writes a file with a plot of a sample vector interpolated with {hermite3},
    and its derivatives. */
    
void generate_test_samples(int32_t *nxP, double **xP, double *tkerP);
  /* Generates a vector {x[0..n-1]} of samples to be interpolated,
    returns {n} in {*nxP} and {x} in {*xP}.
    Also returns in {*tkerP} the central {t} 
    of the kernel plot. */

int32_t main (int32_t argc, char **argv)
  {
    demand(argc == 2, "wrong num of parameters");
    char *prefix = argv[1];
    
    do_deriv_test(0);
    do_deriv_test(1);
    do_deriv_test(2);
    do_deriv_test(3);
    do_plot_test(prefix);

    return 0;
  }

void do_deriv_test(int32_t g)
  {
    demand((g >= 0) && (g <= 3), "invalid test poly degree");
    
    /* Pick an arbitrary polynomial, set {P[0..g]} to its coefs, {D[0..g-1]} to coefs of its deriv: */
    double P[g+1], D[g];
    int32_t i, k;
    for (k = 0; k <= g; k++) 
      { /* P[k] = sin(k); */ /* Pseudorandom. */
        P[k] = (k == g ? 1.0 : 0.0); /* Monomial {r^g}. */
        if (k > 0) { D[k-1] = k*P[k]; }
      }
    
    /* Select number of samples so that we test all special cases of {hermite3_estimate_derivs}. */ 
    int32_t nv = (g < 3 ? g + 1 : 5);
   
    /* Sample {P,D} at unit steps: */
    double v[nv];
    double dv[nv];
    double magmax = 0.0;  /* Max magnitude of internal values, for error checking. */
    for (i = 0; i < nv; i++) 
      { /* Evaluate polynomial and deriv at argument {i}: */
        double p = 0.0;
        double dp = 0.0;
        double mag = 0.0;  /* Magnitude of internal values. */
        for (k = g; k >= 0; k--) 
          { p = p*i + P[k];
            mag = mag*i + fabs(P[k]);
            if (k < g) 
              { dp = dp*i + D[k];
                mag = mag*i + fabs(D[k]);
              }
          }
        v[i] = p;
        dv[i] = dp;
        if (mag > magmax) { magmax = mag; }
      }
    
    /* Estimate derivatives: */
    double hdv[nv];
    hermite3_estimate_derivs(nv, v, hdv);
    
    /* Compare with truth: */
    int32_t nerr = 0;
    for (i = 0; i < nv; i++) 
      { double ei = hdv[i] - dv[i];
        if (fabs(ei) > 1.0e-12*magmax)
          { fprintf(stderr, "** derivative bug: g = %d  i = %3d", g, i);
            fprintf(stderr, "  expected = %24.16e  estimated = %24.16e", dv[i], hdv[i]);
            fprintf(stderr, "  error = %24.16e\n", ei);
            nerr++;
          }
      }
    demand(nerr == 0, "aborted");
  }

void do_plot_test(char *prefix)
  {
    /* Choose the number of test samples {nx} and the samples {x[0..nx-1]}: */
    int32_t nx;
    double *x;
    double tker;
    generate_test_samples(&nx, &x, &tker);

    /* bool_t debug = FALSE; */ /* If TRUE, {interp} will print the weights. */
    
    /* Open the plot file {wr}: */
    char *fname = jsprintf("out/%s.txt", prefix);
    FILE *wr = open_write(fname, TRUE);

    /* Choose subsampling factor {ns} and total number of subsamples {ny}: */
    int32_t hs = 20;    /* Half of subsamples per data sample. */
    int32_t ns = 2*hs;  /* Subssamples per data sample; must be even. */
    int32_t ny = ns*(nx-1) + 1; /* Total samples in subsampled sequnce. */
    
    /* Subsample: */
    double *y = talloc(ny, double);
    hermite3_subsample(nx, x, NULL, ns, ny, y);
    
    /* Plot {x} interpolated on {ny} subsampling points: */
    int32_t k;
    for (k = 0; k < ny; k++)
      { double tk = ((double)k)/((double)ns);
        /* debug = (fabs(t - tker) <= 0.5*(double)nw); */
        double yk = y[k];
        /* Data samples are located at {k = i*ns} for integer {i}: */
        double xj = ((k % ns) == 0 ? x[k/ns] : NAN);
        fprintf(wr, "%10.6f %+12.6f %+12.6f\n", tk, yk, xj);        
      }
    fclose(wr);
    /* Cleanup: */
    free(fname);
    free(x);
    return;
  }  

void generate_test_samples(int32_t *nxP, double **xP, double *tkerP)
  {
    int32_t H_seg = 10;             /* Half-width of each test segment. */
    int32_t H_ker = 6;              /* Max half-width of kernel. */
    int32_t H_out = H_seg + H_ker;  /* Max half-width of interpolated test segment. */

    int32_t W_ker = 2*H_ker + 1;    /* Max total width of kernel. */
    int32_t W_out = 2*H_out + 1;    /* Max total width of interpolated test segment. */

    int32_t DX =  3;       /* Space between interpolated test segments. */
    
    /* The test data consists f a single-sample impulse followed by several
      broad polynomial pulses. Each broad pulse spans {W_seg} data samples
      and its values are defined by a monomial of degree {g}, for {g} from 0
      to {deg+max}; that is, {x[i] = A*((i-c[g])/H_seg)^g} where {x[c[g]]} is the
      central sample of the pulse.  The broad pulses are spaced so that 
      even after interpolation they will be separated by {DX} data steps. */
    
    
    /* Compute the number of data samples {nx}. */
    int32_t deg_max = 4;
    int32_t nx = DX + W_ker + (deg_max+1)*(DX + W_out) + DX;
    double *x = talloc(nx, double);
    
    auto void test_segm(int32_t *kP, int32_t g);
      /* Appends another broad test segment {x[kini..kfin]} with degree {g} 
        to the data sample vector, where {kini} is the input value of {*kP}.
        Also updates {*kP} with {kfin+1}. */
    
    /* Clear all samples: */
    int32_t i;
    for (i = 0; i < nx; i++) { x[i]= 0; }
    
    int32_t ks = 0; /* Next sample to be defined is {x[ks]}. */

    /* Skip some samples: */
    ks += DX;
    
    /* Lay down the single-sample impulse: */
    ks += H_ker;
    x[ks] = 1.0; ks++;
    double tker = (double)ks - 0.5;
    ks += H_ker;
    
    /* Lay down the broad polynomial pulses: */
    int32_t gg;
    for (gg = 0; gg <= deg_max; gg++)
      { ks += DX;
        ks += H_ker;
        test_segm(&ks, gg);
        ks += H_ker;
      }
      
    /* Space after the last pulse: */
    ks += DX;
    
    /* We must be done, return: */
    assert(ks == nx);
    (*nxP) = nx;
    (*xP) = x;
    (*tkerP) = tker;
    return;
    
    /* INTERNAL IMPLEMENTATIONS: */
    
    void test_segm(int32_t *kP, int32_t g)
      { int32_t j;
        for (j = -H_seg; j <= +H_seg; j++)
          { double t = ((double)j)/((double)H_seg);
            x[(*kP)] = pow(t, g);
            (*kP)++;
          }
      }
  }
