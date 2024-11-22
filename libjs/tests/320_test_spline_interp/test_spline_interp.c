#define PROG_NAME "test_interp_spline"
#define PROG_DESC "test of {interp_spline.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-20 06:55:48 by stolfi */ 
/* Created on 2012-03-04 by J. Stolfi, UNICAMP */

#define test_interp_spline_COPYRIGHT \
  "Copyright © 2012  by the State University of Campinas (UNICAMP)"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>


#include <bool.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <ix.h>

#include <interp_spline.h>

int32_t main(int32_t argn, char **argv);

void do_plot_tests(char *prefix);

void do_plot_test(char *prefix, int32_t order, interp_spline_kind_t kind);
  /* Writes a file with a plot of a sample vector interpolated with the specified {order}
    and {kind}, and its derivatives. */
    
void generate_test_samples(uint32_t *nsP, double **sP, double *zkerP);
  /* Generates a vector {s[0..n-1]} of samples to be interpolated,
    returns {n} in {*nsP} and {s} in {*sP}.
    Also returns in {*zkerP} the central {z} 
    of the kernel plot. */

int32_t main (int32_t argc, char **argv)
  {
    demand(argc == 2, "wrong num of parameters");
    char *prefix = argv[1];
    
    do_plot_tests(prefix);

    return 0;
  }
 
#define order_MAX 4
#define deg_MAX (order_MAX + 1)

void do_plot_tests(char *prefix)
  {
    int32_t order, kind;
    for (order = -1; order <= order_MAX; order++)
      { for (kind = 0; kind < interp_spline_kind_NUM; kind++)
          { do_plot_test(prefix, order, kind); }
      }
  }

void do_plot_test(char *prefix, int32_t order, interp_spline_kind_t kind)
  {
    /* Get the number of kernel weigts {nw} : */
    uint32_t nw = interp_spline_compute_num_samples(order, kind);
    if (nw == 0) 
      { /* Invalid {order}+{kind} combination, skip it: */
        return;
      }

   /* Choose the number of test samples {ns} and the samples {s[0..ns-1]}: */
    uint32_t ns;
    double *s;
    double zker;
    generate_test_samples(&ns, &s, &zker);

    /* Allocate storage for the interpolation kernel weights {wt[0..nw-1]} and tap positions {ix[0..nw-1]} : */
    double wt[nw];    /* Interpolation kernel weights. */
    int32_t ix[nw];   /* Tap positions in sample array. */
    
    /* Choose the boundary conditions: */
    ix_reduction_t red = ix_reduction_SINGLE; /* Samples outside {0..ns-1} are "not there". */
    
    bool_t debug = FALSE; /* If TRUE, {interp} will print the weights. */

    auto double interp(double z);
      /* Interpolates the data {s[0..ns-1]} at the fracional index {z}. */
    
    /* Open the plot file {wr}: */
    char redTag = "SERMP"[red];
    char orderTag = (char)(order < 0 ? 'm' : '0'+order);
    char kindTag = "BIO"[kind];
    char *fname = jsprintf("out/%s-r%c-o%c-k%c.txt", prefix, redTag, orderTag, kindTag);
    FILE *wr = open_write(fname, TRUE);

    /* Choose subsampling factor {NU} and total number of subsamples {NP}: */
    uint32_t HU = 64;    /* Half of subsamples per data sample. */
    uint32_t NU = 2*HU;  /* Subssamples per data sample; must be even. */
    uint32_t NP = NU*ns; /* Total subsamples. */
    
    /* Plot {s} interpolated on {NP} subsampling points: */
    for (int32_t k = 0; k <= NP; k++)
      { double z = ((double)k)/((double)NU);
        /* debug = (fabs(z - zker) <= 0.5*(double)nw); */
        double Fz = interp(z);
        /* Data samples are located at {k = HU + i*NU} for integer {i}: */
        double Sz = ((k % NU) == HU ? s[(k-HU)/NU] : NAN);
        fprintf(wr, "%10.6f %+12.6f %+12.6f\n", z, Fz, Sz);        
      }
    fclose(wr);
    /* Cleanup: */
    free(fname);
    free(s);
    return;
    
    /* Internal implementations: */
    
    double interp(double z)
      { /* Get the sample indices {ix[0..nw-1]} that are needed to interpolate at {z}: */
        interp_spline_get_indices(z, ns, red, nw, ix);
        /* Get the weights {wt[0..nw-1]} of those samples in the interpolating kernel: */
        interp_spline_get_weights(z, order, kind, nw, wt);
        /* Interpolate: */
        double sum_ws, sum_wt = 0;
        for (int32_t j = 0; j < nw; j++) 
          { if (ix[j] >= 0) { sum_ws += wt[j]*s[ix[j]]; sum_wt += wt[j]; } }
        double f = sum_ws/sum_wt;
        if (debug)
          { fprintf(stderr, "z = %10.7f  f(z) = %+10.7f\n", z, f);
            for (int32_t k = 0; k < nw; k++) { fprintf(stderr, "  wt[%d] = %10.7f\n", k, wt[k]); }
            fprintf(stderr, "\n");
          }
        return f;
      }
  }  

void generate_test_samples(uint32_t *nsP, double **sP, double *zkerP)
  {
    uint32_t H_seg = 10;             /* Half-width of each test segment. */
    uint32_t H_ker = 6;              /* Max half-width of kernel. */
    uint32_t H_out = H_seg + H_ker;  /* Max half-width of interpolated test segment. */

    uint32_t W_ker = 2*H_ker + 1;    /* Max total width of kernel. */
    uint32_t W_out = 2*H_out + 1;    /* Max total width of interpolated test segment. */

    uint32_t DX =  3;       /* Space between interpolated test segments. */
    
    /* The test data consists f a single-sample impulse followed by several
      broad polynomial pulses. Each broad pulse spans {W_seg} data samples
      and its valuesare defined by a monomial of degree {g}, for {g} from 0
      to {deg+max}; that is, {s[i] = A*((i-c[g])/H_seg)^g} where {s[c[g]]} is the
      central sample of the pulse.  The broad pulses are spaced so that 
      even after interpolation they will be separated by {DX} data steps. */
    
    
    /* Compute the number of data samples {ns}. */
    uint32_t ns = DX + W_ker + (deg_MAX+1)*(DX + W_out) + DX;
    double *s = notnull(malloc(ns*sizeof(double)), "no mem");
    
    auto void test_segm(uint32_t *kP, uint32_t g);
      /* Appends another broad test segment {s[kini..kfin]} with degree {g} 
        to the data sample vector, where {kini} is the input value of {*kP}.
        Also updates {*kP} with {kfin+1}. */
    
    /* Clear all samples: */
    for (int32_t i = 0; i < ns; i++) { s[i]= 0; }
    
    uint32_t ks = 0; /* Next sample to be defined is {s[ks]}. */

    /* Skip some samples: */
    ks += DX;
    
    /* Lay down the single-sample impulse: */
    ks += H_ker;
    s[ks] = 1.0; ks++;
    double zker = (double)ks - 0.5;
    ks += H_ker;
    
    /* Lay down the broad polynomial pulses: */
    for (int32_t gg = 0; gg <= deg_MAX; gg++)
      { ks += DX;
        ks += H_ker;
        test_segm(&ks, gg);
        ks += H_ker;
      }
      
    /* Space after the last pulse: */
    ks += DX;
    
    /* We must be done, return: */
    assert(ks == ns);
    (*nsP) = ns;
    (*sP) = s;
    (*zkerP) = zker;
    return;
    
    /* INTERNAL IMPLEMENTATIONS: */
    
    void test_segm(uint32_t *kP, uint32_t g)
      { for (int32_t j = -(int32_t)H_seg; j <= +(int32_t)H_seg; j++)
          { double x = ((double)j)/((double)H_seg);
            s[(*kP)] = pow(x, g);
            (*kP)++;
          }
      }
  }
