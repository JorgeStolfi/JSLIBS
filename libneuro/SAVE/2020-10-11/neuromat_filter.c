/* See {neuromat_filter.h}. */
/* Last edited on 2013-11-21 02:53:19 by stolfilocal */

/* !!! Replace apodizing with trend-removal !!! */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#include <fftw3.h>

#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <neuromat_eeg.h>
#include <neuromat_poly.h>
#include <neuromat_filter.h>
 
double *neuromat_filter_tabulate_gain(int nf, double fsmp, neuromat_filter_t *gain, bool_t verbose);
  /* Builds a table of Hartley filter coefficients {G[0..nf-1]} from a
    complex Fourier filter transfer function {gain}, suitable for
    scaling a Hartley transform with {nf} coefficients.
    
    Namely, let {F[kf]} be the complex number {gain(kf*fsmp/nf,fsmp/2)}
    for {kf} in {0..nf-1}. The procedure assumes (requires) that {F[kf]}
    and {F[nf-kf]} are complex conjugates; in particular, {F[0]} is pure
    real, and so is {F[nf/2]} if {nf} is even. The procedure stores in
    {G[kf],G[nf-kf]} the coefficient pair that in the Hartley basis is
    equivalent to the coefficient pair {F[kf],F[nf-kf]} in the (complex)
    Fourier basis.
    
    !!! Should fold-over the response. !!! */
 
void neuromat_filter_apply
  ( int nt, 
    int ne, 
    double **val, 
    double fsmp, 
    int tdeg, 
    bool_t tkeep, 
    neuromat_filter_t *gain,
    bool_t verbose
  )
  {
    demand(fsmp > 0.0, "invalid sampling frequency {fsmp}");
    
    /* Allocate the FFT work areas and plans: */
    double *in = (double*) fftw_malloc(sizeof(double)*(nt+4));
    double *out = (double*) fftw_malloc(sizeof(double)*(nt+4));
    fftw_plan pd = fftw_plan_r2r_1d(nt, in, out, FFTW_DHT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    fftw_plan pi = fftw_plan_r2r_1d(nt, out, in, FFTW_DHT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

    /* Fourier transfer function: */
    double *G = neuromat_filter_tabulate_gain(nt, fsmp, gain, verbose);
    
    /* Fitted polynomial coeffs and values: */
    double *P = (tdeg < 0 ? NULL : notnull(malloc(sizeof(double)*(tdeg+1)), "no mem"));
    double *s = (tdeg < 0 ? NULL : notnull(malloc(sizeof(double)*nt), "no mem"));
       
    int it, ie;
    for (ie = 0; ie < ne; ie++) 
      { /* Copy the signal {ie} into the FFT buffer, with mirrored apodized ends and zero padding: */
        for (it = 0; it < nt; it++) { in[it] = val[it][ie]; }
        
        if (tdeg >= 0)
          { /* Fit the polynomial {P}: */
            int maxiter = 5;
            neuromat_poly_fit_robust(nt, NULL, in, NULL, maxiter, tdeg, P, NULL);
            if (verbose) 
              { fprintf(stderr, "  trend = ");
                int r; 
                for (r = 0; r <= tdeg; r++) { fprintf(stderr, " %+12.6f", P[r]); }
                fprintf(stderr, "\n");
              }
            /* Evaluate it: */
            neuromat_poly_eval_multi(tdeg, P, nt, NULL, s);
            /* Subtract it from the data: */
            for (it = 0; it < nt; it++) { in[it] -= s[it]; }
          }

        /* Transform to frequency domain: */
        fftw_execute(pd);

        /* Apply frequency filter: */
        int kf0;
        for (kf0 = 0; kf0 <= nt-kf0; kf0++) 
          { int kf1 = (nt - kf0) % nt; /* The other Hartley element with same absolute freq. */
            double h0 = out[kf0];
            double G0 = G[kf0];
            if (kf1 == kf0)
              { /* Coef {kf0} is the only one with that frequency: */
                out[kf0] = h0*G0;
              }
            else
              { /* Coefs {kf0,kf1} have the same frequency and get mixed: */
                double h1 = out[kf1]; 
                double G1 = G[kf1];
                /* Apply filter: */
                double Gp = G0 + G1;
                double Gm = G0 - G1;
                out[kf0] = 0.5*(h0*Gp + h1*Gm);
                out[kf1] = 0.5*(h1*Gp - h0*Gm);
              }
          }

        /* Return to time domain: */
        fftw_execute(pi);
        
        /* Return filtered signal & trend to {val}: */
        for (it = 0; it < nt; it++) 
          { /* Scale to preserve total power: */
            in[it] /= nt;
            /* Restore trend if any: */
            if ((tdeg >= 0) && tkeep) { in[it] += s[it]; }
            /* Save in {val}: */
            val[it][ie] = in[it];
          }
      }
      
    fftw_destroy_plan(pd);
    fftw_destroy_plan(pi);
    fftw_free(out);
    fftw_free(in);
    free(G);
    if (P != NULL) { free(P); }
    if (s != NULL) { free(s); }
    return;
  }
    
double *neuromat_filter_tabulate_gain(int nf, double fsmp, neuromat_filter_t *gain, bool_t verbose)
  { 
    double *G = notnull(malloc(nf*sizeof(double)), "no mem");
    
    int kf0;
    for (kf0 = 0; kf0 <= nf-kf0; kf0++)
      { int kf1 = (nf - kf0) % nf;  /* Index of the Hartley coeff with the same frequency. */
        /* Compute filter gain at frequency {f}: */
        complex F0 = gain(kf0, nf, fsmp);
        if (verbose) { fprintf(stderr, "%6d %10.6f  %+14.10f %+14.10f", kf0, kf0*fsmp/nf, creal(F0), cimag(F0)); }
        /* Convert Fourier coeffs {Fp,Fm} to Hartley coeffs {G[kf0],G[kf1]}: */
        if (kf0 == kf1)
          { /* Fourier coeff same as Hartley coeff: */
            demand(fabs(cimag(F0)) < 1.0e-13, "Fourier coeff should be pure real");
            G[kf0] = creal(F0);
          }
        else
          { complex F1 = gain(kf1,nf, fsmp);
            if (verbose) { fprintf(stderr, "  %6d %10.6f  %+14.10f %+14.10f", kf1, kf1*fsmp/nf, creal(F1), cimag(F1)); }
            demand(fabs(cimag(F0)+cimag(F1)) < 1.0e-13, "Fourier coeffs should be conjugate");
            /* Force conjugation: */
            double FR = 0.5*(creal(F0) + creal(F1));
            double FI = 0.5*(cimag(F0) - cimag(F1));
            /* Convert Fourier coeff pair to Hartley coeff pair: */
            G[kf1] = FR + FI;
            G[kf0] = FR - FI;
          }
        if (verbose) { fprintf(stderr, "\n"); }
      }
    return G;
  }
  
double neuromat_filter_lowpass_sigmoid(double f, double fa, double fb)
  {
    demand(fa > 0, "invalid {fa}");
    demand(fb > fa, "invalid {fb}");
    double ga = log(fa), gb = log(fb), g = log(f);
    double gm = 0.5*(ga + gb);
    double gr = 0.5*fmax(gb - ga, 1.0e-200);
    double C = 2.19; /* Approx {erf^-1(0.998)}. */
    double z = C*(g - gm)/gr;
    if (z >= 6.0) 
      { return 0; }
    else if (z <= -6.0)
      { return 1; }
    else
      { return 0.5*(1 - erf(z)); }
  }

double neuromat_filter_lowpass_butterworth(double f, double fc, int n)
  {
    if (fc == 0) { return 0; } 
    double s = f/fc;
    if (s >= 1.0e32)
      { return 0; }
    else 
      { double t = log(s);
        if (t <= (-16*M_LN10)/n)
          { return 1; }
        else if (s >= 1.0e16)
          { return exp((0.5*n)*t); }
        else
          { return sqrt(1/(1 + exp(n*t))); }
      }
  }
  
complex neuromat_filter_lowpass_cbutterworth(double f, double fc, int n)
  {
    if (fc == 0) { return 0; } 
    int odd = (n % 2);
    complex s = f/fc*I;
    complex B = (odd ? 1 + s : 1);
    int n2 = n/2;
    double u = M_PI/n/2;
    int k;
    for (k = 1; k <= n2; k++)
      { complex t = s*s - 2*s*cos((2*k+n-1)*u) + 1;
        B *= t;
      }
    return 1/B;
  }

double neuromat_filter_lowpass_gauss(double f, double fc, double gc, double fsmp)
  { 
    demand(fc >= 0, "invalid {fc}");
    demand(fsmp > 0, "invalid {fsmp}");
    /* Compute the Gaussian reference frequency {fd}: */
    demand((gc > 0) && (gc < 1), "invalid {gc}");
    double fd = fc/sqrt(-log(gc));
    demand(fd < 10*fsmp, "invalid {fc,gc}");
    /* Reduce frequency {f} to {[-fmax _ +fmax]}: */
    f = f - fsmp*floor(f/fsmp + 0.5);
    /* Dominating Gaussian term {gf}: */
    double zf = f/fd;
    double gf = exp(-zf*zf);
    /* Range of summation for folded-over tails: */
    int M = (int)ceil(8.5*fd/fsmp);
    int k;
    double gtf = 0, gt0 = 0;
    for (k = M; k > 0; k--)
      { double zkfp = (f + k*fsmp)/fd;
        double gkfp = exp(-zkfp*zkfp);
        double zkfm = (f - k*fsmp)/fd;
        double gkfm = exp(-zkfm*zkfm);
        double zk0 = k*fsmp/fd;
        double gk0 = exp(-zk0*zk0);
        gtf += (gkfp+gkfm);
        gt0 += 2*gk0;
      }
    return (gf + gtf)/(1 + gt0);
  } 

double neuromat_filter_lowpass_biquadratic(double f, double fmax)
  { demand(fmax >= 0, "invalid {fmax}");
    if (fabs(f) >= fmax) { return 0; }
    if (f == 0) { return 1; }
    double z = f/fmax;
    double y = 1 - z*z;
    double g = y*y;
    return g;
  } 

double **neuromat_filter_compute_spectra(int nt, int ne, double **val, int kfmax, bool_t verbose)  
  {
    demand(kfmax <= nt/2, "parameter {kfmax} exceeds the Nyquist frequency");
    
    /* Allocate the output: */
    double **pwr = notnull(malloc(ne*sizeof(double*)),"no mem"); 
    
    /* Allocate the work area: */
    double *in = (double*) fftw_malloc(sizeof(double) * nt);
    double *out = (double*) fftw_malloc(sizeof(double) * nt);
    
    /* Apodizing window: */
    double *W = notnull(malloc(nt*sizeof(double)), "no mem");
    int it;
    for (it = 0; it < nt; it++) { W[it] = 0.5*(1 - cos(2*M_PI*(it+0.5)/nt)); }
    
    /* Do the row transforms: */
    fftw_plan px = fftw_plan_r2r_1d(nt, in, out, FFTW_DHT, FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
    int ie;
    for (ie = 0; ie < ne; ie++)
      { /* Allocate and compute the spectrum of electrode {ie}: */
        pwr[ie] = notnull(malloc((kfmax+1)*sizeof(double)), "no mem"); 
        double sum2_in = 0;
        for (it = 0; it < nt; it++) 
          { double vt = val[it][ie] * W[it];
            sum2_in += vt*vt;
            in[it] =  vt;
          }
        double rms_in = sqrt(sum2_in / nt);
        fftw_execute(px);
        if (verbose) { fprintf(stderr, "  electrode %d\n", ie); }
        int kf0;
        double sum2_out = 0;
        for (kf0 = 0; kf0 <= kfmax; kf0++) 
          { /* Get the index {kf1} of the coefficient with same actual frequency {kf0}: */
            int kf1 = (nt - kf0) % nt;
            /* Compute the total power in that frequency: */
            double c0 = out[kf0];
            double c1 = out[kf1];
            if (verbose) 
              { fprintf(stderr, "    kf0 = %5d out[kf0] = %24.15e", kf0, c0); 
                if (kf1 != kf0) { fprintf(stderr, "  kf1 = %5d out[kf1] = %24.15e", kf1, c1); }
                fprintf(stderr, "\n");
              }
            double pf = (c0*c0 + (kf0 != kf1 ? c1*c1 : 0.0))/nt;
            sum2_out += pf;
            pwr[ie][kf0] = pf/nt;
          }
        double rms_out = sqrt(sum2_out/nt);
        if (verbose) { fprintf(stderr, "rms ratio out/in = %24.15e\n", rms_out/rms_in); }
      }
    fftw_destroy_plan(px);
    
    fftw_free(out);
    fftw_free(in);
    free(W);
    return pwr;
  }

void neuromat_filter_write_spectra
  ( FILE *wr, 
    int nt, 
    int ne, 
    int kfmax, 
    double fsmp, 
    double **pwr
  )
  {
    demand(2*kfmax <= nt, "{kfmax} too high");
    int kf, ie;
    
    for (kf = 0; kf <= kfmax; kf++) 
      { double f = kf*fsmp/nt;
        double flo = (kf == 0 ? 0.0 : kf - 0.5)*fsmp/nt;
        double fhi = (2*kf == nt ? (double)kf : kf + 0.5)*fsmp/nt;
        fprintf(wr, "%8d  %12.7f  %12.7f %12.7f ", kf, f, flo, fhi);
        for (ie = 0; ie < ne; ie++) { fprintf(wr, " %12.5e", pwr[ie][kf]); }
        fprintf(wr, "\n");
      }
    fflush(wr);
  }

