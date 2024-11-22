/* Tests the Butterworth filter formulas. */
/* Last edited on 2024-01-06 08:25:37 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsfile.h>
#include <jsmath.h>
#include <jsrandom.h>

#include <neuromat_filter.h>
#include <neuromat_filter_tabulate_hartley_gains.h>
#include <neuromat_filter_clear_tiny_gains.h>

#include <neuromat_filter_lowpass_sgerf.h>
#include <neuromat_filter_lowpass_gauss.h>
#include <neuromat_filter_lowpass_biquadratic.h>
#include <neuromat_filter_lowpass_butterworth.h>

#include <neuromat_filter_bandpass_log_gauss.h>
#include <neuromat_filter_bandpass_sgerf.h>
#include <neuromat_filter_bandpass_butterworth.h>
#include <neuromat_filter_bandpass_biquadratic.h>
#include <neuromat_filter_bandpass_gauss.h>

typedef double tfi_function_t(double f, int32_t np);
  /* A function that returns a gain for frequency {f}.  The parameter {np} is the order
    (number of poles), which is relevant for some functions like 
     the Butterworth filters. */

void tfi_test_hartley_basis(void);
  /* Tests {neuromat_filte_hartley_basis_eval} and whether the
    Hartley basis is orthonormal under the ordinary dot product. */
    
void tfi_test_fourier_basis(void);
  /* Tests {neuromat_filte_fourier_basis_eval} and whether the
    Fourier basis is orthonormal under the complex dot product. */ 

void tfi_test_hartley_fourier_conversion(void);
  /* Tests the Hartley <--> Fourier conversion functions
    {neuromat_filter_hartley_to_fourier} and {neuromat_filter_fourier_to_hartley}. */

void tfi_test_clear_tiny_gains(int32_t nf, double eps);
  /* Tests {neuromat_filter_clear_tiny_gains} with {nf} frequencies and
    threshold {eps}. Writes a file "out/clear_tiny.txt". Each line has a
    freq index {kf} then three sets of three numbers for the original
    Fourier coeff {Fraw}, the cleaned Fourier coeff {Fkuk}, and the
    difference {D} between them. Each group has the real part, imaginary
    part, and modulus of the complex number. */

void tfi_test_filters(double fsmp, double flo0, double flo1, double fhi1, double fhi0);
  /* Plots various filters with the given parameters. */

void tfi_test_tabulate_hartley_gains(int32_t nf, char *ftype, double fsmp);
  /* Tests {neuromat_filter_tabulate_hartley_gains} with {nf} frequencies,
    sampling frequency {fsmp}, and filter type {ftype}.  The latter 
    is "UG" for unsymmetric Gaussian, "DL" for delay, "BQ", "ER", "BW" for
    lowpass biquadratic sigmoid, erf-log-sigmoid, and Butterworth.
    Those filters except "DL" will be subjected to folding and small-value cleanup. */

void tfi_plot_function(FILE *wr, int32_t npmax, tfi_function_t *gain); 
  /* Writes to {wr} the plot of the function {gain(f,np)} varying {f} gradually from 0 to 0.5. 
    
     Each plot has actually {npmax+1} graphs, with {np} varying
     from 1 to {npmax}.  if the filters do not have an order parameter, {npmax}
     should be 1. */

void tfi_plot_function_pair(char *tag, int32_t npmax, tfi_function_t *lopa, tfi_function_t *band);
  /* Given two versions of a filter function, {lopa} (low-pass) and {band} (bandpass),
    Plots the gain function {lopa(f,np)} to "out/{tag}_lopa.txt" and the function {band(f,np)} to 
    "out/{tag}_band.txt", using {tfi_plot_function}. */
     
int32_t main(int32_t argc, char **argv);
  /* Command line arguments are "{fsmp} {flo0} {flo1} {fhi0} {fhi1}". */

int32_t main(int32_t argc, char **argv)
  {
    /* Assumed sampling frequency: */
    demand(argc == 6, "wrong number of arguments");
    char *rest;
    double fsmp = strtod(argv[1], &rest); demand((*rest) == 0, "invalid argv[1]");
    double flo0 = strtod(argv[2], &rest); demand((*rest) == 0, "invalid argv[2]");
    double flo1 = strtod(argv[3], &rest); demand((*rest) == 0, "invalid argv[3]");
    double fhi1 = strtod(argv[4], &rest); demand((*rest) == 0, "invalid argv[4]");
    double fhi0 = strtod(argv[5], &rest); demand((*rest) == 0, "invalid argv[5]");
    
    demand(flo0 <= flo1, "invalid {flo0,flo1}");
    demand(flo1 <= fhi1, "invalid {flo1,fhi1}");
    demand(fhi1 <= fhi0, "invalid {fhi1,fhi0}");
    demand(fhi1 <= fsmp/2, "limit {fhi1} too large");
   
    /* Tests {neuromat_filter.h}: */
    tfi_test_hartley_basis();
    tfi_test_fourier_basis();
    tfi_test_hartley_fourier_conversion();
    tfi_test_clear_tiny_gains(300, 1.0e-11);
    tfi_test_filters(fsmp, flo0, flo1, fhi1, fhi0);
    
    for (int32_t bnf = 9; bnf <= 10; bnf++)
      { for (int32_t mnf = 0; mnf < 2; mnf++)
        { int32_t nf = (mnf == 0 ? bnf : 50*bnf + (bnf % 2));
          tfi_test_tabulate_hartley_gains(nf, "UG", fsmp);
          tfi_test_tabulate_hartley_gains(nf, "ER", fsmp);
          tfi_test_tabulate_hartley_gains(nf, "BQ", fsmp);
          tfi_test_tabulate_hartley_gains(nf, "BW", fsmp);
          tfi_test_tabulate_hartley_gains(nf, "DL", fsmp);
        }
      }
    
    return 0;
  }
  
void tfi_test_hartley_basis(void)
  { double tol = 1.0e-14;
    for (int32_t nf = 5; nf <= 6; nf++)
      { fprintf(stderr, "--- %s nf = %d ---\n", __FUNCTION__, nf);
        double *etai = talloc(nf, double);
        double *etaj = talloc(nf, double);
        for (int32_t jf = 0; jf < nf; jf++)
          { for (int32_t t = 0; t < nf; t++)
              { etai[t] = neuromat_filter_hartley_basis_eval(nf, jf, t); }
            for (int32_t kf = jf; kf < nf; kf++)
              { for (int32_t t = 0; t < nf; t++)
                  { etaj[t] = neuromat_filter_hartley_basis_eval(nf, kf, t); }
                double dot = 0.0;
                for (int32_t t = 0; t < nf; t++)
                  { dot += etai[t]*etaj[t]; }
                double dot_exp = (jf == kf ? 1.0 : 0.0);
                if (fabs(dot - dot_exp) > tol)
                  { fprintf(stderr, "dot(eta[%d],eta[%d]) = %24.16e\n", jf, kf, dot);
                    demand(FALSE, "Hartley basis not orthonormal");
                  }
              }
          }
        free(etai);
        free(etaj);
      }
  }
  
void tfi_test_fourier_basis(void)
  { double tol = 1.0e-14;
    for (int32_t nf = 5; nf <= 6; nf++)
      { fprintf(stderr, "--- %s nf = %d ---\n", __FUNCTION__, nf);
        complex *phii = talloc(nf, complex);
        complex *phij = talloc(nf, complex);
        for (int32_t jf = 0; jf < nf; jf++)
          { for (int32_t t = 0; t < nf; t++)
              { phii[t] = neuromat_filter_fourier_basis_eval(nf, jf, t); }
            for (int32_t kf = jf; kf < nf; kf++)
              { for (int32_t t = 0; t < nf; t++)
                  { phij[t] = neuromat_filter_fourier_basis_eval(nf, kf, t); }
                complex dot = 0.0;
                for (int32_t t = 0; t < nf; t++)
                  { dot += phii[t]*conj(phij[t]); }
                double dot_exp = (jf == kf ? 1.0 : 0.0);
                if (cabs(dot - dot_exp) > tol)
                  { fprintf(stderr, "dot(phi[%d],phi[%d] = %24.16e + %24.16e * I\n", jf, kf, creal(dot), cimag(dot));
                    demand(FALSE, "Fourier basis not orthonormal");
                  }
              }
          }
        free(phii);
        free(phij); 
      
      }
  }
     
void tfi_test_hartley_fourier_conversion(void)
  { double tol = 1.0e-14;
    for (int32_t nf = 5; nf <= 6; nf++)
      { fprintf(stderr, "--- %s nf = %d ---\n", __FUNCTION__, nf);
        double *x = talloc(nf, double);      /* A random sample vector. */
        double *eta = talloc(nf, double);    /* A Hartley basis element. */
        double *H = talloc(nf, double);      /* The Hartley transform of {x}. */
        complex *phi = talloc(nf, complex);  /* A Fourier basis element. */     
        complex *F = talloc(nf, complex);    /* The Fourier transform of {x}. */
        for (int32_t k = 0; k < 10; k++)
          { /* Generate a random real signal vector: */
            for (int32_t t = 0; t < nf; t++) { x[t] = 2*drandom() - 1; }
            /* Compute the Hartley and Fourier transforms of {x} by brute force: */
            for (int32_t f = 0; f < nf; f++)
              { for (int32_t t = 0; t < nf; t++)
                  { eta[t] = neuromat_filter_hartley_basis_eval(nf, f, t);
                    phi[t] = neuromat_filter_fourier_basis_eval(nf, f, t);
                  }
                double dotH = 0;
                complex dotF = 0;
                for (int32_t t = 0; t < nf; t++)
                  { dotH += x[t]*eta[t];
                    dotF += x[t]*conj(phi[t]);
                  }
                H[f] = dotH;
                F[f] = dotF;
              }
            /* Check the conversion between coefficients: */
            int32_t nerr_HtoF = 0, nerr_FtoH = 0;
            for (int32_t fa = 0; fa < nf; fa++)
              { int32_t fb = (nf - fa) % nf;
                /* Hartley to Fourier: */
                complex Fa, Fb;
                neuromat_filter_hartley_to_fourier(H[fa], H[fb], &Fa, &Fb);
                double err_Fa = cabs(Fa - F[fa]);
                double err_Fb = cabs(Fb - F[fb]);
                if ((err_Fa > tol) || (err_Fb > tol))
                  { fprintf(stderr, "  H[%d],H[%d] --> F[%d],F[%d]:\n", fa, fb, fa, fb);
                    fprintf(stderr, "    input %24.16e , %24.16e\n", H[fa], H[fb]);
                    fprintf(stderr, "    result:\n");
                    fprintf(stderr, "    ( %24.16e + %24.16e * I )", creal(Fa), cimag(Fa));
                    fprintf(stderr, " != ( %24.16e + %24.16e * I )", creal(F[fa]), cimag(F[fa]));
                    fprintf(stderr, " error %12.4e\n", err_Fa);
                    fprintf(stderr, "    ( %24.16e + %24.16e * I )", creal(Fb), cimag(Fb));
                    fprintf(stderr, " != ( %24.16e + %24.16e * I )", creal(F[fb]), cimag(F[fb]));
                    fprintf(stderr, " error %12.4e\n\n", err_Fb);
                    nerr_HtoF++;
                  }
                /* Fourier to Hartley: */
                double Ha, Hb;
                neuromat_filter_fourier_to_hartley(F[fa], F[fb], &Ha, &Hb);
                double err_Ha = cabs(Ha - H[fa]);
                double err_Hb = cabs(Hb - H[fb]);
                if ((err_Ha > tol) || (err_Hb > tol))
                  { fprintf(stderr, "  F[%d],F[%d] --> H[%d],H[%d]:\n", fa, fb, fa, fb);
                    fprintf(stderr, "    input");
                    fprintf(stderr, " ( %24.16e + %24.16e * I )", creal(F[fa]), cimag(F[fa]));
                    fprintf(stderr, " , ( %24.16e + %24.16e * I )\n", creal(F[fb]), cimag(F[fb]));
                    fprintf(stderr, "    result:\n");
                    fprintf(stderr, "    %24.16e", Ha);
                    fprintf(stderr, " != %24.16e", H[fa]);
                    fprintf(stderr, " error %12.4e\n", err_Ha);
                    fprintf(stderr, "    %24.16e", Hb);
                    fprintf(stderr, " != %24.16e", H[fb]);
                    fprintf(stderr, " error %12.4e\n\n", err_Hb);
                    nerr_FtoH++;
                  }
              }
            demand(nerr_HtoF == 0, "Hartley to Fourier failed");
            demand(nerr_FtoH == 0, "Fourier to Hartley failed");
          }
          
        free(x);
        free(eta);
        free(H);
        free(phi);
        free(F);
      }
  }
 
void tfi_test_filters(double fsmp, double flo0, double flo1, double fhi1, double fhi0)
  { 
    fprintf(stderr, "--- %s fsmp = %12.8f parms = %12.8f %12.8f  %12.8f %12.8f ---\n", __FUNCTION__, fsmp, flo0, flo1, fhi1, fhi0);
    demand((0 < flo0) && (flo0 < flo1) && (flo1 < fhi1) && (fhi1 < +INF), "invalid key freqs");
    
    /* Filter gain functions (real): */

    auto double band_lgaus(double f, int32_t np);

    auto double lopa_sgerf(double f, int32_t np);
    auto double band_sgerf(double f, int32_t np);

    auto double lopa_gauss(double f, int32_t np);
    auto double band_gauss(double f, int32_t np);

    auto double lopa_biqua(double f, int32_t np);
    auto double band_biqua(double f, int32_t np);

    auto double lopa_buttw(double f, int32_t np);
    auto double band_buttw(double f, int32_t np);
    
    /* Parameters for the filter gain functions: */
    
    double fm_lgaus, sigma_lgaus, fsup_lgaus, mag_lgaus; int32_t np_lgaus;
    neuromat_filter_bandpass_log_gauss_compute_parms
      ( flo0, flo1, fhi1, fhi0, 
        &fm_lgaus, &np_lgaus, &sigma_lgaus, &mag_lgaus,
        &fsup_lgaus,
        TRUE 
      );
    
    double fm_lo_sgerf, sigma_lo_sgerf, fm_hi_sgerf, sigma_hi_sgerf, fsup_sgerf;
    neuromat_filter_bandpass_sgerf_compute_parms
      ( flo0, flo1, fhi1, fhi0, 
        &fm_lo_sgerf, &sigma_lo_sgerf, &fm_hi_sgerf, &sigma_hi_sgerf,
        &fsup_sgerf,
        TRUE 
      );
   
    double fs_lo_buttw, fs_hi_buttw, fsup_buttw;
    int32_t ord_buttw;
    neuromat_filter_bandpass_butterworth_compute_parms
      ( flo0, flo1, fhi1, fhi0, 
        &fs_lo_buttw, &fs_hi_buttw, &ord_buttw, 
        &fsup_buttw,
        TRUE 
      );
    
    double sigma_lo_gauss, sigma_hi_gauss, mag_gauss, fsup_gauss;
    neuromat_filter_bandpass_gauss_compute_parms
      ( flo0, flo1, fhi1, fhi0, 
        &sigma_lo_gauss, &sigma_hi_gauss, &mag_gauss, 
        &fsup_gauss,
        TRUE 
      );

    tfi_plot_function_pair("LG", np_lgaus,  NULL,         band_lgaus);
    tfi_plot_function_pair("ER", 1,         lopa_sgerf,   band_sgerf);
    tfi_plot_function_pair("BQ", 1,         lopa_biqua,   band_biqua);
    tfi_plot_function_pair("GA", 1,         lopa_gauss,   band_gauss);
    tfi_plot_function_pair("BW", ord_buttw, lopa_buttw,   band_buttw);

    return;

    double band_lgaus(double f, int32_t np)
      { double g = neuromat_filter_bandpass_log_gauss_eval(f, fm_lgaus, np, sigma_lgaus, mag_lgaus);
        return g;
      }

    double lopa_sgerf(double f, int32_t np)
      { double g = neuromat_filter_lowpass_sgerf_eval(f, fm_hi_sgerf, sigma_hi_sgerf);
        return g;
      }

    double band_sgerf(double f, int32_t np)
      { double g = neuromat_filter_bandpass_sgerf_eval(f, fm_lo_sgerf, sigma_lo_sgerf, fm_hi_sgerf, sigma_hi_sgerf);
        return g;
      }

    double lopa_biqua(double f, int32_t np)
      { double g = neuromat_filter_lowpass_biquadratic_eval(f, fhi1, fhi0);
        return g;
      }

    double band_biqua(double f, int32_t np)
      { double g = neuromat_filter_bandpass_biquadratic_eval(f, flo0, flo1, fhi1, fhi0);
        return g;
      }

    double lopa_gauss(double f, int32_t np)
      { double g = neuromat_filter_lowpass_gauss_eval(f, sigma_hi_gauss);
        return g;
      }

    double band_gauss(double f, int32_t np)
      { double g = neuromat_filter_bandpass_gauss_eval(f, sigma_lo_gauss, sigma_hi_gauss, mag_gauss);
        return g;
      }

    double lopa_buttw(double f, int32_t np)
      { double g = neuromat_filter_lowpass_butterworth_eval(f, fs_hi_buttw, np);
        return g;
      }

    double band_buttw(double f, int32_t np)
      { double g = neuromat_filter_bandpass_butterworth_eval(f, fs_lo_buttw, fs_hi_buttw, np);
        return g;
      }
  }
  
void tfi_plot_function_pair(char *tag, int32_t npmax, tfi_function_t *lopa, tfi_function_t *band)
  { 
    for (int32_t which = 0; which < 2; which++)
      { tfi_function_t *gain = (which == 0 ? lopa : band); /* version to plot: low-pass or bandpass. */
        char *xwhich = (which == 0 ? "low" : "band"); 
        fprintf(stderr, "    --- %s tag = %s npmax = %d func = %s pass ---\n", __FUNCTION__, tag, npmax, xwhich);
        if (gain != NULL)
          { char *fname = NULL;
            char *fname = jsprintf("out/%s_%s.txt", tag, (which == 0 ? "lopa" : "band"));
            FILE *wr = open_write(fname, TRUE);
            free(fname);
            tfi_plot_function(wr, npmax, gain);
            fclose(wr);
          }
      }
  }
  
void tfi_plot_function(FILE *wr, int32_t npmax, tfi_function_t *gain) 
  {
    demand(npmax >= 1, "invalid {npmax}");
    double g_tiny = 1.0e-8; /* A small gain. */
    int32_t nf = 30000;  /* Number of Hartley/Fourier terms. */

    int32_t kpmax = (npmax <= 8 ? npmax : 8);
    int32_t np[kpmax+1];     /* Value of {np} for each {kp}. */
    double fg_lo_half[kpmax+1]; /* Freq where the lo shoulder of filter of order {np[kp]} has gain 0.5. */
    double fg_lo_tiny[kpmax+1]; /* Freq where the lo shoulder of filter of order {np[kp]} has gain 0.001. */
    double fg_hi_half[kpmax+1]; /* Freq where the hi shoulder of filter of order {np[kp]} has gain 0.5. */
    double fg_hi_tiny[kpmax+1]; /* Freq where the hi shoulder of filter of order {np[kp]} has gain 0.001. */
    double G_prev[kpmax+1]; /* Prev gain for order {np[kp]}. */
    
    int32_t kfmax = nf/2;
    for (int32_t kp = 0; kp <= kpmax; kp++) 
      { fg_lo_half[kp] = fg_lo_tiny[kp] = NAN;
        fg_hi_half[kp] = fg_hi_tiny[kp] = NAN;
        G_prev[kp] = -INF;
      }
    fprintf(wr, "# npmax = %d\n", npmax); /* For the plot script. */
    for (int32_t kf = 1; kf < kfmax; kf++)
      { double f = ((double)kf)/nf;
        fprintf(wr, "%15.12f", f);
        for (int32_t kp = 1; kp <= kpmax; kp++)
          { np[kp] = ((kp == kpmax) && (npmax > kpmax) ? npmax : kp);
            double G = gain(f, np[kp]);
            fprintf(wr, " %15.12f", G);
            if (G > G_prev[kp])
              { /* Increasing half: */
                if (fabs(G) < 0.5) { fg_lo_half[kp] = f; }
                if (fabs(G) < g_tiny) { fg_lo_tiny[kp] = f; }
              }
            else if (G < G_prev[kp])
              { /* Decreasing half: */
                if ((fabs(G) < 0.5) && (isnan(fg_hi_half[kp]))) { fg_hi_half[kp] = f; }
                if ((fabs(G) < g_tiny) && (isnan(fg_hi_tiny[kp]))) { fg_hi_tiny[kp] = f; }
              }
            G_prev[kp] = G;
          }
        fprintf(wr, "\n");
      }

    for (int32_t kp = 1; kp <= kpmax; kp++)
      { fprintf(stderr, " LO ramp order = %3d:", np[kp]);
        if (! isnan(fg_lo_half[kp])) { fprintf(stderr, " gain is 0.5 at f = %10.6f", fg_lo_half[kp]); }
        if (! isnan(fg_lo_tiny[kp])) { fprintf(stderr, " gain is %12.4e at f = %10.6f", g_tiny, fg_lo_tiny[kp]); }
        fprintf(stderr, "\n");
      }
    fprintf(stderr, "\n");
    for (int32_t kp = 1; kp <= kpmax; kp++)
      { fprintf(stderr, " HI ramp order = %3d:", np[kp]);
        if (! isnan(fg_hi_half[kp])) { fprintf(stderr, " gain is 0.5 at f = %10.6f", fg_hi_half[kp]); }
        if (! isnan(fg_hi_tiny[kp])) { fprintf(stderr, " gain is %12.4e at f = %10.6f", g_tiny, fg_hi_tiny[kp]); }
        fprintf(stderr, "\n");
      }
  }
  
void tfi_test_clear_tiny_gains(int32_t nf, double eps)
  { fprintf(stderr, "--- %s  nf = %d eps = %24.16e ---\n", __FUNCTION__, nf, eps);
    demand(nf >= 4, "bad {NF}");
    
    char *fname =  "out/clear_tiny.txt";
    FILE *wr = open_write(fname, TRUE);
    
    /* Create a vector of conjugate Fourier coeffs of all sizes with random phases, with some zeros: */
    complex Fraw[nf];
    for (int32_t kf = 0; kf < nf; kf++) { Fraw[kf] = 0; }
    double Gmax = 2.0, Gmin = 1.0e-20;
    int32_t nftest = nf/2 - 3;
    for (int32_t kfa = 0; kfa <= nftest; kfa++)
      { int32_t kfb = (nf - kfa) % nf;
        double phk = (kfa == kfb ? 0 : dabrandom(0, M_2_PI));
        complex Fk = Gmax*exp(kfa*(log(Gmin) - log(Gmax))/nftest)*cexp(I*phk);
        Fraw[kfa] = Fk; Fraw[kfb] = conj(Fk);
      }
   
    /* Convert the Fourier coeffs to Hartley form: */
    double Hraw[nf];
    for (int32_t kfa = 0; kfa <= nf/2; kfa++) 
      { int32_t kfb = (nf - kfa) % nf;
        neuromat_filter_fourier_to_hartley(Fraw[kfa], Fraw[kfb], &(Hraw[kfa]), &(Hraw[kfb]));
      }
      
    /* Cleanup the Hartley coeffs: */
    double Hkuk[nf];
    for (int32_t kf = 0; kf < nf; kf++) { Hkuk[kf] = Hraw[kf]; }
    double fsmp = 1.0;
    bool_t verbose = TRUE;
    neuromat_filter_clear_tiny_gains(nf, Hkuk, eps, fsmp, verbose);
    
    /* Convert the cleaned coeffs back to Fourier and compare: */
    complex Fkuk[nf];
    double Emax = -INF;
    double Gtop = NAN;
    for (int32_t kfa = 0; kfa <= nf/2; kfa++) 
      { int32_t kfb = (nf - kfa) % nf;
        neuromat_filter_hartley_to_fourier(Hkuk[kfa], Hkuk[kfb], &(Fkuk[kfa]), &(Fkuk[kfb]));
        assert(fabs(creal(Fkuk[kfa]) - creal(Fkuk[kfb])) < 1.0e-15);
        assert(fabs(cimag(Fkuk[kfa]) + cimag(Fkuk[kfb])) < 1.0e-15);
        double Grawk = cabs(Fraw[kfa]);
        double Gkukk = cabs(Fkuk[kfa]);
        complex Dk = Fkuk[kfa] - Fraw[kfa];
        double Ek = cabs(Dk);
        fprintf(wr, "%5d", kfa);
        fprintf(wr, "  %24.16e %24.16e %24.16e", creal(Fraw[kfa]), cimag(Fraw[kfa]), Grawk);
        fprintf(wr, "  %24.16e %24.16e %24.16e", creal(Fkuk[kfa]), cimag(Fkuk[kfa]), Gkukk);
        fprintf(wr, "  %24.16e %24.16e %24.16e", creal(Dk), cimag(Dk), Ek);
        fprintf(wr, "\n");
        if (Ek > Emax) { Emax = Ek; Gtop = Grawk; }
      }
    fclose(wr);
    fprintf(stderr, "  max correction = %24.16e for original gain %24.16e\n", Emax, Gtop);
  }

void tfi_test_tabulate_hartley_gains(int32_t nf, char *ftype, double fsmp)
  {
    fprintf(stderr, "--- %s  nf = %d  filter = %s  fsmp = %12.8f ---\n", __FUNCTION__, nf, ftype, fsmp);
    double tol = 4.0e-8;

    /* Tabulating parameters: */
    neuromat_filter_gain_t *gain = NULL; /* Filter function to use. */
    double fsup = NAN;           /* Max freq with nonzero gain for folding, or 0. */
    bool_t normalize;    /* True normalizes to unit max gain. */

    /* Filter gain functions and their specific parameters: */

    double delay_amount = NAN;
    auto complex gain_delay(double f);
    
    double ungau_fm = NAN, ungau_sigma = NAN;
    auto complex gain_ungau(double f);
    
    double biqua_fa = NAN, biqua_fb = NAN;
    auto complex gain_biqua(double f);
    
    double sgerf_fm = NAN, sgerf_sigma = NAN;
    auto complex gain_sgerf(double f);
    
    double buttw_fs = NAN; int32_t buttw_ord = INT32_MIN;
    auto complex gain_buttw(double f);

    if (strcmp(ftype, "DL") == 0)
      { /* Delay filter: */
        gain = &gain_delay;
        int32_t kdelay = nf/100 + 3;
        delay_amount = kdelay/fsmp; /* Delay in seconds. */
        fsup = 0.0; /* No {fsmp}-folding */
        normalize = FALSE;  /* No unit-max normalization. */
        fprintf(stderr, "  delay filter, amount = %12.8f s\n", delay_amount);
      }
    else 
      { /* Real filters: */ 
        normalize = TRUE;
        double gain_tiny = 1.0e-6;
        if (strcmp(ftype, "UG") == 0)
          { /* Unsymmetric Gaussian bump on lin freq space with mean below {fsmp/2} but tail well over {fsmp}: */
            gain = &gain_ungau;
            ungau_fm = fsmp/4;
            ungau_sigma = fsmp/3;
            fsup = ungau_fm + 9.0*ungau_sigma;
            fprintf(stderr, "  gauss hump filter  f_avg = %12.8f  f_dev = %12.8f\n", ungau_fm, ungau_sigma);
          }
        else if (strcmp(ftype, "ER") == 0)
          { /* Erf-sigmoid bandpass with midlevel below {fsmp/2} but tail well over {fsmp}: */
            gain = &gain_sgerf;
            sgerf_fm = fsmp/3;
            double fc = 2*fsmp/3;
            sgerf_sigma = neuromat_filter_lowpass_sgerf_compute_sigma(sgerf_fm, fc, gain_tiny);
            fsup = neuromat_filter_lowpass_sgerf_compute_fsup(sgerf_fm, sgerf_sigma);
            fprintf(stderr, "  erf-log-sigmoid filter  fm = %12.8f  sigma = %12.8f\n", sgerf_fm, sgerf_sigma);
          }
        else if (strcmp(ftype, "BQ") == 0)
          { /* Biquadratic sigmoid bandpass with midlevel below {fsmp/2} but tail well over {fsmp}: */
            gain = &gain_biqua;
            biqua_fa = fsmp/3;
            biqua_fb = 2*fsmp/3;
            fsup = biqua_fb;
            fprintf(stderr, "  biquadratic sigmoid filter  fa = %12.8f  fb = %12.8f\n", biqua_fa, biqua_fb);
          }
        else if (strcmp(ftype, "BW") == 0)
          { /* Butterworth bandpass with knee freq below {fsmp/2} but tail well over {fsmp}: */
            gain = &gain_buttw;
            buttw_fs = fsmp/3;
            double fc = 2*fsmp/3;
            buttw_ord = neuromat_filter_lowpass_butterworth_compute_order(buttw_fs, fc, gain_tiny);
            fsup = neuromat_filter_lowpass_butterworth_compute_fsup(buttw_fs, buttw_ord);
            fprintf(stderr, "  Butterworth lowpass filter  fs = %12.8f  ord = %d\n", buttw_fs, buttw_ord);
          }
        else
          { demand(FALSE, "invalid filter type"); }
      }
      
    /* Get the freq index range for plotting/folding: */
    int32_t kf_min, kf_max;
    if (fsup == 0)
      { /* No folding: */
        kf_min = 0; 
        kf_max = nf - 1;
      }
    else
      { kf_max = (int32_t)ceil(fsup/fsmp + 1.0e-8)*nf - 1; 
        kf_min = -kf_max;
      }
    assert(kf_max > 0);
    double f_max = kf_max*fsmp/nf;
    double f_min = kf_min*fsmp/nf;
      
    if ((kf_min != 0) || (kf_max != nf-1))
      { fprintf(stderr, "  folding kf in {%d..%d)", kf_min, kf_max);
        fprintf(stderr, "  = [%12.8f _ %12.8f]", f_min, f_max);
        fprintf(stderr, "  fsup = %12.8f\n", fsup);
      }
    else
      { fprintf(stderr, "  no folding\n"); }

    /* Evaluate the gain without folding or conj-mirroring: */
    int32_t nf_eval = (kf_max - kf_min + 1);
    complex *F_org = talloc(nf_eval, complex); /* Original Fourier filter, unfolded. */
    for (int32_t kf = kf_min; kf <= kf_max; kf++)
      { double f = ((double)kf)*fsmp/nf; 
        F_org[kf - kf_min] = gain(f);
      }

    /* Tabulate the gains in Hartley form without normalization: */
    double *H_raw = talloc(nf, double);   /* Tabulated Hartley before normalization. */
    bool_t tab_verbose = (nf <= 20);
    neuromat_filter_tabulate_hartley_gains(nf, fsmp, gain, fsup, FALSE, H_raw, tab_verbose);

    /* Compare with the original gains: */
    int32_t nerr = 0;
    for (int32_t kfa = 0; kfa <= nf/2; kfa++)
      { int32_t kfb = (nf - kfa) % nf;
      
        double fa = ((double)kfa)*fsmp/nf;
        double fb = ((double)kfb)*fsmp/nf;
      
        complex Fa_sum = 0.0;
        int32_t ifa = kfa, nfa = 0; while (ifa - nf >= kf_min) { ifa = ifa - nf; }
        while (ifa <= kf_max) { Fa_sum += F_org[ifa - kf_min]; ifa = ifa + nf; nfa++; }
        
        complex Fb_sum = 0;
        int32_t ifb = kfb, nfb = 0; while (ifb - nf >= kf_min) { ifb = ifb - nf; }
        while (ifb <= kf_max) { Fb_sum += F_org[ifb - kf_min]; ifb = ifb + nf; nfb++; }
        
        if (fsmp == 0) { assert((nfa == 1) && (nfb == 1)); }
        
        complex Fa_exp = 0.5*(Fa_sum + conj(Fb_sum));
        complex Fb_exp = conj(Fa_exp);
        
        complex Fa_tab, Fb_tab;
        neuromat_filter_hartley_to_fourier(H_raw[kfa], H_raw[kfb], &Fa_tab, &Fb_tab);
        double ea = cabs(Fa_exp - Fa_tab);
        double eb = cabs(Fb_exp - Fb_tab);
        if ((ea > tol) || (eb > tol))
          { fprintf(stderr, "  !! error\n"); 
            fprintf(stderr, "    fa = %14.8f (%d) fb = %14.8f (%d)\n", fa, kfa, fb, kfb); 
            fprintf(stderr, "    Fa_exp = %24.16e + %24.16e * I\n", creal(Fa_exp), cimag(Fa_exp));
            fprintf(stderr, "    Fb_exp = %24.16e + %24.16e * I\n", creal(Fb_exp), cimag(Fb_exp));
            fprintf(stderr, "    H[fa] = %24.16e  H[fb] = %24.16e\n", H_raw[kfa], H_raw[kfb]);
            fprintf(stderr, "    Fa_tab = %24.16e + %24.16e * I  ea = %10.2e\n", creal(Fa_tab), cimag(Fa_tab), ea);
            fprintf(stderr, "    Fb_tab = %24.16e + %24.16e * I  eb = %10.2e\n", creal(Fb_tab), cimag(Fb_tab), eb);
            fprintf(stderr, "\n");
            nerr++;
          }
      }
    demand(nerr == 0, "{neuromat_filter_tabulate_hartley_gains} failed");
    
    double *H_kuk = NULL;
    if (normalize)
      { /* Tabulate again with unit-max norm and eps-cleaning: */
        H_kuk = talloc(nf, double);
        neuromat_filter_tabulate_hartley_gains(nf, fsmp, gain, fsup, normalize, H_kuk, FALSE);
      }

    if (nf >= 100)
      { /* Write file for plot: */
        char *fname = jsprintf("out/%s_tab_%04d.txt", ftype, nf);
        FILE *wr = open_write(fname, TRUE);
        free(fname);
        fprintf(wr, "# fmin = %16.12f\n", ((double)kf_min)*fsmp/nf);
        fprintf(wr, "# fmax = %16.12f\n", ((double)kf_max)*fsmp/nf);
        for (int32_t kf = kf_min; kf <= kf_max; kf++)
          { double fk = ((double)kf)*fsmp/nf; 
            complex Fk_org = F_org[kf - kf_min];
            double Hk_raw, Hk_kuk;
            complex Fk_rec;
            if ((kf >= 0) && (kf < nf))
              { int32_t jf = (nf - kf) % nf;
                Hk_raw = H_raw[kf]; 
                complex Fj_rec;
                if (H_kuk != NULL)
                  { Hk_kuk = H_kuk[kf];
                    double Hj_kuk = H_kuk[jf];
                    neuromat_filter_hartley_to_fourier(Hk_kuk, Hj_kuk, &Fk_rec, &Fj_rec);
                  }
                else
                  { Hk_kuk = -100;
                    double Hj_raw = H_raw[jf];
                    neuromat_filter_hartley_to_fourier(Hk_raw, Hj_raw, &Fk_rec, &Fj_rec);
                  }
              }
            else
              { Hk_raw = Hk_kuk = -100;
                Fk_rec = -100 - 100*I;
              }
            fprintf(wr, "%+4d %+20.15f", kf, fk);
            fprintf(wr, "  %+20.15f %+20.15f", creal(Fk_org), cimag(Fk_org));
            fprintf(wr, "  %+20.15f", Hk_raw);
            fprintf(wr, "  %+20.15f", Hk_kuk);
            fprintf(wr, "  %+20.15f %+20.15f", creal(Fk_rec), cimag(Fk_rec));
            fprintf(wr, "\n");
          }
        fclose(wr);
      }
    
    free(F_org);
    free(H_raw);
    if (H_kuk != NULL) { free(H_kuk); }
    
    return;
    
    complex gain_delay(double f)
      { return cexp(-2*M_PI*I*f*delay_amount);
      }
      
    complex gain_ungau(double f)
      { double z = (f - ungau_fm)/ungau_sigma;
        return (complex)(exp(-0.5*z*z));
      }

    complex gain_sgerf(double f)
      { double g = neuromat_filter_lowpass_sgerf_eval(f, sgerf_fm, sgerf_sigma);
        return (complex)g;
      }

    complex gain_biqua(double f)
      { double g = neuromat_filter_lowpass_biquadratic_eval(f, biqua_fa, biqua_fb);
        return (complex)g;
      }

    complex gain_buttw(double f)
      { double g = neuromat_filter_lowpass_butterworth_eval(f, buttw_fs, buttw_ord);
        double z = 2.0*buttw_ord*(log(f) - log(buttw_fs) + log(M_2_PI));
        if (f > fsmp/2) { fprintf(stderr, "      gain_buttw: f = %12.8f  fs = %12.8f  z = %12.8f  g = %12.8f\n", f, buttw_fs, z, g); }
        return (complex)g;
      }

  }
