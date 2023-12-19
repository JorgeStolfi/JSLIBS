/* Tests the Butterworth filter formulas. */
/* Last edited on 2023-12-18 21:49:31 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
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

void tfi_test_filters(double fsmp, double flo0, double flo1, double fhi1, double fhi0);
  /* Plots various filters with the given parameters. */

void tfi_test_tabulate_filter(int32_t nf, int32_t kfil, double fsmp);
  /* Tests {neuromat_filter_tabulate_hartley_gains} with {nf} frequencies,
    sampling frequency {fsmp}, and filter type {kfil}.  The latter 
    is {0} for a delay filter, {1} for a gaussoid bandpass filter
    that needs folding and tiny-gain cleanup. */

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
    tfi_test_filters(fsmp, flo0, flo1, fhi1, fhi0);
    
    for (int32_t nf = 9; nf <= 10; nf++)
      { for (int32_t kfil = 0; kfil <= 1; kfil++)
          { tfi_test_tabulate_filter(nf, kfil, fsmp);
            tfi_test_tabulate_filter(50*nf, kfil, fsmp);
          }
      }
    
    return 0;
  }
  
void tfi_test_hartley_basis(void)
  { double tol = 1.0e-14;
    for (int32_t nf = 5; nf <= 6; nf++)
      { fprintf(stderr, "  --- %s nf = %d ---\n", __FUNCTION__, nf);
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
      { fprintf(stderr, "  --- %s nf = %d ---\n", __FUNCTION__, nf);
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
      { fprintf(stderr, "  --- %s nf = %d ---\n", __FUNCTION__, nf);
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
    double tiny_sgerf = 1.0e-8;
    double fm_lo_sgerf = sqrt(flo0*flo1);
    double sigma_lo_sgerf = neuromat_filter_lowpass_sgerf_compute_sigma(fm_lo_sgerf, flo0, tiny_sgerf);
    double fm_hi_sgerf = sqrt(fhi0*fhi1);
    double sigma_hi_sgerf = neuromat_filter_lowpass_sgerf_compute_sigma(fm_hi_sgerf, fhi0, tiny_sgerf);
    fprintf(stderr, "erf: fm_lo = %16.12f  sigma_lo = %16.12f", fm_lo_sgerf, sigma_lo_sgerf);
    fprintf(stderr, "  fm_hi = %16.12f  sigma_hi = %16.12f\n", fm_hi_sgerf, sigma_hi_sgerf);
   
    /* double fm_lgauss = sqrt(sqrt(flo0*flo1*fhi1*fhi0)); */ /* Mean frequency. */
    double tiny_lgauss = 1.0e-8;
    double fm_lgauss = fhi1; /* Mean frequency of highest hump. */
    double sigma_lgauss = fabs(neuromat_filter_bandpass_log_gauss_compute_sigma(fm_lgauss, fhi0, tiny_lgauss));
    int32_t np_lgauss = (int32_t)ceil(fmax(0, (log(fhi1) - log(flo1))/(2*sigma_lgauss))) + 1; /* Number of humps. */
    fprintf(stderr, "log_gauss: fm = %16.8f  sigma = %16.8f  np = %d\n", fm_lgauss, sigma_lgauss, np_lgauss);
    
    double tiny_gauss = 1.0e-8;
    double fc_lo_gauss = flo1;
    double sigma_lo_gauss = neuromat_filter_lowpass_gauss_compute_sigma(fc_lo_gauss, tiny_gauss);
    double fc_hi_gauss = fhi0;
    double sigma_hi_gauss = neuromat_filter_lowpass_gauss_compute_sigma(fc_hi_gauss, tiny_gauss);
    fprintf(stderr, "gauss: sigma_lo = %16.12f  sigma_hi = %16.12f\n", sigma_lo_gauss, sigma_hi_gauss);
    
    double tiny_buttw = 1.0e-8;
    double fs_lo_buttw = (flo0 <= 0 ? 0: exp(log(flo0) + 0.05*(log(flo1) - log(flo0))));
    double fs_hi_buttw = (fhi1 >= +INF ? +INF : exp(log(fhi1) + 0.05*(log(fhi0) - log(fhi1))));
    int32_t ord_lo_buttw = neuromat_filter_lowpass_butterworth_compute_order(fs_lo_buttw, flo1, tiny_buttw);
    int32_t ord_hi_buttw = neuromat_filter_lowpass_butterworth_compute_order(fs_hi_buttw, fhi0, tiny_buttw);
    int32_t ord_buttw = (int32_t)imax(ord_lo_buttw, ord_hi_buttw);
    fprintf(stderr, "butterworth: fs_lo = %16.12f  fs_hi = %16.12f  ord = %d\n", fs_lo_buttw, fs_hi_buttw, ord_buttw);

    auto double lopa_gauss(double f, int32_t np);
    auto double band_gauss(double f, int32_t np);

    auto double lopa_sgerf(double f, int32_t np);
    auto double band_sgerf(double f, int32_t np);

    auto double lopa_biqua(double f, int32_t np);
    auto double band_biqua(double f, int32_t np);

    auto double lopa_buttw(double f, int32_t np);
    auto double band_buttw(double f, int32_t np);

    auto double band_lgauss(double f, int32_t np);

    tfi_plot_function_pair("LG", np_lgauss, NULL,         band_lgauss);

    tfi_plot_function_pair("BQ", 1,         lopa_biqua,   band_biqua);
    tfi_plot_function_pair("ER", 1,         lopa_sgerf,   band_sgerf);
    tfi_plot_function_pair("GA", 1,         lopa_gauss,   band_gauss);

    tfi_plot_function_pair("BU", ord_buttw, lopa_buttw,   band_buttw);

    return;

    double band_lgauss(double f, int32_t np)
      { double g = neuromat_filter_bandpass_log_gauss(f, fm_lgauss, np, -sigma_lgauss);
        return g;
      }

    double lopa_sgerf(double f, int32_t np)
      { double g = neuromat_filter_lowpass_sgerf(f, fm_hi_sgerf, sigma_hi_sgerf);
        return g;
      }

    double band_sgerf(double f, int32_t np)
      { double ghi = neuromat_filter_lowpass_sgerf(f, fm_hi_sgerf, sigma_hi_sgerf);
        double glo = neuromat_filter_lowpass_sgerf(f, fm_lo_sgerf, sigma_lo_sgerf);
        double g = ghi*(1 - glo);
        return g;
      }

    double lopa_biqua(double f, int32_t np)
      { double g = neuromat_filter_lowpass_biquadratic(f, fhi1, fhi0);
        return g;
      }

    double band_biqua(double f, int32_t np)
      { double ghi = neuromat_filter_lowpass_biquadratic(f, fhi1, fhi0);
        double glo = neuromat_filter_lowpass_biquadratic(f, flo0, flo1);
        double g = ghi*(1 - glo);
        return g;
      }

    double lopa_gauss(double f, int32_t np)
      { double g = neuromat_filter_lowpass_gauss(f, sigma_hi_gauss);
        return g;
      }

    double band_gauss(double f, int32_t np)
      { double ghi = neuromat_filter_lowpass_gauss(f, sigma_hi_gauss);
        double glo = neuromat_filter_lowpass_gauss(f, sigma_lo_gauss);
        double g = ghi*(1 - glo);
        return g;
      }

    double lopa_buttw(double f, int32_t np)
      { double g = neuromat_filter_lowpass_butterworth(f, fs_hi_buttw, np);
        return g;
      }

    double band_buttw(double f, int32_t np)
      { double ghi = neuromat_filter_lowpass_butterworth(f, fs_hi_buttw, np);
        double glo = neuromat_filter_lowpass_butterworth(f, fs_lo_buttw, np);
        double g = ghi*(1 - glo);
        return g;
      }

  }
  
void tfi_plot_function_pair(char *tag, int32_t npmax, tfi_function_t *lopa, tfi_function_t *band)
  { 
    for (int32_t which = 0; which < 2; which++)
      { tfi_function_t *gain = (which == 0 ? lopa : band); /* version to plot: low-pass or bandpass. */
        if (gain != NULL)
          { char *fname = NULL;
            asprintf(&fname, "out/%s-%s.txt", tag, (which == 0 ? "lopa" : "band"));
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
    double fg_half[kpmax+1]; /* Freq where the filter of order {np[kp]} has gain 0.5. */
    double fg_tiny[kpmax+1]; /* Freq where the filter of order {np[kp]} has gain 0.001. */

    int32_t kfmax = nf/2;
    for (int32_t kp = 0; kp <= kpmax; kp++) { fg_half[kp] = fg_tiny[kp] = NAN; }
    fprintf(wr, "# npmax = %d\n", npmax); /* For the plot script. */
    for (int32_t kf = 1; kf < kfmax; kf++)
      { double f = ((double)kf)/nf;
        fprintf(wr, "%15.12f", f);
        for (int32_t kp = 1; kp <= kpmax; kp++)
          { np[kp] = ((kp == kpmax) && (npmax > kpmax) ? npmax : kp);
            double G = gain(f, np[kp]);
            fprintf(wr, " %15.12f", G);
            if ((fabs(G) < 0.5) && (isnan(fg_half[kp]))) { fg_half[kp] = f; }
            if ((fabs(G) < g_tiny) && (isnan(fg_tiny[kp]))) { fg_tiny[kp] = f; }
          }
        fprintf(wr, "\n");
      }

    for (int32_t kp = 1; kp <= kpmax; kp++)
      { fprintf(stderr, " order = %3d", np[kp]);
        fprintf(stderr, " gain is 0.5 at f = %10.6f", fg_half[kp]);
        fprintf(stderr, " gain is %12.4e at f = %10.6f", g_tiny, fg_tiny[kp]);
        fprintf(stderr, "\n");
      }
  }

void tfi_test_tabulate_filter(int32_t nf, int32_t kfil, double fsmp)
  {
    /* Filter is delay with attenuation if {kfil=0}, gaussoid band if {kfil=1}: */
    demand((kfil >= 0) && (kfil <= 1), "invalid filter type");
    char *ftype = (kfil == 0 ? "DL" : "GB");
    fprintf(stderr, "--- %s  nf = %d  filter = %s  fsmp = %12.8f ---\n", __FUNCTION__, nf, ftype, fsmp);
    double tol = 1.0e-8;

    /* Tabulating parameters: */
    neuromat_filter_gain_t *gain = NULL; /* Filter function to use. */
    int32_t kf_min = INT32_MAX;  /* Min freq index to evaluate. */
    int32_t kf_max = INT32_MIN;  /* Max freq index to eavluat. */
    double fsup;         /* Max freq with nonzero gain for folding, or 0. */
    bool_t normalize;    /* True normalizes to unit max gain. */
    double eps_clean;    /* Threshold for small elem cleanup. */

    /* Filter-specific parameters: */
    double delay_amount = NAN;     /* Delay amount (sec) for delay filter. */
    double gauss_fmean = NAN;   /* Mean freq of gaussoid filter. */
    double gauss_fsigma = NAN;  /* Deviation of gaussoid filter. */

    auto complex gain_delay(double f);
    auto complex gain_gauss(double f);

    if (kfil == 0)
      { /* Delay filter: */
        gain = &gain_delay;
        int32_t kdelay = nf/100 + 3;
        delay_amount = kdelay/fsmp; /* Delay in seconds. */
        kf_min = 0;
        kf_max = nf-1;
        fsup = 0.0; /* No {fsmp}-folding */
        normalize = FALSE;  /* No unit-max normalization. */
        eps_clean = 0; /* No small elem cleanup. */
        fprintf(stderr, "  delay filter, amount = %12.8f s\n", delay_amount);
      }
    else
      { /* Gaussian bump with mean below {fsmp/2} but tail well over {fsmp}: */
        gain = &gain_gauss;
        kf_min = -3*nf + 1;
        kf_max = +3*nf - 1;
        gauss_fmean = fsmp/6;
        double f_min = kf_min*fsmp/nf;
        double f_max = kf_max*fsmp/nf;
        double df = fmin(fabs(f_min - gauss_fmean), fabs(f_max - gauss_fmean));
        double gain_tiny = 1.0e-15;
        gauss_fsigma = df/sqrt(-2*log(gain_tiny));
        fsup = fmax(fabs(f_min), fabs(f_max));
        normalize = TRUE;    /* Apply unit-max normalization. */
        eps_clean = 1.0e-13; /* Clean gains less than this. */
        fprintf(stderr, "  gauss hump filter  f_avg = %12.8f  f_dev = %12.8f\n", gauss_fmean, gauss_fsigma);
      }
    if ((kf_min != 0) || (kf_max != nf-1))
      { fprintf(stderr, "  folding kf in {%d..%d) fsup = %12.8f\n", kf_min, kf_max, fsup); }
    else
      { fprintf(stderr, "  no folding\n"); }
      
    if (eps_clean > 0)
      { fprintf(stderr, "  cleaning gains below %12.4e\n", eps_clean); }
    else
      { fprintf(stderr, "  no small gain cleaning\n"); }

    /* Evaluate the gain without folding or conj-mirroring: */
    int32_t nf_eval = (kf_max - kf_min + 1);
    complex *F_org = talloc(nf_eval, complex); /* Original Fourier filter, unfolded. */
    for (int32_t kf = kf_min; kf <= kf_max; kf++)
      { double f = ((double)kf)*fsmp/nf; 
        F_org[kf - kf_min] = gain(f);
      }

    /* Tabulate the gains in Hartley form without normalization and cleanup: */
    double *H_raw = talloc(nf, double);   /* Tabulated Hartley before normalization, cleanup. */
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
    if (normalize || (eps_clean > 0))
      { /* Tabulate again with unit-max norm and eps-cleaning: */
        H_kuk = talloc(nf, double);
        neuromat_filter_tabulate_hartley_gains(nf, fsmp, gain, fsup, normalize, H_kuk, FALSE);
        if (eps_clean > 0)
          { neuromat_filter_clear_tiny_gains(nf, H_kuk, eps_clean, fsmp, FALSE); }
      }
      
    if (nf >= 100)
      { /* Write file for plot: */
        char *fname = NULL;
        asprintf(&fname, "out/%s_tab_%04d.txt", ftype, nf);
        FILE *wr = open_write(fname, TRUE);
        free(fname);
        fprintf(wr, "# fmin = %16.12f\n", ((double)kf_min)*fsmp/nf);
        fprintf(wr, "# fmax = %16.12f\n", ((double)kf_max)*fsmp/nf);
        for (int32_t kf = kf_min; kf <= kf_max; kf++)
          { double fk = ((double)kf)*fsmp/nf; 
            complex Fk_org = F_org[kf - kf_min];
            double Hk_raw = ((kf >= 0) && (kf < nf) ? H_raw[kf] : -100);
            double Hk_kuk = ((H_kuk != NULL) && (kf >= 0) && (kf < nf) ? H_kuk[kf] : -100);
            fprintf(wr, "%+4d %+20.15f", kf, fk);
            fprintf(wr, "  %+20.15f %+20.15f", creal(Fk_org), cimag(Fk_org));
            fprintf(wr, "  %+20.15f", Hk_raw);
            fprintf(wr, "  %+20.15f", Hk_kuk);
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
      
    complex gain_gauss(double f)
      { double z = (f - gauss_fmean)/gauss_fsigma;
        return exp(-0.5*z*z);
      }
  }
