/* test_sve_blpulse --- test of {sve_minn.h} for bandlimited pulse design */
/* Last edited on 2025-02-16 20:23:58 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <limits.h>

#include <fftw3.h>

#include <bool.h>
#include <sign.h>
#include <argparser.h>
#include <jsrandom.h>
#include <jsstring.h>
#include <affirm.h>
#include <rmxn.h>
#include <rn.h>
#include <jsfile.h>
#include <vec.h>

#include <sve_minn.h>
#include <minn_plot.h>

/* GENERAL PARAMETERS */

typedef struct options_t
  { int32_t ns; /* Number of samples in signal. */
  } options_t;

#define MAXNS 256
  /* Max number of samples. */

/* INTERNAL PROTOTYPES */

int32_t main (int32_t argc, char **argv);

options_t *get_options(int32_t argc, char **argv);
  /* Parses the command-line options. */

void find_compact_pulse(int32_t ns);
  /* Tries to find a pulse that is compact in time and frequency 
    by nonlinear optimization. */ 

void write_solution
  ( char *prefix,
    char *tag,
    sve_goal_t *F,
    int32_t ns,
    double x[],
    double Fx,
    int32_t nf,
    int32_t fMax,
    double P[],
    double W[]
  );
  /* Writes the pulse and its spectrum. Namely, writes to
    "{prefix}-{tag}-a.dat" the pulse samples {x[0..ns-1]}; and to
    "{prefix}-tag-p.dat" its spectrum {P[0..fMax]} and the badness
    weight function {W[0..fmax]}. Also prints the function value {Fx}
    to {stderr} */

void write_pulse(char *prefix, char *tag, int32_t ns, int32_t nf, double x[]);
  /* Writes the pulse samples {x[0..ns-1]} to a file named
    "{fname}". Each line contains a sample index and 
    the sample's value. If {fname} is NULL, writes to {stderr}
    instead. 
    
    The file actually starts with sample index {-2} and continues to
    sample index {nf+1}, inclusive. The samples are assumed to be
    periodic with period {nf}; samples with indices in {ns..nf-1} are
    assumed to be zero. */

void write_spectrum(char *prefix, char *tag, int32_t fMax, double P[], double W[]);
  /* Writes the power spectrum {P[0..fMax]} to a file named
    "{fname}". Each line contains an absolute frequency {f},
    the power {P[f]}, and the badness weight {W[f]}.
    Note that {P} and {F} must have {fMax+1} elements.
    If {fname} is NULL, writes to {stderr}
    instead. */

void pulse_spectrum
  ( int32_t ns,
    int32_t nf,
    double x[], 
    fftw_complex in[], 
    fftw_complex ot[], 
    fftw_plan plan,
    double P[]
  );
  /* Pads the pulse {x[0..ns-1]} with zeros to width {nf}, then
    computes its FFT and reduces it to the extended power spectrum
    {P[0..fMax]} where {fMax == nf/2}. Note that {P} must have
    {fMax+1} elements. Assumes that {in} and {ot} have {nf} complex
    elements. */

double pulse_badness(int32_t fMax, double P[], double W[]);
  /* A badness score that measures the high freq contents of the power
    spectrum {P[0..fMax]}, according to the weights {W[0..fMax]}. Note
    that {P} and {W} must have {fMax+1} elements. */

void fourier_to_spectrum(int32_t nf, fftw_complex ot[], double P[]);
  /* Reduces a complex Fourier transform {ot[0..nf-1][0..1]}
    to a power spectrum {P[0..fMax]} where {fMax == nf/2}. Note
    that {P} must have {fMax+1} elements. */

void plot_spillover
  ( char *prefix,
    char *tag,
    sve_goal_t *F,
    int32_t ns,
    double x[]
  );
  /* Writes a FNI file "{prefix}-{tag}-plt.fni"
    with a plot of the spillover energy on an arbitrary 2D
    subspace through the point {x[0..ns-1]}. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { options_t *o = get_options(argc, argv);
    find_compact_pulse(o->ns);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }
  
void find_compact_pulse(int32_t ns)
  { 
    /* Output file names: */
    char *prefix = NULL; /* Prefix for output file names. */
    char *prefix = jsprintf("out/%03d", ns);
    
    /* Working storage for the goal function: */
    int32_t nf = 2*ns;   /* Number of terms in FFT of expanded pulse. */
    int32_t fMax = nf/2; /* Max frequency in power spectrum of expanded pulse. */
    fftw_complex *in = fftw_malloc(sizeof(fftw_complex)*nf); 
    fftw_complex *ot = fftw_malloc(sizeof(fftw_complex)*nf);
    fftw_plan plan = fftw_plan_dft_1d(nf, in, ot, FFTW_FORWARD, FFTW_ESTIMATE);
    double P[fMax+1];  /* Power spectrum of pulse, indexed by absolute frequency. */
    
    /* Frequency-domain penalty weights {W[0..fMaxEx]}: */
    int32_t nfOk = ns;       /* Max desired nonzero terms in FFT. */
    int32_t fMaxOk = nfOk/2; /* Max desired freq. in nonzero spectrum. */
    double W[fMax+1];
    int32_t f;
    for (f = 0; f <= fMax; f++) 
      { /* double df = (f - fMaxOk)/(fMax - fMaxOk);  */
        /* W[f] = (f <= fMaxOk ? 0 : exp(df*df)-1); */
        W[f] = (f <= fMaxOk ? 0 : 1);
      }
    
    auto double sve_goal(int32_t n, double x[]); 
      /* The goal function for uptimization. */
      
    double Xprev[ns]; /* Guess in previous call of {sve_OK} function. */
    int32_t nok = 0;      /* Counts iterations (actually, calls to {sve_OK}). */
    bool_t sve_debug = FALSE;
    bool_t sve_debug_probes = FALSE;
    
    auto bool_t sve_OK(int32_t iter, int32_t n, double x[], double Fx, double dist, double step, double radius); 
      /* Acceptance criterion function. */

    double x[ns];     /* Initial guess and final solution. */
    
    /* Initialize {x} with a Hann pulse: */
    srand(4615);  srandom(4615);
    int32_t i;
    for (i = 0; i < ns; i++) 
      { double z = 2*M_PI*((double)i)/((double)ns);
        x[i] = 1.0 - cos(z);
      }
    
    /* Print and write the initial guess: */
    fprintf(stderr, "initial guess:\n");
    pulse_spectrum(ns, nf, x, in, ot, plan, P);
    double Fx = sve_goal(ns, x);
    write_solution(prefix, "ini", &sve_goal, ns, x, Fx, nf, fMax, P, W);
    
    /* Optimize iteratively: */
    double *ctr = NULL;
    double dMax = +INFINITY;
    double dBox = FALSE;
    double rMin = 0.000001;
    double rMax = 1000.0;
    double rIni = 0.01;
    double minStep = 0.01*rMin;
    sign_t dir = -1;
    int32_t maxIters = 300;
    
    sve_minn_iterate
      ( ns, &sve_goal, &sve_OK, NULL, 
        x, &Fx, dir, 
        ctr, dMax, dBox, rIni, rMin, rMax, 
        minStep, maxIters, sve_debug, sve_debug_probes
      );
    
    /* Print and write final solution: */
    fprintf(stderr, "final solution:\n");
    pulse_spectrum(ns, nf, x, in, ot, plan, P);
    write_solution(prefix, "fin", &sve_goal, ns, x, Fx, nf, fMax, P, W);
    return;
      
    double sve_goal(int32_t n, double x[])
      { assert(n == ns);
        /* Compute the pulse's power spectrum after padding to {2*n} samples: */
        pulse_spectrum(ns, nf, x, in, ot, plan, P);
        double Fx = pulse_badness(fMax, P, W);
        return Fx;
      }
      
    bool_t sve_OK(int32_t iter, int32_t n, double x[], double Fx, double dist, double step, double radius)
      { assert(n == ns);
        fprintf(stderr, "iteration %d\n", nok);
        if (nok > 0)
          { double d = rn_dist(n, Xprev, x);
            fprintf(stderr, "displacement = %16.10f\n", d);
          }
        if (sve_debug)
          { pulse_spectrum(ns, nf, x, in, ot, plan, P);
            write_solution(NULL, "tmp", &sve_goal, ns, x, Fx, nf, fMax, P, W);
          }
        /* Save guess in {Xprev} for next call: */
        rn_copy(ns, x, Xprev); 
        nok++;
        return FALSE;
      }
  }

void write_solution
  ( char *prefix,
    char *tag,
    sve_goal_t *sve_goal,
    int32_t ns,
    double x[],
    double Fx,
    int32_t nf,
    int32_t fMax,
    double P[],
    double W[]
  )
  {
    /* Print and write the pulse: */
    fprintf(stderr, "pulse = \n");
    write_pulse(NULL, tag, ns, nf, x);
    fprintf(stderr, "\n");

    fprintf(stderr, "spectrum = \n");
    write_spectrum(NULL, tag, fMax, P, W);
    fprintf(stderr, "\n");

    double S = pulse_badness(fMax, P, W);
    fprintf(stderr, "  spillover = %+24.16e\n", S);
    fprintf(stderr, "\n");
    
    double FxN = sve_goal(ns, x);
    fprintf(stderr, "goal function = %+24.16e %+24.16e\n", Fx, FxN);
    demand(Fx == FxN, "inconsistent function value");

    if (prefix != NULL)
      { write_pulse(prefix, tag, ns, nf, x);
        write_spectrum(prefix, tag, fMax, P, W);
        plot_spillover(prefix, tag, sve_goal, ns, x);
      }
  }

void write_pulse(char *prefix, char *tag, int32_t ns, int32_t nf, double x[])
  { char *fname = NULL;
    FILE *wr;
    if (prefix != NULL) 
      { char *fname = jsprintf("%s-%s-a.dat", prefix, tag);
        wr = open_write(fname, TRUE);
      }
    else
      { wr = stderr; }
    int32_t i;
    for (i = -2; i < nf+2; i++)
      { int32_t k = ((i % nf) + nf) % nf;
        double Xi = (k < ns ? x[k] : 0);
        fprintf(wr, "%5d %12.8f\n", i, Xi);
      }
    if ((wr != stderr) && (wr != stdout)) { fclose(wr); }
  }

void write_spectrum(char *prefix, char *tag, int32_t fMax, double P[], double W[])
  { char *fname = NULL;
    FILE *wr;
    if (prefix != NULL) 
      { char *fname = jsprintf("%s-%s-p.dat", prefix, tag);
        wr = open_write(fname, TRUE);
      }
    else
      { wr = stderr; }
    int32_t f;
    for (f = 0; f <= fMax; f++)
      { fprintf(wr, "%5d %12.8f %12.8f\n", f, P[f], W[f]); }
    if ((wr != stderr) && (wr != stdout)) { fclose(wr); }
  }

void plot_spillover
  ( char *prefix,
    char *tag,
    sve_goal_t *sve_goal,
    int32_t ns,
    double x[]
  )
  { /* Choose two orthogonal deformation modes {ua,ub} with mags {ra,rb}: */
    double va[ns], ua[ns];
    for (uint32_t i = 0;  i < ns; i++) { va[i] = 0.05; }
    double ra = rn_dir(ns, va, ua);
    double vb[ns], ub[ns];
    for (uint32_t i = 0;  i < ns; i++) { vb[i] = i*0.05/ns; }
    double dba = rn_dot(ns, vb, ua);
    rn_mix_in(ns, -dba, ua, vb);
    double rb = rn_dir(ns, vb, ub);
    /* Plot the energy as image: */
    double step = fmax(ra, rb)/30;
    float_image_t *img = minn_plot_2D_float_image(ns, x, ua, ra, ub, rb, TRUE, step, sve_goal);
    char *fname = jsprintf("%s-%s-plt.fni", prefix, tag);
    FILE *wr = open_write(fname, TRUE);
    float_image_write(wr, img);
    fclose(wr);
    free(fname);
  }

void pulse_spectrum
  ( int32_t ns, 
    int32_t nf,
    double x[], 
    fftw_complex in[], 
    fftw_complex ot[], 
    fftw_plan plan, 
    double P[]
  )
  { /* The input pulse is {x[0..ns-1]} padded with {nf-ns} zeros: */
    int32_t i;
    for (i = 0; i < nf; i++) 
      { in[i][0] = (i < ns ? x[i] : 0);
        in[i][1] = 0;
      }
    /* Compute the Fourier transform: */
    fftw_execute(plan);
    /* Extract the power spectrum from FFT. */
    fourier_to_spectrum(nf, ot, P);
  }

double pulse_badness(int32_t fMax, double P[], double W[])
  { /* Compute the total power {sumP} and weighted total power {sumPW}: */
    double sumP = 0;
    double sumPW = 0;
    int32_t f;
    for (f = 0; f <= fMax; f++) { sumP += P[f]; sumPW += P[f]*W[f]; }
    /* Compute high-frequency content penalty {hfcp}, independent of absolute power: */
    double hfcp = sumPW / sumP;
    /* Compute the penalty {nupp} for non-unitary power: */
    double nupp = (1 - sumP)*(1 - sumP);
    return hfcp + nupp;
  }

void fourier_to_spectrum(int32_t nf, fftw_complex ot[], double P[])
  {  int32_t fMax = nf/2;  /* Max frequency */
     
     auto double cmod2(/* fftw_complex */ double c[]); /* Complex modulus of {c}, squared. */
     
     double cmod2(/* fftw_complex */ double c[]) { return c[0]*c[0] + c[1]*c[1]; }
     
     /* We must combine Fourier terms with same freq. */
     /* They are {ot[i]} and {ot[n-i]}, modulo {n}. */
     /* Frequency 0 has a single term: */
     P[0] = cmod2(ot[0])/nf;
     /* Other freqs have two terms, except possibly for {fMax}: */
     int32_t f;
     for (f = 1; f <= fMax; f++) 
       { double Pi = cmod2(ot[f]); 
         int32_t j = nf - f;
         if (j != f) { Pi += cmod2(ot[j]); }
         P[f] = Pi/nf;
       }
  }
  
options_t *get_options(int32_t argc, char **argv)
  { argparser_t *pp = argparser_new(stderr, argc, argv);
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem"); 
    o->ns = (int32_t)argparser_get_next_int(pp, 1, MAXNS);
    argparser_finish(pp);
    return o;
  }
