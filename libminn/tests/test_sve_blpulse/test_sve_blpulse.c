/* test_sve_blpulse --- test of {sve_minn.h} for bandlimited pulse design */
/* Last edited on 2017-03-13 22:12:51 by stolfilocal */

#define _GNU_SOURCE
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
#include <rmxn_extra.h>
#include <rn.h>
#include <jsfile.h>
#include <vec.h>

#include <sve_minn.h>
#include <minn_plot.h>

/* GENERAL PARAMETERS */

typedef struct options_t
  { int ns; /* Number of samples in signal. */
  } options_t;

#define MAXNS 256
  /* Max number of samples. */

/* INTERNAL PROTOTYPES */

int main (int argc, char **argv);

options_t *get_options(int argc, char **argv);
  /* Parses the command-line options. */

void find_compact_pulse(int ns);
  /* Tries to find a pulse that is compact in time and frequency 
    by nonlinear optimization. */ 

void write_solution
  ( char *prefix,
    char *tag,
    sve_goal_t *F,
    int ns,
    double x[],
    double Fx,
    int nf,
    int fMax,
    double P[],
    double W[]
  );
  /* Writes the pulse and its spectrum. Namely, writes to
    "{prefix}-{tag}-a.dat" the pulse samples {x[0..ns-1]}; and to
    "{prefix}-tag-p.dat" its spectrum {P[0..fMax]} and the badness
    weight function {W[0..fmax]}. Also prints the function value {Fx}
    to {stderr} */

void write_pulse(char *prefix, char *tag, int ns, int nf, double x[]);
  /* Writes the pulse samples {x[0..ns-1]} to a file named
    "{fname}". Each line contains a sample index and 
    the sample's value. If {fname} is NULL, writes to {stderr}
    instead. 
    
    The file actually starts with sample index {-2} and continues to
    sample index {nf+1}, inclusive. The samples are assumed to be
    periodic with period {nf}; samples with indices in {ns..nf-1} are
    assumed to be zero. */

void write_spectrum(char *prefix, char *tag, int fMax, double P[], double W[]);
  /* Writes the power spectrum {P[0..fMax]} to a file named
    "{fname}". Each line contains an absolute frequency {f},
    the power {P[f]}, and the badness weight {W[f]}.
    Note that {P} and {F} must have {fMax+1} elements.
    If {fname} is NULL, writes to {stderr}
    instead. */

void pulse_spectrum
  ( int ns,
    int nf,
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

double pulse_badness(int fMax, double P[], double W[]);
  /* A badness score that measures the high freq contents of the power
    spectrum {P[0..fMax]}, according to the weights {W[0..fMax]}. Note
    that {P} and {W} must have {fMax+1} elements. */

void fourier_to_spectrum(int nf, fftw_complex ot[], double P[]);
  /* Reduces a complex Fourier transform {ot[0..nf-1][0..1]}
    to a power spectrum {P[0..fMax]} where {fMax == nf/2}. Note
    that {P} must have {fMax+1} elements. */

void plot_spillover
  ( char *prefix,
    char *tag,
    sve_goal_t *F,
    int ns,
    double x[]
  );
  /* Writes a FNI file "{prefix}-{tag}-plt.fni"
    with a plot of the spillover energy on an arbitrary 2D
    subspace through the point {x[0..ns-1]}. */

/* IMPLEMENTATIONS */

int main (int argc, char **argv)
  { options_t *o = get_options(argc, argv);
    find_compact_pulse(o->ns);
    fclose(stderr);
    fclose(stdout);
    return (0);
  }
  
void find_compact_pulse(int ns)
  { 
    /* Debugging printout: */
    bool_t debug = TRUE;

    /* Output file names: */
    char *prefix = NULL; /* Prefix for output file names. */
    asprintf(&prefix, "out/%03d", ns);
    
    /* Working storage for the goal function: */
    int nf = 2*ns;   /* Number of terms in FFT of expanded pulse. */
    int fMax = nf/2; /* Max frequency in power spectrum of expanded pulse. */
    fftw_complex *in = fftw_malloc(sizeof(fftw_complex)*nf); 
    fftw_complex *ot = fftw_malloc(sizeof(fftw_complex)*nf);
    fftw_plan plan = fftw_plan_dft_1d(nf, in, ot, FFTW_FORWARD, FFTW_ESTIMATE);
    double P[fMax+1];  /* Power spectrum of pulse, indexed by absolute frequency. */
    
    /* Frequency-domain penalty weights {W[0..fMaxEx]}: */
    int nfOk = ns;       /* Max desired nonzero terms in FFT. */
    int fMaxOk = nfOk/2; /* Max desired freq. in nonzero spectrum. */
    double W[fMax+1];
    int f;
    for (f = 0; f <= fMax; f++) 
      { /* double df = (f - fMaxOk)/(fMax - fMaxOk);  */
        /* W[f] = (f <= fMaxOk ? 0 : exp(df*df)-1); */
        W[f] = (f <= fMaxOk ? 0 : 1);
      }
    
    auto double F(int n, double x[]); 
      /* The goal function for uptimization. */
      
    double F(int n, double x[])
      { assert(n == ns);
        /* Compute the pulse's power spectrum after padding to {2*n} samples: */
        pulse_spectrum(ns, nf, x, in, ot, plan, P);
        double Fx = pulse_badness(fMax, P, W);
        return Fx;
      }
      
    double Xprev[ns]; /* Guess in previous call of {OK} function. */
    int nok = 0;      /* Counts iterations (actually, calls to {OK}). */
    
    auto bool_t OK(int n, double x[], double Fx); 
      /* Acceptance criterion function. */
      
    bool_t OK(int n, double x[], double Fx)
      { assert(n == ns);
        fprintf(stderr, "iteration %d\n", nok);
        if (nok > 0)
          { double d = rn_dist(n, Xprev, x);
            fprintf(stderr, "displacement = %16.10f\n", d);
          }
        if (debug)
          { pulse_spectrum(ns, nf, x, in, ot, plan, P);
            write_solution(NULL, "tmp", &F, ns, x, Fx, nf, fMax, P, W);
          }
        /* Save guess in {Xprev} for next call: */
        rn_copy(ns, x, Xprev); 
        nok++;
        return FALSE;
      }

    double x[ns];     /* Initial guess and final solution. */
    
    /* Initialize {x} with a Hann pulse: */
    srand(4615);  srandom(4615);
    int i;
    for (i = 0; i < ns; i++) 
      { double z = 2*M_PI*((double)i)/((double)ns);
        x[i] = 1.0 - cos(z);
      }
    
    /* Print and write the initial guess: */
    fprintf(stderr, "initial guess:\n");
    pulse_spectrum(ns, nf, x, in, ot, plan, P);
    double Fx = F(ns, x);
    write_solution(prefix, "ini", &F, ns, x, Fx, nf, fMax, P, W);
    
    /* Optimize iteratively: */
    double dMax = +INFINITY;
    double dBox = FALSE;
    double rMin = 0.000001;
    double rMax = 1000.0;
    double rIni = 0.01;
    double stop = 0.1*rMin;
    sign_t dir = -1;
    int maxIters = 200;
    
    sve_minn_iterate(ns, &F, &OK, x, &Fx, dir, dMax, dBox, rIni, rMin, rMax, stop, maxIters, debug);
    
    /* Print and write final solution: */
    fprintf(stderr, "final solution:\n");
    pulse_spectrum(ns, nf, x, in, ot, plan, P);
    write_solution(prefix, "fin", &F, ns, x, Fx, nf, fMax, P, W);
  }

void write_solution
  ( char *prefix,
    char *tag,
    sve_goal_t *F,
    int ns,
    double x[],
    double Fx,
    int nf,
    int fMax,
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
    
    double FxN = F(ns, x);
    fprintf(stderr, "goal function = %+24.16e %+24.16e\n", Fx, FxN);
    demand(Fx == FxN, "inconsistent function value");

    if (prefix != NULL)
      { write_pulse(prefix, tag, ns, nf, x);
        write_spectrum(prefix, tag, fMax, P, W);
        plot_spillover(prefix, tag, F, ns, x);
      }
  }

void write_pulse(char *prefix, char *tag, int ns, int nf, double x[])
  { char *fname = NULL;
    FILE *wr;
    if (prefix != NULL) 
      { asprintf(&fname, "%s-%s-a.dat", prefix, tag);
        wr = open_write(fname, TRUE);
      }
    else
      { wr = stderr; }
    int i;
    for (i = -2; i < nf+2; i++)
      { int k = ((i % nf) + nf) % nf;
        double Xi = (k < ns ? x[k] : 0);
        fprintf(wr, "%5d %12.8f\n", i, Xi);
      }
    if ((wr != stderr) && (wr != stdout)) { fclose(wr); }
  }

void write_spectrum(char *prefix, char *tag, int fMax, double P[], double W[])
  { char *fname = NULL;
    FILE *wr;
    if (prefix != NULL) 
      { asprintf(&fname, "%s-%s-p.dat", prefix, tag);
        wr = open_write(fname, TRUE);
      }
    else
      { wr = stderr; }
    int f;
    for (f = 0; f <= fMax; f++)
      { fprintf(wr, "%5d %12.8f %12.8f\n", f, P[f], W[f]); }
    if ((wr != stderr) && (wr != stdout)) { fclose(wr); }
  }

void plot_spillover
  ( char *prefix,
    char *tag,
    sve_goal_t *F,
    int ns,
    double x[]
  )
  { int NS = 30; /* Number of pixels on each side of optimum. */
    int i;
    double xa[ns], xb[ns];
    for (i = 0; i < ns; i++) 
      { xa[i] = x[i] + 0.2*(1 - cos((i+0.5)*2*M_PI/ns));
        xb[i] = x[i] + 0.2*(1 - cos((i+0.5)*3*M_PI/ns));
      }
    char *fname = NULL;
    asprintf(&fname, "%s-%s-plt.fni", prefix, tag);
    FILE *wr = open_write(fname, TRUE);
    minn_plot_2D_float_image(wr, ns, F, x, xa, xb, NS);
    fclose(wr);
    free(fname);
  }

void pulse_spectrum
  ( int ns, 
    int nf,
    double x[], 
    fftw_complex in[], 
    fftw_complex ot[], 
    fftw_plan plan, 
    double P[]
  )
  { /* The input pulse is {x[0..ns-1]} padded with {nf-ns} zeros: */
    int i;
    for (i = 0; i < nf; i++) 
      { in[i][0] = (i < ns ? x[i] : 0);
        in[i][1] = 0;
      }
    /* Compute the Fourier transform: */
    fftw_execute(plan);
    /* Extract the power spectrum from FFT. */
    fourier_to_spectrum(nf, ot, P);
  }

double pulse_badness(int fMax, double P[], double W[])
  { /* Compute the total power {sumP} and weighted total power {sumPW}: */
    double sumP = 0;
    double sumPW = 0;
    int f;
    for (f = 0; f <= fMax; f++) { sumP += P[f]; sumPW += P[f]*W[f]; }
    /* Compute high-frequency content penalty {hfcp}, independent of absolute power: */
    double hfcp = sumPW / sumP;
    /* Compute the penalty {nupp} for non-unitary power: */
    double nupp = (1 - sumP)*(1 - sumP);
    return hfcp + nupp;
  }

void fourier_to_spectrum(int nf, fftw_complex ot[], double P[])
  {  int fMax = nf/2;  /* Max frequency */
     auto double cmod2(fftw_complex c); /* Complex modulus of {c}, squared. */
     double cmod2(fftw_complex c) { return c[0]*c[0] + c[1]*c[1]; }
     /* We must combine Fourier terms with same freq. */
     /* They are {ot[i]} and {ot[n-i]}, modulo {n}. */
     /* Frequency 0 has a single term: */
     P[0] = cmod2(ot[0])/nf;
     /* Other freqs have two terms, except possibly for {fMax}: */
     int f;
     for (f = 1; f <= fMax; f++) 
       { double Pi = cmod2(ot[f]); 
         int j = nf - f;
         if (j != f) { Pi += cmod2(ot[j]); }
         P[f] = Pi/nf;
       }
  }
  
options_t *get_options(int argc, char **argv)
  { argparser_t *pp = argparser_new(stderr, argc, argv);
    options_t *o = notnull(malloc(sizeof(options_t)), "no mem"); 
    o->ns = (int)argparser_get_next_int(pp, 1, MAXNS);
    argparser_finish(pp);
    return o;
  }
