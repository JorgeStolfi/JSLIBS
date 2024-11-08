/* test_polygauss --- test of {sve_minn.h} for flat-topped sum of gaussians */
/* Last edited on 2024-11-08 09:52:39 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <limits.h>

#include <bool.h>
#include <sign.h>
#include <argparser.h>
#include <jsrandom.h>
#include <affirm.h>
#include <rn.h>
#include <jsfile.h>
#include <vec.h>

#include <sve_minn.h>
#include <minn_plot.h>

/* INTERNAL PROTOTYPES */

int32_t main (int32_t argc, char **argv);

double tpg_eval_hump(double x, double avg, double dev, double mag);
  /* Evaluates at argument {x} the Gaussian hump with mean {avg}, deviation {dev}, 
    and max value {mag}. */

double tpg_eval_train_simple(double x, int32_t np, double avg[], double dev[], double mag[]);
double tpg_eval_train_serial(double x, int32_t np, double avg[], double dev[], double mag[]);
double tpg_eval_train(double x, int32_t np, bool_t serial, double avg[], double dev[], double mag[]);
  /* Evaluates at argument {x} a train of {np} Gaussian pulses with
    means {avg0..np-1]}, deviations {dev[0..np-1]}, and magnitudes
    {mag[0..np-1]}, sampled at {ns} samples.
    
    If {serial} is false, the result is  simply the sum of the humps.  If {serial} is true,
    each pair of symmetric humps is multiplied by the complement of the combined previous humps. */

double tpg_train_badness(int32_t np, bool_t serial, double avg[], double dev[], double mag[], bool_t verbose);
  /* A badness score that measures the non-uniformity 
    of {F(x)=tpg_eval_train(np,delta,x)} where it is supposed to be flat.
    Specifically, it is the average the curvature squared, weigthed
    {w(x)=tpg_weight(x,np)}. */

double tpg_weight(double x, int32_t np);
  /* A weight function that is mostly 1, with soft shoulders
    between 0 and {h} and {1-h} and 1, where {h = 0.5/(np-1)}. */

typedef struct tpg_sym_t
  { bool_t fixavg;   /* Hump means are fixed. */
    bool_t samedev;  /* All deviations are equal. */
    bool_t fixmag;   /* All magnitudes are fixed. */
  } tpg_sym_t;
  /* Simmetries and other constraints on hump parameters
    {{avg,dev,mag}[0..np-1]} for {tpg_find_parms}.
  
    If {fixavg} is true, each mean {avg[ip]} is fixed at {ip/(np-1)}.
    If false, the first and last means are fixed at 0 and 1, and
    generally {avg[np-1-ip] = 1 - avg[ip]}.
    
    If {samedev} is true, all deviations are equal. If false,
    the deviations {dev[ip]} and {dev[np-1-ip]} are equal for all {ip}.
    
    If {fixmag} is true, all humps have magnitude 1.0.  If false,
    the first and last are 1.0, and generally {mag[np-1-ip]}
    is equal to {mag[ip]}. */

void tpg_find_parms(int32_t np, bool_t serial, double avg[], double dev[], double mag[], tpg_sym_t *sym);
  /* Finds numerically the parameters {{avg,dev,mag}[0..np-1]} so that
    the function {tpg_train_badness(np,serial,avg,dev,mag)} has the
    the flattest top.  The means {avg[0]} and {avg[np-1]} are
    fixed at 0 and 1, respectively.  Also if {samedev} is true the deviations 
    {dev[0..np-1]} are all equal.  Also if {fixmag} is true 
    the amplitudes {mag[0..np-1]} are all fixed at 1. */ 

void tpg_initialize(int32_t np, bool_t serial, double avg[], double dev[], double mag[], tpg_sym_t *sym);
  /* Initializes the vectors {{avg,dev,mag}[0..np-1]} for {tpg_find_parms}. */ 

void tpg_print_parms(char *tag, int32_t np, bool_t serial, double avg[], double dev[], double mag[]);
  /* Prints the parameters {serial} and {{avg,dev,mag}[0..np-1]} to {stderr}. */ 

void tpg_print_packed_parms(char *tag, int32_t nv, double v[]);
  /* Prints{v[0..nv-1]}, the parameters in packed form, to {stderr}. */ 

void tpg_write_plot(char *tag, int32_t iter, int32_t np, bool_t serial, double avg[], double dev[], double mag[], tpg_sym_t *sym);
  /* Writes a file containing the plot of
    {tpg_eval_train(x,np,serial,avg,dev,mag)} and the {np} individual
    humps.
    
    The file name "out/train_{PP}_{tag}_{NNNN}_s{S}_a{A}_d{D}_m{M}.txt", where
    {PP} is the 2-digit {np}, {NNNN} is the 4-digit {iter} (both
    zero-padded); and {S}, {A}, {D}, and {M} are the values of {serial},
    {sym.fixavg}, {sym.samedev}, and {sym.fixmag} (each either 'F' or 'T'). If {iter} is
    negative, the part "_{NNNN}" of the name is omitted. */
    
void tpg_num_variables(int32_t np, tpg_sym_t *sym, int32_t *nv_avg_P, int32_t *nv_dev_P, int32_t *nv_mag_P);
  /* Number of parameters in {avg[0..np-1]}, {dev[0..np-1]}, and {mag[0..np-1]} in a train of {np} humps that 
    are variable and independent, given the constraints {sym}. */
    
#define tpg_avg_d_MIN         (-0.45)
#define tpg_avg_d_MAX         (+0.45)
  /* The difference between the mean {avg[ip]} of a hump and its nominal
    position {ip/(np-1)} may vary in this range times the nominal hump 
    spacing {s=1/(np-1)}. */

#define tpg_dev_serial_MIN  0.25
#define tpg_dev_serial_MAX  4.00
  /* In the serial model, the deviation of a hump may vary in this
    range, times the nominal hump spacing {s=1/(np-1)}. */

#define tpg_dev_simple_MIN  0.25
#define tpg_dev_simple_MAX  0.75
  /* In the simple (non-serial) model, the deviation of a hump may vary
    in this range, times the nominal hump spacing {s=1/(np-1)}. */

#define tpg_mag_MIN 0.20
#define tpg_mag_MAX 1.40
  /* The magnitude of a hump may vary in this range */

void tpg_pack(int32_t np, bool_t serial, double avg[], double dev[], double mag[], tpg_sym_t *sym, int32_t nv, double v[]);
  /* Copies the variable parameters in {{avg,dev,mag}[0..np-1]} to the
    optimization variables {v[0..nv-1]}. The size {nv} must be the sum
    of the counts returned by {tpg_num_variables(...)}.
    
    The parameters are encoded so that their normal range becoems the range {{-1 _+1]}.  
     
    Let {s} be the mean spacing between averages, namely {1/(np-1)} The
    expected range of the average {avg[ip]}, when it is variable, is
    about {(ip±tpg_avg_RAD)*s}; it is encoded as
    {v[k]=(avg[ip]/s-ip)/tpg_avg_RAD}.
     
    The expected range of the deviation {dev[ip]}, when variable, is
    {[tpg_dev_serial_MIN*s _ tpg_dev_serial_MAX*s]} when {serial} is
    true, and {[tpg_dev_simple_MIN*s _ tpg_dev_simple_MAX*s]} when
    {serial} is false. It is therefore encoded as {v[k] =
    (dev[ip]-dmid)/s/drad} where {dmid} and {drad} are the center and
    half-width of that range.
     
    The normal range of the magnitude {mag[ip]}, when variable, is
    assumed to be {[tpg_mag_MIN _ tpg_mag_MAX]}, is enconded as {v[k] =
    (mag[ip]-mmid)/mrad} where {mmid} and {mrad} are the center and
    half-width of that range. */
    
void tpg_unpack(int32_t nv, double v[], int32_t np, bool_t serial, double avg[], double dev[], double mag[], tpg_sym_t *sym);
  /* Copies the optimization variables {v[0..nv-1]} to the corresponding
    parameters in {{avg,dev,mag}[0..np-1]}, undoing the encoding,
    supplying the fixed ones, and duplicating the symmetric ones, as
    determined by the options {serial,samedev,fixmag}. The size {nv}
    must be the sum of the counts returned by {tpg_num_variables()}. */

void tpg_test_pack_unpack(int32_t np);
  /* Tests {tpg_pack,tpg_unpack} for the given parameters. */

void tpg_test_weight(int32_t np);
  /* Prints out the values of {tpg_weight(x,np)} for a few values between 0 and 1. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { 
    int32_t npmin = 2;
    int32_t npmax = 5;
    for (int32_t np = npmin; np <= npmax; np++)
      { tpg_test_pack_unpack(np);
        tpg_test_weight(np);
    
        double avg[np], dev[np], mag[np];
      
        tpg_sym_t sym = (tpg_sym_t){ .fixavg = TRUE, .samedev = TRUE, .fixmag = TRUE };
        bool_t serial = TRUE;

        /* Initialize {avg,dev,mag} with unidormly spaced train: */
        tpg_initialize(np, serial, avg, dev, mag, &sym);
        tpg_print_parms("ini", np, serial, avg, dev, mag);
        tpg_write_plot("ini", -1, np, serial, avg, dev, mag, &sym);
      
        tpg_find_parms(np, serial, avg, dev, mag, &sym);
        
        tpg_print_parms("fin", np, serial, avg, dev, mag);
        tpg_write_plot("fin", -1, np, serial, avg, dev, mag, &sym);
      }
    return 0;
  }

double tpg_eval_train(double x, int32_t np, bool_t serial, double avg[], double dev[], double mag[])
  { assert(avg[0] == 0.0);
    assert(avg[np-1] == 1.0);
    assert(mag[0] == 1.0);
    assert(mag[np-1] == 1.0);
    double F;
    if (serial)
      { F = tpg_eval_train_serial(x, np, avg, dev, mag); }
    else
      { F = tpg_eval_train_simple(x, np, avg, dev, mag); }
    return F;
  }

double tpg_eval_train_simple(double x, int32_t np, double avg[], double dev[], double mag[])
  { double F = 0;
    for (int32_t ip = 0; ip < np; ip++)
      { double G = tpg_eval_hump(x, avg[ip], dev[ip], mag[ip]);
        F += G;
      }
    return F;
  }

double tpg_eval_train_serial(double x, int32_t np, double avg[], double dev[], double mag[])
  { int32_t hp = np/2;  /* Number of humps minus the central one, if any. */
    double F = ((np%2) == 0 ? 0 : tpg_eval_hump(x, avg[hp], dev[hp], mag[hp]));
    for (int32_t ip = hp-1; ip >= 0; ip--)
      { double G0 = tpg_eval_hump(x, avg[ip], dev[ip], mag[ip]);
        double G1 = tpg_eval_hump(x, avg[np-1-ip], dev[np-1-ip], mag[np-1-ip]);
        double G = (1 - F)*(G0 + G1);
        F += G;
      }
    return F;
  }
   
double tpg_eval_hump(double x, double avg, double dev, double mag)
  { double z = (x - avg)/dev;
    if (fabs(z) < 9) 
      { return mag*exp(-0.5*z*z); }
    else
      { return 0.0; }
  }

double tpg_train_badness(int32_t np, bool_t serial, double avg[], double dev[], double mag[], bool_t verbose)
  { 
    /* Evaluate the filter over the region of possible interest: */
    int32_t nx_gap = 30; /* Intervals per gap between humps. */
    int32_t nx = (np - 1)*nx_gap; /* Sampling intervals in the region. */
    double step = 1.0/(nx - 1); /* Sampling step. */
    double xv[nx+1];
    double Fv[nx+1];
    double Wv[nx+1];
    
    double sum_WF = 0, sum_W = 0;
    for (int32_t ix = 0; ix <= nx; ix++)
      { xv[ix] = ((double)ix)/nx;
        Fv[ix] = tpg_eval_train(xv[ix], np, serial, avg, dev, mag);
        Wv[ix] = tpg_weight(xv[ix], np);
        
        sum_WF += Wv[ix]*Fv[ix];
        sum_W += Wv[ix];
      }
    double avg_F = sum_WF/sum_W;
      
    /* Compute weighted average square curvature {msc_F} of {F[ix_min..ix_max]}: */
    double sum_WK2 = 0;
    for (int32_t ix = 1; ix < nx; ix++)
      { double W = Wv[ix];
        double K = ((Fv[ix+1] - Fv[ix]) - (Fv[ix] - Fv[ix-1]))/(step*step);
        double WK2 = W*K*K;
        if (verbose) 
          { fprintf(stderr, "  x = %12.8f  F(x) = %16.14f W(x) = %12.8f", xv[ix], Fv[ix], Wv[ix]);
            fprintf(stderr, "  K(x) = %16.12f  W(x)*K(x)^2 = %16.11f\n", K, WK2);
          }
        sum_WK2 += WK2;
      }
    double msc_F = sum_WK2/sum_W;

    double bad = msc_F/(avg_F*avg_F);
    return bad;
  }
    
double tpg_weight(double x, int32_t np)
  { double F;
    double h0 = 0.40/(np-1);
    double h1 = 0.45/(np-1);
    if ((x < h0) || (x >= 1.0-h0)) 
      { F = 0.0; }
    else
      { if ((x >= h1) && (x <= 1.0-h1))
          { F = 1.0; }
        else
          { double z = fmin(x-h0, 1-h0-x)/(h1 - h0);
            F = 0.5*(1 - cos(z*M_PI));
          }
      }
    return F*F;
  }
  
void tpg_find_parms(int32_t np, bool_t serial, double avg[], double dev[], double mag[], tpg_sym_t *sym)
  { 
    fprintf(stderr, "--- optimizing parms for {np} = %d", np);
    fprintf(stderr, " serial = %c", "FT"[serial]);
    fprintf(stderr, " fixavg = %c samedev = %c fixmag = %c", "FT"[sym->fixavg], "FT"[sym->samedev], "FT"[sym->fixmag]);
    fprintf(stderr, " ----------------------------------------\n");

    demand(np >= 2, "invalid {np}");
    
    /* Debugging printout: */
    bool_t debug = FALSE;
    
    auto double F(int32_t n, double v[]); 
      /* The goal function for uptimization. */

    /* Working storage for the goal function: */
    int32_t nv_avg, nv_dev, nv_mag;
    tpg_num_variables(np, sym, &nv_avg, &nv_dev, &nv_mag);
    int32_t nv = nv_avg + nv_dev + nv_mag; /* Number of optimization variables. */
    double v[nv];
      
    double vprev[nv];  /* Guess in previous call of {OK} function. */
    int32_t nok = 0;   /* Counts iterations (actually, calls to {OK}). */
    
    auto bool_t OK(int32_t n, double v[], double Fv); 
      /* Acceptance criterion function. */
    
    srand(4615);  srandom(4615);

    /* Convert the initial guess paramaters to the optimization variables {v[0..nv-1]}: */
    tpg_pack(np, serial, avg, dev, mag, sym, nv, v);
    
    /* Print the initial solution, packed: */
    tpg_print_packed_parms("ini", nv, v);
    
    double Fv = F(nv, v);
    
    /* Optimize iteratively: */
    double *ctr = NULL;
    double dMax = 1.0;
    double dBox = TRUE;
    double rMin = 1.0e-8;
    double rMax = 1.0/np;
    double rIni = 0.1/np;
    double minStep = 0.1*rMin;
    sign_t dir = -1;
    int32_t maxIters = 300;
    
    sve_minn_iterate
      ( nv, &F, &OK, NULL, 
        v, &Fv, 
        dir, ctr, dMax, dBox, rIni, rMin, rMax, 
        minStep, maxIters, debug
      );
    
    /* Print the optimum solution, packed: */
    tpg_print_packed_parms("fin", nv, v);
    
    /* Unpack optimum solution: */
    tpg_unpack(nv, v, np, serial, avg, dev, mag, sym);

    return;
      
    double F(int32_t n, double v[])
      { assert(n == nv);
        tpg_unpack(nv, v, np, serial, avg, dev, mag, sym);
        double Fv = tpg_train_badness(np, serial, avg, dev, mag, FALSE);
        return Fv;
      }
      
    bool_t OK(int32_t n, double v[], double Fv)
      { assert(n == nv);
        fprintf(stderr, "  iteration %4d", nok);
        if (nok > 0)
          { double d = rn_dist(n, vprev, v);
            fprintf(stderr, "  F = %16.12f displacement = %16.10f ( ", Fv, d);
            for (int32_t iv = 0; iv < nv; iv++) { fprintf(stderr, " %16.10f", v[iv] - vprev[iv]); }
            fprintf(stderr, " )");
          }
        fprintf(stderr, "\n");
        if (debug)
          { tpg_unpack(nv, v, np, serial, avg, dev, mag, sym);
            tpg_write_plot("tmp", nok, np, serial, avg, dev, mag, sym);
          }
        /* Save guess in {vprev} for next call: */
        rn_copy(nv, v, vprev); 
        nok++;
        return FALSE;
      }
  }

void tpg_initialize(int32_t np, bool_t serial, double avg[], double dev[], double mag[], tpg_sym_t *sym)
  {
    for (int32_t ia = 0; ia < np; ia++) { avg[ia] = ((double)ia)/(np-1); }
    for (int32_t id = 0; id < np; id++) { dev[id] = 0.5/(np-1); }
    for (int32_t im = 0; im < np; im++) { mag[im] = 1.0 ; }
  }

void tpg_num_variables(int32_t np, tpg_sym_t *sym, int32_t *nv_avg_P, int32_t *nv_dev_P, int32_t *nv_mag_P)
  { 
    (*nv_avg_P) = (sym->fixavg ? 0 : np/2 - 1);  /* The means are fixed, or symmetric and variable, minus two. */
    (*nv_dev_P) = (sym->samedev ? 1 : (np+1)/2); /* Deviations are symmetrical, maybe all equal. */
    (*nv_mag_P) = (sym->fixmag ? 0 : (np-1)/2);  /* Magnitudes are all 1 or symmatrical, with two fixed. */
  }

void tpg_pack(int32_t np, bool_t serial, double avg[], double dev[], double mag[], tpg_sym_t *sym, int32_t nv, double v[])
  { 
    int32_t nv_avg, nv_dev, nv_mag;
    tpg_num_variables(np, sym, &nv_avg, &nv_dev, &nv_mag);
    assert(nv == nv_avg + nv_dev + nv_mag);
    
    double s = 1.0/(np-1);

    double avg_d_MIN = s*(sym->fixavg ? 0.0 : tpg_avg_d_MIN);
    double avg_d_MAX = s*(sym->fixavg ? 0.0 : tpg_avg_d_MAX);

    double dev_MIN = s*(serial ? tpg_dev_serial_MIN : tpg_dev_simple_MIN);
    double dev_MAX = s*(serial ? tpg_dev_serial_MAX : tpg_dev_simple_MAX);

    double mag_MIN = (sym->fixmag ? 1.0 : tpg_mag_MIN);
    double mag_MAX = (sym->fixmag ? 1.0 : tpg_mag_MAX);
    
    for (int32_t ip = 0; ip < np; ip++)  
      { assert((dev_MIN - 1.0e-6 <= dev[ip]) && (dev[ip] <= dev_MAX + 1.0e-6));
        assert((mag_MIN - 1.0e-6 <= mag[ip]) && (mag[ip] <= mag_MAX + 1.0e-6));
      }
    
    /* Consistency checks: */
    if (sym->fixavg)
      { /* Means should be fixed and equally spaced: */
        assert(nv_avg == 0);
        for (int32_t ia = 0; ia < np; ia++)
          { assert(fabs(avg[ia] - ((double)ia)/(np-1)) < 1.0e-12); }
     }
    else
      { /* Means should be symmetrical about {0.5} and independent except {avg[0],avg[np-1]}. */
        assert(nv_avg == np/2 - 1);
        assert(avg[0] == 0.0);
        assert(avg[np-1] == 1.0);
        for (int32_t ia = 1; ia <= nv_avg; ia++)
          { assert(fabs(avg[ia] + avg[np-1-ia] - 1) < 1.0e-12); }
     }

    if (sym->samedev)
      { /* All deviations should be equal: */ 
        assert(nv_dev == 1);
        for (int32_t id = 1; id < np; id++) { assert(dev[id] == dev[0]); }
      }
    else
      { /* Deviations should be symmetrical, otherwise independent. */
        assert(nv_dev == (np+1)/2);
        for (int32_t id = 0; id < nv_dev; id++) { assert(dev[id] == dev[np-1-id]); }
      }

    if (sym->fixmag)
      { /* All magnitudes should be 1: */ 
        assert(nv_mag == 0);
        for (int32_t im = 0; im < np; im++) { assert(mag[im] == 1.0); }
      }
    else
      { /* Magnitudes should be symmetrical, else independent except {mag[0],mag[np-1]}. */
        assert(nv_mag == (np-1)/2);
        assert(mag[0] == 1.0);
        assert(mag[np-1] == 1.0);
        for (int32_t im = 1; im <= nv_mag ; im++) { assert(mag[im] == mag[np-1-im]); } 
      }

    /* Copy the independent variables: */
    int32_t iv = 0;  /* Counts optimization variables. */
    double amid = 0.5*(avg_d_MIN + avg_d_MAX);
    double arad = 0.5*(avg_d_MAX - avg_d_MIN);
    for (int32_t ia = 1; ia <= nv_avg; ia++) { assert(iv < nv); v[iv] = (avg[ia] - ia*s - amid)/arad; iv++; }
    double dmid = 0.5*(dev_MIN + dev_MAX);
    double drad = 0.5*(dev_MAX - dev_MIN);
    for (int32_t id = 0; id < nv_dev; id++) { assert(iv < nv); v[iv] = (dev[id] - dmid)/drad; iv++; }
    double mmid = 0.5*(tpg_mag_MIN + tpg_mag_MAX);
    double mrad = 0.5*(tpg_mag_MAX - tpg_mag_MIN);
    for (int32_t im = 1; im <= nv_mag; im++) { assert(iv < nv); v[iv] = (mag[im] - mmid)/mrad; iv++; }
    assert(iv == nv);
  }

void tpg_unpack(int32_t nv, double v[], int32_t np, bool_t serial, double avg[], double dev[], double mag[], tpg_sym_t *sym)
  { 
    bool_t debug = FALSE;
    
    int32_t nv_avg, nv_dev, nv_mag;
    tpg_num_variables(np, sym, &nv_avg, &nv_dev, &nv_mag);
    assert(nv == nv_avg + nv_dev + nv_mag);

    if (debug) { for (int32_t iv = 0; iv < nv; iv++) { fprintf(stderr, "    v[%3d] = %16.12f\n", iv, v[iv]); } }

    double s = 1.0/(np-1);

    double avg_d_MIN = s*(sym->fixavg ? 0.0 : tpg_avg_d_MIN);
    double avg_d_MAX = s*(sym->fixavg ? 0.0 : tpg_avg_d_MAX);
    
    double dev_MIN = s*(serial ? tpg_dev_serial_MIN : tpg_dev_simple_MIN);
    double dev_MAX = s*(serial ? tpg_dev_serial_MAX : tpg_dev_simple_MAX);
    if (debug) { fprintf(stderr, "dev MIN = %16.12f MAX = %16.12f\n", dev_MIN, dev_MAX); }

    double mag_MIN = (sym->fixmag ? 1.0 : tpg_mag_MIN);
    double mag_MAX = (sym->fixmag ? 1.0 : tpg_mag_MAX);

    /* Copy the independent variables: */
    int32_t iv = 0;  /* Counts optimization variables. */
    double amid = 0.5*(avg_d_MIN + avg_d_MAX);
    double arad = 0.5*(avg_d_MAX - avg_d_MIN);
    for (int32_t ia = 1; ia <= nv_avg; ia++) { assert(iv < nv); avg[ia] = ia*s + v[iv]*arad + amid; iv++; }
    double dmid = 0.5*(dev_MIN + dev_MAX);
    double drad = 0.5*(dev_MAX - dev_MIN);
    for (int32_t id = 0; id <  nv_dev; id++) { assert(iv < nv); dev[id] = v[iv]*drad + dmid; iv++; }
    double mmid = 0.5*(mag_MIN + mag_MAX);
    double mrad = 0.5*(mag_MAX - mag_MIN);
    for (int32_t im = 1; im <= nv_mag; im++) { assert(iv < nv); mag[im] = v[iv]*mrad + mmid; iv++; }
    assert(iv == nv);
 
    if (sym->fixavg)
      { /* Means should be fixed and equally spaced. */
        assert(nv_avg == 0);
        for (int32_t ia = 0; ia < np; ia++) { avg[ia] = ((double)ia)/(np-1); }
      }
    else
      { /* Means should be symmetrical about {0.5}, otherwise independent except {avg[0],avg[np-1]}. */
        assert(nv_avg == np/2-1);
        avg[0] = 0.0;
        avg[np-1] = 1.0;
        if ((np%2) == 1) { avg[np/2] = 0.5; }
        for (int32_t ia = 1; ia <= nv_avg; ia++) { avg[np-1-ia] = 1.0 - avg[ia]; }
      }

    if (sym->samedev)
      { /* All deviations should be equal. */
        assert(nv_dev == 1);
        for (int32_t id = 1; id < np; id++) { dev[id] = dev[0]; }
      }
    else
      { /* Deviations should be symmetrical, otherwise independent: */
        assert(nv_dev == (np+1)/2);
        for (int32_t id = 0; id < nv_dev; id++) { dev[np-1-id] = dev[id]; }
      }

    if (sym->fixmag) 
      { /* All magnitudes should be 1: */
        assert(nv_mag == 0);
        for (int32_t im = 0; im < np; im++) { mag[im] = 1.0; }
      }
    else
      { /* Magnitudes should be symmetrical, else independent except {mag[0]}. */
        assert(nv_mag == (np-1)/2);
        mag[0] = 1.0;
        mag[np-1] = 1.0;
        for (int32_t im = 1; im <= nv_mag; im++) { mag[np-1-im] = mag[im]; } 
      }
    
    for (int32_t ip = 0; ip < np; ip++)  
      { if (debug) { fprintf(stderr, "  ip = %3d dev = %16.12f\n", ip, dev[ip]); }
        assert((dev_MIN - 1.0e-6 <= dev[ip]) && (dev[ip] <= dev_MAX + 1.0e-6));
        assert((mag_MIN - 1.0e-6 <= mag[ip]) && (mag[ip] <= mag_MAX + 1.0e-6));
      }
  }

void tpg_print_parms(char *tag, int32_t np, bool_t serial, double avg[], double dev[], double mag[])
  { fprintf(stderr, "  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
    fprintf(stderr, "  (%s) %d parameters, {serial} = %c\n", tag, np, "FT"[serial]);
    for (int32_t ip = 0; ip < np; ip++)
      { fprintf(stderr,  "    %2d %16.12f %16.12f %16.12f\n", ip, avg[ip], dev[ip], mag[ip]); }
    double bad = tpg_train_badness(np, serial, avg, dev, mag, TRUE);
    fprintf(stderr, "  badness = %18.15f\n", bad);
    fprintf(stderr, "  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  }

void tpg_print_packed_parms(char *tag, int32_t nv, double v[])
  { fprintf(stderr, "  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
    fprintf(stderr, "  (%s) %d packed parameters parameters\n", tag, nv);
    for (int32_t iv = 0; iv < nv; iv++)
      { fprintf(stderr,  "    %2d %16.12f\n", iv, v[iv]); }
     fprintf(stderr, "  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
  }

void tpg_write_plot(char *tag, int32_t iter, int32_t np, bool_t serial, double avg[], double dev[], double mag[], tpg_sym_t *sym)
  { 
    /* Compute the sample args {xv[0..nx]} and train values {Fv[0..nx]}, and find maximum {F_max}: */
    int32_t nx = 50*(np - 1 + 14);
    double x_min = avg[0] - 7*dev[0];
    double x_max = avg[np-1] + 7*dev[np-1];
    double xv[nx+1];
    double Fv[nx+1];
    double F_max = -INF;
    for (int32_t ix = 0; ix <= nx; ix++)
      { double r = ((double)ix)/nx;
        xv[ix] = (1-r)*x_min + r*x_max;
        Fv[ix] = tpg_eval_train(xv[ix], np, serial, avg, dev, mag);
        if (Fv[ix] > F_max) { F_max = Fv[ix]; }
      }

    /* Open output file: */
    char *itag = NULL, *fname = NULL;
    if (iter > 0) 
      { asprintf(&itag, "_%04d", iter); }
    else
      { asprintf(&itag, "%s", ""); }
    asprintf
      ( &fname, "out/train_%02d_%s%s_s%c_a%c_d%c_m%c.txt", 
        np, tag, itag, 
        "FT"[serial], 
        "FT"[sym->fixavg], "FT"[sym->samedev], "FT"[sym->fixmag]
      );
    FILE *wr = open_write(fname, TRUE);
    free(itag);
    free(fname);
    
    /* Write train and component values: */
    double scale = 1/F_max;
    for (int32_t ix = 0; ix <= nx; ix++)
      { fprintf(wr, "%3d %16.12f %16.12f", ix, xv[ix], scale*Fv[ix]);
        for (int32_t ip = 0; ip < np; ip++)
          { double g = tpg_eval_hump(xv[ix], avg[ip], dev[ip], mag[ip]);
            fprintf(wr, " %16.12f", scale*g);
          }
        fprintf(wr, "\n");
      }
    fclose(wr);
  }

void tpg_test_pack_unpack(int32_t np)
  { 
    fprintf(stderr, "--- testing pack/unpack for {np} = %d", np);
    fprintf(stderr, " ----------------------------------------\n");
    double avg0[np], dev0[np], mag0[np];
    double avg1[np], dev1[np], mag1[np];
    for (bool_t fixavg = FALSE; fixavg <= TRUE; fixavg++)
      {  for (bool_t samedev = FALSE; samedev <= TRUE; samedev++)
          { for (bool_t fixmag = FALSE; fixmag <= TRUE; fixmag++)
              { tpg_sym_t sym = (tpg_sym_t){ .fixavg = fixavg, .samedev = samedev, .fixmag = fixmag };
                fprintf(stderr, "  fixavg = %c samedev = %c fixmag = %c\n", "FT"[sym.fixavg], "FT"[sym.samedev], "FT"[sym.fixmag]);
                for (bool_t serial = FALSE; serial <= TRUE; serial++)
                  { fprintf(stderr, "    serial = %c", "FT"[serial]);
                    int32_t nv_avg, nv_dev, nv_mag;
                    tpg_num_variables(np, &sym, &nv_avg, &nv_dev, &nv_mag);
                    int32_t nv = nv_avg + nv_dev + nv_mag; /* Number of optimization variables. */
                    fprintf(stderr, " %d indep variables avg = %d dev = %d mag = %d\n", nv, nv_avg, nv_dev, nv_mag);
                    double v0[nv], v1[nv];
                    for (int32_t iv = 0; iv < nv; iv++) { v0[iv] = v1[iv] = dabrandom(-1.0, +1.0); }
                    tpg_unpack(nv, v0, np, serial, avg0, dev0, mag0, &sym);
                    tpg_pack(np, serial, avg0, dev0, mag0, &sym, nv, v1);
                    for (int32_t iv = 0; iv < nv; iv++) 
                      { demand(fabs(v0[iv] - v1[iv]) < 1.0e-12, "pack/unpack failed"); }
                    for (int32_t ip = 0; ip < np; ip++)
                      { avg1[ip] = drandom(); dev1[ip] = drandom(); mag1[ip] = drandom(); }
                    tpg_unpack(nv, v1, np, serial, avg1, dev1, mag1, &sym);
                    for (int32_t ip = 0; ip < np; ip++)
                      { assert(fabs(avg0[ip] - avg1[ip]) < 1.0e-12);
                        assert(fabs(dev0[ip] - dev1[ip]) < 1.0e-12);
                        assert(fabs(mag0[ip] - mag1[ip]) < 1.0e-12);
                      }
                  }
              }
          }
      }
  }

void tpg_test_weight(int32_t np)
  { 
    fprintf(stderr, "--- testing weight function for {np} = %d", np);
    fprintf(stderr, " ----------------------------------------\n");
    int32_t nx = 10*(np-1);
    for (int32_t ix = -1; ix <= nx+1; ix++)
      { double x = ((double)ix)/nx;
        double W = tpg_weight(x, np);
        fprintf(stderr, "  x = %+14.10f W(x) = %16.14f\n", x, W);
      }
    fprintf(stderr, "\n");
  }
  
