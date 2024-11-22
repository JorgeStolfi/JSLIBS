/* text_intg - test ODE integrators. */
/* Last edited on 2023-03-29 19:32:02 by stolfi */

#define _GNU_SOURCE
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <jsfile.h>
#include <jsstring.h>
#include <rn.h>
#include <vec.h>

#include <intg_problem.h>
#include <intg_gen.h>
#include <intg_RKF4.h>
#include <intg_RKF2.h>
#include <intg_Euler.h>

typedef intg_prob_t Problem;
  
typedef struct Grater
  { char *tag;
    Intg_T *i; 
  } Grater;
  
typedef struct TestResults
  { int32_t nEvals;     /* Number of RHS evaluations */
    int32_t nSteps;     /* Number of steps really taken */
    int32_t nRedos;     /* Number of steps discarded because of large error */
    double tf;      /* Stopping time */
    double ec, ek;  /* Estimated integration error was approximately {c*dt^k} */
    double mc, mk;  /* Measured integration error was approximately {c*dt^k} */
  } TestResults;

typedef struct Stats
  { double xx, xy, xu, yu, uu;  /* Error/stepsize sums */
  } Stats;
    /* Statistics used for the least-squares fit of integration error.
      Let {dt[i]} be the length of step {i}, {error[i]} be the 
      corresponding error, {x[i]} be {log(dt[i])},
      and {y[i]} be {log(error[i])}.  Then 
        {xx == sum(x[i]*x[i])},
        {xy == sum(x[i]*y[i])},
        {xu == sum(x[i])},
        {yu == sum(y[i])},
        {uu == sum(1)},
      summed over all steps actually taken. */
  
/* INTERNAL PROTOTYPES */

int32_t main(int32_t argn, char **argc);
void DoErrorTest(Grater *g, intg_prob_vec_t *p);
void ErrorTest
  ( Grater *g,       /* Integrator */
    intg_prob_t *p, /* Problem */
    FILE *wr,        /* Plot file, or NULL */
    Time dtSkip,     /* Time step for integration. */
    bool_t ms,         /* TRUE for measured errors, FALSE for estimated ones. */
    double *ek       /* Error grows proportionally to {dt^ek}. */
  );
void WriteErrorTestPoint(FILE *wr, double x, double y);
void WriteErrorTestPlotHeader(FILE *wr, char *iName, intg_prob_vec_t *p, bool_t ms);
void WriteErrorTestIdealPoints(FILE *wr, double ek);
void Trace(char c);

/* IMPLEMENTATIONS */

int32_t main(int32_t argn, char **argc)
  { intg_prob_vec_t p = intg_prob_t_Sample();
    int32_t ng = 3;
    Grater g[ng];
    g[0] = (Grater){"euler", (Intg_T *)Intg_Euler_new()};
    g[1] = (Grater){"rkfo2", (Intg_T *)Intg_RKF2_new()};
    g[2] = (Grater){"rkfo4", (Intg_T *)Intg_RKF4_new()};
    /* Plot true solutions of all problems: */
    int32_t j;
    /* Make plot of error variation with step size: */
    for (j = 0; j < ng; j++)
      { DoErrorTest(&(g[j]), &p);
        fprintf(stdout, "\n");
      }
    return 0;
  }

void DoErrorTest(Grater *g, intg_prob_vec_t *p)
  { double dtSkip = 1.0e-2;
    bool_t ms; /* TRUE for measured error, FALSE for estimated error. */
    for (ms = FALSE; ms <= TRUE; ms++)
      { char *msTag = (ms ? "merr" : "eerr") ;
        int32_t i;
        for (i = 0; i < p->ne; i++)
          { intg_prob_t *pi = &(p->e[i]);
            char *fname = jsprintf("out/%s-%s-%s.plot", pi->tag, g->tag, msTag);
            FILE *wr = open_write(fname, TRUE);
            free(fname);
            double ek;
            WriteErrorTestPlotHeader(wr, g->i->descr, p, ms);
            ErrorTest(g, pi, wr, dtSkip, ms, &ek); 
            WriteErrorTestIdealPoints(wr, ek);
            fclose(wr);
          }
      }
  }

void ErrorTest
  ( Grater *g,       /* Integrator */
    intg_prob_t *p, /* Problem */
    FILE *wr,        /* Plot file, or NULL */
    Time dtSkip,     /* Time step for integration. */
    bool_t ms,         /* TRUE for measured errors, FALSE for estimated ones. */
    double *ek       /* Error grows proportionally to {dt^ek}. */
  )
  { int32_t NSizes = 11;
    double dtMin = 1.0e-6;
    double dtMax = 1.0e-1;
    double dtRef = 1.0e-4;
    
    int32_t n = p->n; /* Dimension of state vector. */
    
    double sumSlope = 0.0;
    int32_t nSlope = 0;
    
    auto void GatherErrorStats(double dt, double error, double errorRef);
    void GatherErrorStats(double dt, double error, double errorRef)
      { if((errorRef > 1.0e-8) && (error > 1.0e-8))
          { double rdt = dt/dtRef;
            double rer = error/errorRef;
            double x = 0.5 * log(rdt*rdt + 1.0e-200)/log(10.0);
            double y = 0.5 * log(rer*rer + 1.0e-200)/log(10.0);
            WriteErrorTestPoint(wr, x, y);
            if (fabs(x) > 0.1) { sumSlope += y/x; nSlope++; }
          }
      }

    /* Work areas for integration: */
    double_vec_t sa = double_vec_new(n);
    double_vec_t va = double_vec_new(n);
    double_vec_t sr = double_vec_new(n);
    double_vec_t er = double_vec_new(n);
    double_vec_t st = double_vec_new(n);
    
    /* Informative header: */
    fprintf
      ( stderr, 
        "=== %s %s [%8.4f + %8.4f] dt = %12.4e ===\n", 
        g->tag, p->tag, p->t0, p->t1, p->dt
      );
    
    auto double GetError(Time ta, State *sa, Velocity *va, Time dt);
      /* Gets error (measured or estimated) for step {dt} statring 
        at time {ta}, state {sa}, vel {va}. */
      
    double GetError(Time ta, State *sa, Velocity *va, Time dt)
      { /* Take reference step, generate state {sr}: */
        Time tr = ta + dt; 
        g->i->step(g->i, p->rhs, ta,sa,va, tr,&sr, &er);
        if (ms)
          { /* Compute state {s1} on true solution: */
            p->sol(tr, &st);
            return rn_dist(n, sr.e, st.e);
          }
        else
          { /* Return estimated error: */
            return rn_norm(n, er.e);
          }
      }
    
    /* Fixed-step integration: */
    Time ta = p->t0;
    while(ta < p->t1)
      { /* Force current state to true solution: */
        p->sol(ta, &sa);
        
        /* Compute current velocity: */
        p->rhs(ta, &sa, &va);
        
        /* Compute errors for reference stepsize: */
        double errorRef = GetError(ta,&sa,&va, dtRef);
        
        /* Compute and tally errors for various stepsizes: */
        int32_t k;
        for (k = 0; k < NSizes; k++)
          { double e = (NSizes == 1 ? 0.5 : ((double)k)/((double)NSizes-1));
            double dt = dtMin * exp(log(dtMax/dtMin)*e);
            if (ta+dt <= p->t1)
              { double error = GetError(ta,&sa,&va, dt);
                GatherErrorStats(dt, error, errorRef);
              }
          }
        
        /* Prepare for next step: */
        ta += dtSkip;
      }

    (*ek) = sumSlope/((double)nSlope);
  }

void WriteErrorTestPoint(FILE *wr, double x, double y)
  { fprintf(wr, "  %12.8f %12.8f\n", x, y); }

void WriteErrorTestPlotHeader(FILE *wr, char *gName, intg_prob_vec_t *p, bool_t ms)
  { fprintf(wr, "# integrator = %s\n", gName);
    fprintf(wr, "# problems =");
    int32_t i;
    for (i = 0; i < p->ne; i++) { fprintf(wr, " %s", p->e[i].tag); }
    fprintf(wr, "\n");
    fprintf(wr, "# %12s %12s\n", "dt/dt0", "err/err0");
  }

void WriteErrorTestIdealPoints(FILE *wr, double ek)
  { int32_t NPoints = 51;
    double rdtMin = 1.0e-3;
    double rdtMax = 1.0e+5;
    int32_t k;
    for (k = 0; k <= NPoints; k++) 
      { double e = (NPoints == 1 ? 0.5 : ((double)k)/((double)NPoints-1));
        double rdt = rdtMin * exp(log(rdtMax/rdtMin)*e);
        double rer = exp(ek*log(rdt));
        double x = 0.5 * log(rdt*rdt + 1.0e-200)/log(10.0);
        double y = 0.5 * log(rer*rer + 1.0e-200)/log(10.0);
        WriteErrorTestPoint(wr, x, y);
      }
    fflush(wr);
  }

void Trace(char c)
  { fputc(c, stderr); }
