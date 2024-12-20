/* text_cases - test ODE integrators. */
/* Last edited on 2024-12-05 10:31:33 by stolfi */

#include <intg_Euler.h>
#include <intg_RKF2.h>
#include <intg_RKF4.h>
#include <intg_gen.h>
#include <intg_problem.h>

#include <vec.h>
#include <affirm.h>
#include <rn.h>
#include <jsstring.h>
#include <jsfile.h>

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef intg_prob_t Problem;
  
typedef struct Grater
  { char *tag;
    Intg_T *i; 
  } Grater;
  
typedef struct TestResults
  { int nEvals;     /* Number of RHS evaluations */
    int nSteps;     /* Number of steps really taken */
    int nRedos;     /* Number of steps discarded because of large error */
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

int main(int argn, char **argc);
void DoErrorTest(char *outDir, Grater *g, intg_prob_vec_t *p);
void AllStandardTests
  ( char *outDir, 
    Grater *g, 
    int ng, 
    intg_prob_vec_t *p, 
    char *tag, 
    Dist tol
  );
void DoAdaptiveTest(char *outDir, Grater *g, intg_prob_t *p, char *tag, Dist tol);
void DoStepTest(char *outDir, Grater *g, intg_prob_vec_t *p, Dist tol);
void PlotTrueSolution(char *outDir, intg_prob_t *p, int nSteps);
TestResults AdaptiveTest
  ( char *outDir, 
    Grater *g,       /* Integrator */
    intg_prob_t *p,  /* Problem */
    FILE *wr,        /* Plot file, or NULL */
    Dist tol,        /* Error tolerance */
    bool_t force,    /* Forces path to follow the true solution */
    double rStep     /* Extra stepsize factor after {g.i.adjustStepSize} */
  );
void ErrorTest
  ( Grater *g,       /* Integrator */
    intg_prob_t *p,  /* Problem */
    FILE *wr,        /* Plot file, or NULL */
    Time dtSkip,     /* Time step for integration. */
    bool_t ms,       /* TRUE for measured errors, FALSE for estimated ones. */
    double *ek       /* Error grows proportionally to {dt^ek}. */
  );
void WriteTestResultsLine(FILE *wr, char *name, TestResults *tr);
void WriteTestResultsHeader(FILE *wr);
void WriteAdaptiveTestPlotLine(FILE *wr, Time t, State *s, Velocity *v);
void WriteAdaptiveTestPlotHeader
  ( FILE *wr, 
    char *iName, 
    Time dt, 
    Time dtMin,
    Time dtMax, 
    Dist tol
  );
void WriteStepTestPlotLine
  ( FILE *wr, 
    double scale,
    TestResults *tr,
    TestResults *trRef,
    int nt  
  );
void WriteStepTestPlotHeader(FILE *wr, char *iName, intg_prob_vec_t *p);
void WriteErrorTestPoint(FILE *wr, double x, double y);
void WriteErrorTestPlotHeader(FILE *wr, char *iName, intg_prob_vec_t *p, bool_t ms);
void WriteErrorTestIdealPoints(FILE *wr, double ek);
void StatsClear(Stats *s);
  /* Initializes the least-squares statistics. */

void StatsGather(Stats *s, double x, double y);
  /* Accumulates statistics for least-squares linear fit {y == a*x + b}. */

void StatsFit(Stats *s, double *a, double *b);
  /* Computes coeffs {a,b} of least-squares linear fit {y == a*x + b}. */

void Trace(char c);
void NL(FILE *wr);
void PVec(FILE *wr, char *fmt, double_vec_t *vt);

/* IMPLEMENTATIONS */

int main(int argn, char **argc)
  { 
    char *outDir = "out";
    intg_prob_vec_t p = intg_prob_t_Sample();
    int ng = 3;
    Grater g[ng];
    g[0] = (Grater){"euler", (Intg_T *)Intg_Euler_new()};
    g[1] = (Grater){"rkfo2", (Intg_T *)Intg_RKF2_new()};
    g[2] = (Grater){"rkfo4", (Intg_T *)Intg_RKF4_new()};
    /* Plot true solutions of all problems: */
    int i, j;
    for (i = 0; i < p.ne; i++) { PlotTrueSolution(outDir, &(p.e[i]), 512); }
    /* Make plot of error variation with step size: */
    for (j = 0; j < ng; j++)
      { DoErrorTest(outDir, &(g[j]), &p);
        fprintf(stdout, "\n");
      }
    /* Test each integrator on all problems: */
    WriteTestResultsHeader(stdout);
    AllStandardTests(outDir, g, ng, &p, "1", 1.0e-1);
    AllStandardTests(outDir, g, ng, &p, "2", 1.0e-2);
    AllStandardTests(outDir, g, ng, &p, "3", 1.0e-3);
    AllStandardTests(outDir, g, ng, &p, "4", 1.0e-6);
    /* Make plot of number of steps as function of reduction factor: */
    for (j = 0; j < ng; j++)
      { DoStepTest(outDir, &(g[j]), &p, 1.0e-3);
        fprintf(stdout, "\n");
      }
    return 0;
  }

void DoErrorTest(char *outDir, Grater *g, intg_prob_vec_t *p)
  { double dtSkip = 1.0e-2;
    bool_t ms; /* TRUE for measured error, FALSE for estimated error. */
    for (ms = FALSE; ms <= TRUE; ms++)
      { char *msTag = (ms ? "merr" : "eerr") ;
        int i;
        for (i = 0; i < p->ne; i++)
          { intg_prob_t *pi = &(p->e[i]);
            char *fName = jsprintf("%s/%s-%s-%s.plot", outDir, pi->tag, g->tag, msTag);
            FILE *wr = open_write(fName, TRUE);
            double ek;
            WriteErrorTestPlotHeader(wr, g->i->descr, p, ms);
            ErrorTest(g, pi, wr, dtSkip, ms, &ek); 
            WriteErrorTestIdealPoints(wr, ek);
            fclose(wr);
            free(fName);
          }
      }
  }
  
void AllStandardTests
  ( char *outDir, 
    Grater *g, 
    int ng, 
    intg_prob_vec_t *p, 
    char *tag, 
    Dist tol
  )
  { int i, j;
    bool_t fixedStep = FALSE;
    for (i = 0; i < p->ne; i++)
      { intg_prob_t *pi = &(p->e[i]);
        /* Choose fixed or adaptive stepsize: */
        if (fixedStep) { pi->dtMin = pi->dt; pi->dtMax = pi->dt; }
        for (j = 0; j < ng; j++)
          { DoAdaptiveTest(outDir, &(g[j]), pi, tag, tol); }
        fprintf(stdout, "\n");
      }
  }

void DoAdaptiveTest(char *outDir, Grater *g, intg_prob_t *p, char *tag, Dist tol)
  { char *fName = NULL;
    char *fName = jsprintf("%s/%s-%s-%s.plot", outDir, p->tag, g->tag, tag);
    FILE *wr = open_write(fName, TRUE);
    TestResults trFree =  AdaptiveTest(outDir, g, p, wr,   tol, FALSE, 1.0);
    TestResults trBound = AdaptiveTest(outDir, g, p, NULL, tol, TRUE,  1.0);
    trFree.mc = trBound.mc;
    trFree.mk = trBound.mk;
    WriteTestResultsLine(stdout, fName, &trFree);
    fclose(wr);
    free(fName);
  }
 
void DoStepTest(char *outDir, Grater *g, intg_prob_vec_t *p, Dist tol)
  { int NScales = 11;
    double MinScale = 0.75;
    double MaxScale = 1.50;
    char *fName = jsprintf("%s/%s-step.plot", outDir, g->tag);
    FILE *wr = open_write(fName, TRUE);
    TestResults *tr = (TestResults *)malloc(p->ne * sizeof(TestResults));
    TestResults *trRef = (TestResults *)malloc(p->ne * sizeof(TestResults));
    WriteStepTestPlotHeader(wr, g->i->descr, p);
    int k, i;
    for (k = 0; k < NScales; k++)
      { double ex = (NScales == 1 ? 0.5 : ((double)k)/((double)NScales-1));
        double scale = MinScale * exp(log(MaxScale/MinScale)*ex);
        for (i = 0; i < p->ne; i++)
          { intg_prob_t *pi = &(p->e[i]);
            tr[i] = AdaptiveTest(outDir, g, pi, NULL, tol, TRUE, scale);
            if (k == 0){trRef[i] = tr[i]; }
          }
        WriteStepTestPlotLine(wr, scale, tr, trRef, p->ne);
      }        
    fclose(wr);
    free(fName);
  }
  
void PlotTrueSolution(char *outDir, intg_prob_t *p, int nSteps)
  { 
    char *fName = jsprintf("%s/%s-true.plot", outDir, p->tag);
    FILE *wr = open_write(fName, TRUE);
    double dt = (p->t1 - p->t0)/((double)nSteps);
    WriteAdaptiveTestPlotHeader(wr, "true solution", p->dt, p->dt, p->dt, 0.0);
    /* Plot true solution: */
    State s = double_vec_new(p->n);
    Velocity v = double_vec_new(p->n);
    int i;
    for (i = 0; i <= nSteps; i++)
      { Time t = p->t0 + dt * ((double)i);
        p->sol(t, &s);
        p->rhs(t, &s, &v);
        WriteAdaptiveTestPlotLine(wr, t, &s, &v);
      }
    fclose(wr);
    free(fName);
  }

#define MAXEVALS (50000)

TestResults AdaptiveTest
  ( char *outDir, 
    Grater *g,       /* Integrator */
    intg_prob_t *p,  /* Problem */
    FILE *wr,        /* Plot file, or NULL */
    Dist tol,        /* Error tolerance */
    bool_t force,    /* Forces path to follow the true solution */
    double rStep     /* Extra stepsize factor after {g.i.adjustStepSize} */
  )
  { int n = p->n; /* Dimension of state vector. */
    
    /* {CallRHS} just calls {p->rhs} and counts the calls: */
    int nEvals = 0;
    
    auto bool_t CallRhs(Time t, State *s, Velocity *v);
    bool_t CallRhs(Time t, State *s, Velocity *v)
      { nEvals++;
        if (nEvals > MAXEVALS) { return TRUE; }
        return p->rhs(t, s, v);
      }

    /* Error statistics: */
    Stats eStats; StatsClear(&eStats);  /* Estimated error × dt statistics. */
    Stats mStats; StatsClear(&mStats);  /* Measured error × dt statistics. */
  
    auto void GatherErrorStats(Stats *s, double dt, double error);
    void GatherErrorStats(Stats *s, double dt, double error)
      { if ((error > 1.0e-8) && (dt > 1.0e-5))
          { double x = 0.5 * log(dt*dt + 1.0e-200);
            double y = 0.5 * log(error*error + 1.0e-200);
            StatsGather(s, x, y);
          }
      }

    bool_t trace = FALSE;
    bool_t debug = TRUE;
    
    int nRedos = 0;  /* Steps redone because of large error. */
    int nSteps = 0;  /* Steps actually taken. */
    
    /* Work areas for integration: */
    double_vec_t s0 = double_vec_new(n);
    double_vec_t s1 = double_vec_new(n);
    double_vec_t v0 = double_vec_new(n);
    double_vec_t er = double_vec_new(n);
    
    fprintf(stderr, "=== AdaptiveTest ========================================\n"); 
    fprintf(stderr, "integrator = %s\n", g->tag); 
    fprintf(stderr, "problem =    %s\n", p->tag); 
    fprintf(stderr, "time interval =   [ %12.8f _ %12.8f ]\n", p->t0, p->t1); 
    fprintf(stderr, "initial time step = %12.4e\n", p->dt); 
    fprintf(stderr, "time step range = [ %12.4e _ %12.4e ]\n", p->dtMin, p->dtMax); 
    fprintf(stderr, "error tolerance = %12.4e\n", tol); 
    fprintf(stderr, "\n");
    
    /* Integration headers: */
    if (wr!=NULL)
      { WriteAdaptiveTestPlotHeader(wr, g->i->descr, p->dt, p->dtMin, p->dtMax, tol); }
    
    /* Initial state and velocity: */
    Time t0 = p->t0;
    p->sol(t0, &s0);
    
    /* Adaptive-step integration: */
    Time dt = p->dt;
    Time tPrev, dtNext, t1;
    while (TRUE)
      { /* Compute the velocity at {t0,s0}: */
        CallRhs(t0, &s0, &v0);
        if (wr!=NULL) { WriteAdaptiveTestPlotLine(wr, t0, &s0, &v0); }
        /* Are we done? */
        if (t0 >= p->t1) { break; }
        if (force && (t0 > p->t0))
          { /* Client wants every step to start from the true curve. */
            /* Save the current state temporarily in {v0}: */
            int i;
            for (i = 0; i < n; i++) { v0.e[i] = s0.e[i]; }
            /* Force the state {s0} to be on the true solution: */
            p->sol(t0, &s0);
            /* Compute the error w.r.t. the true solution: */
            double dt1 = t0 - tPrev;
            double error = rn_dist(n, v0.e, s0.e);
            /* Gather that error in the accumulators {mStats}: */
            GatherErrorStats(&mStats, dt1, error);
            /* Recompute the velocity {v0} for the forced state {s0}: */
            p->rhs(t0, &s0, &v0);  /* Don't count this call in {nEvals}. */
          }
        /* Choose the end of this step: */
        tPrev = t0;
        t1 = t0 + dt; if (t1 > p->t1) { t1 = p->t1; }
        /* Take the step from {t0,s0,v0} to {t1,s1}, adaptively: */
        while (TRUE)
          { /* Try taking the step from {t0,s0,v0} to {t1,s1}: */
            if (debug)
              { fprintf(stderr, "t0 = %14.8f", t0);
                fprintf(stderr, "  dt = %14.4e", dt);
                fprintf(stderr, "  t1 = %14.8f", t1);
                fprintf(stderr, "  nE = %14d", nEvals);
                NL(stderr);
                fprintf(stderr, "  s0 = "); PVec(stderr, "%14.8f", &s0); NL(stderr);
                fprintf(stderr, "  v0 = "); PVec(stderr, "%14.8f", &v0); NL(stderr);
              }
            bool_t failed = g->i->step(g->i, &CallRhs, t0,&s0,&v0, t1,&s1, &er)
            if (failed) 
              { if (debug) { fprintf(stderr, "step failed", t0); }
                break;
              }
           
            /* Step succeeded: */
            if (debug)
              { fprintf(stderr, "  s1 = "); PVec(stderr, "%14.8f", &s1); NL(stderr);
                fprintf(stderr, "  er = "); PVec(stderr, "%14.8f", &er); NL(stderr);
              }
        
            /* Check whether the error is OK: */
            double error = rn_norm(n, er.e);
            if (debug)
              { fprintf(stderr, "  error = %14.4e", error);
                fprintf(stderr, "  tol = %14.4e", tol);
              }
            if ((error < tol) || (nEvals >= MAXEVALS))
              { /* Error is OK, accept it: */
                NL(stderr);
                GatherErrorStats(&eStats, dt, error);
                break;
              }
            /* Ops, last step was too big; reduce the step: */
            if (trace && (nSteps < 200)){ Trace('*'); }
            nRedos++;
            /* Tell {g->i->adjust} that the step was {rStep*dt}. */
            /* Usually {rStep} is 1 except when testing {adjust}. */
            dtNext = g->i->adjust(g, rStep * dt, error, tol, p->dtMin, p->dtMax);
            if (debug)
              { fprintf(stderr, "  dtNExt = %14.4e", dtNext);
                fprintf(stderr, "  dtMin = %14.4e", p->dtMin);
                fprintf(stderr, "  dtMax = %14.4e", p->dtMin);
              }
            /* If we must change {dt}, make sure that it is a significant change: */
            if (dtNext < dt
            
            /* Check whether we are stuck at an endpoint of the range: */
            if (
              { 
            bool_t stuckAtMin = (dt <= p->dtMin) && (dtNext <= p->dtMin);
            if (stuckAtMin) 
              { if (debug) { fprintf(stderr, " stuck at dtMin"); } stuck = TRUE; }
            bool_t stuckAtMax = (dt >= p->dtMax) && (dtNext >= p->dtMax);
            if (stuckAtMax) 
              { if (debug) { fprintf(stderr, " stuck at dtMax"); } stuck = TRUE; }
            if (debug) { NL(stderr); }
            dt = dtNext;
            t1 = t0 + dt;
          }
        /* The step was finally taken, phew. */
        nSteps++;
        if (debug) 
          { fprintf(stderr, "  s1 = "); PVec(stderr, "%14.8f", &s1); NL(stderr);
            NL(stderr);
          }
        if (trace && (nSteps < 200)){ Trace('-'); }
        if (nEvals >= MAXEVALS)
          { fprintf(stderr, "*** {maxEvals} exceeded, integration aborted ***\n");
            break;
          }
        /* Prepare for the next step: */
        t0 = t1;
        int i;
        for (i = 0; i < n; i++) { s0.e[i] = s1.e[i]; }
      }

    if (trace){ fprintf(stderr, "\n"); }
    if (wr!=NULL){ fflush(wr); }
    double ek, ec, mk, mc = 0.0;
    StatsFit(&eStats, &ek, &ec); ec = exp(ec);
    if (force){ StatsFit(&mStats, &mk, &mc); mc = exp(mc); }
    return 
      (TestResults)
        { /*nSteps*/ nSteps, /*nRedos*/ nRedos, /*nEvals*/ nEvals,
          /*ec*/ ec, /*ek*/ ek,
          /*mc*/ mc, /*mk*/ mk,
          /*tf*/ t1
        };
    fprintf(stderr, "=========================================================\n"); 
    fprintf(stderr, "\n"); 
  }

void ErrorTest
  ( Grater *g,       /* Integrator */
    intg_prob_t *p,  /* Problem */
    FILE *wr,        /* Plot file, or NULL */
    Time dtSkip,     /* Time step for integration. */
    bool_t ms,       /* TRUE for measured errors, FALSE for estimated ones. */
    double *ek       /* Error grows proportionally to {dt^ek}. */
  )
  { int NSizes = 11;
    double dtMin = 1.0e-6;
    double dtMax = 1.0e-1;
    double dtRef = 1.0e-4;
    
    int n = p->n; /* Dimension of state vector. */
    
    double sumSlope = 0.0;
    int nSlope = 0;
    
    Stats stats; StatsClear(&stats);  /* Error × dt statistics. */
  
    auto void GatherErrorStats(Stats *s, double dt, double error, double errorRef);
    void GatherErrorStats(Stats *s, double dt, double error, double errorRef)
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
        int k;
        for (k = 0; k < NSizes; k++)
          { double ex = (NSizes == 1 ? 0.5 : ((double)k)/((double)NSizes-1));
            double dt = dtMin * exp(log(dtMax/dtMin)*ex);
            if (ta+dt <= p->t1)
              { double error = GetError(ta,&sa,&va, dt);
                GatherErrorStats(&stats, dt, error, errorRef);
              }
          }
        
        /* Prepare for next step: */
        ta += dtSkip;
      }

    (*ek) = sumSlope/((double)nSlope);
  }

void WriteTestResultsLine(FILE *wr, char *name, TestResults *tr)
  { fprintf(wr, "%-15s", name);
    fprintf(wr, " %7d", tr->nSteps);
    fprintf(wr, " %7d", tr->nRedos);
    fprintf(wr, " %7d", tr->nEvals);
    fprintf(wr, "  %8.4f * dt ^%5.2f", tr->ec, tr->ek);
    fprintf(wr, " ");
    if (tr->mc == 0.0) 
      { fprintf(wr, "%20s", ""); }
    else
      { fprintf(wr, "%8.4f * dt ^%5.2f", tr->mc, tr->mk); }
    fprintf(wr, "\n");
    fflush(wr);
  }

void WriteTestResultsHeader(FILE *wr)
  { fprintf(wr, "%-15s", "TEST");
    fprintf(wr, "  %7s", "nSteps");
    fprintf(wr, "  %7s", "nRedos");
    fprintf(wr, "  %7s", "nEvals");
    fprintf(wr, "   %20s", "estimated error");
    fprintf(wr, "  %20s", "actual error");
    fprintf(wr, "\n") ;
  }

void WriteAdaptiveTestPlotLine(FILE *wr, Time t, State *s, Velocity *v)
  { int i;
    fprintf(wr, "  %12.4f ", t);
    for (i = 0; i < s->ne; i++) { fprintf(wr, " %12.3e", s->e[i]); }
    for (i = 0; i < s->ne; i++) { fprintf(wr, " %12.3e", v->e[i]); }
    fprintf(wr, "\n");;
  }

void WriteAdaptiveTestPlotHeader
  ( FILE *wr, 
    char *iName, 
    Time dt, 
    Time dtMin,
    Time dtMax, 
    Dist tol
  )
  { fprintf(wr, "# integrator = %s\n", iName);
    fprintf(wr, "# dt = %12.4e dtMin = %12.4e dtMax = %12.4e\n", dt, dtMin, dtMax);
    fprintf(wr, "# tol = %12.4e\n", tol);
  }

void WriteStepTestPlotLine
  ( FILE *wr, 
    double scale,
    TestResults *tr,
    TestResults *trRef,
    int nt  
  )
  { fprintf(wr, "  %12.8f ", scale);
    int i;
    for (i = 0; i < nt; i++)
      { int n = tr[i].nEvals;
        int nRef = trRef[i].nEvals;
        double r = ((double)n)/((double)nRef);
        fprintf(wr, " %8.4f", r);
      }
    fprintf(wr, "\n");
  }

void WriteStepTestPlotHeader(FILE *wr, char *gName, intg_prob_vec_t *p)
  { fprintf(wr, "# integrator = %s\n", gName);
    fprintf(wr, "# %12s ", "scale");
    int i;
    for (i = 0; i < p->ne; i++) { fprintf(wr, " %8s", p->e[i].tag); }
    fprintf(wr, "\n");
  }

void WriteErrorTestPoint(FILE *wr, double x, double y)
  { fprintf(wr, "  %12.8f %12.8f\n", x, y); }

void WriteErrorTestPlotHeader(FILE *wr, char *gName, intg_prob_vec_t *p, bool_t ms)
  { fprintf(wr, "# integrator = %s\n", gName);
    fprintf(wr, "# problems =");
    int i;
    for (i = 0; i < p->ne; i++) { fprintf(wr, " %s", p->e[i].tag); }
    fprintf(wr, "\n");
    fprintf(wr, "# %12s %12s\n", "dt/dt0", "err/err0");
  }

void WriteErrorTestIdealPoints(FILE *wr, double ek)
  { int NPoints = 51;
    double rdtMin = 1.0e-3;
    double rdtMax = 1.0e+5;
    int k;
    for (k = 0; k <= NPoints; k++) 
      { double ex = (NPoints == 1 ? 0.5 : ((double)k)/((double)NPoints-1));
        double rdt = rdtMin * exp(log(rdtMax/rdtMin)*ex);
        double rer = exp(ek*log(rdt));
        double x = 0.5 * log(rdt*rdt + 1.0e-200)/log(10.0);
        double y = 0.5 * log(rer*rer + 1.0e-200)/log(10.0);
        WriteErrorTestPoint(wr, x, y);
      }
    fflush(wr);
  }

void StatsClear(Stats *s)
  { s->xx = s->xy = s->xu = s->yu = s->uu = 0.0; }

void StatsGather(Stats *s, double x, double y)
  { s->xx = s->xx + x*x; 
    s->xy = s->xy + x*y;
    s->xu = s->xu + x;   
    s->yu = s->yu + y;
    s->uu = s->uu + 1.0;
  }

void StatsFit(Stats *s, double *a, double *b)
  { double det = s->xx*s->uu - s->xu*s->xu;
    if (fabs(det) < 1.0e-200) { det = 1.0e-200; }
    (*a) = (s->xy*s->uu - s->xu*s->yu)/det;
    (*b) = (s->xx*s->yu - s->xu*s->xy)/det;
  }

void Trace(char c)
  { fputc(c, stderr); }

void NL(FILE *wr)
  { fputc('\n', wr); }
  
void PVec(FILE *wr, char *fmt, double_vec_t *vt)
  { fputc('(', wr);
    int i;
    for (i = 0; i < vt->ne; i++)
      { fputc(' ', wr);
        fprintf(wr, fmt, vt->e[i]);
      }
    fputc(' ', wr);
    fputc(')', wr);
  } 
