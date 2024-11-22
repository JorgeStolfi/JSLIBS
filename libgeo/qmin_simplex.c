/* See qmin_simplex.h */
/* Last edited on 2024-11-22 02:13:45 by stolfi */

#include <stdint.h>
#include <math.h>
#include <limits.h>

#include <rmxn.h>
#include <assert.h>
#include <affirm.h>
#include <bool.h>
#include <gauss_elim.h>

#include <qmin_simplex.h>

void qms_quadratic_min(uint32_t n, double A[], double b[], double x[])
  { bool_t debug = FALSE;

    /* Estimate roundoff error in residual computation: */
    double bmax =  rmxn_max_abs_elem(n, 1, b);
    double tiny = 1.0e-14 * bmax * sqrt(n);  /* Est. roundoff error in residuals. */

    /* Work array: */
    int32_t n1 = n+1; /* Total columns in {A} and {b}. */
    double Ab[n*n1]; /* Active sub-matrices of {A} and {b}, side by side. */
    
    /* Work vector: */
    double y[n]; /* Candidate solution to replace {x}. */

    /* Subset of active (nonzero) variables: */
    uint32_t sit[n]; /* Situation of variable {i}: {x[i]} is active iff {sit[i] < na}. */ 
    uint32_t pos[n]; /* Invariant: {pos[sit[k]] = k = sit[pos[k]]} for all {k} in {0..n-1}. */ 
    /* The active variables are {x[pos[ia]]} for {ia} in {0..na-1}. */
    /* The inactive variables are {x[pos[ii]]} for {ii} in {na..n-1}. */

    auto void activate(uint32_t i); 
      /* Add inactive variable {i} to the active set. Expects {x[i]}
        to be zero, and sets it to a tiny positive value. */  

    auto void inactivate(uint32_t i); 
      /* Deletes active variable {i} from the active set. Expects {x[i]}
        to be essentially zero, and sets it to zero. */  

    auto double compute_residual(uint32_t ieq, double *u);
      /* Computes the resudual {(A u - b)[ieq]} of equation {ieq}.
        Assumes that inactive variables are all zero. */
    
    auto void solve_active(double *u);
      /* Solves the eqsystem {A u = b}, keeping the equations that
        correspond to active variables, and assuming that inactive
        variables are zero. */
    
    auto void find_first_obstacle(double u[], double v[], int32_t *job_P, double *sob_P);
      /* Finds the first coordinate, if any, that becomes negative as
        one moves from {u} to {v}.
        
        Specifically, let {X[ia]} be {u[pos[ia]]} and {Y[ia]} be
        {v[pos[ia]]} for {ia} in {0..na-1}.  The coordinates
        {X[j]} must be strictly positive for all {j}.  Let {P(s)[0..na-1]} be the
        point that is {s} of the way between {X} and {Y}; that
        is, {P(s)[ia] = (1-s)*X[ia] + s*Y[ia]}.
        
        The procedure finds the first index {iob} in {0..na-1} such that
        {P(s)[pos[iob]]} becomes zero when {s} varies from 0 to 1. It
        then returns {*sob_P} the value of {s} when {P(s)[iob]} becomes
        zero, and in {*job_P} the index {pos[iob]}. If no coordinate
        becomes zero along that segment (that is, if {Y[ia]} is positive
        for all {ia}), returns {*job_P = INT32_MAX} and {*sob_P =
        1.0}. */
    
    auto void print_sets(FILE *wr);
      /* Writes to {wr} the current active and inactive lists. */
    
    /* Start with all variables inactive, turn them on one at a time. */
    for (int32_t i = 0; i < n; i++) { sit[i] = (uint32_t)i; pos[i] = (uint32_t)i; x[i] = 0.0; }
    int32_t na = 0; /* Number of active variables. */

    /* Iterate until conditions (1)-(3) are satisfied for column {k}: */
    int32_t jcur = 0; /* Next variable to check for condition (3). */
    int32_t nok = 0; /* Count of variables that checked OK for (3) since last change. */
    while (nok < n)
      { /* Loop invariant: Conditions (1) and (2) are satisfied. */
        /* Also, a variable is active iff it is nonzero. */
        /* Also, the previous {nok} variables before {x[iscan]} (mod {n}) satisfy (3). */

        if (debug) { fprintf(stderr, "  - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"); }
        if (debug) { fprintf(stderr, "  nok = %d\n", nok); }
        if (debug) { gsel_print_array(stderr, 4, "%9.5f", "current solution:",  n, 1,"x",x, ""); }
        if (debug) { print_sets(stderr); }

        /* Check variable {jcur} for  (3), if not satisfied do something about it: */
        if (debug) { fprintf(stderr, "  checking constraints on x[%d]\n", jcur); }
        /* Compute residual {(A x - b)[jcur]}: */
        double res = compute_residual(jcur, x);
        if (debug) { fprintf(stderr, "    residual = %22.14f\n", res); }
        if (sit[jcur] < na)
          { /* Active, it is OK for (3): */
            if (debug) { fprintf(stderr, "    variable %d is active.\n", jcur); }
            assert(fabs(res) <= tiny);
            if (debug) { fprintf(stderr, "    variable %d is OK.\n", jcur); }
            nok++;
          }
        else
          { if (debug) { fprintf(stderr, "    variable %d is inactive.\n", jcur); }
            assert(x[jcur] == 0.0);
            if (res >= -tiny)
              { /* Residual for inactive variable {jcur} is OK, passed (3): */
                if (debug) { fprintf(stderr, "    residual is essentially positive.\n"); }
                if (debug) { fprintf(stderr, "    variable %d is OK.\n", jcur); }
                nok++;
              }
            else
              { /* Residual for inactive variable {jcur} has wrong sign, fails (3): */
                if (debug) { fprintf(stderr, "    residual is negative, possibly not OK.\n"); }
                /* Activate it: */
                activate(jcur);
                /* Make condition (2) true again: */
                /* Solve system with new set of variables: */
                if (debug) { print_sets(stderr); }
                solve_active(y);
                if (y[jcur] <= tiny)
                  { /* Should be positive. Roundoff error, perhaps? Anyway: */
                    if (debug) { fprintf(stderr, "    residual seems to be roundoff.\n"); }
                    inactivate(jcur);
                    /* Consider it OK, keep searching... */
                    if (debug) { fprintf(stderr, "    let's pretend that variable %d is OK.\n", jcur); }
                    nok++;
                  }
                else
                  { /* Variable {x[jcur]} definitely fails (3). */
                    if (debug) { fprintf(stderr, "    variable %d is NOT OK.\n", jcur); }
                    bool_t ok2 = FALSE; /* TRUE when condition (2) is OK. */
                    while (! ok2) 
                      { /* Find first variable in path {x} to {y} to become inactive, if any: */
                        int32_t jmin; double sjmin; 
                        find_first_obstacle(x, y, &jmin, &sjmin);
                        assert((sjmin > 0.0) && (sjmin <= 1.0));
                        /* Advance from {x} towards {y} by the ratio {sjmin}: */
                        for (int32_t j = 0; j < n; j++)
                          { if (sit[j] < na) { x[j] = (1-sjmin)*x[j] + sjmin*y[j]; } }

                        /* Inactivate variables that have become zero: */
                        for (int32_t j = 0; j < n; j++)
                          { if (j == jcur)
                              { /* Forced active, so let's make that quite clear: */ 
                                assert(sit[j] < na);
                                if (x[j] < tiny) { x[j] = tiny; }
                              }
                            else if (sit[j] >= na) 
                              { /* Already inactive, let it stand: */
                                assert(x[j] == 0.0);
                              }
                            else if (j == jmin)
                              { /* Ran into this constraint, force it inactive: */ 
                                inactivate(j);
                              }
                            else if (x[j] < tiny)
                              { /* Coincidentally ran into this constraint too: */ 
                                inactivate(j);
                              }
                          }
                        if (sjmin == 1.0)  
                          { /* Condition (2) should be OK: */ 
                            ok2 = TRUE;
                          }
                        else
                          { /* We ran into a constraint while trying to get to {y}. */
                            /* Therefore condition (2) is not guaranteed yet. */
                            /* Solve system with new set of variables: */
                            if (debug) { print_sets(stderr); }
                            solve_active(y);
                          }
                      }
                    /* We must check condition (3) all over again. */
                    nok = 0;
                  }
              }
          }
        /* Go to next variable: */
        jcur = (jcur + 1) % n;
      }

    void activate(uint32_t jina) 
      { if (debug) { fprintf(stderr, "  activating %d\n", jina); }
        assert(x[jina] == 0);
        int32_t ii = (int32_t)sit[jina];  /* Index such that {pos[ii] = jina}. */
        assert(pos[ii] == jina);
        assert(ii >= na);
        /* Make sure that {jina} is in the *first* inactive list slot: */
        int32_t ii1 = (int32_t)na; /* Index of first inactive list slot. */
        if (ii != ii1)
          { /* Swap {pos[ii]} with {pos[ii1]}: */
            uint32_t jina1 = pos[ii1];
            pos[ii1] = jina; sit[jina] = ii1;
            pos[ii] = jina1; sit[jina1] = ii;
          }
        /* Drop the first slot from the inactive list: */
        na++;
        /* Since {x[jina]} is now active, why not make that quite clear: */
        x[jina] = tiny;
      }

    void inactivate(uint32_t jact)
      { if (debug) { fprintf(stderr, "  deactivating %d\n", jact); }
        assert(x[jact] >= -tiny);  /* Condition (1) (modulo roundoff). */
        assert(x[jact] <= tiny);   /* We can inactivate only if {x} is on constraint. */
        int32_t ia = (int32_t)sit[jact];    /* Index such that {pos[ia] = jact}. */
        assert(ia < na);
        /* Make sure that {jact} is in the *last* active list slot: */
        int32_t ia1 = (int32_t)na-1; /* Index of last active list slot. */
        if (ia != ia1)
          { /* Swap {pos[ia]} with {pos[ia1]}: */
            uint32_t jact1 = pos[ia1];
            pos[ia1] = jact; sit[jact] = ia1;
            pos[ia] = jact1; sit[jact1] = ia;
          }
        /* Drop the last slot from active list: */
        assert(na > 0);
        na--;
        /* Since {x[jact]} is now inactive: */
        x[jact] = 0.0;
      }

    void print_sets(FILE *wr)
      { char *sep; 
        fprintf(wr, "  active = (");
        sep = ""; 
        for (int32_t ia = 0; ia < na; ia++) 
          { fprintf(wr, "%s%d", sep, pos[ia]); sep = ","; 
            assert(pos[ia] >= 0); 
            assert(pos[ia] < n); 
            assert(sit[pos[ia]] == ia); 
          }
        fprintf(wr, ")");
        fprintf(wr, "  inactive = (");
        sep = ""; 
        for (int32_t ia = na; ia < n; ia++) 
          { fprintf(wr, "%s%d", sep, pos[ia]); sep = ","; 
            assert(pos[ia] >= 0); 
            assert(pos[ia] < n); 
            assert(sit[pos[ia]] == ia); 
          }
        fprintf(wr, ")");
        fprintf(wr, "\n");
      }
    
    double compute_residual(uint32_t ieq, double *u)
      { double sum = 0.0;
        /* We need to add only the nonzero variables: */
        for (int32_t ia = 0; ia < na; ia++) 
          { uint32_t j = pos[ia]; sum += A[ieq*n + j]*u[j]; } 
        return sum - b[ieq];
      }
          
    void solve_active(double *u)
      {
        /* Extract active subsystem of {A u = b}, store in {Ab}: */
        for (int32_t ia = 0; ia < na; ia++)
          { int32_t i = (int32_t)pos[ia];
            for (int32_t ja = 0; ja < na; ja++)
              { int32_t j = (int32_t)pos[ja]; Ab[ia*(na+1) + ja] = A[i*n + j]; }
            Ab[ia*(na+1) + na] = b[i];
          }
        /* Solve subsystem: */
        double ua[na];
        if (debug) { gsel_print_array(stderr, 4, "%9.5f", "subsystem:",  na, na+1,"Ab",Ab, ""); }
        gsel_triangularize(na, na+1, Ab, TRUE, 0.0);
        gsel_diagonalize(na, na+1, Ab);
        gsel_normalize(na, na+1, Ab);
        uint32_t rank_ext = gsel_extract_solution(na, na+1, Ab, 1, ua);
        assert(rank_ext <= na);
        if (debug) { gsel_print_array(stderr, 4, "%9.5f", "subsystem's solution:",  na, 1,"ua",ua, ""); }
        /* Unpack {ua} into {u}: */
        for (int32_t i = 0; i < n; i++)
          { int32_t ia =(int32_t) sit[i]; u[i] = (ia < na ? ua[ia] : 0.0); }
        if (debug) { gsel_print_array(stderr, 4, "%9.5f", "new solution:",  n, 1,"u",u, ""); }
      }
    
    void find_first_obstacle(double u[], double v[], int32_t *job_P, double *sob_P)
      { double sob = 1.0; 
        int32_t job = INT32_MAX;
        for (int32_t ia = 0; ia < na; ia++)
          { int32_t j = (int32_t)pos[ia]; 
            assert(u[j] > 0.0);
            if (v[j] <= 0.0)
              { double sj = u[j]/(u[j]-v[j]);
                if (sj > 1.0) { sj = 1.0; }
                if (sj < sob) { sob = sj; job = j; }
              }
          }
        if (debug) { fprintf(stderr, "  pivoting to %d at s = %10.7f\n", job, sob); }
        (*job_P) = job; (*sob_P) = sob;
      }
  }


