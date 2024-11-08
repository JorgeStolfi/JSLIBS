/* Last edited on 2024-11-07 15:25:31 by stolfi */
/* TO FINISH !!! */

#include <stdint.h>

void qms_project_simplex(int32_t m, int32_t n, double M[], double p[], double x[], double q[])
  /* Given an {m × n} matrix {M}, with {m <= n}, and an {n}-vector {p}, 
    finds an {m}-vector {x} such that {x[i] >= 0} for {i in 0..m-1}, 
    {SUM{x[i] : i in 0..m-1} = 1}, and {|p - xM|} is minimum.
    Also returns {x M} in {q}.
    
    In other words, finds the point {q} that is closest to {p} in the
    {(m-1)}-dimensional simplex of {R^n} whose vertices are the 
    rows of {M}, and returns in {x} the barycentric coordinates
    of {p} relative to that simplex. */
  { 
    
    /* Build matrix {H} for point location relative to the simplex: */
    double H[m

    find m[0..7] such that:

      m[0..7] >= 0;
      m H = [f.R, f.G, f.B, 1]
      minimi

       { v0.R v0.G v0.B  1 }
       { v1.R v1.G v1.B  1 }
       { v2.R v2.G v2.B  1 }
   H = { v3.R v3.G v3.B  1 }
       { n0.R n0.G n0.B  0 }
       { n1.R n1.G n1.B  0 }
       { n2.R n2.G n2.B  0 }
       { n3.R n3.G n3.B  0 }
  
    bool_t debug = TRUE;

    /* Estimate roundoff error in residual computation: */
    double bmax =  rmxn_max_abs_elem(n, 1, b);
    double tiny = 1.0e-14 * bmax * sqrt(n);  /* Est. roundoff error in residuals. */

    /* Work array: */
    int32_t n1 = n+1; /* Total columns in {A} and {b}. */
    double Ab[n*n1]; /* Active sub-matrices of {A} and {b}, side by side. */
    
    /* Work vector: */
    double y[n]; /* Candidate solution to replace {x}. */

    /* Subset of active (nonzero) variables: */
    int32_t na; /* Number of active variables. */
    int32_t sit[n]; /* Situation of variable {i}: {x[i]} is active iff {sit[i] < na}. */ 
    int32_t pos[n]; /* Invariant: {pos[sit[k]] = k = sit[pos[k]]} for all {k} in {0..n-1}. */ 
    /* The active variables are {x[pos[ia]]} for {ia} in {0..na-1}. */
    /* The inactive variables are {x[pos[ii]]} for {ii} in {na..n-1}. */

    auto void activate(int32_t i); 
      /* Add inactive variable {i} to the active set. Expects {x[i]}
        to be zero, and sets it to a tiny positive value. */  

    auto void inactivate(int32_t i); 
      /* Deletes active variable {i} from the active set. Expects {x[i]}
        to be essentially zero, and sets it to zero. */  

    auto double compute_residual(int32_t i, double *u);
      /* Computes the resudual {(A u - b)[i]}.  Assumes that inactive
        variables are all zero. */
    
    auto void solve_active(double *u);
      /* Solves the eqsystem {A u = 0}, considering the equations that
        correspond to active variables, and assuming that inactive
        variables are zero. */
    
    auto void find_first_obstacle(double *x, double *y, int32_t *job, double *sob);
      /* Assumes that {x[j]} is non-negative for all {j}. Returns in
        {*job} the index of the first coordinate that becomes zero
        when one moves along the segment from {x} to {y}. Also returns
        in {*sob} the relative distance from {x} up to that point
        (scaled so that {dist(x,y) == 1}). If no coordinate becomes
        zero along that segment (that is, if {y} is positive), returns
        {*job = -1} and {*sob = 1.0}. */
    
    auto void print_sets(FILE *wr);
      /* Writes to {wr} the current active and inactive lists. */
    
    /* Start with all variables inactive, turn them on one at a time. */
    int32_t i;
    for (i = 0; i < n; i++) { sit[i] = i; pos[i] = i; x[i] = 0.0; }
    na = 0;

    /* Iterate until conditions (1)-(3) are satisfied for column {k}: */
    i = 0; /* Next variable to check for condition (3). */
    int32_t nok = 0; /* Number of variables that checked OK for (3) since last change. */
    while (nok < n)
      { /* Loop invariant: Conditions (1) and (2) are satisfied. */
        /* Also, a variable is active iff it is nonzero. */
        /* Also, the previous {nok} variables before {x[iscan]} (mod {n}) satisfy (3). */

        if (debug) { fprintf(stderr, "  - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"); }
        if (debug) { fprintf(stderr, "  nok = %d\n", nok); }
        if (debug) { gsel_print_array(stderr, 4, "%9.5f", "current solution:",  n, 1,"x",x, ""); }
        if (debug) { print_sets(stderr); }

        /* Check variable {i} for  (3), if not satisfied do something about it: */
        if (debug) { fprintf(stderr, "  checking constraints on x[%d]\n", i); }
        /* Compute residual {(A x - b)[i]}: */
        double res = compute_residual(i, x);
        if (debug) { fprintf(stderr, "    residual = %22.14f\n", res); }
        if (sit[i] < na)
          { /* Active, it is OK for (3): */
            if (debug) { fprintf(stderr, "    variable %d is active.\n", i); }
            assert(fabs(res) <= tiny);
            if (debug) { fprintf(stderr, "    variable %d is OK.\n", i); }
            nok++;
          }
        else
          { if (debug) { fprintf(stderr, "    variable %d is inactive.\n", i); }
            assert(x[i] == 0.0);
            if (res >= -tiny)
              { /* Residual for inactive variable {i} is OK, passed (3): */
                if (debug) { fprintf(stderr, "    residual is essentially positive.\n"); }
                if (debug) { fprintf(stderr, "    variable %d is OK.\n", i); }
                nok++;
              }
            else
              { /* Residual for inactive variable {i} has wrong sign, fails (3): */
                if (debug) { fprintf(stderr, "    residual is negative, possibly not OK.\n"); }
                /* Activate it: */
                activate(i);
                /* Make condition (2) true again: */
                /* Solve system with new set of variables: */
                if (debug) { print_sets(stderr); }
                solve_active(y);
                if (y[i] <= tiny)
                  { /* Should be positive. Roundoff error, perhaps? Anyway: */
                    if (debug) { fprintf(stderr, "    residual seems to be roundoff.\n"); }
                    inactivate(i);
                    /* Consider it OK, keep searching... */
                    if (debug) { fprintf(stderr, "    let's pretend that variable %d is OK.\n", i); }
                    nok++;
                  }
                else
                  { /* Variable {x[i]} definitely fails (3). */
                    if (debug) { fprintf(stderr, "    variable %d is NOT OK.\n", i); }
                    bool_t ok2 = FALSE; /* TRUE when condition (2) is OK. */
                    while (! ok2) 
                      { /* Find first variable in path {x} to {y} to become inactive, if any: */
                        int32_t jmin; double sjmin; 
                        find_first_obstacle(x, y, &jmin, &sjmin);
                        /* Advance from {x} towards {y} by the ratio {sjmin}: */
                        int32_t j;
                        for (j = 0; j < n; j++)
                          { if (sit[j] < na) { x[j] = (1-sjmin)*x[j] + sjmin*y[j]; } }

                        /* Inactivate variables that have become zero: */
                        for (j = 0; j < n; j++)
                          { if (j == i)
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
        i = (i + 1) % n;
      }

    void activate(int32_t i) 
      { if (debug) { fprintf(stderr, "  activating %d\n", i); }
        assert(x[i] == 0);
        int32_t ia = sit[i];  /* Index of active list slot that contains {i}. */
        assert(ia >= na);
        /* Make sure that {i} is in the *first* inactive list slot: */
        int32_t ia1 = na; /* Index of first inactive list slot. */
        if (ia != ia1)
          { /* Swap {pos[ia]} with {pos[ia1]}: */
            int32_t i1 = pos[ia1];
            pos[ia1] = i;
            pos[ia] = i1;
            /* Fix {iax[i]}, {iax[j]}: */
            sit[i1] = ia;
            sit[i] = ia1;
          }
        /* Drop the first slot from the inactive list: */
        na++;
        /* Since {x[i]} is now active, why not make that quite clear: */
        x[i] = tiny;
      }

    void inactivate(int32_t i)
      { if (debug) { fprintf(stderr, "  deactivating %d\n", i); }
        assert(x[i] >= -tiny);     /* Condition (1) (modulo roundoff). */
        assert(x[i] <= tiny);  /* We can inactivate only if {x} is on constraint. */
        int32_t ia = sit[i];  /* Index of active list slot that contains {i}. */
        assert(ia < na);
        /* Make sure that {i} is in the *last* active list slot: */
        int32_t ia1 = na-1; /* Index of last active list slot. */
        if (ia != ia1)
          { /* Swap {pos[ia]} with {pos[ia1]}: */
            int32_t i1 = pos[ia1];
            pos[ia1] = i;
            pos[ia] = i1;
            /* Fix {iax[i]}, {iax[j]}: */
            sit[i1] = ia;
            sit[i] = ia1;
          }
        /* Drop the last slot from active list: */
        na--;
        /* Since {x[i]} is now inactive: */
        x[i] = 0.0;
      }

    void print_sets(FILE *wr)
      { int32_t ia;
        char *sep; 
        fprintf(wr, "  active = (");
        sep = ""; 
        for (ia = 0; ia < na; ia++) 
          { fprintf(wr, "%s%d", sep, pos[ia]); sep = ","; 
            assert(pos[ia] >= 0); 
            assert(pos[ia] < n); 
            assert(sit[pos[ia]] == ia); 
          }
        fprintf(wr, ")");
        fprintf(wr, "  inactive = (");
        sep = ""; 
        for (ia = na; ia < n; ia++) 
          { fprintf(wr, "%s%d", sep, pos[ia]); sep = ","; 
            assert(pos[ia] >= 0); 
            assert(pos[ia] < n); 
            assert(sit[pos[ia]] == ia); 
          }
        fprintf(wr, ")");
        fprintf(wr, "\n");
      }
    
    double compute_residual(int32_t i, double *u)
      { double sum = 0.0;
        int32_t ja;
        /* We need to add only the nonzero variables: */
        for (ja = 0; ja < na; ja++) 
          { int32_t j = pos[ja]; sum += A[i*n + j]*u[j]; } 
        return sum - b[i];
      }
          
    void solve_active(double *u)
      {
        /* Extract active subsystem of {A u = b}, store in {Ab}: */
        int32_t ia, ja;
        for (ia = 0; ia < na; ia++)
          { int32_t i = pos[ia];
            for (ja = 0; ja < na; ja++)
              { int32_t j = pos[ja]; Ab[ia*(na+1) + ja] = A[i*n + j]; }
            Ab[ia*(na+1) + na] = b[i];
          }
        /* Solve subsystem: */
        double ua[na];
        if (debug) { gsel_print_array(stderr, 4, "%9.5f", "subsystem:",  na, na+1,"Ab",Ab, ""); }
        gsel_triangularize(na, na+1, Ab, TRUE, 0.0);
        gsel_diagonalize(na, na+1, Ab);
        gsel_normalize(na, na+1, Ab);
        int32_t rank_ext = gsel_extract_solution(na, na+1, Ab, 1, ua);
        if (debug) { gsel_print_array(stderr, 4, "%9.5f", "subsystem's solution:",  na, 1,"ua",ua, ""); }
        arrsert(rank_ext <= na);
        /* Unpack {ua} into {u}: */
        int32_t i;
        for (i = 0; i < n; i++)
          { int32_t ia = sit[i]; u[i] = (ia < na ? ua[ia] : 0.0); }
        if (debug) { gsel_print_array(stderr, 4, "%9.5f", "new solution:",  n, 1,"u",u, ""); }
      }
    
    void find_first_obstacle(double *u, double *v, int32_t *job, double *sob)
      { (*sob) = +INFINITY; 
        (*job) = -1;
        int32_t ja;
        for (ja = 0; ja < na; ja++)
          { int32_t j = pos[ja]; 
            assert(u[j] > 0.0);
            if (v[j] <= 0.0)
              { double sj = u[j]/(u[j]-v[j]);
                if (sj > 1.0) { sj = 1.0; }
                if (sj < (*sob)) { (*sob) = sj; (*job) = j; }
              }
          }
        if ((*job) < 0) { (*sob) = 1.0; }
        if (debug) { fprintf(stderr, "  pivoting to %d at s = %10.7f\n", (*job), (*sob)); }
      }
  }


  {

    bool_t s[4]; /* {s[i] = TRUE} means vertex {i} belongs to the face containing  the projection. */
    /* Assumes {m} is the barycentric coords of {f} relative to {v}. */

    "let {s[0..3] = FALSE}, {n[0..3] = UNDEF}."

    "loop"

      { "compute mixing coefs {m} of {u[0..3]} to give {f}, where {u[i] = (s[i] ? n[i] : v[i])}"

        if ("{m} is all >= 0") { return; }

        "let {i} be such that {m[i]} is the most negative coord of {m}";

        "set {s[i] = TRUE}, {n[i] =} unit vector normal to face opposite to {v[i]}";

      }

