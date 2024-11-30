/* See qmin_simplex.h */
/* Last edited on 2024-11-29 22:25:11 by stolfi */

#include <stdint.h>
#include <math.h>
#include <limits.h>

#include <rmxn.h>
#include <assert.h>
#include <affirm.h>
#include <bool.h>
#include <sign.h>

#include <gauss_elim.h>
#include <gauss_elim_solve.h>

#include <qmin_simplex.h>

#define ER(...) fprintf(stderr, ##__VA_ARGS__)

void qmin_simplex(uint32_t n, double A[], double b[], double x[])
  { bool_t debug = FALSE;

    /* Estimate roundoff error in residual computation: */
    double bmax =  rmxn_max_abs_elem(n, 1, b);
    double tiny = 1.0e-14 * bmax * sqrt(n);  /* Est. roundoff error in residuals. */

    /* Work array: */
    uint32_t n1 = n+1; /* Total columns in {A} and {b}. */
    double Ab[n*n1]; /* Active sub-matrices of {A} and {b}, side by side. */
    
    /* Work vector: */
    double y[n]; /* Candidate solution to replace {x}. */
    
    /* Total number of variables, active or inactive: */
    uint32_t nx = n;

    /* Subset of active (nonzero) variables from {0..nx-1}: */
    uint32_t na; /* Number of active variables. */
    uint32_t ia_from_ix[nx]; /* Situation of variable {i}: {x[i]} is active iff {ia_from_ix[i] < na}. */ 
    uint32_t ix_from_ia[nx]; /* Invariant: {ix_from_ia[ia_from_ix[k]] = k = ia_from_ix[ix_from_ia[k]]} for all {k} in {0..n-1}. */ 
    /* The active variables are {x[ix_from_ia[i]]} for {ia} in {0..na-1}. */
    /* The inactive variables are {x[ix_from_ia[i]]} for {i} in {na..n-1}. */
    
    /* The highest set of active variables seen in lex order, to detect loops: */
    bool_t hiset[nx]; /* Bit {hiset[ix]} is true if {x[ix]} was part of the high set. */
    for (uint32_t ix = 0;  ix < n; ix++) { hiset[ix] = FALSE; }
    
    /* Start with all variables inactive, turn them on one at a time. */
    for (uint32_t ix = 0;  ix < n; ix++) { ia_from_ix[ix] = (uint32_t)ix; ix_from_ia[ix] = (uint32_t)ix; x[ix] = 0.0; }
    na = 0;
    
    auto void move_towards(double y[]);
      /* Called when some variable {x[ix_act]} was tentatively changed from 
        inactive to active, and, as a result, the new solution {y[0..nx-1]} of the
        reduced system yielded {y[ix_cur]} positive. 
        
        Tries to set {y} as the new current point {x}, but, if that
        would make some active variable {x[ix_hit]} other than
        {x[ix_act]} zero or negative, moves only part of the way, until
        the first of those variables hits 0. Adjusts the active set if
        necessary.
        
        Specifically, if all active variables in {y} are positive, just
        sets {x} to {y} and leaves the active set unchanged.
        
        Otherwise, computes {z} along the line from {x} to {y} when the
        first active variable {z[ix_hit]} would become zero, and sets
        {x} to {z}. See {find_first_obstacle}. Then inactivates all
        active variables whose {x} become essentially 0. At least one of
        them will remain active. */

    auto void activate(uint32_t ix); 
      /* Add inactive variable {x[ix]} to the active set. Expects {x[ix]}
        to be zero, and sets it to a tiny positive value. */  

    auto void inactivate(uint32_t ix); 
      /* Deletes active variable {ix} from the active set. Expects {x[ix]}
        to be essentially zero, and sets it to zero. */  

    auto double gradient_component(uint32_t ix);
      /* Computes component {ix} of the gradient of {Q} at {x}, namely
        the {dQ/dx[ix] = (A x - b)[ix]}, the residual of equation {ix}
        in the system that defines the optimum of {Q}.
        
        Assumes that all inactive variables are zero; that is,
        {x[ix_from_ia[ia]] = 0} for {ia} in {na..nx-1}. */
    
    auto void solve_active(double y[]);
      /* Solves the equation system {U u = v}, where {u[0..na-1]} is the
        values of the active variables, and {U} is those rows of and
        columns of {A} that refers to those variables, assuming that
        inactive variables are zero. Then returns the full solution
        vector {y[0..nx-1]} with active variables copied from {u}, and
        inactive ones set to zero. That is, sets {y[ix_from_ia[ia]] =
        u[ia]} for {ia} in {0..na-1}, and {y[ix_from_ia[ia]] = 0} for
        {ia} in {na..nx-1}. */

    auto void find_first_obstacle(double y[], int32_t *ix_hit_P, double *frac_hit_P);
      /* Assuming that the current solution {x[0..nx-1} satisfies all constraints,
        and {y[0..nx-1]} is a new candidate solution which may violate some of them, 
        finds the first active variable, if any, that violates some constraint
        as one moves from {x} to {y} in a straight line.
        
        Specifically, let {P(s)} be the point of {\RR^na} that is {s} of
        the way between {x} and {y}; that is, {P(s)[ia] = (1-s)*x[ia] +
        s*y[ia]}. The procedure finds the index {ix_hit} in {0..nx-1}
        such that {P(s)[ix_hit]} is the first coordinate of {P(s)} to
        become zero when {s} varies from 0 to 1. It then sets {*ix_hit_P}
        to that index, and {*frac_hit_P} the value of {s} when
        {P(s)[ix_hit]} becomes zero. If no coordinate becomes zero along
        that segment (that is, if {v[ia]} is positive for all {ia}),
        returns {*ix_hit_P = -1} and {*frac_hit_P = 1.0}.
        
        The procedure expects that, for every {ix} in the inactive set,
        both {x[ix]} and {y[ix]} are exactly zero, so {P(s)} will be
        zero too. For every {ix} in the active set, the procedure
        requires that {x[ix]} be definitely positive. 
        
        The procedure also requires that {y[ix_pos]} be definitely
        positive for at least one {ix_pos} in the active set. Thus
        {P(s)[ix_pos]} will be definitely positive for any {s} in
        {[0_1]}, the returned {ix_hit} will not be {ix_pos}, and at least
        {x[ix_pos]} will remain active. */
        
    auto bool_t detect_loop(void);
      /* Check the lex order of the current active variable set against {hiset}.
        If the current set is less than {hiset}, return false.
        If the current set is equal to {hiset}, return true.
        If the current set is greater than {hiset}, set it as the new {hiset},
        and return false.. */
    
    auto void print_sets(void);
      /* Writes to {stderr} the current active and inactive lists. */

    /* Iterate until conditions (1)-(3) are satisfied for every variable: */
    uint32_t ix_cur = 0; /* Next variable to check for condition (3). */
    uint32_t nok = 0; /* Count of variables that checked OK for (3) since last change. */
    uint32_t max_iters = nx*(nx+5); /* Let's hope. But the worst case may be exponential, it seems. */
    uint32_t niter = 0; /* Number of main iterations done. */
    bool_t looping = FALSE;
    while ((nok < nx) && (niter < max_iters) && (! looping))
      { /* Loop invariant: Conditions (1) and (2) are satisfied. */
        /* Also, a variable is active iff it is nonzero. */
        /* Also, the previous {nok} variables before {x[ix_cur]} (mod {n}) satisfy (3). */

        if (debug) { ER("  - - - - - - - - - - - - - - - - - - - - - - - - - - -\n"); }
        if (debug) { ER("  nok = %d\n", nok); }
        if (debug) { gauss_elim_print_array(stderr, 4, "%9.5f", "current solution:",  nx, 1,"x",x, ""); }
        if (debug) { ER("  current active set:\n"); }
        if (debug) { print_sets(); }
        if (debug) { ER("  checking variable %d.\n", ix_cur); }
      
        if (ia_from_ix[ix_cur] < na)
          { /* Active, it is OK for (3): */
            if (debug) { ER("    variable %d is active, retained active.\n", ix_cur); }
            nok++;
          }
        else
          { /* Inactive. Check (3): */
            if (debug) { ER("    variable %d is inactive.\n", ix_cur); }
            assert(x[ix_cur] == 0.0);
            /* Check condition (3) for variable {x[ix_cur]}: */
            if (debug) { ER("    checking constraint {dQ/dx[%d] > 0}\n", ix_cur); }
            /* Compute component {ix_cur} of the gradient: */
            double grad_ix_cur = gradient_component(ix_cur);
            if (debug) { ER("    gradient {dQ/dx[%d] = (A x - b)[%d]} = %22.14f\n", ix_cur, ix_cur, grad_ix_cur); }
            if (grad_ix_cur >= -tiny)
              { if (debug) { ER("    negative of gradient {-dQ/dx[%d]} pushes into constraint {x[%d] >= 0}.\n", ix_cur, ix_cur); }
                if (debug) { ER("    variable %d is OK as inactuve.\n", ix_cur); }
                nok++;
              }
             else
              { if (debug) { ER("    negative gradient {-dQ/dx[%d]} pushes away from constraint {x[%d] >= 0}.\n", ix_cur, ix_cur); }
                if (debug) { ER("    activating variable %d.\n", ix_cur); }
                activate(ix_cur);
                /* Try to make condition (2) true again. */
                /* Solve system for new set of active variables: */
                if (debug) { print_sets(); }
                solve_active(y);
                if (y[ix_cur] <= tiny)
                  { /* Should be positive. Roundoff error, perhaps? Anyway: */
                    if (debug) { ER("    solving for active vars yields {x[%d]} =  roundoff (%24.16e).\n", ix_cur, y[ix_cur]); }
                    if (debug) { ER("    iactivating {x[%d]} again and leaving it zero.\n", ix_cur); }
                    inactivate(ix_cur);
                    /* Consider it OK, keep searching... */
                    if (debug) { ER("    let's pretend that variable %d is OK.\n", ix_cur); }
                    nok++;
                  }
                else
                  { if (debug) { ER("    solving for active variables yields {y[%d] = %24.16e}.\n", ix_cur, y[ix_cur]); }
                    if (debug) { ER("    so variable %d must defintely remain active.\n", ix_cur); }
                    if (debug) { ER("    move to {y} unless hits another constraint first.\n"); }
                    move_towards(y);
                    if (debug) { ER("    new active and inactive sets:\n"); }
                    if (debug) { print_sets(); }
                    /* Since the active set and/or current solution {x} changed, we must check all variables again: */
                    nok = 0;
                    /* Count as another main iteration, although incomplete: */
                    niter++;
                  }
              }
          }
          
        /* Check current active variable set against {hiset}, update as needed: */
        looping = detect_loop();
        
        /* Go to next variable: */
        ix_cur = (ix_cur + 1) % nx;
      }
    if (nok >= nx)
      { if (debug) { ER("  seems to have converged in %d iterations...\n", niter); } }
    else if (niter >= max_iters)
      { if (debug) { ER("  too many iterations (%d), gave up...\n", niter); } }
    else if (looping)
      { if (debug) { ER("  looping detected during iteration (%d)...\n", niter); } }
    else
      { assert(FALSE); }
    
    return;
      
    void move_towards(double y[])
      { /* Find first variable in path {x} to {y} to become inactive, if any: */
        int32_t ix_hit; double frac_hit; 
        find_first_obstacle(y, &ix_hit, &frac_hit);
        if (debug) { ER("    moving {x} to {%24.16e} of the way towards {y}\n", frac_hit); }
        assert((frac_hit > 0.0) && (frac_hit <= 1.0));
        if (debug) { 
          { if (ix_hit < 0) 
              { ER("    moving {x} all the way to {y}\n"); }
            else
              { ER("    moving {x} to {%24.16e} of the way towards {y}\n", frac_hit); }
          }}
             
        int32_t ix_pos = -1; /* Some variable that remained positive. */
        for (uint32_t ia = 0; ia < nx; ia++)
          { uint32_t ix = ix_from_ia[ia];
            if (ia > na)
              { /* {x[ix]} is inactive: */
                assert(x[ix] == 0); 
                assert(y[ix] == 0);
                assert(ix != ix_hit);
              }
            else
              { if (ix_hit < 0)
                  { x[ix] = y[ix]; }
                else if (ix == ix_hit)
                  { /* Ran into this constraint, force it inactive: */ 
                    inactivate(ix);
                  }
                else
                  { x[ix] = (1-frac_hit)*x[ix] + frac_hit*y[ix];
                    if (x[ix] < tiny)
                      { /* Coincidentally ran into this constraint too: */ 
                        inactivate(ix);
                      }
                  }
                if (x[ix] >= tiny) { ix_pos = (int32_t)ix; }
              }
          }
        /* At least SOME variable must have remained positive: */
        assert(ix_pos != -1);
      }
    
    void find_first_obstacle(double y[], int32_t *ix_hit_P, double *frac_hit_P)
      { double frac_hit = 1.0; 
        int32_t ix_hit = -1;
        for (uint32_t ia = 0;  ia < na; ia++)
          { int32_t ix = (int32_t)ix_from_ia[ia]; 
            assert(x[ix] > 0.0);
            if (y[ix] <= 0.0)
              { /* Could be this one: */
                double sj = x[ix]/(x[ix]-y[ix]);
                assert(sj >= 0);
                if (sj > 1.0) { sj = 1.0; }
                if (sj < frac_hit) { frac_hit = sj; ix_hit = ix; }
              }
          }
        if (debug) { ER("    hit constraint on {x[%d]} at {s = %24.16e}\n", ix_hit, frac_hit); }
        (*ix_hit_P) = ix_hit;
        (*frac_hit_P) = frac_hit;
      }

    void activate(uint32_t ix) 
      { if (debug) { ER("  activating %d\n", ix); }
        assert(x[ix] == 0);
        uint32_t ia = ia_from_ix[ix];  /* Index such that {ix_from_ia[ia] = ix}. */
        assert(ix_from_ia[ia] == ix);
        assert(ia >= na);
        /* Make sure that {ix} is in the *first* inactive list slot: */
        uint32_t ia1 = na; /* Index of first inactive list slot. */
        if (ia != ia1)
          { /* Swap {ix_from_ia[ia]} with {ix_from_ia[ia1]}: */
            uint32_t ix1 = ix_from_ia[ia1];
            ix_from_ia[ia1] = ix; ia_from_ix[ix] = ia1;
            ix_from_ia[ia] = ix1; ia_from_ix[ix1] = ia;
          }
        /* Drop the first slot from the inactive list: */
        na++;
        /* Since {x[ix]} is now active, why not make that quite clear: */
        x[ix] = tiny;
      }

    void inactivate(uint32_t ix)
      { if (debug) { ER("  deactivating %d\n", ix); }
        assert(x[ix] <= tiny);   /* We can inactivate only if {x} is on constraint. */
        assert(x[ix] >= -tiny);  /* Condition (1) (modulo roundoff). */
        uint32_t ia = ia_from_ix[ix];    /* Index of {ix} in active list. */
        assert(ia < na);
        /* Make sure that {ix} is in the *last* active list slot: */
        uint32_t ia1 = (uint32_t)(na-1); /* Index of last active list slot. */
        if (ia != ia1)
          { /* Swap {ix_from_ia[ia]} with {ix_from_ia[ia1]}: */
            uint32_t ix1 = ix_from_ia[ia1];
            ix_from_ia[ia1] = ix; ia_from_ix[ix] = ia1;
            ix_from_ia[ia] = ix1; ia_from_ix[ix1] = ia;
          }
        /* Drop the last slot from active list: */
        assert(na > 0);
        na--;
        /* Since {x[ix]} is now inactive: */
        x[ix] = 0.0;
      }

    void print_sets(void)
      { char *sep; 
        sep = ""; 
        ER("      active = ( ");
        for (uint32_t ia = 0;  ia < na; ia++) 
          { ER("%s%d", sep, ix_from_ia[ia]); sep = ","; 
            assert(ix_from_ia[ia] >= 0); 
            assert(ix_from_ia[ia] < n); 
            assert(ia_from_ix[ix_from_ia[ia]] == ia); 
          }
        ER(" ) inactive = (");
        sep = ""; 
        for (int32_t ia = (int32_t)na; ia < n; ia++) 
          { ER("%s%d", sep, ix_from_ia[ia]); sep = ","; 
            assert(ix_from_ia[ia] < n); 
            assert(ia_from_ix[ix_from_ia[ia]] == ia); 
          }
        ER(" )\n");
      }
    
    double gradient_component(uint32_t ix)
      { double sum = 0.0;
        /* We need to add only the nonzero variables: */
        for (uint32_t ia = 0;  ia < na; ia++) 
          { uint32_t jx = ix_from_ia[ia]; sum += A[ix*n + jx]*x[jx]; } 
        return sum - b[ix];
      }
          
    void solve_active(double y[])
      {
        /* Extract active subsystem of {A u = b}, store in {Ab}: */
        for (uint32_t ia = 0;  ia < na; ia++)
          { uint32_t ix = (uint32_t)ix_from_ia[ia];
            for (uint32_t ja = 0;  ja < na; ja++)
              { uint32_t j = ix_from_ia[ja]; 
                Ab[ia*(na+1) + ja] = A[ix*n + j];
              }
            Ab[ia*(na+1) + na] = b[ix];
          }
        /* Solve subsystem: */
        double ua[na];
        if (debug) { gauss_elim_print_array(stderr, 4, "%9.5f", "subsystem:",  na, na+1,"Ab",Ab, ""); }
        uint32_t rank_ext;
        gauss_elim_solve_packed(na, na, 1, Ab, ua, 0.0, &rank_ext, NULL);
        assert(rank_ext <= na);
        if (debug) { gauss_elim_print_array(stderr, 4, "%9.5f", "subsystem's solution:",  na, 1,"ua",ua, ""); }
        /* Unpack {ua[0..na-1]} into {y[0..nx-1]}: */
        for (uint32_t ix = 0;  ix < nx; ix++)
          { uint32_t ia = ia_from_ix[ix]; 
            y[ix] = (ia < na ? ua[ia] : 0.0);
          }
        if (debug) { gauss_elim_print_array(stderr, 4, "%9.5f", "new solution:",  nx, 1,"y",y, ""); }
      }
       
    bool_t detect_loop(void)
      { sign_t sgn = 0; /* Sign of comparison of current set with {hiset}, equal for now. */
        for (uint32_t ix = 0; (ix < nx) && (sgn >= 0); ix++)
          { /* Is {x[ix]} active? */
            bool_t actix = (ia_from_ix[ix] < na);
            if (sgn == 0)
              { /* Current set matched {hiset[0..ix-1]} so far. Check {ix}: */
                if (hiset[ix] != actix)
                  { if (actix) 
                      { /* Current set is lex greater, note and start updating {hiset}: */
                        sgn = +1; hiset[ix] = TRUE;
                      }
                    else
                      { /* Current set is lex less than {hiset}, note but leave {histe} unchanged: */
                        sgn = -1;
                      }
                   }
                 else
                   { /* Still matching, do nothing: */ }
               }
             else if (sgn > 0)
               { /* Continue updating {hiset}: */
                 hiset[ix] = actix;
               }
             else
               { /* Already excluded in {for}: */
                 assert(FALSE);
               }
           }
        return (sgn == 0);
      }

   }


