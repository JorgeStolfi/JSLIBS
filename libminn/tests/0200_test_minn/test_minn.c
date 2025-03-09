#define PROG_NAME "test_minn"
#define PROG_DESC "test of {minn.h}"
#define PROG_VERS "1.0"

/* Last edited on 2025-02-16 20:31:49 by stolfi */ 
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define tmnn_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <rn.h>
#include <rmxn.h>
#include <rmxn_throw.h>
#include <rmxn_ellipsoid.h>
#include <jsmath.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <affirm.h>

#include <minn.h>

void tmnn_do_one_test(int32_t it, char *fname, uint32_t n, bool_t verbose);
  /* Tests the functions of {minn.h}, iteration {it}. */
  
typedef double tmnn_goal_func_t(uint32_t n, double v[], double vOpt[]);
  /* Type of a "raw" function for optimization, that has a minimum at
    {v == vOpt}. */

void tmnn_test_minn_uniform
  ( uint32_t n, 
    tmnn_goal_func_t *F2,
    tmnn_goal_func_t *B2,
    bool_t box, 
    minn_method_t meth, 
    int32_t it,
    bool_t verbose
  );
  /* Tests {minn_uniform} with goal function {F2},
    bias function {B2}, and parameters {dMax,box,meth}. */

void tmnn_test_minn_subspace
  ( uint32_t n, 
    tmnn_goal_func_t *F2,
    tmnn_goal_func_t *B2,
    bool_t box, 
    minn_method_t meth, 
    int32_t it,
    bool_t verbose
  );
  /* Tests {minn_subspace} with goal function {F2},
    bias function {B2}, and parameters {box,meth}. */

void tmnn_test_minn_ellipsoid_constrained
  ( uint32_t n, 
    tmnn_goal_func_t *F2,
    tmnn_goal_func_t *B2,
    minn_method_t meth, 
    int32_t it,
    bool_t verbose
  );
  /* Tests {minn_ellipsoid_constrained} with goal function {F},
    bias function {B}, and parameters {box,meth}. */

void tmnn_choose_radii(uint32_t n, double zeros, double arad[]);
  /* Chooses finite non-negative search domain radii {arad[0..n-1]}.
    If {zeros} is true, some of the radii may be zero; else they are all positive. */

void tmnn_choose_tolerances
  ( uint32_t n,
    double arad[],
    bool_t uniform,
    minn_method_t meth,
    double atol[]
  );
  /* Chooses the tolerance(s) for the {minn.h} functions.
    
    The parameters {arad[0..n-1]} are assumed to be the radii
    of the search domain in each coordinate axis, all
    non-negative. Some may be 0, some may be {+INF}.
    
    The procedure chooses a suitable tolerance {atol[i]}
    for each coordinate axis {i} in {0..n-1}. 
    If {arad[i]} is zero, the toterance {atol[i]} will be zero.
    If {uniform} is true, all elements will be equal. */
    
void tmnn_choose_optimum
  ( uint32_t n, 
    double arad[], 
    double atol[],
    bool_t box,
    double vOpt[]
  );
  /* Chooses and stores into {vOpt[0..n-1]} the coefficiensts of the
    optimum solution, assuming that he search domain radii are
    {arad[0..n-1]} and the tolerances are {atol[0..n-1]}.
    
    If {arad} is {NULL}, then {rad} should be finite and positive and
    is assumed to be the domain radius {arad[i]} for all axes. 
    Otherwise {rad} should be {NAN}.
    
    If {atol} is {NULL}, then {tol} should be finite and positive and
    is assumed to be the tolerance {atol[i]} for all axes.
    Otherwise {tol} should be {NAN}. */
    
void tmnn_choose_constraints(uint32_t n, uint32_t q, double A[]);
  /* The parameter {A} must have {q*n} elements, and is interpreted
    as a {q} by {n} array, stored by rows.
    
    The procedure chooses a set of {q} random {n}-vectors that
    define some linear constraints of the search, and stores them as the
    rows of {A}.  Some of the constraints may be redundant and they may not 
    include the constraints implied by zero radii {arad[i]}. */

void tmnn_check_solution
  ( uint32_t n, 
    minn_goal_t *F,
    double vSol[],
    double FSol,
    double vOpt[], 
    bool_t box,
    double arad[],
    double atol[]
  );
  /* Checks the final solution {vSol[0..n-1]}. */
  
void tmnn_debug_points
  ( char *name, 
    int32_t k,
    uint32_t n, 
    double v[], 
    double vOpt[], 
    bool_t box,
    double arad[], 
    double atol[], 
    double Fv, 
    double Bv
  );
  /* Prints the probe vector {v}, its domain inclusion function,
    and its relation to {vOpt}, accounting for the domain radii {arad}
    and the toleances {atol}.  These two may be {NULL}.
    
    If {Bv} is not {NAN}, then {Fv} is assumed to be the "raw" function value, without
    the bias, and {Bv} the bias.  If {Bv} is {NAN}, then {Fv} is assumed to be the
    total function value, including the bias.
    
    The printout is prefixed with "{name][{k}] = " if {k >= 0},
    else with "{name} final = ". */

bool_t tmnn_too_many_points(uint32_t n, double arad[], double atol[]);
  /* Returns {TRUE} if the {minn_enum} method would make too many 
    function evaluations.  Assumes {box=TRUE}. */

int32_t main(int32_t argn, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {
    srandom(4615*417);
    
    demand(argc == 3, "invalid command line args");
    char *fname = argv[1];
    uint32_t n = (uint32_t)atoi(argv[2]);
    
    uint32_t nt = 100; /* Number of tests, */
    for (int32_t it = 0;  it < nt; it++)
      { bool_t verbose = (it < 10);
        tmnn_do_one_test(it, fname, n, verbose);
      }
    return 0;
  }

void tmnn_do_one_test (int32_t it, char *fname, uint32_t n, bool_t verbose)
  { 
    /* Goal functions for {r2_opt}: */

    tmnn_goal_func_t *F = NULL; /* */
      /* The raw goal function for optimization, without bias.
        Selected by {fname}. */

    double bias_mag = 1.0e-6;
      /* The magnitude of the bias term. */
    
   auto double B(uint32_t n, double v[], double vOpt[]);
      /* The bias term, a small quadratic function with 
        minimum at the origin. */
   
    /* Alternatives for {F_raw}: */

    auto double F_indiff(uint32_t n, double v[], double vOpt[]);
      /* A goal function that returns 1.0 always.  With bias,
        the minimum of this function should be at {(0,..0)}. */
    
    auto double F_optdst(uint32_t n, double v[], double vOpt[]);
      /* A goal function that retrurns the squared distance 
        from {v[0..n-1]} to {vOpt[0..n-1]} in a slightly distorted metric.
        With zero or small enough bias, the minimum of this function
        should be close to {vOpt}. */
    
    auto double F_optnoq(uint32_t n, double v[], double vOpt[]);
      /* A goal function similar to {F_optdst(n,v,vOpt)} but with terms
        of higher degree on {v-vOpt} that still has a quadratic minimum
        at {vOpt}. With zero or small enough bias, the minimum of this
        function should be close to {vOpt}. */
    
    /* Choose the target function: */
    if (strcmp(fname, "indiff") == 0)
      { F = &F_indiff; }
    else if (strcmp(fname, "optdst") == 0)
      { F = &F_optdst; }
    else if (strcmp(fname, "optnoq") == 0)
      { F = &F_optnoq; }
    else 
      { demand(FALSE, "invalid goal function name"); }
    
    for (int32_t ibox = 0;  ibox < 2; ibox++)
      { bool_t box = (ibox != 0);
        for (int32_t imeth = 0;  imeth < 2; imeth++)
          { minn_method_t meth = (imeth == 0 ? minn_method_ENUM : minn_method_QUAD);
            tmnn_test_minn_uniform(n, F, B, box, meth, it, verbose);
            tmnn_test_minn_subspace(n, F, B, box, meth, it, verbose);
            if (! box)
              { /* Has no {box} option, the domain is always an ellipsoid: */
                tmnn_test_minn_ellipsoid_constrained(n, F, B, meth, it, verbose);
              }
          }
      }
    return;
    
    /* INTERNAL IMPLS */
    
    double B(uint32_t n1, double v[], double vOpt[])
      { assert(n1 == n);
        double Bv = 0;
        if (bias_mag > 0)
          { double d2 = rn_norm_sqr(n, v); 
            Bv = bias_mag * d2;
          }
        return Bv;
      }
    
    double F_indiff(uint32_t n1, double v[], double vOpt[])
      { assert(n1 == n);
        double Fv = 1.0;
        return Fv;
      }
    
    double F_optdst(uint32_t n1, double v[], double vOpt[])
      { assert(n1 == n);
        /* Compute sum of relative square distances from actual optimum {opt}: */
        double Fv = 0.0;
        if (n > 0)
          { for (int32_t i = 0;  i < n; i++)
              { double di = v[i] - vOpt[i];
                int32_t j = (i + 1) % (int32_t)n;
                double dj = v[j] - vOpt[j];
                double dij = (1 + 0.25*i)*di + 0.33*dj;
                Fv += dij*dij;
              }
          }
        return Fv;
      }
      
    double F_optnoq(uint32_t n1, double v[], double vOpt[])
      { assert(n1 == n);
        double Fv = 0;
        if (n > 0)
          { for (int32_t i = 0;  i < n; i++)
              { double di = v[i] - vOpt[i];
                int32_t j = (i + 1) % (int32_t)n;
                double dj = v[j] - vOpt[j];
                double dij = (1 + 0.25*i)*di - 0.67*dj;
                double eij = exp(dij) - dij;
                Fv += eij;
              }
          }
        return Fv;
      }
  }
  
void tmnn_test_minn_uniform
  ( uint32_t n, 
    tmnn_goal_func_t *F,
    tmnn_goal_func_t *B,
    bool_t box, 
    minn_method_t meth, 
    int32_t it,
    bool_t verbose
  )
  { 
    bool_t debug_points = FALSE;

    if (verbose) 
      { fprintf(stderr, "--- testing {minn_uniform}");
        fprintf(stderr, " n = %d", n);
        fprintf(stderr, " box = %c meth = %c\n", "FT"[box], "EQ"[meth]);
      }
      
    /* Create an {arad} vector for convenience: */
    double arad[n];
    for (int32_t i = 0;  i < n; i++) { arad[i] = 1.0; }
    
    /* Choose the tolerance on each axis: */
    double atol[n]; /* Tolerance on each axis. */
    bool_t unif_tol = FALSE; /* Different tolerance on each axis. */
    tmnn_choose_tolerances(n, arad, unif_tol, meth, atol);
      
    if ((meth == minn_method_ENUM) && tmnn_too_many_points(n, arad, atol))
      { fprintf(stderr, "too many points in enum grid, skipped");
        return; 
      }
    
    /* Choose the optimum: */
    double vOpt[n]; /* Expected optimum point. */
    tmnn_choose_optimum(n, arad, atol, box, vOpt);
      
    uint32_t nF = 0; /* Counts calls to {nF}. */
    
    auto double F_minn(uint32_t n, double v[]);
      /* Goal function: just {F(n,v,vOpt)}. */
    
    /* Call the minimizer: */
    double vSol[n];   /* Computed minimum. */
    double FSol = NAN;
    minn_uniform(n, F_minn, box, atol, meth, vSol, &FSol);
    if (verbose)
      { fprintf(stderr, "did %d function evaluations\n", nF);
        tmnn_debug_points("v", -1, n, vSol, vOpt, box, arad, atol, FSol, NAN);
      }
    
    /* Check the solution: */
    tmnn_check_solution(n, F_minn, vSol, FSol, vOpt, box, arad, atol);
    
    return;
    
    /* INTERNAL IMPS */
    
    double F_minn(uint32_t n, double v[])
      { double Fv = F(n, v, vOpt);
        double Bv = B(n, v, vOpt);
        if (debug_points) 
          { tmnn_debug_points("v", (int32_t)nF, n, v, vOpt, box, arad, atol, Fv, Bv); }
        nF++;
        return Fv + Bv;
      }
  }
  
void tmnn_test_minn_subspace
  ( uint32_t n, 
    tmnn_goal_func_t *F,
    tmnn_goal_func_t *B,
    bool_t box, 
    minn_method_t meth, 
    int32_t it,
    bool_t verbose
  )
  { 
    bool_t debug_points = FALSE;

    if (verbose) 
      { fprintf(stderr, "--- testing {minn_subspace}");
        fprintf(stderr, " n = %d box = %c meth = %c\n", n, "FT"[box], "EQ"[meth]);
      }

    /* Choose the subspace: */
    uint32_t d = uint32_abrandom(0, n);
    double U[d*n];
    double urad[d];
    rmxn_throw_ortho_complement(n, 0, NULL, d, U);
    bool_t rad_zeros = FALSE; /* All radii must be positive. */
    tmnn_choose_radii(d, rad_zeros, urad);
    
    /* Choose the tolerance on each axis: */
    double utol[d]; /* Tolerance on each axis. */
    bool_t unif_tol = FALSE; /* Different tolerance on each axis. */
    tmnn_choose_tolerances(d, urad, unif_tol, meth, utol);
      
    if ((meth == minn_method_ENUM) && tmnn_too_many_points(d, urad, utol))
      { fprintf(stderr, "too many points in enum grid, skipped");
        return; 
      }
    
    /* Choose the optimum: */
    double uOpt[d]; /* Expected optimum point in the {U} basis. */
    tmnn_choose_optimum(d, urad, utol, box, uOpt);
    double vOpt[n]; /* Expected optimum point in {\RR^n}. */
    rmxn_map_row(d, n, uOpt, U, vOpt);
      
    uint32_t nF = 0; /* Counts calls to {nF}. */
    
    auto double F_minn(uint32_t n, double v[]);
      /* Goal function: just {F(n,v,vOpt)}. */
    
    /* Call the minimizer: */
    double vSol[n];   /* Computed minimum. */
    double FSol = NAN;
    minn_subspace(n, F_minn, d, U, urad, box, utol, meth, vSol, &FSol);
    if (verbose)
      { fprintf(stderr, "did %d function evaluations\n", nF);
        tmnn_debug_points("v", -1, n, vSol, vOpt, box, NULL, NULL, FSol, NAN);
      }

    /* Check the solution in the subspace {<U>}: */
    double uSol[d]; /* Solution coords in the {U} basis. */
    rmxn_map_col(d, n, U, vSol, uSol); 
    tmnn_check_solution(d, NULL, uSol, FSol, uOpt, box, urad, utol);

    /* Check the solution in {\RR^N}: */
    double vPro[n]; /* Solution {vSol} projected onto {<U>}. */
    rmxn_map_row(d, n, uSol, U, vPro);
    double dp = rn_dist(n, vSol, vPro);
    demand(dp < 1.0e-8, "solution is not in subspace {<U>}");
    tmnn_check_solution(n, F_minn, vSol, FSol, vOpt, box, NULL, NULL);
    
    return;
    
    /* INTERNAL IMPS */
    
    double F_minn(uint32_t n, double v[])
      { double Fv = F(n, v, vOpt);
        double Bv = B(n, v, vOpt);
        if (debug_points) 
          { tmnn_debug_points("v", (int32_t)nF, n, v, vOpt, box, NULL, NULL, Fv, Bv); }
        nF++;
        return Fv + Bv;
      }
  }    
    
void tmnn_test_minn_ellipsoid_constrained
  ( uint32_t n, 
    tmnn_goal_func_t *F,
    tmnn_goal_func_t *B,
    minn_method_t meth, 
    int32_t it,
    bool_t verbose
  )
  { 
    bool_t debug_points = FALSE;

    if (verbose) 
      { fprintf(stderr, "--- testing {minn_ellipsoid_constrained}");
        fprintf(stderr, " n = %d meth = %c\n", n, "EQ"[meth]);
      }

    bool_t box = FALSE; /* Domain is always an ellipsoid. */
    
    /* Choose the radii of the base ellipsoid: */
    double arad[n];
    bool_t rad_zeros = TRUE; /* Some radii may be zero. */
    tmnn_choose_radii(n, rad_zeros, arad);
    
    /* Choose the tolerances for nonzero axes: */
    double atol[n]; /* Tolerance on each axis. */
    bool_t unif_tol = TRUE; /* Same tolerance on every axis. */
    tmnn_choose_tolerances(n, arad, unif_tol, meth, atol);
    double tol = (n == 0 ? NAN : atol[0]);
    
    /* Choose the explicit constraints: */
    uint32_t q = uint32_abrandom(0, n+2);
    double A[q*n];
    tmnn_choose_constraints(n, q, A);
 
    /* Normalize the constraints and add the implicit ones: */
    uint32_t m;       /* Number of normalized constraints. */
    double *C = NULL; /* Normalized constraint matrix. */
    bool_t debug_norm = FALSE;
    rmxn_ellipsoid_normalize_constraints(n, arad, q, A, debug_norm, &m, &C);
    assert((m >= 0) && (m <= n));
   
    /* Compute the search ellipsoid {\RF(U,urad)}: */
    uint32_t d = n - m;
    double U[d*n];
    double urad[d];
    rmxn_ellipsoid_cut(n, arad, m, C, d, U, urad);
    
    /* Create a vector {utol[0..d-1]} for convenience: */
    double utol[d];
    for (int32_t k = 0;  k < d; k++) { utol[k] = tol; }

    /* Choose the optimum: */
    double uOpt[d]; /* Expected optimum point in the {U} basis. */
    tmnn_choose_optimum(d, urad, utol, box, uOpt);
    double vOpt[n]; /* Expected optimum point in {\RR^n}. */
    rmxn_map_row(d, n, uOpt, U, vOpt);
      
    uint32_t nF = 0; /* Counts calls to {nF}. */
    
    auto double F_minn(uint32_t n, double v[]);
      /* Goal function: just {F(n,v,vOpt)}. */
    
    /* Call the minimizer: */
    double vSol[n];   /* Computed minimum. */
    double FSol = NAN;
    minn_ellipsoid_constrained(n, F_minn, arad, q, A, tol, meth, vSol, &FSol);
    if (verbose)
      { fprintf(stderr, "did %d function evaluations\n", nF);
        tmnn_debug_points("v", -1, n, vSol, vOpt, box, arad, atol, FSol, NAN);
      }

    /* Check the solution in the subspace {<U>}: */
    double uSol[d]; /* Solution coords in the {U} basis. */
    rmxn_map_col(d, n, U, vSol, uSol); 
    tmnn_check_solution(d, NULL, uSol, FSol, uOpt, box, urad, NULL);

    /* Check the solution in the space {\RR^n}: */
    double vPro[n]; /* Solution {vSol} projected onto {<U>}. */
    rmxn_map_row(d, n, uSol, U, vPro);
    double dp = rn_dist(n, vSol, vPro);
    demand(dp < 1.0e-8, "solution is not in subspace {<U>}");
    tmnn_check_solution(n, F_minn, vSol, FSol, vOpt, box, arad, atol);
    
    return;
    
    /* INTERNAL IMPS */
    
    double F_minn(uint32_t n, double v[])
      { double Fv = F(n, v, vOpt);
        double Bv = B(n, v, vOpt);
        if (debug_points) 
          { tmnn_debug_points("v", (int32_t)nF, n, v, vOpt, box, arad, atol, Fv, Bv); }
        nF++;
        return Fv + Bv;
      }
  }    
    
void tmnn_choose_radii(uint32_t n, double zeros, double arad[])
  {
    double radMin = 0.50, radMax = 2.50; /* Range for {urad[i]}. */ 
    for (int32_t i = 0;  i < n; i++)
      { if (zeros && (i == 0 || (drandom() < 0.20)))
          { arad[i] = 0.0; }
        else
          { arad[i] = dabrandom(radMin, radMax); }
      }
  }

void tmnn_choose_tolerances
  ( uint32_t n,
    double arad[],
    bool_t uniform,
    minn_method_t meth,
    double atol[]
  )
  { double tolMin = +INF; /* Min {atol[i]} for axes with non-zero radius. */
    for (int32_t i = 0;  i < n; i++) 
      { double radi = arad[i];
        if (radi == 0)
          { atol[i] = 0.0; }
        else if (isfinite(radi))
          { /* Bounded search. If {meth} is enum, don't let {atol[i]} be too small: */
            double tolMin = (meth == minn_method_ENUM ? radi/4 : 0.01*radi);
            double tolMax = radi/3;
            atol[i] = dabrandom(tolMin, tolMax);
          }
        else
          { /* Unbounded search. Tolerance is irrelevant: */
            atol[i] = dabrandom(0.01, 0.10);
          }
        if (radi > 0) { tolMin = fmin(tolMin, atol[i]); }
      }
    if (uniform) 
      { assert(tolMin > 0);
        if (tolMin == +INF) { tolMin = 0.1; }
        for (int32_t i = 0;  i < n; i++) { atol[i] = tolMin; } 
      }
  }

void tmnn_choose_optimum
  ( uint32_t n,
    double arad[],
    double atol[],
    bool_t box,
    double vOpt[]
  )
  { if (box)
      { /* Try to put near a corner: */
        for (int32_t i = 0;  i < n; i++)
          { double radi = arad[i];
            double toli = atol[i];
            vOpt[i] = radi - 2*toli;
            assert(fabs(vOpt[i]) < radi - toli);
          }
      }
    else
      { /* Put anywhere inside the ball/ellipsoid: */
        rn_throw_ball(n, vOpt);
        rn_weigh(n, arad, vOpt, vOpt);
      }
  }   
  
#define tmnn_MAX_PTS 10000
  /* Max function evaluations allowd in enum search. */

bool_t tmnn_too_many_points(uint32_t n, double arad[], double atol[])
  { double np = 1;
    for (int32_t i = 0;  i < n; i++)
      { double radi = (arad == NULL ? +INF : arad[i]);
        assert(radi >= 0);
        double toli = (atol == NULL ? 0.0 : atol[i]);
        assert(isfinite(toli) && (toli >= 0.0));
        if (radi == 0)
          { continue; }
        else if ((toli == 0) || (radi == +INF))
          { return TRUE; }
        else 
          { double npi = floor(radi/toli + 1.0e-8);
            if (npi > tmnn_MAX_PTS/np) { return TRUE; }
            np = np*npi;
          }
      }
    return FALSE;
  }

void tmnn_debug_points
  ( char *name,
    int32_t k,
    uint32_t n, 
    double v[], 
    double vOpt[], 
    bool_t box,
    double arad[], 
    double atol[], 
    double Fv, 
    double Bv
  )
  { if (k >= 0)
      { fprintf(stderr, "  %s[%4d] = ", name, k); }
    else
      { fprintf(stderr, "  %s final = ", name); }
    rn_print(stderr, n, v);
    fprintf(stderr, "\n");
    
    /* Print function and bias values: */
    fprintf(stderr, "  function value");
    if (! isnan(Bv))
      { fprintf(stderr, " raw = %14.9f bias = %14.9f", Fv, Bv);
        fprintf(stderr, " total = %14.9f\n", Fv+Bv);
      }
    else
      { fprintf(stderr, " = %14.9f\n", Fv); }
    
    if (arad != NULL)
      { /* Check belonging to search domain: */
        double vr[n];
        for (int32_t i = 0;  i < n; i++)
          { double radi = arad[i];
            if (radi == +INF)
              { vr[i] = 0.0; }
            else if (radi == 0)
              { vr[i] = (v[i] == 0 ? 0.0 : +INF); }
            else
              { vr[i] = v[i]/radi; }
          }
        double vrnorm;
        if (box)
          { vrnorm = rn_L_inf_norm(n, vr); }
        else
          { vrnorm = rn_norm(n, vr); }
        fprintf(stderr, "  domain inclusion norm = %14.9f\n", vrnorm);
      }
    
    /* Print absolute distance from optimum: */
    double dopt = rn_dist(n, v, vOpt);
    fprintf(stderr, "  dist from optimum = %14.9f\n", dopt);
    
    if (atol != NULL)
      { /* Print distance from optimum relative to tolerances: */
        double vd[n];
        rn_sub(n, v, vOpt, vd);
        double vdt[n];
        rn_unweigh(n, atol, vd, vdt);
        double vdtnorm  = rn_norm(n, vdt);
        fprintf(stderr, "  tolerance-relative distance = %14.9f\n", vdtnorm);
      }
    fprintf(stderr, "\n");
  }
        
void tmnn_check_solution
  ( uint32_t n, 
    minn_goal_t *F,
    double vSol[],
    double FSol,
    double vOpt[], 
    bool_t box,
    double arad[],
    double atol[]
  )
  {
    fprintf(stderr, "!! %s NOT IMPLEMENTED !!\n", __FUNCTION__); 
  }

void tmnn_choose_constraints(uint32_t n, uint32_t q, double A[])
  { 
    for (int32_t i = 0;  i < q; i++)
      { for (int32_t j = 0;  j < n; j++)
          { A[i*(int32_t)n + j] = dabrandom(-2.0, +2.0); }
      }
  }
