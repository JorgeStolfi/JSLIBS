/* test_sve_step --- tests functions of {sve_minn_h} */
/* Last edited on 2025-04-02 08:25:43 by stolfi */

#include <values.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>

#include <vec.h>
#include <rn.h>
#include <rmxn.h>
#include <rmxn_throw.h>
#include <rmxn_spin.h>
#include <rmxn_shift.h>
#include <rmxn_canonical_simplex.h>
#include <rmxn_regular_simplex.h>
#include <sym_eigen.h>
#include <affirm.h>
#include <jsrandom.h>
#include <bool.h>

#include <sve_minn.h>

/* GENERAL PARAMETERS */

#define MAX_RUNS 100
  /* Max number of trials per test. */

#define MAX_DIM 10
  /* Max number of corners in simplex. */

#define MAX_VARS 10
  /* Max number of variables in goal function. */

/* INTERNAL PROTOTYPES */

int32_t main (int32_t argc, char **argv);

void test_all(uint32_t n, uint32_t trial, bool_t nice, bool_t verbose);
  /* Tests the routines of {sve_minn.h} with some function of {n}
    variables. If {nice} is TRUE, uses a nice test function and simplex, otherwise
    uses random ones. */

void test_sve_sample_function(uint32_t n, sve_goal_t *goal, bool_t nice, bool_t verbose);
  /* Tests {sve_sample_function} on the goal function {F} of {n} variables
    and a random simplex in {\RR^n}.  If {nice}, the simplex will be regular and 
    centered at the origin. */

void test_sve_minn_quadratic_optimum(uint32_t n, sve_goal_t *goal, double xRef[], bool_t nice, bool_t verbose);
  /* Tests the routine {sve_minn_quadratic_optimum} with the values of {goal}
    evaluated at the corners and mid-edges of a random simplex of {\RR^n}.
    Checks whether the result is {xRef[0..n-1]}. 
    
    If {nice} is true, the simplex will be regular and centered at {xRef[0..n-1],
    otherwise it will be irregular and centered some distance away from {xRef}. */

void test_sve_minn_single_step(uint32_t n, sve_goal_t *goal, double xRef[], bool_t dBox, sign_t dir, bool_t nice, bool_t verbose);
  /* Tests the routine {sve_minn_single_step} with the given {goal} function of {n}
    variables, assuming its minimum is at {xRef}. 
    
    The domain will be a cube if {dBox} is true, or a ball if {dBox} is false.
    
    If {nice} is true, the search domain will be centered at {xRef[0..n-1],
    otherwise it will be centered some distance away from {xRef}. */

void tsve_pick_problem(uint32_t n, bool_t nice, double A[], double xRef[], double *CP, sign_t *dirP, bool_t verbose);
  /* These procedures pick the parameters {xRef,A,C} of a quadratic
     function with the formula {F(x)=(x-xRef)'*A*(x-xRef)+C}. Uses
     {tsve_pick_nice_problem} if {nice} is true, or
     {tsve_pick_random_problem} if {nice} is false.
     
     In any case, the stationary point {xRef} will be in the cube {U^n}
     where {U==[-1_+1]}. The entries of the coefficient matrix {A} and
     the constant term {C} will be in {U}. The parameter {C} will be
     stored into {*CP}.
     
     Also stores in {*dirP} the type of stationary solution: {-1} for
     minimum, {+1} for maximum, {0} for saddle point. */

void tsve_check_problem_matrix(uint32_t n, double A[], sign_t dir, bool_t verbose);
  /* Checks whether the eigenvalues of {A} are all well away from zero.
    If {dir} is not zero, also checks whether their sign is opposite to {dir}. */

void tsve_print_problem(FILE *wr, uint32_t n, double A[], double xRef[], double C, sign_t dir);
  /* Writes the elements {xREf,A,C,dir} of the quadtratic function to {wr}. */

void tsve_throw_simplex(uint32_t n, double rSim, double xRef[], bool_t nice, double v[], bool_t verbose);
  /* Picks an {n}-dimensional simplex {v} in {\R^n}.
    If {nice} is true, the simplex will  be regular, with radius {rSim} and center at {xRef[0..n-1]}.
    If {nice} is false, the simplexy will be irregular, with radius at most {rSim}, centered at 
    a random point that lies at distance at most {rSim} from {xRef}.  If {xRef} is {NULL},
    assumes it is the origin. */
     
void tsve_throw_domain_center(uint32_t n, double dMax, bool_t dBox, double xRef[], double rMin, bool_t nice, double dCtr[], bool_t verbose);
  /* Chooses the center {dCtr[0..n-1]} of the search domain for {sve_minn_single_step}
    radius {dMax}.  The domain will be a box if {dBox} is true, or a ball if {dBox} is false.
    The domain center will be the optimum {xRef[0..n-1]} if {nice} is true, 
    or some random point such that {xRef} is inside the domain. */
   
void tsve_throw_initial_point(uint32_t n, double dCtr[], double dMax, bool_t dBox, bool_t nice, double xIni[], bool_t verbose);
  /* Chooses an initial point {xIni[0..n-1]} for {sve_minn_single_step} in the domain with center {dCtr} and
    radius {dMax}.  The domain will be a box if {dBox} is true, or a ball if {dBox} is false.
    The initial point will be the domain center if {nice} is true, or a random point in
    the domain if {nice} is false. */

void tsve_check_function_values(uint32_t n, sve_goal_t *goal, double v[], double Fv[]);
  /* Compares the values of {goal} at the vertices and midpoints of a simplex
    with vertices {v[0..(n+1)*n-1]} with the values {Fv[0..nf-1]} where {nf = (n+1)*(n+2)/2}. */

void tsve_check_optimum_position(uint32_t n, double xRef[], double xCmp[]);
  /* Compares the computed optimum {xCmp[0..n-1]} against the 
    expected optimum {xRef[0..n-1]} */

double tsve_eval_quadratic(uint32_t n, double A[], double xRef[], double C, const double x[]);
  /* Evaluates the function {F(x)=(x-xRef)'*A*(x-xRef) + C} at the point {x[0..n-1]}. */

void tsve_convert_to_cartesian(uint32_t n, double cm[], double v[], double xCmp[], bool_t verbose);
  /* Converts the computed stationary point from barycentric coordinates
    {cm[0..n]}, relative to the simplex {v[0..(n+1)*n-1]} of {\RR^n}, to
    Cartesian coordinates {xCmp[0..n-1]}. */

double tsve_max_abs_elem(uint32_t cols, uint32_t rows, double N[]);
  /* Returns the maximum absolute value of the elements in the matrix {N} with {cols} columns and {rows} rows. */

double tsve_max_abs_col_elem(uint32_t cols, uint32_t rows, double N[], uint32_t j);
  /* Returns the maximum absolute value of the elements in column {j} of
     the matrix {N} with {cols} columns and {rows} rows. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  { for (uint32_t trial = 0;  trial < MAX_RUNS; trial++) 
      { /* Choose num of variables {n}, niceness, and verbosity: */
        uint32_t n = (trial < 3 ? trial : uint32_abrandom(1, MAX_VARS));
        bool_t nice = (trial < 5);
        bool_t verbose = (trial < 6);
        test_all(n, trial, nice, verbose);
      }
    fprintf(stderr, "done.\n");
    return 0;
  }
    
void test_all(uint32_t n, uint32_t trial, bool_t nice, bool_t verbose)
  { 
    fprintf(stderr, "entering %s  n = %d trial = %d nice = %c\n", __FUNCTION__, n, trial, "FT"[nice]);
    
    srand(1665 + 2*trial);
    srandom(1665 + 2*trial);

    if (verbose) { fprintf(stderr, "  generating problem...\n"); }
    double A[n*n];  /* Main coefficient matrix of {F}. */
    double xRef[n]; /* True stationary point. */
    double C;       /* Value of function {F} at {xRef}. */
    sign_t dir;
    
    tsve_pick_problem(n, nice, A, xRef, &C, &dir, verbose); 
    bool_t dBox = ((trial % 2) == 1); /* For {sve_minn_single_step}. */
    
    auto double goal(uint32_t n, const double x[]);

    test_sve_sample_function(n, goal, nice, verbose);
    test_sve_minn_quadratic_optimum(n, goal, xRef, nice, verbose);
    test_sve_minn_single_step(n, goal, xRef, dBox, dir, nice, verbose);
    
    fprintf(stderr, "exiting %s\n", __FUNCTION__);

    return;

    double goal(uint32_t n, const double x[])
      {
        return tsve_eval_quadratic(n, A, xRef, C, x);
      }
  }

void test_sve_sample_function(uint32_t n, sve_goal_t *goal, bool_t nice, bool_t verbose)
  { 
    uint32_t nv = n+1;
    if (verbose) { fprintf(stderr, "  --- testing {sve_sample_function} ---\n"); }
    if (verbose) { fprintf(stderr, "  n = %d  nice = %c ---\n", n, "FT"[nice]); }

    if (verbose) { fprintf(stderr, "  choosing the probe simplex...\n"); }
    double v[nv*n]; /* Elem {v[i*n + j]} is coord {j} of simplex vertex {i}. */
    if (nice)
      { rmxn_regular_simplex(n, v); }
    else
      { double rSim = dabrandom(1.0, 9.0);
        tsve_throw_simplex(n, rSim, NULL, nice, v, verbose);
      }

    if (verbose) { fprintf(stderr, "  gathering the function samples...\n"); }
    uint32_t nf = nv*(n+2)/2;
    double Fv[nf]; /* Function values at simplex corners and edge midpoints. */
    sve_sample_function(n, goal, v, Fv);

    if (verbose) { fprintf(stderr, "  checking values...\n"); }
    tsve_check_function_values(n, goal, v, Fv);
    
    if (verbose) { fprintf(stderr, "  ---\n\n"); }
  }

void test_sve_minn_quadratic_optimum(uint32_t n, sve_goal_t *goal, double xRef[], bool_t nice, bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "  --- testing {sve_minn_quadratic_optimum} ---\n"); }
    if (verbose) { fprintf(stderr, "  n = %d  nice = %c ---\n", n, "FT"[nice]); }
    
    uint32_t nv = n+1;
    if (verbose) { fprintf(stderr, "  generating the simplex with %d vertices ...\n", nv); }
    double *v = talloc(nv*n, double);  /* Storage for the simplex. */
    double rSim = 1.0;
    tsve_throw_simplex(n, rSim, xRef, nice, v, verbose);
    
    uint32_t nf = nv*(n+2)/2;
    if (verbose) { fprintf(stderr, "  sampling the function at the %d probe points ...\n", nf); }
    double *Fv = talloc(nf, double);   /* Storage for function samples. */
    sve_sample_function(n, goal, v, Fv);
    
    if (verbose) { fprintf(stderr, "  computing the stationary point...\n"); }
    double cm[nv];
    bool_t sve_debug = verbose;
    bool_t sve_debug_system = verbose;
    sve_minn_quadratic_optimum(n, Fv, cm, sve_debug, sve_debug_system);
    if (verbose) 
      { rn_gen_print(stderr, nv, cm, "%12.8f", "  barycentric coords:\n  [ ", "\n    ", " ]\n"); }
    
    if (verbose) { fprintf(stderr, "  converting barycentric to Cartesian...\n"); }
    double xCmp[n]; 
    tsve_convert_to_cartesian(n, cm, v, xCmp, verbose);

    if (n >= 1)
      { /* Optimum should be unique: */
        if (verbose) { fprintf(stderr, "  checking position of computed optimum...\n"); }
        tsve_check_optimum_position(n, xRef, xCmp);
      }
    
    if (verbose) { fprintf(stderr, "  ---\n\n"); }
  }

void test_sve_minn_single_step(uint32_t n, sve_goal_t *goal, double xRef[], bool_t dBox, sign_t dir, bool_t nice, bool_t verbose)
  { 
    if (verbose) { fprintf(stderr, "  --- testing {sve_minn_single_step} ---\n"); }
    if (verbose) { fprintf(stderr, "  n = %d  dBox = %c  nice = %c ---\n", n, "FT"[dBox], "FT"[nice]); }
    
    uint32_t nv = n+1;
    uint32_t nf = nv*(n+2)/2;
    if (verbose) { fprintf(stderr, "  allocating array for %d vertices and %d sample values ...\n", nv, nf); }
    double *v = talloc(nv*n, double);  /* Storage for the simplex. */
    double *Fv = talloc(nf, double);    /* Storage for function samples. */
    
    double rMin = 1.0e-5;               /* Min simplex radius, for numerical stability. */
    
    if (verbose) { fprintf(stderr, "  choosing the search domain ...\n"); }
    double dMax = dabrandom(1.0, 9.0);  /* Search domain radius. */
    double dCtr[n];                    /* Domain center. */
    tsve_throw_domain_center(n, dMax, dBox, xRef, rMin, nice, dCtr, verbose);

    if (verbose) { fprintf(stderr, "  choosing the starting point ...\n"); }
    double xIni[n];          /* Initial guess. */
    tsve_throw_initial_point(n, dCtr, dMax, dBox, nice, xIni, verbose);
    
    if (verbose) { fprintf(stderr, "  performing the single step ...\n"); }
    double x[n];
    rn_copy(n, xIni, x);  
    double Fx = goal(n, x);
    
    double radius = 0.9*dMax;            /* Initial simplex radius. */
    uint32_t nEvals = 10000;
    double dStep = NAN;
    bool_t debug = verbose;
    bool_t debug_probes = verbose;
    sve_minn_single_step
      ( n, x, &Fx, goal, dir, dCtr, dMax, dBox, &radius, rMin, 
        v, Fv, &nEvals, &dStep, debug, debug_probes
      );
      
    demand(nEvals == 10000 + nf + 1, "wrong {nEvals}");

    if (n >= 1)
      { /* Optimum should be unique: */
        if (verbose) { fprintf(stderr, "  checking position of computed optimum...\n"); }
        tsve_check_optimum_position(n, xRef, x);
      }
    
    if (verbose) { fprintf(stderr, "  ---\n\n"); }
  }
  
void tsve_pick_problem(uint32_t n, bool_t nice, double A[], double xRef[], double *CP, sign_t *dirP, bool_t verbose)
  { 
    /* Pick the stationary point {xRef}: */
    for (uint32_t j = 0;  j < n; j++)
      { xRef[j] = (nice ? ((double)j+1)/((double)n+1) : dabrandom(-1.0, +1.0)); }

    /* Pick a diagonal coefficient matrix {A}: */
    for (uint32_t i = 0;  i < n; i++) 
      { for (uint32_t j = 0;  j < n; j++)
          { A[i*n + j] = (i == j ? dabrandom(0.2, 1.0) : 0.0); }
      }
      
    /* Choose a stationary point type and adjust the matrix: */
    sign_t dir; 
    if (n == 0)
      { /* Function will be constant: */
        dir = 0;
      }
    else
      { /* Random {dir}, but not 0 if {n} is 1. */
        do { dir = (sign_t)int32_abrandom(-1,+1); } while ((dir == 0) && (n == 1)); 
        switch(dir)
          { case +1: /* Function must have a maximum: */
              /* Negate all elements of diagonal {A}: */
              for (uint32_t i = 0; i < n; i++) { A[i*n + i] = -A[i*n + i]; }
              break;
            case -1: /* Function must have a minimum: */
              /* Nothing to do: */
              break;
            case 0: /* Does not matter if min or max, so make it saddle: */
              /* Negate last diag elem and some middle ones: */
              assert(n >= 2);
              A[n*n-1] = - A[n*n-1];
              for (uint32_t i = 1; i < n-1; i++) 
                { if (drandom() < 0.5) { A[i*n + i] = -A[i*n + i]; } }
              break;
            default:
              assert(FALSE);
          }
      }
    
    
    if (! nice)
      { if (verbose)
          { rmxn_gen_print
              ( stderr, n, n, A, "%12.8f", 
                "  coefficient matrix {A} before twisting:\n  [ ", "\n    ", " ]\n",
                "[ ", " ", " ]"
              );
          }
        /* Twist the matrix: */
        double R[n*n], M[n*n];
        rmxn_throw_ortho(n, R);
        rmxn_mul(n, n, n, A, R, M);
        rmxn_tr_mul(n, n, n, R, M, A);
      }
      
    tsve_check_problem_matrix(n, A, dir, verbose);

    /* Choose the function value at {xRef}: */
    double C = (nice ? 0.5 : dabrandom(-1.0, +1.0));

    if (verbose) { tsve_print_problem(stderr, n, A, xRef, C, dir); }
      
    (*dirP) = dir;
    (*CP) = C;
  }
      
void tsve_check_problem_matrix(uint32_t n, double A[], sign_t dir, bool_t verbose)
  {
    double d[n];
    uint32_t nev;
    double R[n*n];
    sym_eigen(n, A, d, R, &nev);
    if (verbose)
      { fprintf(stderr, "    num eigenvalues of {A} = %d\n", nev);
        rn_gen_print(stderr, nev, d, "%12.8f", "    eigenvalues:\n  [ ", "\n    ", " ]\n");
      }
    
    affirm(nev == n, "matrix {A} is not full rank");
    for (int32_t i = 0; i < n; i++)
      { demand(fabs(d[i]) > 1.0e-3, "eigenvalue too small");
        if (dir != 0) { demand(dir*d[i] < 0, "eigenvalue has wrong sign"); }
      }
  }

void tsve_print_problem(FILE *wr, uint32_t n, double A[], double xRef[], double C, sign_t dir)
  { 
    rn_gen_print(wr, n, xRef, "%12.8f", "  true stationary point {xRef}:\n  [ ", "\n    ", " ]\n");
    fprintf(wr, "  value {C} at stationary point = %24.16e\n", C);
    fprintf(wr, "  stationary value type {dir} = %+d\n", dir);
    rmxn_gen_print
      ( wr, n, n, A, "%12.8f", 
        "  coefficient matrix {A}:\n  [ ", "\n    ", " ]\n",
        "[ ", " ", " ]"
      );
  }

void tsve_throw_simplex(uint32_t n, double rSim, double xRef[], bool_t nice, double v[], bool_t verbose)
  { /* Generate a regular {n}-simplex in {R^n}, with unit radius and center at origin: */
    uint32_t nv = n+1;
    rmxn_regular_simplex(n, v);
    
    if (! nice)
      { /* Scale each vertex by a random amount less than 1: */
        for (uint32_t iv = 0; iv < nv; iv++)
          { double *vi = &(v[iv*n]);
            double scale = dabrandom(0.1, 1.0);
            rn_scale(n, scale, vi, vi);
          }
        
        /* Spin the simplex to random orientation: */
        rmxn_spin_rows(nv, n, v, v);
        
        /* Choose a center at most 1 away from origin and shift it there: */
        double dCtr[n];
        rn_throw_ball(n, dCtr);
        rmxn_shift_rows(nv, n, v, dCtr, v);
     }

    /* Scale the simplex (and shift) to radius {rSim}: */
    rmxn_scale(nv, n, rSim, v, v);
    
    if (xRef != NULL)
      { /* Shift everything by {xRef}: */
        rmxn_shift_rows(nv, n, v, xRef, v);
      }

    if (verbose) 
      { rmxn_gen_print
          ( stderr, nv, n, v, "%12.8f", 
            "  simplex vertices:\n  [ ", "\n    ", " ]\n",
            "[ ", " ", " ]"
          );
      }
  }

void tsve_check_function_values(uint32_t n, sve_goal_t *goal, double v[], double Fv[])
  { uint32_t nv = n+1;
    double tol = 1.0e-6;
    double x[n];
    int32_t i01 = 0;
    for (int32_t i0 = 0; i0 < nv; i0++)
      { double *v0 = &(v[i0*(int32_t)n]);
        for (int32_t i1 = 0;  i1 <= i0; i1++)
          { double *v1 = &(v[i1*(int32_t)n]);
            rn_mix(n, 0.5, v0, 0.5, v1, x);
            double Fx_cmp = Fv[i01];
            double Fx_ref = goal(n, x);
            double diff = Fx_cmp - Fx_ref;
            if (fabs(diff) > tol)
              { fprintf(stderr, "**\n");
                rn_gen_print(stderr, n, x, "%12.8f", "  argument vector x:\n  [ ", "\n    ", " ]");
                if (i0 == i1) 
                  { fprintf(stderr, " = v[%d]\n", i0); }
                else
                  { fprintf(stderr, " = (v[%d]+v[%d])/2\n", i0, i1); }
                fprintf(stderr, "F(x) ref =   %22.16e\n", Fx_ref);
                fprintf(stderr, "F(x) cmp =   %22.16e\n", Fx_cmp);
                fprintf(stderr, "difference = %22.16e\n", diff);
                fprintf(stderr, "tolerance =  %22.16e\n", tol);
                demand(FALSE, "** computed sample value {Fv[i01]} does not match {F(x)}");
              }
            i01++;
          }
      }
  }

void tsve_check_optimum_position(uint32_t n, double xRef[], double xCmp[])
  { double tol = 1.0e-6 * sqrt((rn_norm_sqr(n, xRef) + rn_norm_sqr(n, xCmp))/2);
    for (uint32_t j = 0; j < n; j++)
      { double Xcmpj = xCmp[j];
        double Xrefj = xRef[j];
        double diffj = Xcmpj - Xrefj;
        if (fabs(diffj) > tol)
          { fprintf(stderr, "xCmp[%d] =   %22.16e\n", j, xCmp[j]);
            fprintf(stderr, "xRef[%d] =   %22.16e\n", j, xRef[j]);
            fprintf(stderr, "difference = %22.16e\n", diffj);
            fprintf(stderr, "tolerance =  %22.16e\n", tol);
            demand(FALSE, "** computed optimum does not match reference sol");
          }
      }
  }

double tsve_eval_quadratic(uint32_t n, double A[], double xRef[], double C, const double x[])
  { double F = C;
    for (uint32_t i = 0;  i < n; i++)
      { for (uint32_t j = 0;  j < n; j++)
          { F += (x[i] - xRef[i])*A[i*n + j]*(x[j] - xRef[j]); }
      }
    return F;
  }

void tsve_convert_to_cartesian(uint32_t n, double cm[], double v[], double xCmp[], bool_t verbose)
  { uint32_t nv = n+1;
    for (uint32_t j = 0;  j < n; j++)
      { double Xsum = 0;
        for (uint32_t i = 0;  i < nv; i++) { Xsum += cm[i] * v[i*n + j]; }
        xCmp[j] = Xsum;
      }
    if (verbose) 
      { rn_gen_print(stderr, n, xCmp, "%12.8f", "  computed optimum:\n  [ ", "\n    ", " ]\n"); }
  }

double tsve_max_abs_elem(uint32_t cols, uint32_t rows, double N[])
  { double emax = 0.0;
    for (uint32_t i = 0;  i < rows; i++)
      { for (uint32_t j = 0;  j < cols; j++)
          { double Mij = fabs(N[i*cols + j]);
            if (Mij > emax) { emax = Mij; }
          }
      }
    return emax;
  }

double tsve_max_abs_col_elem(uint32_t cols, uint32_t rows, double N[], uint32_t j)
  { double emax = 0.0;
    for (uint32_t i = 0;  i < rows; i++)
      { double Mij = fabs(N[i*cols + j]);
        if (Mij > emax) { emax = Mij; }
      }
    return emax;
  }
  
void tsve_throw_domain_center(uint32_t n, double dMax, bool_t dBox, double xRef[], double rMin, bool_t nice, double dCtr[], bool_t verbose)
  { 
    /* Start with the optimum point: */
    rn_copy(n, xRef, dCtr);
    
    if (! nice)
      { double dx[n]; 
        if (dBox) { rn_throw_cube(n, dx); } else { rn_throw_ball(n, dx); }
        rn_scale(n, 0.9*dMax, dx, dx);
        rn_add(n, dx, xRef, dCtr);
      }

    if (verbose) 
      { rn_gen_print(stderr, n, dCtr, "%12.8f", "  domain center {dCtr}:\n  [ ", "\n    ", " ]\n"); }
  }

void tsve_throw_initial_point(uint32_t n, double dCtr[], double dMax, bool_t dBox, bool_t nice, double xIni[], bool_t verbose)
  { 
    /* Start with the domain center: */
    rn_copy(n, dCtr, xIni);
    
    if (! nice) 
      { double dx[n]; 
        if (dBox) { rn_throw_cube(n, dx); } else { rn_throw_ball(n, dx); }
        rn_scale(n, 0.9*dMax, dx, dx);
        rn_add(n, dx, xIni, xIni);
      }

    if (verbose) 
      { rn_gen_print(stderr, n, xIni, "%12.8f", "  starting guess {xIni}:\n  [ ", "\n    ", " ]\n"); }
   }
