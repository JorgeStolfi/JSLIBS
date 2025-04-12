/* See {sve_minn.h} */
/* Last edited on 2025-04-02 09:44:03 by stolfi */

#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <affirm.h>
#include <vec.h>
#include <jsmath.h>
#include <rmxn.h>
#include <rmxn_spin.h>
#include <rmxn_shift.h>
#include <rmxn_regular_simplex.h>
#include <gausol_print.h>
#include <gausol_solve.h>

#include <sve_minn.h>
#include <sve_minn_iterate.h>

#define Pr fprintf
#define Er stderr

/* INTERNAL PROTOTYPES */

bool_t sve_clip_candidate
  ( uint32_t n,
    double y[],
    double dCtr[],
    double dMax,
    bool_t dBox, 
    bool_t debug
  );
  /* Pulls the point {y[0..n-1]} towards the point {dCtr[0..n-1]} so that 
    it lies within the search domain defined by {dMax} and {dBox}. 
    If {dCtr} is {NULL}, assumes a vector of {n} zeros.
    Returns true iff {y} was changed. */

uint32_t sve_take_best_candidate
  ( sign_t dir,
    uint32_t n, 
    double x[], double *FxP,  /* Previous guess and its value. */
    double y[], double Fy,    /* Computed stationary point and value. */
    uint32_t nv, double v[], uint32_t nf, double Fv[],
    double *dStep,
    bool_t debug
  );
  /* Chooses the best candidate point among the previous solution 
    and the probe points (candidates) examined by {sve_minn_single_step}.

    The candidates for the next probe center are the previous best
    solution {x[0..n-1]}, the estimated stationary point {y[0..n-1]} of
    the approximating quadratic, as computed by
    {sve_minn_quadratic_optimum}, and the probe points of this iteration
    (namely the vertices {V(i)} and edge midpoints {V(i,j)} of the
    current probe simplex).
    
    Upon entry, assumes that the computed stationary point {y} has been
    clipped to the search domain, and that {*FxP} and {Fy} are the
    function values at {x} and {y}, respectively. It also assumes that
    {v[0..nv*n-1]} contains the (unprojected) coordinates of the simplex
    vertices, by rows; and that {Fv[0..nf-1]} contains the function
    values at all {nf} probe points.
    
    If the best candidate is not {x}, its coordinates will be stored in
    {x[0..n-1]}, its function value in {*FxP}, and the Euclidean
    distance between it and the original {x} in in {*dStepP}. Otherwise
    {x} and {*FxP} are not changed and {dStepP} is set to zero.
    
    If {dir} is zero, the chosen point is always the computed stationary
    point {y}. Otherwise it is the candidate with minimum (if {dir==-1})
    or maximum (if {dir==+1}) function value.
    
    The returned value is 0 if the winner is the computed stationary point
    {y}, 1 if it is one of the probe point, and 2 if it is the
    previous best solution {x} (i.e. if the step failed altogether).  */

void sve_print_probes(FILE *wr, uint32_t nv, uint32_t n, double v[], uint32_t nf, double Fv[]);
  /* Prints the probe values {Fv[0..nf-1]} and the probe points.
    Assumes that {v[0..nv*n-1]} are the coordinates of the vertices,
    stored by rows. */

void sve_minn_quadratic_optimum(uint32_t n, double Fv[], double cm[], bool_t debug, bool_t debug_system)
  { uint32_t nv = n+1; /* Number of vertices in simplex. */
    uint32_t rows = nv+1; /* {n+1} stationary eqs and one unit-sum eq. */
    uint32_t cols = nv+1; /* {n+1} barycentric coords and one Lagrange multip. */
    double M[rows*cols]; /* Main systems matrix. */
    double b[rows];  /* RHS vector. */
    /* Fill in the main submatrix: */
    assert(cols >= 1);
    for (uint32_t i = 0; i < nv; i++)
      { for (uint32_t j = 0;  j < nv; j++)
          { uint32_t ij = (j <= i ? i*(i+1)/2 + j : j*(j+1)/2 + i);
            double Fij = Fv[ij];
            M[i*cols + j] = Fij;
            if (i == j) { b[i] = Fij/4; }
          }
      }
    /* Fill in the row {nv} and column {nv} with the unit-sum constraint: */
    uint32_t ije = n+1; /* Index of constraint row & column. */
    for (uint32_t i = 0;  i <= n; i++) 
      { M[i*cols + ije] = 1;
        M[ije*cols + i] = 1;
      }
    M[ije*cols + ije] = 0;
    b[ije] = 1;
    if (debug_system)
      { gausol_print_system
          ( stderr, 4, "%12.7f", "quadratic step system {M x = b}:",
            rows,NULL,n, cols,NULL,n, "M",M, 1,"b",b, 0,NULL,NULL, ""
          );
      }
    /* Solve the system: */
    double x[rows];
    uint32_t rank_ext;
    gausol_solve(rows, cols, M, 1, b, x, TRUE, TRUE, 0.0, NULL, &rank_ext);
    if (debug_system)
      { fprintf(stderr, "    rank = %d  Lagrange mult = %24.16e\n", rank_ext, x[ije]);
        gausol_print_array
          ( stderr, 4, "%17.10f", "quadratic system raw solution:", 
            rows,NULL,n, 1,NULL,0, "x",x, ""
          );
      }
    assert(rank_ext <= rows);
    if (rank_ext < rows) 
      { Pr(Er, "%s: warning - solution with %d degrees of indeterminacy\n", __FUNCTION__, rows - rank_ext); }
    /* Check unit sum condition: */
    double sum = 0.0;
    for (uint32_t i = 0;  i <= n; i++) { cm[i] = x[i]; sum += cm[i]; }
    if (fabs(sum - 1.0) > 0.5e-7) 
      { Pr(Er, "%s: warning - unit-sum constraint violated, sum = %24.16e\n", __FUNCTION__, sum); }
    /* Just to be sure: */
    if ((sum != 0) && (sum != 1)) { for (uint32_t i = 0;  i <= n; i++) { cm[i] /= sum; } }
  }

void sve_sample_function(uint32_t n, sve_goal_t *F, double v[], double Fv[])
  { double x[n];
    uint32_t nv = n + 1;
    for (uint32_t i = 0;  i < nv; i++)
      { for (uint32_t j = 0;  j <= i; j++)
          { /* Set {x[0..n-1]} to the midpoint of simplex corners {i,j}: */
            double *vi = &(v[i*n]);
            double *vj = &(v[j*n]);
            for (uint32_t k = 0;  k < n; k++) { x[k] = (vi[k] + vj[k])/2; }
            /* Get the function's value {F(x)} at {x}, store {F(x)} into {Fv}: */
            uint32_t ij = i*(i+1)/2 + j;
            Fv[ij] = F(n, x);
          }
      }
  }

uint32_t sve_minn_single_step
  ( uint32_t n, 
    double x[], double *FxP,
    sve_goal_t *F, sign_t dir,
    double dCtr[], double dMax, bool_t dBox, 
    double *radiusP, double rMin,
    double v[], double Fv[],
    uint32_t *nEvalsP,
    double *dStepP,
    bool_t debug,
    bool_t debug_probes
  )
  { 
    demand(rMin <= dMax, "invalid {rMin}, larger than {dMax}");

    uint32_t nv = n + 1;          /* Number of simplex vertices. */
    uint32_t nf = (n+1)*(n+2)/2;  /* Number {nf} of probe points. */
    
    double y[n];   /* Simplex center, then computed stationary point. */
    double cm[nv]; /* Barycentric coords of stationary point. */
    /* The probe simplex center {y} is in principle the current optimum {x}: */
    rn_copy(n, x, y);

    /* Compute distance from {dCtr}: */
    double dist;
    if (dCtr == NULL)
      { dist = (dBox ? rn_L_inf_norm(n, y) : rn_norm(n, y)); }
    else
      { dist = (dBox ? rn_L_inf_dist(n, dCtr, y) : rn_dist(n, dCtr, y)); }

    /* Adjust simplex center {y} and {radius} so that it fits in domain: */
    double radius = (*radiusP);
    if (debug)  { Pr(Er, "      dist(x, dCtr) = %16.12f  raw simplex radius = %12.8f\n", dist, radius); }
    demand((radius >= rMin) && (radius <= dMax), "invalid simplex {radius}");
    if (dMax < INFINITY)
      { /* Ensure that the simplex fits in the domain box/ball: */
        double dExtra = dist + radius - dMax; 
        if (dExtra > 0.0)
          { /* Simplex risks falling out of {dCtr,dMax} ball/box. */
            if (debug) { Pr(Er, "      dist plus simplex radius =  %16.12f  exceeds {dMax} by %16.12f\n", dist+radius, dExtra); }
            /* Compute new radius {rSafe} such that the ball 
              with center at the adjusted point will fit inside that ball/box. */
            double rSafe = radius - 0.50000001*dExtra;
            if (rSafe < rMin) { rSafe = rMin; }              
            if (rSafe > dMax) { rSafe = dMax; }
            /* Adjust {y} so that it is at {dMax-rSafe} from {dCtr}: */
            bool_t y_clipped = sve_clip_candidate(n, y, dCtr, dMax - rSafe, dBox, debug);
            if(y_clipped && debug_probes)
              { Pr(Er, "      clipped simplex center y = \n");
                rn_gen_print(Er, n, y, "%20.16f", "    [ ", "\n      ", " ]\n");
              }
            /* Use this adjusted radius: */
            assert(rSafe <= radius); /* Should always be the case. */
            radius = rSafe;
            if (debug) { Pr(Er, "      simplex radius reduced to %12.8f\n", radius); }
            (*radiusP) = radius;
          }
      }
    if(debug) { Pr(Er, "\n"); }

    /* Choose a regular simplex of radius {radius} around {y}: */
    rmxn_regular_simplex(n, v);
    rmxn_spin_rows(nv, n, v, v);
    double scale = radius/rmxn_regular_simplex_radius(n);
    rmxn_scale(nv, n, scale, v, v);
    rmxn_shift_rows(nv, n, v, y, v);
    
    /* Evaluate the goal function at sampling points: */
    if (debug_probes)
      { /* Print the current simplex: */
        Pr(Er, "      simplex vertices:\n");
        rmxn_gen_print(Er, nv, n, v, "%20.16f", "", "\n", "", "        [ ", " ", " ]");
        Pr(Er, "\n");
        for (uint32_t iv = 0;  iv < nv; iv++)
          { double *vi = &(v[iv*n]);
            Pr(Er, "        dist(v[%3d], y) = %16.12f\n", iv, rn_dist(n, vi, y));
          }
      } 
    sve_sample_function(n, F, v, Fv);
    if (debug_probes)
      { /* Print the probe points and values: */
        Pr(Er, "      probe values and points:\n");
        sve_print_probes(Er, nv, n, v, nf, Fv); 
        Pr(Er, "\n");
      }
    if (nEvalsP != NULL) { (*nEvalsP) += nf; }
    
    /* Optimize as a quadratic obtaining barycentric coords {cm[0..nv-1]}: */
    bool_t debug_system = debug_probes;
    sve_minn_quadratic_optimum(n, Fv, cm, debug, debug_system);

    /* Convert barycentric coords {cm} to Cartesian coords {y[0..n-1]}: */
    rmxn_map_row(nv, n, cm, v, y);
    if (debug) 
      { Pr(Er, "      computed quadratic stationary point y = \n");
        rn_gen_print(Er, n, y, "%20.16f", "    [ ", "\n      ", " ]\n");
        Pr(Er, "\n");
      }

    if (dMax < INFINITY)
      { /* Adjust the estimated min {y} to be inside the sphere/box {dCtr,dMax}: */
        bool_t y_clipped = sve_clip_candidate(n, y, dCtr, dMax, dBox, debug);
        if(y_clipped && debug)
          { Pr(Er, "      clipped stationary point y = \n");
            rn_gen_print(Er, n, y, "%20.16f", "    [ ", "\n      ", " ]\n");
          }
      }
      
    /* Evaluate at new point: */
    double Fy = F(n, y);
    if (nEvalsP != NULL) { (*nEvalsP)++; }
    if (debug) { Pr(Er, "      function at stationary point {y} = %22.16e\n", Fy); }

    uint32_t stepKind = sve_take_best_candidate(dir, n, x, FxP, y, Fy, nv, v, nf, Fv, dStepP, debug);
    return stepKind;
  }

bool_t sve_clip_candidate
  ( uint32_t n,
    double y[],
    double dCtr[],
    double dMax,
    bool_t dBox,
    bool_t debug
  )
  {
    double dist;
    if (dCtr == NULL)
      { dist = (dBox ? rn_L_inf_norm(n, y) : rn_norm(n, y)); }
    else
      { dist = (dBox ? rn_L_inf_dist(n, dCtr, y) : rn_dist(n, dCtr, y)); }
    if (debug) { Pr(Er, "        distance from center = %20.16e\n", dist); }
            
    if (dist > dMax) 
      { /* Point {y} is outside domain, curb it: */
        double s = (1.0 - 1.0e-12)*dMax/dist;
        if (debug) 
          { Pr(Er, "        moved too far - contracting towards {dCtr} by s = %20.16f\n", s); }
        for (uint32_t k = 0;  k < n; k++) { y[k] = (1-s)*(dCtr == NULL ? 0.0 : dCtr[k]) + s*y[k]; }
        return TRUE;
      }
    else
      { return FALSE; }
  }

uint32_t sve_take_best_candidate
  ( sign_t dir,
    uint32_t n, 
    double x[], double *FxP,
    double y[], double Fy,
    uint32_t nv, double v[], uint32_t nf, double Fv[],
    double *dStepP,
    bool_t debug
  )
  {
    /* If looking for a stationary point, keep {y} and {*FyP}: */
    if (dir == 0) 
      { /* Ignore simplex vertices and previous guess: */
        rn_copy(n, y, x);
        (*FxP) = Fy;
        return 0;
      }
    
    /* Grab the current function value at {y}: */
    double Fx = (*FxP);
    
    /* Find the optimum value among the values {Fv[0..nf-1]} sampled at simplex nodes: */
    double FOpt = -INF*(double)dir; /* {FOpt} is the optimum sample value. */
    int32_t iOpt = -1, jOpt = -1; /* {V(iOpt,jOpt)} is the optimum sample point. */
    for (int32_t i = 0;  i < nv; i++)
      { for (int32_t j = 0;  j <= i; j++)
          { int32_t ij = i*(i+1)/2 + j;
            if (dir*Fv[ij] >= dir*FOpt) 
              { iOpt = i; jOpt = j; FOpt = Fv[ij]; }
          }
      }
    if (nv > 0) { assert((iOpt >= 0) && (jOpt >= 0)); }

    /* Choose the best guess: */
    uint32_t stepKind; 
    if ((dir*Fx > dir*Fy) && (dir*Fx > dir*FOpt))
      { if (debug) { Pr(Er, "      step failed - the optimum is still the input guess {x}\n"); }
        stepKind = 2;
        (*dStepP) = 0;
      }
    else 
      { if (dir*FOpt > dir*Fy)
          { if (debug) { Pr(Er, "      the optimum is the sampling point V(%d,%d)\n", iOpt, jOpt); }
            double *vi = &(v[iOpt*(int32_t)n]);
            double *vj = &(v[jOpt*(int32_t)n]);
            double vOpt[n];
            for (uint32_t k = 0; k < n; k++) { vOpt[k] = (vi[k] + vj[k])/2; }
            (*dStepP) = rn_dist(n, x, vOpt);
            rn_copy(n, vOpt, x);
            (*FxP) = FOpt;
            stepKind = 1;
          }
        else
          { if (debug) { Pr(Er, "      the optimum is the (clipped) quadratic stationary point {y}\n"); }
            (*dStepP) = rn_dist(n, x, y);
            rn_copy(n, y, x);
            (*FxP) = Fy;
            stepKind = 0;
          }
      }
    return stepKind;
  }

void sve_print_probes(FILE *wr, uint32_t nv, uint32_t n, double v[], uint32_t nf, double Fv[])
  {
    for (uint32_t i = 0;  i < nv; i++)
      { for (uint32_t j = 0;  j <= i; j++)
          { uint32_t ij = i*(i+1)/2 + j;
            fprintf(wr, "        %24.16e", Fv[ij]);
            fprintf(wr, "  ");
            double *vi = &(v[i*n]);
            double *vj = &(v[j*n]);
            for (uint32_t k = 0;  k < n; k++) { fprintf(wr, " %20.16f", (vi[k] + vj[k])/2); }
            fprintf(wr, "\n");
          }
      }
  }
