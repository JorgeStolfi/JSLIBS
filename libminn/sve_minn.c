/* See {sve_minn.h} */
/* Last edited on 2024-09-15 15:40:16 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include <gauss_elim.h>
#include <bool.h>
#include <affirm.h>
#include <vec.h>
#include <jsmath.h>
#include <rmxn.h>
#include <rmxn_extra.h>
#include <rmxn_regular_simplex.h>

#include <sve_minn.h>

#define Pr fprintf
#define Er stderr

/* INTERNAL PROTOTYPES */

void sve_clip_candidate
  ( int32_t n,
    double y[],
    double ctr[],
    double dMax,
    bool_t box, 
    bool_t debug
  );
  /* Pulls the point {y[0..n-1]} towards the point {ctr[0..n-1]} so that 
    it lies within the search domain defined by {dMax} and {box}. 
    If {ctr} is {NULL}, assumes a vector of {n} zeros. */

int32_t sve_take_step
  ( int32_t n,
    sign_t dir,
    double y[],
    double *FyP,
    double x[],
    double Fx,
    int32_t nv,
    double v[],
    int32_t nf,
    double Fv[],
    bool_t debug
  );
  /* Chooses the next next probe center for {sve_minn_iterate}.

    The candidates for the next probe center are the estimated
    stationary point {y[0..n-1]} computed by {sve_minn_step}, the
    current probe center {x[0..n-1]}, and the probe points of this
    iteration (namely the vertices {V(i)} and edge midpoints {V(i,j)}
    of the current probe simplex).
    
    Upon entry, assumes that {y} has been clipped as appropriate; that
    {Fx} and {*FyP} are the function values at {x} and {y},
    respectively; that {v} contains the coordinates of the simplex
    vertices, by rows; and that {Fv} contains the function values at
    all probe points.
    
    The chosen next center is stored in {y[0..n-1]} and its function
    value in {*FyP}.
    
    If {dir} is zero, the chosen point is always the quadratic 
    statonary point; that is, the procedure is a no-op.
    Otherwise it is the candidate with minimum (if {dir==-1}) or maximum 
    (if {dir==+1}) function value.
    
    The returned value is 0 if the winner is the estimated minimum
    {y}, 1 if it is one of the probe point, and 2 if it is the
    original center itself (i.e. if the step failed altogether). */

void sve_print_probes(FILE *wr, int32_t nv, int32_t n, double v[], int32_t nf, double Fv[]);
  /* Prints the probe values {Fv[0..nf-1]} and the probe points.
    Assumes that {v[0..nv*n-1]} are the coordinates of the vertices,
    stored by rows. */

void sve_minn_step(int32_t n, double Fv[], double cm[], bool_t debug)
  { int32_t nv = n+1; /* Number of vertices in simplex. */
    int32_t rows = nv+1; /* {n+1} stationary eqs, and one unit-sum eq. */
    int32_t cols = nv+2; /* {n+1} barycentric coords, one Lagrange multip, and the indep term. */
    double M[rows*cols];
    /* Fill in the main submatrix: */
    int32_t jb = cols - 1; /* Column of independent term. */
    for (int32_t i = 0; i < nv; i++)
      { for (int32_t j = 0; j < nv; j++)
          { int32_t ij = (j <= i ? i*(i+1)/2 + j : j*(j+1)/2 + i);
            double Fij = Fv[ij];
            M[i*cols + j] = Fij;
            if (i == j) { M[i*cols + jb] = Fij/4; }
          }
      }
    /* Fill in the row {nv} and column {nv} with the unit-sum constraint: */
    int32_t ije = n+1; /* Index of constraint row & column. */
    for (int32_t i = 0; i <= n; i++) { M[i*cols + ije] = 1; M[ije*cols + i] = 1; }
    M[ije*cols + ije] = 0;
    M[ije*cols + jb] = 1;
    if (debug)
      { rmxn_gen_print(stderr, rows, cols, M, "%12.7f", "  [ ", "\n    ", " ]\n", "[ ", " ", " ]"); }
    /* Solve the system: */
    gsel_triangularize(rows, cols, M, TRUE, 0.0);
    gsel_diagonalize(rows, cols, M);
    gsel_normalize(rows, cols, M);
    double x[n+2];
    int32_t rank_ext = gsel_extract_solution(rows, cols, M, 1, x);
    if (rank_ext < rows) 
      { Pr(Er, "%s: warning - solution with %d degrees of indeterminacy\n", __FUNCTION__, rows - rank_ext); }
    /* Extract the solution: */
    double sum = 0.0;
    for (int32_t i = 0; i <= n; i++) { cm[i] = x[i]; sum += cm[i]; }
    if (fabs(sum - 1.0) > 0.5e-7) 
      { Pr(Er, "%s: warning - normalization failed, sum = %24.16e\n", __FUNCTION__, sum); }
    /* Just to be sure: */
    if ((sum != 0) && (sum != 1)) { for (int32_t i = 0; i <= n; i++) { cm[i] /= sum; } }
  }

void sve_sample_function(int32_t n, sve_goal_t *F, double v[], double Fv[])
  { double x[n];
    int32_t nv = n + 1;
    for (int32_t i = 0; i < nv; i++)
      { for (int32_t j = 0; j <= i; j++)
          { /* Set {x[0..n-1]} to the midpoint of simplex corners {i,j}: */
            double *vi = &(v[i*n]);
            double *vj = &(v[j*n]);
            for (int32_t k = 0; k < n; k++) { x[k] = (vi[k] + vj[k])/2; }
            /* Get the function's value {F(x)} at {x}, store {F(x)} into {Fv}: */
            int32_t ij = i*(i+1)/2 + j;
            Fv[ij] = F(n, x);
          }
      }
  }

void sve_minn_iterate
  ( int32_t n, 
    sve_goal_t *F, 
    sve_pred_t *OK,
    double x[],
    double *FxP,
    sign_t dir, 
    double ctr[],
    double dMax,
    bool_t box,
    double rIni,
    double rMin, 
    double rMax,
    double stop,
    int32_t maxIters,
    bool_t debug
  )
  { 
    if (debug) { Pr(Er, ">> enter %s >>\n", __FUNCTION__); }
    bool_t debug_probes = TRUE;     /* TRUE to print probe points & values. */
    bool_t debug_termination = TRUE; /* TRUE to print iteration termination info. */
    
    demand((rMin >= sve_minn_MIN_RADIUS) && (rMin <= rMax), "invalid {rMin,rMax}");
    demand((rMin <= rIni) && (rIni <= rMax), "invalid {rIni}");
    
    /* Allocate storage for the simplex: */
    int32_t nv = n+1; /* Number of vertices in simplex. */
    double_vec_t vv = double_vec_new(nv*n); 
    double *v = vv.e; /* Cartesian coords of simplex vertices. */
    
    /* Allocate storage for the sample values: */
    int32_t nf = (n+1)*(n+2)/2; /* Number of probe points. */
    double_vec_t Fvv = double_vec_new(nf);
    double *Fv = Fvv.e; /* Sampled function values. */

    /* Iteration variables: */
    double cm[nv]; /* Barycentric coords of next solution. */
    double y[n];    /* Cartesian coords of next solution; also simplex center. */
    double radius = rIni; /* Current probe simplex radius: */
    double dPrev = -1; /* Distance moved in previous iteration (-1 if none). */
    int32_t nIters = 0; /* Counts quadratic step iterations. */
    int32_t nEvals = 0; /* Counts function evaluations. */
    
    /* Get initial function value: */
    double Fx = (*FxP);
    
    while(TRUE) 
      { if (debug) 
          { Pr(Er, "  iteration %4d\n", nIters);
            Pr(Er, "    current x = \n");
            rn_gen_print(Er, n, x, "%20.16f", "    [ ", "\n      ", "   ]\n");
            Pr(Er, "    function = %22.16e\n", Fx);
          }
        /* Check for termination: */
        if ((OK != NULL) && (OK(n, x, Fx))) 
          { if (debug_termination) { Pr(Er, "  client is satisfied\n"); }
            break;
          }
        /* Check budget: */
        if (nIters >= maxIters) 
          { if (debug_termination) { Pr(Er, "  iteration limit exhausted\n"); }
            break;
          }
          
        /* QUADRATIC OPTIMIZATION STEP */
        
        /* Ensure that the probe simplex radius is in {[rMin _ rMax]}: */
        if (debug)  { Pr(Er, "  raw radius = %12.8f\n", radius); }
        if (radius < rMin)
          { radius = rMin;
            if (debug)  { Pr(Er, "  radius augmented to {rMin} = %12.8f\n", radius); }
          }
        if (radius > rMax) 
          { radius = rMax; 
            if (debug)  { Pr(Er, "  radius reduced to {rMax} = %12.8f\n", radius); }
          }
        
        /* The probe simplex center {y} is in principle the current guess {x}: */
        rn_copy(n, x, y);
        
        /* Adjust simplex center {y} and {radius} so that there is room for the probe simplex: */
        if (dMax < INFINITY)
          { /* Get distance {dist} from domain center {ctr} to current center candiate {y}: */
            double dist;
            if (ctr == NULL)
              { dist = (box ? rn_L_inf_norm(n, y) : rn_norm(n, y)); }
            else
              { dist = (box ? rn_L_inf_dist(n, ctr, y) : rn_dist(n, ctr, y)); }
            if (debug) { Pr(Er, "  dist(y, ctr) = %16.12f radius = %16.12f sum = %16.12f\n", dist, radius, dist+radius); }
            /* Ensure that the simplex fits in the domain box/ball: */
            double dExtra = dist + radius - dMax; 
            if (dExtra > 0.0)
              { /* Simplex risks falling out of {ctr,dMax} ball/box. */
                if (debug) { Pr(Er, "  simplex sphere overshoots box/ball by %16.12f\n", dExtra); }
                /* Compute new radius {rSafe} such that the ball 
                  with center at the adjusted point will fit inside that ball/box. */
                double rSafe = radius - 0.50000001*dExtra;
                if (rSafe < rMin) { rSafe = rMin; }              
                if (rSafe > dMax) { rSafe = dMax; }
                assert((rMin <= rSafe) && (rSafe <= rMax)); /* Should be always the case. */
                /* Adjust {y} so that it is at {dMax-rSafe} from {ctr}: */
                sve_clip_candidate(n, y, ctr, dMax - rSafe, box, debug);
                /* Use this adjusted radius: */
                assert(rSafe <= radius); /* Should always be the case. */
                radius = rSafe;
                if (debug) { Pr(Er, "  radius adjusted for {dMax} = %12.8f\n", radius); }
              }
            else
              { Pr(Er, "  simplex is inside domain box/ball\n"); }
          }
        if(debug) { Pr(Er, "\n"); }

        /* Choose a regular simplex around {y}: */
        rmxn_regular_simplex(n, v);
        rmxn_spin_rows(nv, n, v, v);
        /* Compute the scaling factor {scale} for the regular simplex: */
        double scale = radius/rmxn_regular_simplex_radius(n);
        rmxn_scale(nv, n, scale, v, v);
        rmxn_shift_rows(nv, n, v, y, v);
        /* Evaluate the goal function at sample points: */
        if (debug && debug_probes)
          { /* Print the current simplex: */
            Pr(Er, "    simplex vertices:\n");
            rmxn_gen_print(Er, nv, n, v, "%20.16f", "", "\n", "", "    [ ", " ", " ]");
            Pr(Er, "\n");
            for (int32_t iv = 0; iv < nv; iv++)
              { double *vi = &(v[iv*n]);
                Pr(Er, "    dist(v[%3d], y) = %16.12f\n", iv, rn_dist(n, vi, y));
              }
          } 
        sve_sample_function(n, F, v, Fv);
        if (debug && debug_probes)
          { /* Print the probe points and values: */
            Pr(Er, "    probe values and points:\n");
            sve_print_probes(Er, nv, n, v, nf, Fv); 
            Pr(Er, "\n");
          }
        nEvals += nf;
        /* Optimize as a quadratic obtaining barycentric coords {cm[0..nv-1]}: */
        sve_minn_step(n, Fv, cm, debug);
        nIters++;
        /* Convert solution {cm} to Cartesian coordinates {y[0..n-1]}: */
        rmxn_map_row(nv, n, cm, v, y);
        if (debug) 
          { Pr(Er, "    raw y = \n");
            rn_gen_print(Er, n, y, "%20.16f", "    [ ", "\n      ", " ]\n");
            Pr(Er, "\n");
          }
        
        /* STEP ADJUSTMENT */
        
        if (dMax < INFINITY)
          { /* Adjust the estimated min {y} to be inside the sphere/box {ctr,dMax}: */
            sve_clip_candidate(n, y, ctr, dMax, box, debug);
          }
        /* Evaluate at new point: */
        double Fy = F(n, y);
        nEvals++;
        if (debug) 
          { Pr(Er, "    clipped y = \n");
            rn_gen_print(Er, n, y, "%20.16f", "    [ ", "\n      ", " ]\n");
            Pr(Er, "    function = %22.16e\n", Fy);
            Pr(Er, "\n");
          }
        /* Set {y} to the best of all points seen so far: */
        int32_t stepKind = sve_take_step(n, dir, y, &Fy, x, Fx, nv, v, nf, Fv, debug);
        double dStep = rn_dist(n, x, y); /* Length of this step: */

        /* Update the point {x}: */
        rn_copy(n, y, x);
        Fx = Fy;
        if (debug) 
          { Pr(Er, "    choice: %s\n", ((char*[3]){"quadratic min", "sample point", "old center"})[stepKind]);
            Pr(Er, "    new guess = \n");
            rn_gen_print(Er, n, x, "%20.16f", "    [ ", "\n      ", " ]\n");
            Pr(Er, "    function = %22.16e\n", Fx);
            Pr(Er, "    step len = %22.16e\n", dStep);
            Pr(Er, "\n");
          }
          
        /* Estimate the distance {dEst} to the minimum */
        double dEst;
        if (stepKind == 0)
          { /* Quadratic minimization succeeded, assume geometric convergence: */
            dEst = (dPrev <= dStep ? +INF : dStep*dStep/(dPrev - dStep));
          }
        else if (stepKind == 1)
          { /* Moved to a simplex vertex: */
            dEst = +INF;
          }
        else if (stepKind == 2)
          { /* Current guess is less than simplex: */
            dEst = radius;
          }
        else
          { assert(FALSE); }
        /* Check for convergence: */
        bool_t dEst_is_small = (dEst < stop);
        bool_t no_progress = ((dStep < stop) && (radius <= rMin));
        if (dEst_is_small || no_progress) 
          { if (debug_termination) 
              { Pr(Er, "  seems to have converged:\n");
                Pr(Er, "    dPrev = %22.16e\n", dPrev);
                Pr(Er, "    dStep = %22.16e\n", dStep);
                Pr(Er, "    dEst =  %22.16e\n", dEst);
              }
            break;
          }
        /* The motion is at least the center-to-edge radius of the simplex: */
        /* Remember how much we moved in this iteration: */
        dPrev = dStep;
        /* Adjust simplex radius for next iteration: */
        double rNew = dStep/2;
        if (rNew < 0.250*radius) { rNew = 0.250*radius; }
        if (rNew > 2.0*radius) { rNew = 2.0*radius; }
        radius = rNew;
      }
    if (debug) { Pr(Er, "  did %d iterations and %d function evaluations.\n", nIters, nEvals); }
    free(v);
    free(Fv);
    
    /* Return final function value: */
    (*FxP) = Fx;
    if (debug) { Pr(Er, "<< leave %s <<\n", __FUNCTION__); }
  }

void sve_clip_candidate
  ( int32_t n,
    double y[],
    double ctr[],
    double dMax,
    bool_t box,
    bool_t debug
  )
  {
    double dist;
    if (ctr == NULL)
      { dist = (box ? rn_L_inf_norm(n, y) : rn_norm(n, y)); }
    else
      { dist = (box ? rn_L_inf_dist(n, ctr, y) : rn_dist(n, ctr, y)); }
    if (dist > dMax) 
      { /* Moved too much, curb it: */
        double s = (1.0 - 1.0e-12)*dMax/dist;
        if (debug) 
          { Pr(Er, "  moved too far from initial guess d = %20.16e\n", dist);
            Pr(Er, "  contracting towards initial guess by s = %20.16f\n", s);
          }
        for (int32_t k = 0; k < n; k++) { y[k] = (1-s)*(ctr == NULL ? 0.0 : ctr[k]) + s*y[k]; }
      }
  }

int32_t sve_take_step
  ( int32_t n,
    sign_t dir,
    double y[],
    double *FyP,
    double x[],
    double Fx,
    int32_t nv,
    double v[],
    int32_t nf,
    double Fv[],
    bool_t debug
  )
  {
    /* If looking for a stationary point, keep {y} and {*FyP}: */
    if (dir == 0) { return 0; }
    
    /* Grab the current function value at {y}: */
    double Fy = (*FyP);
    
    /* Find the optimum value among {Fv[0..nf-1]}: */
    double FOpt = -INF*(double)dir; /* {FOpt} is the optimum sample value. */
    int32_t iOpt = -1, jOpt = -1; /* {V(iOpt,jOpt)} is the optimum sample point. */
    for (int32_t i = 0; i < nv; i++)
      { for (int32_t j = 0; j <= i; j++)
          { int32_t ij = i*(i+1)/2 + j;
            if (dir*Fv[ij] >= dir*FOpt) { iOpt = i; jOpt = j; FOpt = Fv[ij]; }
          }
      }
    if (nv > 0) { assert((iOpt >= 0) && (jOpt >= 0)); }

    /* Choose the next center {x} and set the nominal step {d}: */
    if ((dir*Fx > dir*Fy) || (dir*FOpt > dir*Fy))
      { /* Quadratic step failed: */
        if (debug) { Pr(Er, "    quadratic step failed - not optimum\n"); }
        if (dir*Fx >= dir*FOpt)
          { /* The minimum is still {x}: */
            if (debug) { Pr(Er, "    the optimum is still the simplex center\n"); }
            rn_copy(n, x, y);
            (*FyP) = Fx;
            return 2;
          }
        else
          { /* The minimum is {V(iOpt,jOpt)}: */
            if (debug) { Pr(Er, "    the optimum is V(%d,%d)\n", iOpt, jOpt); }
            double *vi = &(v[iOpt*n]);
            double *vj = &(v[jOpt*n]);
            for (int32_t k = 0; k < n; k++) { y[k] = (vi[k] + vj[k])/2; }
            (*FyP) = FOpt;
            return 1;
          }
      }
    else
      { /* Quadratic step seems to have succeeded: */ 
        if (debug) { Pr(Er, "    quadratic step yielded a new optimum\n"); }
        return 0;
      }
  }

void sve_print_probes(FILE *wr, int32_t nv, int32_t n, double v[], int32_t nf, double Fv[])
  {
    for (int32_t i = 0; i < nv; i++)
      { for (int32_t j = 0; j <= i; j++)
          { int32_t ij = i*(i+1)/2 + j;
            fprintf(wr, "    %24.16e", Fv[ij]);
            fprintf(wr, "  ");
            double *vi = &(v[i*n]);
            double *vj = &(v[j*n]);
            for (int32_t k = 0; k < n; k++) { fprintf(wr, " %20.16f", (vi[k] + vj[k])/2); }
            fprintf(wr, "\n");
          }
      }
  }
