/* See {r2_opt.h}. */
/* Last edited on 2024-11-08 09:53:51 by stolfi */

#define _GNU_SOURCE
#include <math.h>
#include <limits.h>
#include <assert.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
  
#include <bool.h>
#include <sign.h>
#include <r2.h>
#include <i2.h>
#include <rn.h>
#include <jsmath.h>
#include <affirm.h>

#include <sve_minn.h>

#include <r2_opt.h>

/* INTERNAL PROTOTYPES */

bool_t r2_opt_coord_is_variable(double arij, double asij);
  /* Returns {TRUE} iff a coordinate with search radius {arij}
    and search step {asij} is variable; {FALSE} if it is fixed. */

/* IMPLEMENTATIONS */

void r2_opt_single_scale_enum
  ( int32_t ni,                   /* Number of points to optimize. */
    i2_t iscale,              /* Object scaling exponent along each axis. */  
    r2_opt_goal_func_t *f2,   /* Function that evaluates the goal function. */
    r2_t arad[],              /* Max coordinate adjustment for each point. */
    r2_t astp[],              /* Adjustment step for each point. */
    r2_t p[],                 /* (IN/OUT) Corresponding points in each object. */
    double *f2p,              /* (OUT) Goal function value at the best solution. */
    bool_t debug              /* Prints diagnostics if true. */
  )
  {
    if (debug) { fprintf(stderr, ">> enter %s >>\n", __FUNCTION__); }
    r2_t p0[ni]; /* Saved initial guess. */

     /* The enumeration variables {v[0..nv-1]} are the
       displacements from {p0[0..ni-1]} to {p[0..ni-1]},
       divided by the corresponding steps.  Each {v[k]} 
       is at most {r[k]} in absolute value. */

    int32_t nc = 2*ni; /* Number of coordinates (including fixed ones). */
    int32_t v[nc];     /* Enumeration variables. */
    int32_t r[nc];     /* Limits for the enumeration variables. */
    int32_t nv = 0;    /* Number of variable coordinates. */

    /* Save {p0}, compute the integer radii {r[k]} and initialize {v[k]}: */
    for (int32_t i = 0; i < ni; i++)
      { p0[i] = p[i];
        for (int32_t j = 0; j < 2; j++)
          { int32_t k = 2*i + j;
            double arij = arad[i].c[j];
            double asij = astp[i].c[j];
            if (r2_opt_coord_is_variable(arij, asij))
              { r[k] = (int32_t)floor(arij/asij);
                affirm(r[k] > 0, "bug in r[k]");
                nv++;
              }
            else
              { r[k] = 0; }
            v[k] = -r[k];
          }
      }
    if (debug) { fprintf(stderr, "%d coordinates (%d variable)\n", nc, nv); }
    
    /* Enumerate all valid vectors {v[0..nc-1]}: */
    (*f2p) = +INF;   /* Minimum mismatch found so far. */
    r2_t pv[nc];     /* Trial alignment vector. */
    while (TRUE) 
      { /* Compute {pv} from {v}: */
        for (int32_t i = 0; i < ni; i++)
          { for (int32_t j = 0; j < 2; j++)
              { int32_t k = 2*i + j;
                pv[i].c[j] = p0[i].c[j] + v[k]*astp[i].c[j];
              }
          }
        /* Evaluate the function at {pv}: */
        double f2v = f2(ni,pv,iscale);
        if (f2v < (*f2p))
          { /* Update the current optimum: */
            if (debug) { fprintf(stderr, "found better solution f2 = %24.15e\n", f2v); }
            for (int32_t i = 0; i < ni; i++) { p[i] = pv[i]; }
            (*f2p) = f2v;
          }
        /* Increment the next {v[k]} that can be incremented, reset previous ones to min: */
        { int32_t k = 0;
          while ((k < nc) && (v[k] >= r[k])) { v[k] = -r[k]; k++; }
          if (k >= nc){ /* Done: */ return; }
          v[k]++;
        }
      }
    if (debug) { fprintf(stderr, "<< leave %s <<\n", __FUNCTION__); }
  }

void r2_opt_single_scale_quadopt
  ( int32_t ni,                   /* Number of points to optimize. */
    i2_t iscale,              /* Object scaling exponent along each axis. */  
    r2_opt_goal_func_t *f2,   /* Function that evaluates the goal function. */
    r2_t arad[],              /* Max coordinate adjustment for each point. */
    r2_t astp[],              /* Desired adjustment precision for each point. */
    r2_t p[],                 /* (IN) Initial guesses; (OUT) Best solution. */
    double *f2p,              /* (OUT) Goal function value at the best solution. */
    bool_t debug              /* Prints diagnostics if true. */
  )
  {
    if (debug) { fprintf(stderr, ">> enter %s >>\n", __FUNCTION__); }
    int32_t maxIters = 10;
    
    /* Find the number {nv} of variables to optimize and the relative precision {tol}: */
    int32_t nv = 0;
    double tol = +INF;
    for (int32_t i = 0; i < ni; i++) 
      { for (int32_t j = 0; j < 2; j++)
          { double arij = arad[i].c[j];
            double asij = astp[i].c[j];
            if (r2_opt_coord_is_variable(arij, asij))
              { /* This coordinate is variable: */
                nv++;
                assert(arij > 0.0);
                double tij = asij/arij; /* Rel precision required in this coordinate. */
                if (tij < tol) { tol = tij; }
              }
          }
      }
    
    if (nv > 0)
      { 
        /* Apply nonlinear optimization. */

        /* The optimization variables {x[0..nv-1]} are the coordinates
          of the displacements from {p0[0..ni-1].c[j]} to {p[0..ni-1].c[j]},
          except the fixed ones, divided by the corresponding radii {arad[i].c[j]}.
          Their range is therefore the box  {[-1 _ +1]^nv}. */

        /* Save initial guess {p} in {p0}: */
        r2_t p0[ni]; /* Saved initial guess. */
        for (int32_t i = 0; i < ni; i++) { p0[i] = p[i]; }
        
        /* These functions assume that the initial guess was saved in {p0[0..ni-1]}: */ 

        auto void points_to_vars(r2_t q[], double y[]);
          /* Stores the candidate displacements {q[0..ni-1]-p0[0..ni-1]} into {y[0..nv-1]}. */

        auto void vars_to_points(double y[], r2_t q[]);
          /* Stores the optimization variables {y[0..nv-1]} into the candidate {q[0..ni-1]}, adding {p0[0..ni-1]}. */

        auto double f2_for_sve(int32_t nx, double x[]);
          /* Computes the minimization goal function from the given argument {x[0..nv-1]}.
            Expects {nx == nv}. Also sets {p[0..ni-1]} from {x[0..nx-1]}. */

        /* Compute the initial goal function value: */
        double z[nv];
        points_to_vars(p, z);
        double Fz = f2_for_sve(nv, z);
        
        /* Optimize: */
        sign_t dir = -1;                   /* Look for minimum. */
        /* All these parameters are realative to the search radius in each coord: */
        double *ctr = NULL;                /* Center of search domain is the origin. */
        double dMax = 1.0;                 /* Max deviation from initial guess. */
        bool_t dBox = TRUE;                /* Search in box, not ball. */
        double rIni = 0.5;                 /* Initial probe simplex radius. */
        double rMin = fmin(tol, 0.25);     /* Minimum probe simplex radius. */
        double rMax = 0.70;                /* Maximum probe simplex radius. */
        double minStep = 0.25*tol;         /* Stop when {x} moves less than this. */
        sve_minn_iterate
          ( nv, 
            &f2_for_sve, NULL, NULL, 
            z, &Fz,
            dir, ctr, dMax, dBox, rIni, rMin, rMax, minStep,
            maxIters,
            debug
          );

        /* Return the optimal vector: */
        vars_to_points(z, p);
        (*f2p) = Fz;

        /* Local implementations: */

        void points_to_vars(r2_t q[], double y[])
          { int32_t k = 0;
            for (int32_t i = 0; i < ni; i++)
              { for (int32_t j = 0; j < 2; j++)
                  { double asij = astp[i].c[j];
                    double arij = arad[i].c[j];
                    if (r2_opt_coord_is_variable(arij, asij)) 
                      { assert(arij > 0.0);
                        y[k] = (q[i].c[j] - p0[i].c[j])/arij;
                        k++;
                      }
                  }
              }
            assert(k == nv);
          }

        void vars_to_points(double y[], r2_t q[])
          { int32_t k = 0;
            for (int32_t i = 0; i < ni; i++)
              { for (int32_t j = 0; j < 2; j++)
                  { double asij = astp[i].c[j];
                    double arij = arad[i].c[j];
                    if (r2_opt_coord_is_variable(arij, asij))
                      { q[i].c[j] = p0[i].c[j] + y[k]*arij; 
                        k++;
                      }
                  }
              }
            assert(k == nv);
          }

        double f2_for_sve(int32_t nx, double x[])
          { assert(nx == nv);
            /* Convert variables {x[0..nx-1]} to displacements {p[0..ni-1]}: */
            vars_to_points(x, p);
            /* Evaluate the client function: */
            double Q2 = f2(ni, p, iscale);
            return Q2;
          }

      }
    else
      { /* Just compute the mismatch for the given vector: */
        (*f2p) = f2(ni, p, iscale);
      }
    
    if (debug) { fprintf(stderr, "<< leave %s <<\n", __FUNCTION__); }
    return;
  }

void r2_opt_multi_scale
  ( int32_t ni,                  /* Number of points to optimize. */
    r2_opt_goal_func_t *f2,  /* Function that evaluates the goal function. */
    bool_t quadopt,          /* Use quadratic optimization? */
    r2_t arad[],             /* Max coordinate adjustment for each point. */
    r2_t astp[],             /* Adjustment step or desired precision for each point. */
    r2_t p[],                /* (IN) Initial guesses; (OUT) Best solution. */
    double *f2p,             /* (OUT) Goal function value at the best solution. */
    bool_t debug              /* Prints diagnostics if true. */
  )
  {
    if (debug) { fprintf(stderr, ">> enter %s >>\n", __FUNCTION__); }
    /* Find the num of variables {nv} and max relative search radius {umax.c[j]} on each axis: */
    r2_t umax = (r2_t){{ -INF, -INF }};
    int32_t nv = 0;
    for (int32_t i = 0; i < ni; i++)
      { for (int32_t j = 0; j < 2; j++)
          { double arij = arad[i].c[j];
            double asij = astp[i].c[j];
            if (r2_opt_coord_is_variable(arij, asij))
              { nv++;
                assert(asij > 0.0);
                double uij = arij/asij; /* Rel uncertainty this coordinate. */
                if (uij > umax.c[j]) { umax.c[j] = uij; }
              }
          }
      }
      
    if (nv == 0)
      { /* Nothing to do: */ 
        i2_t zscale = (i2_t){{ 0, 0 }};
        (*f2p) = f2(ni, p, zscale);
      }
    else
      { /* Save the initial guess as it defines the center of the domain: */
        r2_t p0[ni];
        for (int32_t i = 0; i < ni; i++) { p0[i] = p[i]; }
      
        /* Compute the initial scale {mscale}: */
        i2_t mscale = (i2_t){{ 0, 0 }};
        while ((umax.c[0] >= 0.5) || (umax.c[1] >= 0.5)) 
          { for (int32_t j = 0; j < 2; j++) 
              { if (umax.c[j] >= 0.5) 
                  { mscale.c[j]++;
                    umax.c[j] = umax.c[j]/2;
                  }
              }
          }
          
        /* Compute the search radius {rad} and step {stp} for scale {mscale}: */
        r2_t rad[ni];  /* Search radius at current scale. */
        r2_t stp[ni];  /* Search step at current scale (if {quadopt} is false). */
        
        auto void initialize_stp(void);
          /* Sets {stp} of variable coordiantes to be a little more than
            {arad}, so that the search at the max scale assumes that the
            search interval is the whole of the original interval. */

        auto void update_rad_stp(i2_t cscale);
          /* Assumes that the current candidate {p} for scale {cscale}
            has uncertainty {±stp}. Sets {rad} to the 
            proper search radius for that scale.
            Then sets {stp} to the proper search step;
            which is {0.5*rad}, or {astp} if the
            scale is zero and {astp} is smaller. */

        /* Initialize {stp[i].c[j]} to be a bit more than {arad[i]c[j]}: */
        initialize_stp();
        
        /* Now solve the problem at increasing scales: */
        i2_t iscale = mscale;
        while(TRUE)
          { 
            /* Choose the proper search radius {rad} and step {stp} for this scale: */
            update_rad_stp(iscale);
          
            /* Solve the problem at scale {iscale}: */
            if (quadopt)
              { r2_opt_single_scale_quadopt(ni, iscale, f2, rad, stp, p, f2p, debug); }
            else
              { r2_opt_single_scale_enum(ni, iscale, f2, rad, stp, p, f2p, debug); }

            /* Are we done? */
            if ((iscale.c[0] == 0) && (iscale.c[1] == 0)) { break; }

            /* Reduce the lasgest elements of {iscale}: */
            int32_t xscale = (int32_t)imax(iscale.c[0], iscale.c[1]);
            for (int32_t j = 0; j < 2; j++) 
              { if (iscale.c[j] == xscale) 
                  { iscale.c[j]--; }
              }
          }

        /* INTERNAL IMPLEMENTATIONS */

        void initialize_stp(void)
          { int32_t nv1 = 0;
            for (int32_t i = 0; i < ni; i++)
              { for (int32_t j = 0; j < 2; j++) 
                  { /* Get original search radius and step: */
                    double arij = arad[i].c[j];
                    double asij = astp[i].c[j];
                    if (r2_opt_coord_is_variable(arij, asij))
                      { stp[i].c[j] = 1.0001*arij; }
                    else
                      { stp[i].c[j] = 0.0; }
                  }
              }
            affirm(nv1 == nv, "bug nv1");
          }

        void update_rad_stp(i2_t cscale)
          { int32_t nv1 = 0;
            for (int32_t i = 0; i < ni; i++)
              { for (int32_t j = 0; j < 2; j++) 
                  { /* Get original search radius and step: */
                    double arij = arad[i].c[j];
                    double asij = astp[i].c[j];
                    if (r2_opt_coord_is_variable(arij, asij))
                      { /* The search radius in principle is the previous step: */
                        assert(stp[i].c[j] > 0.0);
                        /* Compute the new search region, clip to the original search region: */
                        double pij = p[i].c[j]; /* Current candidate. */
                        double p0ij = p0[i].c[j]; /* Initial guess. */
                        double pij_min = fmax(pij-stp[i].c[j], p0ij-arij);
                        double pij_max = fmin(pij+stp[i].c[j], p0ij+arij);
                        assert(pij_max >= pij_min);
                        /* Set {p[i].c[j]} to the center of the uncertainty region: */
                        p[i].c[j] = (pij_min + pij_max)/2;
                        /* Set the search radius to that the uncertainty region: */
                        rad[i].c[j] = (pij_max - pij_min)/2;
                        if (iscale.c[j] == 0) 
                          { /* Use the original step at scale 0: */
                            stp[i].c[j] = asij;
                          }
                        else 
                          { /* At middle scales, use 2/3 of the radius as step: */
                            stp[i].c[j] = 2.0/3.0*rad[i].c[j];
                          }
                      }
                    else 
                      { rad[i].c[j] = 0.0; stp[i].c[j] = 0.0; }
                  }
              }
            affirm(nv1 == nv, "bug nv1");
          }
      }

    if (debug) { fprintf(stderr, "<< leave %s <<\n", __FUNCTION__); }
    return;
  }
  
bool_t r2_opt_coord_is_variable(double arij, double asij)
  {
    demand(arij >= 0.0, "invalid search radius");
    if (arij > 0.0) { demand(asij > 0, "invalid search step"); }
    return ((arij > 0.0) && (arij >= asij));  
  }

double r2_opt_rel_disp_sqr(int32_t ni, r2_t p[], r2_t q[], r2_t arad[], r2_t astp[])
  {
    double d2 = 0.0;
    for (int32_t j = 0; j < 2; j++)
      { for (int32_t i = 0; i < ni; i++)
          { double arij = arad[i].c[j];
            double asij = astp[i].c[j];
            if (r2_opt_coord_is_variable(arij, asij))
              { double dij = (p[i].c[j] - q[i].c[j])/arij;
                d2 += dij*dij;
              }
          }
      }
    return d2;
  }
    
            
    

