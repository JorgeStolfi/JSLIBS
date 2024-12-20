/* See {minn.h}. */
/* Last edited on 2024-12-05 13:10:49 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
 
#include <bool.h>
#include <rn.h>
#include <rmxn.h>
#include <rmxn_ellipsoid.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <sym_eigen.h>
#include <affirm.h>

#include <minn.h>
#include <minn_enum.h>
#include <minn_quad.h>
    
void minn_uniform
  ( uint32_t n,          /* Dimension of search space. */
    minn_goal_t *F,     /* Function to be minimized. */
    bool_t box,         /* True to search in the unit cube, false in the unit ball. */
    double atol[],      /* Desired precision along each coordinate. */
    minn_method_t meth, /* Minimizaton method do use.*/
    double v[],         /* (OUT) Minimum vector found. */
    double *Fval_P      /* (OUT) Goal function value at the minimum. */
  )
  {
    demand(n >= 0, "invalid {n}");
    switch (meth)
      {
        case minn_method_ENUM:
          minn_enum(n, F, box, atol, v, Fval_P);
          break;
          
        case minn_method_QUAD:
          /* Compute the general tolerance {tol}: */
          double tol = +INF;
          for (uint32_t i = 0;  i < n; i++) { tol = fmin(tol, atol[i]); }
          minn_quad(n, F, box, tol, v, Fval_P);
          break;
          
        default:
          assert(FALSE);
      }
  }

void minn_subspace
  ( uint32_t n,          /* Dimension of search space. */
    minn_goal_t *F,     /* Function to be minimized. */
    uint32_t d,          /* Dimension of search domain. */
    double U[],         /* Main axis directions of the search domain. */
    double urad[],      /* Radii of the search domain. */
    bool_t box,         /* True the search domain is a box, false it is an ellipsoid. */
    double utol[],      /* Desired precision along each {U} row. */
    minn_method_t meth, /* Minimizaton method do use.*/
    double v[],         /* (OUT) Minimum vector found. */
    double *Fval_P      /* (OUT) Goal function value at the minimum. */
  )
  {
    demand(n >= 0, "invalid {n}");
    demand((d >= 0) && (d <= n), "invalid {d}");
    /* Map the domain to the unit {d}-ball or {d}-cube. */
    /* Compute the tolerances {xtol[0..d-1]} in each mapped axis: */
    double xtol[d];
    for (uint32_t k = 0;  k < d; k++)
      { double urk = urad[k];
        demand(isfinite(urk) && (urk > 0), "invalid {urad[k]}");
        xtol[k] = utol[k]/urk;
      }

    double ys[d], vt[n]; /* Work vectors for {unmap,F_unit}. */
    
    auto void unmap_vec(double x[], double v[]);
      /* Unmaps {x[0..d-1]} from the unit cube/ball to a vector {v[0..n-1]}
        in the original domain {D}. */
        
    auto double F_unit(uint32_t nx, double x[]);
      /* Unmaps {x[0..d-1]} to vector {v[0..n-1]} and evaluates
        the given goal function {F} on it. Expects {nx} to be {d}. */
        
    double x[d]; /* Minumum in the unit ball/cube. */
    minn_uniform(d, &F_unit, box, xtol, meth, x, Fval_P);
    unmap_vec(x, v);
    return;
    
    /* INTERNAL IMPS */
    
    void unmap_vec(double xs[], double vs[])
      { if (U == NULL)
          { assert(d == n);
            rn_weigh(n, xs, urad, vs);
          }
        else
          { rn_weigh(n, xs, urad, ys);
            rmxn_map_row(d, n, ys, U, vs);
          }
      }
      
    double F_unit(uint32_t nx, double xt[])
      { assert(nx == d);
        unmap_vec(xt, vt);
        return F(n, vt);
      }
  }

void minn_ellipsoid_constrained
  ( uint32_t n,          /* Dimension of search space. */
    minn_goal_t *F,     /* Function to be minimized. */
    double arad[],      /* Readii of the ellipsoid. */
    uint32_t q,          /* Number of explicit constraints. */
    double A[],         /* Constraint matrix. */
    double tol,         /* Desired precision. */
    minn_method_t meth, /* Minimizaton method do use.*/
    double v[],         /* (OUT) Minimum vector found. */
    double *Fval_P      /* (OUT) Goal function value at the minimum. */
  )
  {
    bool_t debug = FALSE;
    demand(n >= 0, "invalid {n}");
    demand(q >= 0, "invalid {q}");

    /* Normalize the constraints and add the implicit ones: */
    uint32_t m;        /* Number of normalized constraints. */
    double *C = NULL; /* Normalized constraint matrix. */
    rmxn_ellipsoid_normalize_constraints(n, arad, q, A, debug, &m, &C);
    assert((m >= 0) && (m <= n));

    /* Compute the search ellipsoid {\RF(U,urad)}: */
    uint32_t d = (uint32_t)(n - m);
    double U[d*n];
    double urad[d];
    rmxn_ellipsoid_cut(n, arad, m, C, d, U, urad);
    
    /* Minimize over the ellipsoid {\RF(U,urad)}: */
    bool_t box = FALSE; /* Cutting a box is messy. */
    double utol[d];
    for (uint32_t k = 0;  k < d; k++) { utol[k] = tol; }
    minn_subspace(n, F, d, U, urad, box, utol, meth, v, Fval_P);
    free(C);
  }
