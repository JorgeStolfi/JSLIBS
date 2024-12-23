/* See pst_proc_map_fractal.h */
/* Last edited on 2024-12-22 12:38:49 by stolfi */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <values.h>

#include <r2.h>
#include <r2_extra.h>
#include <r2x2.h>

#include <pst_proc_map.h>
#include <pst_proc_map_fractal.h>

/* INTERNAL PROTOTYPES */

void pst_proc_map_fractal_find_octant(r2_t *p, double seed, r2_t *a, double *rbc, r2_t *b, double *rca, r2_t *c, double *rab);
  /* Assumes that {p} is in the signed unit square. Determines the
    octant of the square where {*p} lies, and returns the corners of
    that octant (a right isosceles triangle) in {*a,*b,*c}. 
    
    The point {*a} will always be the origin. The point {*b} will be a corner of
    the signed unit square.  The point {*c} will be on the X or Y axis, with 
    coordinate {+1} or {-1}.
    
    The procedure also chooses random-looking keys {*rbc,*rca,*rab} for
    the three sides, depending on the given {seed}, so that shared edges
    of the octants get the same keys. */

void pst_proc_map_fractal_split_edge(r2_t *a, double fa, r2_t *b, double fb, double rab, r2_t *m, double *fm);
  /* Splits the edge {ab} at the midpoint and chooses a value for it,
    using the shape key {rab}.  The midpoint is returned in {*m}, and its 
    chosen value in {*fm}.  
    
    If {a,fa} and {b,fb} are swapped, and {rab} is replaced by {1-rab},
    the chosen value {fm} does not change. */
    
/* IMPLEMENTATIONS */

void pst_proc_map_fractal_cone(r2_t *p, double eps, double seed, double *z, r2_t *dz)
  { 
    r2_t debug_point = (r2_t){{ 0.3, 0.3 }};
    double debug_rad = 0.01;
    bool_t debug = (r2_dist(p, &debug_point) <= debug_rad);
    
    auto void fract
      ( r2_t *a, double fa, double rbc,
        r2_t *b, double fb, double rca,
        r2_t *c, double fc, double rab,
        int32_t ind
      );
      /* Given a triangle {a,b,c} that contains {p} (possibly at its border),
        computes the height and derivatives at {p}. Stores the result 
        in {*z,*dz}.  
        
        Assumes that the triangle is right isoceles, entirely
        contained in one octant, with sides at 0, 45, or 90 degrees.
        Assumes also that {a,b} is the longest edge of the triangle, 
        and that at leat one vertex is inside the unit circle. 
        
        The height at each corner is {fa,fb,fc}, respectively; these
        must be positive if the vertex is inside the circle, and zero
        otherwise. For any other point inside the triangle, the function
        will be in the range spanned by {fa,fb,fc}. 
        
        The parameters {rbc,rca,rab} are random key values in {(0 _ 1)}
        that completely determine the shape of the function along the
        edge opposite to {a,b,c}, respectively.  
        
        Swapping {a,fa} with {b,fb} and replacing {rbc,rca,rab} by
        {1-rca,1-rbc,1-rab} will yield the same function values at every
        point in the triangle. Ditto for . 
        
        Uses affine interpolation of {fa,fb,fc} if the side {ab} is
        shorter than {eps}.
        
        The {ind} parameter is the indenttation for debugging printouts. */
    
    auto void fract_split
      ( r2_t *a, double fa, double rbc,
        r2_t *b, double fb, double rca,
        r2_t *c, double fc, double rab,
        int32_t ind
      );
      /* Same as {fract}, but assumes that the triangle {a,b,c} is not
        small enough to use affine interpolation. */

    /* Initialize to zero by default: */
    (*z) = 0.0;
    dz->c[0] = 0.0; dz->c[1] = 0.0;
    
    /* Grab coordinates {x,y} of {p}: */
    double x = p->c[0];
    double y = p->c[1];
    
    /* Outside the unit circle the height is zero: */
    if (x*x + y*y >= 1.0) { return; }
    
    /* Decide which octant {p} is in, and set {a0,b0,c0} and {rbc0,rca0,rab0} accordingly: */
    r2_t a0, b0, c0;
    double rbc0, rca0, rab0;
    pst_proc_map_fractal_find_octant(p, seed, &a0, &rbc0, &b0, &rca0, &c0, &rab0);
    
    /* Choose corner heights {fa0,fb0,fc0}: */
    double fa0 = 1.0; /* At origin. */
    double fb0 = 0.0; /* Outside unit circle. */
    double fc0 = 0.0; /* On edge of unit circle. */
   
    /* Call recursion to define {*z,*dz} if not zero: */
    if (debug) { fprintf(stderr, "----------------------------------------------------\n"); }
    if (debug) { fprintf(stderr, "p = ( %8.5f, %8.5f )\n", p->c[0], p->c[1]); }
    fract(&a0, fa0, rbc0, &b0, fb0, rca0, &c0, fc0, rab0, 0);
    if (debug) { fprintf(stderr, "z = %8.5f  dz = ( %8.5f, %8.5f )\n", (*z), dz->c[0], dz->c[1]); }
    if (debug) { fprintf(stderr, "----------------------------------------------------\n"); }
    
    return;
    
    void fract
      ( r2_t *a, double fa, double rbc,
        r2_t *b, double fb, double rca,
        r2_t *c, double fc, double rab,
        int32_t ind
      )
      { 
        if (debug) { fprintf(stderr, "%*s fa = %8.5f  fb = %8.5f  fc = %8.5f\n", ind, "", fa, fb, fc); }
        assert ((fa > 0) || (fb > 0) || (fc > 0));
          
        if (r2_dist_sqr(a,b) < eps*eps)
          { /* Triangle is small enough, use an affine approximation: */
            pst_proc_map_function_affine(p, a, fa, b, fb, c, fc, z, dz);
          }
        else  
          { /* Triangle is not small enough; divide and recurse. */
            fract_split(a, fa, rbc, b, fb, rca, c, fc, rab, ind);
          }
      }
        
    void fract_split
      ( r2_t *a, double fa, double rbc,
        r2_t *b, double fb, double rca,
        r2_t *c, double fc, double rab,
        int32_t ind
      )    
      {
        assert((fa > 0) || (fb > 0) || (fc > 0));
        /* Since {ab} is the longest side, split it at the midpoint {m}: */
        r2_t m;
        double fm;
        pst_proc_map_fractal_split_edge(a, fa, b, fb, rab, &m, &fm);
        
        /* Determine which side of the line {c--m} contains {p}: */
        double xcm = m.c[0] - c->c[0];
        double ycm = m.c[1] - c->c[1];
        double xcp = x - c->c[0];
        double ycp = y - c->c[1];
        double xca = a->c[0] - c->c[0];
        double yca = a->c[1] - c->c[1];
        double Dp = xcm*ycp - ycm*xcp;
        double Da = xcm*yca - ycm*xca;
        /* The recursion requires computing keys {ram,rmb,rmc,rcm} 
          for the  edges {a--m} and {m--b}, {m--c} and {c--m}. 
          
          Define a "flip" operation as exchanging {a<-->b}, {fa<-->fb},
          and {rca<-->rbc}, and complementing all three input keys {rbc,rca,rab}.
          Such a flip must not change the value at any point along those four 
          edges.  
          
          Therefore, {rmc} and {rcm} must obviously be complements of each other,
          and must be unaffected by a flip.  Also {ram} and {rmb} must be 
          swapped and complemented by a flip, and must depend only on the global
          seed and on {rab}, because it is the only key shared with the adjacent
          triangle. */
          
        /* Choose 'random' key {rmc} for edge {m} to {c}:  */
        /* Note that it is invariant under a flip. */
        double sab = 2*fabs(rab - 0.5); /* A number in {0_1} invariant under flip. */
        double rmc = 0.5*(1 + sin(33*(rbc-rca) + 27*sab + 19*seed + 0.5)); 
        if (Dp*Da >= 0)
          { /* Point {p} is on the same side as corner {a}. Recurse with {c,a,m}, if not all zero: */
            if ((fc > 0) || (fa > 0) || (fm > 0))
              { /* Compute new 'random' key for the edge {am}: */
                double ram = 0.5*(1 + sin(17*rab + 43*seed + 0.5)); 
                fract(c, fc, ram, a, fa, rmc, &m, fm, rca, ind+1);
              }
          }
        else
          { /* Point {p} is on the same side as corner {b}. Recurse with {b,c,m}, if not all zero: */
            if ((fb > 0) || (fc > 0) || (fm > 0))
              { /* Compute new 'random' key for the edge {mb}: */
                double rmb = 0.5*(1 - sin(17*(1-rab) + 43*seed + 0.5)); 
                fract(b, fb, 1-rmc, c, fc, rmb, &m, fm, rbc, ind+1);
              }
          }
      }
  }

void pst_proc_map_fractal_find_octant(r2_t *p, double seed, r2_t *a, double *rbc, r2_t *b, double *rca, r2_t *c, double *rab)
  { 
    (*a) = (r2_t){{ 0.0, 0.0 }}; /* Always at the origin. */
    double x = p->c[0];
    double y = p->c[1];
    int32_t koc; /* Index of octant. Clockwise, {0,4,5,2,3,7,6,1}, with 0 being {x>y>0}. */
    if (fabs(x) >= fabs(y))
      { /* Octants 0..3: */
        c->c[1] = 0.0; koc = 0;
        if (x >= 0)
          { b->c[0] = c->c[0] = +1.0; }
        else
          { b->c[0] = c->c[0] = -1.0; koc = koc + 2; }
        if (y >= 0)
          { b->c[1] = +1.0; }
        else
          { b->c[1] = -1.0; koc = koc + 1; }
      }
    else
      { /* Octants 4..7: */
        c->c[0] = 0.0; koc = 4;
        if (y >= 0)
          { b->c[1] = c->c[1] = +1.0; }
        else
          { b->c[1] = c->c[1] = -1.0; koc = koc + 2; }
        if (x >= 0)
          { b->c[0] = +1.0; }
        else
          { b->c[0] = -1.0; koc = koc + 1; }
      }

    /* Choose {*rbc,*rca,*rab}: */
    int32_t kaa = koc; /* Edge {b0,c0} is different for each octant. */
    int32_t kbb = (koc >> 1); /* Edge {a0,c0} is shared 0=0+1, 1=2+3, 2=4+5, 3=6+7. */
    int32_t kcc = (koc == 5 ? 2 : (koc == 6 ? 1 : koc & 3)); /* Edge {a0,c0} is shared 0=0+4, 1=1+6, 2=2+5, 3=3+7. */
    (*rbc) = 0.5*(1 + sin(kaa + 0.5)); 
    (*rca) = 0.5*(1 + sin(kbb + 8.5)); 
    (*rab) = 0.5*(1 + sin(kcc + 16.5));
  }

void pst_proc_map_fractal_split_edge(r2_t *a, double fa, r2_t *b, double fb, double rab, r2_t *m, double *fm)
  { 
    double tm = 0.5; /* Fractional position of split point {m} on {a,b}. */
    r2_mix(1-tm, a, tm, b, m); /* The split point on {a,b}. */

    assert((fa >= 0) && (fb >= 0));
    assert((0 <= rab) && (rab <= 1));
    
    if ((fa == 0) && (fb == 0))
      { /* Both points outside. Assume that whole edge is outside: */
        (*fm) = 0.0;
        assert(((*fm) > 0) || ((fa <= 0) && (fb <= 0)));
      }
    else
      { /* Estimate the part {a1,b1} of {a,b} inside the circle: */
        double ta, tb; /* Fractional positions of {a1,b1} in {a,b}. */
        double eps = 1.0e-8; /* A fudge factor to account for rounding errors. */
        r2_clip_seg_to_unit_disk(a, b, &ta, &tb);
        assert((0 <= ta) && (ta <= tb) && (tb <= 1.0)); /* Otherwise {fa} and {fb} should be zero. */
        assert((fa == 0) || (ta < 0 + eps));
        assert((fb == 0) || (tb > 1 - eps));
        /* Choose the value {fm} at {m}: */
        if ((tm <= ta + eps) || (tm >= tb - eps))
          { /* Seems that {m} is practically outside the unit circle: */
            (*fm) = 0.0;
          }
        else
          { /* The point {m} is inside the circle. */
            /* Estimate the range of values {fa1,fb1} for {fm}: */
            double fa1, fb1;
            double th = (ta + tb)/2; /* Fractional position where the range is {[fa _ fb]}. */
            if (fabs(tm - th) < eps)
              { /* The range is just {[fa_fb]}: */
                fa1 = fa; fb1 = fb;
              }
            else
              { double tr = (tb - ta)/2 + 2*eps; /* Fractional half-length of the clipped segment. */
                if (tm < th)
                  { /* Point {m} is closer to {a1} than to {b1}: */
                    double s = (th - tm + eps)/tr;
                    assert((0 < s) && (s < 1));
                    fa1 = fa;
                    fb1 = fb + s*(fa - fb);
                  }
                else
                  { /* Point {m} is closer to {b1} than to {a1}: */
                    double s = (tm - th + eps)/tr;
                    assert((0 < s) && (s < 1));
                    fa1 = fa + s*(fb - fa);
                    fb1 = fb;
                  }
              }
            /* Compute the {fa1,fb1} interpolation factor {q}: */
            double q = 2*rab - 1; /* Maps {rab} from {[0_1]} to {[-1_+1]}. */
            q = 0.5*q*q*q; /* To bias towards 0. */
            q = (q + 1)/2; /* Maps back to {[0_1]}. */
            q = eps + (1-2*eps)*q; /* Keep away from extremes. */
            (*fm) = fa1 + q*(fb1 - fa1);
          }
      }
  }
