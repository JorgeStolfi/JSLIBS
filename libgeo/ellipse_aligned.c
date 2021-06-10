/* See ellipse_aligned.h */
/* Last edited on 2021-06-09 19:48:23 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <r2x2.h>
#include <interval.h>
#include <affirm.h>

#include <ellipse_aligned.h>

#define Pr fprintf
#define Er stderr
  
void ellipse_aligned_bbox(double rx, double ry, interval_t bbox[])
  { LO(bbox[0]) = - rx;
    HI(bbox[0]) = + rx;
    LO(bbox[1]) = - ry;
    HI(bbox[1]) = + ry;
  }

void ellipse_aligned_int_bbox
  ( double cx, double cy,  /* Center of ellipse. */
    double rx, double ry, 
    double mrg, /* Extra margin. */
    int32_t *xLoP,  /* (OUT) Min X of clip area. */
    int32_t *xHiP,  /* (OUT) Max X of clip area. */
    int32_t *yLoP,  /* (OUT) Min Y of clip area. */
    int32_t *yHiP   /* (OUT) Max Y of clip area. */
  )
  {
    demand(mrg >= 0, "the extra margin {mrg} must be non-negative");
    (*xLoP) = (int32_t)floor(cx - rx - mrg);
    (*xHiP) = (int32_t)ceil (cx + rx + mrg);
    (*yLoP) = (int32_t)floor(cy - ry - mrg);
    (*yHiP) = (int32_t)ceil (cy + ry + mrg);
  }

double ellipse_aligned_position(double xp, double yp, double rx, double ry, double *cosP, double *sinP)
  {
    demand(rx >= 0, "invalid rx"); 
    demand(ry >= 0, "invalid ry");
    bool_t get_cs = ((cosP != NULL) || (sinP != NULL));
    
    if (rx < ry)
      { /* Ellipse is mostly vertical; swap axes: */
        return ellipse_aligned_position(yp, xp, ry, rx, sinP, cosP);
      }
    
    double cost = NAN;
    double sint = NAN;
    double pos = NAN;
    
    if (rx == 0)
      { /* Ellipse is a point: */
        assert(ry == 0);
        if ((xp == 0) && (yp == 0))
          { /* Point is at origin, assume it is on the boundary: */
            pos = 1.0;
          }
        else
          { if (get_cs) 
              { double rp = hypot(xp,yp); 
                cost = xp/rp; sint = yp/rp;
              }
            pos = +INF;
          }
      }
    else if (ry == 0)
      { /* Ellipse is a horizontal line segment: */
        if (yp == 0)
          { /* Point is on the major axis: */
            if (fabs(xp) <= rx)
              { /* Assume point is on the boundary: */
                cost = xp/rx; sint = NAN;
                pos = 1.0;
              }
            else
              { /* Point is outside assume it is on the X axis: */
                cost = (xp > 0 ? +1.0 : -1.0); sint = 0.0;
                pos = xp/fabs(rx);
              }
          }
        else 
          { /* Point is not on the major axis: */
            cost = 0.0; sint = (yp > 0 ? +1.0 : -1.0);
            pos = +INF;
          }
      }
    else
      { /* Compute the canonical coords {xe,ye} of {p}: */
        double xe = xp/rx;
        double ye = yp/ry;
        /* Compute the canonical radius {rp} of {p}: */
        double rp = hypot(xe, ye);
        cost = xe/rp;
        sint = ye/rp;
        pos = rp;
      }
    if (cosP != NULL) { (*cosP) = cost; }
    if (sinP != NULL) { (*sinP) = sint; }
    return pos;
  }
  
double ellipse_aligned_nearest_point(double rx, double ry, double xp, double yp, double *xqP, double *yqP)
  {
    bool_t debug = FALSE;
    /* !!! Delete excess debugging printouts. !!! */

    if (debug) { Pr(Er, "  rx = %13.8f  ry = %13.8f\n", rx, ry); }
    if (debug) { Pr(Er, "  p rel = [ %13.8f %13.8f ]\n", xp, yp); }

    assert(rx >= 0); assert(ry >= 0);
    
    if (rx < ry)
      { /* Swap axes and try again: */
        return ellipse_aligned_nearest_point(ry, rx, yp, xp, yqP, xqP);
      }
    
    double xq = NAN;
    double yq = NAN;
    double dist = NAN;
    if (rx == 0)
      { assert(ry == 0); /* Since {0 <= ry <= rx}. */
        /* Ellipse is essentially a point: */
        xq = 0; yq = 0;
        dist = hypot(xp - xq, yp - yq);
      }
    else if ((rx - ry) == 0)
      { /* Ellipse is essentially a circle: */
        if ((xp == 0) && (yp == 0))
          { /* {p} is at the center: */
            xq = rx; yq = 0; 
            dist = -rx;
          }
        else
          { /* {p} is eccentric: */
            double rp = hypot(xp, yp);
            double scale = rx/rp;
            xq = xp * scale;
            yq = yp * scale;
            dist = rp - rx;
          }
      }
    else if (xp == 0)
      { /* Point {p} is on the Y axis; nearest pt is {sgn(yp)*v*ry}. */
        xq = 0; yq = (yp >= 0 ? ry : -ry);
        dist = fabs(yp) - ry;
      }
    else
      { /* Compute the aspect ratio {S} of the ellipse: */  
        assert(rx > 0);
        double S = ry/rx;

        if (debug) { Pr(Er, "  S     = %13.8f\n", S); }

        if (S == 0)
          { /* Ellipse is essentially a segment in the {u} direction: */
            double dpq;
            if (xp > +rx)
              { xq = +rx; dpq = hypot(xp - xq, yp); }
            else if (xp < -rx) 
              { xq = -rx; dpq = hypot(xp - xq, yp); }
            else
              { xq = xp; dpq = fabs(yp); }
            yq = 0;
            dist = dpq;
          }

        /* Compute the slope {T} of {p}: */
        double T = yp/xp;

        if (debug) { Pr(Er, "  T     = %13.8f\n", T); }

        if (T == 0)
          { /* Point {p} is essentially on the major axis. */
            double e = 1 - S*S;
            xq = xp/e;
            if (debug) { Pr(Er, "  xq    = %13.8f\n", xq); }
            if (fabs(xq) >= rx)
              { /* Point {p} is beyond the focus. */
                /* The nearest point is {ctr ± rx*u}: */
                xq = (xp >= 0 ? +rx : -rx); yq = 0;
                dist = fabs(xp) - rx;    
              }
            else
              { /* Point {p} is between the foci. */
                double cq = xq/rx;
                assert(cq <= 1);
                yq = ry*sqrt(1 - cq*cq);
                if (debug) { Pr(Er, "  yq    = %13.8f\n", yq); }
                dist = - hypot(xq - xp, yq);
              }
          }
        else
          { 
            /* Ellipse has non-empty interior, and {p} is not on the axes. */

            /* Compute the adimensional parameters {A,B} for 1st quadrant: */
            assert(S > 0);
            assert(S <= 1);
            assert(isfinite(T));
            double A = S*fabs(T); 
            assert(! isnan(A));

            double d2 = ((rx - ry)/rx)*(rx + ry);
            assert(! isnan(d2));
            assert(isfinite(xp));
            assert(xp != 0);
            double B = d2/fabs(xp);
            assert(! isnan(B));

            if (debug) { Pr(Er, "  A     = %13.8f\n", A); }
            if (debug) { Pr(Er, "  B     = %13.8f\n", B); }

            /* Solve the characteristic polynomial for {t}: */
            double t = ellipse_aligned_compute_t(A, B);

            if (debug) { Pr(Er, "  t sol = %13.8f\n", t); }

            assert(t >= 0.0);
            assert(t <= 1.0);

            /* Compute the cosine and sine {ct,st} of the angular argument of {q}: */
            double dt = 1 + t*t; 
            double ct = (1 - t*t)/dt; 
            double st = 2*t/dt;

            /* Compute the {u,v} coords of {q}: */
            xq = (xp >= 0 ? +1 : -1)*rx*ct; 
            yq = (yp >= 0 ? +1 : -1)*ry*st; 

            /* Compute {q} and the distance: */
            double dpq = hypot(xp - xq, yp - yq);

            /* Return with correct sign: */
            if ((xp > rx) || (yp > ry)) 
              { /* {p} is definitely outside: */ dist = dpq; }
            else
              { assert(rx > 0);
                assert(ry > 0);

                /* Compute the canonical coords {xe,ye} of {p}: */
                double xe = xp/rx;
                double ye = yp/ry;
                dist = (xe*xe + ye*ye < 1.0 ? -dpq : +dpq);
              }
          }
      }
    if (xqP != NULL) { (*xqP) = xq; }
    if (yqP != NULL) { (*yqP) = yq; }
    return dist;
  }
  
double ellipse_aligned_compute_t(double A, double B)
  {
    bool_t debug = FALSE;
    /* !!! Delete excess debugging printouts. !!! */

    /* The polynomial and its derivatives: */
    /* {P(t)    = A*(1+t**2)*(1-t**2) - 2*t*((B+1)*t**2 - (B-1))} */
    /* {P'(t)   = -2*((2*A*t + 3*(B+1))*t**2 - (B-1))} */
    /* {P''(t)  = -12*t*(A*t + (B+1))} */
    /* {P'''(t) = -12*(2*A*t + (B+1))} */
    
    demand(A >= 0, "{A} must be non-negative");
    demand(B >= 0, "{B} must be non-negative");
    
    if (A == 0)
      { /* {P(t)} reduces to {2*t*((B+1)*t**2 - (B-1))}. */
        /* The roots in {[0:1]} are {t=0} and {t=sqrt((B-1)/(B+1))}. */
        if (B-1 <= 0)
          { return 0; }
        else
          { return sqrt((B-1)/(B+1)); }
      }
    
    /* Since {P(0) == A > 0} and {P(1) == -4 < 0}, there is a root in {[0:1]}. */
    /* Solve {P(t)==0} by Newton-Raphson. */
    /* Since {P''(t) <= 0} in {[0:1]}, we can use {t0==1} as starting guess. */
    /* !!! Should use a quadratic approximation to get a good {t0}. !!! */
    /* !!! Or even use a quadratic Newton, since {P'' < 0} and {P''' < 0}. !!! */
    
    double t = 1.0;
    double Bp = B + 1;
    double Bm = B - 1;
    int32_t nIts = 0;  /* Counts iterations. */
    if (debug) { Pr(Er, "  t[%d] = %13.8f\n", nIts, t); }
    while (TRUE)
      { /* Compute {P(t)}: */
        double t2 = t*t;
        double Pt = A*(1+t2)*(1-t2) - 2*t*(Bp*t2 - Bm);
        if (debug) { Pr(Er, "  P[%d] = %13.8f\n", nIts, Pt); }
        if (Pt >= 0) { /* We must have reached or passed the root: */ break; }
        
        /* Compute {P'(t)}: */
        double Dt = -2*((2*A*t + 3*Bp)*t2 - Bm);
        if (debug) { Pr(Er, "  D[%d] = %13.8f\n", nIts, Dt); }
        if (Dt >= 0) { /* We must have reached or passed the root: */ break; }
        
        /* Compute the new {t}: */
        double dt = -Pt/Dt;
        if (debug) { Pr(Er, "  d[%d] = %13.8f\n", nIts, dt); }
        
        assert(dt <= 0);
        double ot = t;
        t += dt;
        nIts++;
        if (debug) { Pr(Er, "  t[%d] = %13.8f\n", nIts, t); }
        if (t < 0) { /* Must be a very small positive root: */ t = 0; break; }
        if (t >= ot) { /* We cannot progress any further: */ break; }
        demand(nIts <= 1000, "did not converge in 1000 iterations");
        if (debug) { Pr(Er, "\n"); }
      }
    if (debug) { Pr(Er, "  t fin = %13.8f\n", t); }
    if (nIts > 25) { Pr(Er, "%s required %d iterations\n", __FUNCTION__, nIts); }
    if (debug) { Pr(Er, "\n"); }
    return t;    
  }
  
bool_t ellipse_aligned_inside(double rx, double ry, double xp, double yp)
  {
    assert(rx >= 0); assert(ry >= 0);
    bool_t debug = FALSE;

    if (debug) { Pr(Er, "  rx = %13.8f  ry = %13.8f\n", rx, ry); }
    if (debug) { Pr(Er, "  p rel = [ %13.8f %13.8f ]\n", xp, yp); }

    if ((rx == 0) || (ry == 0))
      { /* Ellipse is a point or line segment: */
        return FALSE;
      }
    
    /* Test the implicit equation: */
    double xd = xp/rx;
    double yd = yp/ry;
    return xd*xd + yd*yd <= 1.0;
  }

double ellipse_aligned_border_position(double rx, double ry, double hwd, double xp, double yp)
  {
    bool_t debug = FALSE;
    assert(rx >= ry); assert(ry >= 0);
    
    if (debug) { Pr(Er, "  rx = %13.8f  ry = %13.8f\n", rx, ry); }
    if (debug) { Pr(Er, "  p rel = [ %13.8f %13.8f ]\n", xp, yp); }
    
    if (rx == ry)
      { /* Ellipse is a circle: */
        return fmin(+1.0, fmax(-1.0, (hypot(xp,yp) - rx)/hwd));
      }
    
    /* Check bounding rectangle: */
    if ((fabs(xp) >= rx + hwd) || (fabs(yp) >= ry + hwd)) { return +1.0; }
    
    if (ry == 0)
      { /* Ellipse is a line segment: */
        if (fabs(xp) <= rx)
          { /* Straight part of stroke: */
            return fmin(1.0, fabs(yp)/hwd);
          }
        else
          { /* Round end caps: */
            double xq = (xp > 0 ? +rx : -rx);
            return fmin(1.0, hypot(xp - xq, yp)/hwd);
          }
      }

    /* Let {Q} be the upright square centered at {p} with side {2*hwd}. */
    /* Test the outermost corner of {Q}: */
    double aoc = (fabs(xp) + hwd)/rx;
    double boc = (fabs(yp) + hwd)/ry;
    if (aoc*aoc + boc*boc < 1.0) 
      { /* {Q} is inside the ellipse, so {p} is well inside: */
        return -1.0;
      }
    /* Test the innermost corner of {Q}, clipped to the quadrant: */
    double aic = fmax(0.0, (fabs(xp) - hwd))/rx;
    double bic = fmax(0.0, (fabs(yp) - hwd))/ry;
    if (aic*aic + bic*bic > 1.0) 
      { /* {Q} is outside the ellipse, so {p} is well outside: */
        return +1.0;
      }
    /* We are close enough to the boundary. */
    /* Cmpute the exact signe distance {d} to the boundary: */
    double d = ellipse_aligned_nearest_point(rx, ry, xp, yp, NULL, NULL);
    return fmax(-1.0, fmin(+1.0, d/hwd));
  }

void ellipse_aligned_print(FILE *wr, double rx, double ry, char *fmt)
  { fprintf(wr, "{ rx: ");
    fprintf(wr, fmt, rx);
    fprintf(wr, " ry: ");
    fprintf(wr, fmt, ry);
    fprintf(wr, " }");
    fflush(wr);
  }
