/* See rn_classif_test.h */
/* Last edited on 2024-12-05 10:24:29 by stolfi */

#include <math.h>
#include <stdio.h>
#include <assert.h>

#include <r2.h>
#include <rn.h>
#include <jsrandom.h>
#include <affirm.h>
#include <jsmath.h>

#include <rn_classif_test.h>

double rn_classif_test_dist_ellipse(double x, double y, double A, double B, double eps);
  /* Approximate signed distance from {(x,y)} to the ellipse with radii {A,B}.
    Barely accurate only within {eps} from the ellipse where {eps << min(A,B)}.
    In particular, returns {±eps} when the distance is about {eps}. */
  
double rn_classif_test_dist_ellipse(double x, double y, double A, double B, double eps)
  {
    assert(A > 2*eps);
    assert(B > 2*eps);
    double ro = hypot(x/(A+eps), y/(B+eps));
    double ri = hypot(x/(A-eps), y/(B-eps)) + 1e-100;
    double d = eps*(2*ri*ro - ri - ro)/(ri - ro);
    return d;
  }

void rn_classif_test_check_saturn(int *NA, int *NC)
  { if ((*NA) == 0) { (*NA) = 2; }
    if ((*NC) == 0) { (*NC) = 2; }
    demand((*NA) == 2, "invalid number of attributes");
    demand((*NC) == 2, "invalid number of classes");
  }

#define saturn_A1 0.875
#define saturn_B1 0.375
#define saturn_A2 0.500
#define saturn_B2 0.500

int rn_classif_test_label_saturn(int NA, int NC, double p[])
  { assert(NA == 2);
    assert(NC == 2);
    /* Get coordinates: */
    double x = p[0];
    double y = p[1];
    /* Rescale {p} so that each ellipse becomes the unit circle: */
    double eps = 0.01; /* Boundary tolerance. */
    double e1 = fabs(rn_classif_test_dist_ellipse(x, y, saturn_A1, saturn_B1, eps));
    double e2 = fabs(rn_classif_test_dist_ellipse(x, y, saturn_A2, saturn_B2, eps));
    int cl;
    if (e1 < e2)
      { cl = (e1 > eps ? 0 : 1); }
    else
      { cl = (e2 > eps ? 0 : 2); }
    return cl;
  }

void rn_classif_test_throw_saturn(int i, int NA, int NC, double p[], int *class)
  { assert(NA == 2);
    assert(NC == 2);
    assert(i >= 0);
    int M = NC; /* Sampling period. */
    int cl = (i % M) + 1;
    double A = (cl == 1 ? saturn_A1 : saturn_A2);
    double B = (cl == 1 ? saturn_B1 : saturn_B2);
    double t = drandom();
    t = t + 0.25*sin(4*M_PI*t);
    double x = A*cos(2*M_PI*t);
    double y = B*sin(2*M_PI*t);
    p[0] = x;
    p[1] = y;
    (*class) = cl;
  }

void rn_classif_test_check_petals(int *NA, int *NC)
  { if ((*NA) == 0) { (*NA) = 2; }
    if ((*NC) == 0) { (*NC) = 4; }
    demand((*NA) == 2, "invalid number of attributes");
    demand((*NC) == 4, "invalid number of classes");
  }

#define petals_Rmax 0.80
#define petals_Rmin 0.70
#define petals_H    0.40
#define petals_Rint 0.01

int rn_classif_test_label_petals(int NA, int NC, double p[])
  { assert(NA == 2);
    assert(NC == 4);
    /* Get coordinates: */
    double x = p[0];
    double y = p[1];
    /* Rotate petal to X-axis, assign tentative class: */
    int cl = 1;
    if (y > +x) { double t = x; x = +y; y = t; cl = 2; }
    if (y < -x) { double t = x; x = -y; y = t; cl = 5 - cl; }
    assert((cl >= 1) && (cl <= 4));
    /* Shift point left by {Rint}: */
    x = x - petals_Rint;
    double eps = 0.01; /* Boundary tolerance. */
    if (x > 0) 
      { /* Expand {y} by {1/(H*x)}: */
        y = y/(petals_H*x); 
        /* Distance from origin should be in range {Rmin,Rmax}: */
        double Reps = eps/(petals_H*x); /* Boundary tolerance. */
        double r = hypot(x,y);
        if ((r > petals_Rmax + Reps) || (r < petals_Rmin - Reps)) { cl = 0; }
      }
    else
      { /* Distance from X-axis must be small: */
        if (fabs(y) > eps) { cl = 0; } 
      }
    return cl;
  }

void rn_classif_test_throw_petals(int i, int NA, int NC, double p[], int *class)
  { assert(NA == 2);
    assert(NC == 4);
    assert(i >= 0);
    int M = NC; /* Sampling period. */
    int cl = (i % M) + 1;
    /* Throw a point in the half-ring with radii {Rmin,Rmax} and positive {x}: */
    double t = M_PI*drandom();
    double R = petals_Rmin + (petals_Rmax - petals_Rmin)*drandom();
    double xt = R*sin(t);
    double yt = R*cos(t);
    /* Squash {y} by {H*x}: */
    yt = petals_H*yt*xt;
    /* Shift petal right by {Rint}: */
    xt = xt + petals_Rint;
    /* Rotate according to class: */
    double s = 0.5*M_PI*(cl-1);
    double x = xt*cos(s) - yt*sin(s);
    double y = xt*sin(s) + yt*cos(s);
    p[0] = x;
    p[1] = y;
    (*class) = cl;
  }

void rn_classif_test_check_vessel(int *NA, int *NC)
  { if ((*NA) == 0) { (*NA) = 2; }
    if ((*NC) == 0) { (*NC) = 3; }
    demand((*NA) == 2, "invalid number of attributes");
    demand((*NC) == 3, "invalid number of classes");
  }

#define vessel_A    1.00
#define vessel_B    0.70
#define vessel_xmin 0.10
#define vessel_xmax 0.50
#define vessel_ymax 0.40

int rn_classif_test_label_vessel(int NA, int NC, double p[])
  {
    assert(NA == 2);
    assert(NC == 3);
    double xc = (vessel_xmax + vessel_xmin)/2;
    double yc = 0;
    double xh = (vessel_xmax - vessel_xmin)/2;
    double yh = vessel_ymax;
    /* Get coordinates: */
    double x = p[0];
    double y = p[1];
    /* Check against ellipse: */
    double eps = 0.01; /* Boundary tolerance. */
    double e1 = fabs(rn_classif_test_dist_ellipse(x, y, vessel_A, vessel_B, eps));
    /* Check against rectangle: */
    double xr = fabs(x) - xc;
    double yr = y - yc;
    double ex23 = fabs(xr) - xh;
    double ey23 = fabs(yr) - yh;
    double e23 = fmax(ex23, ey23);
    int cl;
    if (e1 < e23)
      { cl = (e1 > eps ? 0 : 1); }
    else
      { cl = (e23 > eps ? 0 : (x < 0 ? 2 : 3)); }
    return cl;
  }

void rn_classif_test_throw_vessel(int i, int NA, int NC, double p[], int *class)
  { assert(NA == 2);
    assert(NC == 3);
    double xc = (vessel_xmax + vessel_xmin)/2;
    double yc = 0;
    double xh = (vessel_xmax - vessel_xmin)/2;
    double yh = vessel_ymax;
    assert(i >= 0);
    int M = NC; /* Sampling period. */
    int cl = (i % M) + 1;
    double x, y;
    if (cl == 1)
      { double t = drandom();
        t = t + 0.25*sin(4*M_PI*t);
        x = vessel_A*cos(2*M_PI*t);
        y = vessel_B*sin(2*M_PI*t);
      }
    else
      { double xt = 2*drandom() - 1; xt = (3 - xt*xt)*xt/2.0;
        double yt = 2*drandom() - 1; yt = (3 - yt*yt)*yt/2.0;
        x = xc + xt*xh;
        y = yc + yh*yt;
        int xsgn = 2*cl - 5;
        x = x*xsgn;
      }
    p[0] = x;
    p[1] = y;
    (*class) = cl;
  }

void rn_classif_test_check_mballs(int *NA, int *NC)
  { if ((*NA) == 0) { (*NA) = 2; }
    if ((*NC) == 0) { (*NC) = 2; }
    demand(((*NA) >= 1) && ((*NA) <= 8), "invalid number of attributes");
    demand((*NC) >= 2, "invalid number of classes");
    int64_t NT = ipow((*NC)-1,(*NA)); /* Total numberof balls/cubelets. */
    demand(NT < 256, "too many domain components"); 
  }

#define mballs_r_rel (2.0/3.0)

int rn_classif_test_label_mballs(int NA, int NC, double p[])
  {
    assert((NA >= 1) && (NA <= 8));
    assert(NC >= 2);
    int NB = NC - 1;          /* Number of balls/cubelets along each axis, and of ball domains. */
    
    /* Compute cubelet and ball geometry: */
    double R_abs = 1.0/NB;       /* Inradius of each cubelet. */
    double S = 2.0*R_abs;        /* Spacing between ball centers, and side of each cubelet. */
    double r_rel = mballs_r_rel; /* Radius of ball relative to inradius of cubelet. */

    /* Get distance from {p} to the cube {U^NA}: */
    int k;
    double eps = 0.01;
    double du = 0.0;
    for (k = 0; k < NA; k++) { du = fmax(du, fabs(p[k]) - 1); }
    int cl;
    if (du > eps)
      { /* Outside {U^NA}: */
        cl = 0;
      }
    else
      { /* Find the cubelet indices {ix[0..NA-1]} and scale cubelet to {U^NA}: */
        int ix[NA]; 
        double u[NA];  /* Point {p} relative to cubelet. */
        for (k = 0; k < NA; k++) 
          { /* Get cubelet index along axis {k} */
            double uk = (p[k] + 1)/S;
            ix[k] = (int)floor(uk);
            /* Fix roundoff errors: */
            if (ix[k] < 0) { ix[k] = 0; }
            if (ix[k] >= NB) { ix[k] = NB-1; }
            /* Scale {p} from absolute to cubelet-relative coords: */
            u[k] = 2*(uk - ix[k]) - 1;
          } 
        /* Check against ball: */
        double r = rn_norm(NA, u);
        if (r > r_rel) 
          { /* Background class: */
            cl = 1;
          }
        else
          { /* Get class of ball: */
            /* !!! Should use a more irregular assignment. !!! */
            int ixsum = 0;
            for (k = 0; k < NA; k++) { ixsum += ix[k]; } 
            cl = 2 + (ixsum % NB);
          }
      }
    return cl;
  }

void rn_classif_test_throw_mballs(int i, int NA, int NC, double p[], int *class)
  { assert((NA >= 1) && (NA <= 8));
    assert(NC >= 2);
    int NB = NC - 1;          /* Number of balls/cubelets along each axis, and of ball domains. */
    
    /* Compute cubelet and ball geometry: */
    double R_abs = 1.0/NB;       /* Inradius of each cubelet. */
    double S = 2.0*R_abs;        /* Spacing between ball centers, and side of each cubelet. */
    double r_rel = mballs_r_rel; /* Radius of ball relative to inradius of cubelet. */
    double r_abs = r_rel*R_abs;  /* Radius of each ball. */
    
    assert(i >= 0);
    /* In each sampling subperiod we generate one point in one ball and {NH} points in the rest of the cubelet: */
    /* Compute {NH} so that the sampling density is about the same: */
    double H = pow(S,NA);                  /* Measure of one cubelet. */
    double B = pow(r_abs,NA)*ball_vol(NA); /* Measure of one ball. */
    assert(B <= 0.5001*H);
    int NH = (int)((H - B)/B + 0.5); 
    if (NH < 1) { NH = 1; }
    int NS = NH + 1; /* Samples in each cubelet (sampling subperiod). */
    int M = NS*(int)ipow(NB,NA); /* Sampling period: */

    /* Determine the cubelet {icub} and the sample index {ismp} in cubelet: */
    int phase = i % M; /* Phase in period. */
    int icub = phase / NS;  /* Cubelet index, in {0..NB^NA-1}. */
    int ismp = phase % NS;  /* Sample index in cubelet. */

    /* Break the cubelet index {icub} into separate indices {ix[0..NA-1]}: */
    int ix[NA];
    int k;
    int tmp = icub;
    for (k = 0; k < NA; k++) { ix[k] = tmp % NB; tmp = tmp/NB; } 

    /* Compute the class {cl} of the sample: */
    /* !!! Should use a more irregular assignment. !!! */
    int ixsum = 0;
    for (k = 0; k < NA; k++) { ixsum += ix[k]; } 
    int cl = (ismp == 0 ? 1 : 2 + (ixsum % NB));

    /* Pick a point in {U^NA} on the proper side of the ball with radius {R = 2*r_abs/S}: */
    if (cl == 1)
      { /* Pick a point in {U^NA} outside the ball: */
        /* This should be fast because {B <= H/2}. */
        do 
          { rn_throw_cube (NA, p); }
        while (rn_norm_sqr(NA, p) < r_rel*r_rel);
      }
    else
      { /* Pick a point inside the ball: */
        /* Since we may have {B << H} we can't do trial ane error. */
        rn_throw_ball (NA, p);
        rn_scale(NA, r_rel, p, p);
      }

    /* Now scale and shift to the proper cubelet: */
    for (k = 0; k < NA; k++) { p[k] = S*(ix[k] + (p[k] + 1)/2) - 1; } 
    (*class) = cl;
  }

void rn_classif_test_check_shells(int *NA, int *NC)
  { if ((*NA) == 0) { (*NA) = 2; }
    if ((*NC) == 0) { (*NC) = 2; }
    demand((*NA) == 2, "invalid number of attributes");
    demand((*NC) >= 2, "invalid number of classes");
  }

int rn_classif_test_label_shells(int NA, int NC, double p[])
  {
    assert(NA == 2);
    assert(NC >= 2);
    int cl;
    double rp = rn_norm(NA, p);
    if (rp >= 1)
      { /* Outside the unit ball: */ 
        cl = 0;
      }
    else
      { /* Find the shell containing {p}: */ 
        double vp = pow(rp, NA); /* Measure of ball with radius {rp}, in {[0_NC]} */
        cl = (int)ceil(NC*vp + 0.5);
        assert(cl >= 1);
        if (cl > NC) { cl = 1; }
      }
    return cl;
  }

void rn_classif_test_throw_shells(int i, int NA, int NC, double p[], int *class)
  { assert(NA == 2);
    assert(NC >= 2);
    assert(i >= 0);
    /* Cycle among classes: */
    int cl = (i % NC) + 1;
    /* Now throw a point until it falls in the right class: */
    /* !!! There is a better way but I am too lazy to work out the math... !! */
    int clp;
    do
      { /* Throw point {p} in the unit ball: */ 
        rn_throw_ball(NA, p);
         /* Make sure {p} is inside {U^NA}: */
        rn_scale(NA, 0.999, p, p);
        /* Classify that point: */ 
        clp = rn_classif_test_label_shells(NA, NC, p);
      }
    while (clp != cl);
    (*class) = cl;
  }

void rn_classif_test_check_staoli(int *NA, int *NC)
  { if ((*NA) == 0) { (*NA) = 2; }
    if ((*NC) == 0) { (*NC) = 2; }
    demand((*NA) >= 1, "invalid number of attributes");
    demand((*NC) == 3, "invalid number of classes");
  }

void rn_classif_test_box_widths_staoli(int NA, double *wd2, double *wd3);
  /* Computes the side {*wd2} of the fat box and the smaller size {*wd3}
    of the thin box in the {staoli} problem. */

void rn_classif_test_get_box_staoli(int NA, int cl, double *xctr, double *xrad, double *yrad);
  /* Computes the dimension and position of the box with class {cl}
    (either 2 or 3) in the {staoli} problem with dimension {NA}. The
    box center is at coordinate {*xctr} on axis zero; the box extends
    plus or minus {*xrad} along that axis and plus or minus {*yrad}
    along the other axes. */

void rn_classif_test_box_widths_staoli(int NA, double *wd2, double *wd3)
  { 
    assert(NA >= 1);
    double g = pow(4, ((double)NA-1)/((double)NA));
    /* Choose the sizes so that the boxes are {2*wd3} apart */
    /* and at least {wd3} away from the edges of {U^NA}. */
    (*wd3) = 2/(g + 5);
    (*wd2) = g*(*wd3);
  }

void rn_classif_test_get_box_staoli(int NA, int cl, double *xctr, double *xrad, double *yrad)
  { 
    assert(NA >= 1);
    assert((cl == 2) || (cl == 3));
    double wd2, wd3;
    rn_classif_test_box_widths_staoli(NA, &wd2, &wd3);
    double b = 1 - (wd2 + 3*wd3)/2; /* Margin around boxes. */
    if (cl == 2)
      { /* Uniform box: */
        (*xctr) = - 1.0 + b + wd2/2;
        (*xrad) = wd2/2;
        (*yrad) = wd2/2;
      }
    else
      { /* Thin box: */
        (*xctr) = + 1.0 - b - wd3/2;
        (*xrad) = wd3/2;
        (*yrad) = 2*wd3;
      }
  }

int rn_classif_test_label_staoli(int NA, int NC, double p[])
  {
    assert(NA >= 1);
    assert(NC == 3);
    /* Check against the domains {2..NC}: */
    int cl;
    for (cl = 2; cl <= NC; cl++) 
      { /* Get the geometry of the {cl} domain: */
        double xctr, xrad, yrad;
        rn_classif_test_get_box_staoli(NA, cl, &xctr, &xrad, &yrad);
        /* Check {p} against that box: */
        bool_t inside = TRUE;
        int t;
        for (t = 0; t < NA; t++)
          { double zctr = (t == 0 ? xctr : 0);
            double zrad = (t == 0 ? xrad : yrad);
            if (fabs(p[t] - zctr) > zrad) { inside = FALSE; }
          }
        if (inside) { return cl; }
      }
    /* Not in any domain {2..NC}: */
    return 1;
  }

void rn_classif_test_throw_staoli(int i, int NA, int NC, double p[], int *class)
  { assert(NA >= 1);
    assert(NC == 3);
    assert(i >= 0);
    /* Compute the measures of the three domains: */
    double wd2, wd3;
    rn_classif_test_box_widths_staoli(NA, &wd2, &wd3);
    double v2 = pow(wd2,NA);
    double v3 = wd3*pow(4*wd3,NA-1);
    assert(fabs(v2 - v3) < 0.001*(v2 + v3)); 
    double vu = pow(2, NA); /* Measure of {U^NA}. */
    assert(v2 + v3 < vu);
    /* Choose the period {NP} so that the densities are about the same: */
    int NP = (int)(2*vu/(v2 + v3) + 0.5);
    if (NP < 3) { NP = 3; }
    /* Cycle among domains: */
    int phase = i % NP;
    int cl;
    if (phase == 0)
      { cl = 2; }
    else if (phase == 1)
      { cl = 3; }
    else
      { cl = 1; }
    /* Now throw a point: */
    if (cl == 1)
      { /* Throw a point in {U^NA} until it is in class 1: */
        int clp;
        do
          { /* Throw point {p} in the unit ball: */ 
            rn_throw_cube(NA, p);
            /* Classify that point: */ 
            clp = rn_classif_test_label_staoli(NA, NC, p);
          }
        while (clp != cl);
      }
    else
      { /* Throw a point in {U^NA}: */
        rn_throw_cube(NA, p);
        /* Scale and shift to the domain {cl}: */
        double xctr, xrad, yrad;
        rn_classif_test_get_box_staoli(NA, cl, &xctr, &xrad, &yrad);
        p[0] = p[0]*xrad + xctr;
        int t;
        for (t = 1; t < NA; t++) { p[t] = p[t]*yrad; }
      }
    (*class) = cl;
  }
