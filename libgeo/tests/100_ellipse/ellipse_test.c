/* ellipse_test --- test program for ellipse_crs.h, ellipse_ouv.h  */
/* Last edited on 2013-05-24 20:28:54 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <affirm.h>
#include <bool.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <r2.h>

#include <ellipse_crs.h>
#include <ellipse_ouv.h>

/* GENERAL PARAMETERS */

#define N_RUNS 1000
  /* Number of test runs. */

#define Pr fprintf
#define Er stderr

/* INTERNAL PROTOTYPES */

int main (int argc, char **argv);

void test_ellipse_tools(int trial, bool_t verbose);
  /* Tests the function {ellipse_ouv_nearest_point}. */

ellipse_crs_t throw_ellipse_crs(double R);
  /* Generates a random ellipse in the square {[-R_+R]×[-R_+R]}. */

r2_t throw_point(ellipse_crs_t *E, double R);
  /* Generates a random point in the square {[-R_+R]×[-R_+R]}
    that is an interesting test case for the ellipse {E}. */

/* CHECKING LSQ FITTING */

void check_ellipse_nearest_point(ellipse_crs_t *E, r2_t *p, r2_t *q, double dpq);
  /* Checks whether {q,dpq} are the point of {E} nearest to {p} and the signed
    distance from {p} to [E}, respectively. */

void check_ellipse_border_position(ellipse_crs_t *E, r2_t *p, double hwd, double fp);
  /* Checks whether {fp} is the position of {p} relative to the boundary
    of {E} fattened by {hwd}.  Assumes that {ellipse_crs_nearest_point}
    has been checked. */

/* IMPLEMENTATIONS */

int main (int argc, char **argv)
  { int i;
    for (i = 0; i < N_RUNS; i++) 
      { test_ellipse_tools(i, i < 5); }
    fclose(Er);
    fclose(stdout);
    return (0);
  }

void test_ellipse_tools(int trial, bool_t verbose)
  { 
    srand(1665 + 2*trial);
    srandom(1665 + 2*trial);

    Pr(Er, "\n");
    Pr(Er, "======================================================================\n");
    Pr(Er, "%s (%d)\n", __FUNCTION__, trial);
    
    /* Generate the test data. */
    double R = 500;
    ellipse_crs_t E;    /* Test ellipse. */
    r2_t p;              /* Test point. */
    if (trial == 0)
      { /* Simple case --- a centered circle: */
        E.ctr = (r2_t) {{ 0,0 }};
        E.rad = 100;
        E.str = (r2_t) {{ 0,0 }};
        p = (r2_t) {{ 100, 100 }};
      }
    else if (trial == 1)
      { E.ctr = (r2_t) {{ 0,0 }};
        E.rad = 100;
        E.str = (r2_t) {{ 50,0 }};
        p = (r2_t) {{ 100, 100 }};
      }
    else if (trial == 2)
      { E.ctr = (r2_t) {{ 8.508132, -244.675234 }};
        E.rad = 7.793289;
        E.str = (r2_t) {{ 251.915949, -482.314121 }};
        p = (r2_t) {{ -86.2103454, -143.890590 }};
      }
    else if (trial == 3)
      { E.ctr = (r2_t) {{ 224.00000, 413.00000 }};
        E.rad = 29.97611;
        E.str = (r2_t) {{ -2.27656, 3.91519 }};
        p = (r2_t) {{ -0.72380, -2.23947 }};
      }
    else
      { E = throw_ellipse_crs(R);
        p = throw_point(&E, R);
      }
      
    if (verbose) 
      { Pr(Er, "test case:\n");
        Pr(Er, "  ellipse (crs) = "); ellipse_crs_print(Er, &E, "%24.16e"); Pr(Er, "\n");
      }
    
    ellipse_ouv_t F;
    ellipse_crs_to_ouv(&E, &F);
    
    if (verbose) 
      { Pr(Er, "test case:\n");
        Pr(Er, "  ellipse (ouv) = "); ellipse_ouv_print(Er, &F, "%24.16e"); Pr(Er, "\n");
      }
   
    if (verbose) { Pr(Er, "  point = "); r2_print(Er, &p); Pr(Er, "\n"); }
    
    /* Test {ellipse_crs_nearest_point}: */ 
    r2_t q;
    double dpq = ellipse_crs_nearest_point(&E, &p, &q);
    if (verbose) 
      { Pr(Er, "  found = "); r2_print(Er, &q);  
        Pr(Er, "  dist = %+12.6f", dpq);
        Pr(Er, "\n");
      }
    check_ellipse_nearest_point(&E, &p, &q, dpq);
    
    /* Test {ellipse_crs_border_position}: */ 
    double hwd = ( drandom() < 0.250 ? 1.5+drandom() : 0.20*drandom()) * E.rad;
    double fp = ellipse_crs_border_position(&E, hwd, &p);
    if (verbose) 
      { Pr(Er, "  bpos = %+10.6f", fp);
        Pr(Er, "\n");
      }
    check_ellipse_border_position(&E, &p, hwd, fp);
    
    Pr(Er, "done.\n");
    Pr(Er, "======================================================================\n");
  }

void check_ellipse_nearest_point(ellipse_crs_t *E, r2_t *p, r2_t *q, double dpq) 
  { 
    /* Grab vectrs ans semidiameters: */
    ellipse_ouv_t F;
    ellipse_crs_to_ouv(E, &F);
    double a = F.a;
    double b = F.b;
    
    /* General magnitude of dimensions: */
    double mag = r2_dist(p, &(E->ctr)) + a;
    
    /* Compute the UV coordinates of {p}: */
    r2_t hp; r2_sub(p, &(E->ctr), &hp);
    double up = r2_dot(&hp, &(F.u));
    double vp = r2_dot(&hp, &(F.v));

    /* Compute the UV coordinates of {q}: */
    r2_t hq; r2_sub(q, &(E->ctr), &hq);
    double uq = r2_dot(&hq, &(F.u));
    double vq = r2_dot(&hq, &(F.v));

    /* Compute rel radial positions of {p,q}: */
    double fp, fq;
    if (a <= 0)
      { /* Ellipse is a point: */
        fp = ((up == 0) && (vp == 0) ? 00 : +1);
        fq = ((uq == 0) && (vq == 0) ? 00 : +1);
      }
    else if (b <= 0)
      { /* Ellipse is a segment: */
        fp = ((vp == 0) && (fabs(up) <= b) ? 00 : +1);
        fq = ((vq == 0) && (fabs(uq) <= b) ? 00 : +1);
      }
    else
      { /* Ellipse has non-empty interior: */
        fp = hypot(up/a, vp/b);
        fq = hypot(uq/a, vq/b);
      }
      
    /* Compute the UV gradient {ug,vg} at {q}: */
    double ug = (a == 0 ? 0.0 : 2*uq/(a*a));
    double vg = (b == 0 ? 0.0 : 2*vq/(b*b));
    
    /* Compute the angle between {q-p} and the gradient: */
    double G = ug*(vq - vp) - vg*(uq - up);
    double mg = hypot(ug,vg);
    double mpq = hypot(uq - up, vq - vp);
    double sinG = (mg*mpq == 0 ? 0.0 : G/mg/mpq);
      
    auto void dump_all(void);
    
    void dump_all(void)
      {
        Pr(Er, "----------------------------------------\n");
        Pr(Er, "  ellipse ="); ellipse_crs_print(Er, E, "%24.16e"); Pr(Er, "\n");
        Pr(Er, "  p ori ="); r2_print(Er, p); Pr(Er, "\n");
        Pr(Er, "  q fin ="); r2_print(Er, q); Pr(Er, "\n");
        Pr(Er, "\n");
        Pr(Er, "  p rel = [ %24.16e %24.16e ]\n", up, vp);
        Pr(Er, "  q rel = [ %24.16e %24.16e ]\n", uq, vq);
        Pr(Er, "\n");
        Pr(Er, "  fp = %24.16e\n", fp);
        Pr(Er, "  fq = %24.16e\n", fq);
        Pr(Er, "\n");
        Pr(Er, "  grad  = ( %24.16e %24.16e )\n", ug, vg);
        Pr(Er, "  disp  = ( %24.16e %24.16e )\n", uq - up, vq - vp); 
        Pr(Er, "  cros  = %24.16e\n", G);  
        Pr(Er, "  sinG  = %24.16e\n", sinG); 
        Pr(Er, "----------------------------------------\n");
      }
    
    if (a <= 0)
      { /* Ellipse is a point: */
        double tol_q = 1.0e-7*mag;
        if (r2_dist(&(E->ctr), q) > tol_q)
          { dump_all(); demand(FALSE, "{q} should be the center"); }
      }
    else if (b <= 0)
      { /* Ellipse is a segment: */
        double tol_q = 1.0e-7*mag;
        if (fabs(vq) > tol_q)
          { dump_all(); demand(FALSE, "{q} should be the on {ctr--u} line"); }
        if (fabs(uq) - a > tol_q)
          { dump_all(); demand(FALSE, "{q} should be on major diameter"); }
      }
    else
      {
        /* Check whether {q} is on the ellipse: */
        double tol_fq = 1.0e-7;
        if (fabs(fq - 1.0) > tol_fq)
          { dump_all(); demand(FALSE, "point {q} not on ellipse"); }

        /* Check whether {p-q} is collinear with {gradF(q)}: */
        double tol_G = 1.0e-6*mag/b;
        if (fabs(G) > tol_G)
          { dump_all(); demand(FALSE, "vector {q-p} not normal to {E}"); }
      }
    
    /* Check the distance and sign: */
    double xpq = r2_dist(p, q);
    if (fp < 1.0) { xpq = -xpq; }
    double tol_dist = 1.0e-7*mag;
    if (fabs(xpq - dpq) > tol_dist)
      { dump_all(); demand(FALSE, "distance is incorrect"); }
  }

void check_ellipse_border_position(ellipse_crs_t *E, r2_t *p, double hwd, double fp) 
  { 
    /* Grab vectrs ans semidiameters: */
    ellipse_ouv_t F;
    ellipse_crs_to_ouv(E, &F);
    /* double a = F.a; */
    /* double b = F.b; */
    
    /* General magnitude of dimensions: */
    /* double mag = r2_dist(p, &(E->ctr)) + a; */
    
    /* Compute the UV coordinates of {p}: */
    r2_t hp; r2_sub(p, &(E->ctr), &hp);
    double up = r2_dot(&hp, &(F.u));
    double vp = r2_dot(&hp, &(F.v));

    /* Compute the signed distance from the border: */
    double dpE = ellipse_crs_nearest_point(E, p, NULL);
    
    /* Compute the expected border-relative position {fp_exp}: */
    double fp_exp = fmax(-hwd, fmin(+hwd, dpE))/hwd;
    
    auto void dump_all(void);
    
    void dump_all(void)
      {
        Pr(Er, "ellipse ="); ellipse_crs_print(Er, E, "%24.16e"); Pr(Er, "\n");
        Pr(Er, "p ori ="); r2_print(Er, p); Pr(Er, "\n");
        Pr(Er, "\n");
        Pr(Er, "p rel = [ %24.16e %24.16e ]\n", up, vp);
        Pr(Er, "\n");
        Pr(Er, "f cmp = %24.16e\n", fp);
        Pr(Er, "f exp = %24.16e\n", fp_exp);
        Pr(Er, "\n");
      }
    
    /* Check border relative position: */
    double tol_fp = 1.0e-7;
    if (fabs(fp - fp_exp) > tol_fp)
      { dump_all(); demand(FALSE, "border position is incorrect"); }
  }

ellipse_crs_t throw_ellipse_crs(double R)
  {
    ellipse_crs_t E;
    
    /* Half small, half big: */
    E.rad = (drandom() < 0.500 ? 8 : R) * drandom();
    
    /* Half roundish, half not: */
    double K = (drandom() < 0.500 ? 24 : R);
    E.str.c[0] = K*(2*drandom() - 1);
    E.str.c[1] = K*(2*drandom() - 1);
    
    double xmrg = fabs(E.str.c[0]) + E.rad;
    double ymrg = fabs(E.str.c[1]) + E.rad;
    E.ctr.c[0] = (R - 2*xmrg)*(2*drandom() - 1);
    E.ctr.c[1] = (R - 2*ymrg)*(2*drandom() - 1);
    
    return E;
  }

r2_t throw_point(ellipse_crs_t *E, double R)
  {
    /* Get the main vectors and semi-diameters: */
    ellipse_ouv_t F;
    ellipse_crs_to_ouv(E, &F);
    double a = F.a;
    double b = F.b;

    r2_t p; /* The chosen point. */
    
    double coin = drandom();
    if ((coin -= 0.100) < 0)
      { p = E->ctr; }
    else if ((coin -= 0.100) < 0)
      { r2_mix(1.0, &(E->ctr), 2*a*drandom(), &(F.u), &p); }
    else if ((coin -= 0.100) < 0)
      { r2_mix(1.0, &(E->ctr), 2*b*drandom(), &(F.v), &p); }
    else
      { p = (r2_t){{ R*(2*drandom() - 1), R*(2*drandom() - 1) }}; }
    return p;
  }

