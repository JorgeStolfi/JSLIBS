/* See pst_geom.h */
/* Last edited on 2009-03-03 13:51:22 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <r3.h> 
#include <hr3.h> 
#include <r3x3.h> 
#include <jsrandom.h>
#include <affirm.h>

#include <pst_geom.h>
#include <pst_sphere_basic.h>

#define Pr fprintf
#define Er stderr

#define X c[0]
#define Y c[1]
#define Z c[2]
  /* Cartesian coordinates of {r2_t} or {r3_t}. */

#define hm c.c[0]
#define hx c.c[1]
#define hy c.c[2]
#define hz c.c[3]
  /* Homogeneous coordinates of an {hr2_point_t} or {hr3_point_t}. */

/* IMPLEMENTATIONS */

void pst_sphere_G_D_R_to_dst_rad_maj
  ( double G,      /* Angular spread {1/F} of camera (possibly zero). */
    double D,      /* Dist from optical center to sphere center (finite). */
    double R,      /* The radius of the sphere. */
    double *dst,   /* (OUT) Dist from {Q} to center {ctr} of ellipse, or {+INF}. */
    double *rad,   /* (OUT) Length of minor semidiameter of ellipse. */
    double *maj    /* (OUT) Length of major semidiameter of ellipse. */
  )
  {
    bool_t debug = TRUE;
    demand(isfinite(G), "invalid {G}");
    demand(isfinite(R) && (R >= 0), "invalid {R}");
    demand(isfinite(D) && (D >= 0), "invalid {D}");
    /* Compute {beta2}, the minor radius over the radius of the sphere, squared: */
    double t = R*G, t2 = t*t;
    demand(t < 1, "camera is lower than top of sphere");
    double beta2 = 1/(1 - t2); 
    /* assert(beta2 >= 1); */
    /* Compute the apparent transverse radius (minor semidiameter): */
    double radT = R*sqrt(beta2); /* Worth computing, for {*rad} and {*str} */
    if (dst != NULL) { (*dst) = beta2*D; }
    if (rad != NULL) { (*rad) = radT; }
    if (maj != NULL)
      { /* Compute {alfa2}, the major radius over the minor radius, squared: */
        double q = D*G, q2 = q*q;
        double alfa2 = 1 + beta2*q2;
        /* Not worth testing for {alfa2 == 1}: */
        (*maj) = sqrt(alfa2)*radT;
      }
    if (debug) 
      { if (dst != NULL) { Pr(Er, "  dst = %18.12f\n", (*dst)); }
        if (rad != NULL) { Pr(Er, "  rad = %18.12f\n", (*rad)); }
        if (maj != NULL) { Pr(Er, "  maj = %18.12f\n", (*maj)); }
      }
  }

void pst_sphere_G_dst_rad_to_D_R_maj
  ( double G,      /* Angular spread {1/F} of camera. */
    double dst,    /* Dist from {Q} to center {ctr} of ellipse (finite). */
    double rad,    /* Length of minor semidiameter of ellipse. */
    double *D,     /* (OUT) Dist from optical center to sphere center (finite). */
    double *R,     /* (OUT) The radius of the sphere. */
    double *maj    /* (OUT) Length of major semidiameter. */
  )
  {
    bool_t debug = TRUE;
    demand(isfinite(G), "invalid {G}");
    demand(isfinite(dst) && (dst >= 0), "invalid {dst}");
    demand(isfinite(rad) && (rad >= 0), "invalid {rad}");
    /* Compute the the ratio {u=rad/F} where {F} is the focal length: */
    double u = rad*G, u2 = u*u;
    /* Compute {beta2}, the minor radius over the radius of the sphere, squared: */
    double beta2 = 1 + u2;
    if (debug) { Pr(Er, "  beta2 = %16.12f beta = %16.12f\n", beta2, sqrt(beta2)); }
    /* assert(beta2 >= 1); */
    /* Compute the radius of the sphere {*R}: */
    if (D != NULL) { (*D) = dst/beta2; }
    if (R != NULL) { (*R) = rad/sqrt(beta2); }
    if (maj != NULL)
      { /* Compute {alfa2}, the major radius over the minor radius, squared: */
        double dG = dst*G, d2G2 = dG*dG;
        double alfa2 = 1 + beta2*d2G2;
        if (debug) { Pr(Er, "  alfa2 = %16.12f alfa = %16.12f\n", alfa2, sqrt(alfa2)); }
        (*maj) = sqrt(alfa2)*rad;
      }
    if (debug) 
      { if (D != NULL)   { Pr(Er, "  D =   %18.12f\n", (*D)); }
        if (R != NULL)   { Pr(Er, "  R =   %18.12f\n", (*R)); }
        if (maj != NULL) { Pr(Er, "  maj = %18.12f\n", (*maj)); }
      }
  }

void pst_sphere_G_rad_maj_to_R_D_dst
  ( double G,      /* Angular spread {1/F} of camera. */
    double rad,    /* Length of minor semidiameter of ellipse. */
    double maj,    /* Length of major semidiameter of ellipse. */
    double *R,     /* (OUT) The radius of the sphere. */
    double *D,     /* (OUT) dist from {Q} to center {K} of sphere, or {+INF}. */
    double *dst    /* (OUT) Dist from {Q} to center {ctr} of ellipse, or {+INF}. */
  )
  {
    demand(FALSE, "!!! not implemented yet !!!");
  }

void pst_sphere_dst_rad_maj_to_G_D_R
  ( double dst,    /* Dist from {Q} to center {ctr} of ellipse (finite). */
    double rad,    /* Length of minor semidiameter of ellipse. */
    double maj,    /* Length of major semidiameter of ellipse. */
    double *G,     /* (OUT) Angular spread {1/F} of camera. */
    double *D,     /* (OUT) dist from {Q} to center {K} of sphere. */
    double *R      /* (OUT) The radius of the sphere. */
  )
  {
    bool_t debug = TRUE;
    demand(isfinite(dst) && (dst >= 0),   "invalid {DST}");
    demand(isfinite(rad) && (rad > 0),    "invalid {rad}");
    demand(isfinite(maj) && (maj >= rad), "invalid {maj}");
    if (rad == maj)
      { /* Top view, assume view from infinity: */
        if (G != NULL) { (*G) = 0; /* Arbitrary choice. */ }
        if (D != NULL) { (*D) = dst; /* Actually meaningless, but... */ }
        if (R != NULL) { (*R) = rad; }
      }
    else
      { /* Oblique view from finite distance: */
        demand(FALSE, "!!! not implemented yet !!!");
      }
    if (debug) 
      { if (G != NULL) { Pr(Er, "  G =   %18.12f\n", (*G)); }
        if (D != NULL) { Pr(Er, "  D =   %18.12f\n", (*D)); }
        if (R != NULL) { Pr(Er, "  R =   %18.12f\n", (*R)); }
      }
  }
  
hr3_pmap_t pst_sphere_perspective_view_map(hr3_point_t *O, r2_t *K, double R)
  { 
    /* Compute the sphere center to camera direction {wdr} and distance {D}: */
    r3_t wdr; double D;
    double Om = O->hm;
    if (Om == 0)
      { double x = O->hx;
        double y = O->hy;
        double z = O->hz;
        wdr = (r3_t) {{ x, y, z }};
        double dd = r3_dir(&wdr, &wdr);
        demand(dd > 0, "invalid camera viewpoint");
        D = +INF;
      }
    else
      { /* Get the Cartesian NCS coordinates {O} of the viewpoint: */
        r3_t OC = r3_from_hr3(O);
        wdr = (r3_t) {{ K->X - OC.X, K->Y - OC.Y, OC.Z }};
        D = r3_dir(&wdr, &wdr);
      }
    
    /* Compute the other two cardinal directions {udr,vdr} of the VCS: */
    r3_t udr = (r3_t) {{ -wdr.Y, +wdr.X, 0.0 }};
    double du = r3_dir(&udr, &udr);
    if (du == 0) { udr = (r3_t) {{ 1.0, 0.0, 0.0 }}; }
    r3_t vdr; r3_cross(&wdr, &udr, &vdr);
    (void)r3_dir(&vdr, &vdr); /* Just to be sure... */
    
    /* Assemble the rotation part of the VCS-to-NCS matrix {M.inv}: */
    hr3_pmap_t M;
    int j;
    for (j = 0; j < 4; j++) 
      { M.inv.c[0][j] = (j == 0 ? 1 : 0); /* For now. */
        M.inv.c[1][j] = (j == 0 ? 0 : R*udr.c[j-1]);
        M.inv.c[2][j] = (j == 0 ? 0 : R*vdr.c[j-1]);
        M.inv.c[3][j] = (j == 0 ? 0 : R*wdr.c[j-1]);
      }

    if (D != +INF)
      { /* We must pre-multiply {M.inv} by a cylindrical-to-conical matrix {A}, */
        /* and post-multiply by a translation from {K'} to {K}. */
        
        /* Compute the distance {L} from {K} to {K'}: */
        double R = R;                     /* radius of the sphere of sphere: */
        double T = sqrt((D - R)*(D + R));   /* Dist from viewpoint to horizon. */
        double S = R*(T/D);               /* Radius of horizon. */
        double L = R*R/T;               /* Dist from sphere ctr to horizon ctr. */

        /* Append to {M.inv} the translation from {K'} to {K}: */
        for (j = 1; j < 4; j++) { M.inv.c[0][j] = -L*wdr.c[j-1]; }

        /* The matrix {A} must take the camera from {Z==+INF} to {Z==D-L}. */
        /* The unit sphere remains a unit sphere, but squeezed down along the */
        /* {W} axis, so that the equator becomes the horizon of the new camera. */
        r4x4_t A; r4x4_ident(&A);
        double g = (R-L)/(D-R);
        A.c[0][0] = R;
        A.c[1][1] = S;
        A.c[2][2] = S;
        A.c[3][0] = g*R;
        A.c[3][3] = g*(D-L);

        /* Debugging the code: */
        bool_t debug_A = TRUE;
        if (debug_A)
          { int i, j, k;
            for (i = 0; i < 4; i++)
              for (k = +1; k>= -1; k -= 2)
                { r4_t p = (r4_t) {{ 0, 0, 0, 0 }};
                  p.c[i] = k; 
                  for (j = 0; j < 2; j++)
                    { r4_t q;
                      r4x4_map_row(&p, &A, &q);
                      double qt = q.c[0];
                      double qi = q.c[i];
                      if (qt != 0)
                        { r4_scale(1/qt, &q, &q); }
                      else if (qi != 0)
                        { r4_scale(1/qi, &q, &q); }
                      r4_print(Er, &p); 
                      Pr(Er, " --> ");
                      r4_print(Er, &q); 
                      Pr(Er, "\n");
                      p.c[0] = 1;
                    }
                }
          }
        /* Combine {A} with {M.inv}: */
        r4x4_mul(&A, &(M.inv), &(M.inv));
      }
    
    /* Compute the NCS-to-VCS matrix: */
    r4x4_inv(&(M.inv), &(M.dir));
    
    return M;
  }
    

