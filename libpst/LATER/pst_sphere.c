/* See pst_geom.h */
/* Last edited on 2009-03-03 14:06:25 by stolfi */

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
#include <pst_argparser.h>

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

void pst_sphere_compute_rad_maj_dst_from_G_R_D
  ( pst_sphere_t *S, /* The sphere. */
    hr3_point_t *O,  /* Camera's viepoint. */
    r2_t *ctr,       /* (OUT) The center of the sphere's projection on the image. */
    double *rad,     /* (OUT) The minor semidiameter of that projection. */
    r2_t *str        /* (OUT) The stretch vector of the projection. */
  )
  {
    demand(S->R > 0, "bad sphere radius");
    double Om = O->hm; demand((! isnan(Om)) && (Om >= 0), "bad {O.m}");
    if (Om == 0)
      { /* Parallel projection (a rather common case): */
        if (rad != NULL) { (*rad) = S->R; }
        if (ctr != NULL) { (*ctr) = S->K; }
        if (str != NULL) { (*str) = (r2_t) {{ 0, 0 }}; }
        return;
      }
    double Oz = O->hz; demand(! isnan(Oz), "{O.z} is NAN");
    demand(Oz > 0, "bad {O.z}");
    /* Get the center and spread of the camera: */
    double G = Om/Oz;
    /* Compute {beta2}, the minor radius over the radius of the sphere, squared: */
    double t = S->R*G, t2 = t*t;
    demand(t < 1, "camera is lower than top of sphere");
    double beta2 = 1/(1 - t2); 
    /* assert(beta2 >= 1); */
    /* Compute the apparent transverse radius (minor semidiameter): */
    double radT = S->R*sqrt(beta2); /* Worth computing, for {*rad} and {*str} */
    if (radT != NULL) { (*rad) = radT; }
    if ((ctr != NULL) || (str != NULL))
      { /* Get the optical center {Q}: */
        double Ox = O->hx; demand(! isnan(Ox), "{O.x} is NAN");
        double Oy = O->hy; demand(! isnan(Oy), "{O.y} is NAN");
        r2_t Q = (r2_t) {{ Ox/Om, Oy/Om }};
        /* Get the displacement {Dsp} from {Q} to {K}: */
        r2_t Dsp; r2_sub(S->K, &Q, &Dsp);
        if (ctr != NULL) 
          { /* Compute the apparent center {*ctr}: */
            /* Not worth testing for {beta2 == 1}: */
            r2_mix(1.0, &Q, beta2, &Dsp, ctr);
          }
        if (str != NULL)
          { /* Compute {D}, the distance from {Q} to {K}: */
            double D = r2_norm(&Dsp);
            /* Compute {alfa2}, the major radius over the minor radius, squared: */
            double q = D*G, q2 = q*q;
            double alfa2 = 1 + beta2*q2;
            /* Compute the stretch vector {*str}: */
            /* Not worth testing for {alfa2 == 1}: */
            r2_scale((sqrt(alfa2)- 1)*radT/D, &Dsp, str);
          }
      }
  }

void pst_sphere_compute_ctr_rad_str_from_S_O
  ( pst_sphere_t *S, /* The sphere. */
    hr3_point_t *O,  /* Camera's viepoint. */
    r2_t *ctr,       /* (OUT) The center of the sphere's projection on the image. */
    double *rad,     /* (OUT) The minor semidiameter of that projection. */
    r2_t *str        /* (OUT) The stretch vector of the projection. */
  )
  {
    demand(S->R > 0, "bad sphere radius");
    double Om = O->hm; demand((! isnan(Om)) && (Om >= 0), "bad {O.m}");
    if (Om == 0)
      { /* Parallel projection (a rather common case): */
        if (rad != NULL) { (*rad) = S->R; }
        if (ctr != NULL) { (*ctr) = S->K; }
        if (str != NULL) { (*str) = (r2_t) {{ 0, 0 }}; }
        return;
      }
    double Oz = O->hz; demand(! isnan(Oz), "{O.z} is NAN");
    demand(Oz > 0, "bad {O.z}");
    /* Get the center and spread of the camera: */
    double G = Om/Oz;
    /* Compute {beta2}, the minor radius over the radius of the sphere, squared: */
    double t = S->R*G, t2 = t*t;
    demand(t < 1, "camera is lower than top of sphere");
    double beta2 = 1/(1 - t2); 
    /* assert(beta2 >= 1); */
    /* Compute the apparent transverse radius (minor semidiameter): */
    double radT = S->R*sqrt(beta2); /* Worth computing, for {*rad} and {*str} */
    if (radT != NULL) { (*rad) = radT; }
    if ((ctr != NULL) || (str != NULL))
      { /* Get the optical center {Q}: */
        double Ox = O->hx; demand(! isnan(Ox), "{O.x} is NAN");
        double Oy = O->hy; demand(! isnan(Oy), "{O.y} is NAN");
        r2_t Q = (r2_t) {{ Ox/Om, Oy/Om }};
        /* Get the displacement {Dsp} from {Q} to {K}: */
        r2_t Dsp; r2_sub(S->K, &Q, &Dsp);
        if (ctr != NULL) 
          { /* Compute the apparent center {*ctr}: */
            /* Not worth testing for {beta2 == 1}: */
            r2_mix(1.0, &Q, beta2, &Dsp, ctr);
          }
        if (str != NULL)
          { /* Compute {D}, the distance from {Q} to {K}: */
            double D = r2_norm(&Dsp);
            /* Compute {alfa2}, the major radius over the minor radius, squared: */
            double q = D*G, q2 = q*q;
            double alfa2 = 1 + beta2*q2;
            /* Compute the stretch vector {*str}: */
            /* Not worth testing for {alfa2 == 1}: */
            r2_scale((sqrt(alfa2)- 1)*radT/D, &Dsp, str);
          }
      }
  }

void pst_sphere_compute_R_D_maj_from_G_dst_rad
  ( double G,          /* Camera resolution {1/F}. */
    double rad,        /* Length of minor semidiameter of ellipse. */
    double dst,        /* Dist from {Q} to center {ctr} of ellipse, or {+INF}. */
    double *R,         /* (OUT) The radius of the sphere. */
    double *D,         /* (OUT) dist from {Q} to center {K} of sphere, or {+INF}. */
    double *maj        /* (OUT) Length of major semidiameter. */
  )
  {
    bool_t debug = TRUE;
    demand(isfinite(G), "invalid {G}");
    demand(! isnan(dst) && (dst >= 0), "invalid {dst}");
    demand(! isnan(rad) && rad > 0, "invalid {rad}");
    /* Compute the the ratio {u=rad/F} where {F} is the focal length: */
    double u = rad*G, u2 = u*u;
    /* Compute {beta2}, the minor radius over the radius of the sphere, squared: */
    double beta2 = 1 + u2;
    if (debug) { Pr(Er, "  beta2 = %16.12f beta = %16.12f\n", beta2, sqrt(beta2)); }
    /* assert(beta2 >= 1); */
    /* Compute the radius of the sphere {*R}: */
    if (R != NULL) { (*R) = rad/sqrt(beta2); }
    if (D != NULL) { (*D) = dst/beta2; }
    if (maj != NULL)
      { /* Compute {alfa2}, the major radius over the minor radius, squared: */
        double dG = dst*G, d2G2 = dG*dG;
        double alfa2 = 1 + beta2*d2G2;
        if (debug) { Pr(Er, "  alfa2 = %16.12f alfa = %16.12f\n", alfa2, sqrt(alfa2)); }
        (*maj) = sqrt(alfa2)*rad;
      }
    if (debug) 
      { if (R != NULL)   { Pr(Er, "  R =   %18.12f\n", (*R)); }
        if (D != NULL)   { Pr(Er, "  D =   %18.12f\n", (*D)); }
        if (maj != NULL) { Pr(Er, "  maj = %18.12f\n", (*maj)); }
      }
  }

void pst_sphere_compute_R_D_G_c_s_from_rad_maj_dst
  ( double rad,        /* Length of minor semidiameter of ellipse. */
    double maj,        /* Length of major semidiameter of ellipse. */
    double dst,        /* Dist from {Q} to center {ctr} of ellipse, or {+INF}. */
    double *R,         /* (OUT) The radius of the sphere. */
    double *D,         /* (OUT) dist from {Q} to center {K} of sphere, or {+INF}. */
    double *G,         /* (OUT) Camera resolution inferred from {E,Q}. */
    double *c,         /* (OUT) Co-sine of elevation of {O} seen from {K}. */
    double *s          /* (OUT) Sine of elevation of {O} seen from {K}. */
  )
  {
    bool_t debug = TRUE;
    demand(! isnan(dst) && (dst >= 0), "invalid {DST}");
    demand(! isnan(rad) && (rad > 0), "invalid {rad}");
    demand(! isnan(maj) && (maj >= rad), "invalid {maj}");
    if (rad == maj)
      { /* Top view, assume view from infinity: */
        if (R != NULL) { (*R) = rad; }
        if (D != NULL) { (*D) = 0; /* Actually meaningless, but... */ }
        if (G != NULL) { (*G) = 0; /* Arbitrary choice. */ }
      }
    else if (dst == +INF)
      { /* Oblique view from infinite distance: */
        if (R != NULL) { (*R) = rad; }
        if (D != NULL) { (*D) = 0; }
        if (G != NULL) { (*G) = 0; }
      }
    else
      { /* Oblique view from finite distance: */
        double v = sqrt((alfa - 1)*(alfa2 + 1)/beta2);
        double G = v/dst;
        if (debug)
          { Pr(Er, "  alfa2 = %16.12f  alfa = %16.12f\n", alfa2, sqrt(alfa2)); }
                /* Compute the stretch vector {*str}: */
                /* Not worth testing for {alfa2 == 1}: */
                r2_scale((sqrt(alfa2)- 1)*rad, &dir, str);
      }
  
ellipse_crs_t sphere_to_ellipse(pst_sphere_t *S, hr3_point_t *O)
  {
    ellipse_crs_t E;
    pst_sphere_and_viewpoint_to_ellipse_params
      ( S, O, &(E.ctr), &(E.rad), &(E.str) );
    return E;
  }

void pst_sphere_from_ellipse_A2(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S);
void pst_sphere_from_ellipse_A3(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S);
void pst_sphere_from_ellipse_A4(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S);
void pst_sphere_from_ellipse_A5(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S);

void pst_sphere_from_ellipse_B1(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S);
void pst_sphere_from_ellipse_B2(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S);
void pst_sphere_from_ellipse_B3(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S);
void pst_sphere_from_ellipse_B4(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S);
void pst_sphere_from_ellipse_B6(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S);
  /* The subcases of {pst_sphere_from_ellipse}. */

void pst_sphere_from_ellipse
  ( ellipse_crs_t *E, /* Projection of a sphere. */
    hr3_point_t *O,   /* (IN) Incomplete viewpoint (OUT) Inferred viewpoint. */
    pst_sphere_t *S   /* (OUT) Reconstructed sphere. */
  )
  {
    /* Paranoia: */
    demand(E->rad > 0, "bad {E.rad}");
    demand((! isnan(O->hm)) && (O->hm >= 0), "bad {O.m}");
    if (! isnan(O->hz)) { demand(O->hz > 0, "bad {O.z}"); }

    if (isnan(O->hx) && isnan(O->hy) && isnan(O->hz)) 
      { /* Case A5 or B5; push {O} to infinity and handle as A5: */
        O->hm = 0.0; 
        /* Since {Z} has to be positive, we may as well: */
        O->hz = 1.0;
      }
    
    /* Split by case: */
    if (O->hm == 0)
      { /* Case A: Viewpoint is at infinity. */
        
        if ((O->hx == 0) && (O->hy == 0))
          { /* Set {O} to {[0,0,0,1]} and handle in subcase A2: */
            O->hz = 1.0;
          }
            
        if (! (isnan(O->hx) || isnan(O->hy) || isnan(O->hz)))
          { /* Subcase A2: {O} infinite at known direction: */
            pst_sphere_from_ellipse_A2(E, O, S);
            return;
          }
        
        if 
          ( ((! isnan(O->hx)) || (O->hy == 0)) && 
            ((! isnan(O->hy)) || (O->hx == 0))
          )
          { /* Subcase A3: cylindrical from known azimuth, orthogonal in limit. */
            pst_sphere_from_ellipse_A3(E, O, S);
            return;
          }
        
        if ((isnan(O->hx) && isnan(O->hy)) || isnan(O->hz))
          { /* Subcase A5: cylindrical from unknown direction. */
            pst_sphere_from_ellipse_A5(E, O, S);
            return;
          }
          
        if ((! isnan(O->hz)) && (isnan(O->hx) != isnan(O->hy)))
          { /* Subcase A4: cylindrical from source on known infinite line. */
            pst_sphere_from_ellipse_A4(E, O, S);
            return;
          }
          
        /* We should have covered all 'A' cases: */
        assert(FALSE); 
      }
    else
      { /* Case B: finite viewpoint. */
        
        if ((! isnan(O->hx)) && (! isnan(O->hy)) && (! isnan(O->hz)))
          { /* Subcase B1: {O} is finite and known. */
            pst_sphere_from_ellipse_B1(E, O, S);
            return;
          }
          
        if (((! isnan(O->hx)) || (! isnan(O->hy))) && (! isnan(O->hz)))
          { /* Subcase B2: {O} is finite on a {X}- or {Y}-parallel line. */
            pst_sphere_from_ellipse_B2(E, O, S);
            return;
          }
        
        if ((isnan(O->hx) && isnan(O->hy)) && (! isnan(O->hz)))
          { /* Subcase B3: Known finite {F}, unknown {Q} (necessarily finite). */
            /* That is, {O} lies on an {X,Y}-parallel plane. */
            pst_sphere_from_ellipse_B3(E, O, S);
            return;
          }
        
        if ((isnan(O->hx) && isnan(O->hy)) && (! isnan(O->hz)))
          { /* Subcase B4: Known finite {Q}, unknown {F} (possibly infinite). */
            /* That is, {O} lies on a {Z}-parallel line. */
            pst_sphere_from_ellipse_B4(E, O, S);
            return;
          }
        
        /* Subcase B5 was handled as A5. */
        
        if ((isnan(O->hx) != isnan(O->hy)) && isnan(O->hz))
          { /* Subcase B6: {O} lies on a {X,Z}- or {Y,Z}-parallel plane. */
            pst_sphere_from_ellipse_B6(E, O, S);
            return;
          }
          
        /* We should have covered all 'B' cases: */
        assert(FALSE); 
      }
  }

void pst_sphere_from_ellipse_A2(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S)
  {
    /* Subcase A2: cyclindrical projection from known direction. */
    /* Compute the unit-length direction vector {N}: */
    r3_t N = (r3_t){{ O->hx, O->hy, O->hz }};
    double len = r3_dir(&N, &N); assert(len != 0);
    demand(N.Z > 0, "direction is too horizontal");
    /* Choose a radius {radA} of {E} to preserve the projected area from {N}: */
    double radE = E->rad;
    double majE = radE + r2_norm(&(E->str)); 
    double radA = sqrt(N.Z*majE*radE)
    /* The radius of {S} is the minor semidiameter of the projection: */
    S->R = radA;
    /* The center of {S} projects to the center of {E}: */
    S->K = E->ctr;
    /* Normalize {O} to unit Euclidean length: */
    O->hm = 0; O->hx = N.X; O->hy = N.Y; O->hz = N.Z;     
  }

void pst_sphere_from_ellipse_A3(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S)
  {
    /* Subcase A3: fixed azimuth {{H}, orthogonal in limit. */
    /* Get the 2D unit vector {H} of the azimuth: */
    r2_t H = (r2_t) {{ O->hx, O->hy }};
    /* If {H} is {(0,NAN)} or {(NAN,0)}, we can fix the {NAN} to 1: */ 
    if (isnan(H.X) && (! isnan(H.Y)) && (H.Y == 0)) { H.X = 1; }
    if (isnan(H.Y) && (! isnan(H.X)) && (H.X == 0)) { H.Y = 1; }
    double len = r2_dir(&H, &H); 
    assert(len > 0); /* Else we sould have gone to case A2. */
    /* Choose an aspect ratio {czA = rad/maj} for {E} from {E.rad,E.str,H}: */
    double TH = fabs(r2_dot(&(E->str), &H));
    double czA = E->rad/(E->rad + TH);
    /* Choose a radius {radA} for {E} that preserves area with aspect {czA}: */
    double radE = E->rad;
    double majE = radE + r2_norm(&(E->str));
    double radA = sqrt(czA*majE*radE)
    /* The radius of {S} is the minor semidiameter of the projection: */
    S->R = radA;
    /* The center of {S} projects to the center of {E}: */
    S->K = E->ctr;
    /* Compute the elevation sine {szA} from the aspect {czA}: */
    double szA = sqrt(fmax(0, 1 - czA*czA));
    /* Get the direction {N} of the viewpoint from {H} (or {-H}): */
    r3_t N = (r3_t){{ szA*H.X, szA*H.Y, czA }};
    (void)r3_dir(&N, &N); /* Just to be sure... */
    O->hm = 0; O->hx = N.X; O->hy = N.Y; O->hz = N.Z; 
  }

void pst_sphere_from_ellipse_A4(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S)
  {
    /* Subcase A4: cylindrical from source on known infinite line. */
    /* The center is preserved: */
    S->K = E->ctr;
    demand(FALSE, "!!! not implemented yet !!!");
  }

void pst_sphere_from_ellipse_A5(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S)
  {
    /* Subcase A5: cylindrical from unknown direction. */
    /* The radius is the minor semidiameter: */
    S->R = E->rad;
    /* The center projects to the center: */
    S->K = E->ctr;
    /* Get the azimuth {H} from the stretch of {E}: */
    r2_t H;
    double TE = r2_dir(&(E->str), &H);
    if (TE == 0) { H = (r2_t) {{ 0, 0 }}; }
    /* Compute the elevation trigs {czE,szE} from the aspect of {E}: */
    double czE = E->rad/(E->rad + TE);
    double szE = sqrt(fmax(0, 1 - czE*czE));
    /* Assemble the unit direction {N} of the viewpoint: */
    r3_t N = (r3_t){{ sz*H.X, sz*H.Y, cz }};
    (void)r3_dir(&N, &N); /* Just to be sure... */
    /* Set {O} normalized to unit Euclidean length: */
    O->hm = 0; O->hx = N.X; O->hy = N.Y; O->hz = N.Z; 
  }


void pst_sphere_from_ellipse_B1(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S)
  {
    /* Subcase B1: Conical, from finite known {O}. */
    /* Get the optical center {Q} and resolution {G}: */
    r3_t Q = (r2_t) {{ O->hx/O->hm, O->hy/O->hm }};
    double G = O->hz/O->hm;
    /* Get the displacement {dsp} from {Q} to {ctr}: */
    r2_t dsp; r2_sub(&(E->ctr), &Q, &dsp);
    /* Choose a radius {radA} for {E}: */
    /* !!! This is a hack, rethink: !!! */
    double radA = E->rad;
    /* Compute the sphere from {Q,G,dsp,radA}: */
    pst_sphere_compute_R_K_str_from_Q_G_dsp_rad
      ( &Q, G, &dsp, radA, &(S->K), &(S-R), NULL );
  }

void pst_sphere_from_ellipse_B2(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S)
  {
    /* Subcase B2: {O} is finite on a {X}- or {Y}-parallel line. */
    demand(FALSE, "!!! not implemented yet !!!");
  }

void pst_sphere_from_ellipse_B3(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S)
  {
    /* Subcase B3: Known finite {F}, unknown {Q} (necessarily finite). */
    /* That is, {O} lies on an {X,Y}-parallel plane. */
    /* Get the resolution {G}: */
    double G = O->hz/O->hm;
    /* Get the view azimuth {H} from the stretch of {E}: */
    r2_t H;
    double TE = r2_dir(&(E->str), &H);
    if (TE == 0) { H = (r2_t) {{ 0, 0 }}; }
    /* Get the major radius {maj} of {E}: */
    double maj = E->rad + TE;
    /* Compute the radius {S->R} and {dst = |Q - E.ctr|} from {G,rad,maj}: */
    double dst;
    pst_sphere_compute_R_dst_from_G_rad_maj(G,E->rad,maj,&(S-R),&D);
    /* Get the optical center {Q} from {H} and the aspect of {E}: */
    r3_t Q = (r2_t) {{ O->hx/O->hm, O->hy/O->hm }};
    /* Compute the sphere from {Q,G,dsp,E-rad}: */
    pst_sphere_compute_R_K_str_from_Q_G_dsp_rad
      ( &Q, G, &dsp, E->rad, &(S->K), &(S-R), NULL );
    /* Compute the 
    demand(FALSE, "!!! not implemented yet !!!");
  }

void pst_sphere_from_ellipse_B4(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S)
  {
    /* Subcase B4: Known finite {Q}, unknown {F} (possibly infinite). */
    /* That is, {O} lies on a {Z}-parallel line. */
    /* Get the optical center {Q}: */
    r2_t Q = (r2_t) {{ O->hx/O->hm, O->hy/O->hm }};
    /* Get the displacement {dsp} from {Q} to {ctr}: */
    r2_t dsp; r2_sub(&(E->ctr), &Q, &dsp);
    /* Compute the major radius {maj} of the ellipse: */
    double maj = E->rad + r2_norm(&(E->str));
    /* Compute {S.R}, {S.K}, and the resolution {G}: */
    double G = NAN;
    pst_sphere_compute_R_K_G_from_Q_dsp_rad_maj
      ( &Q, &dsp, &(E->rad), maj, &(S-R), &(S->K), &G ); 
    assert(isfinite(G)); /* Note: may be zero or very small. */
    /* Update {O} from {Q,G}: */
    double F = fabs(1/G);
    if (F > 9.9e+99) 
      { /* Practically at infinity, so put it there: */
        (*O) = (hr3_point) {{ 0, 0, 0, 1 }};
      }
    else
      { /* Finite enough: */
        (*O) = (hr3_point) {{ 1, Q.X, Q.Y, F }}; 
      }
  }

void pst_sphere_from_ellipse_B6(ellipse_crs_t *E, hr3_point_t *O, pst_sphere_t *S)
  {
    /* Subcase B6: {O} lies on a {X,Z}- or {Y,Z}-parallel plane. */
    demand(FALSE, "!!! not implemented yet !!!");
  }

/* INTERNAL TOOLS */

void pst_sphere_compute_stretch_from_apparent_params
  ( hr3_point_t *O,    /* Camera's viepoint. */
    r2_t *ctr,         /* The center of the sphere's projection on the image. */
    double rad,        /* The minor semidiameter of the projection. */
    r2_t *str          /* (OUT) The major semidiameter of the projection. */
  );
  /* Computes the expected stretch vector {*str} of a sphere's
    projection, given the center {ctr} of the projection and the
    transverse radius {rad} (minor semidiameter) of the projection.
    Note that {ctr} and {rad} are apparent parameters of the
    projection, not the nominal parameters of the sphere. */

/* SPHERE PROJECTION GEOMETRY */

void pst_geom_sphere_compute_projection_from_nominal_params
  ( hr3_point_t *O, /* Camera projection. */
    r2_t *K,         /* The nominal center of the sphere. */
    double R,        /* The radius of the sphere. */
    ellipse_crs_t *EP  /* (OUT) The sphere's projection on the image. */
  )
  { pst_geom_sphere_compute_apparent_params_from_nominal_params
      ( C, K, R, &(EP->ctr), &(EP->rad), &(EP->str) );
  }  

void pst_geom_sphere_compute_stretch_from_apparent_params
  ( hr3_point_t *O , /* Camera projection. */
    r2_t *ctr,   /* The center of the sphere's projection on the image. */
    double rad,  /* The minor semidiameter of the projection. */
    r2_t *str   /* (OUT) The major semidiameter of the projection. */
  )
  { 
    pst_geom_sphere_compute_nominal_params_and_stretch_from_apparent_params
      ( C, ctr, rad, /*K:*/ NULL, /*R:*/ NULL, /*str:*/ str );
  }
  
/* MAPPING BETWEEN 2D AND 3D COORDINATES */

void pst_geom_sphere_perspective_compute_point_and_normal
  ( hr3_pmap_t *M, /* NCS-to-VCS coordinate mapping. */
    r2_t *pim,     /* Image coords of the projection {P'} of a sphere point {P}. */
    r2_t *K,     /* Image coords of nominal center {K} of the sphere. */
    r3_t *psp,     /* (OUT) NCS coordinates of the point {P} on the sphere. */
    r3_t *nrm      /* (OUT) NCS coordinates of the sphere's normal at {P}. */
  )
  { 
    /* Outline of computation:

        (1) Assemble the homogeneous NCS coordinates of {P'},
        namely {[m',x',y',z'] = [1,X',Y',0]};

        (2) Multiply that 4-vector by the map M, to obtain the 
        homogeneous VCS coordinates {[t',u',v',w']} of {P'};

        (3) Compute the Cartesian VCS coordinates {U'} and {V'} of
        that point, that is, {U' = u'/t'} and {V' = V'/T'}.
        
        (4) If {U'^2 + V'^2 > 1}, then {P'} was outside
        the sphere's projection.  Set {P} and {N} to all {NAN}s and
        return.

        (5) Compute the Cartesian VCS coordinates of {P}
        {(U,V,W) = (U', V', sqrt(1 - U^2 - V^2)}.

        (6) Convert the Cartesian VCS coordinates {(U,V,W)}
        to the homogeneous VCS coordinates {[t,u,v,w] = [1,U,V,W]};

        (7) Multiply that 4-vector by the inverse of map {M},
        to obtain the homogeneous NCS coordinates {[m,x,y,z]} of {P};

        (8) Compute the homogeneous Cartesian coordinates of {P} as
        {(X,Y,Z) = (x/m,y/m,z/m)}.

        (9) Compute the normal direction {N} by subtracting the NCS
        coordinates of {K} from those of {P}, and normalizing to unit
        length. */
      
    /* (1): */
    r4_t hpim_i = (r4_t){{ 1, pim->X, pim->Y, 0 }}; 
    /* (2): */
    r4_t hpim_v; r4x4_map_row(&hpim_i, &(M->dir), &hpim_v); 
    /* (3): */
    double t1 = hpim_v.c[0], u1 = hpim_v.c[1], v1 = hpim_v.c[2];
    double U1 = u1/t1, V1 = v1/t1;    
    /* (4): */
    double rUV2 = U1*U1 + V1*V1;
    if (rUV2 > 1.0) 
      { (*psp) = (*nrm) = (r3_t) {{ NAN, NAN, NAN }}; return; } 
    /* (5): */
    double U = u1, V = V1, W = sqrt(fmax(0, 1 - U*U - V*V));            
    /* (6): */
    r4_t hpst_v = (r4_t){{ 1, U, V, W }};                               
    /* (7): */ 
    r4_t hpst_i; r4x4_map_row(&hpst_v, &(M->inv), &hpst_i);               
    /* (8): */
    double m =  hpst_i.c[0]; 
    (*psp) = (r3_t) {{ hpst_i.c[1]/m, hpst_i.c[2]/m, hpst_i.c[3]/m }};  
    /* (9): */
    r3_t pds = (r3_t) {{ psp->X - K->X, psp->Y - K->Y, psp->Z }};  
    r3_dir(&pds, nrm);                                                  
  }

hr3_pmap_t pst_geom_sphere_perspective_view_map(hr3_point_t *O , ellipse_crs_t *E)
  { /* Get the sphere's nominal parameters: */
    r2_t K;
    double R;
    r2_t str;
    pst_geom_sphere_compute_nominal_params_and_stretch_from_apparent_params
      ( C, &(E->ctr), E->rad, &K, &R, &str ); 
  
    /* Compute the sphere center to camera direction {wdr} and distance {D}: */
    r3_t wdr; double D;
    double m = O->hm;
    if (m == 0)
      { /* !!! Fix after changing the camera's viewpoint to be {hr3_point_t}. */
        double x = O->hx;
        double y = O->hy;
        double z = O->hz;
        wdr = (r3_t) {{ x, y, z }};
        double dd = r3_dir(&wdr, &wdr);
        demand(dd > 0, "invalid camera viewpoint");
        D = +INF;
      }
    else
      { /* Get the Cartesian NCS coordinates {O} of the viewpoint: */
        r3_t OC = r3_from_hr3(&(C->O));
        wdr = (r3_t) {{ OC.X - K.X, OC.Y - K.Y, OC.Z }}; 
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
    

