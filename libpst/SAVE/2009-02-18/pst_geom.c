/* See pst_geom.h */
/* Last edited on 2009-02-28 19:44:13 by stolfi */

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

/* SPHERE PROJECTION GEOMETRY */

void pst_geom_sphere_compute_apparent_params_from_nominal_params
  ( pst_camera_t *C, /* Camera projection. */
    r2_t *tct,       /* The nominal center of the sphere. */
    double trd,      /* The nominal radius of the sphere. */
    r2_t *ctrP,      /* (OUT) The center of the sphere's projection on the image. */
    double *radP,    /* (OUT) The minor semidiameter of that projection. */
    r2_t *strP       /* (OUT) The stretch vector of the projection. */
  )
  {
    assert(trd > 0);
    double m = C->vpt.c.c[0];
    if (m == 0)
      { /* Parallel projection (a rather common case): */
        if (radP != NULL) { (*radP) = trd; }
        if (ctrP != NULL) { (*ctrP) = *tct; }
        if (strP != NULL) { (*strP) = (r2_t) {{ 0, 0 }}; }
        return;
      }
    /* Get the center and spread of the camera: */
    double spr = pst_camera_spread(C);
    r2_t opc = pst_camera_center(C);
    /* Compute {beta2}, the minor radius over the nominal radius, squared: */
    double t = trd*spr, t2 = t*t;
    demand(t < 1, "camera is lower than top of sphere");
    double beta2 = 1/(1 - t2); 
    /* assert(beta2 >= 1); */
    /* Compute the apparent transverse radius (minor semidiameter): */
    double rad = trd*sqrt(beta2); /* Worth computing, for {*radP} and {*strP} */
    if (radP != NULL) { (*radP) = rad; }
    if ((ctrP != NULL) || (strP != NULL))
      { /* Get the displacement {Dsp} from {opc} to {tct}: */
        r2_t Dsp; r2_sub(tct, &opc, &Dsp);
        if (ctrP != NULL) 
          { /* Compute the apparent center {*ctrP}: */
            /* Not worth testing for {beta2 == 1}: */
            r2_mix(1.0, &opc, beta2, &Dsp, ctrP);
          }
        if (strP != NULL)
          { /* Compute {D}, the distance from {opc} to {tct}: */
            double D = r2_norm(&Dsp);
            /* Compute {alfa2}, the major radius over the minor radius, squared: */
            double q = D*spr, q2 = q*q;
            double alfa2 = 1 + beta2*q2;
            /* Compute the stretch vector {*strP}: */
            /* Not worth testing for {alfa2 == 1}: */
            r2_scale((sqrt(alfa2)- 1)*rad/D, &Dsp, strP);
          }
      }
  }

void pst_geom_sphere_compute_projection_from_nominal_params
  ( pst_camera_t *C, /* Camera projection. */
    r2_t *tct,         /* The nominal center of the sphere. */
    double trd,        /* The nominal radius of the sphere. */
    ellipse_crs_t *EP  /* (OUT) The sphere's projection on the image. */
  )
  { pst_geom_sphere_compute_apparent_params_from_nominal_params
      ( C, tct, trd, &(EP->ctr), &(EP->rad), &(EP->str) );
  }  

void pst_geom_sphere_compute_nominal_params_and_stretch_from_apparent_params
  ( pst_camera_t *C, /* Camera projection. */
    r2_t *ctr,    /* The center of the sphere's projection on the image. */
    double rad,   /* The minor semidiameter of that projection. */
    r2_t *tctP,   /* (OUT) The nominal center of the sphere. */
    double *trdP, /* (OUT) The nominal radius of the sphere. */
    r2_t *strP    /* (OUT) The correct stetch vector of the projection. */
  )
  {
    bool_t debug = TRUE;
    assert(rad > 0);
    double m = C->vpt.c.c[0];
    if (m == 0)
      { /* Parallel projection (a rather common case): */
        if (trdP != NULL) { (*trdP) = rad; }
        if (tctP != NULL) { (*tctP) = *ctr; }
        if (strP != NULL) { (*strP) = (r2_t) {{ 0, 0 }}; }
        return;
      }
    /* Get the center and spread of the camera: */
    double spr = pst_camera_spread(C);
    r2_t opc = pst_camera_center(C);
    /* Compute the the ratio {u=rad/F} where {F} is the focal length: */
    double u = rad*spr, u2 = u*u;
    /* Compute {beta2}, the minor radius over the nominal radius, squared: */
    double beta2 = 1 + u2;
    if (debug) { Pr(Er, "  beta2 = %16.12f  beta = %16.12f\n", beta2, sqrt(beta2)); }
    /* assert(beta2 >= 1); */
    /* Compute the nominal radius {*trdP}: */
    if (trdP != NULL) { (*trdP) = rad/sqrt(beta2); }
    if ((tctP != NULL) || (strP != NULL))
      { /* Get the displacement {dsp} from {opc} to {ctr}: */
        r2_t dsp; r2_sub(ctr, &opc, &dsp);
        if (tctP != NULL) 
          { /* Compute the nominal center {*tctP}: */
            /* Not worth testing for {beta2 == 1}: */
            r2_mix(1.0, &opc, 1/beta2, &dsp, tctP);
          }
        if (strP != NULL)
          { /* Get the direction {dir} and distance {d} from {opc} to {ctr}: */
            r2_t dir;
            double d = r2_dir(&dsp, &dir);
            if (d == 0)
              { /* Sphere on optical axis: */
                (*strP) = (r2_t) {{ 0, 0 }};
              }
            else
              { /* Compute {alfa2}, the major radius over the minor radius, squared: */
                double v = d*spr, v2 = v*v;
                double alfa2 = 1 + beta2*v2;
                if (debug)
                  { Pr(Er, "  alfa2 = %16.12f  alfa = %16.12f\n", alfa2, sqrt(alfa2)); }
                /* Compute the stretch vector {*strP}: */
                /* Not worth testing for {alfa2 == 1}: */
                r2_scale((sqrt(alfa2)- 1)*rad, &dir, strP);
              }
          }
      }
    if (debug) 
      { if (tctP != NULL)  { Pr(Er, "  tct = "); r2_print(Er, tctP); Pr(Er, "\n"); }
        if (trdP != NULL)  { Pr(Er, "  trd = %18.12f\n", (*trdP)); }
        if (strP != NULL)  { Pr(Er, "  str = "); r2_print(Er, strP); Pr(Er, "\n"); }
      }
  }

void pst_geom_sphere_compute_stretch_from_apparent_params
  ( pst_camera_t *C, /* Camera projection. */
    r2_t *ctr,   /* The center of the sphere's projection on the image. */
    double rad,  /* The minor semidiameter of the projection. */
    r2_t *strP   /* (OUT) The major semidiameter of the projection. */
  )
  { 
    pst_geom_sphere_compute_nominal_params_and_stretch_from_apparent_params
      ( C, ctr, rad, /*tctP:*/ NULL, /*trdP:*/ NULL, /*strP:*/ strP );
  }
  
/* MAPPING BETWEEN 2D AND 3D COORDINATES */

void pst_geom_sphere_perspective_compute_point_and_normal
  ( hr3_pmap_t *M, /* NCS-to-VCS coordinate mapping. */
    r2_t *pim,     /* Image coords of the projection {P'} of a sphere point {P}. */
    r2_t *tct,     /* Image coords of nominal center {K} of the sphere. */
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
    r3_t pds = (r3_t) {{ psp->X - tct->X, psp->Y - tct->Y, psp->Z }};  
    r3_dir(&pds, nrm);                                                  
  }

hr3_pmap_t pst_geom_sphere_perspective_view_map(pst_camera_t *C, ellipse_crs_t *E)
  { /* Get the sphere's nominal parameters: */
    r2_t tct;
    double trd;
    r2_t str;
    pst_geom_sphere_compute_nominal_params_and_stretch_from_apparent_params
      ( C, &(E->ctr), E->rad, &tct, &trd, &str ); 
  
    /* Compute the sphere center to camera direction {wdr} and distance {D}: */
    r3_t wdr; double D;
    double m = C->vpt.c.c[0];
    if (m == 0)
      { /* !!! Fix after changing the camera's viewpoint to be {hr3_point_t}. */
        double x = C->vpt.c.c[1];
        double y = C->vpt.c.c[2];
        double z = C->vpt.c.c[3];
        wdr = (r3_t) {{ x, y, z }};
        double dd = r3_dir(&wdr, &wdr);
        demand(dd > 0, "invalid camera viewpoint");
        D = +INF;
      }
    else
      { /* Get the Cartesian NCS coordinates {vpt} of the viewpoint: */
        r3_t vpt = r3_from_hr3(&(C->vpt));
        wdr = (r3_t) {{ vpt.X - tct.X, vpt.Y - tct.Y, vpt.Z }}; 
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
        M.inv.c[1][j] = (j == 0 ? 0 : trd*udr.c[j-1]);
        M.inv.c[2][j] = (j == 0 ? 0 : trd*vdr.c[j-1]);
        M.inv.c[3][j] = (j == 0 ? 0 : trd*wdr.c[j-1]);
      }

    if (D != +INF)
      { /* We must pre-multiply {M.inv} by a parallel-to-conical matrix {A}, */
        /* and post-multiply by a translation from {K'} to {K}. */
        
        /* Compute the distance {L} from {K} to {K'}: */
        double R = trd;                     /* Nominal radius of sphere: */
        double T = sqrt((D - R)*(D + R));   /* Dist from viewpoint to horizon. */
        double S = trd*(T/D);               /* Radius of horizon. */
        double L = trd*trd/T;               /* Dist from sphere ctr to horizon ctr. */

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
    
/* UTILITIES */

void pst_geom_clip_dir(r3_t *udir, r3_t *sdir, double ard)      
  { double ucos = r3_dot(sdir, udir);
    if (ucos < cos(ard))
      { r3_t para, perp; /* Comps. of {udir} parallel and perpendic. to {sdir}. */
        r3_decomp(udir, sdir, &para, &perp);
        while (r3_L_inf_norm(&perp) == 0.0)
          { /* Set {perp} to a random direction perpendicular to {sdir}. */
            r3_throw_ball(&perp);
            r3_decomp(&perp, sdir, &para, &perp);
          }
        /* Normalize {perp} to unit length. */
        (void)r3_dir(&perp, &perp);
        /* Set {udir} to point {ard} radians away from {sdir}: */
        double cr = cos(ard), sr = sin(ard);
        r3_mix(cr, sdir, sr, &perp, udir);
      }
  }

ellipse_crs_t *pst_geom_sphere_new(void)
  { ellipse_crs_t *E = 
      (ellipse_crs_t *)notnull(malloc(sizeof(ellipse_crs_t)), "no mem");
    E->ctr = (r2_t){{ 0, 0 }};
    E->rad = 0;
    E->str = (r2_t){{ 0, 0 }};
    return E;
  }

