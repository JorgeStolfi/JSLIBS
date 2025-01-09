/* See pst_geom.h */
/* Last edited on 2025-01-01 15:10:45 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <values.h>
#include <stdlib.h>
#include <assert.h>

#include <float_image.h>
#include <r2.h>
#include <r3.h> 
#include <r3x3.h> 
#include <jsrandom.h>
#include <affirm.h>

#include <pst_geom.h>
#include <pst_argparser.h>

void pst_geom_clip_dir(r3_t *udir, r3_t *sdir, double rad)      
  { double ucos = r3_dot(sdir, udir);
    if (ucos < cos(rad))
      { r3_t para, perp; /* Comps. of {udir} parallel and perpendic. to {sdir}. */
        r3_decomp(udir, sdir, &para, &perp);
        while (r3_L_inf_norm(&perp) == 0.0)
          { /* Set {perp} to a random direction perpendicular to {sdir}. */
            r3_throw_ball(&perp);
            r3_decomp(&perp, sdir, &para, &perp);
          }
        /* Normalize {perp} to unit length. */
        (void)r3_dir(&perp, &perp);
        /* Set {udir} to point {rad} radians away from {sdir}: */
        double cr = cos(rad), sr = sin(rad);
        r3_mix(cr, sdir, sr, &perp, udir);
      }
  }

r3_t pst_geom_sphere_compute_normal(r2_t *uv)
  { double u = uv->c[0];
    double v = uv->c[1];
    /* Check if inside unit {u,v} disk: */
    double uv2 = u*u + v*v;
    if (uv2 >= 0.999999)
      { /* Outside - return null normal: */
        return (r3_t) {{ 0.0, 0.0, 0.0 }};
      }
    else
      { /* Inside - compute {w} coordinate: */
        double w = sqrt(1.000001 - uv2);
        /* The normal is the points's position: */
        return (r3_t) {{ u, v, w }};
      }
  }

void pst_geom_sphere_view_matrices
  ( ellipse_crs_t *geo, /* Geometry of sphere's projection. */
    r3x3_t *xym_to_uvm, 
    r3x3_t *uvw_to_xyz
  )
  { /* Initialize {xym_to_uvm} with translation matrix that subtracts {ctr}: */
    r3x3_ident(xym_to_uvm);
    xym_to_uvm->c[2][0] = -geo->ctr.c[0];
    xym_to_uvm->c[2][1] = -geo->ctr.c[1];
    
    /* Initialize {uvw_to_xyz} with the identity: */
    r3x3_ident(uvw_to_xyz);
    /* Compute min and max radii: */
    double mstr = r2_norm(&(geo->str));
    double rmin = geo->rad;
    double rmax = geo->rad + mstr;
    if (rmin != rmax) 
      { /* Elliptic shape, must change basis: */
        /* Get directions {du,dv} of shortest and longest axes: */
        r2_t dmin, dmax;     
        (void)r2_dir(&(geo->str), &dmax);
        dmin = (r2_t){{ +dmax.c[1], -dmax.c[0] }};
        /* Append to {xym_to_uvm} a rotation that maps {x,y} to {u,v} */
        r3x3_t rot;
        r3x3_ident(&rot);
        rot.c[0][0] = dmin.c[0];
        rot.c[1][0] = dmin.c[1];
        rot.c[0][1] = dmax.c[0];
        rot.c[1][1] = dmax.c[1];
        r3x3_mul(xym_to_uvm, &rot, xym_to_uvm);
        /* Set {uvw_to_urz} to rotate around {u} axis, changing {u,v,w} into {u,r,z}: */
        r3x3_t uvw_to_urz;
        r3x3_ident(&uvw_to_urz);
        double ct = rmin/rmax;         /* Cosine of rotation angle. */
        double st = sqrt(1.0 - ct*ct); /* Sine of rotation angle. */
        uvw_to_urz.c[1][1] = +ct;
        uvw_to_urz.c[1][2] = +st;
        uvw_to_urz.c[2][1] = -st;
        uvw_to_urz.c[2][2] = +ct;
        /* Append to {uvw_to_urz} the inverse rotation from {u,r,z} to {x,y,z}: */
        r3x3_t tor;
        r3x3_ident(&tor);
        tor.c[0][0] = dmin.c[0];
        tor.c[0][1] = dmin.c[1];
        tor.c[1][0] = dmax.c[0];
        tor.c[1][1] = dmax.c[1];
        r3x3_mul(&uvw_to_urz, &tor, uvw_to_xyz);
      }
    /* Append to {xym_to_uvm} an unscale by {rmin,rmax}: */
    r3x3_t unrad;
    r3x3_ident(&unrad);
    unrad.c[0][0] = 1/rmin;
    unrad.c[1][1] = 1/rmax;
    r3x3_mul(xym_to_uvm, &unrad, xym_to_uvm);
    
    /* fprintf(stderr, "xym_to_uvm:\n"); */
    /* r3x3_print(stderr, xym_to_uvm); */
    /* fprintf(stderr, "\n"); */
    /* fprintf(stderr, "uvw_to_xyz:\n"); */
    /* r3x3_print(stderr, uvw_to_xyz); */
    /* fprintf(stderr, "\n"); */
  }

ellipse_crs_t *pst_geom_sphere_parse
  ( argparser_t *pp,
    bool_t next,
    r2_t *ctrdef, 
    double *ctrAdj, 
    double *radAdj, 
    double *strAdj,
    double *adjStep
  )
  { 
    /* If there is no "-sphere" keyword, return NULL: */
    if (! pst_keyword_present(pp, "-sphere", next)) { return NULL; }
  
    /* Keyword is present; create and return a new descriptor: */
    ellipse_crs_t *geo = pst_geom_sphere_new();
    bool_t center_given = FALSE;
    bool_t radius_given = FALSE;
    bool_t stretch_given = FALSE;
    if (ctrAdj != NULL) { *ctrAdj = 0.0; }
    if (radAdj != NULL) { *radAdj = 0.0; }
    if (strAdj != NULL) { *strAdj = 0.0; }
    if (adjStep != NULL) { *adjStep = +INF; }
    while (TRUE)
      { if (argparser_keyword_present_next(pp, "center"))
          { if (center_given) { argparser_error(pp, "duplicate sphere's \"center\""); }
            geo->ctr.c[0] = argparser_get_next_double(pp, 0.0, 1.0e+20);
            geo->ctr.c[1] = argparser_get_next_double(pp, 0.0, 1.0e+20);
            if ((ctrAdj != NULL) && argparser_keyword_present_next(pp, "adjust"))
              { *ctrAdj = argparser_get_next_double(pp, 0.0, +DBL_MAX); }
            center_given = TRUE;
          }
        else if (argparser_keyword_present_next(pp, "radius"))
          { if (radius_given) { argparser_error(pp, "duplicate sphere's \"radius\""); }
            geo->rad = argparser_get_next_double(pp, 0.01, 1.0e+20);
            if ((radAdj != NULL) && argparser_keyword_present_next(pp, "adjust"))
              { *radAdj = argparser_get_next_double(pp, 0.0, +DBL_MAX); }
            radius_given = TRUE;
          }
        else if (argparser_keyword_present_next(pp, "stretch"))
          { if (stretch_given) { argparser_error(pp, "duplicate sphere's \"stretch\""); }
            geo->str.c[0] = argparser_get_next_double(pp, -1.0e+10, +1.0e+10);
            geo->str.c[1] = argparser_get_next_double(pp, -1.0e+10, +1.0e+10);
            if ((strAdj != NULL) && argparser_keyword_present_next(pp, "adjust"))
              { *strAdj = argparser_get_next_double(pp, 0.0, +DBL_MAX); }
            stretch_given = TRUE;
          }
        else
          { /* No valid parameter keyword next -- assume end of sphere's geometry spec: */
            break;
          }
      }
    if ((adjStep != NULL) && argparser_keyword_present_next(pp, "step"))
      { *adjStep = argparser_get_next_double(pp, 1.0e-20, 1.0e+20); }

    /* Provide defaults: */
    if (! center_given)
      { if (ctrdef == NULL)
          { argparser_error(pp, "must specify sphere's \"center\""); }
        else
          { geo->ctr = *ctrdef; }
      }

    if (! radius_given)
      { argparser_error(pp, "must specify sphere's \"radius\""); }

    if (! stretch_given) 
      { geo->str = (r2_t){{ 0.0, 0.0 }}; }

    return geo;
  }

void pst_geom_sphere_write(FILE *wr, ellipse_crs_t *geo)
  { r2_t *ctr = &(geo->ctr);
    double rad = geo->rad;
    r2_t *str = &(geo->str);
    fprintf(wr, "  -sphere\n");
    fprintf(wr, "    center %7.2f %7.2f\n", ctr->c[0], ctr->c[1]);
    fprintf(wr, "    radius %7.2f\n", rad);
    fprintf(wr, "    stretch %7.2f %7.2f\n", str->c[0], str->c[1]);
    fflush(wr);
  }
