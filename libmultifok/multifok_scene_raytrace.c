/* See {multifok_scene_raytrace.h}. */
/* Last edited on 2024-10-25 07:16:12 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <interval.h>
#include <r3.h>
#include <frgb.h>
#include <r2.h>
#include <frgb_ops.h>
#include <jsrandom.h>

#include <multifok_scene.h>
#include <multifok_scene_raytrace.h>
#include <multifok_scene_tree.h>

#define ot_FLAT multifok_scene_object_type_FLAT
#define ot_RAMP multifok_scene_object_type_RAMP
#define ot_DISK multifok_scene_object_type_DISK
#define ot_BALL multifok_scene_object_type_BALL
 
#define ZMIN multifok_scene_ZMIN
#define ZMAX multifok_scene_ZMAX
  /* Shorter name. */

#define FUDGE (1.0e-6)
  /* Fudge amount to expand bounding boxes to account for roundoff. */

double multifok_scene_raytrace_tree
  ( multifok_scene_tree_t *tr, 
    r3_t *p, 
    r3_t *d, 
    double tMin,
    double tMax,
    bool_t verbose, 
    int32_t level,
    multifok_scene_object_t **oHit_P
  );
  /* Same as {multifok_scene_raytrace} but for a subtree {tr} of the whole
    scene tree, at depth {level}, with the ray clipped to the parameter
    range {[tMin _ tMax]}. */

double multifok_scene_raytrace_object
  ( multifok_scene_object_t *obj,
    r3_t *p,
    r3_t *d,
    double tMin,
    double tMax,
    bool_t verbose
  );
  /* Checks whether the ray {R} defined by point {p} and unit direction vector {d}
    intersects the upper surface of the object {obj} at 
    a parameter value {tHit} in the range {[tMin _ tMax]}. If it does,
    returns {tHit}, otherwise returns {+INF}. */

void multifok_scene_raytrace_get_ray_bbox
  ( r3_t *p,
    r3_t *d,
    double tMin,
    double tMax,
    interval_t bbox[]
  );
  /* Sets {bbox[0..2]} to the bounding box of the part of the ray
    with parameter value in the range {[tMin _ tMax]}. */
    
/* TRACING SPECIFI OBJECTS */
    
double multifok_scene_raytrace_FLAT(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose);   
  /* Like {multifok_scene_raytrace_object}, but specialized 
    for a {FLAT} floor object. Returns the ray's {t} parameter 
    if hit, {+INF} otherwise. */
      
double multifok_scene_raytrace_RAMP(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose);   
  /* Like {multifok_scene_raytrace_object}, but specialized for a {RAMP}
    object. Returns the ray's {t} parameter if hit, {+INF} otherwise. */
      
double multifok_scene_raytrace_DISK(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose);   
  /* Like {multifok_scene_raytrace_object}, but specialized for a {DISK}
    object. Returns the ray's {t} parameter if hit, {+INF} otherwise. */
      
double multifok_scene_raytrace_BALL(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose);   
  /* Like {multifok_scene_raytrace_object}, but specialized for a {BALL} 
    object. Returns the ray's {t} parameter if hit, {+INF} otherwise. */

void multifok_scene_raytrace
  ( multifok_scene_t *scene,
    multifok_scene_tree_t *tree,
    r3_t *p, 
    r3_t *d, 
    bool_t verbose,
    multifok_scene_object_t **oHit_P, 
    r3_t *pHit_P
  )
  {
    if (verbose || isnan(p->c[2]) || isnan(d->c[2])) 
      { r3_gen_print(stderr, p, "%12.8f", "        tracing ray p = ( ", " ", " )");
        r3_gen_print(stderr, d, "%12.8f", "  d = ( ", " ", " )\n");
      }
      
    demand(d->c[2] < 0, "invalid direction vector");
    demand(p->c[2] >= 0, "invalid focus plane {Z}");
   
    /* Assumes that the {Z} of the scene's visible surface is in {[ZMIN _ ZMAX]}: */
    double tMin = (ZMAX + 2*FUDGE - p->c[2])/d->c[2];
    double tMax = (ZMIN + 2*FUDGE - p->c[2])/d->c[2]; 
    assert(tMin < tMax);
    int32_t level = 0;
    double tHit = multifok_scene_raytrace_tree(tree, p, d, tMin, tMax, verbose, level, oHit_P);
    if (tHit != +INF)
      { assert((tHit >= tMin) && (tHit <= tMax));
        assert((*oHit_P) != NULL);
        double zHit = p->c[2] + tHit*d->c[2];
        double xHit = p->c[0] + tHit*d->c[0];
        double yHit = p->c[1] + tHit*d->c[1];
        (*pHit_P) = (r3_t) {{ xHit, yHit, zHit }};
      }
    else
      { (*pHit_P) = (r3_t) {{ NAN, NAN, NAN }}; }
  }

double multifok_scene_raytrace_tree
  ( multifok_scene_tree_t *tr, 
    r3_t *p, 
    r3_t *d, 
    double tMin,
    double tMax,
    bool_t verbose, 
    int32_t level,
    multifok_scene_object_t **oHit_P
  )
  { 
    multifok_scene_object_t *oHit = NULL;
    double tHit = +INF;
    
    if (tr != NULL)
      { if (verbose) { fprintf(stderr, "    %*sray tMin = %12.8f tMax = %12.8f\n", 2*level, "", tMin, tMax); }
        assert(tMin <= tMax);
        /* Get the ray's bounding box {ray_bbox[0..2]}: */
        interval_t ray_bbox[3];
        multifok_scene_raytrace_get_ray_bbox(p, d, tMin, tMax, ray_bbox);

        /* Check if it it intersects the tree's bbox: */
        bool_t ok = TRUE; /* Set to false if ray's bbox is disjoint from tree's bbox. */
        for (int32_t j = 0; (j < 3) & ok; j++) 
          { if (ray_bbox[j].end[1] < tr->bbox[j].end[0]) { ok = FALSE; }
            if (ray_bbox[j].end[0] > tr->bbox[j].end[1]) { ok = FALSE; }
          }
        if (ok)
          { /* Bounding boxes intersect, we have a chance. */
            /* Ray-trace the root object: */
            multifok_scene_object_t *oRoot = tr->obj;
            if (verbose) { fprintf(stderr, "        %*strying root obj %d\n", 2*level, "", oRoot->ID); }
            double tHit_root = multifok_scene_raytrace_object(oRoot, p, d, tMin, tMax, verbose);
            if (tHit_root != +INF)
              { /* Root object hit replaces previous hit: */
                assert((tHit_root >= tMin) && (tHit_root <= tMax));
                if (verbose) { fprintf(stderr, "        %*shit root obj %d at t = %12.8f tMax = %12.8f\n", 2*level, "", oRoot->ID, tHit_root, tMax); }
                /* Save object, bring down {tMax}: */
                oHit = oRoot;
                tHit = tHit_root;
                tMax = tHit_root;
              }
            else
              { if (verbose) { fprintf(stderr, "        %*smissed/rejected root obj %d\n", 2*level, "", oRoot->ID); } 
              }

            /* Ray trace the subtrees. We start with the subtree
              {sub[ic]} that seems closest to the ray in direction
              {axis}. Then we ray-trace {sub[ic]} and {sub[1-ic]},
              in that order, remembering the highest hit as we
              go. */
            int8_t axis = tr->axis;
            int32_t ic; /* First subtree to try. */
            if ((tr->sub[0] != NULL) && (tr->sub[1] != NULL))
              { double obc = interval_mid(&(tr->obj->bbox[axis]));
                double chd[2]; /* Distance to coord {axis} of child box centers. */
                for (int32_t jc = 0; jc < 2; jc++) 
                  { chd[jc] = fabs(obc - interval_mid(&(tr->sub[jc]->bbox[axis]))); }
                ic = (chd[0] < chd[1] ? 0 : 1);
              }
            else
              { /* No subtrees, {ic} is irrelevant. */
                ic = 0;
              }

            for (int32_t kc = 0; kc < 2; kc++) 
              { /* Ray-trace {sub[ic]} */
                if (verbose) { fprintf(stderr, "        %*strying subtree %d\n", 2*level, "", ic); }
                multifok_scene_object_t *oHit_sub;
                double tHit_sub = multifok_scene_raytrace_tree(tr->sub[ic], p, d, tMin, tMax, verbose, level+1, &oHit_sub);
                if (tHit_sub != +INF)
                  { /* We got a hit: */
                    assert((tHit_sub >= tMin) && (tHit_sub <= tMax));
                    assert(oHit_sub != NULL);
                    if (verbose) 
                      { fprintf(stderr, "        %*sgot hit with subtree %d obj %d at t = %12.8f tMax = %12.8f\n", 2*level, "", ic, oHit_sub->ID, tHit_sub, tMax); }
                    oHit = oHit_sub;
                    tHit = tHit_sub;
                    tMax = tHit_sub;
                  }
                else
                  { if (verbose) { fprintf(stderr, "        %*smissed hit with subtree %d\n", 2*level, "", ic); } 
                  }
                /* Try the other subtree: */
                ic = 1 - ic;
              }
          }
      }
    /* Return to caller: */
    if (verbose) 
      { if (oHit != NULL) 
          { fprintf(stderr, "        %*sreturning hit with obj %d at t = %12.8f\n", 2*level, "", oHit->ID, tHit);
          }
        else
          { fprintf(stderr, "        %*sreturning with no hit\n", 2*level, ""); }
      }
    (*oHit_P) = oHit;
    return tHit;
  }
        
void multifok_scene_raytrace_get_ray_bbox
  ( r3_t *p,
    r3_t *d,
    double tMin,
    double tMax,
    interval_t bbox[]
  )
  { 
    for (int32_t j = 0; j < 3; j++) 
      { double cLo = p->c[j] + d->c[j]*tMin; /* Coordinate at bottom of ray. */
        double cHi = p->c[j] + d->c[j]*tMax; /* Coordinate at top of ray. */
        bbox[j] = (interval_t){{ fmin(cLo, cHi) - 0.5*FUDGE, fmax(cLo, cHi) + 0.5*FUDGE }};
      }
  }

double multifok_scene_raytrace_object
  ( multifok_scene_object_t *obj,
    r3_t *p,
    r3_t *d,
    double tMin,
    double tMax,
    bool_t verbose
  )
  { 
    /* Get the ray's bounding box */
    interval_t ray_bbox[3];
    multifok_scene_raytrace_get_ray_bbox(p, d, tMin, tMax, ray_bbox);
    
    /* Check bbox intersection: */
    for (int32_t j = 0; j < 3; j++)
      { if (ray_bbox[j].end[1] < obj->bbox[j].end[0]) { return +INF; }
        if (ray_bbox[j].end[0] > obj->bbox[j].end[1]) { return +INF; }
      }
    /* Bounding boxes intersect.  Try th eobject itself: */
    double tHit = +INF;
    switch(obj->type)
      { case ot_FLAT: 
          tHit = multifok_scene_raytrace_FLAT(obj->bbox, p, d, verbose); break;
        case ot_RAMP:
          tHit = multifok_scene_raytrace_RAMP(obj->bbox, p, d, verbose); break;
        case ot_DISK:
          tHit = multifok_scene_raytrace_DISK(obj->bbox, p, d, verbose); break;
        case ot_BALL:
          tHit = multifok_scene_raytrace_BALL(obj->bbox, p, d, verbose); break;
        default: demand(FALSE, "unrecogneized object type");
      }
    /* If we still have a hit, check its {Z} against {zMIn}: */
    if (isfinite(tHit)) 
      { /* Ray hits object, but maybe below {zMin}: */
        char *typeX = multifok_scene_object_type_to_string(obj->type);
        if (verbose) { fprintf(stderr, "    hit %s object %d at t = %+12.8f", typeX, obj->ID, tHit); } 
        if (tHit < tMin)
          { /* Must be rounding error? */
            if (verbose) { fprintf(stderr, " (!! clipped, tMin = %12.8f)\n", tMin); }
            tHit = tMin;
          }
        else if (tHit > tMax)
          { if (verbose) { fprintf(stderr, " (rejected, tMax = %12.8f)\n", tMax); }
            tHit = +INF;
          }
        else
          { if (verbose) { fprintf(stderr, " (accepted)\n"); }
          }
      }
    return tHit;
  }

double multifok_scene_raytrace_DISK(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose)
  { double pX = p->c[0], pY = p->c[1], pZ = p->c[2];
    double dX = d->c[0], dY = d->c[1], dZ = d->c[2];
    demand(dZ < 0.0, "invalid direction vector");
    
    r3_t ctr = multifok_scene_box_center(bbox);
    double cX = ctr.c[0], cY = ctr.c[1], cZ = ctr.c[2];
    r3_t rad = multifok_scene_box_radius(bbox);
    double r_obj = fmin(rad.c[0], rad.c[1]);

    /* Compute the parameter {tHit} where {ray(tHit)} hits the disk's plane: */
    double tHit = (cZ - pZ)/dZ; 
    /* Compute horizontal displacement {sX,sY} from object's center to {ray(tHit)}: */ 
    double sX = pX + dX*tHit - cX; /* {X} position of ray hit rel to object ctr */
    double sY = pY + dY*tHit - cY; /* {Y} position of ray hit rel to object ctr */
    /* Check if ray hits object: */
    double r2_ray = sX*sX + sY*sY;
    if (r2_ray <= r_obj*r_obj)
      { return tHit; }
    else
      { return +INF; }
  }
        
double multifok_scene_raytrace_BALL(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose)
  { double pX = p->c[0], pY = p->c[1], pZ = p->c[2];
    double dX = d->c[0], dY = d->c[1], dZ = d->c[2];
    demand(dZ < 0.0, "invalid direction vector");
    r3_t ctr = multifok_scene_box_center(bbox);
    r3_t rad = multifok_scene_box_radius(bbox);
    double cX = ctr.c[0], cY = ctr.c[1], cZ = ctr.c[2];
    double r_obj = fmin(fmin(rad.c[0], rad.c[1]), rad.c[2]);

    /* Form the coefficients {A*t^2 + B*t + C} of */
    /* {d2_obj(t) = sX(t)^2 + sY(t)^2 + sZ(t)^2 - r2_obj} as function of {t}: */

    double LX = pX - cX;
    double LY = pY - cY;
    double LZ = pZ - cZ;

    /* sX(t) = pX + dX*t - cX  = LX + t*dX; */
    /* sY(t) = pY + dY*t - cY  = LY + t*dY; */
    /* sZ(t) = pZ + dZ*t - cZ  = LZ + t*dZ;    */

    double A = dX*dX + dY*dY + dZ*dZ;
    double B = 2*(LX*dX + LY*dY + LZ*dZ);
    double C = LX*LX + LY*LY + LZ*LZ - r_obj*r_obj;
    assert(A > 1.0e-12);

    double Delta = B*B - 4*A*C;
    if (Delta > 0) 
      { return (sqrt(Delta) - B)/(2*A); }
    else
      { return +INF; }
  }
        
double multifok_scene_raytrace_FLAT(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose)
  { 
    /* Get the ray parameters: */
    double pX = p->c[0], pY = p->c[1], pZ = p->c[2];
    double dX = d->c[0], dY = d->c[1], dZ = d->c[2];
    demand(dZ < 0, "invalid ray direction");

    /* Get the {Z} of the plane:: */
    double zObj = interval_mid(&(bbox[2])); /* {Z} of object. */

    /* Dertermine the ray time and coordinates of hit point: */
    double tHit = (zObj - pZ)/dZ;
    
    /* Check if hit inside the {XY} bounding box: */
    double xHit = pX + dX*tHit;
    if ((xHit < bbox[0].end[0]) || (xHit > bbox[0].end[1])) { return +INF; }
    double yHit = pY + dY*tHit;
    if ((yHit < bbox[1].end[0]) || (xHit > bbox[1].end[1])) { return +INF; }
    return tHit;
  }

double multifok_scene_raytrace_RAMP(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose)
  { 
    /* Get the ray parameters: */
    double pX = p->c[0], pY = p->c[1], pZ = p->c[2];
    double dX = d->c[0], dY = d->c[1], dZ = d->c[2];
    demand(dZ < 0, "invalid ray direction");

    /* Grab the {X} and {Z} ranges of the plane strip: */
    double bXlo = (2*bbox[0].end[0] + bbox[0].end[1])/3;
    double bXhi = (bbox[0].end[1] + 2*bbox[0].end[1])/3;
    double bZlo = bbox[2].end[0] + FUDGE; /* In case of roundoff errors. */
    double bZhi = bbox[2].end[1] - FUDGE; /* In case of roundoff errors. */

    /* Determine the equation {A*X + B*Z + C = 0} of the tilted */
    /* strip's plane, which goes through {(bXlo,0,bZlo)} and {(bXhi,0,bZhi)}: */
    double A = -(bZhi - bZlo);
    double B = +(bXhi - bXlo);
    double C = -(A*bXlo + B*bZlo);
    /* Substitute the ray's {X(t)} and {Z(t)} into that equation: */
    /* Namely, {X(t) = pX + dX*t}, {Z(t) = pZ + dZ*t} to get {AR*t + BR = 0}. */
    double AR = A*dX + B*dZ;
    double BR = A*pX + B*pZ + C;
    /* Get {tHit} of strip by solving the equation {AR*t + BR = 0}: */
    double tHit_ramp;
    if (fabs(AR) < 1.0e-12)
      { /* Ray is almost parallel to the strip: */
        tHit_ramp = +INF;
      }
    else
      { /* Ray is oblique to strip: */
        tHit_ramp = -BR/AR;
      }

    /* Get {tHit} of lower and upper flat parts: */
    double tHit_lo = (bZlo - pZ)/dZ;
    double tHit_hi = (bZhi - pZ)/dZ;
    assert(tHit_lo < tHit_hi);
    
    /* Select the proper hit: */
    double tHit = fmax(tHit_hi, fmin(tHit_lo, tHit_ramp));
    
    /* Check if hit inside the {XY} bounding box: */
    double xHit = pX + dX*tHit;
    if ((xHit < bbox[0].end[0]) || (xHit > bbox[0].end[1])) { return +INF; }
    double yHit = pY + dY*tHit;
    if ((yHit < bbox[1].end[0]) || (xHit > bbox[1].end[1])) { return +INF; }
    return tHit;
  }
 
/* SCENE TO IMAGE COORDINATE MAPPING */

r3_t multifok_scene_coords_to_image_coords
  ( r3_t *p_scene,
    double zFoc,
    interval_t scene_dom[],
    int32_t NX,
    int32_t NY
  )
  {
    double pixSize = multifok_scene_pixel_size(NX, NY, scene_dom);
    
    double x_scene_ctr = interval_mid(&(scene_dom[0]));
    double y_scene_ctr = interval_mid(&(scene_dom[1]));
    double x_img_ctr = 0.5*NX;
    double y_img_ctr = 0.5*NY;

    double x_img = y_img_ctr + (p_scene->c[0] - x_scene_ctr)/pixSize;
    double y_img = x_img_ctr + (p_scene->c[1] - y_scene_ctr)/pixSize;
    double z_img = (p_scene->c[2] - zFoc)/pixSize;
    return (r3_t){{ x_img, y_img, z_img }};    
  }
  
r3_t multifok_scene_coords_from_image_coords
  ( r3_t *p3_img,
    int32_t NX,
    int32_t NY,
    interval_t scene_dom[],
    double zFoc
  )
  {
    double pixSize = multifok_scene_pixel_size(NX, NY, scene_dom);

    double x_scene_ctr = interval_mid(&(scene_dom[0]));
    double y_scene_ctr = interval_mid(&(scene_dom[1]));
    double x_img_ctr = 0.5*NX;
    double y_img_ctr = 0.5*NY;

    double x_scene = x_scene_ctr + (p3_img->c[0] - x_img_ctr)*pixSize;
    double y_scene = y_scene_ctr + (p3_img->c[1] - y_img_ctr)*pixSize;
    double z_scene = zFoc + p3_img->c[2]*pixSize;
    return (r3_t){{ x_scene, y_scene, z_scene }};
  }
    
double multifok_scene_pixel_size(int32_t NX, int32_t NY, interval_t scene_dom[])
  { demand((NX > 0) && (NY > 0), "invalid {NX,NY}");
    double WX = 2*interval_rad(&(scene_dom[0]));
    double WY = 2*interval_rad(&(scene_dom[1]));
    demand((WX > 0) && (WY > 0), "invalid scene domain");
    double pixSize = fmin(WX/((double)NX), WY/((double)NY));
    return pixSize;
  }

frgb_t multifok_scene_raytrace_compute_hit_color
  ( multifok_scene_object_t *obj,
    r3_t *q,
    multifok_scene_raytrace_pattern_t *pattern
  )
  { if (obj == NULL)
      { return (frgb_t){{ 0.500, 0.500, 0.500 }}; }
    else
      { interval_t *bbox = obj->bbox; 
        /* Express {q} relative to object: */
        r3_t u; 
        for (int32_t j = 0; j < 3; j++) 
          { u.c[j] = q->c[j] - interval_mid(&(bbox[j])); }
        int32_t ID = obj->ID;
        if (ID != -1)
          { /* Rotate {u} about the {Z} axis as a function of {ID}: */
            double ang = 3*ID + 1, ca = cos(ang), sa = sin(ang);
            double ux = + ca*u.c[0] + sa*u.c[1];
            double uy = - sa*u.c[0] + ca*u.c[1];
            u = (r3_t){{ ux, uy, u.c[2] }};
          }
        /* Evaluate the 3D pattern at {u}: */
        double r = pattern(u.c[0],u.c[1],u.c[1]);
        /* Use it to mix the object's {fg} and {bg} colors: */
        frgb_t clr = frgb_mix((1-r), &(obj->bg), r, &(obj->fg));
        return clr;
      }
  }

#define multifok_scene_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

