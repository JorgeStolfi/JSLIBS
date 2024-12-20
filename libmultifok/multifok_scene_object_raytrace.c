/* See {multifok_scene_object_raytrace.h}. */
/* Last edited on 2024-12-05 10:36:40 by stolfi */

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
#include <jsrandom.h>
#include <jsqroots.h>

#include <multifok_scene.h>
#include <multifok_scene_raytrace.h>
#include <multifok_scene_object.h>

#include <multifok_scene_object_raytrace.h>

#define ot_FLAT multifok_scene_object_type_FLAT
#define ot_RAMP multifok_scene_object_type_RAMP
#define ot_DISK multifok_scene_object_type_DISK
#define ot_BALL multifok_scene_object_type_BALL
#define ot_CONE multifok_scene_object_type_CONE

#define FUDGE (1.0e-6)
  /* Fudge amount to expand bounding boxes to account for roundoff. */

double multifok_scene_object_raytrace
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
    bool_t ok = TRUE;
    for (uint32_t j = 0;  (j < 3) && ok; j++)
      { if (ray_bbox[j].end[1] < obj->bbox[j].end[0]) { ok = FALSE; }
        if (ray_bbox[j].end[0] > obj->bbox[j].end[1]) { ok = FALSE; }
      }
    if (! ok)
      { if (verbose) 
          { fprintf(stderr, "    ray and object bboxes don't intersect");
            multifok_scene_print_box(stderr, " obj = ", obj->bbox, "");
            multifok_scene_print_box(stderr, " ray = ", ray_bbox, "\n");
          } 
        return +INF;
      }
    /* Bounding boxes intersect.  Try the object itself: */
    double tHit = +INF;
    switch(obj->type)
      { case ot_FLAT: 
          tHit = multifok_scene_object_raytrace_FLAT(obj->bbox, p, d, verbose); break;
        case ot_RAMP:
          tHit = multifok_scene_object_raytrace_RAMP(obj->bbox, p, d, verbose); break;
        case ot_DISK:
          tHit = multifok_scene_object_raytrace_DISK(obj->bbox, p, d, verbose); break;
        case ot_BALL:
          tHit = multifok_scene_object_raytrace_BALL(obj->bbox, p, d, verbose); break;
        case ot_CONE:
          tHit = multifok_scene_object_raytrace_CONE(obj->bbox, p, d, verbose); break;
        default: demand(FALSE, "unrecognized object type");
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

double multifok_scene_object_raytrace_DISK(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose)
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
        
double multifok_scene_object_raytrace_BALL(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose)
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
      { /* Return the SMALLER root: */
        return -(B + sqrt(Delta))/(2*A);
      }
    else
      { return +INF; }
  }
        
double multifok_scene_object_raytrace_CONE(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose)
  { /* Grab important cordinates: */
    double pX = p->c[0], pY = p->c[1], pZ = p->c[2];
    double dX = d->c[0], dY = d->c[1], dZ = d->c[2];
    demand(dZ < 0.0, "invalid direction vector");
    double cX = interval_mid(&(bbox[0]));
    double cY = interval_mid(&(bbox[1]));
    double loZ = bbox[2].end[0];
    double hiZ = bbox[2].end[1];
    
    double tMin = (hiZ - pZ)/dZ;
    double tMax = (loZ - pZ)/dZ;
    
    /* The cone equation is {(X-cX)^2 + (Y-cY)^2 = ((Z-hiZ)/2)^2}. */
    /* Substitute {X = pX+t*dX}, {Y=pY+t*dY}, {Z=pZ+t*DZ}. */
    /* {(t*dX+pX-cX)^2 + (t*dY+pY-cY)^2 = ((t*dZ+pZ-hiZ)/2)^2} */

    double LX = pX - cX;  
    double LY = pY - cY;  
    double LZ = pZ - hiZ; 
    
    /* Collect terms to get {A*t^2 + B^t + C = 0}: */

    double A = dX*dX + dY*dY - dZ*dZ/4;
    double B = 2*LX*dX + 2*LY*dY - LZ*dZ/2;
    double C = LX*LX + LY*LY - LZ*LZ/4;
    
    double r1, r2, im;
    int32_t sgn = roots_quadratic(A, B, C, &r1, &r2, &im);
    if (im != 0)
      { /* No real root: */ return +INF; }
    else if (isnan(r1))
      { /* Degenerate, no roots: */ return +INF; }
    else if (isnan(r2))
      { /* Only one root {r1}: */
        if ((r1 >= tMin) && (r1 <= tMax)) 
          { return r1; }
        else
          { return +INF; }
      }
    else
      { assert(r1 <= r2);
        assert(sgn > 0);
        if ((r1 >= tMin) && (r1 <= tMax)) 
          { return r1; }
        else if ((r2 >= tMin) && (r2 <= tMax)) 
          { return r2; }
        else
          { return +INF; }
      }
  }
        
double multifok_scene_object_raytrace_FLAT(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose)
  { 
    /* Get the ray parameters: */
    double pX = p->c[0], pY = p->c[1], pZ = p->c[2];
    double dX = d->c[0], dY = d->c[1], dZ = d->c[2];
    demand(dZ < 0, "invalid ray direction");

    /* Get the {Z} of the plane:: */
    double zObj = interval_mid(&(bbox[2])); /* {Z} of object. */

    /* Dertermine the ray time and {XY} coordinates of hit point: */
    double tHit = (zObj - pZ)/dZ;
    double xHit = pX + dX*tHit;
    double yHit = pY + dY*tHit;
    if (verbose) { fprintf(stderr, "    tHit = %.16e xHit = %.16e yHit = %.16e zObj = %.16e\n", tHit, xHit, yHit, zObj); }
    
    /* Check if hit inside the {XY} bounding box: */
    if ((xHit < bbox[0].end[0]) || (xHit > bbox[0].end[1])) { return +INF; }
    if ((yHit < bbox[1].end[0]) || (yHit > bbox[1].end[1])) { return +INF; }
    return tHit;
  }

double multifok_scene_object_raytrace_RAMP(interval_t bbox[], r3_t *p, r3_t *d, bool_t verbose)
  { 
    /* Get the ray parameters: */
    double pX = p->c[0], pY = p->c[1], pZ = p->c[2];
    double dX = d->c[0], dY = d->c[1], dZ = d->c[2];
    demand(dZ < 0, "invalid ray direction");

    /* Grab the {X} range of the plane strip:  */
    /* It should be a bit less than the scene's {X} range: */
    double bXrad = interval_rad(&(bbox[0]));
    double bXlo = bbox[0].end[0] + WD_RAMP*bXrad;
    double bXhi = bbox[0].end[1] - WD_RAMP*bXrad;

    /* Grab the {Z} range of the plane strip: */
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
    assert(tHit_lo > tHit_hi);  /* Note! Because {dZ < 0}. */
    
    /* Select the proper hit: */
    double tHit = fmax(tHit_hi, fmin(tHit_lo, tHit_ramp));
    
    /* Compute {XY} coordinates of hit point: */ 
    double xHit = pX + dX*tHit;
    double yHit = pY + dY*tHit;
    
    /* Check if hit inside the {XY} bounding box: */
    if ((xHit < bbox[0].end[0]) || (xHit > bbox[0].end[1])) { return +INF; }
    if ((yHit < bbox[1].end[0]) || (yHit > bbox[1].end[1])) { return +INF; }
    return tHit;
  }
