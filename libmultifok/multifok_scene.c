/* See {multifok_scene.h}. */
/* Last edited on 2023-01-22 19:10:11 by stolfi */

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
#include <r2.h>
#include <r3.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <jsrandom.h>
#include <jsqroots.h>

#include <multifok_scene.h>

int32_t multifok_scene_choose_obj_count(interval_t box[], double rMin, double rMax, double minSep);
  /* Chooses the ideal number of objects (disks or balls) that {multifok_scene_throw} 
    should try to put in a scene.  The parameter {minSep} has the meaning described 
    under {multifok_scene_throw}. */

void multifok_scene_throw_colors(frgb_t *tone, frgb_t *bg, frgb_t *fg);
  /* Generates two random contrasting colors, with general hue depending on {warm}. */

void multifok_scene_ray_trace
  ( multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    r3_t *p, 
    r3_t *d, 
    bool_t debug,
    frgb_t *ray_clr_P, 
    r3_t *ray_hit_P
  );
  /* Traces one ray {R} that goes through the point {p} with the
    direction parallel to {d}, assumed to be not horizontal. Finds the
    object (disk, spher, or backplane) with max {Z} that hits that ray.
    Returns its color {clr(R)} in {*ray_clr_P} and the hit point
    {hit(R)} in {*ray_hit_P}.
    
    The color {clr(R)} is the pixel {fs[0..2]} returned by
    {pattern(x,y,z,kd,3,fs)}. If {hit(R)} is on a disk or sphere with
    center {ctr}, {kd} is the index if the object that was hit, the
    pattern coordinates {(x,y,z)} are {hit(R)-ctr}. If {hit(R)} is on
    the backplane, {kd} is set to {-1}, {x} and {y} are the coords of
    {hit(R)}, and {z} is set to 0, ignoring the actual {Z} of the
    backplane. */

bool_t multifok_scene_ray_trace_object(multifok_scene_object_t *object, r3_t *p, r3_t *d, bool_t debug, r3_t *hit_P);
  /* Traces one ray {R} that goes through the point {p} with the direction parallel to {d}, assumed to be 
    not horizontal.  If the ray hits a disk or sphere, returns {TRUE} and stores the hit point in {*hit_P}.
    Otherwise returns {FALSE} and leaves {*hit_P} unchanged. */
    
void multifok_scene_ray_trace_backplane(interval_t box[], bool_t back_tilt, r3_t *p, r3_t *d, bool_t debug, r3_t *back_hit_P);   
  /* Traces one ray {R} that goes through the point {p} with the direction parallel to {d}, assumed to be 
    not horizontal.  Computes the intersection of the ray with the backplane and returns that point in {*back_hit_P}.
    The backplane is defined by the box {box[0..2]} and the flag {back_tilt}, as described under {multifok_scene_t}. */
    
/* IMPLEMENTATIONS */

multifok_scene_t *multifok_scene_throw
  ( interval_t box[],
    bool_t back_tilt,
    double rMin, 
    double rMax,
    double minSep, 
    bool_t verbose
  )
  {
    demand((box[2].end[0] >= 0.0) && (box[2].end[1] <= multifok_scene_ZMAX), "bad {box[2]}");
    srandom(4615*417);

    /* Generate the objects {obj[0..ND-1]}. If overlapping, every try is valid, otherwise we need many more tries:*/
    multifok_scene_object_t *obj = NULL;
    int32_t ND = 0; /* Number of objects actually generated. */
    if (! back_tilt)
      { int32_t ND_max = multifok_scene_choose_obj_count(box, rMin, rMax, minSep);
        if (verbose) { fprintf(stderr, "trying to generate %d objects\n", ND_max); }
    
        obj = (multifok_scene_object_t*)notnull(malloc(ND_max*sizeof(multifok_scene_object_t)), "no mem");
        int32_t NT = (minSep >= 0 ? 50 : 1)*ND_max; /* Number of tries. */
        for (int32_t kt = 0; (kt < NT) && (ND < ND_max); kt++)
          { /* Generate a random object {objk}: */
            multifok_scene_object_t objk = multifok_scene_object_throw(box, minSep, rMin, rMax);
            bool_t ok = TRUE; /* False if object overlaps the previous ones in {XY}. */
            if (minSep >= 0)
              { /* Check for {XY} overlaps: */
                for (int32_t jd = 0; (jd < ND) && ok; jd++)
                  { multifok_scene_object_t *dj = &(obj[jd]);
                    /* Check for overlap between {objk} and {dj} plus min separation: */
                    double dxy = hypot(objk.ctr.c[0] - dj->ctr.c[0], objk.ctr.c[1] - dj->ctr.c[1]);
                    if (dxy < objk.rad + minSep + dj->rad) { ok = FALSE; }
                  }
              }
            if (ok) 
              { obj[ND] = objk;
                if (verbose) 
                  { r3_gen_print(stderr, &(objk.ctr), "%12.6f", "  ctr = ( ", " ", " )");
                    fprintf(stderr, " rad = %12.6f", objk.rad);
                    fprintf(stderr, " %s\n", (objk.flat ? "disk" : "ball"));
                  }
                for (int32_t j = 0; j < 3; j++)
                  { if ((j == 2) || (minSep >= 0))
                      { /* Check for full containment in box: */
                        double radj = (objk.flat && (j == 2) ? 0.0 : objk.rad);
                        double xmin = objk.ctr.c[j] - radj;
                        double xmax = objk.ctr.c[j] + radj;
                        double bmin = box[j].end[0];
                        double bmax = box[j].end[1];
                        if (j != 2)
                          { bmin += minSep;
                            bmax -= minSep; 
                          }
                        affirm(xmin > bmin, "object to close to box lo");
                        affirm(xmax < bmax, "object to close to box hi");
                      }
                  }
                ND++;
              }
          }
        assert(ND <= ND_max);
        if (ND < ND_max)
          { if (verbose) { fprintf(stderr, "generate only %d objects\n", ND); }
            /* Trim array: */
            obj = (multifok_scene_object_t*)notnull(realloc(obj, ND*sizeof(multifok_scene_object_t)), "no mem");
          }
      }
      
    /* Choose the background colors: */
    frgb_t tone = (frgb_t){{ 0.0, 0.0, 0.0 }};
    frgb_t bg, fg;
    multifok_scene_throw_colors(&tone, &bg, &fg);
     
    /* Store data in scene: */
    multifok_scene_t *scene = (multifok_scene_t *)notnull(malloc(sizeof(multifok_scene_t)), "no mem");
    for (int32_t j = 0; j < 3; j++) { scene->box[j] = box[j]; }
    scene->back_tilt = back_tilt;
    scene->ND = ND;
    scene->obj = obj;
    scene->bg = bg;
    scene->fg = fg;
    return scene;
  }

int32_t multifok_scene_choose_obj_count(interval_t box[], double rMin, double rMax, double minSep)
  {
    /* Compute the average area {aObj} of an object with radius in {[rMin _ rMax]}, accounting for min sep: */
    double rfat = (minSep >= 0 ? 0.5*minSep : 0);
    double rr0 = rMin + rfat; /* Min radius object occup. */
    double rr1 = rMax + rfat; /* Max radius object occup. */
    double aObj = M_PI*(rr1*rr1*rr1 - rr0*rr0*rr0)/(rr1 - rr0)/3; /* Average disk area. */
    double rObj = sqrt(aObj/M_PI); /* RMS average object radius. */
    
    /* Compute the {XY} area {aBox} available for those disks: */
    double wd[2];
    for (int32_t j = 0; j < 2; j++)
      { wd[j] = box[j].end[1] - box[j].end[0]; 
        if (minSep >= 0)
          { wd[j] -= 2*minSep; }
        else
          { wd[j] += 2*rObj; }
        demand (wd[j] >= 0, "box too tight");
      }
    double aBox = wd[0]*wd[1]; /* Total useful area of box. */

    /* Decide max number of objects. If non-overlapping, limited by area ratio, else more than that: */
    int32_t ND_max = (minSep >= 0 ? 1 : 5)*(int32_t)ceil(aBox/aObj); /* Max number of objects. */
    return ND_max;
  }

multifok_scene_object_t multifok_scene_object_throw(interval_t box[], double minSep, double rMin, double rMax)
  {
    demand(rMin <= rMax, "invalid radius interval");

    /* Choose disk or ball: */
    bool_t flat = (drandom() < 0.5);
    
    /* Get the preliminary range of positions {[blo[0..2] _ bho[0..2]} for the object: */
    double blo[3], bhi[3];
    for (int32_t j = 0; j < 3; j++) 
      { blo[j] = box[j].end[0]; 
        bhi[j] = box[j].end[1];
        /* Shrink the box by a small amount: */
        double eps = 0.0001*(bhi[j] - blo[j]);
        blo[j] += eps;
        bhi[j] -= eps;
        /* If {minSep} is positive, shrink it by {minSep} in {X} and {Y}: */
        if ((minSep > 0) && ((j == 0) || (j == 1)))
          { blo[j] += minSep;
            bhi[j] -= minSep;
          }
      }

    /* Choose the object's radius: */
    double rad = rMin + drandom()*drandom()*(rMax - rMin);

    /* Reduce {rad} as needed so that the ball can be placed at all: */
    for (int32_t j = 0; j < 3; j++) 
      { double bw = bhi[j] - blo[j];
        if (((minSep >= 0) || (j == 2)) && (2*rad > bw))
          { rad = 0.4999*bw; }
      }

    /* Choose the center {ctr}: */
    r3_t ctr;
    for (int32_t j = 0; j < 3; j++)
      { double xlo = blo[j];
        double xhi = bhi[j];
        /* Decide the object's size along axis {j}: */
        double radj = ((j == 2) && flat ? 0.0 : rad);
        /* Decide if the object's {j}-projection must be fully inside {[blo[j] _ bhi[j]]}: */
        bool_t fully_inside = (j == 2) || (minSep >= 0);
        /* If it has to be fully inside, shrink the range for {ctr.c[j]}: */
        if (fully_inside) { xlo += radj; xhi -= radj; }
        demand (xlo < xhi, "box is too tight");
        /* Choose the center's coordinate: */
        ctr.c[j] = xlo + drandom()*(xhi - xlo);
      }
      
    /* Choose the background colors: */
    frgb_t tone_disk = (frgb_t){{ 1.0f, 0.3f, 0.0f }};
    frgb_t tone_ball = (frgb_t){{ 0.0f, 0.3f, 1.0f }};
    frgb_t bg, fg;
    multifok_scene_throw_colors((flat ? &tone_disk : &tone_ball), &bg, &fg);

    /* Assemble the object: */
    multifok_scene_object_t obj;
    obj.ctr = ctr;
    obj.rad = rad;
    obj.flat = flat;
    obj.bg = bg;
    obj.fg = fg;
    return obj;
  }
  
void multifok_scene_throw_colors(frgb_t *tone, frgb_t *bg, frgb_t *fg)
  {
    double wmax = +INF;
    for (int32_t ic = 0; ic < 3; ic++) 
      { double vi = tone->c[ic];
        bg->c[ic] = (float)(0.25*vi);
        wmax = fmin(wmax, 1.0 - bg->c[ic]);
      }
    for (int32_t ic = 0; ic < 3; ic++) 
      { fg->c[ic] = (float)(bg->c[ic] + wmax);
        assert(fg->c[ic] <= 1.00001);
        fg->c[ic] = (float)fmin(fg->c[ic], 1.0);
      }
  }

void multifok_scene_ray_trace
  ( multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    r3_t *p, 
    r3_t *d, 
    bool_t debug,
    frgb_t *ray_clr_P, 
    r3_t *ray_hit_P
  )
  {
    bool_t verbose = FALSE;
    
    if (verbose) 
      { r3_gen_print(stderr, p, "%12.6f", "  ray p = ( ", " ", " )");
        r3_gen_print(stderr, d, "%12.6f", "  d = ( ", " ", " )\n");
      }
    demand(d->c[2] == 1.0, "invalid direction vector");
    demand(p->c[2] >= 0, "invalid focus plane {Z}");
    
    int32_t ND = scene->ND;
    int32_t NC = 3; /* Colors have 3 components. */
    
    r3_t ray_hit = (r3_t){{ NAN, NAN, -INF }};  /* Highest point where ray hits scene. */
    
    /* Check ray hits with the objs: */
    int32_t kdMax = -1;                  /* Index of highest object hit, or {-1}. */
    for (int32_t kd = 0; kd < ND; kd++)
      { multifok_scene_object_t *obj = &(scene->obj[kd]);
        r3_t obj_hit;
        bool_t hit = multifok_scene_ray_trace_object(obj, p, d, debug, &obj_hit);
        if (hit)
          { /* We have a hit: */
            if (obj_hit.c[2] > ray_hit.c[2])
              { ray_hit = obj_hit;
                kdMax = kd;
              }
          }
      }
      
    /* Check ray hit with backplane: */
    r3_t back_hit;
    multifok_scene_ray_trace_backplane(scene->box, scene->back_tilt, p, d, debug, &back_hit);
 
    /* Decide if ray hit backplane or objects first; set {fg,bg} colors of hit point. */
    if ((kdMax < 0) || (back_hit.c[2] > ray_hit.c[2]))
      { /* Replace {ray_hit}: */
        ray_hit = back_hit;
        kdMax = -1; /* Hit object is no longer a object. */
      }
  
    /* Get hit point colors (backplane color is indep of {Z}, even if tilted): */
    double zPat = (kdMax >= 0 ? ray_hit.c[2] : 0.0);
    frgb_t ray_clr;
    pattern(ray_hit.c[0], ray_hit.c[1], zPat, kdMax, NC, ray_clr.c);

    (*ray_clr_P) = ray_clr;
    (*ray_hit_P) = ray_hit;
  }

bool_t multifok_scene_ray_trace_object(multifok_scene_object_t *obj, r3_t *p, r3_t *d, bool_t debug, r3_t *hit_P)
  { 
    bool_t verbose = debug;
    
    /* Grab the relevant coordinates: */
    double pX = p->c[0], pY = p->c[1], pZ = p->c[2];
    double dX = d->c[0], dY = d->c[1], dZ = d->c[2];
    demand(dZ = 1.0, "invalid direction vector");
    r3_t *c = &(obj->ctr);
    double cX = c->c[0], cY = c->c[1], cZ = c->c[2];
    double r_obj = obj->rad, r2_obj = r_obj*r_obj;
    
    /* Let {ray(Z) = p + d*(Z - pZ)} be the point on the ray at height {Z}. */
    /* For generic {Z}, let {eZ} be {Z - pZ}, then {ray(Z) = p + d*eZ}. */
    /* For generic {Z}, let {s} be {ray(Z) - c}, that is {(p - c) + d*eZ}. */
        
    bool_t hit;     /* True if ray hit {obj}. */
    double eZ_hit;  /* Value of {eZ} at hit point, if {hit} is true. */

    if (obj->flat)
      { /* Compute the point {h = ray(cZ)} where the ray hits the object's plane: */
        eZ_hit = cZ - pZ; /* {Z} distance from object plane to {p}. */
        /* Compute displacement {sX,sY} from object's center to {h}: */ 
        double sX = pX + dX*eZ_hit - cX; /* {X} position of ray hit rel to object ctr */
        double sY = pY + dY*eZ_hit - cY; /* {Y} position of ray hit rel to object ctr */
        /* Check if ray hits object: */
        double r2_ray = sX*sX + sY*sY;
        hit = (r2_ray <= r2_obj);
        if (hit && verbose) { fprintf(stderr, "    hit disk eZ = %+12.6f\n", eZ_hit); }
      }
    else
      { /* Form the coefficients {A*eZ*eZ + B*eZ + C} of */
        /* {r2_obj = sX*sX + sY*sY + sZ*sZ - r2_obj} as function of {eZ}: */
       
        double LX = pX - cX;
        double LY = pY - cY;
        double LZ = pZ - cZ;
       
        /* sX = pX + dX*eZ - cX  = eZ*dX + LX; */
        /* sY = pY + dY*eZ - cY  = eZ*dY + LY; */
        /* sZ = Z - cZ           = eZ + LZ;    */
        
        double A = dX*dX + dY*dY + 1.0;
        double B = 2*(LX*dX + LY*dY + LZ);
        double C = LX*LX + LY*LY + LZ*LZ - r2_obj;
        
        double Delta = B*B - 4*A*C;
        hit = (Delta > 0);
        if (hit) 
          { eZ_hit = (sqrt(Delta) - B)/(2*A); 
            if (verbose) { fprintf(stderr, "    hit ball eZ = %+12.6f\n", eZ_hit); }
          }
      }
    
    if (hit) { (*hit_P) = (r3_t) {{ pX + dX*eZ_hit, pY + dY*eZ_hit, pZ + eZ_hit }}; }
    return hit;
  }
   
void multifok_scene_ray_trace_backplane(interval_t box[], bool_t back_tilt, r3_t *p, r3_t *d, bool_t debug, r3_t *back_hit_P)
  {     
    bool_t verbose = debug;
    
    /* Let {ray(Z) = p + d*(Z - pZ)} be the point on the ray at height {Z}. */

    /* Grab the relevant coordinates: */
    double pX = p->c[0], pY = p->c[1], pZ = p->c[2];
    double dX = d->c[0], dY = d->c[1], dZ = d->c[2];

    /* Grab the relevant box coordinates: */
    double bXlo = box[0].end[0];
    double bXhi = box[0].end[1];
    double bZlo = box[2].end[0];
    double bZhi = box[2].end[1];
    
    double Z_hit; /* {Z} of point where ray hits backplane. */
    if (back_tilt)
      { demand(dZ = 1.0, "invalid direction vector");
        /* Determine the backplane equation {(Z - bZ) = A*(X - bX)}: */
        /* The plane contains the bottom left and top right edges of the box. */
        /* But we must fudge it a bit because rays may be cast a {skosh} oustde the box. */
        double skosh = 1.0; /* Limit to box overflow due to pixel subsampling. */
        assert(pX >= bXlo - skosh);
        assert(pX <= bXhi + skosh);
        double bX = bXlo - skosh;
        double box_dx = bXhi - bXlo + 2*skosh;
        double bZ = bZlo;
        double box_dz = bZhi - bZlo;
        double A = box_dz/box_dx; /* Backplane tilt, {dZ/dX}. */
        /* From the ray equation, {(X - pX) = dX*(Z - pZ)}. */
        /* Hence at the hit point, {Z - bZ + A*bX = A*X = A*(dX*Z - dX*pZ + pX)}. */
        /* That is, {Z*(1 - A*dX) = bZ + A*(pX - dX*pZ) - A*bX: */
        /* That is, {Z = (bZ + A*(pX - bX - dX*pZ))/(1 - A*dX): */
        double den = 1 - A*dX;
        demand(den > 1.0e-4, "backplane too tilted for aperture");
        Z_hit = (bZ + A*(pX - bX - dX*pZ))/den;
        /* We must clip the backplane to the box {Z} range anyway because of tilted rays. */
        if (Z_hit < bZlo) { Z_hit = bZlo; }
        if (Z_hit > bZhi) { Z_hit = bZhi; }
      }
    else
      { Z_hit = bZlo; }
    
    if (verbose) { fprintf(stderr, "    hit backplane Z = %+12.6f\n", Z_hit); }

    /* For generic {Z}, let {eZ} be {Z - pZ}, then {ray(Z) = p + d*eZ}. */
    double eZ = Z_hit - pZ; /* {Z} distance from hit point to {p}. */
    (*back_hit_P) = (r3_t){{ pX + dX*eZ, pY + dY*eZ, Z_hit }};
  }

#define multifok_scene_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

