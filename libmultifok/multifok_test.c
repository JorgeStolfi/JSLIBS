/* See {multifok_test.h}. */
/* Last edited on 2023-01-10 20:44:46 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <vec.h>
#include <wt_table.h>
#include <affirm.h>
#include <interval.h>
#include <r2.h>
#include <r3.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <jsrandom.h>

#include <float_image.h>

#include <multifok_test.h>

void multifok_test_throw_colors(frgb_t *bg, frgb_t *fg);
  /* Generates two random contrasting colors. */

void multifok_test_ray_trace
  ( multifok_test_scene_t *scene,
    multifok_test_pattern_t *pattern,
    r3_t *p, 
    r3_t *d, 
    frgb_t *clrHit_P, 
    double *zHit_P
  );
  /* Traces one ray that goes through the point {p} with the direction parallel to {d}, assumed to be 
    not horizontal.  Finds the object (disk or backplane) with max {Z} that hits that ray.  Returns its color in
    {*clrHit_P} and the disk's {Z} coordinate in {*zHit_P}.
    
    The color of each object is a mix of its two colors with the
    proportion returned by {pattern(x,y,z)}. */

void multifok_test_compute_pixel_color_and_blur_radius
  ( multifok_test_scene_t *scene,
    multifok_test_pattern_t *pattern,
    int32_t ix, 
    int32_t iy, 
    double rPix, 
    double zFoc, 
    double zDep, 
    int32_t NR,
    frgb_t *clr_P,
    double *rBlur_P
  );
  /* Computes the color {clr} and nominal focal blur radius {rBlur} ay
    the pixel on column {ix} and row {iy}, by ray-tracing the given {scene}
    with {NR} rays. Returns these valus in {*clr_P} and {*rBlur_P}.
    
    The rays will be generated randomly. Each ray {R} will pass through a
    random point {p} in a fuzzy neighborood of the pixel center
    {(ix+0.5,iy+0.5,zFoc)} with radius {rPix}. Its direction {d} will be
    tilted from the vertical by random amounts proportional to {1/zDep}
    (or strictly vertical if {zDep=+INF}).
    
    Each ray is traced with
    {multifok_test_ray_trace(scene,pattern,p,d,&clrHit,&zHit)}. Its
    color {clr{R}) is the returned {clrHit}, and its nominal focal blur
    radius {rBlur(R)} is {rPix} plus (in the {hypot} sense) the absolute
    difference between {zFoc} and the {Z} coordinate {zHit} of the hit
    point, divided by {zDep}. If {zDep} is {+INF} (vertical rays) then
    {rBlur(R)} is just {rPix}.

    The computed pixel color {clr} and its blur radius {rBlur} will be the average of the ray
    colors {clr(R}} and blur radii {rBlur(R)}. */
    
/* IMPLEMENTATIONS */

multifok_test_scene_t *multifok_test_scene_throw
  ( interval_t box[],
    int32_t ND,
    double rMin, 
    double rMax,
    double minSep, 
    bool_t verbose
  )
  {
    demand(box[2].end[0] >= 0, "box extends below {Z=0}");
    
    multifok_test_scene_t *scene = (multifok_test_scene_t *)notnull(malloc(sizeof(multifok_test_scene_t)), "no mem");
    
    if (verbose) { fprintf(stderr, "trying to generate %d disks\n", ND); }
    multifok_test_disk_t *dsk = (multifok_test_disk_t*)notnull(malloc(ND*sizeof(multifok_test_disk_t)), "no mem");
    int32_t ND_out = 0; /* Number of disks actually generated. */
    int32_t NT = (minSep >= 0 ? ND : 5*ND); /* Number of tries. */
    for (int32_t kt = 0; (kt < NT) && (ND_out < ND); kt++)
      { /* Generate a random disk {dk}: */
        multifok_test_disk_t dk = multifok_test_disk_throw(box, rMin, rMax);
        bool_t ok = TRUE; /* False if disk overlaps the previous ones. */
        if (minSep >= 0)
          { /* Check for {XY} overlaps: */
            for (int32_t jd = 0; (jd < ND_out) && ok; jd++)
              { multifok_test_disk_t *dj = &(dsk[jd]);
                /* Check for overlap between {dk} and {dj} plus min separation: */
                double dxy = hypot(dk.ctr.c[0] - dj->ctr.c[0], dk.ctr.c[1] - dj->ctr.c[1]);
                if (dxy < dk.rad + minSep + dj->rad) { ok = FALSE; }
              }
          }
        if (ok) 
          { dsk[ND_out] = dk;
            if (verbose) 
              { r3_gen_print(stderr, &(dk.ctr), "%12.6f", "  ctr = ( ", " ", " )");
                fprintf(stderr, " rad = %12.6f\n", dk.rad);
              }
            ND_out++;
          }
      }
    assert(ND_out <= ND);
    if (ND_out < ND)
      { if (verbose) { fprintf(stderr, "generate only %d disks\n", ND_out); }
        /* Trim array: */
        dsk = (multifok_test_disk_t*)notnull(realloc(dsk, ND_out*sizeof(multifok_test_disk_t)), "no mem");
      }
      
    /* Store data in scene: */
    for (int32_t j = 0; j < 3; j++) { scene->box[j] = box[j]; }
    scene->ND = ND_out;
    scene->dsk = dsk;
    multifok_test_throw_colors(&(scene->bg), &(scene->fg));
    return scene;
  }

multifok_test_disk_t multifok_test_disk_throw(interval_t box[], double rMin, double rMax)
  {
    /* Ensure that the box is large enough in {X} and {Y} to fit at least one disk: */
    for (int32_t j = 0; j < 2; j++)
      { demand(box[j].end[1] - box[j].end[0] > 2*rMax, "{box} too small"); }
    demand(rMin <= rMax, "invarid radius interval");
  
    /* Gnereate a disk fully inside the box. */
    double rad = rMin + drandom()*drandom()*(rMax - rMin);
    r3_t ctr;
    for (int32_t j = 0; j < 3; j++)
      { double xlo = box[j].end[0];
        double xhi = box[j].end[1];
        if (j < 2)
          { /* Shrink the box in {X} and {Y} by {rad}: */
            xlo += rad; xhi -= rad;
          }
        assert(xlo < xhi);
        ctr.c[j] = xlo + drandom()*(xhi - xlo);
      }

    multifok_test_disk_t dsk;
    dsk.ctr = ctr;
    dsk.rad = rad;
    multifok_test_throw_colors(&(dsk.bg), &(dsk.fg));
    return dsk;
  }
  
void multifok_test_throw_colors(frgb_t *bg, frgb_t *fg)
  {
    for (int32_t ic = 0; ic < 3; ic++) 
      { double vi = drandom();
        bg->c[ic] = (float)vi;
        fg->c[ic] = (float)(vi < 0.5 ? 0.5 + vi : vi - 0.5);
      }
  }

void multifok_test_images_make
  ( int32_t NX, 
    int32_t NY, 
    multifok_test_scene_t *scene,
    multifok_test_pattern_t *pattern,
    double rPix,
    double zFoc, 
    double zDep, 
    int32_t NR,
    float_image_t **cimg_P,
    float_image_t **rimg_P
  )
  {
    demand((NX >= 10) && (NY >= 10), "image is too small");
    int32_t NC = 3;
    float_image_t *cimg = float_image_new(NC, NX, NY);
    float_image_t *rimg = float_image_new(1, NX, NY);
    
    for (int32_t iy = 0; iy < NY; iy++)
      { for (int32_t ix = 0; ix < NX; ix++)
          { /* Compute the color and nominal focus of pixel {ix,iy} by ray-tracing, assuming image plane at {Z=zFoc}: */
            frgb_t clr;
            double rBlur;
            multifok_test_compute_pixel_color_and_blur_radius(scene, pattern, ix, iy, rPix, zFoc, zDep, NR, &clr, &rBlur);
            for (int32_t ic = 0; ic < NC; ic++)
              { float_image_set_sample(cimg, ic, ix, iy, clr.c[ic]); }
            float_image_set_sample(rimg, 0, ix, iy, (float)rBlur);         
          }
      }
    (*cimg_P) = cimg;
    (*rimg_P) = rimg;
  }

void multifok_test_compute_pixel_color_and_blur_radius
  ( multifok_test_scene_t *scene,
    multifok_test_pattern_t *pattern,
    int32_t ix, 
    int32_t iy, 
    double rPix, 
    double zFoc, 
    double zDep, 
    int32_t NR,
    frgb_t *clr_P,
    double *rBlur_P
  )
  { 
    int32_t NC = 3;
    double clrTot[NC]; /* Average color of rays. */
    for (int32_t ic = 0; ic < NC; ic++) { clrTot[ic] = 0.0; }
    double rbTot = 0.0; /* Sum of {rb(R)}. */
    for (int32_t ir = 0; ir < NR; ir++)
      { /* Generate point in pixel {íx,iy} at {Z=zFoc}: */
        r3_t p;
        p.c[0] = (double)ix + 0.5 + rPix*dgaussrand();
        p.c[1] = (double)iy + 0.5 + rPix*dgaussrand();
        p.c[2] = zFoc;
        /* Generate a ray with gaussian slope inversely prop to {zDep}: */
        r3_t d;
        d.c[0] = dgaussrand()/zDep;
        d.c[1] = dgaussrand()/zDep;
        d.c[2] = 1.0;
        /* Raytrace the disks: */
        double zHit;
        frgb_t clrHit;
        multifok_test_ray_trace(scene, pattern, &p, &d, &clrHit, &zHit);
        assert(zHit >= 0);
        double dz = (zHit - zFoc)/zDep;
        double rb = hypot(dz, rPix);
        rbTot += rb;
        for (int32_t ic = 0; ic < NC; ic++) { clrTot[ic] += (double)clrHit.c[ic]; }
      }
    /* Compute average ray color: */
    frgb_t clrPix;
    for (int32_t ic = 0; ic < NC; ic++) { clrPix.c[ic] = (float)(clrTot[ic]/NR); }
    /* Compute the average focus radius: */
    double rBlur = rbTot/NR;
    (*clr_P) = clrPix;
    (*rBlur_P) = rBlur;
  }

void multifok_test_ray_trace
  ( multifok_test_scene_t *scene,
    multifok_test_pattern_t *pattern,
    r3_t *p, 
    r3_t *d, 
    frgb_t *clrHit_P, 
    double *zHit_P
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
    
    int32_t kdMax = -1;   /* Index of highest disk hit by ray, or {-1} if none. */
    r3_t hitMax;          /* Point where the ray hits the highest disk. */
    for (int32_t kd = 0; kd < ND; kd++)
      { multifok_test_disk_t *dk = &(scene->dsk[kd]);
        /* Compute the point {hk} where the ray hits the disk's plane: */
        double zhk = dk->ctr.c[2];
        double ezk = zhk - p->c[2]; /* {Z} distance from disk plane to p. */
        r3_t hk = (r3_t){{ p->c[0] + d->c[0]*ezk, p->c[1] + d->c[1]*ezk, zhk }};
        /* Compute displacement {exk,eyk} from disk's center to {hk}: */ 
        double exk = hk.c[0] - dk->ctr.c[0]; /* {X} position of ray hit rel to disk ctr */
        double eyk = hk.c[1] - dk->ctr.c[1]; /* {Y} position of ray hit rel to disk ctr */
        /* Check if ray hits disk: */
        double r2_dsk = dk->rad*dk->rad;
        double r2_ray = exk*exk + eyk*eyk;
        if (r2_ray <= r2_dsk)
          { /* We have a hit: */
            if (zhk > hitMax.c[2])
              { hitMax = hk;
                kdMax = kd;
              }
          }
      }

    /* Get the topmost hit disk, if any: */
    multifok_test_disk_t *dskMax = (kdMax < 0 ? NULL : &(scene->dsk[kdMax]));
    frgb_t *bg;
    frgb_t *fg;
    if (dskMax != NULL)
      { /* Get disk colors: */
        bg = &(dskMax->bg);
        fg = &(dskMax->fg);
      }
    else
      { /* Simulate hit with a backplane at {Z = 0}: */
        double zHit = 0.0;
        double ezb = zHit - p->c[2]; /* {Z} distance from bg plane to {p}. */
        hitMax = (r3_t){{ p->c[0] + d->c[0]*ezb, p->c[1] + d->c[1]*ezb, zHit }};
        /* Get backplane colors: */
        bg = &(scene->bg);
        fg = &(scene->fg);
      }
    double pv = pattern(hitMax.c[0], hitMax.c[1], hitMax.c[2]);
    frgb_t clrHit = frgb_mix(pv, fg, 1-pv, bg);

    (*clrHit_P) = clrHit;
    (*zHit_P) = hitMax.c[2];
  }

#define multifok_test_C_COPYRIGHT \
    "© 2017 by the State University of Campinas (UNICAMP)"

