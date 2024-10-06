/* See {multifok_test.h}. */
/* Last edited on 2024-10-01 14:37:06 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <wt_table.h>
#include <wt_table_hann.h>
#include <affirm.h>
#include <interval.h>
#include <fget.h>
#include <i2.h>
#include <r2.h>
#include <r3.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <jsqroots.h>
#include <rn.h>

#include <float_image.h>
#include <float_image_paint.h>
#include <float_image_read_pnm.h>
#include <float_image_write_pnm.h>

#include <multifok_window.h>
#include <multifok_window.h>
#include <multifok_basis.h>
#include <multifok_score.h>
#include <multifok_term.h>
#include <multifok_scene.h>

#include <multifok_test.h>

#define ZMAX multifok_scene_ZMAX
  /* Max scene {Z} coordinate. */

#define EQUALS "============================================================"
#define SHARPS "############################################################"

void multifok_test_choose_ray_tilts_and_weights
  ( int32_t NR_min, 
    double zDep, 
    int32_t *NR_P,
    r2_t **tilt_P, 
    double **wr_P
  );
  /* Chooses a number {NR} of rays, their deviations from the vertical {tilt[0..NR-1},
    and their weights {wr[0..NR-1]}. Allocates the arrays and returns the choices
    in {*NR_P}, {*tilt_P}, {*wr_P}.
    
    The number {NR} of rays will be at least {NR_min}, but will be rounded up as needed to
    provide good coverage.
    
    The tilt of each ray {kr} will be such that for each unit of displacement in {Z}
    the ray will deviate from the vertical by {tilt[kr].c[0]} in {X} and {tilt[kr].c[0]} in {Y}.
    The maximum tilt will be {2/zDep}.  The first ray will always be vertical.
    
    As a special case, if {zDep} is {+INF}, then {NR_min} must be 1, and there
    will be only the first vertical ray. Otherwise {NR} will be at least 7. */

void multifok_test_compute_focus_plane_point_color_sharpness_and_zave
  ( multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    r3_t *p, 
    double zDep, 
    int32_t NR,
    r2_t tilt[],
    double wr[],
    bool_t debug,
    frgb_t *pt_colr_P,
    double *pt_shrp_P,
    double *pt_zave_P,
    double *pt_zvar_P
  );
  /* Computes the color {colr(p)}, nominal sharpness {shrp(p)}, scene {Z} average
    {zave(p)}, and scene {Z} variance {zvar(p)} at the point {p}, assumed to be a point
    on the focus plane, by ray-tracing the given {scene} with about
    {NR} rays through {p}. Returns these valus in {*pt_colr_P},
    {*pt_shrp_P}, {*pt_zave_P}, and {*pt_zvar_P}.
    
    The rays will go through the point {p} and will have 
    directions {(tilt[ir].c[0, tilt[ir].c[1], 1.0)}, for {ir}
    in {0..NR-1}.
    
    Each ray {R} is traced with {multifok_scene_ray_trace} with
    parameters {(scene,pattern,p,d,&ray_htob,&ray_htpt)}. Its color
    {colr(R)} and hit point {hit(R)} are the restults returned in {ray_colr} and
    {ray_htpt}. Its sharpness {shrp(R)} is the depth of focus {zDep}
    divided by the absolute difference between {zFoc} and {zval(R)=hit(R).Z}. In
    particular, if {zDep} is {+INF} (vertical rays), or {zval(R)} is
    equal to {zFoc} then {shrp(R)} is {+INF}.

    The computed point color {colr(p)}, its sharpness {shrp(p)}, and average
    height {zave(p)} will be be the average of the ray colors {colr(R}},
    sharpnesses {shrp(R)}, and {zval(R)} over all rays {R}, with the 
    respective weights {wr[0..NR-1]}. The variance 
    {zvar(p)} is the variance of {zval(R)} over those rays.*/

void multifok_test_compute_pixel_color_sharpness_and_zave
  ( int32_t ix,
    int32_t iy, 
    multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    int32_t NX,
    int32_t NY,
    double zFoc, 
    double zDep, 
    int32_t NP,
    double wp[],
    int32_t NR,
    r2_t tilt[],
    double wr[],
    bool_t debug,
    frgb_t *pix_colr_P,
    double *pix_shrp_P,
    double *pix_zave_P,
    double *pix_zvar_P
  );
  /* Computes the color {colr(pix)}, nominal sharpness {shrp(pix)}, scene {Z}
    average {zave(pix)}, and scene {Z} variance at the pixel {pix} on column {ix}
    and row {iy}, by ray-tracing the given {scene} with multiple rays.
    Returns these valus in {*pix_colr_P} {*pix_shrp_P}, {*pix_zave_P}, and {*pix_zvar_P}.
    
    The image-to-scene scale  implicitly chosen so that the image domain fits snugly in
    the scene domain
    
    Specifically, the procedure generates an array of {NP} by {NP}
    sample point in and around the pixel {pix}, as explained in
    {multifok_test_images_make}. For eeach sample ppoint {p}, procedure
    {multifok_test_compute_focus_plane_point_color_sharpness_and_zave} is
    used ith parameters {(scene, pattern, &p, zDep, NR, tilt, wr, &pt_colr,
    &pt_shrp, &pt_zave, &pt_zvar)} to compute the average point color
    {colr(p)=pt_colr}, the average sharpness {shrp(p)=pt_shrp}, the
    average {Z} coordinate {zave(p)=pt_zave}, and its variance {zvar(p)=pt_zvar}.

    The computed pixel color {colr(pix)}, sharpness {shrp(pix)}, and scene {Z} height {zave(pix)} will be 
    averages of {colr(p)}, {shrp(p)}, and {zave(p)} over those points, with 2D Hann weights.
    The varainces are appropriately combined into {zvar(p)} */
    
/* IMPLEMENTATIONS */

void multifok_test_images_make
  ( int32_t NX, 
    int32_t NY, 
    multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    double zFoc, 
    double zDep, 
    int32_t NP,
    int32_t NR_min,
    float_image_t **csimg_P,
    float_image_t **shimg_P,
    float_image_t **azimg_P,
    float_image_t **dzimg_P
  )
  {
    demand((NX >= 10) && (NY >= 10), "image is too small");
    int32_t NC = 3;
    
    float_image_t *csimg = float_image_new(NC, NX, NY);
    float_image_t *shimg = float_image_new(1, NX, NY);
    float_image_t *azimg = float_image_new(1, NX, NY);
    float_image_t *dzimg = float_image_new(1, NX, NY);
   
    /* Pixel sampling weights: */
    double wp[NP];
    wt_table_hann_fill(NP, 0.0, wp, NULL);
    wt_table_normalize_sum(NP, wp);
    
    /* Ray tilts and weights: */
    int32_t NR; /* Actual number of aperture rays. */
    r2_t *tilt;
    double *wr;
    multifok_test_choose_ray_tilts_and_weights(NR_min, zDep, &NR, &tilt, &wr);
    
    for (int32_t iy = 0; iy < NY; iy++)
      { for (int32_t ix = 0; ix < NX; ix++)
          { /* fprintf(stderr, "(%d,%d)", ix, iy); */
            bool_t debug = FALSE;
            
            if (debug)
              { /* Restrict debugging to one pixel */ 
                if (scene->NO == 0)
                  { debug = ((ix == 10) && (iy == 10)); }
                else
                  { /* Lets debug a pixel right over object 0: */
                    r3_t *ctr_obj = &(scene->objs[0].ctr);
                    r2_t p_img_obj = multifok_test_scene_to_pixel(ctr_obj->c[0], ctr_obj->c[1], scene->dom, NX, NY);
                    debug = ((fabs(ix - p_img_obj.c[0]) < 0.55) && (fabs(iy - p_img_obj.c[1]) < 0.55)); 
                  }
              }
            
            /* Compute the color and nominal focus of pixel {ix,iy} by ray-tracing, assuming image plane at {Z=zFoc}: */
            frgb_t pix_colr;
            double pix_shrp;
            double pix_zave;
            double pix_zvar;
            
            multifok_test_compute_pixel_color_sharpness_and_zave
              ( ix, iy, scene, pattern, NX, NY, zFoc, zDep, NP, wp, NR, tilt, wr, debug,
                &pix_colr, &pix_shrp, &pix_zave, &pix_zvar
              );
            for (int32_t ic = 0; ic < NC; ic++)
              { float_image_set_sample(csimg, ic, ix, iy, pix_colr.c[ic]); }
            float_image_set_sample(shimg, 0, ix, iy, (float)pix_shrp);         
            float_image_set_sample(azimg, 0, ix, iy, (float)pix_zave);         
            float_image_set_sample(dzimg, 0, ix, iy, (float)sqrt(pix_zvar));         
          }
      }
    (*csimg_P) = csimg;
    (*shimg_P) = shimg;
    (*azimg_P) = azimg;
    (*dzimg_P) = dzimg;
  }
    
void multifok_test_estimate_pixel_zave_zdev
  ( int32_t ix,
    int32_t iy,
    multifok_scene_t *scene,
    int32_t NX,
    int32_t NY,
    double *zave_P,
    double *zdev_P
  )
  { 
    r3_t d = (r3_t){{ 0.0, 0.0, 1.0 }}; /* Ray direction. */
    int32_t NR = 9; /* Number of rays to trace. */
    double zval[NR]; /* {Z} surface coordinates at the subpixel rays. */
    int32_t kz = 0;
    for (int32_t ex = -1; ex <= +1; ex++)
      { for (int32_t ey = -1; ey <= +1; ey++)
          { /* Compute a point {p}directly above center of pixel {(ix,iy)} perturbed by {0.5*(ex,ey)}: */
            r3_t p2 = multifok_test_pixel_to_scene(ix + 0.5*(ex + 1), iy + 0.5*(ey + 1), NX, NY, scene->dom);
            r3_t p = (r3_t){{ p2.c[0], p2.c[1], ZMAX + 1.0 }};
            /* Ray-trace the scene from {p} straight down: */
            multifok_scene_object_t *ray_htob; /* Object hit, or {NULL} if floor. */
            r3_t ray_htpt; /* Hit point. */
            multifok_scene_ray_trace(scene, &p, &d, FALSE, &ray_htob, &ray_htpt);
            assert(kz < NR);
            zval[kz] = ray_htpt.c[2]; /* {Z} of hit points. */
            kz++; 
          }
      }
    assert(kz == NR);
    double sum_zval = 0.0;
    for (int32_t kz = 0; kz < NR; kz++) { sum_zval += zval[kz]; }
    double zave = sum_zval/NR;
    double sum_dz2 = 1.0e-200;
    for (int32_t kz = 0; kz < NR; kz++) { double dz = zval[kz] - zave; sum_dz2 += dz*dz; }
    double zdev = sqrt(sum_dz2/NR);
    (*zave_P) = zave;
    (*zdev_P) = zdev;
  }

void multifok_test_compute_pixel_color_sharpness_and_zave
  ( int32_t ix,
    int32_t iy, 
    multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    int32_t NX,
    int32_t NY,
    double zFoc, 
    double zDep, 
    int32_t NP,
    double wp[],
    int32_t NR,
    r2_t tilt[],
    double wr[],
    bool_t debug,
    frgb_t *pix_colr_P,
    double *pix_shrp_P,
    double *pix_zave_P,
    double *pix_zvar_P
  )
  { 
    bool_t verbose = debug;
    
    if (verbose) 
      { fprintf(stderr, SHARPS "\n");
        fprintf(stderr, "computing pixel ix = %d iy = %d\n", ix, iy);
      }
    
    int32_t NC = 3;
    /* Scan sampling points, store sample point averages and variances, and sample weights: */
    double sum_wp_colr[NC]; /* Sum of {w(p)*colr(p)}. */
    for (int32_t ic = 0; ic < NC; ic++) { sum_wp_colr[ic] = 0.0; }
    double sum_wp_shrp = 0.0; /* Sum of {w(p)*shrp(R)}. */
    double sum_wp_zave = 0.0; /* Sum of {w(p)*zave(p)}. */
    double sum_wp_dp2 = 0.0; /* Sum of {w(p)*dp(p)^2}, {dp(p)=|(p.x-ctr.x,p.y-ctr.y)|}. */
    double sum_wp = 0.0; /* Sum of {w(p)}. */
    int32_t NS = NP*NP; /* Number of sampling points. */
    double pt_zaves[NS];
    double pt_zvars[NS];
    double pt_wps[NS];
    int32_t ks = 0; /* Sample index in {0..NS-1}. */
    for (int32_t jx = 0; jx < NP; jx++)
      { double dpx = 2*(jx + 0.5)/NP - 1.0;
        double x_img = ix + 0.5 + dpx;
        for (int32_t jy = 0; jy < NP; jy++)
          { double dpy = 2*(jy + 0.5)/NP - 1.0;
            double y_img = iy + 0.5 + dpy;
            /* Map point of image to scene cordinates: */
            r3_t p_scene = multifok_test_pixel_to_scene(x_img, y_img, NX, NY, scene->dom);
            /* Set the {Z} coordinate to {zFoc} instead of zero: */
            p_scene.c[2] = zFoc;
            double pt_wp = wp[jx]*wp[jy];
            frgb_t pt_colr;
            double pt_shrp;
            double pt_zave;
            double pt_zvar;
            if (verbose) 
              { fprintf(stderr, "  sample point jx = %d jy = %d", jx, jy);
                r3_gen_print(stderr, &p_scene, "%12.6f", "  = ( ", " ", " )\n");
              }
            multifok_test_compute_focus_plane_point_color_sharpness_and_zave
              ( scene, pattern, &p_scene, zDep, NR, tilt, wr, debug,
                &pt_colr, &pt_shrp, &pt_zave, &pt_zvar
              );
            if (verbose) { fprintf(stderr, "  shrp = %12.8f zave = %12.6f zvar = %16.12f\n", pt_shrp, pt_zave, pt_zvar); }
            /* Accumulate: */
            for (int32_t ic = 0; ic < NC; ic++) 
              { sum_wp_colr[ic] += pt_wp*(double)pt_colr.c[ic]; }
            assert(pt_shrp > 0.0);
            sum_wp_shrp += pt_wp*pt_shrp;
            sum_wp_zave += pt_wp*pt_zave;
            sum_wp_dp2 += pt_wp*(dpx*dpx + dpy*dpy);
            sum_wp += pt_wp;
            /* Save averages and variances, and sample weights: */
        
            pt_zaves[ks] = pt_zave;
            pt_zvars[ks] = pt_zvar;
            pt_wps[ks] = pt_wp;
            ks++;
          }
      }
    assert(ks == NS);
    
    /* Compute average pixel color, radius, sharpness, height: */
    frgb_t pix_colr;
    for (int32_t ic = 0; ic < NC; ic++) { pix_colr.c[ic] = (float)(sum_wp_colr[ic]/sum_wp); }
    double pix_srad = sqrt(sum_wp_dp2/sum_wp);
    double pix_shrp = sum_wp_shrp/sum_wp;
    assert(pix_shrp > 0);
     /* Adjust sharpness to account for pixel averaging: */
    double pix_krad = (pix_shrp == +INF ? 0.0 : 1/pix_shrp); /* Avg sampling kernel radius. */
    pix_krad = hypot(pix_srad, pix_krad);
    pix_shrp = (pix_krad == 0 ? +INF : 1.0/pix_krad);
    double pix_zave = sum_wp_zave/sum_wp;
    if (verbose) { fprintf(stderr, "  srad = %12.8f  krad = %12.8f shrp = %12.8f\n", pix_srad, pix_krad, pix_shrp); }
    
    /* Compute pixel {Z} variance: */
    double sum_wp_dz2 = 0.0;
    for (int32_t ks = 0; ks < NS; ks++)
      { double dz = pt_zaves[ks] - pix_zave; 
        sum_wp_dz2 += pt_wps[ks]*(dz*dz + pt_zvars[ks]);
      }
    double pix_zvar = sum_wp_dz2/sum_wp;
    
    if (verbose) 
      { fprintf(stderr, "pixel average srad = %12.8f shrp = %12.8f zave = %12.6f zvar = %16.12f\n", pix_srad, pix_shrp, pix_zave, pix_zvar); }

    /* Return pixel data: */
    (*pix_colr_P) = pix_colr;
    (*pix_shrp_P) = pix_shrp;
    (*pix_zave_P) = pix_zave;
    (*pix_zvar_P) = pix_zvar;
    
    if (verbose) { fprintf(stderr, SHARPS "\n"); }
  }
 
void multifok_test_compute_focus_plane_point_color_sharpness_and_zave
  ( multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    r3_t *p, 
    double zDep, 
    int32_t NR,
    r2_t tilt[],
    double wr[],
    bool_t debug,
    frgb_t *pt_colr_P,
    double *pt_shrp_P,
    double *pt_zave_P,
    double *pt_zvar_P
  )
  {
    bool_t verbose = debug;
    
    if (verbose) 
      { fprintf(stderr, "  " EQUALS "\n");
        r3_gen_print(stderr, p, "%12.6f", "  computing point ( ", " ", " )\n");
      }
    
    demand(zDep > 0.0, "invalid {zDep}");
    double zFoc = p->c[2];
    int32_t NC = 3;
    interval_t zRange = scene->dom[2];
    double sum_wr_colr[NC]; /* Sum of {wr(R)*colr(R)}. */
    for (int32_t ic = 0; ic < NC; ic++) { sum_wr_colr[ic] = 0.0; }
    double sum_wr_shrp = 0.0; /* Sum of {wr(R)*shrp(R)}. */
    double sum_wr_zval = 0.0; /* Sum of {wr(R)*zval(R)}. */
    double sum_wr = 0.0; /* Sum of {wr(R)}. */
    double ray_zvals[NR]; /* Stores ray heights for {zdv} computation. */
    for (int32_t ir = 0; ir < NR; ir++)
      { /* Generate a ray with predefined deviation from vertical: */
        r3_t d;
        d.c[0] = tilt[ir].c[0];
        d.c[1] = tilt[ir].c[1];
        d.c[2] = 1.0;
        /* Raytrace the objects: */
        multifok_scene_object_t *ray_htob; /* Object hit, or {NULL} if floor. */
        r3_t ray_htpt; /* Hit point. */
        multifok_scene_ray_trace(scene, p, &d, debug, &ray_htob, &ray_htpt);
        double ray_zval = ray_htpt.c[2];
        /* The scene surface must be entirely contained in the scene's box: */
        if ((ray_zval < zRange.end[0]) || (ray_zval > zRange.end[1]))
          { fprintf(stderr, "    ray hit at Z = %16.12f, outside range {%16.12f _ %16.12f]\n", ray_zval, zRange.end[0], zRange.end[1]);
            assert(FALSE);
          }
        double dz = ray_zval - zFoc;
        double ray_shrp = ((zDep == +INF) || (dz == 0) ? +INF : fabs(zDep/dz));/* Sharpness at hit point. */
        
        /* Compute the hit point's color: */
        frgb_t ray_colr = multifok_scene_compute_hit_color(ray_htob, &(ray_htpt), pattern);

        /* Accumulate point properties: */
        double wri = wr[ir];
        if (verbose) {fprintf(stderr, "    ray %3d shrp = %12.8f zval = %12.6f wri = %16.12f\n", ir, ray_shrp, ray_zval, wri); }
        for (int32_t ic = 0; ic < NC; ic++) { sum_wr_colr[ic] += wri*(double)ray_colr.c[ic]; }
        sum_wr_shrp += wri*ray_shrp;
        sum_wr_zval += wri*ray_zval;
        sum_wr += wri;
        ray_zvals[ir] = ray_zval;
      }

    /* Compute average color, sharpnes, and height at point {p}: */
    assert(sum_wr > 0);
    if (verbose) {fprintf(stderr, "  sum_wr_shrp = %12.8f sum_wr_zval = %12.6f sum_wr = %16.12f\n", sum_wr_shrp, sum_wr_zval, sum_wr); }
    frgb_t pt_colr;
    for (int32_t ic = 0; ic < NC; ic++) { pt_colr.c[ic] = (float)(sum_wr_colr[ic]/sum_wr); }
    double pt_shrp = sum_wr_shrp/sum_wr;
    assert(pt_shrp > 0.0);
    double pt_zave = sum_wr_zval/sum_wr;
    
    /* Compute Z height variance: */
    double sum_wr_dz2 = 0;
    for (int32_t ir = 0; ir < NR; ir++)
      { double dz = ray_zvals[ir] - pt_zave;
        sum_wr_dz2 += wr[ir]*dz*dz;
      }
    double pt_zvar = sum_wr_dz2/sum_wr;
    
    /* Return results: */
    (*pt_colr_P) = pt_colr;
    (*pt_shrp_P) = pt_shrp;
    (*pt_zave_P) = pt_zave;
    (*pt_zvar_P) = pt_zvar;
    
    if (verbose) { fprintf(stderr, "  " EQUALS "\n"); }
  }
  
r3_t multifok_test_pixel_to_scene
  ( double x_img,
    double y_img,
    int32_t NX,
    int32_t NY,
    interval_t scene_dom[]
  )
  {
    /* Scale from image pixel {XY} to scene {XY}: */
    double WX_scene = 2*interval_rad(&(scene_dom[0]));
    double WY_scene = 2*interval_rad(&(scene_dom[1]));
    double scale_img = fmax(((double)NX)/WX_scene, ((double)NY)/WY_scene);

    double x_scene_ctr = interval_mid(&(scene_dom[0]));
    double y_scene_ctr = interval_mid(&(scene_dom[1]));
    double x_img_ctr = 0.5*NX;
    double y_img_ctr = 0.5*NY;
    double x_scene = x_scene_ctr + (x_img - x_img_ctr)/scale_img;
    double y_scene = y_scene_ctr + (y_img - y_img_ctr)/scale_img;
    return (r3_t){{ x_scene, y_scene, 0 }};
  }
  
r2_t multifok_test_scene_to_pixel
  ( double x_scene,
    double y_scene,
    interval_t scene_dom[],
    int32_t NX,
    int32_t NY
  )
  {
    /* Scale from image pixel {XY} to scene {XY}: */
    double WX_scene = 2*interval_rad(&(scene_dom[0]));
    double WY_scene = 2*interval_rad(&(scene_dom[1]));
    double scale_img = fmax(((double)NX)/WX_scene, ((double)NY)/WY_scene);

    double x_scene_ctr = interval_mid(&(scene_dom[0]));
    double y_scene_ctr = interval_mid(&(scene_dom[1]));
    double x_img_ctr = 0.5*NX;
    double y_img_ctr = 0.5*NY;
    double x_img = y_img_ctr + (x_scene - x_scene_ctr)*scale_img;
    double y_img = x_img_ctr + (y_scene - y_scene_ctr)*scale_img;
    return (r2_t){{ x_img, y_img }};    
  }
  
void multifok_test_choose_ray_tilts_and_weights
  ( int32_t NR_min, 
    double zDep, 
    int32_t *NR_P,
    r2_t **tilt_P, 
    double **wr_P
  )
  {
    bool_t verbose = TRUE;
    
    demand(NR_min >= 1, "invalid {NR_min}");
    demand(zDep > 0.0, "invalid {zDep}");
    
    /* Determine the number of ray layers (besides the central ray): */
    int32_t NL;
    if (isfinite(zDep)) 
      { /* Determine the number {NL} of layers beyond the central ray: */
        NL = (int32_t)ceil(0.5*(1 + sqrt(1 + 4.0*(NR_min-1)/3.0)));
        assert(NL >= 1);
      }
    else
      { demand(NR_min == 1, "{NR_min} should be 1 if {zDep} is infinite");
        NL = 0;
      }
    /* Compute the total number of rays {NR}: */
    int32_t NR = 1 + 6*(NL-1)*NL/2;
    assert(NR >= NR_min);
    /* Allocate the arrays: */
    r2_t *tilt = (r2_t*)notnull(malloc(NR*sizeof(r2_t)), "no mem"); 
    double *wr = rn_alloc(NR); 
    /* Fill the central ray: */
    tilt[0] = (r2_t){{ 0.0, 0.0 }};
    wr[0] = 1.0;
    if (NL > 0)
      { /* Generate the rays by layers: */
        double rMax = 2.0/zDep; /* Limit of ray tilt from vertical. */
        double rStep = rMax/(NL + 0.5); /* Radius increment between layers. */
        double rad = rStep; /* Radius of next layer. */
        int32_t kr = 1; /* Number of rays generated. */
        for (int32_t kl = 0; kl < NL; kl++)
          { /* Compute the weight for all rays in this layer: */
            double cr = 0.5*(1 + cos(M_PI*rad/rMax)); /* Apodizing weight. */
            cr = 1 - (1 - cr)*(1 - cr);
            double wrl = 2.0*cr/3.0;
            /* Generate {6*kl} rays on layer {kl} */
            double dang = 2*M_PI/(6*kl); /* Angular step in each layer. */
            double ang = 0.0;
            for (int32_t ks = 0; ks < 6; ks++)
              { for (int32_t ir = 0; ir < kl; ir++) 
                  { tilt[kr].c[0] = rad*cos(ang);
                    tilt[kr].c[1] = rad*sin(ang);
                    wr[kr] = wrl;
                    kr++;
                    ang += dang;
                  }
              }
            rad += rStep;
          }
        assert(kr == NR);
      }
      
    if (verbose)
      { fprintf(stderr, "using %d aperture rays (central plus %d layers)\n", NR, NL); 
        for (int32_t kr = 0; kr < NR; kr++)
          { fprintf(stderr, "  tilt[%3d] = ( %+12.6f %+12.6f )", kr, tilt[kr].c[0], tilt[kr].c[1]); 
            fprintf(stderr, "  wr = %12.8f\n", wr[kr]);
          }
      }
    /* Return the results: */
    (*NR_P) = NR;
    (*tilt_P) = tilt;
    (*wr_P) = wr;
  }

/* Terms data I/O */

/* IMAGE I/O */

float_image_t *multifok_test_read_color_image(char *inPrefix, char *tag, char *code, float vMin, float vMax)
  {
    char *fname = NULL;
    asprintf(&fname, "%s%s-%s.ppm", inPrefix, tag, code);
    bool_t isMask = FALSE;
    double gama = multifok_test_image_gamma;
    double bias = multifok_test_image_bias;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = FALSE;
    float_image_t *img = float_image_read_pnm_named(fname, isMask, gama, bias, yup, warn, verbose);
    
    int32_t NC = (int32_t)img->sz[0];
    demand(NC == 3, "color image should have 3 channels");
    
    if ((vMin != 0.0) || (vMax != 1.0))
      { for (int32_t ic = 0; ic < NC; ic++)
          { float_image_rescale_samples(img, ic, 0.0,1.0, vMin,vMax); }
      }
    
    free(fname);
    return img;
  }

float_image_t *multifok_test_read_grayscale_image(char *inPrefix, char *tag, char *code, float vMin, float vMax)
  {
    char *fname = NULL;
    asprintf(&fname, "%s%s-%s.pgm", inPrefix, tag, code);
    bool_t isMask = FALSE;
    double gama = multifok_test_image_gamma;
    double bias = multifok_test_image_bias;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = FALSE;
    float_image_t *img = float_image_read_pnm_named(fname, isMask, gama, bias, yup, warn, verbose);

    int32_t NC = (int32_t)img->sz[0];
    demand(NC == 1, "grayscale image should have 1 channels");

    if ((vMin != 0.0) || (vMax != 1.0))
      { float_image_rescale_samples(img, 0, 0.0,1.0, vMin,vMax); }

    free(fname);
    return img;
  }
      
float_image_t *multifok_test_read_scene_color_image(char *inPrefix, char *tag)
  {
    return multifok_test_read_color_image(inPrefix, tag, "cs", 0.0, 1.0);
  } 

float_image_t *multifok_test_read_sharpness_image(char *inPrefix, char *tag)
  {
    return multifok_test_read_grayscale_image(inPrefix, tag, "sh", 0.0, 1.0);
  }
  
float_image_t *multifok_test_read_zave_image(char *inPrefix, char *tag)
  {
    return multifok_test_read_grayscale_image(inPrefix, tag, "az", 0.0, ZMAX);
  }

float_image_t *multifok_test_read_zdev_image(char *inPrefix, char *tag)
  {
    return multifok_test_read_grayscale_image(inPrefix, tag, "dz", 0.0, ZMAX);
  }


void multifok_test_write_color_image
  ( float_image_t *img, 
    char *outPrefix, 
    char *tag, 
    char *code, 
    float vMin, 
    float vMax
  )
  {
    int32_t NC = (int32_t)img->sz[0];
    demand(NC == 3, "color image should have 3 channels");
    
    float_image_t *cpy = float_image_copy(img);

    if ((vMin != 0.0) || (vMax != 1.0))
      { for (int32_t ic = 0; ic < NC; ic++)
          { float_image_rescale_samples(cpy, ic, vMin, vMax, 0.0, 1.0); }
      }

    char *fname = NULL;
    asprintf(&fname, "%s%s-%s.ppm", outPrefix, tag, code);
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    double gamma = multifok_test_image_gamma;
    double bias = multifok_test_image_bias;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = TRUE;
    float_image_write_pnm_named(fname, cpy, isMask, gamma, bias, yup, warn, verbose);
    float_image_free(cpy);
    free(fname);
  }
  
void multifok_test_write_grayscale_image(float_image_t *img, char *outPrefix, char *tag, char *code, float vMin, float vMax)
  {
    int32_t NC = (int32_t)img->sz[0];
    demand(NC == 1, "grayscale image should have 1 channel");
    
    float_image_t *cpy = float_image_copy(img);
    
    if ((vMin != 0.0) || (vMax != 1.0))
      { float_image_rescale_samples(cpy, 0, vMin, vMax, 0.0, 1.0); }

    char *fname = NULL;
    asprintf(&fname, "%s%s-%s.pgm", outPrefix, tag, code);
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    double gamma = multifok_test_image_gamma;
    double bias = multifok_test_image_bias;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = TRUE;
    float_image_write_pnm_named(fname, cpy, isMask, gamma, bias, yup, warn, verbose);
    float_image_free(cpy);
    free(fname);
  }

void multifok_test_draw_crosses(float_image_t *img, int32_t ch, i2_vec_t pix, float val)
  { 
    int32_t NC = (int32_t)img->sz[0];
    
    /* Draw crosses over the image: */
    double rad = 6.0;
    double hwd = 0.5;
    bool_t empty = TRUE;
    bool_t diagonal = TRUE;
    for (int32_t kq = 0; kq < pix.ne; kq++)
      { double xctr = pix.e[kq].c[0] + 0.5;
        double yctr = pix.e[kq].c[1] + 0.5;
        for (int32_t kc = 0; kc < NC; kc++)
          { float_image_paint_cross(img, kc, xctr, yctr, rad, empty, 2.0*hwd, diagonal, 0.0f, 3);
            float valk = (kc == ch ? val : 0.0f);
            float_image_paint_cross(img, kc, xctr, yctr, rad, empty, hwd, diagonal, valk, 3);
          }
      }
  }

/* SPECIALIZED IMAGE OUTPUT */

void multifok_test_write_scene_color_image(float_image_t *csimg, char *outPrefix, char *tag)
  {
    multifok_test_write_color_image(csimg, outPrefix, tag, "cs", 0.0, 1.0);
  }

void multifok_test_write_sharpness_image(float_image_t *shimg, char *outPrefix, char *tag)
  {
    float vMin = 0.0f; 
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(shimg, 0, &vMin, &vMax);
    demand(vMin == 0.0f, "negative samples in sharpness image");

    multifok_test_write_grayscale_image(shimg, outPrefix, tag, "sh", 0.0, vMax);
  }

void multifok_test_write_zave_image(float_image_t *azimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(azimg, outPrefix, tag, "az", 0.0, ZMAX);
  }

void multifok_test_write_zdev_image(float_image_t *dzimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(dzimg, outPrefix, tag, "dz", 0.0, ZMAX);
  }

void multifok_test_write_score_image(float_image_t *scimg, char *outPrefix, char *tag)
  {
    float vMin = 0.0f; 
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(scimg, 0, &vMin, &vMax);
    demand(vMin == 0.0f, "negative samples in score image");

    multifok_test_write_grayscale_image(scimg, outPrefix, tag, "sc", 0.0, 1.0);
  }

void multifok_test_write_score_error_image(float_image_t *esimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(esimg, outPrefix, tag, "es", -1.0, +1.0);
  }
    
void multifok_test_write_estimated_Z_image(float_image_t *czimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(czimg, outPrefix, tag, "cz", 0.0, ZMAX);
  }
 
void multifok_test_write_Z_error_image(float_image_t *ezimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(ezimg, outPrefix, tag, "ez", -ZMAX, +ZMAX);
  }
      
void multifok_test_write_reconstructed_color_image(float_image_t *crimg, char *outPrefix, char *tag)
  {
    multifok_test_write_color_image(crimg, outPrefix, tag, "cr", 0.0, 1.0);
  }

void multifok_test_write_color_error_image(float_image_t *erimg, char *outPrefix, char *tag)
  { 
    multifok_test_write_color_image(erimg, outPrefix, tag, "er", -1.0f, +1.0f);
  }

void multifok_test_write_pix_sel_image(float_image_t *bgimg, i2_vec_t pix, char *outPrefix, char *tag)
  { demand(bgimg->sz[0] == 3, "background must be a color image");
    float_image_t *psimg = float_image_copy(bgimg);
    multifok_test_draw_crosses(psimg, 0, pix, 1.0f);
    multifok_test_write_color_image(psimg, outPrefix, tag, "ps", 0.0f, 1.0f);
    float_image_free(psimg);
  }

void multifok_test_write_window_average_image(float_image_t *avimg, char *outPrefix, char *tag)
  { 
    multifok_test_write_grayscale_image(avimg, outPrefix, tag, "av", 0.0f, 1.0f);
  }

void multifok_test_write_window_gradient_image(float_image_t *gvimg, char *outPrefix, char *tag)
  { float vMin = 0.0f; 
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(gvimg, 0, &vMin, &vMax);
    demand(vMin == 0.0f, "negative samples in deviation image");
    
    multifok_test_write_grayscale_image(gvimg, outPrefix, tag, "gv", vMin, vMax);
  };

void multifok_test_write_window_deviation_image(float_image_t *dvimg, char *outPrefix, char *tag)
  { float vMin = 0.0f; 
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(dvimg, 0, &vMin, &vMax);
    demand(vMin == 0.0f, "negative samples in deviation image");

    multifok_test_write_grayscale_image(dvimg, outPrefix, tag, "dv", vMin, vMax);
  };

void multifok_test_write_normalized_image(float_image_t *nrimg, char *outPrefix, char *tag)
  { float vMin = -1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(nrimg, 0, &vMin, &vMax);
    double vR = fmax(fabs(vMin), fabs(vMax));
    vMin = (float)-vR;
    vMax = (float)+vR;

    multifok_test_write_grayscale_image(nrimg, outPrefix, tag, "nr", vMin, vMax);
  };

void multifok_test_write_sample_weights_image(float_image_t *wsimg, char *outPrefix, char *tag)
  { float vMin = 0.0f; 
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(wsimg, 0, &vMin, &vMax);
    demand(vMin == 0.0f, "image has negative sample weights");

    multifok_test_write_grayscale_image(wsimg, outPrefix, tag, "ws", vMin, vMax);
  }
       
void multifok_test_write_basis_elem_image(float_image_t *beimg, char *outPrefix, char *tag)       
  {
    float vMin = -1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(beimg, 0, &vMin, &vMax);
    double vR = fmax(fabs(vMin), fabs(vMax));
    vMin = (float)-vR;
    vMax = (float)+vR;

    multifok_test_write_grayscale_image(beimg, outPrefix, tag, "be", vMin, vMax);
  }

void multifok_test_write_basis_coeff_image(float_image_t *bcimg, char *outPrefix, char *tag)
  {
    float vMin = -1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(bcimg, 0, &vMin, &vMax); /* Ignores infinities. */
    float vR = (float)fmax(fabs(vMin), fabs(vMax));
    vMin = -vR;
    vMax = +vR;

    multifok_test_write_grayscale_image(bcimg, outPrefix, tag, "bc", vMin, vMax);
  }

void multifok_test_write_basis_coeff_squared_image(float_image_t *bqimg, char *outPrefix, char *tag)
  {
    double avg, dev;
    float_image_compute_sample_avg_dev(bqimg, 0, &avg, &dev);
    
    float vMin = 0.0; 
    float vMax = (float)fmax(1.0e-38, avg + 2.5*dev);

    multifok_test_write_grayscale_image(bqimg, outPrefix, tag, "bq", vMin, vMax);
  }

void multifok_test_write_quadratic_term_image(float_image_t *tmimg, char *outPrefix, char *tag)
  {
    float vMin = -1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(tmimg, 0, &vMin, &vMax); /* Ignores infinities. */
    float vR = (float)fmax(fabs(vMin), fabs(vMax));
    vMin = -vR;
    vMax = +vR;
    multifok_test_write_grayscale_image(tmimg, outPrefix, tag, "tm", vMin, vMax);
  }

void multifok_test_write_pixel_mask_image(float_image_t *mkimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(mkimg, outPrefix, tag, "mk", 0.0, 1.0);
  }

FILE* multifok_test_open_text_file(char *outPrefix, char *tag)
  { 
    char *fname = NULL;
    asprintf(&fname, "%s%s.txt", outPrefix, tag);
    FILE *wr = open_write(fname, TRUE);
    free(fname);
    return wr;
  }

void multifok_test_write_basis_elem_names(char *outPrefix, int32_t NB, char *belName[])  
  { char *fname = NULL;
    asprintf(&fname, "%s-bnames.txt", outPrefix);
    FILE *wr = open_write(fname, TRUE);
    multifok_basis_write_elem_names(wr, NB, belName);
    fclose(wr);
  }
 
void multifok_test_write_term_names
  ( char *outPrefix, 
    int32_t NT, 
    char *termName[]
  )
  { char *fname = NULL;
    asprintf(&fname, "%s-tnames.txt", outPrefix);
    FILE *wr = open_write(fname, TRUE);
    multifok_term_write_names(wr, NT, termName);
    fclose(wr);
  }
  
void multifok_test_read_term_weights_and_names
  ( char *fname, 
    int32_t *NT_P,
    double **wt_P,
    char ***termName_P,
    bool_t verbose
  )
  {
    FILE *rd = open_read(fname, TRUE);
    multifok_term_read_weights_and_names(rd, NT_P, wt_P, termName_P, verbose);
    fclose(rd);
  }

void multifok_test_write_term_weights_and_names
  ( char *outPrefix, 
    int32_t NT, 
    double wt[], 
    char *termName[]
  )
  { char *fname = NULL;
    asprintf(&fname, "%s-twts.txt", outPrefix);
    FILE *wr = open_write(fname, TRUE);
    multifok_term_write_weights_and_names(wr, NT, wt, termName);
    fclose(wr);
  }
   
void multifok_test_read_term_weights_names_get_indices
  ( char *fname, 
    int32_t NB, 
    char *belName[],
    int32_t *NT_P,
    double **wt_P, 
    char ***termName_P, 
    int32_t *NP_P, 
    multifok_term_prod_t **prix_P,
    bool_t verbose
  )
  {
    FILE *rd = open_read(fname, TRUE);
    multifok_score_read_term_weights_names_get_indices
      ( rd, 
        NB, belName, 
        NT_P, wt_P, termName_P, 
        NP_P, prix_P, 
        verbose
      );
    fclose(rd);
  }
   
void multifok_test_read_term_index_table
  ( char *fname, 
    int32_t NB,
    char *belName[],
    int32_t *NP_P,
    multifok_term_prod_t **prix_P,
    int32_t *NT_P,
    char ***termName_P,
    bool_t verbose
  )
  {
    FILE *rd = open_read(fname, TRUE);
    multifok_term_read_index_table(rd, NB, belName, NP_P, prix_P, NT_P, termName_P, verbose);
    fclose(rd);
  }

void multifok_test_write_term_index_table
  ( char *outPrefix, 
    int32_t NP, 
    multifok_term_prod_t prix[],
    bool_t verbose
  )
  {
    char *fname = NULL;
    asprintf(&fname, "%s-prix.txt", outPrefix);
    FILE *wr = open_write(fname, TRUE);
    multifok_term_write_index_table(wr, NP, prix, verbose);
    fclose(wr);
    free(fname);
  }

#define multifok_test_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

