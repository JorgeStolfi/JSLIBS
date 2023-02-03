/* See {multifok_test.h}. */
/* Last edited on 2023-01-31 21:45:44 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <wt_table.h>
#include <affirm.h>
#include <interval.h>
#include <fget.h>
#include <r2.h>
#include <r3.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <jsqroots.h>

#include <float_image.h>
#include <float_image_read_pnm.h>
#include <float_image_write_pnm.h>

#include <multifok_window.h>
#include <multifok_window.h>
#include <multifok_basis.h>
#include <multifok_score.h>
#include <multifok_term.h>
#include <multifok_scene.h>

#include <multifok_test.h>

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

void multifok_test_compute_focus_plane_point_color_sharpness_and_zav
  ( multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    r3_t *p, 
    double zDep, 
    int32_t NR,
    r2_t tilt[],
    double wr[],
    bool_t debug,
    frgb_t *pt_clr_P,
    double *pt_shr_P,
    double *pt_zav_P,
    double *pt_zvr_P
  );
  /* Computes the color {clr(p)}, nominal sharpness {shr(p)}, scene {Z} average
    {zav(p)}, and scene {Z} variance {zvr(p)} at the point {p}, assumed to be a point
    on the focus plane, by ray-tracing the given {scene} with about
    {NR} rays through {p}. Returns these valus in {*pt_clr_P},
    {*pt_shr_P}, {*pt_zav_P}, and {*pt_zvr_P}.
    
    The rays will go through the point {p} and will have 
    directions {(tilt[ir].c[0, tilt[ir].c[1], 1.0)}, for {ir}
    in {0..NR-1}.
    
    Each ray {R} is traced with {multifok_scene_ray_trace} with
    parameters {(scene,pattern,p,d,&ray_hob,&ray_hpt)}. Its color
    {clr(R)} and hit point {hit(R)} are the restults returned in {ray_clr} and
    {ray_hpt}. Its sharpness {shr(R)} is the depth of focus {zDep}
    divided by the absolute difference between {zFoc} and {zht(R)=hit(R).Z}. In
    particular, if {zDep} is {+INF} (vertical rays), or {zht(R)} is
    equal to {zFoc} then {shr(R)} is {+INF}.

    The computed point color {clr(p)}, its sharpness {shr(p)}, and average
    height {zav(p)} will be be the average of the ray colors {clr(R}},
    sharpnesses {shr(R)}, and {zht(R)} over all rays {R}, with the 
    respective weights {wr[0..NR-1]}. The variance 
    {zvr(p)} is the variance of {zht(R)} over those rays.*/

void multifok_test_compute_pixel_color_sharpness_and_zav
  ( multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    int32_t ix, 
    int32_t iy, 
    double zFoc, 
    double zDep, 
    int32_t NP,
    double wp[],
    int32_t NR,
    r2_t tilt[],
    double wr[],
    bool_t debug,
    frgb_t *pix_clr_P,
    double *pix_shr_P,
    double *pix_zav_P,
    double *pix_zvr_P
  );
  /* Computes the color {clr(pix)}, nominal sharpness {shr(pix)}, scene {Z}
    average {zav(pix)}, and scene {Z} variance at the pixel {pix} on column {ix}
    and row {iy}, by ray-tracing the given {scene} with multiple rays.
    Returns these valus in {*pix_clr_P} {*pix_shr_P}, {*pix_zav_P}, and {*pix_zvr_P}.
    
    Specifically, the procedure generates an array of {NP} by {NP}
    sample point in and around the pixel {pix}, as explained in
    {multifok_test_images_make}. For eeach sample ppoint {p}, procedure
    {multifok_test_compute_focus_plane_point_color_sharpness_and_zav} is
    used ith parameters {(scene, pattern, &p, zDep, NR, tilt, wr, &pt_clr,
    &pt_shr, &pt_zav, &pt_zvr)} to compute the average point color
    {clr(p)=pt_clr}, the average sharpness {shr(p)=pt_shr}, the
    average {Z} coordinate {zav(p)=pt_zav}, and its variance {zvr(p)=pt_zvr}.

    The computed pixel color {clr(pix)}, sharpness {shr(pix)}, and scene {Z} height {zav(pix)} will be 
    averages of {clr(p)}, {shr(p)}, and {zav(p)} over those points, with 2D Hann weights.
    The varainces are appropriately combined into {zvr(p)} */
    
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
    float_image_t **cimg_P,
    float_image_t **simg_P,
    float_image_t **zimg_P,
    float_image_t **dimg_P
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
    bool_t normalize_wp = FALSE;
    wt_table_fill_hann(NP, wp, normalize_wp);
    
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
                    r3_t *ctr = &(scene->objs[0].ctr);
                    debug = ((ix == (int32_t)(ctr->c[0] + 0.5)) && (iy == (int32_t)(ctr->c[1] + 0.5))); 
                  }
              }
            
            /* Compute the color and nominal focus of pixel {ix,iy} by ray-tracing, assuming image plane at {Z=zFoc}: */
            frgb_t pix_clr;
            double pix_shr;
            double pix_zav;
            double pix_zvr;
            
            multifok_test_compute_pixel_color_sharpness_and_zav
              ( scene, pattern, ix, iy, zFoc, zDep, NP, wp, NR, tilt, wr, debug,
                &pix_clr, &pix_shr, &pix_zav, &pix_zvr
              );
            for (int32_t ic = 0; ic < NC; ic++)
              { float_image_set_sample(csimg, ic, ix, iy, pix_clr.c[ic]); }
            float_image_set_sample(shimg, 0, ix, iy, (float)pix_shr);         
            float_image_set_sample(azimg, 0, ix, iy, (float)pix_zav);         
            float_image_set_sample(dzimg, 0, ix, iy, (float)sqrt(pix_zvr));         
          }
      }
    (*cimg_P) = csimg;
    (*simg_P) = shimg;
    (*zimg_P) = azimg;
    (*dimg_P) = dzimg;
  }

void multifok_test_compute_pixel_color_sharpness_and_zav
  ( multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    int32_t ix, 
    int32_t iy, 
    double zFoc, 
    double zDep, 
    int32_t NP,
    double wp[],
    int32_t NR,
    r2_t tilt[],
    double wr[],
    bool_t debug,
    frgb_t *pix_clr_P,
    double *pix_shr_P,
    double *pix_zav_P,
    double *pix_zvr_P
  )
  { 
    bool_t verbose = debug;
    
    if (verbose) 
      { fprintf(stderr, SHARPS "\n");
        fprintf(stderr, "computing pixel ix = %d iy = %d\n", ix, iy);
      }

    int32_t NC = 3;
    /* Scan sampling points, store sample point averages and variances, and sample weights: */
    double sum_wp_clr[NC]; /* Sum of {w(p)*clr(p)}. */
    for (int32_t ic = 0; ic < NC; ic++) { sum_wp_clr[ic] = 0.0; }
    double sum_wp_shr = 0.0; /* Sum of {w(p)*rPix/rBlur(R)}. */
    double sum_wp_zav = 0.0; /* Sum of {w(p)*zav(p)}. */
    double sum_wp_dp2 = 0.0; /* Sum of {w(p)*dp(p)^2}, {dp(p)=|(p.x-ctr.x,p.y-ctr.y)|}. */
    double sum_wp = 0.0; /* Sum of {w(p)}. */
    r3_t p;
    p.c[2] = zFoc;
    int32_t NS = NP*NP; /* Number of sampling points. */
    double pt_zavs[NS];
    double pt_zvrs[NS];
    double pt_wps[NS];
    int32_t ks = 0; /* Sample index in {0..NS-1}. */
    for (int32_t jx = 0; jx < NP; jx++)
      { double dpx = 2*(jx + 0.5)/NP - 1.0;
        p.c[0] = ix + 0.5 + dpx;
        for (int32_t jy = 0; jy < NP; jy++)
          { double dpy = 2*(jy + 0.5)/NP - 1.0;
            p.c[1] = iy + 0.5 + dpy;
            double pt_wp = wp[jx]*wp[jy];
            frgb_t pt_clr;
            double pt_shr;
            double pt_zav;
            double pt_zvr;
            if (verbose) 
              { fprintf(stderr, "  sample point jx = %d jy = %d", jx, jy);
                r3_gen_print(stderr, &p, "%12.6f", "  = ( ", " ", " )\n");
              }
            multifok_test_compute_focus_plane_point_color_sharpness_and_zav
              ( scene, pattern, &p, zDep, NR, tilt, wr, debug,
                &pt_clr, &pt_shr, &pt_zav, &pt_zvr
              );
            if (verbose) { fprintf(stderr, "  shr = %12.8f zav = %12.6f zvr = %16.12f\n", pt_shr, pt_zav, pt_zvr); }
            /* Accumulate: */
            for (int32_t ic = 0; ic < NC; ic++) 
              { sum_wp_clr[ic] += pt_wp*(double)pt_clr.c[ic]; }
            sum_wp_shr += pt_wp*pt_shr;
            sum_wp_zav += pt_wp*pt_zav;
            sum_wp_dp2 += pt_wp*(dpx*dpx + dpy*dpy);
            sum_wp += pt_wp;
            /* Save averages and variances, and sample weights: */
        
            pt_zavs[ks] = pt_zav;
            pt_zvrs[ks] = pt_zvr;
            pt_wps[ks] = pt_wp;
            ks++;
          }
      }
    assert(ks == NS);
    
    /* Compute average pixel color, radius, sharpness, height: */
    frgb_t pix_clr;
    for (int32_t ic = 0; ic < NC; ic++) { pix_clr.c[ic] = (float)(sum_wp_clr[ic]/sum_wp); }
    double pix_rad = sqrt(sum_wp_dp2/sum_wp);
    double pix_shr = sum_wp_shr/sum_wp;
    double pix_zav = sum_wp_zav/sum_wp;

    /* Adjust sharpness to account for pixel averaging: */
    if (verbose) { fprintf(stderr, "adjusting pixel average shr: before = %12.8f", pix_shr); }
    double pix_rbl = (isfinite(pix_shr) ? 1/pix_shr : 0.0); /* Avg blur radius. */
    pix_rbl = hypot(pix_rad, pix_rbl);
    pix_shr = pix_rad/pix_rbl;
    if (verbose) { fprintf(stderr, "  rad = %12.8f  after = %12.8f\n", pix_rad, pix_shr); }
    
    /* Compute pixel {Z} variance: */
    double sum_wp_dz2 = 0.0;
    for (int32_t ks = 0; ks < NS; ks++)
      { double dz = pt_zavs[ks] - pix_zav; 
        sum_wp_dz2 += pt_wps[ks]*(dz*dz + pt_zvrs[ks]);
      }
    double pix_zvr = sum_wp_dz2/sum_wp;
    
    if (verbose) 
      { fprintf(stderr, "pixel average rad = %12.8f shr = %12.8f zav = %12.6f zvr = %16.12f\n", pix_rad, pix_shr, pix_zav, pix_zvr); }

    /* Return pixel data: */
    (*pix_clr_P) = pix_clr;
    (*pix_shr_P) = pix_shr;
    (*pix_zav_P) = pix_zav;
    (*pix_zvr_P) = pix_zvr;
    
    if (verbose) { fprintf(stderr, SHARPS "\n"); }
  }
 
void multifok_test_compute_focus_plane_point_color_sharpness_and_zav
  ( multifok_scene_t *scene,
    multifok_scene_pattern_t *pattern,
    r3_t *p, 
    double zDep, 
    int32_t NR,
    r2_t tilt[],
    double wr[],
    bool_t debug,
    frgb_t *pt_clr_P,
    double *pt_shr_P,
    double *pt_zav_P,
    double *pt_zvr_P
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
    double sum_wr_clr[NC]; /* Sum of {wr(R)*clr(R)}. */
    for (int32_t ic = 0; ic < NC; ic++) { sum_wr_clr[ic] = 0.0; }
    double sum_wr_shr = 0.0; /* Sum of {wr(R)*shr(R)}. */
    double sum_wr_zht = 0.0; /* Sum of {wr(R)*zht(R)}. */
    double sum_wr = 0.0; /* Sum of {wr(R)}. */
    double ray_zhts[NR]; /* Stores ray heights for {zdv} computation. */
    for (int32_t ir = 0; ir < NR; ir++)
      { /* Generate a ray with predefined deviation from vertical: */
        r3_t d;
        d.c[0] = tilt[ir].c[0];
        d.c[1] = tilt[ir].c[1];
        d.c[2] = 1.0;
        /* Raytrace the objects: */
        multifok_scene_object_t *ray_hob; /* Object hit, or {NULL} if floor. */
        r3_t ray_hpt; /* Hit point. */
        multifok_scene_ray_trace(scene, p, &d, debug, &ray_hob, &ray_hpt);
        double ray_zht = ray_hpt.c[2];
        /* The scene surface must be entirely contained in the scene's box: */
        if ((ray_zht < zRange.end[0]) || (ray_zht > zRange.end[1]))
          { fprintf(stderr, "    ray hit at Z = %16.12f, outside range {%16.12f _ %16.12f]\n", ray_zht, zRange.end[0], zRange.end[1]);
            assert(FALSE);
          }
        double dz = ray_zht - zFoc;
        double ray_shr = ((dz != 0.0) && isfinite(zDep) ? fabs(zDep/dz) : +INF);/* Sharpness at hit point. */
        
        /* Compute the hit point's color: */
        frgb_t ray_clr = multifok_scene_compute_hit_color(ray_hob, &(ray_hpt), pattern);

        /* Accumulate point properties: */
        double wri = wr[ir];
        if (verbose) {fprintf(stderr, "    ray %3d shr = %12.8f zht = %12.6f wri = %16.12f\n", ir, ray_shr, ray_zht, wri); }
        for (int32_t ic = 0; ic < NC; ic++) { sum_wr_clr[ic] += wri*(double)ray_clr.c[ic]; }
        sum_wr_shr += wri*ray_shr;
        sum_wr_zht += wri*ray_zht;
        sum_wr += wri;
        ray_zhts[ir] = ray_zht;
      }

    /* Compute average color, sharpnes, and height at point {p}: */
    assert(sum_wr > 0);
    if (verbose) {fprintf(stderr, "  sum_wr_shr = %12.8f sum_wr_zht = %12.6f sum_wr = %16.12f\n", sum_wr_shr, sum_wr_zht, sum_wr); }
    frgb_t pt_clr;
    for (int32_t ic = 0; ic < NC; ic++) { pt_clr.c[ic] = (float)(sum_wr_clr[ic]/sum_wr); }
    double pt_shr = sum_wr_shr/sum_wr;
    double pt_zav = sum_wr_zht/sum_wr;
    
    /* Compute Z height variance: */
    double sum_wr_dz2 = 0;
    for (int32_t ir = 0; ir < NR; ir++)
      { double dz = ray_zhts[ir] - pt_zav;
        sum_wr_dz2 += wr[ir]*dz*dz;
      }
    double pt_zvr = sum_wr_dz2/sum_wr;
    
    /* Return results: */
    (*pt_clr_P) = pt_clr;
    (*pt_shr_P) = pt_shr;
    (*pt_zav_P) = pt_zav;
    (*pt_zvr_P) = pt_zvr;
    
    if (verbose) { fprintf(stderr, "  " EQUALS "\n"); }
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
    double *wr = (double*)notnull(malloc(NR*sizeof(double)), "no mem"); 
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

#define ZMAX multifok_scene_ZMAX

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


void multifok_test_write_color_image(float_image_t *img, char *outPrefix, char *tag, char *code, float vMin, float vMax)
  {
    int32_t NC = (int32_t)img->sz[0];
    demand(NC == 3, "color image should have 3 channels");
    
    if ((vMin != 0.0) || (vMax != 1.0))
      { for (int32_t ic = 0; ic < NC; ic++)
          { float_image_rescale_samples(img, ic, vMin, vMax, 0.0, 1.0); }
      }

    char *fname = NULL;
    asprintf(&fname, "%s%s-%s.ppm", outPrefix, tag, code);
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    double gamma = multifok_test_image_gamma;
    double bias = multifok_test_image_bias;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = TRUE;
    float_image_write_pnm_named(fname, img, isMask, gamma, bias, yup, warn, verbose);
    free(fname);
  }
  
void multifok_test_write_grayscale_image(float_image_t *img, char *outPrefix, char *tag, char *code, float vMin, float vMax)
  {
    int32_t NC = (int32_t)img->sz[0];
    demand(NC == 1, "grayscale image should have 1 channel");
    
    if ((vMin != 0.0) || (vMax != 1.0))
      { float_image_rescale_samples(img, 0, vMin, vMax, 0.0, 1.0); }

    char *fname = NULL;
    asprintf(&fname, "%s%s-%s.pgm", outPrefix, tag, code);
    bool_t isMask = FALSE; /* Assume uniform distr. of pixel values in encoding/decoding. */
    double gamma = multifok_test_image_gamma;
    double bias = multifok_test_image_bias;
    bool_t yup = TRUE;
    bool_t warn = TRUE;
    bool_t verbose = TRUE;
    float_image_write_pnm_named(fname, img, isMask, gamma, bias, yup, warn, verbose);
    free(fname);
  }

void multifok_test_write_scene_color_image(float_image_t *csimg, char *outPrefix, char *tag)
  {
    multifok_test_write_color_image(csimg, outPrefix, tag, "cs", 0.0, 1.0);
  }

void multifok_test_write_sharpness_image(float_image_t *shimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(shimg, outPrefix, tag, "sh", 0.0, 1.0);
  }

void multifok_test_write_zavg_image(float_image_t *azimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(azimg, outPrefix, tag, "az", 0.0, ZMAX);
  }

void multifok_test_write_zdev_image(float_image_t *dzimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(dzimg, outPrefix, tag, "dz", 0.0, ZMAX);
  }

void multifok_test_write_score_image(float_image_t *shimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(shimg, outPrefix, tag, "sc", 0.0, 1.0);
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
  { float vMin = -1.0f;
    float vMax = +1.0f;
    multifok_test_write_color_image(erimg, outPrefix, tag, "er", vMin, vMax);
  }

void multifok_test_write_window_average_image(float_image_t *avimg, char *outPrefix, char *tag)
  { 
    multifok_test_write_grayscale_image(avimg, outPrefix, tag, "av", 0.0f, 1.0f);
  }

void multifok_test_write_window_deviation_image(float_image_t *dvimg, char *outPrefix, char *tag)
  { float vMin = 0.0f;  /* To void {NAN} values if {img} is all zeros. */
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
  { float vMin = 0.0f;  /* To void {NAN} values if {img} is all zeros. */
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

void multifok_test_read_term_indices_names_and_weights
  ( char *twfname, 
    int32_t NB, 
    char *belName[], 
    int32_t *NP_P, 
    multifok_term_prod_t **prix_P,
    int32_t *NT_P,
    double **wt_P, 
    char ***termName_P,
    bool_t verbose
  )
  {
    /* Read term names and weights from file: */
    FILE *rd = open_read(twfname, TRUE);
    multifok_score_read_term_names_and_weights(rd, NB, belName, NP_P, prix_P, NT_P, wt_P, termName_P, verbose);
    fclose(rd);
  }

void multifok_test_write_basis_elem_names(char *outPrefix, int32_t NB, char *belName[])  
  { char *fname = NULL;
    asprintf(&fname, "%s-bnames.txt", outPrefix);
    FILE *wr = open_write(fname, TRUE);
    multifok_basis_write_elem_names(wr, NB, belName);
    fclose(wr);
  }

void multifok_test_write_term_names(char *outPrefix, int32_t NT, char *termName[])  
  { char *fname = NULL;
    asprintf(&fname, "%s-tnames.txt", outPrefix);
    FILE *wr = open_write(fname, TRUE);
    multifok_term_write_names(wr, NT, termName);
    fclose(wr);
  }
       
#define multifok_test_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

