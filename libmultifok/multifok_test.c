/* See {multifok_test.h}. */
/* Last edited on 2023-01-25 07:47:54 by stolfi */

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

#include <multifok_test.h>
#include <multifok_focus_op.h>
#include <multifok_scene.h>

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
    parameters {(scene,pattern,p,d,&ray_clr,&ray_hit)}. Its color
    {clr(R)} and hit point {hit(R)} are the restults returned in {ray_clr} and
    {ray_hit}. Its sharpness {shr(R)} is the depth of focus {zDep}
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
    float_image_t *cimg = float_image_new(NC, NX, NY);
    float_image_t *simg = float_image_new(1, NX, NY);
    float_image_t *zimg = float_image_new(1, NX, NY);
    float_image_t *dimg = float_image_new(1, NX, NY);
   
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
          { /* Compute the color and nominal focus of pixel {ix,iy} by ray-tracing, assuming image plane at {Z=zFoc}: */
            frgb_t pix_clr;
            double pix_shr;
            double pix_zav;
            double pix_zvr;
            /* bool_t debug = ((ix == 137) && (iy == NY-111)); */
            bool_t debug = FALSE;
            multifok_test_compute_pixel_color_sharpness_and_zav
              ( scene, pattern, ix, iy, zFoc, zDep, NP, wp, NR, tilt, wr, debug,
                &pix_clr, &pix_shr, &pix_zav, &pix_zvr
              );
            for (int32_t ic = 0; ic < NC; ic++)
              { float_image_set_sample(cimg, ic, ix, iy, pix_clr.c[ic]); }
            float_image_set_sample(simg, 0, ix, iy, (float)pix_shr);         
            float_image_set_sample(zimg, 0, ix, iy, (float)pix_zav);         
            float_image_set_sample(dimg, 0, ix, iy, (float)sqrt(pix_zvr));         
          }
      }
    (*cimg_P) = cimg;
    (*simg_P) = simg;
    (*zimg_P) = zimg;
    (*dimg_P) = dimg;
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
    
    if (verbose) { fprintf(stderr, "=== computing pixel ix = %d iy = %d ===\n", ix, iy); }

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
    
    demand(zDep > 0.0, "invalid {zDep}");
    double zFoc = p->c[2];
    int32_t NC = 3;
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
        frgb_t ray_clr;
        r3_t ray_hit;
        multifok_scene_ray_trace(scene, pattern, p, &d, debug, &ray_clr, &ray_hit);
        double ray_zht = ray_hit.c[2];
        /* The scene surface must be entirely contained in the scene's box: */
        assert(ray_zht >= scene->box[2].end[0]);
        assert(ray_zht <= scene->box[2].end[1]);
        double dz = ray_zht - zFoc;
        double ray_shr = ((dz != 0.0) && isfinite(zDep) ? fabs(zDep/dz) : +INF);/* Sharpness at hit point. */
        
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
    return multifok_test_read_color_image(inPrefix, tag, "c", 0.0, 1.0);
  } 

float_image_t *multifok_test_read_sharpness_image(char *inPrefix, char *tag)
  {
    return multifok_test_read_grayscale_image(inPrefix, tag, "s", 0.0, 1.0);
  }
  
float_image_t *multifok_test_read_zave_image(char *inPrefix, char *tag)
  {
    return multifok_test_read_grayscale_image(inPrefix, tag, "z", 0.0, ZMAX);
  }

float_image_t *multifok_test_read_zdev_image(char *inPrefix, char *tag)
  {
    return multifok_test_read_grayscale_image(inPrefix, tag, "d", 0.0, ZMAX);
  }

void multifok_test_read_term_names_and_weights
  ( char *twfname, 
    int32_t NB, 
    char *belName[],
    int32_t *NT_P, 
    char ***tname_P, 
    double **wt_P, 
    i2_t **tix_P,
    bool_t verbose
  )
  {
    /* Read term names and weights from file: */
    FILE *rd = open_read(twfname, TRUE);
    int32_t NT = 0;
    string_vec_t tname = string_vec_new(50);
    double_vec_t wt = double_vec_new(50);
    while (TRUE)
      { fget_skip_spaces(rd);
        int32_t c = fgetc(rd);
        if (c == EOF) 
          { break; }
        else if (c == '#')
          { do { c = fgetc(rd); } while ((c != '\n') && (c != EOF)); 
            if (c == EOF) { break; }
            continue;
          }
        else if (c == '\n')
          { continue; }
        else
          { ungetc(c, rd); 
            string_vec_expand(&tname,NT);
            double_vec_expand(&wt,NT);
            tname.e[NT] = fget_string(rd);
            wt.e[NT] = fget_double(rd);
            NT++;
            fget_comment_or_eol(rd, '#');
          }
      }
    string_vec_trim(&tname,NT);
    double_vec_trim(&wt,NT);
    fclose(rd);
    
    /* Convert term names to indices: */
    i2_t *tix = (i2_t*)notnull(malloc(NT*sizeof(i2_t)), "no mem");
    multifok_focus_op_term_indices_from_names(NB, belName, NT, tname.e, tix);
    
    if (verbose)
      { for (int32_t kt = 0; kt < NT; kt++)
          { fprintf(stderr, "%3d  %+16.8f", kt, wt.e[kt]);
            char *tnk = tname.e[kt];
            fprintf(stderr, " * %-12s", tnk);
            fprintf(stderr, " = %3d %3d\n", tix[kt].c[0], tix[kt].c[1]);
          }
      }
    
    (*NT_P) = NT;
    (*tname_P) = tname.e;
    (*tix_P) = tix;
    (*wt_P) = wt.e;
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

void multifok_test_write_scene_color_image(float_image_t *cimg, char *outPrefix, char *tag)
  {
    multifok_test_write_color_image(cimg, outPrefix, tag, "c", 0.0, 1.0);
  }

void multifok_test_write_sharpness_image(float_image_t *simg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(simg, outPrefix, tag, "s", 0.0, 1.0);
  }

void multifok_test_write_zavg_image(float_image_t *zimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(zimg, outPrefix, tag, "z", 0.0, ZMAX);
  }

void multifok_test_write_zdev_image(float_image_t *dimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(dimg, outPrefix, tag, "d", 0.0, ZMAX);
  }

void multifok_test_write_score_image(float_image_t *simg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(simg, outPrefix, tag, "f", 0.0, 1.0);
  }

void multifok_test_write_score_error_image(float_image_t *eimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(eimg, outPrefix, tag, "e", -1.0, +1.0);
  }
    
void multifok_test_write_estimated_Z_image(float_image_t *uimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(uimg, outPrefix, tag, "u", 0.0, ZMAX);
  }
 
void multifok_test_write_Z_error_image(float_image_t *vimg, char *outPrefix, char *tag)
  {
    multifok_test_write_grayscale_image(vimg, outPrefix, tag, "v", -ZMAX, +ZMAX);
  }
      
void multifok_test_write_reconstructed_color_image(float_image_t *rimg, char *outPrefix, char *tag)
  {
    multifok_test_write_color_image(rimg, outPrefix, tag, "r", 0.0, 1.0);
  }

void multifok_test_write_color_error_image(float_image_t *qimg, char *outPrefix, char *tag)
  { float vMin = -1.0f;
    float vMax = +1.0f;
    multifok_test_write_color_image(qimg, outPrefix, tag, "q", vMin, vMax);
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

void multifok_test_write_normalized_image(float_image_t *nimg, char *outPrefix, char *tag)
  { float vMin = -1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(nimg, 0, &vMin, &vMax);
    double vR = fmax(fabs(vMin), fabs(vMax));
    vMin = (float)-vR;
    vMax = (float)+vR;

    multifok_test_write_grayscale_image(nimg, outPrefix, tag, "n", vMin, vMax);
  };

void multifok_test_write_sample_weights_image(float_image_t *wimg, char *outPrefix, char *tag)
  { float vMin = 0.0f;  /* To void {NAN} values if {img} is all zeros. */
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(wimg, 0, &vMin, &vMax);
    demand(vMin == 0.0f, "image has negative sample weights");

    multifok_test_write_grayscale_image(wimg, outPrefix, tag, "w", vMin, vMax);
  }
       
void multifok_test_write_basis_elem_image(float_image_t *aimg, char *outPrefix, char *tag)       
  {
    float vMin = -1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(aimg, 0, &vMin, &vMax);
    double vR = fmax(fabs(vMin), fabs(vMax));
    vMin = (float)-vR;
    vMax = (float)+vR;

    multifok_test_write_grayscale_image(aimg, outPrefix, tag, "a", vMin, vMax);
  }

void multifok_test_write_basis_coeff_image(float_image_t *bimg, char *outPrefix, char *tag)
  {
    float vMin = -1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float vMax = +1.0e-38f;  /* To void {NAN} values if {img} is all zeros. */
    float_image_update_sample_range(bimg, 0, &vMin, &vMax); /* Ignores infinities. */
    float vR = (float)fmax(fabs(vMin), fabs(vMax));
    vMin = -vR;
    vMax = +vR;

    multifok_test_write_grayscale_image(bimg, outPrefix, tag, "b", vMin, vMax);
  }

void multifok_test_write_basis_coeff_squared_image(float_image_t *bqimg, char *outPrefix, char *tag)
  {
    double avg, dev;
    float_image_compute_sample_avg_dev(bqimg, 0, &avg, &dev);
    
    float vMin = 0.0; 
    float vMax = (float)fmax(1.0e-38, avg + 2*dev);

    multifok_test_write_grayscale_image(bqimg, outPrefix, tag, "bq", vMin, vMax);
  }

void multifok_test_write_basis_elem_names(char *outPrefix, int32_t NB, char *belName[])  
  { char *fname = NULL;
    asprintf(&fname, "%s-belnames.txt", outPrefix);
    FILE *wr = open_write(fname, TRUE);
    for (int32_t kb = 0; kb < NB; kb++)
      { fprintf(wr, "%s\n", belName[kb]); }
    fclose(wr);
  }
       
void multifok_test_write_term_names_and_weights(char *outPrefix, int32_t NT, char *tname[], double wt[])  
  { char *fname = NULL;
    asprintf(&fname, "%s-terms.txt", outPrefix);
    FILE *wr = open_write(fname, TRUE);
    for (int32_t kt = 0; kt < NT; kt++)
      { fprintf(wr, "%-12s %+12.4f\n", tname[kt], wt[kt]); }
    fclose(wr);
  }

#define multifok_test_C_COPYRIGHT \
    "© 2023 by the State University of Campinas (UNICAMP)"

