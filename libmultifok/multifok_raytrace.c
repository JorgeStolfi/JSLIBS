/* See {multifok_raytrace.h}. */
/* Last edited on 2024-10-24 15:34:08 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <bool.h>
#include <r2.h>
#include <r3.h>
#include <i2.h>
#include <affirm.h>
#include <float_image.h>

#include <multifok_frame.h>

#include <multifok_raytrace.h>
  
#define EQUALS "============================================================"
#define SHARPS "############################################################"

void multifok_raytrace_compute_pixel_properties
  ( i2_t *iPix, 
    multifok_raytrace_proc_t *trace_ray,
    multifok_raytrace_img_to_scene_map_t *map_point,
    int32_t NX,
    int32_t NY,
    int32_t NS,
    r2_t uSmp[],
    double wSmp[],
    int32_t NR,
    r2_t tRay[],
    double wRay[],
    r3_t *dRef,
    double zFoc,
    double zDep,
    bool_t debug,
    multifok_raytrace_report_ray_proc_t report_ray,
    int32_t NC,
    float colr_pix[],
    double *vBlr_pix_P,
    double *hAvg_pix_P,
    double *hVar_pix_P
  );
  /* Computes the color {colr(pix)}, blurring indicator {vBlr(pix)},
    height average {hAvg(pix)}, and height variance {hVar(pix)} at
    the pixel {pix} on column {iPix.c[0]} and row {iPix.c[1]}, by ray-tracing 
    some scene or object with multiple rays. Returns these valus in
    {colr_pix[0..NC-1]}, {*vBlr_pix_P}, {*hAvg_pix_P}, and {*hVar_pix_P}.
    
    Specifically, the procedure generates a set of {NS} saub-sampling
    points around the center {pix}, as explained in
    {multifok_raytrace_make_frame}. For each sample point {p}, the
    procedure {multifok_raytrace_compute_image_point_properties} is
    called to compute the average sampling point color {colr(p)}, the
    average height {hAvg(p)}, and its variance {hVar(p)}. These values
    are then combined with the weights {wSmp[0..NS-1]} to obtain the
    pixel color {colr(pix)} and the average and variance
    {hAvg(pix),hVar(pix)} of the height.
    
    The blurring indicator {vBlr} is the weighted mean square radius of the 
    {X,Y} deviations of the ray hit points from the pixel center's {X,Y}.
    This value will always be 1.0 or more on account of the spread of the 
    sampling points {p}. 
    
    If {debug} is true, prints to {stderr} info on the computation of
    the pixel's values. In that case, if {report_ray} is not {NULL},
    calls it for each ray that contributed to those values. */

void multifok_raytrace_compute_sample_point_properties
  ( multifok_raytrace_proc_t *trace_ray,
    r3_t *pRay, 
    int32_t NR,
    r2_t tRay[],
    double wRay[],
    int32_t iniRay,
    int32_t stpRay,
    r3_t *dRef,
    double zFoc,
    double zDep,
    bool_t debug,
    i2_t *iPix,
    r2_t *pSmp,
    double wSmp,
    multifok_raytrace_report_ray_proc_t report_ray,
    int32_t NC,
    float colr_pt[],
    double *vBlr_pt_P,
    double *hAvg_pt_P,
    double *hVar_pt_P
  );
  /* Computes the color {colr(p)}, blurring indicator {vBlr(p)}, height
    average {hAvg(p)}, and height variance {hVar(p)} at the point {p}
    with scene coordinates {pRay}, by ray-tracing the scene or object
    with a set of {KR = NR/stpRay} rays through {pRay}. Returns these
    valus in {colr_pt[0..NC-1]}, {*vBlr_pt_P}, {*hAvg_pt_P}, and
    {*hVar_pt_P}.
    
    As a special case, if {NR} is 1, then {iniRay} and {stpRay} must be
    0, and that single ray ({KR=1}), assumed to be vertical, ({tRay[0] =
    (0,0)}) will be used.
    
    The rays will have weights {wRay[iniRay + jr*stpRay]} directions
    defined by {tRay[iniRay + jr*stpRay]}, for {jr} in {0..KR-1}, as
    described in {multifok_sampling_compute_ray_direction}.
    Currently the ray rotation angle {aRay} will be zero for all 
    
    Each ray {R} with direction {dRay} is traced with the {trace_ray}
    procedure with parameters {pRay,dRay,NC} and should return the
    values of {colr(R)}, {vBlr(R)}, and {hHit(R)}.

    The computed point color {colr(p)}, its blurring indicator {vBlr(p)}, and
    average height {hAvg(p)} will be be the average of the ray colors
    {colr(R}}, blurring indicator {vBlr(R)}, and {hHit(R)} over all rays {R},
    with the respective weights {wRay[0..NR-1]}. The variance {hVar(p)}
    is the weighted variance of {hHit(R)} over those rays. 
    
    The burring indicator {vBlr(p)} of the point will be zero if {zDep} is
    {+INF} or the sampling point {p} lies on the scene's surface.
    In this latter case, {hAvg} will be {zFoc} and {hDev} will be zero.
    
    If {debug} is true, also prints debugging infrormation, and, if
    {report_ray} is not {NULL}, also calls it with the ray data.
    The parameters {iPix,pSmp,wSmp} are used only for this purpose. */

/* IMPLEMENTATIONS */
     
r3_t multifok_raytrace_compute_ray_direction(r2_t *tRay, double aRay, double zDep)
  { demand(zDep >= 1.0e-6, "invalid {zDep}");
    r2_t tRot; r2_rot(tRay, aRay, &tRot);
    double dzRay = (isfinite(zDep) && (r2_norm_sqr(&tRot) != 0) ? -zDep/2 : -1.0);
    r3_t dRay = (r3_t){{ tRot.c[0], tRot.c[1], dzRay }};
    double norm = r3_dir(&dRay, &dRay);
    assert(norm > 1.0e-6);
    return dRay;
  }

multifok_frame_t *multifok_raytrace_make_frame
  ( int32_t NC, 
    int32_t NX, 
    int32_t NY, 
    multifok_raytrace_proc_t *trace_ray, 
    multifok_raytrace_img_to_scene_map_t *map_point,
    r3_t *dRef,
    double zFoc,
    double zDep,
    int32_t NS,
    r2_t uSmp[],
    double wSmp[],
    int32_t NR,
    r2_t tRay[],
    double wRay[],
    bool_t verbose,
    multifok_raytrace_debug_pred_t *debug_pix,
    multifok_raytrace_report_ray_proc_t *report_ray
  )
  {
    float_image_t *sVal = float_image_new(NC, NX, NY);
    float_image_t *shrp = float_image_new(1, NX, NY);
    float_image_t *hAvg = float_image_new(1, NX, NY);
    float_image_t *hDev = float_image_new(1, NX, NY);

    for (int32_t yPix = 0; yPix < NY; yPix++)
      { for (int32_t xPix = 0; xPix < NX; xPix++)
          { /* fprintf(stderr, "(%d,%d)", xPix, yPix); */
            i2_t iPix = (i2_t){{ xPix, yPix }};
            bool_t deb_pix = (debug_pix == NULL ? FALSE : debug_pix(&iPix));
            
            /* Compute the properties of pixel {iPix}: */
            float colr_pix[NC];
            double vBlr_pix;
            double hAvg_pix;
            double hVar_pix;
 
            multifok_raytrace_compute_pixel_properties
              ( &iPix, trace_ray, map_point,
                NX, NY, NS, uSmp, wSmp, NR, tRay, wRay, dRef, zFoc, zDep, 
                deb_pix, report_ray,
                NC, colr_pix, &vBlr_pix, &hAvg_pix, &hVar_pix
              );
            assert(vBlr_pix >= 1.0);
            double shrp_pix = 1/vBlr_pix;
              
            for (int32_t ic = 0; ic < NC; ic++)
              { float_image_set_sample(sVal, ic, xPix, yPix, colr_pix[ic]); }
            float_image_set_sample(hAvg, 0, xPix, yPix, (float)hAvg_pix);         
            float_image_set_sample(hDev, 0, xPix, yPix, (float)sqrt(hVar_pix)); 
            float_image_set_sample(shrp, 0, xPix, yPix, (float)shrp_pix);         
            
          }
      }
      
    multifok_frame_t *frame = multifok_frame_from_images(sVal,shrp,hAvg,hDev, zFoc,zDep);

    return frame;
  }

void multifok_raytrace_compute_pixel_properties
  ( i2_t *iPix, 
    multifok_raytrace_proc_t *trace_ray,
    multifok_raytrace_img_to_scene_map_t *map_point,
    int32_t NX,
    int32_t NY,
    int32_t NS,
    r2_t uSmp[],
    double wSmp[],
    int32_t NR,
    r2_t tRay[],
    double wRay[],
    r3_t *dRef,
    double zFoc,
    double zDep,
    bool_t debug,
    multifok_raytrace_report_ray_proc_t report_ray,
    int32_t NC,
    float colr_pix[],
    double *vBlr_pix_P,
    double *hAvg_pix_P,
    double *hVar_pix_P
  )
  { 
    if (debug) 
      { fprintf(stderr, "  " SHARPS "\n");
        fprintf(stderr, "  computing pixel iPix = ");
        i2_gen_print(stderr, iPix, "%d", "[ ", " ", " ]");
        fprintf(stderr, " with NS = %d NR = %d\n", NS, NR);
      }

    int32_t KR, stpRay;
    if (NR == 1)
      { /* Only one ray to use for all sampling points: */
        KR = 1;
        stpRay = 0;
      }
    else
      { demand((NR % NS) == 0, "ray count {N} must be multiple of {NS}");
        KR = NR/NS;  /* Rays per sampling point. */
        assert(KR >= 1);
        stpRay = NR/KR;  /* Rain index increment for same sampling point. */
      }
      
    /* Scan sampling points, accumulate: */
    double sum_w_colr[NC]; /* Sum of {w(p)*colr(p)}. */
    for (int32_t ic = 0; ic < NC; ic++) { sum_w_colr[ic] = 0.0; }
    double sum_w_vBlr = 0.0; /* Sum of {w(p)*vBlr(R)}. */
    double sum_w_hAvg = 0.0; /* Sum of {w(p)*hAvg(p)}. */
    double sum_w_r2u = 0.0; /* Sum of {w(p)*|p-ctr|^2}. */
    double sum_w = 0.0; /* Sum of {w(p)}. */
    double hAvg_pt[NS]; /* Saved {zAve(p)} for all sampling points. */
    double hVar_pt[NS]; /* Saved {hVar(p)} for all sampling points. */
    for (int32_t ks = 0; ks < NS; ks++)
      { double x_img = iPix->c[0] + 0.5 + uSmp[ks].c[0];
        double y_img = iPix->c[1] + 0.5 + uSmp[ks].c[1];
        double wk = wSmp[ks];
        r2_t pSmp = (r2_t){{ x_img, y_img }};
        r3_t pRay; map_point(&pSmp, &pRay);
        if (debug) 
          { fprintf(stderr, "    sample point %d", ks);
            r2_gen_print(stderr, &pSmp, "%12.6f", "  = image ( ", " ", " )");
            r3_gen_print(stderr, &pRay, "%12.6f", "  = scene ( ", " ", " )");
            fprintf(stderr, " weight = %10.8f\n", wk);
          }
        float colr_pt[NC];
        double vBlr_pt;
        int32_t iniRay = (NR == 1 ? 0 : ks); /* First ray index for next sampling point. */
        multifok_raytrace_compute_sample_point_properties
          ( trace_ray, &pRay, NR, tRay, wRay, iniRay, stpRay,  
            dRef, zFoc, zDep, 
            debug, iPix, &pSmp, wk, report_ray, 
            NC, colr_pt, &vBlr_pt, &(hAvg_pt[ks]), &(hVar_pt[ks])
          );
        if (debug) { fprintf(stderr, "    sample point vBlr = %12.8f hAvg = %12.6f hVar = %16.12f\n", vBlr_pt, hAvg_pt[ks], hVar_pt[ks]); }
        /* Accumulate: */
        for (int32_t ic = 0; ic < NC; ic++) 
          { sum_w_colr[ic] += wk*(double)colr_pt[ic]; }
        assert(vBlr_pt >= 0.0);
        sum_w_vBlr += wk*vBlr_pt;
        sum_w_hAvg += wk*hAvg_pt[ks];
        sum_w_r2u += wk*r2_norm_sqr(&(uSmp[ks]));
        sum_w += wk;
      }
    
    /* Compute average pixel color {colr_pix[0..NC-1]}: */
    for (int32_t ic = 0; ic < NC; ic++) 
      { colr_pix[ic] = (float)(sum_w_colr[ic]/sum_w); }
    
    /* Compute average blurring indicator {vBlr_pix}: */
    double vBlr_pix = sum_w_vBlr/sum_w;
    assert((! isnan(vBlr_pix)) && (vBlr_pix >= 0));
    
    /* Adjust blurring indicator {vBlr_pix} to account for pixel averaging: */
    vBlr_pix = hypot(vBlr_pix, 1.0);
    
    /* Compute pixel height average {hAvg_pix}: */
    double hAvg_pix = sum_w_hAvg/sum_w;
    
    /* Compute pixel height variance {hVar_pix}: */
    double sum_w_dz2 = 0.0;
    for (int32_t ks = 0; ks < NS; ks++)
      { double dz = hAvg_pt[ks] - hAvg_pix; 
        sum_w_dz2 += wSmp[ks]*(dz*dz + hVar_pt[ks]);
      }
    double hVar_pix = sum_w_dz2/sum_w;
    
    if (debug) 
      { double ruAvg = sqrt(sum_w_r2u/sum_w);
        fprintf(stderr, "  pixel [%4d,%4d]", iPix->c[0], iPix->c[1]);
        fprintf(stderr, "  vBlr = %12.8f hAvg = %12.6f hVar = %16.12f", vBlr_pix, hAvg_pix, hVar_pix);
        fprintf(stderr, "  avg |uSmp|^2 = %12.6f\n", ruAvg);
      }

    /* Return pixel data: */
    (*vBlr_pix_P) = vBlr_pix;
    (*hAvg_pix_P) = hAvg_pix;
    (*hVar_pix_P) = hVar_pix;
    
    if (debug) { fprintf(stderr, "  " SHARPS "\n"); }
  }
 
void multifok_raytrace_compute_sample_point_properties
  ( multifok_raytrace_proc_t *trace_ray,
    r3_t *pRay, 
    int32_t NR,
    r2_t tRay[],
    double wRay[],
    int32_t iniRay,
    int32_t stpRay,
    r3_t *dRef,
    double zFoc,
    double zDep,
    bool_t debug,
    i2_t *iPix,
    r2_t *pSmp,
    double wSmp,
    multifok_raytrace_report_ray_proc_t report_ray,
    int32_t NC,
    float colr_pt[],
    double *vBlr_pt_P,
    double *hAvg_pt_P,
    double *hVar_pt_P
  )
  {
    if (debug) 
      { fprintf(stderr, "      " EQUALS "\n");
        r3_gen_print(stderr, pRay, "%12.6f", "      computing properties of scene point ( ", " ", " )");
        fprintf(stderr, " point sampling weight = %.8f\n", wSmp);
      }

    int32_t KR;
    if (NR == 1)
      { demand(stpRay == 0, "invalid {stpRay} for {NR == 1}");
        KR = 1;
      }
    else
      { demand((NR % stpRay) == 0, "invalid {stpRay}");
        KR = NR/stpRay; /* Number of rays for this sampling point. */
      }
    demand(iniRay < NR, "invalid {iniRay}");
    
    double sum_w_colr[NC]; /* Sum of {wRay(R)*colr(R)}. */
    for (int32_t ic = 0; ic < NC; ic++) { sum_w_colr[ic] = 0.0; }
    double sum_w_vBlr = 0.0; /* Sum of {wRay(R)*vBlr(R)}. */
    double sum_w_hHit = 0.0; /* Sum of {wRay(R)*hHit(R)}. */
    double sum_w = 0.0; /* Sum of {wRay(R)}. */
    double hHit_ray[KR]; /* Stores heights of rays used. */
    double w_ray[KR];    /* Stores weights of rays used. */
    for (int32_t jr = 0; jr < KR; jr++)
      { /* INdex of ray in {tRay,wRay}: */
        int32_t ir = iniRay + jr*stpRay;
        assert((ir >= 0) && (ir < NR));
        
        /* Save weight of this ray: */
        w_ray[jr] = wRay[ir];
        
        /* Get the ray's direction {dRay}: */
        double aRay = 0.0;
        r3_t dRay = multifok_raytrace_compute_ray_direction(&(tRay[ir]), aRay, zDep);
        if (debug)
          { fprintf(stderr, "        ray %4d direction = ", ir);
            r3_gen_print(stderr, &(dRay), "%+9.6f", "( ", " ", " )\n");
          }

        /* Trace a ray throug {pRay} with direction {dRay}: */
        float colr_ray[NC];
        r3_t pHit;
        trace_ray(pRay, &dRay, debug, &pHit, NC, colr_ray);
        assert(isfinite(pHit.c[2]));
        assert(isfinite(colr_ray[0]));

        /* Compute and save {hHit_ray} = height of {pHit} from the {zFoc=0} image plane: */
        r3_t u; r3_sub(&pHit, pRay, &u);
        double dHit = r3_dot(&u, dRef); /* Depth of {pHit} away from curr image plane. */
        assert(isfinite(dHit));
        hHit_ray[jr] = zFoc - dHit;
        
        /* Compute {vBlr_ray}: */
        r3_t qRef;  /* Point on ray {pRay,dRef} with same depth as {pHit}. */
        r3_mix(1.0, pRay, dHit, dRef, &qRef);
        double vBlr_ray = r3_dist_sqr(&pHit, &qRef);
        assert(isfinite(vBlr_ray));
        
        if (debug)
          { if (report_ray != NULL)
              { report_ray(iPix, pSmp, wSmp, pRay, &dRay, wRay[ir], &pHit, hHit_ray[jr], vBlr_ray); }
            else
              { fprintf(stderr, "        vBlr = %12.8f hHit = %12.6f w = %16.12f\n", vBlr_ray, hHit_ray[jr], wRay[ir]); }
          }

        /* Accumulate the ray properties onto the point properties: */
        for (int32_t ic = 0; ic < NC; ic++)
          { double colri = fmax(0.0, fmin(1.0, colr_ray[ic]));
            sum_w_colr[ic] += wRay[ir]*(double)colri;
          }
        sum_w_vBlr += w_ray[jr]*vBlr_ray;
        sum_w_hHit += w_ray[jr]*hHit_ray[jr];
        sum_w += w_ray[jr];
      }

    /* Compute average color, blurring indicator, and height at point {p}: */
    assert(sum_w > 0);
    if (debug) {fprintf(stderr, "      sum_w_vBlr = %12.8f sum_w_hHit = %12.6f sum_w = %16.12f\n", sum_w_vBlr, sum_w_hHit, sum_w); }
    for (int32_t ic = 0; ic < NC; ic++) { colr_pt[ic] = (float)(sum_w_colr[ic]/sum_w); }
    double vBlr_pt = sum_w_vBlr/sum_w;
    assert((! isnan(vBlr_pt)) && (vBlr_pt >= 0.0));
    double hAvg_pt = sum_w_hHit/sum_w;
    
    /* Compute Z height variance: */
    double sum_w_dz2 = 0;
    for (int32_t jr = 0; jr < KR; jr++)
      { double dz = hHit_ray[jr] - hAvg_pt;
        sum_w_dz2 += wRay[jr]*dz*dz;
      }
    double hVar_pt = sum_w_dz2/sum_w;
    
    /* Return results: */
    (*vBlr_pt_P) = vBlr_pt;
    (*hAvg_pt_P) = hAvg_pt;
    (*hVar_pt_P) = hVar_pt;
    
    if (debug) { fprintf(stderr, "      " EQUALS "\n"); }
  }

void multifok_raytrace_show_ray_data
  ( FILE *wr,
    i2_t *iPix,
    double pixSize,
    r2_t *pSmp,
    double wSmp,
    r3_t *pRay,
    r3_t *dRay, 
    double wRay,
    r3_t *pHit, 
    double hHit,
    double vBlr
  )
  {
    i2_gen_print(wr, iPix, "%4d", "pixel [", ",", "]");
    fprintf(wr, " size %.6f × %.6f", pixSize, pixSize);
    fprintf(wr, " subsampling point ( %.3f %.3f )", pSmp->c[0], pSmp->c[1]);
    fprintf(wr, " weight = %.8f", wSmp);
    r3_gen_print(wr, pRay, "%10.4f", " = scene ( ", " ", " )");
    r3_gen_print(wr, dRay, "%+9.6f", " dir = ( ", " ", " )");
    fprintf(wr, " weight = %.8f", wRay);
    r3_gen_print(wr, pHit, "%10.4f", " hit = ( ", " ", " )");
    fprintf(wr, " vBlr = %10.4f hHit = %10.4f", vBlr, hHit);
    fprintf(wr, "\n");
    fflush(wr);
  }
    
void multifok_raytrace_write_ray_data
  ( FILE *wr,
    i2_t *iPix,
    double pixSize,
    r2_t *pSmp,
    double wSmp,
    r3_t *pRay,
    r3_t *dRay, 
    double wRay,
    r3_t *pHit, 
    double hHit,
    double vBlr
  )
  {
    fprintf(wr, "%4d %4d", iPix->c[0], iPix->c[1]);
    fprintf(wr, " %9.6f", pixSize);
    fprintf(wr, "  %11.6f %11.6f", pSmp->c[0], pSmp->c[1]);
    fprintf(wr, "  %10.8f", wSmp);
    r3_gen_print(wr, pRay, "%12.6f", "  ", " ", " ");
    r3_gen_print(wr, dRay, "%+9.6f", "  ", " ", " ");
    fprintf(wr, "   %10.8f", wRay);
    r3_gen_print(wr, pHit, "%12.6f", "  ", " ", " ");
    fprintf(wr, "  %12.8f %12.6f", vBlr, hHit);
    fprintf(wr, "\n");
    fflush(wr);
  }
