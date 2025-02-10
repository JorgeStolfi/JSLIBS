/* See {multifok_raytrace.h}. */
/* Last edited on 2025-02-08 11:10:49 by stolfi */

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
#include <jsrandom.h>
#include <affirm.h>
#include <float_image.h>

#include <multifok_frame.h>
#include <multifok_sampling.h>

#include <multifok_raytrace.h>
  
#define PINGOS "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - "
#define DASHES "----------------------------------------------------------------------"
#define EQUALS "======================================================================"
#define SHARPS "######################################################################"

void multifok_raytrace_compute_pixel_properties
  ( i2_t *iPix, 
    multifok_raytrace_proc_t *trace_ray,
    multifok_raytrace_img_to_scene_map_t *map_point,
    int32_t xLo,
    int32_t xHi,
    int32_t yLo,
    int32_t yHi,
    multifok_sampling_t *samp,
    r3_t *dRef,
    double zFoc,
    double zDep,
    multifok_raytrace_debug_pred_t *debug_pixel,
    multifok_raytrace_report_ray_proc_t report_ray,
    multifok_raytrace_report_pixel_proc_t report_pixel,
    r3_t *sNrm_pix_P,
    int32_t NC,
    float sVal_pix[],
    double *vBlr_pix_P,
    double *hAvg_pix_P,
    double *hVar_pix_P
  );
  /* Computes the color {sVal(pix)}, blurring indicator {vBlr(pix)},
    height average {hAvg(pix)}, average normal {sNrm(pix)}, 
    and height variance {hVar(pix)} at
    the pixel {pix} on column {iPix.c[0]} and row {iPix.c[1]}, by ray-tracing 
    some scene or object with multiple rays. Returns these valus in
    {sVal_pix[0..NC-1]}, {*vBlr_pix_P}, {*hAvg_pix_P}, and {*hVar_pix_P}.
    
    Specifically, the procedure generates a set of {NS=samp.NS}
    sampoints (sub-sampling points) around the center of {pix} and
    throws a set {RAYS(p,iPix) of {KR=samp.KR} rays from each sampoint
    {p}, as explained in {multifok_raytrace_compute_point_properties}.
    The results of tracing the rays in {RAYS(p,iPix)} for each sampling
    point {p} are averaged yielding the average sampoint color
    {sVal(p)}, the average normal {sNrm(p)}, the average height
    {hAvg(p)}, and its variance {hVar(p)}. Then these values for all
    sampoints {p} in in the pixel are combined with the weights
    {wSmp[0..NS-1]} to obtain the pixel color {sVal(pix)}, the pixel
    normal {sNrm(pix)}, and the average and variance
    {hAvg(pix),hVar)(pix)} of the height.
    
    The blurring indicator {vBlr} is the weighted mean square radius of the 
    {X,Y} deviations of the ray hit points from the pixel center's {X,Y}.
    This value will always be 1.0 or more on account of the spread of the 
    sampoints {p}. 
    
    If {debug_pixel} procedure is called upon entry. If it returns true,
    the procedure prints to {stderr} info on the computation of the
    pixel's values. In that case, if {report_pixel} is not {NULL}, it
    calls it once before returning. In that case also, if {report_ray}
    is not {NULL}, calls it for each ray that contributed to the pixel's
    values. */

r3_t multifok_raytrace_choose_ray_direction(i2_t *iPix, i2_t *iSmp, r3_t *dRef, double zDep);
  /* Returns a unit ray direction vector {dRay} to ray-trace a scene from 
    some sampoint {pObs}.

    The ray will be pointing in the general diretion {dRef} (that is,
    {dot(dRay,dRef)>0}) but tilted relative to it by a random amount.

    The image plane is assumed to be coincide in scene space with the
    in-focus plane {F}, that passes through {pObs} and is perpendicular to
    the direction vector {dRef}. Let {Q} be a plane parallel to {F} and
    located {zDep} away from {F} in the direction {dRef}. Let {u} be the
    point {pObs + zDep*dRef} on {Q}, and {q[r]} be the point where the ray
    {R[r]} hits {Q}. The direction {dir(R[r])} of each ray {R[r]} will
    be chosen so that {q[r]} has a 2D Gaussian distribution on {Q}, with
    mean {u} and root-mean-square distance 1.0 from {u}.
    
    As a special case, if {zDep} is {+oo}, the result will be
    just {dRef}.
    
    The ray direction will be a pseudorandom function of the pixel
    indices {iPix} and the sampoint indices {iSmp} relative to that pixel.
    NOTE: this procefure re-seeds the random number generator with some
    function of {iPix} and {iSmp}.  */

/* IMPLEMENTATIONS */
     
r3_t multifok_raytrace_choose_ray_direction(i2_t *iPix, i2_t *iSmp, r3_t *dRef, double zDep)
  { if (zDep == +INF) { return *dRef; }
  
    demand(isfinite(zDep) && (zDep >= 1.0e-6), "invalid {zDep}");

    /* Choose two vectors perpendicular to {dRef}: */
    r3_t r0, r1;
    r3_throw_ortho_pair(dRef, &r0, &r1);
    
    /* Generate a  random vector {v[r]} from {pObs} to {q[r]}: */
    double s0 = dgaussrand(), s1 = dgaussrand();
    r3_t v; r3_scale(zDep, dRef, &v);
    r3_mix_in(s0, &r0, &v);
    r3_mix_in(s1, &r1, &v);
    
    /* Normalize and return: */
    r3_t dRay;
    double m = r3_dir(&v, &dRay);
    assert(m >= 0.9999*zDep);
    return dRay;
  }

void multifok_raytrace_paint_frame_rectangle
  ( multifok_frame_t *fr,
    int32_t xLo,
    int32_t xHi,
    int32_t yLo,
    int32_t yHi,
    multifok_raytrace_proc_t *trace_ray,
    multifok_raytrace_img_to_scene_map_t *map_point,
    r3_t *dRef,
    double zFoc,
    double zDep,
    multifok_sampling_t *samp,
    bool_t verbose,
    multifok_raytrace_debug_pred_t *debug_pixel,
    multifok_raytrace_report_ray_proc_t *report_ray,
    multifok_raytrace_report_pixel_proc_t report_pixel
  )
  {
    int32_t NC, NXF, NYF;
    float_image_get_size(fr->sVal, &NC, &NXF, &NYF);
    demand((0 <= xLo) && (xLo <= xHi) && (xHi < NXF), "invalid {xLo..xHi}");
    demand((0 <= yLo) && (yLo <= yHi) && (yHi < NYF), "invalid {yLo..yHi}");
    
    for (int32_t yPix = yLo; yPix <= yHi; yPix++)
      { for (int32_t xPix = xLo; xPix <= xHi; xPix++)
          { /* fprintf(stderr, "(%d,%d)", xPix, yPix); */
            i2_t iPix = (i2_t){{ xPix, yPix }};
            
            /* Compute the properties of pixel {iPix}: */
            r3_t sNrm_pix;
            float sVal_pix[NC];
            double vBlr_pix;
            double hAvg_pix;
            double hVar_pix;
 
            multifok_raytrace_compute_pixel_properties
              ( &iPix, trace_ray, map_point,
                xLo, xHi, yLo, yHi, samp, dRef, zFoc, zDep, 
                debug_pixel, report_ray, report_pixel,
                &sNrm_pix, NC, sVal_pix, &vBlr_pix, &hAvg_pix, &hVar_pix
              );
            double shrp_pix = 1/vBlr_pix;
            double hDev_pix = sqrt(hVar_pix);
              
            for (int32_t ic = 0;  ic < NC; ic++)
              { float_image_set_sample(fr->sVal, ic, xPix, yPix, sVal_pix[ic]); }
            float_image_set_sample(fr->hAvg, 0, xPix, yPix, (float)hAvg_pix);         
            float_image_set_sample(fr->hDev, 0, xPix, yPix, (float)hDev_pix); 
            float_image_set_sample(fr->shrp, 0, xPix, yPix, (float)shrp_pix); 
            for (int32_t j = 0; j < 3; j++)
              { float_image_set_sample(fr->sNrm, j, xPix, yPix, (float)sNrm_pix.c[j]); }        
          }
      }
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
    multifok_sampling_t *samp,
    bool_t verbose,
    multifok_raytrace_debug_pred_t *debug_pixel,
    multifok_raytrace_report_ray_proc_t *report_ray,
    multifok_raytrace_report_pixel_proc_t *report_pixel
  )
  {
    float_image_t *sVal = float_image_new(NC, NX, NY);
    float_image_t *shrp = float_image_new(1, NX, NY);
    float_image_t *hAvg = float_image_new(1, NX, NY);
    float_image_t *hDev = float_image_new(1, NX, NY);
    float_image_t *sNrm = float_image_new(3, NX, NY);
  
    multifok_frame_t *fr = multifok_frame_from_images(sVal,shrp,hAvg,hDev,sNrm, zFoc,zDep);
    
    multifok_raytrace_paint_frame_rectangle
      ( fr, 0, NX-1, 0, NY-1, 
        trace_ray, map_point, dRef, zFoc, zDep, samp,
        verbose, debug_pixel, report_ray, report_pixel
      );
      
    return fr;

  }

void multifok_raytrace_compute_pixel_properties
  ( i2_t *iPix, 
    multifok_raytrace_proc_t *trace_ray,
    multifok_raytrace_img_to_scene_map_t *map_point,
    int32_t xLo,
    int32_t xHi,
    int32_t yLo,
    int32_t yHi,
    multifok_sampling_t *samp,
    r3_t *dRef,
    double zFoc,
    double zDep,
    multifok_raytrace_debug_pred_t *debug_pixel,
    multifok_raytrace_report_ray_proc_t report_ray,
    multifok_raytrace_report_pixel_proc_t report_pixel,
    r3_t *sNrm_pix_P,
    int32_t NC,
    float sVal_pix[],
    double *vBlr_pix_P,
    double *hAvg_pix_P,
    double *hVar_pix_P
  )
  { 
    uint32_t NS = samp->NS; /* Count of sampoints per pixel. */
    uint32_t KR = samp->KR; /* Rays per sampoint. */
    demand(NS >= 1, "invalid pixel sampoint count {NS}.");
    demand(KR >= 1, "invalid ray count per sampoint {KR}");
    bool_t deb_pix = (debug_pixel == NULL ? FALSE : debug_pixel(iPix));
    if (deb_pix) 
      { fprintf(stderr, "  " SHARPS "\n");
        fprintf(stderr, "  computing properties of pixel iPix = ");
        i2_gen_print(stderr, iPix, "%d", "[ ", " ", " ]");
        fprintf(stderr, " with NS = %d KR = %d\n", NS, KR);
      }
      
    int32_t xPix = iPix->c[0];
    int32_t yPix = iPix->c[1];
      
    /* Scan sampoints, accumulate: */
    r3_t sum_w_sNrm = (r3_t){{0,0,0}}; /* Sum of {w(p)*sNrm(p)}. */
    double sum_w_sVal[NC]; /* Sum of {w(p)*sVal(p)}. */
    for (int32_t ic = 0;  ic < NC; ic++) { sum_w_sVal[ic] = 0.0; }
    double sum_w_vBlr = 0.0; /* Sum of {w(p)*vBlr(R)}. */
    double sum_w_hAvg = 0.0; /* Sum of {w(p)*hAvg(p)}. */
    double sum_w = 0.0; /* Sum of {w(p)}. */
    double hAvg_pt[NS]; /* Saved {zAve(p)} for all sampoints. */
    double hVar_pt[NS]; /* Saved {hVar(p)} for all sampoints. */
    double step = samp->step; /* Spacing between sampoints. */
    for (uint32_t ks = 0;  ks < NS; ks++)
      { bool_t deb_pt = deb_pix;
        i2_t *iSmp = &(samp->iSmp[ks]);
        double wk = samp->wSmp[ks];
        double x_img = xPix + 0.5 + iSmp->c[0]*step;
        double y_img = yPix + 0.5 + iSmp->c[1]*step;
        r2_t pImg = (r2_t){{ x_img, y_img }};
        r3_t pSmp; map_point(&pImg, &pSmp);
        r3_t pRay = pSmp;
        if (deb_pt) 
          { fprintf(stderr, "    " EQUALS "\n");
            fprintf(stderr, "    computing properties of sample point %d", ks);
            i2_gen_print(stderr, iSmp, "%+3d", "  = ( ", " ", " )\n");
            r2_gen_print(stderr, &pImg, "%12.6f", "    image coords ( ", " ", " )\n");
            r3_gen_print(stderr, &pRay, "%12.6f", "    scene coords ( ", " ", " )\n");
            fprintf(stderr, "    weight = %10.8f\n", wk);
          }

        r3_t sNrm_pt;
        float sVal_pt[NC];
        double vBlr_pt;
        multifok_raytrace_compute_point_properties
          ( trace_ray, &pRay, KR,
            dRef, zFoc, zDep, 
            deb_pt, iPix, iSmp, step, wk, report_ray, 
            &sNrm_pt, NC, sVal_pt, &vBlr_pt, &(hAvg_pt[ks]), &(hVar_pt[ks])
          );
        double hDev_pt = sqrt(hVar_pt[ks]);
        if (deb_pt) 
          { fprintf(stderr, "    sample point vBlr = %12.8f hAvg = %12.6f  hDev = %16.12f", vBlr_pt, hAvg_pt[ks], hDev_pt);
            r3_gen_print(stderr, &sNrm_pt, "%+7.4f", " normal = ( ", " ", " )\n");
            fprintf(stderr, "    " EQUALS "\n");
          }
        /* Accumulate: */
        for (int32_t j = 0; j < 3; j++) { sum_w_sNrm.c[j] += wk*sNrm_pt.c[j]; }
        for (int32_t ic = 0;  ic < NC; ic++) 
          { sum_w_sVal[ic] += wk*(double)sVal_pt[ic]; }
        assert(vBlr_pt >= 0.0);
        double r2Smp = step*step*(double)i2_norm_sqr(iSmp); /* Contrib of sampoint pos to {vBlr}. */
        sum_w_vBlr += wk*(vBlr_pt + r2Smp);
        sum_w_hAvg += wk*hAvg_pt[ks];
        sum_w += wk;
      }
    
    assert(sum_w > 0);
    if (deb_pix) 
      { fprintf(stderr, "  sum_w_vBlr = %12.8f", sum_w_vBlr); 
        fprintf(stderr, " sum_w_hAvg = %12.6f", sum_w_hAvg); 
        r3_gen_print(stderr, &sum_w_sNrm, "%+8.5f", " sum_w_sNrm = ( ", " ", " )");
        fprintf(stderr, " sum_w = %16.12f\n", sum_w); 
      }

    /* Compute average pixel normal {sNrm_pix[0..NC-1]}: */
    r3_t sNrm_pix;
    for (int32_t j = 0; j < 3; j++) { sNrm_pix.c[j] = sum_w_sNrm.c[j]/sum_w; }
    (void)r3_dir(&(sNrm_pix), &(sNrm_pix));
    if ((! isfinite(sNrm_pix.c[0])) || (! isfinite(sNrm_pix.c[1])) || (! isfinite(sNrm_pix.c[2])))
      { sNrm_pix = (r3_t){{0,0,0}}; }

    /* Compute average pixel color {sVal_pix[0..NC-1]}: */
    for (int32_t ic = 0;  ic < NC; ic++) 
      { sVal_pix[ic] = (float)(sum_w_sVal[ic]/sum_w); }
    
    /* Compute average blurring indicator {vBlr_pix}: */
    double vBlr_pix = sum_w_vBlr/sum_w;
    assert((! isnan(vBlr_pix)) && (vBlr_pix >= 0));
    
    /* Sharpness indicator {shrp} is {1/vBlr}: */
    double shrp_pix = 1.0/vBlr_pix;
    
    /* Compute pixel height average {hAvg_pix}: */
    double hAvg_pix = sum_w_hAvg/sum_w;
    
    /* Compute pixel height variance {hVar_pix}: */
    double sum_w_dz2 = 0.0;
    for (uint32_t ks = 0;  ks < NS; ks++)
      { double dz = hAvg_pt[ks] - hAvg_pix; 
        double wk = samp->wSmp[ks];
        sum_w_dz2 += wk*(dz*dz + hVar_pt[ks]);
      }
    double hVar_pix = sum_w_dz2/sum_w;
    double hDev_pix = sqrt(hVar_pix);
    
    if (deb_pix) 
      { fprintf(stderr, "  pixel [%4d,%4d]", iPix->c[0], iPix->c[1]);
        fprintf(stderr, "  vBlr = %12.8f hAvg = %12.6f hDev = %16.12f", vBlr_pix, hAvg_pix, hDev_pix);
        r3_gen_print(stderr, &sNrm_pix, "%+7.4f", " sNrm = ( ", " ", " )\n");
        if (report_pixel != NULL)
          { r2_t ctr2 = (r2_t){{ xPix + 0.5, yPix + 0.5 }};
            r3_t pCtr; map_point(&ctr2, &pCtr); 
            report_pixel
              ( iPix, &pCtr, zFoc, zDep, 
                shrp_pix, hAvg_pix, hDev_pix, &sNrm_pix, NC, sVal_pix 
              );
          }
      }

    /* Return pixel data: */
    (*sNrm_pix_P) = sNrm_pix;
    (*vBlr_pix_P) = vBlr_pix;
    (*hAvg_pix_P) = hAvg_pix;
    (*hVar_pix_P) = hVar_pix;
    
    if (deb_pix) { fprintf(stderr, "  " SHARPS "\n"); }
  }
 
void multifok_raytrace_compute_point_properties
  ( multifok_raytrace_proc_t *trace_ray,
    r3_t *pObs,
    uint32_t KR,
    r3_t *dRef,
    double zFoc,
    double zDep,
    bool_t deb_pt,
    i2_t *iPix,
    i2_t *iSmp,
    double step,
    double wSmp,
    multifok_raytrace_report_ray_proc_t report_ray,
    r3_t *sNrm_pt_P,
    int32_t NC,
    float sVal_pt[],
    double *vBlr_pt_P,
    double *hAvg_pt_P,
    double *hVar_pt_P
  )
  {
    if (deb_pt) { fprintf(stderr, "    entering {multifok_raytrace_compute_point_properties} ...\n"); }

    demand(KR >= 1, "invalid {KR}"); 
    
    r3_t sum_w_sNrm = (r3_t){{0,0,0}}; /* Sum of {wRay(R)*sNrm(R)}. */
    double sum_w_sVal[NC]; /* Sum of {wRay(R)*sVal(R)}. */
    for (int32_t ic = 0;  ic < NC; ic++) { sum_w_sVal[ic] = 0.0; }
    double sum_w_vBlr = 0.0; /* Sum of {wRay(R)*vBlr(R)}. */
    double sum_w_hHit = 0.0; /* Sum of {wRay(R)*hHit(R)}. */
    double sum_w = 1.0e-200; /* Sum of {wRay(R)}, biased to avoid {0/0}. */
    double hHit_ray[KR]; /* Stores heights of rays used. */
    double w_ray[KR];    /* Stores weights of rays used. */
    /* Cast up to {2*KR} rays trying to get {KR} hits: */
    
    /* Restart the random number generator: */
    int32_t rseed = 417 + 4615*iPix->c[0] + 17*iPix->c[1];
    rseed = 121*rseed + 11*iSmp->c[0] + iSmp->c[1];
    srandom((uint32_t)rseed);

    for (uint32_t ir = 0;  ir < KR; ir++)
      { bool_t deb_ray = deb_pt;
        /* Get the ray's direction {dRay}: */
        r3_t dRay = multifok_raytrace_choose_ray_direction(iPix, iSmp, dRef, zDep);
        if (deb_ray)
          { fprintf(stderr, "      " DASHES "\n");
            fprintf(stderr, "      computing properties of ray %4d", ir);
            r3_gen_print(stderr, &(dRay), "%+9.6f", "      direction = ( ", " ", " )\n");
          }

        /* Trace a ray throug {pObs} with direction {dRay}: */
        float sVal_ray[NC];
        r3_t pHit, sNrm;
        trace_ray(pObs, &dRay, deb_ray, &pHit, &sNrm, NC, sVal_ray);
        assert(isfinite(pHit.c[0]) && isfinite(pHit.c[1]) && isfinite(pHit.c[2]));

        /* Save weight of this ray: */
        w_ray[ir] = 1.0;

        /* Compute and save {hHit_ray} = height of {pHit} from the {zFoc=0} image plane: */
        r3_t u; r3_sub(&pHit, pObs, &u);
        double dHit = r3_dot(&u, dRef); /* Depth of {pHit} away from curr image plane. */
        assert(isfinite(dHit));
        hHit_ray[ir] = zFoc - dHit;

        /* Compute {vBlr_ray}: */
        r3_t qRef;  /* Point on ray {pObs,dRef} with same depth as {pHit}. */
        r3_mix(1.0, pObs, dHit, dRef, &qRef);
        double vBlr_ray = r3_dist_sqr(&pHit, &qRef);
        assert(isfinite(vBlr_ray));

        if (deb_ray)
          { fprintf(stderr, "        ray vBlr = %12.8f hHit = %12.6f weight = %16.12f", vBlr_ray, hHit_ray[ir], w_ray[ir]);
            r3_gen_print(stderr, &sNrm, "%+8.5f", " sNrm = ( ", " ", " )\n");
            if (report_ray != NULL)
              { report_ray(iPix, iSmp, step, wSmp, pObs, &dRay, w_ray[ir], &pHit, hHit_ray[ir], vBlr_ray, &sNrm, NC, sVal_ray); }
            fprintf(stderr, "      " DASHES "\n");
          }

        /* Accumulate the ray properties onto the point properties: */
        for (int32_t j = 0; j < 3; j++) { sum_w_sNrm.c[j] += w_ray[ir]*sNrm.c[j]; }
        for (int32_t ic = 0;  ic < NC; ic++)
          { double sVali = fmax(0.0, fmin(1.0, sVal_ray[ic]));
            assert(isfinite(sVali));
            sum_w_sVal[ic] += w_ray[ir]*(double)sVali;
          }
        sum_w_vBlr += w_ray[ir]*vBlr_ray;
        sum_w_hHit += w_ray[ir]*hHit_ray[ir];
        sum_w += w_ray[ir];
      }

    assert(sum_w > 0);
    if (deb_pt) 
      { fprintf(stderr, "    sum_w_vBlr = %12.8f", sum_w_vBlr); 
        fprintf(stderr, " sum_w_hHit = %12.6f", sum_w_hHit); 
        r3_gen_print(stderr, &sum_w_sNrm, "%+8.5f", " sum_w_sNrm = ( ", " ", " )");
        fprintf(stderr, " sum_w = %16.12f\n", sum_w); 
      }
    
    /* Compute the mean normal {sNrm_pt} at sampoint: */
    r3_t sNrm_pt;
    for (int32_t j = 0; j < 3; j++) { sNrm_pt.c[j] = sum_w_sNrm.c[j]/sum_w; }
    (void)r3_dir(&(sNrm_pt), &(sNrm_pt));
    if ((! isfinite(sNrm_pt.c[0])) || (! isfinite(sNrm_pt.c[1])) || (! isfinite(sNrm_pt.c[2])))
      { sNrm_pt = (r3_t){{0,0,0}}; }
     
    /* Compute the mean color {sVal_pt} at sampoint: */
    for (int32_t ic = 0;  ic < NC; ic++) { sVal_pt[ic] = (float)(sum_w_sVal[ic]/sum_w); }
    
    /* Compute the mean blur {vBlr_pt} at sampoint: */
    double vBlr_pt = sum_w_vBlr/sum_w;
    assert((! isnan(vBlr_pt)) && (vBlr_pt >= 0.0));

    /* Compute the mean height {hAvg_pt} at sampoint: */
    double hAvg_pt = sum_w_hHit/sum_w;
    
    /* Compute the height variance {hVar_pt} at sampoint: */
    double sum_w_dz2 = 0;
    for (uint32_t ir = 0;  ir < KR; ir++)
      { double dz = hHit_ray[ir] - hAvg_pt;
        sum_w_dz2 += w_ray[ir]*dz*dz;
      }
    double hVar_pt = sum_w_dz2/sum_w;
    
    if (deb_pt) 
      { fprintf(stderr, "    vBlr_pt = %12.8f hAvg_pt = %12.6f hVar_pt = %12.6f\n", vBlr_pt, hAvg_pt, hVar_pt); 
        r3_gen_print(stderr, &sNrm_pt, "%+8.5f", "    sNrm_pt = ( ", " ", " )\n");
        fprintf(stderr, "    exiting {multifok_raytrace_compute_point_properties}\n");
      }

    /* Return results: */
    (*sNrm_pt_P) = sNrm_pt;
    (*vBlr_pt_P) = vBlr_pt;
    (*hAvg_pt_P) = hAvg_pt;
    (*hVar_pt_P) = hVar_pt;
  }

void multifok_raytrace_show_ray_data
  ( FILE *wr,
    int32_t indent,
    i2_t *iPix,
    double pixSize,
    i2_t *iSmp,
    double step,
    double wSmp,
    r3_t *pRay,
    r3_t *dRay, 
    double wRay,
    r3_t *pHit, 
    double hHit,
    double vBlr,
    r3_t *sNrm,
    int32_t NC,
    float sVal[]
  )
  {
    fprintf(wr, "%*s", indent, "");
    i2_gen_print(wr, iPix, "%4d", "pixel [", ",", "]");
    fprintf(wr, " size %.6f × %.6f", pixSize, pixSize);
    fprintf(wr, " sub-sampling point ( %+3d %+3d )", iSmp->c[0], iSmp->c[1]);
    fprintf(wr, " step %10.8f", step);
    fprintf(wr, " weight = %.8f\n", wSmp);
    fprintf(wr, "%*s", indent, "");
    r3_gen_print(wr, pRay, "%10.4f", "ray start = ( ", " ", " )");
    r3_gen_print(wr, dRay, "%+9.6f", " dir = ( ", " ", " )");
    fprintf(wr, " weight = %.8f\n", wRay);
    fprintf(wr, "%*s", indent, "");
    r3_gen_print(wr, pHit, "%10.4f", "pHit = ( ", " ", " )");
    fprintf(wr, " vBlr = %10.4f hHit = %10.4f\n", vBlr, hHit);
    fprintf(wr, "%*s", indent, "");
    r3_gen_print(wr, sNrm, "%+8.5f", "sNrm = ( ", " ", " )\n");
    fprintf(wr, "%*s", indent, "");
    fprintf(wr, "sVal = ( ");
    for (uint32_t c = 0; c < 3; c++)
      { fprintf(wr, " %7.5f", (c < NC ? sVal[c] : 0.0)); }
    fprintf(wr, " )\n");
    fflush(wr);
  }
    
void multifok_raytrace_write_ray_data
  ( FILE *wr,
    i2_t *iPix,
    double pixSize,
    i2_t *iSmp,
    double step,
    double wSmp,
    r3_t *pRay,
    r3_t *dRay, 
    double wRay,
    r3_t *pHit, 
    double hHit,
    double vBlr,
    r3_t *sNrm,
    int32_t NC,
    float sVal[]
  )
  {
    fprintf(wr, "%4d %4d", iPix->c[0], iPix->c[1]);
    fprintf(wr, " %9.6f", pixSize);
    fprintf(wr, "  %+3d %+3d", iSmp->c[0], iSmp->c[1]);
    fprintf(wr, "  %8.5f", step);
    fprintf(wr, "  %10.8f", wSmp);
    r3_gen_print(wr, pRay, "%12.6f", "  ", " ", " ");
    r3_gen_print(wr, dRay, "%+9.6f", "  ", " ", " ");
    fprintf(wr, "   %10.8f", wRay);
    r3_gen_print(wr, pHit, "%12.6f", "  ", " ", " ");
    fprintf(wr, "  %12.8f %12.6f", vBlr, hHit);
    r3_gen_print(wr, sNrm, "%+8.5f", "  ", " ", " ");
    for (uint32_t c = 0; c < 3; c++)
      { fprintf(wr, " %7.5f", (c < NC ? sVal[c] : 0.0)); }
    fprintf(wr, "\n");
    fflush(wr);
  }

void multifok_raytrace_show_pixel_data
  ( FILE *wr,
    int32_t indent,
    i2_t *iPix,
    double pixSize,
    r3_t *pCtr,
    double zFoc,
    double zDep,
    double shrp,
    double hAvg,
    double hDev,
    r3_t *sNrm,
    int32_t NC,
    float sVal[]
  )
  {
    fprintf(wr, "%*s", indent, "");
    i2_gen_print(wr, iPix, "%4d", "pixel [", ",", "]");
    fprintf(wr, " size %.6f × %.6f", pixSize, pixSize);
    r3_gen_print(wr, pCtr, "%10.4f", " center ( ", " ", " )\n");
    fprintf(wr, "%*s", indent, "");
    fprintf(wr, "zFoc = %12.8f  zDep = %12.6f\n", zFoc, zDep);
    fprintf(wr, "shrp = %10.4f  hAvg = %10.4f  hDev = %10.4f", shrp, hAvg, hDev);
    r3_gen_print(wr, sNrm, "%+8.5f", " sNrm = ( ", " ", " )");
    fprintf(wr, "  color = (");
    for (uint32_t c = 0; c < NC; c++) { fprintf(wr, " %7.5f", sVal[c]); }
    fprintf(wr, "\n\n");
    fflush(wr);
  }

void multifok_raytrace_write_pixel_data
  ( FILE *wr,
    i2_t *iPix,
    double pixSize,
    r3_t *pCtr,
    double zFoc,
    double zDep,
    double shrp,
    double hAvg,
    double hDev,
    r3_t *sNrm,
    int32_t NC,
    float sVal[]
  )
  {
    fprintf(wr, "%4d %4d", iPix->c[0], iPix->c[1]);
    fprintf(wr, " %9.6f", pixSize);
    r3_gen_print(wr, pCtr, "%12.6f", "  ", " ", " ");
    fprintf(wr, "  %12.8f %12.6f", zFoc, zDep);
    fprintf(wr, "  %12.8f", shrp);
    fprintf(wr, "  %12.8f %12.6f", hAvg, hDev);
    fprintf(wr, "  %1d", NC);
    r3_gen_print(wr, sNrm, "%+8.5f", "  ", " ", " ");
    for (uint32_t c = 0; c < 3; c++)
      { fprintf(wr, " %7.5f", (c < NC ? sVal[c] : 0.0)); }
    fprintf(wr, "\n");
    fflush(wr);
  }
