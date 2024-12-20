/* See {multifok_raytrace.h}. */
/* Last edited on 2024-12-15 20:36:32 by stolfi */

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
  
#define EQUALS "============================================================"
#define SHARPS "############################################################"

void multifok_raytrace_compute_pixel_properties
  ( i2_t *iPix, 
    multifok_raytrace_proc_t *trace_ray,
    multifok_raytrace_img_to_scene_map_t *map_point,
    int32_t NX,
    int32_t NY,
    multifok_sampling_t *samp,
    r3_t *dRef,
    double zFoc,
    double zDep,
    multifok_raytrace_debug_pred_t *debug_pixel,
    multifok_raytrace_report_ray_proc_t report_ray,
    multifok_raytrace_report_pixel_proc_t report_pixel,
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
    
    Specifically, the procedure generates a set of {NS=samp.NS}
    sampoints (sub-sampling points) around the center of {pix} and throws
    a set {RAYS(p,iPix) of {KR=samp.KR} rays from each sampoint {p}, as explained in
    {multifok_raytrace_compute_point_properties}.
    The results of tracing the rays in {RAYS(p,iPix)} for each sampling
    point {p} are averaged yielding the average sampoint color
    {colr(p)}, the average height {hAvg(p)}, and its variance {hVar(p)}.
    Then these values for all sampoints {p} in in the ]pixel are
    combined with the weights {wSmp[0..NS-1]} to obtain the pixel color
    {colr(pix)} and the average and variance {hAvg(pix),hVar)(pix)} of
    the height.
    
    The blurring indicator {vBlr} is the weighted mean square radius of the 
    {X,Y} deviations of the ray hit points from the pixel center's {X,Y}.
    This value will always be 1.0 or more on account of the spread of the 
    sampoints {p}. 
    
    If {deb_pix} is true, prints to {stderr} info on the computation of
    the pixel's values. In that case, if {report_ray} is not {NULL},
    calls it for each ray that contributed to those values. */

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
    
    for (int32_t yPix = 0;  yPix < NY; yPix++)
      { for (int32_t xPix = 0;  xPix < NX; xPix++)
          { /* fprintf(stderr, "(%d,%d)", xPix, yPix); */
            i2_t iPix = (i2_t){{ xPix, yPix }};
            
            /* Compute the properties of pixel {iPix}: */
            float colr_pix[NC];
            double vBlr_pix;
            double hAvg_pix;
            double hVar_pix;
 
            multifok_raytrace_compute_pixel_properties
              ( &iPix, trace_ray, map_point,
                NX, NY, samp, dRef, zFoc, zDep, 
                debug_pixel, report_ray, report_pixel,
                NC, colr_pix, &vBlr_pix, &hAvg_pix, &hVar_pix
              );
            double shrp_pix = 1/vBlr_pix;
            double hDev_pix = sqrt(hVar_pix);
              
            for (int32_t ic = 0;  ic < NC; ic++)
              { float_image_set_sample(sVal, ic, xPix, yPix, colr_pix[ic]); }
            float_image_set_sample(hAvg, 0, xPix, yPix, (float)hAvg_pix);         
            float_image_set_sample(hDev, 0, xPix, yPix, (float)hDev_pix); 
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
    multifok_sampling_t *samp,
    r3_t *dRef,
    double zFoc,
    double zDep,
    multifok_raytrace_debug_pred_t *debug_pixel,
    multifok_raytrace_report_ray_proc_t report_ray,
    multifok_raytrace_report_pixel_proc_t report_pixel,
    int32_t NC,
    float colr_pix[],
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
        fprintf(stderr, "  computing pixel iPix = ");
        i2_gen_print(stderr, iPix, "%d", "[ ", " ", " ]");
        fprintf(stderr, " with NS = %d KR = %d\n", NS, KR);
      }
      
    int32_t xPix = iPix->c[0];
    int32_t yPix = iPix->c[1];
      
    /* Scan sampoints, accumulate: */
    double sum_w_colr[NC]; /* Sum of {w(p)*colr(p)}. */
    for (int32_t ic = 0;  ic < NC; ic++) { sum_w_colr[ic] = 0.0; }
    double sum_w_vBlr = 0.0; /* Sum of {w(p)*vBlr(R)}. */
    double sum_w_hAvg = 0.0; /* Sum of {w(p)*hAvg(p)}. */
    double sum_w = 0.0; /* Sum of {w(p)}. */
    double hAvg_pt[NS]; /* Saved {zAve(p)} for all sampoints. */
    double hVar_pt[NS]; /* Saved {hVar(p)} for all sampoints. */
    double step = samp->step; /* Spacing between sampoints. */
    for (uint32_t ks = 0;  ks < NS; ks++)
      { i2_t *iSmp = &(samp->iSmp[ks]);
        double wk = samp->wSmp[ks];
        double x_img = xPix + 0.5 + iSmp->c[0]*step;
        double y_img = yPix + 0.5 + iSmp->c[1]*step;
        r2_t pSmp = (r2_t){{ x_img, y_img }};
        r3_t pRay; map_point(&pSmp, &pRay);
        if (deb_pix) 
          { fprintf(stderr, "    sample point %d", ks);
            i2_gen_print(stderr, iSmp, "%+3d", "  = ( ", " ", " )");
            r2_gen_print(stderr, &pSmp, "%12.6f", "  = image ( ", " ", " )");
            r3_gen_print(stderr, &pRay, "%12.6f", "  = scene ( ", " ", " )");
            fprintf(stderr, " weight = %10.8f\n", wk);
          }
        float colr_pt[NC];
        
        double vBlr_pt;
        multifok_raytrace_compute_point_properties
          ( trace_ray, &pRay, KR,
            dRef, zFoc, zDep, 
            deb_pix, iPix, iSmp, step, wk, report_ray, 
            NC, colr_pt, &vBlr_pt, &(hAvg_pt[ks]), &(hVar_pt[ks])
          );
        double hDev_pt = sqrt(hVar_pt[ks]);
        if (deb_pix) { fprintf(stderr, "    sample point vBlr = %12.8f hAvg = %12.6f  hDev = %16.12f\n", vBlr_pt, hAvg_pt[ks], hDev_pt); }
        /* Accumulate: */
        double r2Smp = step*step*(double)i2_norm_sqr(iSmp); /* Contrib of sampoint pos to {vBlr}. */
        for (int32_t ic = 0;  ic < NC; ic++) 
          { sum_w_colr[ic] += wk*(double)colr_pt[ic]; }
        assert(vBlr_pt >= 0.0);
        sum_w_vBlr += wk*(vBlr_pt + r2Smp);
        sum_w_hAvg += wk*hAvg_pt[ks];
        sum_w += wk;
      }
    
    /* Compute average pixel color {colr_pix[0..NC-1]}: */
    for (int32_t ic = 0;  ic < NC; ic++) 
      { colr_pix[ic] = (float)(sum_w_colr[ic]/sum_w); }
    
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
        if (report_pixel != NULL)
          { r2_t ctr2 = (r2_t){{ xPix + 0.5, yPix + 0.5 }};
            r3_t pCtr; map_point(&ctr2, &pCtr); 
            report_pixel
              ( iPix, &pCtr, zFoc, zDep, 
                shrp_pix, hAvg_pix, hDev_pix, NC, colr_pix 
              );
          }
      }

    /* Return pixel data: */
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
    bool_t debug,
    i2_t *iPix,
    i2_t *iSmp,
    double step,
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
        r3_gen_print(stderr, pObs, "%12.6f", "      computing properties of scene point ( ", " ", " )");
        fprintf(stderr, " point sampling weight = %.8f\n", wSmp);
      }

    demand(KR >= 1, "invalid {KR}"); 
    
    double sum_w_colr[NC]; /* Sum of {wRay(R)*colr(R)}. */
    for (int32_t ic = 0;  ic < NC; ic++) { sum_w_colr[ic] = 0.0; }
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

    for (uint32_t jr = 0;  jr < KR; jr++)
      { /* Get the ray's direction {dRay}: */
        r3_t dRay = multifok_raytrace_choose_ray_direction(iPix, iSmp, dRef, zDep);
        if (debug)
          { fprintf(stderr, "        ray %4d direction = ", jr);
            r3_gen_print(stderr, &(dRay), "%+9.6f", "( ", " ", " )\n");
          }

        /* Trace a ray throug {pObs} with direction {dRay}: */
        float colr_ray[NC];
        r3_t pHit;
        trace_ray(pObs, &dRay, debug, &pHit, NC, colr_ray);
        assert(isfinite(pHit.c[0]) && isfinite(pHit.c[1]) && isfinite(pHit.c[2]));

        /* Save weight of this ray: */
        w_ray[jr] = 1.0;

        /* Compute and save {hHit_ray} = height of {pHit} from the {zFoc=0} image plane: */
        r3_t u; r3_sub(&pHit, pObs, &u);
        double dHit = r3_dot(&u, dRef); /* Depth of {pHit} away from curr image plane. */
        assert(isfinite(dHit));
        hHit_ray[jr] = zFoc - dHit;

        /* Compute {vBlr_ray}: */
        r3_t qRef;  /* Point on ray {pObs,dRef} with same depth as {pHit}. */
        r3_mix(1.0, pObs, dHit, dRef, &qRef);
        double vBlr_ray = r3_dist_sqr(&pHit, &qRef);
        assert(isfinite(vBlr_ray));

        if (debug)
          { if (report_ray != NULL)
              { report_ray(iPix, iSmp, step, wSmp, pObs, &dRay, w_ray[jr], &pHit, hHit_ray[jr], vBlr_ray); }
            else
              { fprintf(stderr, "        vBlr = %12.8f hHit = %12.6f w = %16.12f\n", vBlr_ray, hHit_ray[jr], w_ray[jr]); }
          }

        /* Accumulate the ray properties onto the point properties: */
        for (int32_t ic = 0;  ic < NC; ic++)
          { double colri = fmax(0.0, fmin(1.0, colr_ray[ic]));
            assert(isfinite(colri));
            sum_w_colr[ic] += w_ray[jr]*(double)colri;
          }
        sum_w_vBlr += w_ray[jr]*vBlr_ray;
        sum_w_hHit += w_ray[jr]*hHit_ray[jr];
        sum_w += w_ray[jr];
      }

    assert(sum_w > 0);
    if (debug) {fprintf(stderr, "      sum_w_vBlr = %12.8f sum_w_hHit = %12.6f sum_w = %16.12f\n", sum_w_vBlr, sum_w_hHit, sum_w); }
    for (int32_t ic = 0;  ic < NC; ic++) { colr_pt[ic] = (float)(sum_w_colr[ic]/sum_w); }
    double vBlr_pt = sum_w_vBlr/sum_w;
    assert((! isnan(vBlr_pt)) && (vBlr_pt >= 0.0));
    double hAvg_pt = sum_w_hHit/sum_w;
    
    /* Compute Z height variance: */
    double sum_w_dz2 = 0;
    for (uint32_t jr = 0;  jr < KR; jr++)
      { double dz = hHit_ray[jr] - hAvg_pt;
        sum_w_dz2 += w_ray[jr]*dz*dz;
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
    double vBlr
  )
  {
    fprintf(wr, "%*s", indent, "");
    i2_gen_print(wr, iPix, "%4d", "pixel [", ",", "]");
    fprintf(wr, " size %.6f × %.6f", pixSize, pixSize);
    fprintf(wr, " sub-sampling point ( %+3d %+3d )", iSmp->c[0], iSmp->c[1]);
    fprintf(wr, " step %10.8f", step);
    fprintf(wr, " weight = %.8f\n", wSmp);
    fprintf(wr, "%*s", indent, "");
    r3_gen_print(wr, pRay, "%10.4f", " = scene ( ", " ", " )");
    r3_gen_print(wr, dRay, "%+9.6f", " dir = ( ", " ", " )");
    fprintf(wr, " weight = %.8f", wRay);
    r3_gen_print(wr, pHit, "%10.4f", " hit = ( ", " ", " )");
    fprintf(wr, " vBlr = %10.4f hHit = %10.4f\n", vBlr, hHit);
    fprintf(wr, "\n");
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
    double vBlr
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
    int32_t NC,
    float colr[]
  )
  {
    fprintf(wr, "%*s", indent, "");
    i2_gen_print(wr, iPix, "%4d", "pixel [", ",", "]");
    fprintf(wr, " size %.6f × %.6f", pixSize, pixSize);
    r3_gen_print(wr, pCtr, "%10.4f", " center ( ", " ", " )\n");
    fprintf(wr, "%*s", indent, "");
    fprintf(wr, "zFoc = %12.8f  zDep = %12.6f\n", zFoc, zDep);
    fprintf(wr, "shrp = %10.4f  hAvg = %10.4f  hDev = %10.4f", shrp, hAvg, hDev);
    fprintf(wr, "  color = (");
    for (uint32_t c = 0; c < NC; c++) { fprintf(wr, " %7.5f", colr[c]); }
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
    int32_t NC,
    float colr[]
  )
  {
    fprintf(wr, "%4d %4d", iPix->c[0], iPix->c[1]);
    fprintf(wr, " %9.6f", pixSize);
    r3_gen_print(wr, pCtr, "%12.6f", "  ", " ", " ");
    fprintf(wr, "  %12.8f %12.6f", zFoc, zDep);
    fprintf(wr, "  %12.8f", shrp);
    fprintf(wr, "  %12.8f %12.6f", hAvg, hDev);
    fprintf(wr, "  %1d", NC);
    for (uint32_t c = 0; c < 3; c++)
      { fprintf(wr, " %7.5f", (c < NC ? colr[c] : 0.0)); }
    fprintf(wr, "\n");
    fflush(wr);
  }
