/* See {msm_image_tools.h} */
/* Last edited on 2022-10-20 06:38:32 by stolfi */

#define msm_image_tools_C_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>

#include <vec.h>
#include <float_image.h>
#include <frgb.h>
#include <jsmath.h>
#include <affirm.h>

#include <msm_basic.h>
#include <msm_rung.h>
#include <msm_seq_desc.h>
#include <msm_pairing.h>
#include <msm_cand.h>
#include <msm_cand_vec.h>
#include <msm_dyn.h>
#include <msm_image.h>
#include <msm_image_paint.h>

#include <msm_image_tools.h>

void msm_image_seq_seq_score_write_named
  ( msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1, 
    msm_rung_score_proc_t *score, 
    msm_pairing_t *pr, 
    char *name, 
    char *tag
  )
  { int32_t n0 = seq0->size;
    int32_t n1 = seq1->size;
    msm_image_t *gim = msm_image_alloc(1, n0, n1, TRUE);
    fprintf(stderr, "  image size = %d × %d  scale factor = %d\n", n0, n1, gim->scale);
    /* Paint the sample-to-sample distances: */
    msm_image_seq_seq_score_paint(seq0, seq1, score, gim, 0, 0);
    /* Find pixel value range: */
    double lov, hiv;
    msm_image_compute_min_max_range(gim, 0, &lov, &hiv);
    /* Compute the maximum abs value of data: */
    double maxf = (fabs(lov) > fabs(hiv) ? fabs(lov) : fabs(hiv));
    /* Reduce that value to better show the interesting part of the plot: */
    maxf = 0.25*maxf;
    /* Colorize the image: */
    msm_image_t *cim = msm_image_colorize(gim, maxf);
    if (pr != NULL)
      { frgb_t clr = (frgb_t){{ 0.000f, 0.700f, 0.000f }}; /* Dark green. */
        msm_image_pairing_paint(pr, n0, n1, clr.c, FALSE, cim, 0, 0);
      }
    /* Write image: */
    msm_image_write_as_pnm(cim, 0, 1, name, tag);
    /* Reclaim image space: */
    float_image_free(gim->fim); free(gim);
    float_image_free(cim->fim); free(cim);
  }

void msm_image_cand_write_named(msm_cand_t *cd, float v, char *name, char *tag)
  { /* Get sequence sizes {n0,n1}: */
    int32_t n0 = cd->seq[0].size;
    int32_t n1 = cd->seq[1].size;
    /* Choose the unscaled image size {nv0,nv1 = n0,n1}, for folded plot: */
    int32_t nv0 = n0, nv1 = n1;
    /* Allocate the virtual image {img}: */
    msm_image_t *img = msm_image_alloc(1, nv0, nv1, TRUE);
    /* Fill the image with "no value": */
    float_image_fill(img->fim, (float)(-INF));
    /* Paint candidate on image: */
    msm_image_cand_paint(cd, &v, FALSE, img, 0, 0);
    /* Normalize image and write it out: */
    msm_image_normalize_and_write_as_pgm(img, name, tag);
    /* Reclaim image space: */
    float_image_free(img->fim); free(img);
  }
    
void msm_image_cand_vec_write_named
  ( msm_cand_vec_t *cdv,
    msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1,
    bool_t scale, 
    char *name, 
    char *tag
  )
  { /* Get sequence sizes {n0,n1}: */
    int32_t n0 = seq0->size;
    int32_t n1 = seq1->size;
    /* Choose the unscaled image size {nv0,nv1 = n0,n1}: */
    int32_t nv0 = n0, nv1 = n1;
    /* Allocate the virtual image {img}: */
    msm_image_t *img = msm_image_alloc(1, nv0, nv1, scale);
    /* Paint candidates: */
    msm_image_cand_vec_paint(cdv, seq0, seq1, img, 0, 0);
    /* Normalize image and write it out: */
    msm_image_normalize_and_write_as_pgm(img, name, tag);
    /* Reclaim image space: */
    float_image_free(img->fim); free(img);
  }

void msm_image_dyn_tableau_write_named
  ( int32_t n0,                /* Length of first sequence. */
    int32_t n1,                /* Length of second sequence. */
    msm_dyn_tableau_t *tb, /* Tableau to plot. */
    msm_rung_t gopt,       /* Optimal end-rung or {msm_rung_none} */
    bool_t scale,          /* TRUE to adjust the scale to avoid large images. */   
    bool_t colorize,       /* TRUE to map scores to RGB colors instead of gray values. */   
    char *name,            /* Prefix for file name. */
    char *tag              /* Tag for file name. */
  )
  { /* Find out the X and Y ranges covered by {tb}: */
    int32_t i0Min, i0Max;  /* Min and max X in {tb}. */
    int32_t i1Min, i1Max;  /* Min and max Y in {tb}. */
    msm_dyn_tableau_get_i0_i1_ranges(tb, &i0Min, &i0Max, &i1Min, &i1Max);
    fprintf(stderr, "plotting tableau\n");
    fprintf(stderr, "  coordinate ranges");
    fprintf(stderr, "[%4d .. %4d]×[%4d .. %4d]\n", i0Min, i0Max, i1Min, i1Max);
    /* Choose the unscaled image size {nv0,nv1 = n0,n1}, for UNfolded plot: */
    int32_t nv0 = i0Max - i0Min + 1;
    int32_t nv1 = i1Max - i1Min + 1;
    /* Allocate a virtual greyscale image {gim}: */
    msm_image_t *gim = msm_image_alloc(1, nv0, nv1, scale);
    /* Initialize the image with {-INF}: */
    float_image_fill(gim->fim, -INF);
    /* Paint tableau {tb} into float image {gim}, note maximum abs value {maxf}: */
    msm_image_dyn_tableau_scores_paint(n0, n1, tb, gim, i0Min, i1Min);
    /* Find the pixel value range {[-maxf _ +maxf]}: */
    double lov, hiv;
    msm_image_compute_min_max_range(gim, 0, &lov, &hiv);
    double maxf = (fabs(lov) > fabs(hiv) ? fabs(lov) : fabs(hiv));
    frgb_t opath_clr; /* Color for the optimum path. */
    msm_image_t *cim = NULL;
    double vmin,vmax;
    if (colorize)
      { /* Convert {gim} to a color pnm image. Maps 0 to gray, {-maxf} to blue, {+maxf} to red. */
        cim = msm_image_colorize(gim, maxf);
        assert(cim->fim->sz[0] == 3); /* Hope so. */
        /* Paint the indices that are multiples of {n0} or {n1} in green: */
        opath_clr = (frgb_t){{ 0.000f, 0.800f, 0.000f }}; /* Slightly dark reddish yellow */
        vmin = 0;  vmax = 1;
      }
    else
      { cim = gim;
        assert(cim->fim->sz[0] == 1); /* Hope so. */
        opath_clr = (frgb_t){{ 0.000f, 0.000f, 0.000f }}; /* Black. */
        vmin = lov; vmax = hiv;
      }

    /* Paint the optimum path: */
    msm_image_dyn_tableau_pairing_paint(n0, n1, tb, gopt, opath_clr.c, cim, i0Min, i1Min);
    /* Write {cim} to PPM image. */
    msm_image_write_as_pnm(cim, vmin, vmax, name, tag);
    /* Reclaim the image's storage: */
    if(cim != gim) { float_image_free(cim->fim); free(cim); }
    float_image_free(gim->fim); free(gim);
  }
