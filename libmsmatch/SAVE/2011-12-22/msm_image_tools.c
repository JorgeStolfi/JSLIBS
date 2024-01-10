/* See {msm_image_tools.h} */
/* Last edited on 2008-01-29 17:21:29 by hcgl */

#define msm_image_tools_C_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>

#include <vec.h>
#include <float_image.h>
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
  ( msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp, 
    msm_rung_score_proc_t *score, 
    msm_pairing_t *pr, 
    char *name, 
    char *tag
  )
  { int nx = xp->npos;
    int ny = yp->npos;
    msm_image_t *gim = msm_image_alloc(1, nx, ny);
    fprintf(stderr, "  image size = %d × %d  scale factor = %d\n", nx, ny, gim->scale);
    /* Paint the sample-to-sample distances: */
    msm_image_seq_seq_score_paint(xp, yp, score, gim, 0, 0);
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
      { float clr[3] = { 0.000, 0.700, 0.000 }; /* Dark green. */
        msm_image_pairing_paint(pr, nx, ny, clr, FALSE, cim, 0, 0);
      }
    /* Write image: */
    msm_image_write_as_pnm(cim, 0, 1, name, tag);
    /* Reclaim image space: */
    float_image_free(gim->fim); free(gim);
    float_image_free(cim->fim); free(cim);
  }

void msm_image_cand_write_named(msm_cand_t *cd, float v, char *name, char *tag)
  { /* Get sequence sizes {nx,ny}: */
    int nx = cd->seq[0].npos;
    int ny = cd->seq[1].npos;
    /* Choose the unscaled image size {nix,niy = nx,ny}, for folded plot: */
    int nix = nx, niy = ny;
    /* Allocate the virtual image {img}: */
    msm_image_t *img = msm_image_alloc(1, nix, niy);
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
    msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp,
    char *name, 
    char *tag
  )
  { /* Get sequence sizes {nx,ny}: */
    int nx = xp->npos;
    int ny = yp->npos;
    /* Choose the unscaled image size {nix,niy = nx,ny}, for folded plot: */
    int nix = nx, niy = ny;
    /* Allocate the virtual image {img}: */
    msm_image_t *img = msm_image_alloc(1, nix, niy);
    /* Paint candidates: */
    msm_image_cand_vec_paint(cdv, xp, yp, img, 0, 0);
    /* Normalize image and write it out: */
    msm_image_normalize_and_write_as_pgm(img, name, tag);
    /* Reclaim image space: */
    float_image_free(img->fim); free(img);
  }

void msm_image_dyn_tableau_write_named
  ( int nx,                /* Length of first sequence. */
    int ny,                /* Length of second sequence. */
    msm_dyn_tableau_t *tb, /* Tableau to plot. */
    msm_rung_t gopt,       /* Optimal end-rung or {msm_rung_none} */
    char *name,            /* Prefix for file name. */
    char *tag              /* Tag for file name. */
  )
  { /* Find out the X and Y ranges covered by {tb}: */
    int xMin, xMax;  /* Min and max X in {tb}. */
    int yMin, yMax;  /* Min and max Y in {tb}. */
    msm_dyn_tableau_get_X_Y_ranges(tb, &xMin, &xMax, &yMin, &yMax);
    fprintf(stderr, "plotting tableau\n");
    fprintf(stderr, "  original coordinate ranges");
    fprintf(stderr, "[%4d .. %4d]×[%4d .. %4d]\n", xMin, xMax, yMin, yMax);
    /* Round {xMin,yMin} down to multiples of the sequence lengths: */
    xMin = ifloor(xMin, nx);
    yMin = ifloor(yMin, ny);
    /* Round {xMax,yMax} up to multiples of the sequence lengths, minus 1: */
    xMax = iceil(xMax + 1, nx) - 1;
    yMax = iceil(yMax + 1, ny) - 1;
    fprintf(stderr, "  adjusted coordinate ranges");
    fprintf(stderr, "[%4d .. %4d]×[%4d .. %4d]\n", xMin, xMax, yMin, yMax);
    /* Choose the unscaled image size {nix,niy = nx,ny}, for UNfolded plot: */
    int NVX = xMax - xMin + 1;
    int NVY = yMax - yMin + 1;
    /* Allocate a virtual greyscale image {gim}: */
    msm_image_t *gim = msm_image_alloc(1, NVX, NVY);
    /* Initialize the image with {-INF}: */
    float_image_fill(gim->fim, -INF);
    /* Paint tableau {tb} into float image {gim}, note maximum abs value {maxf}: */
    msm_image_dyn_tableau_scores_paint(nx, ny, tb, gim, xMin, yMin);
    /* Find the pixel value range {[-maxf _ +maxf]}: */
    double lov, hiv;
    msm_image_compute_min_max_range(gim, 0, &lov, &hiv);
    double maxf = (fabs(lov) > fabs(hiv) ? fabs(lov) : fabs(hiv));
    /* Convert {gim} to a color pnm image. Maps 0 to gray, {-maxf} to blue, {+maxf} to red. */
    msm_image_t *cim = msm_image_colorize(gim, maxf);
    int NC = cim->fim->sz[0]; /* Channels in color image. */
    assert(NC == 3); /* Hope so. */
    /* Paint the indices that are multiples of {nx} or {ny} in green: */
    float coord_clr[3] = { 0.100, 0.700, 0.000 }; /* Slightly yellowish dark green */
    msm_image_seq_periods_paint(nx, ny, coord_clr, cim, xMin, yMin);
    /* Paint the optimum path in yellow: */
    float opath_clr[3] = { 0.900, 0.800, 0.000 }; /* Slightly dark reddish yellow */
    msm_image_dyn_tableau_pairing_paint(nx, ny, tb, gopt, opath_clr, cim, xMin, yMin);
    /* Write {cim} to PPM image. */
    msm_image_write_as_pnm(cim, 0, 1, name, tag);
    /* Reclaim the image's storage: */
    float_image_free(gim->fim); free(gim);
    float_image_free(cim->fim); free(cim);
  }
