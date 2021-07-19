/* See {msm_image.h} */
/* Last edited on 2021-07-18 00:42:32 by jstolfi */

#define msm_image_C_COPYRIGHT \
  "Copyright © 2005 by the State University of Campinas (UNICAMP)"

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include <affirm.h>
#include <float_image.h>
#include <float_image_to_uint16_image.h>
#include <jsmath.h>
#include <uint16_image.h>
#include <uint16_image_write_pnm.h>

#include <msm_basic.h>

#include <msm_image.h>

void msm_image_reduce_coords(int *x, int *y, int NVX, int NVY);
  /* Reduces the (unscaled) coordinate {*x} to the range {0..NVX-1}
    modulo {NVX}, and the coordinate {*y} to the range {0..NVY-1} modulo
    {NVY}. */

float_image_t *float_image_colorize(float_image_t *gfim, double maxf);
/* Returns a color image with same size as {gfim}, where negative samples 
  are mapped to bluish tones, positive samples to reddish tones, 0 to white,
  {±INF} to black, and {NAN} to grey. */

msm_image_t *msm_image_alloc(int NC, int NVX, int NVY, bool_t scale)
  { /* Choose the image downscaling factor {scale}: */
    int nbig = (NVX > NVY ? NVX : NVY); /* Length of largest sequence. */
    int nmax = 1024; /* Maximum alowed image size. */
    int ppp = (scale ? (nbig + nmax - 1)/nmax : 1); /* Positions per pixel. */
    if (ppp == 0) { ppp = 1; }
    /* Compute actual image size {NPX,NPY}: */
    int NPX = (NVX - 1)/ppp + 1;
    int NPY = (NVY - 1)/ppp + 1;
    fprintf(stderr, "  virtual image %4d × %4d", NVX, NVY);
    fprintf(stderr, "  actual image %4d × %4d", NPX, NPY);
    fprintf(stderr, "  scale =  %4d positions per pixel\n", ppp);
    /* Allocate image: */
    float_image_t *fim = float_image_new(NC, NPX, NPY);
    /* Return results: */
    msm_image_t *img = (msm_image_t *)notnull(malloc(sizeof(msm_image_t)), "out of mem");
    (*img) = (msm_image_t) 
      { .NV = { NVX, NVY },
        .scale = ppp,
        .fim = fim
      };
    return img;
  }

void msm_image_compute_avg_dev_range(msm_image_t *img, int c, double *lo, double *hi, double ns)
  { double avg, dev;
    float_image_compute_sample_avg_dev(img->fim, c, &avg, &dev);
    /* Choose black and white values: */
    double vlo = avg - ns*dev;
    double vhi = avg + ns*dev;
    assert(vlo <= vhi);
    if (vlo == vhi)
      { vlo = nextafter(vlo, -INF);
        vhi = nextafter(vhi, +INF);
      }
    /* Return them: */
    (*lo) = vlo; (*hi) = vhi;
  }

void msm_image_compute_min_max_range(msm_image_t *img, int c, double *lo, double *hi)
  { /* Find range {vMin,vMax} of samples: */
    float vMin = +INF, vMax = -INF;
    float_image_update_sample_range(img->fim, c, &vMin, &vMax);
    /* fprintf(stderr, "  float_image score range = [ %10.6f _ %10.6f ]\n", vMin, vMax); */
    /* Provide default range {[0_0]} if there are no finite pixels: */
    if (vMin > vMax) { vMin = vMax = 0; }
    /* Make sure that {vMin < vMax}: */
    if (vMin == vMax) 
      { /* Fudge the bounds to make them different: */
        float eps = (float)(1e-6 * fabs(vMin));
        if (eps < 1e-30f) { eps = 1e-30f; }
        vMin -= eps; vMax += eps;
      }
    /* Return them: */
    (*lo) = vMin; (*hi) = vMax;
    fprintf(stderr, "  msm_image   score range = [ %10.6f _ %10.6f ]\n", (*lo), (*hi));
  }

float_image_t *float_image_colorize(float_image_t *gfim, double maxf)
  { /* Check the channel count and get image dimensions: */
    demand(gfim->sz[0] == 1, "input image must be grayscale");
    int NC = 3; /* Number of channels in output image. */
    int NX = (int)gfim->sz[1]; /* Number of columns. */
    int NY = (int)gfim->sz[2]; /* Number of rows. */
    float_image_t *cfim = float_image_new(NC, NX, NY);
    /* Limiting colors: */
    float zer[3] = { 1.000f, 1.000f, 1.000f }; /* Lum = 1.000; for {v = 0}. */
    float pos[3] = { 1.000f, 0.333f, 0.000f }; /* Lum = 0.500; for {v=+maxf}. */
    float neg[3] = { 0.000f, 0.667f, 1.000f }; /* Lum = 0.500; for {v=-maxf}. */
    /* Convert monochrome pixels from {gfim} to color pixels of {cfim}: */
    int x, y, c;
    for(y = 0; y < NY; y++)
      { for(x = 0; x < NX; x++)
          { /* Get pixel {v} from monochrome image: */
            float v = float_image_get_sample(gfim, 0, x, y);
            /* Compute the corresponding color {clr}: */
            float clr[NC];
            if (isnan(v))
              { /* Assume pixel as undefined, colorize it with grey: */
                for (c = 0; c < NC; c++) { clr[c] = 0.5; }
              }
            else if (fabsf(v) == fabsf(INF))
              { /* Colorize it with black: */
                for (c = 0; c < NC; c++) { clr[c] = 0.0; }
              }
            else
              { /* Choose the limiting color {lim} appropriate to the sign of {v}: */
                float *lim = (v >= 0 ? pos : neg); /* Limiting color for sign of {v}. */
                /* Compute intensity {s} of {v} relative to {±maxf}, and its complement {t}: */
                double s = fabs(v)/maxf;
                assert(s >= 0); 
                if (s > 1) { s = 1; }
                double t = 1 - s;
                /* Interpolate with {s} between {zer} and {lim} to get the color {clr}: */
                for (c = 0; c < NC; c++) { clr[c] = (float)(t*zer[c] + s*lim[c]); }
              }
            /* Store {clr} in color image: */
            for (c = 0; c < NC; c++) 
              { float_image_set_sample(cfim, c, x, y, clr[c]); }
          }
      }
    return cfim;
  }

msm_image_t *msm_image_colorize(msm_image_t *gim, double maxf)
  { msm_image_t *cim = (msm_image_t *)notnull(malloc(sizeof(msm_image_t)), "out of mem");
    cim->NV[0] = gim->NV[0];
    cim->NV[1] = gim->NV[1]; 
    cim->scale = gim->scale;
    cim->fim = float_image_colorize(gim->fim, maxf);
    return cim;
  }

void msm_image_write_as_pnm
  ( msm_image_t *img, 
    double minv, 
    double maxv, 
    char *name,
    char *tag
  )
  { /* Chck the channel count: */
    int NC = (int)img->fim->sz[0];
    demand((NC == 1) || (NC == 3), "bad number of channels"); 
    /* Prepare the per-channel {lo,hi} parameters from single {minv,maxv}: */
    double lov[NC];
    double hiv[NC];
    int c; 
    for (c = 0; c < NC; c++) { lov[c] = minv; hiv[c] = maxv; }
    /* Convert the float image {img} to an integer image {pnm}: */
    bool_t yup = TRUE, verbose = TRUE;
    bool_t isMask = TRUE; /* Assume 0 and 1 are important in encoding/decoding. */
    uint16_image_t *pnm = float_image_to_uint16_image(img->fim, isMask, NC, lov, hiv, NULL, 255, yup, verbose);
    /* Construct the file name {fileName}: */
    char *ext = (NC == 1 ? ".pgm" : ".ppm");
    /* Open the output file {wr} and write {pnm} there: */
    FILE *wr = msm_open_write(name, tag, ext, TRUE);
    bool_t forceplain = FALSE;
    uint16_image_write_pnm_file(wr, pnm, forceplain, verbose);
    /* Cleanup: */
    fclose(wr);
    uint16_image_free(pnm);
  }

void msm_image_reduce_coords(int *x, int *y, int NVX, int NVY)
  { (*x) = (int)imod((*x), NVX); 
    (*y) = (int)imod((*y), NVY);
  }

void msm_image_set(msm_image_t *img, int ic, int x, int y, float v)
  { /* fprintf(stderr, "    setting %10.5f to virtual pixel ( %1d %4d %4d )\n", v, ic, x, y); */
    msm_image_reduce_coords(&x, &y, img->NV[0], img->NV[1]);
    int ipx = x/img->scale;
    int ipy = y/img->scale;
    float *p = float_image_get_sample_address(img->fim, ic, ipx, ipy);
    (*p) = v;
    /* fprintf(stderr, "    actual pixel ( %1d %4d %4d ) became %10.5f\n", ic, ipx, ipy, (*p)); */
  }

void msm_image_add(msm_image_t *img, int ic, int x, int y, float v)
  { /* fprintf(stderr, "    adding %10.5f to virtual pixel ( %4d %4d )\n", v, ic, x, y); */
    msm_image_reduce_coords(&x, &y, img->NV[0], img->NV[1]);
    int ipx = x/img->scale;
    int ipy = y/img->scale;
    float *p = float_image_get_sample_address(img->fim, ic, ipx, ipy);
    if (isfinite((*p)))
      { (*p) += v; }
    else
      { (*p) = v; }
    /* fprintf(stderr, "    actual pixel ( %4d %4d ) became %10.5f\n", ic, ipx, ipy, (*p)); */
  }

void msm_image_max(msm_image_t *img, int ic, int x, int y, float v)
  { msm_image_reduce_coords(&x, &y, img->NV[0], img->NV[1]);
    int ipx = x/img->scale;
    int ipy = y/img->scale;
    float *p = float_image_get_sample_address(img->fim, ic, ipx, ipy);
    if (v > (*p)) { (*p) = v; }
  }

void msm_image_normalize_and_write_as_pgm(msm_image_t *img, char *name, char *tag)
  { /* Chck the channel count: */
    int NC = (int)img->fim->sz[0];
    demand(NC == 1, "image must be grayscale");
    /* Find actual range of image (excluding infinite pixels): */
    double lov, hiv;
    msm_image_compute_min_max_range(img, 0, &lov, &hiv);
    /* Adjust {lov} so that lowest score maps to gray: */
    lov = lov - (hiv - lov);
    /* Write out as PGM: background white, candidates gray to black: */
    msm_image_write_as_pnm(img, hiv, lov, name, tag);
  }
