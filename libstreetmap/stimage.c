/* See stimage.h */
/* Last edited on 2017-07-01 00:03:58 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <values.h>
#include <assert.h>

#include <fget.h>
#include <nget.h>
#include <affirm.h>
#include <uint16_image.h>
#include <uint16_image_read_pnm.h>
#include <uint16_image_write_pnm.h>

#include <stmap.h>
#include <stimage.h>

/* Index of sample given X and Y indices (assumed to be within bounds): */
#define INDEX(IMG,IX,IY) \
  ((IX)+(IMG)->nx*(IY))

/* Index of sample given X and Y indices (-1 if outside bounds): */
#define CHKINDEX(IMG,IX,IY) \
  ( ((IX) >= 0) && ((IX) < (IMG)->nx) && ((IY) >= 0) && ((IY) < (IMG)->ny) ? \
    ((IX)+(IMG)->nx*(IY)) : -1 \
  )

DataImage *st_image_new
  ( int mmXLo,  /* Lowest X sample coordinate, in mm. */
    int nx,     /* Number of samples in X direction. */
    int mmYLo,  /* Lowest Y sample coordinate, in mm. */
    int ny,     /* Number of samples in Y direction. */
    int mmStep  /* Step between samples, in mm */
  )
  {
    DataImage *img = (DataImage *)malloc(sizeof(DataImage));
    affirm(img != NULL, "out of mem");
    img->nx = nx;
    img->ny = ny;
    img->mmStep = mmStep;
    img->mmXLo = mmXLo; 
    img->mmYLo = mmYLo;
    /* Compute image size: */
    int ns = nx*ny;
    img->d = (double*)malloc(ns*sizeof(double));
    affirm(img->d != NULL, "out of mem");
    int i;
    for (i = 0; i < ns; i++) { img->d[i] = 0; }
    return img;
  }

DataImage *st_image_for_rect(Interval rx, Interval ry, double step)
  { /* Compute number of steps: */
    int nx = (int)((rx.hi - rx.lo)/step + 0.5);
    int ny = (int)((ry.hi - ry.lo)/step + 0.5);
    /* Adjust step and bounds to integer number of millimeters: */
    int mmStep = (int)ceil(step*1000.0);
    int mmXLo = (int)((rx.lo + 0.5*step)*1000 + 0.5); 
    int mmYLo = (int)((ry.lo + 0.5*step)*1000 + 0.5);
    /* Allocate image: */
    return st_image_new(mmXLo, nx, mmYLo, ny, mmStep);
  }

int st_image_locate(DataImage *img, double x, double y)
  { /* Nearest sample: */
    int ix = (int)floor((x*1000.0 - img->mmXLo)/img->mmStep + 0.5);
    int iy = (int)floor((y*1000.0 - img->mmYLo)/img->mmStep + 0.5);
    return CHKINDEX(img,ix,iy);
  }

double st_image_get(DataImage *img, double x, double y)
  { int k = st_image_locate(img, x, y);
    return (k >= 0 ? img->d[k] : 0.0);
  }

void st_image_set(DataImage *img, double x, double y, double val)
  { int k = st_image_locate(img, x, y);
    if (k >= 0) { img->d[k] = val; } 
  }

void st_image_increment(DataImage *img, double x, double y, double val)
  { int k = st_image_locate(img, x, y);
    if (k >= 0) { img->d[k] += val; } 
  }

double st_image_interpolate(DataImage *img, double x, double y)
  { /* Displacements in mm from lower corner of area: */
    double st = (double)(img->mmStep);
    int dx = (int)floor(x*1000.0 - img->mmXLo + 0.5*st);
    int dy = (int)floor(y*1000.0 - img->mmYLo + 0.5*st);
    
    /* Indices of nearest sample: */
    int ix = dx/img->mmStep;
    int iy = dy/img->mmStep;
    
    /* Fraction of step towards next sample: */
    double fx = (double)(dx % img->mmStep)/st;
    double fy = (double)(dy % img->mmStep)/st;

    /* Add values four adjacent pixels, with bilinear weights: */
    int sx, sy;
    double sum = 0.0;
    for (sy = 0; sy <= 1; sy++)
      { double ty = (sy == 0 ? 1.0-fy : fy);
        for (sx = 0; sx <= 1; sx++)
          { double tx = (sx == 0 ? 1.0-fx : fx);
            int k = CHKINDEX(img,ix+sx,iy+sy); 
            if (k >= 0) { sum = tx*ty*img->d[k]; }
          }
      }
    return sum;
  }

DataImage *st_image_blur(DataImage *img, double rad)
  { /* Compute a unidimensional Gaussian with variance {rad^2/2}: */
    int nr = (int)ceil(3.0*rad*1000/img->mmStep);
    double *mask = (double *)malloc((nr+1)*sizeof(double));
    double S = rad*rad/2.0;
    int i;
    for (i = 0; i <= nr; i++) 
      { double x = ((double)i*img->mmStep)/1000.0;
        mask[i] = exp(-x*x/S);
      }
    /* Now convolve bidimensional gaussian (normalized) with given image: */
    DataImage *res = st_image_new
      ( img->mmXLo, img->nx, img->mmYLo, img->ny, img->mmStep );
    int ix, iy;
    for (iy = 0; iy < img->ny; iy++)
      { for (ix = 0; ix < img->nx; ix++)
          { /* Compute weighted average of samples around {[ix,iy]}: */
            int dx, dy;
            double sw = 0, swd = 0;
            for (dy = -nr; dy <= +nr; dy++)
              { int jy = iy+dy;
                if ((jy >= 0) && (jy < img->ny))
                  { for (dx = -nr; dx <= +nr; dx++)
                      { int jx = ix+dx;
                        if ((jx >= 0) && (jx < img->nx))
                          { int kimg = INDEX(img,jx,jy);
                            double mx = mask[abs(dx)];
                            double my = mask[abs(dy)];
                            swd += img->d[kimg]*mx*my;
                            sw += mx*my;
                          }
                      }
                  }
              }
            if (sw > 0.0) { swd /= sw; }
            int kres = INDEX(res,ix,iy);
            res->d[kres] = swd;            
          }
      }
    return res;
  }

Interval st_image_sample_range(DataImage *img)
  { double lo = +INFINITY;
    double hi = -INFINITY;
    int ns = img->nx*img->ny, k;
    for (k = 0; k < ns; k++)
      { double dk = img->d[k];
        if (dk < lo) { lo = dk; }
        if (dk > hi) { hi = dk; }
      }
    return (Interval){lo, hi};
  }

void st_image_write_as_pgm
  ( FILE *wr, 
    DataImage *img,
    double bval,
    double wval
  )
  { /* Allocate image with one pixel per sample; */
    uint16_image_t *qim = uint16_image_new(img->nx, img->ny, 1);
    qim->maxval = 65535; /* Some versions of PGM/PBM can't handle more. */
    double mv = (double)qim->maxval;
    /* Convert data samples to gray tone values: */
    int ix, iy;
    for (iy = 0; iy < img->ny; iy++)
      { uint16_t *smp = qim->smp[img->ny-1-iy]; 
        for (ix = 0; ix < img->nx; ix++)
          { int k = INDEX(img,ix,iy);
            double z = (img->d[k] - bval)/(wval - bval);
            double zs = floor(z*mv + 0.5);
            if (zs < 0.0) { zs = 0.0; }
            if (zs > mv) { zs = mv; }
            smp[ix] = (uint16_t)((int)zs);
          }
      }
    bool_t forceplain = FALSE;
    bool_t verbose = FALSE;
    uint16_image_write_pnm_file(wr, qim, forceplain, verbose);
    fprintf(wr, "xLo(mm) = %d\n", img->mmXLo);
    fprintf(wr, "yLo(mm) = %d\n", img->mmYLo);
    fprintf(wr, "step(mm) = %d\n", img->mmStep);
    fprintf(wr, "bval = %24.16e\n", bval);
    fprintf(wr, "wval = %24.16e\n", wval);
    uint16_image_free(qim);
  }

DataImage *st_image_read_from_pgm(FILE *rd)
  { /* Read pgm image; */
    uint16_image_t *qim = uint16_image_read_pnm_file(rd);
    assert(qim->chns == 1);
    int mmXLo = nget_int(rd, "xLo(mm)"); fget_eol(rd);
    int mmYLo = nget_int(rd, "yLo(mm)"); fget_eol(rd);
    int mmStep = nget_int(rd, "step(mm)"); fget_eol(rd);
    double bval = nget_double(rd, "bval"); fget_eol(rd);
    double wval = nget_double(rd, "wval"); fget_eol(rd);
    DataImage *img = st_image_new(mmXLo, qim->cols, mmYLo, qim->rows, mmStep);
    double mv = (double)qim->maxval;
    /* Convert gray pixels to sample values: */
    int ix, iy;
    for (iy = 0; iy < img->ny; iy++)
      { uint16_t *smp = qim->smp[img->ny-1-iy]; 
        for (ix = 0; ix < img->nx; ix++)
          { int k = INDEX(img,ix,iy);
            double zs = (double)smp[ix];
            double z = bval + (wval - bval)*(zs/mv);
            img->d[k] = z;
          }
      }
    uint16_image_free(qim);
    return img;
  }
