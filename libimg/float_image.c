/* See float_image.h */
/* Last edited on 2023-03-19 08:41:31 by stolfi */ 

#define _GNU_SOURCE
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include <sample_conv.h>
#include <frgb.h>
#include <frgb_ops.h>
#include <float_image_color.h>
#include <ix.h>
#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <affirm.h>
#include <bool.h>

#include <float_image.h>

/* IMPLEMENTATIONS */

float float_image_get_sample(float_image_t *A, int32_t c, int32_t x, int32_t y)
  { assert((c >= 0) && (c < A->sz[0]));
    assert((x >= 0) && (x < A->sz[1]));
    assert((y >= 0) && (y < A->sz[2]));
    ix_index_t ix[3] = {c, x, y};
    ix_pos_t p = ix_position(3, ix, A->bp, A->st);
    return A->sample[p];
  }

float *float_image_get_sample_address(float_image_t *A, int32_t c, int32_t x, int32_t y)
  { assert((c >= 0) && (c < A->sz[0]));
    assert((x >= 0) && (x < A->sz[1]));
    assert((y >= 0) && (y < A->sz[2]));
    ix_index_t ix[3] = {c, x, y};
    ix_pos_t p = ix_position(3, ix, A->bp, A->st);
    return &(A->sample[p]);
  }

void float_image_set_sample(float_image_t *A, int32_t c, int32_t x, int32_t y, float v)
  { assert((c >= 0) && (c < A->sz[0]));
    assert((x >= 0) && (x < A->sz[1]));
    assert((y >= 0) && (y < A->sz[2]));
    ix_index_t ix[3] = {c, x, y};
    ix_pos_t p = ix_position(3, ix, A->bp, A->st);
    A->sample[p] = v;
  }

void float_image_get_pixel(float_image_t *A, int32_t x, int32_t y, float v[])
  { float *sp = float_image_get_sample_address(A, 0, x, y);
    int32_t NC = (int32_t)A->sz[0];
    for (int32_t c = 0; c < NC; c++) { v[c] = (*sp); sp += A->st[0]; }
  }

void float_image_set_pixel(float_image_t *A, int32_t x, int32_t y, float v[])
  { float *sp = float_image_get_sample_address(A, 0, x, y);
    int32_t NC = (int32_t)A->sz[0];
    for (int32_t c = 0; c < NC; c++) { (*sp) = v[c]; sp += A->st[0]; }
  }

void float_image_fill_pixel(float_image_t *A, int32_t x, int32_t y, float v)
  { float *sp = float_image_get_sample_address(A, 0, x, y);
    int32_t NC = (int32_t)A->sz[0];
    for (int32_t c = 0; c < NC; c++) { (*sp) = v; sp += A->st[0]; }
  }

void float_image_get_sample_row
  ( float_image_t *A, 
    int32_t c, 
    int32_t xmin, 
    int32_t xmax, 
    int32_t y, 
    float v[] 
  )
  { if (xmin <= xmax)
      { int32_t NX = (int32_t)A->sz[1];  /* Count of columns of {A}. */
        int32_t nrx = xmax - xmin + 1; /* Count of samples to process. */
        /* Get pointer {sp} to next sample to fetch: */
        int32_t xmin_get = (xmin < 0 ? 0 : xmin); /* First column to actually fetch from {A}. */
        float *sp = ((xmin < NX) && (xmax >= 0) ? float_image_get_sample_address(A, c, xmin_get, y) : NULL);
        for (int32_t ix = 0; ix < nrx; ix++) 
          { if ((ix >= 0) && (ix < NX))
              { v[ix] = (*sp); sp += A->st[1]; }
            else
              { v[ix] = NAN; }
          }
      }
  }

void float_image_set_sample_row
  ( float_image_t *A, 
    int32_t c, 
    int32_t xmin, 
    int32_t xmax, 
    int32_t y, 
    float v[] 
  )
  { if (xmin <= xmax)
      { int32_t nrx = xmax - xmin + 1; /* Count of samples to process. */
        /* Get pointer {sp} to next sample to set: */
        float *sp = float_image_get_sample_address(A, c, xmin, y);
        for (int32_t ix = 0; ix < nrx; ix++) { (*sp) = v[ix]; sp += A->st[1]; }
      }
  }

void float_image_fill_sample_row 
  ( float_image_t *A,  
    int32_t c,  
    int32_t xmin,  
    int32_t xmax,  
    int32_t y,  
    float v 
  )
  { if (xmin <= xmax)
      { int32_t nrx = xmax - xmin + 1; /* Count of samples to process. */
        /* Get pointer {sp} to next sample to set: */
        float *sp = float_image_get_sample_address(A, c, xmin, y);
        for (int32_t ix = 0; ix < nrx; ix++) { (*sp) = v; sp += A->st[1]; }
      }
  }

void float_image_get_pixel_row(float_image_t *A, int32_t xmin, int32_t xmax, int32_t y, float v[])
  { if (xmin <= xmax)
      { int32_t NC = (int32_t)A->sz[0];
        int32_t NX = (int32_t)A->sz[1];  /* Count of columns of {A}. */
        /* Get pointer {pxp} to first sample of next pixel to set: */
        int32_t xmin_get = (xmin < 0 ? 0 : xmin); /* First column to actually fetch from {A}. */
        float *pxp = ((xmin < NX) && (xmax >= 0) ? float_image_get_sample_address(A, 0, xmin_get, y) : NULL);
        int32_t nrx = xmax - xmin + 1; /* Count of pixels to process. */
        int32_t k = 0; /* Index of sample in {v}. */
        for (int32_t ix = 0; ix < nrx; ix++)
          { float *sp = pxp; /* Pointer to next unfetched sample in {A}. */
            for (int32_t c = 0; c < NC; c++) 
              { if ((ix >= 0) && (ix < NX)) 
                  { v[k] = (*sp); sp += A->st[0]; }
                else
                  { v[ix] = NAN; }
                k++; 
              }
            pxp += A->st[1];
          }
      }
  }

void float_image_set_pixel_row(float_image_t *A, int32_t xmin, int32_t xmax, int32_t y, float v[])
  { if (xmin <= xmax)
      { int32_t NC = (int32_t)A->sz[0];
        /* Get pointer {pxp} to first sample of next pixel to set: */
        float *pxp = float_image_get_sample_address(A, 0, xmin, y);
        int32_t nrx = xmax - xmin + 1; /* Count of pixels to process. */
        int32_t k = 0; /* Index of sample in {v}. */
        for (int32_t ix = 0; ix < nrx; ix++)
          { float *sp = pxp; /* Pointer to next unset sample in {A}. */
            for (int32_t c = 0; c < NC; c++) 
              { (*sp) = v[k]; sp += A->st[0]; k++; }
            pxp += A->st[1];
          }
      }
  }

void float_image_fill_pixel_row(float_image_t *A, int32_t xmin, int32_t xmax, int32_t y, float v)
  { if (xmin <= xmax)
      { int32_t NC = (int32_t)A->sz[0];
        /* Get pointer {pxp} to first sample of next pixel to set: */
        float *pxp = float_image_get_sample_address(A, 0, xmin, y);
        int32_t nrx = xmax - xmin + 1; /* Count of pixels to process. */
        for (int32_t ix = 0; ix < nrx; ix++)
          { float *sp = pxp; /* Pointer to to first sample of next unset pixel. */
            for (int32_t c = 0; c < NC; c++) 
              { (*sp) = v; sp += A->st[0]; }
            pxp += A->st[1];
          }
      }
  }

void float_image_fill_channel(float_image_t *A, int32_t c, float v)
  { int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1]; 
    int32_t NY = (int32_t)A->sz[2];
    demand((c >= 0) && (c < NC), "invalid channel");
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float_image_set_sample(A, c, x, y, v); }
      }
  }

void float_image_set_channel(float_image_t *A, int32_t cA, float_image_t *V, int32_t cV)
  { int32_t NCV = (int32_t)V->sz[0];
    demand((cV >= 0) && (cV < NCV), "invalid {V} channel");
    int32_t NCA = (int32_t)A->sz[0];
    demand((cA >= 0) && (cA < NCA), "invalid {A} channel");
    int32_t NX = (int32_t)V->sz[1]; 
    int32_t NY = (int32_t)V->sz[2];
    float_image_check_size(A, -1, NX, NY);
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float v = float_image_get_sample(V, cV, x, y);
            float_image_set_sample(A, cA, x, y, v);
          }
      }
  }

void float_image_get_window_samples
  ( float_image_t *A, 
    int32_t c, 
    int32_t x, 
    int32_t y, 
    int32_t nwx, 
    int32_t nwy, 
    bool_t rep, 
    float v[] 
  )
  {
    int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1]; 
    int32_t NY = (int32_t)A->sz[2];
    demand((c >=0) && (c < NC), "invalid channel");
    
    demand(nwx % 2 == 1, "{nwx} must be odd");
    demand(nwy % 2 == 1, "{nwy} must be odd");
    int32_t rwx = nwx/2; /* Half-width in X. */
    int32_t rwy = nwy/2; /* Half-width in Y. */
    
    /* Scan the window pixels: */
    int32_t kw = 0;
    for (int32_t iy = -rwy; iy <= rwy; iy++)
      { for (int32_t ix = -rwx; ix <= rwx; ix++)
          { /* Compute the sample's indices {xp,yp} in input img, or -1 if non-existant: */
            int32_t xp = x + ix;
            int32_t yp = y + iy;
            if (rep)
              { /* Map invalid {xp,yp} to nearest pixel in domain: */
                if (yp < 0) { yp = 0; } else if (yp >= NY) { yp = NY - 1; }
                if (xp < 0) { xp = 0; } else if (xp >= NX) { xp = NX - 1; }
              }
            /* Set invalid {xp,yp} to -1: */
            if ((xp < 0) || (xp >= NX) || (yp < 0) || (yp >= NY)) { xp = yp = -1; }
            /* Get and set the value: */  
            v[kw] = ( yp < 0 ? NAN : float_image_get_sample(A, c, xp, yp) );
            kw++;
          }
      }
  }

void float_image_get_window_pixels
  ( float_image_t *A, 
    int32_t x, 
    int32_t y, 
    int32_t nwx, 
    int32_t nwy, 
    bool_t rep, 
    float v[] 
  )
  {
    int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1]; 
    int32_t NY = (int32_t)A->sz[2];
    
    demand(nwx % 2 == 1, "{nwx} must be odd");
    demand(nwy % 2 == 1, "{nwy} must be odd");
    int32_t rwx = nwx/2; /* Half-width in X. */
    int32_t rwy = nwy/2; /* Half-width in Y. */
    
    /* Scan the window pixels: */
    int32_t kw = 0;
    for (int32_t iy = -rwy; iy <= rwy; iy++)
      { for (int32_t ix = -rwx; ix <= rwx; ix++)
          { /* Compute the sample's indices {xp,yp} in input img, or -1 if non-existant: */
            int32_t xp = x + ix;
            int32_t yp = y + iy;
            if (rep)
              { /* Map invalid {xp,yp} to nearest pixel in domain: */
                if (yp < 0) { yp = 0; } else if (yp >= NY) { yp = NY - 1; }
                if (xp < 0) { xp = 0; } else if (xp >= NX) { xp = NX - 1; }
              }
            /* Set invalid {xp,yp} to -1: */
            if ((xp < 0) || (xp >= NX) || (yp < 0) || (yp >= NY)) { xp = yp = -1; }

            /* Scan the channels: */
            for (int32_t c = 0; c < NC; c++)
              { /* Get and set the value: */  
                v[kw] = ( yp < 0 ? NAN : float_image_get_sample(A, c, xp, yp) );
                kw++;
              }
          }
      }
    
  }

void float_image_mix_channels
  ( double sA,
    float_image_t *A,
    int32_t cA, 
    double sB,
    float_image_t *B,
    int32_t cB, 
    float_image_t *R,
    int32_t cR
  )
  { int32_t NCR = (int32_t)R->sz[0];
    demand((cR >= 0) && (cR < NCR), "invalid {R} channel");
    int32_t NCA = (int32_t)A->sz[0];
    demand((cA >= 0) && (cA < NCA), "invalid {A} channel");
    int32_t NCB = (int32_t)B->sz[0];
    demand((cB >= 0) && (cB < NCB), "invalid {B} channel");
    int32_t NX = (int32_t)R->sz[1]; 
    int32_t NY = (int32_t)R->sz[2];
    float_image_check_size(A, -1, NX, NY);
    float_image_check_size(B, -1, NX, NY);
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float vA = float_image_get_sample(A, cA, x, y);
            float vB = float_image_get_sample(B, cB, x, y);
            float_image_set_sample(R, cR, x, y, (float)(sA*vA + sB*vB));
          }
      }
  }

void float_image_fill(float_image_t *A, float v)
  { int32_t NC = (int32_t)A->sz[0];
    for (int32_t c = 0; c < NC; c++) { float_image_fill_channel(A, c, v); }
  }

void float_image_fill_pixels(float_image_t *A, float v[])
  { int32_t NC = (int32_t)A->sz[0];
    for (int32_t c = 0; c < NC; c++) { float_image_fill_channel(A, c, v[c]); }
  }

void float_image_assign(float_image_t *A, float_image_t *V)
  { int32_t NC = (int32_t)A->sz[0];
    demand(V->sz[0] == NC, "incompatible channel counts");
    for (int32_t c = 0; c < NC; c++) { float_image_set_channel(A, c, V, c); }
  }

void float_image_fill_rectangle
  ( float_image_t *A, 
    int32_t xmin, 
    int32_t xmax, 
    int32_t ymin, 
    int32_t ymax, 
    float v
  )
  { if ((xmin > xmax) || (ymin > ymax)) { return; }
    demand((xmin >= 0) && (xmax < A->sz[1]), "invalid column range");
    demand((ymin >= 0) && (ymax < A->sz[2]), "invalid row range");
    float *prow = float_image_get_sample_address(A, 0, xmin, ymin);
    for (int32_t y = ymin; y <= ymax; y++)
      { float *ppix = prow;
        for (int32_t x = xmin; x <= xmax; x++)
          { float *psmp = ppix;
            for (int32_t c = 0; c < A->sz[0]; c++) { (*psmp) = v; psmp += A->st[0]; }
            ppix += A->st[1];
          }
        prow += A->st[2];
      }
  }

void float_image_fill_channel_rectangle
  ( float_image_t *A, 
    int32_t c, 
    int32_t xmin, 
    int32_t xmax, 
    int32_t ymin, 
    int32_t ymax,
    float v
  )
  { if ((xmin > xmax) || (ymin > ymax)) { return; }
    demand((c >= 0) && (c < A->sz[0]), "invalid channel");
    demand((xmin >= 0) && (xmax < A->sz[1]), "invalid column range");
    demand((ymin >= 0) && (ymax < A->sz[2]), "invalid row range");
    float *prow = float_image_get_sample_address(A, c, xmin, ymin);
    for (int32_t y = ymin; y <= ymax; y++)
      { float *psmp = prow;
        for (int32_t x = xmin; x <= xmax; x++) { (*psmp) = v; psmp += A->st[1]; }
        prow += A->st[2];
      }
  }

void float_image_fill_rectangle_pixels
  ( float_image_t *A, 
    int32_t xmin, 
    int32_t xmax, 
    int32_t ymin, 
    int32_t ymax, 
    float v[]
  )
  { if ((xmin > xmax) || (ymin > ymax)) { return; }
    demand((xmin >= 0) && (xmax < A->sz[1]), "invalid column range");
    demand((ymin >= 0) && (ymax < A->sz[2]), "invalid row range");
    float *prow = float_image_get_sample_address(A, 0, xmin, ymin);
    for (int32_t y = ymin; y <= ymax; y++)
      { float *ppix = prow;
        for (int32_t x = xmin; x <= xmax; x++)
          { float *psmp = ppix;
            for (int32_t c = 0; c < A->sz[0]; c++) { (*psmp) = v[c]; psmp += A->st[0]; }
            ppix += A->st[1];
          }
        prow += A->st[2];
      }
  }

void float_image_assign_channel_rectangle
  ( float_image_t *A, 
    int32_t cA, 
    int32_t xminA, int32_t xmaxA,  
    int32_t yminA, int32_t ymaxA,  
    float_image_t *V,  
    int32_t cV,  
    int32_t xminV, 
    int32_t yminV 
  )
  { if ((xminA > xmaxA) || (yminA > ymaxA)) { return; }
    int32_t xmaxV = xminV + xmaxA - xminA;
    int32_t ymaxV = yminV + ymaxA - yminA;
    demand((cA >= 0) && (cA < A->sz[0]), "invalid A channel");
    demand((cV >= 0) && (cV < A->sz[0]), "invalid V channel");
    demand((xminA >= 0) && (xmaxA < A->sz[1]) && (xmaxV < V->sz[1]), "invalid column range");
    demand((yminA >= 0) && (ymaxA < A->sz[2]) && (ymaxV < V->sz[2]), "invalid row range");
    float *prowA = float_image_get_sample_address(A, cA, xminA, yminA);
    float *prowV = float_image_get_sample_address(V, cV, xminV, yminV);
    for (int32_t y = yminA; y <= ymaxA; y++)
      { float *psmpA = prowA;
        float *psmpV = prowV;
        for (int32_t x = xminA; x <= xmaxA; x++) 
          { (*psmpA) = (*psmpV); psmpA += A->st[1]; psmpV += V->st[1]; }
        prowA += A->st[2];
        prowV += V->st[2];
      }
  }

void float_image_assign_rectangle
  ( float_image_t *A,  
    int32_t xminA, int32_t xmaxA,  
    int32_t yminA, int32_t ymaxA,  
    float_image_t *V,  
    int32_t xminV,  
    int32_t yminV 
  )
  { if ((xminA > xmaxA) || (yminA > ymaxA)) { return; }
    int32_t xmaxV = xminV + xmaxA - xminA;
    int32_t ymaxV = yminV + ymaxA - yminA;
    demand(A->sz[0] == V->sz[0], "channel counts must agree");
    int32_t NC = (int32_t)A->sz[0];
    demand((xminA >= 0) && (xmaxA < A->sz[1]) && (xmaxV < V->sz[1]), "invalid column range");
    demand((yminA >= 0) && (ymaxA < A->sz[2]) && (ymaxV < V->sz[2]), "invalid row range");
    float *prowA = float_image_get_sample_address(A, 0, xminA, yminA);
    float *prowV = float_image_get_sample_address(V, 0, xminV, yminV);
    for (int32_t y = yminA; y <= ymaxA; y++)
      { float *ppixA = prowA;
        float *ppixV = prowV;
        for (int32_t x = xminA; x <= xmaxA; x++)
          { float *psmpA = ppixA;
            float *psmpV = ppixV;
            for (int32_t c = 0; c < NC; c++) 
              { (*psmpA) = (*psmpV); psmpA += A->st[0]; psmpV += V->st[0]; }
            ppixA += A->st[1];
            ppixV += V->st[1];
          }
        prowA += A->st[2];
        prowV += V->st[2];
      }
  }

void float_image_get_gradient_sobel
  ( float_image_t *A, 
    int32_t c, 
    int32_t x, 
    int32_t y,
    double *fxP, 
    double *fyP
  )
  {
    int32_t NX = (int32_t)A->sz[1]; 
    int32_t NY = (int32_t)A->sz[2];
    /* Compute the horizontal and vertical derivatives in channel {c}: */
    double sum_wy_fx = 0.0;
    double sum_wy = 0.0;
    double sum_wx_fy = 0.0;
    double sum_wx = 0.0;
    for (int32_t d = -1; d <= +1; d++)
      { /* Sobel weight for slice {d}: */
        double wd = (d == 0 ? 2.0 : 1.0);
        if ((x-1 >= 0) && (x+1 < NX))
          { /* Accumulate the horizontal derivative: */
            int32_t yd = y + d;
            double nwy = ((yd < 0) || (yd >= NY) ? 0.0 : wd);
            if (nwy > 0.0)
              { double fxp = float_image_get_sample(A, c, x+1, yd);
                double fxm = float_image_get_sample(A, c, x-1, yd);
                double fx = (fxp - fxm)/2;
                sum_wy_fx += nwy*fx;
                sum_wy += nwy;
              }
          }

        if ((y-1 >= 0) && (y+1 < NY))
          { /* Accumulate the vertical derivative: */
            int32_t xd = x + d;
            double nwx = ((xd < 0) || (xd >= NX) ? 0.0 : wd);
            if (nwx > 0.0)
              { double fyp = float_image_get_sample(A, c, xd, y+1);
                double fym = float_image_get_sample(A, c, xd, y-1);
                double fy = (fyp - fym)/2;
                sum_wx_fy += nwx*fy;
                sum_wx += nwx;
              }
          }
      }

    /* Compute the derivatives for pixel {x,y}: */
    (*fxP) = (sum_wy == 0 ? 0.0 : sum_wy_fx/sum_wy);
    (*fyP) = (sum_wx == 0 ? 0.0 : sum_wx_fy/sum_wx);
  }

void float_image_get_local_avg_var
  ( float_image_t *A, 
    int32_t c, 
    int32_t x, 
    int32_t y,
    int32_t hw,
    double wt[],
    double *avgP, 
    double *varP
  )
  {
    int32_t NX = (int32_t)A->sz[1]; 
    int32_t NY = (int32_t)A->sz[2];
    double sum_w = 0.0;
    double sum_w_f = 0.0;
    double sum_w_f2 = 0.0;
    for (int32_t dx = -hw; dx <= +hw; dx++)
      { int32_t xs = x + dx;
        double wx = ((xs < 0) || (xs >= NX) ? 0 : wt[dx + hw]); 
        for (int32_t dy = -hw; dy <= +hw; dy++)
          { int32_t ys = y + dy;
            double wy = ((ys < 0) || (ys >= NY) ? 0 : wt[dy + hw]); 
            double w = wx*wy;
            if (w > 0)
              { double f = float_image_get_sample(A, c, xs, ys);
                if (! isnan(f))
                  { sum_w += w;
                    sum_w_f += w*f;
                    sum_w_f2 += w*f*f;
                  }
              }
          }
      }
    /* Compute the local average and variance (NAN if total weight is zero): */
    double avg = sum_w_f/sum_w;
    double var = sum_w_f2/sum_w - avg*avg;
    (*avgP) = avg;
    (*varP) = var;
  }

void float_image_local_avg_var
  ( float_image_t *A, 
    int32_t cA, 
    int32_t hw,
    double wt[],
    float_image_t *M, 
    int32_t cM,
    float_image_t *V, 
    int32_t cV
  )
  {
    /* Get {A} size and check chennels: */
    int32_t NCA = (int32_t)A->sz[0];
    demand((cA >= 0) && (cA < NCA), "invalid {A} channel");
    int32_t NXA = (int32_t)A->sz[1]; 
    int32_t NYA = (int32_t)A->sz[2];
    /* Get {M} size and check chennels: */
    int32_t NXM, NYM;
    if (M != NULL)
      { int32_t NCM = (int32_t)M->sz[0];
        demand((cM >= 0) && (cM < NCM), "invalid {M} channel");
        NXM = (int32_t)M->sz[1];
        NYM = (int32_t)M->sz[2]; 
      }
    else
      { NXM = NYM = 0; }
    /* Get {V} size and check chennels: */
    int32_t NXV, NYV;
    if (V != NULL)
      { int32_t NCV = (int32_t)V->sz[0];
        demand((cV >= 0) && (cV < NCV), "invalid {V} channel");
        NXV = (int32_t)V->sz[1]; 
        NYV = (int32_t)V->sz[2]; 
      }
    else
      { NXV = NYV = 0; }
      
    /* Determine computation domain and ofsets from it: */
    int32_t NX = (int32_t)imax(NXM, NXV);
    int32_t NY = (int32_t)imax(NYM, NYV);
    int32_t dXA = (NX - NXA)/2;
    int32_t dYA = (NY - NYA)/2;
    int32_t dXM = (NX - NXM)/2;
    int32_t dYM = (NY - NYM)/2;
    int32_t dXV = (NX - NXV)/2;
    int32_t dYV = (NY - NYV)/2;
    
    for (int32_t x = 0; x < NX; x++)
      { for (int32_t y = 0; y < NY; y++)
          { /* Extract local mean and variance: */
            double avg, var;
            float_image_get_local_avg_var(A, cA, x - dXA, y - dYA, hw, wt, &avg, &var);
            if (M != NULL)
              { int32_t xM = x - dXM;
                int32_t yM = y - dYM;
                if ((xM >= 0) && (xM < NXM) && (yM >= 0) && (yM < NYM))
                  { float_image_set_sample(M, cM, xM, yM, (float)avg); }
              }
            if (V != NULL) 
              { int32_t xV = x - dXV;
                int32_t yV = y - dYV;
                if ((xV >= 0) && (xV < NXV) && (yV >= 0) && (yV < NYV))
                  { float_image_set_sample(V, cV, xV, yV, (float)var); }
              }
          }
      }
  }
  
double float_image_get_dilated
  ( float_image_t *A, 
    int32_t c, 
    int32_t x, 
    int32_t y,
    int32_t hw,
    double wt[]
  )
  {
    int32_t NX = (int32_t)A->sz[1]; 
    int32_t NY = (int32_t)A->sz[2];
    double max_w_a = -INF;
    for (int32_t dx = -hw; dx <= +hw; dx++)
      { int32_t xs = x + dx;
        double wx = ((xs < 0) || (xs >= NX) ? 0 : wt[dx + hw]); 
        for (int32_t dy = -hw; dy <= +hw; dy++)
          { int32_t ys = y + dy;
            double wy = ((ys < 0) || (ys >= NY) ? 0 : wt[dy + hw]); 
            double w = wx*wy;
            if (w > 0)
              { double a = float_image_get_sample(A, c, xs, ys);
                double w_a = w*a;
                if (w_a > max_w_a) { max_w_a = w_a; }
              }
          }
      }
    return max_w_a;
  }

void float_image_make_grayscale(float_image_t *A)
  {
    int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    if (NC == 1) { return; }
    demand(NC == 3, "channel count must be 1 or 3");
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { frgb_t pix = fic_get_frgb_pixel(A, 0, 1, 2, x, y);
            double value = frgb_get_Y(&pix);
            for (int32_t c = 0; c < NC; c++) { pix.c[c] = (float)value; }
            fic_set_frgb_pixel(A, 0, 1, 2, x, y, &pix);
          }
      }
  }

void float_image_apply_gamma(float_image_t *A, int32_t c, double gamma, double bias)
  { if (gamma == 1) { return; }
    demand(gamma > 0, "gamma must be positive"); 
    int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    demand((c >= 0) && (c < NC), "invalid channel index");
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float *pA = float_image_get_sample_address(A, c, x, y);
            (*pA) = sample_conv_gamma((*pA), gamma, bias);
          }
      }
  }

void float_image_log_scale(float_image_t *A, int32_t c, double vref, double base)
  { demand((isfinite(vref)) && (vref > 0), "invalid vref");
    demand((isfinite(base)) && (base != 1) && (base > 0), "invalid base");
    double logBase = (base == M_E ? 1.0 : log(base));
    int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    demand((c >= 0) && (c < NC), "invalid channel index");
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float *pA = float_image_get_sample_address(A, c, x, y);
            (*pA) = sample_conv_log((*pA), vref, logBase);
          }
      }
  }

void float_image_undo_log_scale(float_image_t *A, int32_t c, double vref, double base)
  { demand((isfinite(vref)) && (vref > 0), "invalid vref");
    demand((isfinite(base)) && (base != 1) && (base > 0), "invalid base");
    double logBase = (base == M_E ? 1.0 : log(base));
    int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    demand((c >= 0) && (c < NC), "invalid channel index");
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float *pA = float_image_get_sample_address(A, c, x, y);
            (*pA) = sample_conv_undo_log((*pA), vref, logBase);
          }
      }
  }

void float_image_rescale_samples(float_image_t *A, int32_t c, float a0, float a1, float z0, float z1)
  { int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    demand((c >= 0) && (c < NC), "invalid channel index");
    double scale = (z1 - z0)/(a1 - a0);
    demand(! isnan(scale), "scale factor is NaN");
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float *v = float_image_get_sample_address(A, c, x, y);
            (*v) = (float)(z0 + scale*((*v) - a0));
          }
      }
  }

void float_image_square_samples(float_image_t *A, int32_t c)
  { int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    demand((c >= 0) && (c < NC), "invalid channel index");
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float *smp = float_image_get_sample_address(A, c, x, y);
            float v = *smp;
            if (!isnan(v)) { (*smp) = v*v; }
          }
      }
  }

double float_image_compute_sample_sum(float_image_t *A, int32_t c)
  { int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    if ((c < 0) || (c >= NC))
      { /* Invalid channel, all samples are zero (or there are no samples). */
        return 0.0;
      }

    /* Compute the sum of samples {sum}: */
    double sum = 0;
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { double v = float_image_get_sample(A, c, x, y);
            if (isfinite(v))
              { /* Neither {±INF} nor {NAN}: */ sum += v; }
          }
      }
    return sum;
  }

double float_image_compute_total_energy(float_image_t *A, int32_t c, double avg)
  { int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    if ((c < 0) || (c >= NC))
      { /* Invalid channel, all samples are zero (or there are no samples). */
        return 0.0;
      }

    /* Compute the sum of squares {sum2}: */
    double sum2 = 0;
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { double v = float_image_get_sample(A, c, x, y);
            if (isfinite(v))
              { /* Neither {±INF} nor {NAN}: */ v -= avg; sum2 += v*v; }
          }
      }
    return sum2;
  }

void float_image_replace_nan_samples(float_image_t *A, int32_t c, float v)
  { int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    demand((c >= 0) && (c < NC), "invalid channel index");
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x = 0; x < NX; x++)
          { float *p = float_image_get_sample_address(A, c, x, y);
            if (isnan(*p)) { (*p) = v; }
          }
       }
  }

void float_image_compute_sample_avg_dev(float_image_t *A, int32_t c, double *avg, double *dev)
  { int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    if ((c < 0) || (c >= NC))
      { /* Invalid channel, all samples are zero (or there are no samples). */
        (*avg) = 0.0;
        (*dev) = 0.0;
      }
    else
      { /* Valid channel. */
        /* Compute the sample average {sa}: */
        double sum = 0;
        int32_t tot = 0;
        for (int32_t y = 0; y < NY; y++)
          { for (int32_t x = 0; x < NX; x++)
              { double v = float_image_get_sample(A, c, x, y);
                if (isfinite(v))
                  { /* Neither {±INF} nor {NAN}: */ sum += v; tot++; }
              }
          }
        double sa = (tot <= 0 ? 0.0 : sum/tot);
        /* Compute the sample variance {sv}: */
        sum = 0;
        for (int32_t y = 0; y < NY; y++)
          { for (int32_t x = 0; x < NX; x++)
              { double v = float_image_get_sample(A, c, x, y);
                if (isfinite(v))
                  { /* Neither {±INF} nor {NAN}: */ v -= sa; sum += v*v; }
              }
          }
        double sv = (tot <= 1 ? 0.0 : sum/(tot - 1));
        /* Return the results: */
        if (avg != NULL) { (*avg) = sa; }
        if (dev != NULL) { (*dev) = sqrt(sv); }
      }
  }

void float_image_update_sample_range(float_image_t *A, int32_t c, float *vMin, float *vMax)
  { 
    if ((vMin == NULL) && (vMax == NULL)) { /* Nothing to do: */ return; }
    int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    if ((c < 0) || (c >= NC))
      { /* Invalid channel, all samples are zero. (But there may be no samples!) */
        if ((NX > 0) && (NY > 0))
          { float v = 0.0;
            if (vMin != NULL) { if (v < (*vMin)) { (*vMin) = v; } }
            if (vMax != NULL) { if (v > (*vMax)) { (*vMax) = v; } }
          }
      }
    else
      { /* Valid channel. */
        for (int32_t y = 0; y < NY; y++)
          { for (int32_t x = 0; x < NX; x++)
              { float v = float_image_get_sample(A, c, x, y);
                if (isfinite(v))
                  { /* Neither {±INF} nor {NAN}: */ 
                    if (vMin != NULL) { if (v < (*vMin)) { (*vMin) = v; } }
                    if (vMax != NULL) { if (v > (*vMax)) { (*vMax) = v; } }
                  }
              }
          }
      }
  }

float float_image_spectrum_max_sample(float_image_t *A, int32_t c, bool_t centered)
  { 
    int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    /* Determine the constant term position: */
    int32_t x0 = (centered ? NX/2 : 0);
    int32_t y0 = (centered ? NY/2 : 0);
    /* Collect the max sample value: */
    float vMax = 0.0;
    if ((c >= 0) && (c < NC))
      { /* Valid channel. */
        for (int32_t y = 0; y < NY; y++)
          { for (int32_t x = 0; x < NX; x++)
              { if ((x != x0) || (y != y0))
                  { float v = float_image_get_sample(A, c, x, y);
                    if (isfinite(v))
                      { assert(v >= 0.0);
                        vMax = fmaxf(vMax, v);
                      }
                  }
              }
          }
      }
    return vMax;
  }

void float_image_flip_x(float_image_t *A, int32_t c, int32_t ix)
  { int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    demand((c >= 0) && (c < NC), "invalid channel index");
    
    /* Reduce the index {ix} to the range {0..NX-1}: */
    ix = (ix % NX); if (ix < 0) { ix += NX; }
    
    /* Compute the first and last indices for the X loop: */
    int32_t fst_x = ix / 2 + 1; 
    int32_t lst_x = (ix + NX - 1) / 2;
    for (int32_t y = 0; y < NY; y++)
      { for (int32_t x1 = fst_x; x1 <= lst_x; x1++)
          { /* Compute the partner: */
            int32_t x2 = ix - x1; if (x2 < 0) { x2 += NX; }
            assert(x2 != x1);
            assert(x2 < NX);
            float *smp1 = float_image_get_sample_address(A, c, x1, y);
            float *smp2 = float_image_get_sample_address(A, c, x2, y);
            float tmp = (*smp1); (*smp1) = (*smp2); (*smp2) = tmp;
          }
      }
  }
  
void float_image_flip_y(float_image_t *A, int32_t c, int32_t iy)
  { int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    demand((c >= 0) && (c < NC), "invalid channel index");
    
    /* Reduce the index {iy} to the range {0..NY-1}: */
    iy = (iy % NY); if (iy < 0) { iy += NY; }
    
    /* Compute the first and last indices for the Y loop: */
    int32_t fst_y = iy / 2 + 1; 
    int32_t lst_y = (iy + NY - 1) / 2;
    for (int32_t x = 0; x < NX; x++)
      { for (int32_t y1 = fst_y; y1 <= lst_y; y1++)
          { /* Compute the partner: */
            int32_t y2 = iy - y1; if (y2 < 0) { y2 += NY; }
            assert(y2 < NY);
            assert(y2 != y1);
            float *smp1 = float_image_get_sample_address(A, c, x, y1);
            float *smp2 = float_image_get_sample_address(A, c, x, y2);
            float tmp = (*smp1); (*smp1) = (*smp2); (*smp2) = tmp;
          }
      }
  }

void float_image_shift(float_image_t *A, int32_t c, int32_t dx, int32_t dy)
  { int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    demand((c >= 0) && (c < NC), "invalid channel index");
    
    /* Reduce the shifts to the ranges {0..NX-1}, {0..NY-1}: */
    dx = (dx % NX); if (dx < 0) { dx += NX; }
    dy = (dy % NY); if (dy < 0) { dy += NY; }
    
    /* Uses the two-flip method to obtain a shift. This may be faster
      than decomposition, depending on cache size. */
    
    if (dx != 0)
      { /* Apply two flips in X: */
        float_image_flip_x(A, c, 0);
        float_image_flip_x(A, c, dx);
      }
      
    if (dy != 0)
      { /* Apply two flips in Y: */
        float_image_flip_y(A, c, 0);
        float_image_flip_y(A, c, dy);
      }
  }

#define float_image_max_samples (1024u*1024u*1024u)
  /* Maximum total elements (2^30), for sanity checks. */

float_image_t *float_image_new (int32_t NC, int32_t NX, int32_t NY)
  { /* Sanity checks: */
    demand((NC >= 0) && (NC < float_image_max_size), "too many channels");
    demand((NX >= 0) && (NX < float_image_max_size), "too many columns");
    demand((NY >= 0) && (NY < float_image_max_size), "too many rows");
    ix_count_t NS = ((ix_count_t)NC)*NX*NY;
    demand(NS < float_image_max_samples, "too many samples");
    /* Allocate header: */
    float_image_t *A = (float_image_t *)notnull(malloc(sizeof(float_image_t)), "no mem");
    A->sz[0] = NC; A->st[0] = (NC < 2 ? 0 : 1);
    A->sz[1] = NX; A->st[1] = (NX < 2 ? 0 : NC);
    A->sz[2] = NY; A->st[2] = (NY < 2 ? 0 : NC*NX);
    A->bp = 0;
    if ((NC == 0) || (NX == 0) || (NY == 0))
      { A->sample = (float *)NULL; }
    else
      { A->sample = (float *)notnull(malloc(NS*sizeof(float)), "no mem"); }
    (void)ix_parms_are_valid(3, A->sz, A->bp, A->st, /*die:*/ TRUE);
    return A;
  }
  
float_image_t *float_image_copy (float_image_t *A)
  { int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    float_image_t *B = float_image_new(NC, NX, NY);
    for(int32_t y = 0; y < NY; y++)
      { for(int32_t x = 0; x < NX; x++)
          { for (int32_t c = 0; c < NC; c++) 
              { float v = float_image_get_sample(A, c, x, y);
                float_image_set_sample(B, c, x, y, v);
              }
          }
      }
    return B;
  }

float_image_t *float_image_crop
  ( float_image_t *A, 
    int32_t cLo,
    int32_t cHi,
    int32_t xLo,
    int32_t xHi,
    int32_t yLo,
    int32_t yHi,
    float bg
  )
  { /* Grab the dimensions of {A}: */
    int32_t NCA, NXA, NYA;
    float_image_get_size(A, &NCA, &NXA, &NYA);
    /* Compute the dimensions of the result {R} and allocate it: */
    int32_t NCR = (int32_t)imax(0, cHi - cLo);
    int32_t NXR = (int32_t)imax(0, xHi - xLo);
    int32_t NYR = (int32_t)imax(0, yHi - yLo);
    float_image_t *R = float_image_new(NCR, NXR, NYR);
    /* Fill {R} as appropriate: */
    for (int32_t c = 0; c < NCR; c++)
      { int32_t cA = cLo + c;
        bool_t cOK = ((cA >= 0) && (cA < NCA));
        for (int32_t x = 0; x < NXR; x++)
          { int32_t xA = xLo + x;
            bool_t xOK = ((xA >= 0) && (xA < NXA));
            for (int32_t y = 0; y < NYR; y++)
              { int32_t yA = yLo + y;
                bool_t yOK = ((yA >= 0) && (yA < NYA));
                float v;
                if (cOK && xOK && yOK)
                  { v = float_image_get_sample(A, cA, xA, yA); }
                else
                  { v = bg; }
                float_image_set_sample(R, c, x, y, v);
              }
          }
      }
    return R;
  }

void float_image_free(float_image_t *A)
  { if (A == NULL) return;
    if (A->sample != NULL) { free(A->sample); }
    free(A);
  }
  
void float_image_get_size(float_image_t *A, int32_t *NC, int32_t *NX, int32_t *NY)
  { if (NC != NULL) { (*NC) = (int32_t)A->sz[0]; }
    if (NX != NULL) { (*NX) = (int32_t)A->sz[1]; }
    if (NY != NULL) { (*NY) = (int32_t)A->sz[2]; }
  }
  
void float_image_check_size(float_image_t *A, int32_t NC, int32_t NX, int32_t NY)
  { if (NC >= 0) { demand(((int32_t)A->sz[0]) == NC, "wrong number of channels"); }
    if (NX >= 0) { demand(((int32_t)A->sz[1]) == NX, "wrong number of columns"); }
    if (NY >= 0) { demand(((int32_t)A->sz[2]) == NY, "wrong number of rows"); }
  }

#define float_image_file_version "2006-03-25"

void float_image_write(FILE *wr, float_image_t *A)
  { 
    int32_t NC = (int32_t)A->sz[0];
    int32_t NX = (int32_t)A->sz[1];
    int32_t NY = (int32_t)A->sz[2];
    filefmt_write_header(wr, "float_image_t", float_image_file_version);
    fprintf(wr, "NC = %d\n", NC);
    fprintf(wr, "NX = %d\n", NX);
    fprintf(wr, "NY = %d\n", NY);
    for(int32_t y = 0; y < NY; y++)
      { if (y > 0) { fprintf(wr, "\n"); }
        for(int32_t x = 0; x < NX; x++)
          { fprintf(wr, "%5d %5d", x, y);
            for (int32_t c = 0; c < NC; c++) 
              { float v = float_image_get_sample(A, c, x, y);
                fprintf(wr, " %+14.7e", v);
              }
            fprintf(wr, "\n");
          }
      }
    filefmt_write_footer(wr, "float_image_t");
    fflush(wr);
  }

float_image_t *float_image_read(FILE *rd)
  {
    filefmt_read_header(rd, "float_image_t", float_image_file_version);
    int32_t NC = nget_int32(rd, "NC"); fget_eol(rd);
    int32_t NX = nget_int32(rd, "NX"); fget_eol(rd);
    int32_t NY = nget_int32(rd, "NY"); fget_eol(rd);
    float_image_t *A = float_image_new(NC, NX, NY);
    for (int32_t y = 0; y < NY; y++)
      { if (y > 0) { fget_eol(rd); }
        for (int32_t x = 0; x < NX; x++)
          { int32_t xr = fget_int32(rd);
            int32_t yr = fget_int32(rd);
            demand((xr == x) && (yr == y), "bad pixel indices");
            /* Read channels of pixel {(x,y)}: */
            for (int32_t c = 0; c < NC; c++)
              { double v = fget_double(rd);
                float_image_set_sample(A, c, x, y, (float)v);
              }
            fget_eol(rd);
          }
      }
    filefmt_read_footer(rd, "float_image_t");
    return A;
  }

void float_image_debug_pixel(char *label, double x, double y, int32_t chns, float f[], char *tail)
  { 
    fprintf(stderr, "%s(%9.4f,%9.4f) = (", label, x, y);
    for (int32_t ich = 0; ich < chns; ich++) 
      { fprintf(stderr, " %7.4f", f[ich]); }
    fprintf(stderr, " )%s", tail);
  }

void float_image_debug_double_pixel(char *label, double x, double y, int32_t chns, double v[], char *tail)
  { 
    fprintf(stderr, "%s(%9.4f,%9.4f) = (", label, x, y);
    for (int32_t ich = 0; ich < chns; ich++) 
      { fprintf(stderr, " %7.4f", v[ich]); }
    fprintf(stderr, " )%s", tail);
  }
