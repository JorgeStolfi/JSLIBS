/* See {neuromat_{eeg_|}image.h}. */
/* Last edited on 2020-10-11 01:51:14 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <float_image.h>
#include <float_image_paint.h>
#include <frgb.h>
#include <frgb_path.h>
#include <affirm.h>

#include <jsfile.h>

#include <neuromat_eeg.h>
#include <neuromat_image.h>
#include <neuromat_eeg_image.h>

float_image_t *neuromat_eeg_image_schematic_head_mask(int NX, int NY, r2_t *ctr, r2_t *rad)
  {
    float_image_t *img = float_image_new(1, NX, NY);
    float_image_fill(img, 0.0);
    float vfill = 1.0; /* Ellipse fill color. */
    float vdraw = NAN; /* Outline color. */
    float hwd = 0.0; /* Half-width of drawing pen. */
    int msub = 3; /* Subsampling order. */
    (void)float_image_paint_ellipse_aligned(img, 0, ctr, rad, hwd, vfill, vdraw, msub);
    return img;
  }
  
void neuromat_eeg_image_draw_electrodes
  ( float_image_t *img,
    int ne, 
    r2_t pos[], 
    double drad, 
    int ilo, 
    float vlo[], 
    float vhi[]
  )
  { 
    int NC = (int)img->sz[0];
    int ie;
    float vdraw = NAN;   /* Outline color. (Not used) */
    float hwd = 0.0;     /* Half-width of drawing pen. */
    bool_t round = TRUE; /* Round dots? */
    bool_t diagonal = FALSE; /* Diagonal crosses? (Not used.) */
    int msub = 3; /* Subsampling order. */
    for (ie = 0; ie < ne; ie++) 
      { r2_t *posi = &(pos[ie]);
        double cx = posi->c[0];
        double cy = posi->c[1];
        float *vfill = (ie == ilo ? vlo : vhi); /* Fill color vector. */
        int c;
        for (c = 0; c < NC; c++)
          { (void)float_image_paint_dot
              ( img, c, cx, cy, drad, hwd, 
                round, diagonal, vfill[c], vdraw, msub
              );
          } 
      }
  }

void neuromat_eeg_image_paint_potentials
  ( int ne, 
    double val[], 
    float_image_t *bas[],
    float_image_t *msk, 
    int c, float_image_t *img
  )
  {
    int NC = (int)img->sz[0];
    int NX = (int)img->sz[1];
    int NY = (int)img->sz[2];
    
    assert(msk->sz[0] == 1);
    assert(msk->sz[1] == NX);
    assert(msk->sz[2] == NY);
    
    demand((c >= 0) && (c < NC), "invalid channel");

    int ie;
    int ix, iy;
    for (iy = 0; iy < NY; iy++)
      { for (ix = 0; ix < NX; ix++)
          { double wxy = float_image_get_sample(msk, 0, ix, iy);
            if (wxy != 0)
              { double sum = 0;
                for (ie = 0; ie < ne; ie++) 
                  { double basi = float_image_get_sample(bas[ie], 0, ix, iy);
                    sum += val[ie] * basi;
                  }
                assert(! isnan(sum));
                float *smp = float_image_get_sample_address(img, c, ix, iy);
                double val = wxy*sum + (1-wxy)*(*smp);
                (*smp) = (float)val;
              }
          }
      }
  }

float_image_t *neuromat_eeg_image_make_timeline_bar
  ( int NX, 
    int nt, 
    int nc, 
    double **val, 
    int nm, 
    int ic_mark[], 
    int hw,
    int *track_xloP, 
    int *track_xszP, 
    int track_y[]
  )
  {
    demand(hw >= 1, "invalid {hw}");
    int nk = (nm < 1 ? 1 : nm); /* Number of tracks. */

    /* Choose the horizontal start and  of the tracks: */
    int xmrg = (int)imax(NX/25, 3*hw);  /* Margin pixels on sides. */
    int xlo = xmrg;  /* Left offset to start of track (pixels). */
    int xsz = NX - 2*xmrg; /* Horizontal extent of tracks (pixels). */
    demand(xsz >= 4*hw, "timeline image width is too small");
    (*track_xloP) = xlo;
    (*track_xszP) = xsz;
    
    /* Choose the vertical positions of the tracks and paint them: */
    int ytw = 4*hw; /* Track width including margins. */
    int yts = 1*hw; /* Extra intertrack spacing. */
    int ytot = nk*ytw + (nk-1)*yts; /* Total height of tracks (pixels). */
    int NY = 2*ytw + ytot + 2*ytw; /* Total height of timeline image. */

    /* Build an image with the timeline bar, and get its column span {*bxminP..*bxmaxP}: */
    float_image_t *img = float_image_new(3, NX, NY);

    /* Paint the tracks and marker channel ranges: */
    int ymrg = (NY - ytot)/2;  /* Margin pixels at bottom. */
    frgb_t ftrack = (frgb_t){{ 0.500f, 0.500f, 0.500f }}; /* Color or track line. */
    frgb_t fmark = (frgb_t){{ 1.000f, 0.900f, 0.800f }}; /* Color or tics and marks line. */
    for (int ik = 0; ik < nk; ik++)
      { track_y[ik] = ymrg + ik*(ytw + yts); 
        neuromat_image_paint_time_track(img, hw, xlo, xsz, track_y[ik], &ftrack);
        if (nm > 0)
          { int ic = ic_mark[ik];
            neuromat_image_paint_marker_ranges(img, hw, nt, nc, val, ic, xlo, xsz, track_y[ik], &fmark);
          }
      }
    return img;
  }
    
void neuromat_eeg_image_paint_marker_dots
  ( float_image_t *frame,
    int nc, 
    double valt[],
    int nm,
    int ic_mark[], 
    double vscale[], 
    r2_t ctr_mark[],
    double rad_mark,
    int style
  )
  { int km;
    float vdraw = NAN;   /* Outline color. (Not used) */
    float hwd = 0.0;     /* Half-width of drawing pen. */
    bool_t round = TRUE; /* Round dots? */
    bool_t diagonal = FALSE; /* Diagonal crosses? (Not used.) */
    int msub = 4; /* Subsampling order. */
    for (km = 0; km < nm; km++)
      { int ic = ic_mark[km]; /* Index of channel in data frames: */
        demand((ic >= 0) && (ic < nc), "invalid marker channel");
        /* Convert marker value to {[-1 _ +1]}: */
        float zk = (float)(valt[ic]/vscale[ic]); /* Marker value scaled to {[-1..+1]} */
        frgb_t rgbk = frgb_path_map_signed(zk, 1, style);
        int c;
        for (c = 0; c < 3; c++)
          { r2_t *ctr = &(ctr_mark[km]);
            float vfill = rgbk.c[c];
            (void)float_image_paint_dot
              ( frame, c, ctr->c[0], ctr->c[1], rad_mark, hwd, 
                round, diagonal, vfill, vdraw, msub
              );
          }
      }
  }

 
