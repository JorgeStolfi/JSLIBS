/* See {neuromat_{eeg_|}image.h}. */
/* Last edited on 2021-08-25 02:25:58 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
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

float_image_t *neuromat_eeg_image_make_idealized_scalp_mask
  ( int32_t NX, 
    int32_t NY, 
    r2_t *ctr, 
    r2_t *rad
  )
  {
    float_image_t *img = float_image_new(1, NX, NY);
    float_image_fill(img, 0.0);
    float vfill = 1.0; /* Ellipse fill color. */
    float vdraw = NAN; /* Outline color. */
    float hwd = 0.0; /* Half-width of drawing pen. */
    int32_t msub = 3; /* Subsampling order. */
    (void)float_image_paint_ellipse_aligned(img, 0, ctr, rad, hwd, vfill, vdraw, msub);
    return img;
  }
  
void neuromat_eeg_image_compute_pot_field
  ( int32_t ne, 
    double val[], 
    float_image_t *bas[],
    float_image_t *msk, 
    float_image_t *img
  )
  {
    int32_t NC = (int32_t)img->sz[0];
    int32_t NX = (int32_t)img->sz[1];
    int32_t NY = (int32_t)img->sz[2];
    demand(NC == 1, "field image should be monochrome");
    
    assert(msk->sz[0] == 1);
    assert(msk->sz[1] == NX);
    assert(msk->sz[2] == NY);

    for (int32_t iy = 0; iy < NY; iy++)
      { for (int32_t ix = 0; ix < NX; ix++)
          { double wxy = float_image_get_sample(msk, 0, ix, iy);
            double smp;
            if (wxy == 0)
              { smp = 0.0; }
            else
              { smp = 0.0;
                for (int32_t ie = 0; ie < ne; ie++) 
                  { double basi = float_image_get_sample(bas[ie], 0, ix, iy);
                    smp += val[ie] * basi;
                  }
                assert(! isnan(smp));
              }
            float_image_set_sample(img, 0, ix, iy, (float)smp);
          }
      }
  }

void neuromat_eeg_image_paint_electrodes
  ( int32_t ne,
    r2_t pos[], 
    double drad,
    double hwd,
    int32_t ie_spec, 
    frgb_t *fcfill, 
    frgb_t *fcdraw, 
    float_image_t *img
  )
  { 
    int32_t NC = (int32_t)img->sz[0];
    demand(NC == 4, "image should be RGBA");
    if ((fcfill == NULL) && (fcdraw == NULL)) { return; }
    bool_t round = TRUE; /* Round dots? */
    bool_t diagonal = FALSE; /* Diagonal crosses? (Not used.) */
    int32_t msub = 3; /* Subsampling order. */
    int32_t cop = 3; /* Index of opacity channel. */
    for (int32_t ie = 0; ie < ne; ie++) 
      { r2_t *posi = &(pos[ie]);
        double cx = posi->c[0];
        double cy = posi->c[1];
        frgb_t *cfill = (ie == ie_spec ? fcdraw : fcfill); /* Fill color, or NULL. */
        frgb_t *cdraw = (ie == ie_spec ? fcfill : fcdraw); /* Draw color, or NULL. */
        for (int32_t c = 0; c < NC; c++)
          { if (c == cop) 
              { /* Draw the opacity mask: */
                float vfill = (cfill == NULL ? NAN : 1.000f);
                float vdraw = (cdraw == NULL ? NAN : 1.000f);
                (void)float_image_paint_dot
                  ( img, c, cx, cy, drad, hwd, 
                    round, diagonal, vfill, vdraw, msub
                  );
              }
            else 
              { /* Draw the colors, but don't antialias edged: */
                float vfill = (cfill == NULL ? NAN : cfill->c[c]); /* Fill color. */
                float vdraw = ((cdraw == NULL) || (hwd <= 0) ? NAN : cdraw->c[c]); /* Draw color. */
                double dradx, hwdx; /* Dot parameters for the un-antialiased color: */
                if (isnan(vfill) && (! isnan(vdraw)))
                  { /* Empty outline: */
                    dradx = drad; hwdx = hwd + 1.5;
                  }
                else if ((! isnan(vfill)) && isnan(vdraw))
                  { /* Filled dot without outline: */
                    dradx = drad + 1.5; hwdx = 0.0;
                  }
                else 
                  { /* filled and drawn: */
                    dradx = drad + 1.5; hwdx = hwd + 1.5;
                  }
                (void)float_image_paint_dot
                  ( img, c, cx, cy, dradx, hwdx, 
                    round, diagonal, vfill, vdraw, msub
                  );
              }
          }
      }
  }

void neuromat_eeg_image_paint_timeline_bars  
  ( int32_t nt, 
    int32_t nc, 
    double **val, 
    int32_t nm, 
    int32_t ic_mark[], 
    int32_t hw,
    frgb_t fc[],
    float_image_t *img,
    int32_t *track_xloP, 
    int32_t *track_xszP, 
    int32_t track_y[]
  )
  {
    int32_t NC = (int32_t)img->sz[0];
    int32_t NX = (int32_t)img->sz[1];
    int32_t NY = (int32_t)img->sz[2];
    demand(NC == 4, "image should be RGBA");
    
    demand(hw >= 1, "invalid {hw}");
    int32_t nk = (nm < 1 ? 1 : nm); /* Number of tracks. */

    /* Choose the horizontal start and  of the tracks: */
    int32_t xmrg = (int32_t)imax(NX/25, 3*hw);  /* Margin pixels on sides. */
    int32_t xlo = xmrg;  /* Left offset to start of track (pixels). */
    int32_t xsz = NX - 2*xmrg; /* Horizontal extent of tracks (pixels). */
    demand(xsz >= 4*hw, "timeline image width is too small");
    (*track_xloP) = xlo;
    (*track_xszP) = xsz;
    
    /* Choose the vertical positions of the tracks and paint them: */
    int32_t ymrg = 2*hw;           /* Extra marging at bottom and top. */
    int32_t ythw = 2*hw;           /* Track half-width including highlights and tic marks. */
    int32_t ytw = 2*ythw;          /* Track width including highlights and tic marks. */
    int32_t ytsp = 1*hw;           /* Extra intertrack spacing. */
    int32_t ytstep = ytw + ytsp;   /* Y step between track axes. */
    int32_t y0 = ymrg + ytsp; /* Y coordinate of lowest timeline. */
    int32_t ytot = ymrg + ytsp   + ythw + (nk-1)*ytstep + ythw + ytsp + ymrg; /* Total height of tracks (pixels) incl margins. */
    demand(NY >= ytot, "image not tall enough for timeline bars"); 

    /* Paint the tracks and marker channel ranges: */
    frgb_t ftrack = (frgb_t){{ 0.500f, 0.500f, 0.500f }}; /* Color or track line. */
    for (int32_t ik = 0; ik < nk; ik++)
      { frgb_t fmark = fc[ik]; /* Color or tics and marks line. */
        track_y[ik] = y0 + ik*ytstep; 
        int32_t ic = (nm == 0 ? -1 : ic_mark[ik]); /* EEG channel index. */
        neuromat_image_paint_time_track(img, hw, xlo, xsz, track_y[ik], ftrack);
        if (ic >= 0) { neuromat_image_paint_marker_ranges(img, hw, nt, nc, val, ic, xlo, xsz, track_y[ik], fmark); }
      }
  }
    
void neuromat_eeg_image_paint_marker_dots
  ( int32_t nc, 
    double valt[],
    int32_t nm,
    int32_t ic_mark[], 
    double vscale[], 
    frgb_t fcpos[],
    frgb_t fcneg[],
    r2_t ctr_mark[],
    double rad_mark,
    float_image_t *img
  )
  { 
    int32_t NC = (int32_t)img->sz[0];
    demand(NC == 4, "image should be RGBA");

    float vdraw = NAN;   /* Outline color. (Not used) */
    float hwd = 0.0;     /* Half-width of drawing pen. */
    bool_t round = TRUE; /* Round dots? */
    bool_t diagonal = FALSE; /* Diagonal crosses? (Not used.) */
    int32_t msub = 4; /* Subsampling order. */
    for (int32_t km = 0; km < nm; km++)
      { int32_t ic = ic_mark[km]; /* Index of channel in data frames: */
        demand((ic >= 0) && (ic < nc), "invalid marker channel");
        /* Convert marker value to {[-1 _ +1]}: */
        r2_t *ctrk = &(ctr_mark[km]);
        double zk = valt[ic]/vscale[ic]; /* Marker value scaled to {[-1..+1]} */
        double sk = sqrt(fabs(zk));
        frgb_t *rgbk = (zk > 0 ? &(fcpos[km]) : &(fcneg[km]));
        for (int32_t c = 0; c < NC; c++)
          { double rkc = (c < 3 ? rad_mark + 1.5 : rad_mark);
            float vfill = (float)(c < 3 ? sk*rgbk->c[c] : 1.000);
            (void)float_image_paint_dot
              ( img, c, ctrk->c[0], ctrk->c[1], rkc, hwd, 
                round, diagonal, vfill, vdraw, msub
              );
          }
      }
  }

 
