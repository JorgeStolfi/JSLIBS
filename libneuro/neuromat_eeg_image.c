/* See {neuromat_{eeg_|}image.h}. */
/* Last edited on 2021-08-28 21:11:22 by stolfi */

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

void neuromat_image_paint_marker_ranges
  ( float_image_t *img,
    int32_t hw, 
    int32_t nt, 
    int32_t nc,
    double **val,
    int32_t ic_mark, 
    int32_t xlo, 
    int32_t xsz,
    int32_t y, 
    frgb_t fc
  );
  /* Paints onto {img} tic marks and highlights showing the
    ranges of frames where the marker channel {ic} is nonzero.
    The marks are located on the timeline defined by {xlo,xsz,y}.
    The pixels will be set to color {fc} with  opacity 1.
    
    The procedure assumes that there are {nt} data frames with {nc}
    channels, and that {val[it][ic]} is the value of channel
    {ic} in data frame {it} for {ic} in {0..nc-1} and {it} in {0..nt-1}.

    The procedure scans the samples {val[0..nt-1][ic_mark]} of the given
    channel {ic_mark} and identifies the maximal runs of consecutive frames
    where those samples are non-zero. For each run, spanning frames
    {it_ini} to {it_fin} inclusive, the procedure paints with color {fc}
    two tic marks at the positions corresponding to {it_ini} and
    {it_fin}, and a horizontal line of between those two positions,
    as per {neuromat_image_paint_time_range}.  
    
    The frame with index{it} is assumed to have time {it+0.5}.
    The timeline spans the time interval from {tlo=0.0} 
    to {thi = nt}. */

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

    for (uint32_t iy = 0;  iy < NY; iy++)
      { for (uint32_t ix = 0;  ix < NX; ix++)
          { double wxy = float_image_get_sample(msk, 0, ix, iy);
            double smp;
            if (wxy == 0)
              { smp = 0.0; }
            else
              { smp = 0.0;
                for (uint32_t ie = 0;  ie < ne; ie++) 
                  { double basi = float_image_get_sample(bas[ie], 0, ix, iy);
                    smp += val[ie] * basi;
                  }
                assert(! isnan(smp));
              }
            float_image_set_sample(img, 0, ix, iy, (float)smp);
          }
      }
  }

float_image_t *neuromat_eeg_image_electrodes_overlay
  ( int32_t ne, 
    r2_t pos[], 
    double drad, 
    double hwd,
    int32_t ie_spec, 
    frgb_t *fcfill, 
    frgb_t *fcdraw,
    int32_t NX,
    int32_t NY
  )
  { 
    int32_t NC = 4;
    float_image_t *img = float_image_new(NC, NX, NY);
    float_image_fill(img, 0.000f);
    if ((fcfill == NULL) && (fcdraw == NULL)) { return img; }
    bool_t round = TRUE; /* Round dots? */
    bool_t diagonal = FALSE; /* Diagonal crosses? (Not used.) */
    int32_t msub = 3; /* Subsampling order. */
    int32_t cop = 3; /* Index of opacity channel. */
    for (uint32_t ie = 0;  ie < ne; ie++) 
      { r2_t *posi = &(pos[ie]);
        double cx = posi->c[0];
        double cy = posi->c[1];
        frgb_t *cfill = (ie == ie_spec ? fcdraw : fcfill); /* Fill color, or NULL. */
        frgb_t *cdraw = (ie == ie_spec ? fcfill : fcdraw); /* Draw color, or NULL. */
        for (uint32_t c = 0;  c < NC; c++)
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
    return img;
  }

float_image_t *neuromat_eeg_image_make_time_tracks  
  ( int32_t nt, 
    int32_t nc, 
    double **val, 
    int32_t nm, 
    neuromat_eeg_marker_spec_t marker[], 
    int32_t hw,
    int32_t NX,
    int32_t *mkdots_xctrP,
    double *mkdots_radP,
    int32_t *track_xloP, 
    int32_t *track_xszP, 
    int32_t track_y[],
    int32_t *slider_hhP
  )
  {
    demand(hw >= 1, "invalid {hw}");

    /* Choose the horizontal start and end of the tracks: */
    int32_t xmrg = 4*hw;                /* Extra marging at left and right. */
    int32_t drad = 2*hw;                /* Dot radius. */
    int32_t xspc = 3*hw;                /* Spacing between dots and time tracks. */
    int32_t xdctr = xmrg + drad;        /* X coord of dots center. */
    int32_t xlo = xdctr + drad + xspc;  /* Left offset to start of track (pixels). */
    int32_t xsz = NX - 2*xlo;           /* Horizontal extent of tracks (pixels). */
    demand(xsz >= 4*hw, "timeline image width is too small");
    (*mkdots_xctrP) = xdctr;
    (*mkdots_radP) = (double)drad;
    (*track_xloP) = xlo;
    (*track_xszP) = xsz;
    
    /* Choose the vertical positions of the tracks and paint them: */
    int32_t ymrg = 4*hw;                  /* Extra marging at bottom and top. */
    int32_t slhh = 4*hw;                  /* Max extent of slider from top and bot timelines. */
    int32_t ythw = 2*hw;                  /* Track half-width including highlights and tic marks. */
    int32_t yspc = 2*hw;                  /* Extra intertrack spacing. */
    int32_t ytstep = ythw + yspc + ythw;  /* Y step between track axes. */
    int32_t NY = ymrg + slhh + ythw + (nm-1)*ytstep + ythw + slhh + ymrg; /* Total height of tracks (pixels) incl margins. */
    (*slider_hhP) = slhh;

    /* Allocate the image: */
    int32_t NC = 4;
    float_image_t *img = float_image_new(NC, NX, NY);
    float_image_fill(img, 0.000f);
    
    /* Paint the tracks and marker channel ranges: */
    frgb_t ftrack = (frgb_t){{ 0.500f, 0.500f, 0.500f }}; /* Color or track line. */
    int32_t y0 = ymrg + yspc + ythw;      /* Y coordinate of lowest timeline. */
    for (uint32_t im = 0;  im < nm; im++)
      { int32_t yk = y0 + im*ytstep;
        track_y[im] = yk; 
        neuromat_image_paint_time_track(img, hw, xlo, xsz, yk, ftrack);
        int32_t ic = marker[im].ic; /* EEG channel index. */
        if ((ic >= 0) && (ic < nc))
          { frgb_t fmark = marker[im].color; /* Color or tics and marks line. */
            neuromat_image_paint_marker_ranges(img, hw, nt, nc, val, ic, xlo, xsz, yk, fmark);
          }
      }
    return img;
  }

void neuromat_image_paint_marker_ranges
  ( float_image_t *img,
    int32_t hw, 
    int32_t nt, 
    int32_t nc,
    double **val,
    int32_t ic_mark, 
    int32_t xlo, 
    int32_t xsz,
    int32_t y, 
    frgb_t fc
  )
  { 
    int32_t NC = (int32_t)img->sz[0];
    demand(NC == 4, "result image should be RGBA");
    double tlo = 0.0; /* Time just before first frame. */
    double thi = (double)nt; /* Time just after last frame. */
    /* Scan samples of channel {ic}: */
    /* Pretend that marker channel is zero before and after all frames. */
    int32_t it_ini = -1;  /* Index of start frame of run, or -1 if not started yet. */
    int32_t it_fin = -1;  /* Index of end frame of run, or -1 if not started yet. */
    for (uint32_t it = 0;  it <= nt; it++)
      { double smp = (it >= nt ? 0.0 : val[it][ic_mark]);
        assert(! isnan(smp));
        if (smp != 0.0)
          { /* Start or extend a run: */
            if (it_ini < 0) { it_ini = it; }
            it_fin = it;
          }
        else if (it_ini >= 0)
          { /* End of run: */
            assert((it_ini >= 0) && (it_ini <= it_fin) && (it_fin < nt));
            double tini = it_ini + 0.5;
            double tfin = it_fin + 0.5;
            neuromat_image_paint_time_range_and_tics(img, hw, tlo, tini, tfin, thi, xlo, xsz, y, fc);
            it_ini = -1;
            it_fin = -1;
          }
      }
  }

void neuromat_eeg_image_paint_marker_dots
  ( float_image_t *img,
    int32_t nc, 
    double valt[],
    int32_t nm,
    neuromat_eeg_marker_spec_t marker[], 
    r2_t mkdot_ctr[],
    double mkdot_rad
  )
  { 
    int32_t NC = (int32_t)img->sz[0];
    demand(NC == 4, "image should be RGBA");

    float vdraw = NAN;   /* Outline color. (Not used) */
    float hwd = 0.0;     /* Half-width of drawing pen. */
    bool_t round = TRUE; /* Round dots? */
    bool_t diagonal = FALSE; /* Diagonal crosses? (Not used.) */
    int32_t msub = 4; /* Subsampling order. */
    for (uint32_t im = 0;  im < nm; im++)
      { /* Get the marker's value {vr} scaled and clipped to {[-1..+1]} */
        int32_t ic = marker[im].ic; /* Index of channel in data frames: */
        double vm = ((ic >= 0) && (ic < nc) ? valt[ic] : 0.0);
        double vr = fmax(-1.0, fmin(+1.0, vm/(marker[ic].vref + 1.0e-100)));
        /* Compute the dot's color {fc}: */
        frgb_t fc;
        if (fabs(vr) < 1.0e-6)
          { fc = (frgb_t){{ 0.000f, 0.000f, 0.000f }}; }
        else
          { /* Compute max saturated color {fc}: */
            frgb_t fc = marker[im].color;
            if (vr < 0)
              { /* Complement {fc} relative to its saturation range: */
                double clo = fmin(fc.c[0], fmin(fc.c[1], fc.c[2]));
                double chi = fmax(fc.c[0], fmax(fc.c[1], fc.c[2]));
                for (uint32_t ia = 0;  ia < 3; ia++)  { fc.c[ia] = (float)(clo + (chi - fc.c[ia])/(chi - clo + 1.0e-100)); }
              }
            /* Scale {fc} (non-linearly) by {fabs(vr)}: */
            double vs = sqrt(fabs(vr));
            for (uint32_t ia = 0;  ia < 3; ia++)  { fc.c[ia] = (float)(vs*fc.c[ia]); }
          }
         /* Now paint the dot: */
         r2_t *ctrk = &(mkdot_ctr[im]);
         for (uint32_t c = 0;  c < NC; c++)
          { double rkc = (c < 3 ? mkdot_rad + 1.5 : mkdot_rad);
            float vfill = (float)(c < 3 ? fc.c[c] : 1.000);
            (void)float_image_paint_dot
              ( img, c, ctrk->c[0], ctrk->c[1], rkc, hwd, 
                round, diagonal, vfill, vdraw, msub
              );
          }
      }
  }

 
