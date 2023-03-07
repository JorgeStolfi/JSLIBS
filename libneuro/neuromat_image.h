#ifndef neuromat_image_H
#define neuromat_image_H

/* NeuroMat generic image tools. */
/* Last edited on 2023-03-07 17:17:49 by stolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <sign.h>
#include <frgb.h>
#include <float_image.h>
#include <r2.h>

#include <neuromat_eeg.h>

/* COLORIZING MONOCHROMATIC IMAGES */

void neuromat_image_colorize_field
  ( float_image_t *cim,
    float_image_t *fld, 
    float_image_t *msk, 
    double vmax, 
    int32_t style
  );
  /* Stores into {cim} a false-color version of image {fld}, where the
    value range {[-vmax__+vmax]} is turned into various colors with
    {frgb_path_map_signed(z,1,style)}. 
    
    The image {cim} must have 4 chanels (RGB plus opacity). If {msk} is not null, it
    must be be a monochromatic image with samples in {[0__1]}, no
    {NAN}s. If {msk} is null, it is assumed to be an image with all
    pixels equal to 1. 
    
    Pixels of {fld} outside that range are clipped to that range. Pixels
    where {fld} is {NAN} or where the mask image {msk} is zero are
    set to black with opacity zero. Pixels where the mask is nonzero
    are set to the RGB color with that opacity.
    
    */

/* COLORIZING AN IMAGE OVERLAY */

void neuromat_image_colorize_signed_overlay
  ( float_image_t *cim,
    float_image_t *ovr, 
    sign_t sgn,
    frgb_t fc
  );
  /* Converts a signed monochromatic overlat {ovr} into an RGBA image {cim}.
    
    The image {cim} must have 4 channels. The image {ovr} must have a
    single channel and the same width and height as {cim}.  The sign {sgn}
    should be either {+1} or {-1}.
    
    Each pixel of {ovr} must be {NAN} or in the range {[-1__+1]}. If its
    is not {NAN} and its sign is {sgn}, channels {0..2} of the corresponding pixel of {cim}
    are set to the color {fc}, and channel 3 (opacity) is set to the absolute value of {ovr}.
    Otherwise the pixel of {cim} is set to all zeros (black, transparent). */

/* SHOWING TIMES AND TIME RANGES ON A TIMELINE BAR

  The procedures below draw time tracks, ranges, tics, or sliders with sample value {smp} into
  an RGBA image {img}.
  
  !!! Tic marks and the slider should be antialiased !!!
  
  The range highlights and tic marks are vertically centered on the
  corresponding /time line/, the axis of the coresponding time track.
  The timeline is located {y} pixels above the bottom of the image,
  starts {xlo} pixels to the right of the left image border, and extends
  for exactly {xsz} pixels.
  
  Horizontal lines are drawn vertically centered on that timeline, with
  vertical width {2*hw} pixels.
  
  Each tic mark is {2*hw} pixels wide and {4*hw} pixels tall,
  vertically centered on the timeline. 
  
  The timeline represents some time interval {tlo__thi}. An arbitrary
  time {t} is mapped affinely from that interval to a real displacement
  along the timeline, in the interval {0__xsz}, and then rounded to 
  the nearest integer coordinate.  Thus the time {t = tlo}
  corresponds to the point {xlo} pixels to the right of the image
  border, and the time span {thi-tlo} corresponds to exactly {xsz}
  pixels.
  */

void neuromat_image_paint_time_track
  ( float_image_t *img,
    int32_t hw, 
    int32_t xlo, 
    int32_t xsz, 
    int32_t y, 
    frgb_t fc
  );
  /* Paints onto {img} a /time track/, as a horizontal stripe of of
    half-width {hw} centered on the timeline. The track starts {xlo}
    pixels to the right of the left image border, and extends for
    exactly {xsz} pixels. The pixels will be set to color {fc} with
    opacity 1. */

void neuromat_image_paint_time_range_and_tics
  ( float_image_t *img,
    int32_t hw, 
    double tlo,
    double tini,
    double tfin,
    double thi, 
    int32_t xlo, 
    int32_t xsz,
    int32_t y, 
    frgb_t fc
  );
  /* Paints onto image {img} with color {fc} two tic marks and a
    horizontal line over the timeline, showing the time range
    {tini__tfin} in the global time range {tlo__thi}. The highlights and
    tics are located on the timeline defined by {xlo,xsz,y}.
    
    See {neuromat_image_paint_time_range} and 
    {neuromat_image_paint_tic}. */
    
void neuromat_image_paint_time_range
  ( float_image_t *img,
    int32_t hw, 
    double tlo,
    double tini,
    double tfin,
    double thi, 
    int32_t xlo, 
    int32_t xsz,
    int32_t y, 
    frgb_t fc
  );
  /* Paints a horizontal straight line stroke with color {fc} over the
    timeline, showing the range {tini__tfin} in the global time range
    {tlo__thi}.

    The stroke ends are clipped to the timeline span.  If {tini>=tfin},
    or {tfin-tini} is too small, nothing is drawn.  */

void neuromat_image_paint_tic
  ( float_image_t *img,
    int32_t hw, 
    double tlo,
    double t,
    double thi, 
    int32_t xlo, 
    int32_t xsz,
    int32_t y, 
    frgb_t fc
  );
  /* Paints a vertical tic mark with color {fc} showing
    the position of time {t} in the global time range {tlo__thi}.  

    If {t<tlo} or {t>thi}, the tic mark is omitted. */

void neuromat_image_paint_slider
  ( float_image_t *img,
    int32_t hw,
    int32_t hh,
    double tlo,
    double t,
    double thi, 
    int32_t xlo, 
    int32_t xsz,
    int32_t ylo,
    int32_t yhi, 
    frgb_t fc
  );
 /* Paints a slider with color {fc} showing
    the position of time {t} in the global time range {tlo__thi}. 
    Fails if {t<tlo} or {t>thi}. 
    
    The slider will span the timelines with ordinates {y=ylo}
    and {y=yhi}.  I will extend at least {hh} below and above
    those lines. */

#endif
