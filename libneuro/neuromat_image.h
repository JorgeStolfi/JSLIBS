#ifndef neuromat_image_H
#define neuromat_image_H

/* NeuroMat generic image tools. */
/* Last edited on 2017-09-29 01:23:24 by jstolfi */

#define _GNU_SOURCE
#include <frgb.h>
#include <float_image.h>
#include <r2.h>

#include <neuromat_eeg.h>

/* SHOWING TIMES AND TIME RANGES ON A TIMELINE BAR

  The procedures below draw times and time ranges with color {fc} into
  an image {img}, along a horizontal timeline.
  
  The timeline is an imaginary zero-width line segment located {y}
  pixels above the bottom edge of the image. It starts {xlo} pixels to
  the right of the left image border, and extends for exactly {xsz}
  pixels.
  
  Horizontal lines are drawn vertically centered on that timeline, with
  vertical width {2*hw} pixels.
  
  Each tic mark is {2*hw} pixels wide and {4*hw} pixels tall,
  vertically centered on the timeline. 
  
  The timeline represents some time interval {tlo__thi}. An arbitrary
  time {t} is mapped affinely from that interval to a real displacement
  along the timeline, in the interval {0__xsz}.  Thus the time {t = tlo}
  corresponds to the point {xlo} pixels to the right of the image
  border, and the time span {thi-tlo} corresponds to exactly {xsz}
  pixels. */

void neuromat_image_paint_time_track
  ( float_image_t *img, 
    int hw, 
    int xlo, 
    int xsz, 
    int y, 
    frgb_t *fc
  );
  /* Paints into {img} a horizontal stripe with color {fc},
    spanning the whole length of the timeline ({xsz} pixels). */

void neuromat_image_paint_marker_ranges
  ( float_image_t *img, 
    int hw, 
    int nt, 
    int nc,
    double **val,
    int ic_mark, 
    int xlo, 
    int xsz,
    int y,
    frgb_t *fc
  );
  /* Paints into image {img} tic marks and highlights showing the
    ranges of frames where the marker channel {ic} is nonzero.
    
    The procedure assumes that there are {nt} data frames with {nc}
    channels, and that {val[it][ic]} is the value of channel
    {ic} in data frame {it} for {ic} in {0..nc-1} and {it} in {0..nt-1}.

    The procedure scans the samples {val[0..nt-1][ic_mark]} of the given
    channel {ic_mark} and identifies the maximal runs of consecutive frames
    where those samples are non-zero. For each run, spanning frames
    {it_ini} to {it_fin} inclusive, the procedure paints with color {fc}
    two tic marks at the positions corresponding to {it_ini} and
    {it_fin}, and horizontal line of between those two positions,
    as per {neuromat_image_paint_time_range}. 
    
    The frame with index{it} is assumed to have time {it+0.5}.
    The timeline spans the time interval from {tlo=0.0} 
    to {thi = nt}.

 */

void neuromat_image_paint_time_range_and_tics
  ( float_image_t *img, 
    int hw, 
    double tlo,
    double tini,
    double tfin,
    double thi, 
    int xlo, 
    int xsz,
    int y,
    frgb_t *fc
  );
  /* Paints into image {img} with color {fc} two tic marks and a
    horizontal line over the timeline, showing the time range {tini _
    tfin} in the global time range {tlo__thi}.
    
    See {neuromat_image_paint_time_range} and 
    {neuromat_image_paint_tic}. */
    
void neuromat_image_paint_time_range
  ( float_image_t *img, 
    int hw, 
    double tlo,
    double tini,
    double tfin,
    double thi, 
    int xlo, 
    int xsz,
    int y,
    frgb_t *fc
  );
  /* Paints a horizontal straight line stroke with color {fc} over the
    timeline, showing the range {tini__tfin} in the global time range
    {tlo__thi}.

    The stroke ends are clipped to the timeline span.  If {tini>=tfin},
    or {tfin-tini} is too small, nothing is drawn.  */

void neuromat_image_paint_tic
  ( float_image_t *img, 
    int hw, 
    double tlo,
    double t,
    double thi, 
    int xlo, 
    int xsz,
    int y,
    frgb_t *fc
  );
  /* Paints a vertical tic mark with color {fc} showing
    the position of time {t} in the global time range {tlo__thi}.  

    If {t<tlo} or {t>thi}, the tic mark is omitted. */

void neuromat_image_paint_slider
  ( float_image_t *img, 
    int hw,
    double tlo,
    double t,
    double thi, 
    int xlo, 
    int xsz,
    int ylo,
    int yhi,
    frgb_t *fc
  );
 /* Paints a slider with color {fc} showing
    the position of time {t} in the global time range {tlo__thi}. 
    Fails if {t<tlo} or {t>thi}. 
    
    The slider will span the timelines with ordinates {y=ylo}
    and {y=yhi}.  I will extend at least {6*hw} below and above
    those lines. */

void neuromat_image_colorize
  ( float_image_t *fim, 
    float_image_t *msk, 
    double vmax, 
    int style, 
    frgb_t *bgr,
    float_image_t *cim
  );
  /* Stores into {cim} a false-color version of image {fim}, where the
    range {[-vmax__+vmax]} is turned into various colors with
    {frgb_path_signed(z,1,style)}. Pixels where the image is {NAN} are
    mapped to {bgr}.
    
    If {msk} is not null, it is assumed to be a mask image with samples
    in {[0__1]}, that has been multiplied into {fim}. The colorized
    image is then masked by {msk}, so that pixels where {msk} is zero
    are painted with color {bgr}. */

void neuromat_image_paint_overlay(float_image_t *img, float_image_t *ovr, float col[]);
  /* Paints the overlay image {ovr} over the image {img}, colorizing it with {col}.
    
    The image {img} may have any number {nc} of channels.  The image {ovr} must have a single
    channel and the same width and height as {img}.  Each pixel of {ovr} must be in the 
    range {[0__1]}, and is interpreted as the fraction of the pixel that is covered by the
    overlay.  If it is 0, the corresponding pixel of {img} is not affected.  If it is 1,
    the pixel is replaced by {col[0..nc-1]}.  For intermediate values, the pixel of {img}
    is blended with {col[0..nc-1]}. */

#endif
