#ifndef neuromat_image_H
#define neuromat_image_H

/* NeuroMat generic image tools. */
/* Last edited on 2013-12-06 04:51:18 by stolfilocal */

#define _GNU_SOURCE
#include <frgb.h>
#include <float_image.h>

#include <neuromat_eeg.h>

/* TIME SCALE AND TRIGGERS */

float_image_t *neuromat_image_make_slider_bar(int NC, int NX, int NY, int *bxminP, int *bxmaxP);
  /* Creates an image with {NC} channels of a slider bar, namely a black rectangle with a grey horizontal line,
    2 pixel wide, spanning almost  its whole width.  The height {NY} must be even and at least 4
    pixels. Returns in {(*bxminP)} and {(*bxmaxP)} the leftmost and righmost pixel column indices
    of the grey bar. */

double neuromat_image_slider_bar_position(int bxmin, int bxmax, double t, int nt);
  /* Returns the fractional column index along the slider bar
    that corresponds to a particular time {t} in the range {[0 _ nt]}.
    Assumes that the slider bar spans pixel columns {bxmin..bxmax} inclusive.
    Namely, if {t=0} returns {bxmin}, and if {t=nt} returns {bxmax+1}.
    Other values of {t} are interpolated linearly. */
    
void neuromat_image_paint_slider_dot(float_image_t *img, int bxmin, int bxmax, double yctr, double t, int nt);
 /* Assumes that {img} contains slider bar spanning pixel columns
   {bxmin..bxmax} inclusive, vertically centered at ordinate {yctr}.
   Paints a dot over the slider bar, at an abscissa corresponding to
   time {t} in the range {[0 _ nt]} (see
   {neuromat_image_slider_bar_position}). */

float_image_t *neuromat_image_colorize(float_image_t *fim, float_image_t *msk, double vmax, int style, frgb_t *bgr);
  /* Returns a false-color version of image {fim}, where the range
    {[-vmax _ +vmax]} is turned into various colors with
    {frgb_path_signed(z,1,style)}. Pixels where the image is {NAN} are
    mapped to {bgr}.
    
    If {msk} is not null, it is assumed to be a mask image with samples 
    in {[0 _ 1]}, that has been multiplied into {fim}.
    The colorized image is then masked by {msk}, so that pixels
    where {msk} is zero are painted with color {bgr}. */

void neuromat_image_paint_overlay(float_image_t *fim, float_image_t *ovr, float col[]);
  /* Paints the overlay image {ovr} over the image {fim}, colorizing it with {col}.
    
    The image {fim} may have any number {nc} of channels.  The image {ovr} must have a single
    channel and the same width and height as {fim}.  Each pixel of {ovr} must be in the 
    range {[0 _ 1]}, and is interpreted as the fraction of the pixel that is covered by the
    overlay.  If it is 0, the corresponding pixel of {fim} is not affected.  If it is 1,
    the pixel is replaced by {col[0..nc-1]}.  For intermediate values, the pixel of {fim}
    is blended with {col[0..nc-1]}. */

#endif
