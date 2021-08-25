#ifndef neuromat_eeg_image_H
#define neuromat_eeg_image_H

/* NeuroMat EEG-specific image tools. */
/* Last edited on 2021-08-25 02:25:37 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <r2.h>
#include <frgb.h>
#include <float_image.h>

#include <neuromat_eeg.h>

float_image_t *neuromat_eeg_image_make_idealized_scalp_mask
  ( int32_t NX, 
    int32_t NY, 
    r2_t *ctr, 
    r2_t *rad
  );
  /* Returns a monochromatic mask for the upper half of the idealized 3D scalp on the
    idealized 2D scalp plane, namely 1 inside  the ellipse with center {ctr} and radius {rad},
    0 outside, with antialiased boders. */

void neuromat_eeg_image_compute_pot_field
  ( int32_t ne, 
    double val[], 
    float_image_t *bas[],
    float_image_t *msk, 
    float_image_t *img
  );
  /* Sets the single-channel image {img} to the linear combination of the single-chanel images
    {bas[0..ne-1]} with coefficients {val[0..ne-1]}.  Only sets pixels 
    where {msk} is non-zero; where {msk} is zero, stores 0 instead. */

void neuromat_eeg_image_paint_electrodes
  ( int32_t ne, 
    r2_t pos[], 
    double drad, 
    double hwd,
    int32_t ie_spec, 
    frgb_t *fcfill, 
    frgb_t *fcdraw,
    float_image_t *img
  );
  /* Draws dots onto image of image {img} at the electrode positions {pos[0..ne-1]}.
    The positions are assumed to be in image coordinates (in pixels, from the lower
    left corner).
    
    The image {img} must have 4 channels (RGB + opacity). Each dot is an antialiased disk with radius {drad} pixels, filled
    with {fcfill} (if not {NULL}) with outline drawn with a pen of width {2*hwd} with
    value {fcdraw} (if not {NULL}). The opacity channel will be set accordingly.
    
    If {ie_spec} is in {0..ne-1}, the color schema is reversed for the dot
    with index {ie_spec}. */

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
  );
  /* Paints onto the RGBA image {img} a set of {nk} /time tracks/ for {nm} marker channels in {val}.
    Each track is a gray line  spanning most of the width of the image.
    
    The procedure assumes that {val} has {nt} data frames with {nc}
    channels. Specifically that {val[it][ic]} is the value of channel
    {ic} in data frame {it} for {ic} in {0..nc-1} and {it} in {0..nt-1},
    of which {nm} are assumed to be the marker ones. The procedure
    assumes that the channel index of marker {im} is {ic_mark[im]}, for
    {im} in {0..nm-1}.

    If {nm} is zero, the image will contain {nk=1} tracks.
    
    If {nm} is positive, the image it will contain {nk=nm} tracks, one for
    each marker channel. In the latter case, each track will have tics
    and highlights showing where the corresponding marker channel is
    nonzero.
    
    In any case, each track will be {2*hw} pixels wide, with {hw}
    additional pixels all around reserved for tic marks, and another
    {hw} pixels between adjacent tracks.  
    
    The highlights and tics on each marker channel 
    {ic=ic_mark[im]} are painted with {neuromat_image_paint_marker_ranges}
    using color {fc[im]}. 
    
    The procedure returns in {*track_xloP} and {*track_xszP} the 
    horizontal offset and total width of the tracks, in pixels.
    Note that tic marks may extend {hw} pixels to the left or to the
    right of this range.
    
    The procedure also returns in {track_y[k]} the vertical offset from image
    edge to the centerline of track number {k}.  Note that 
    {track_y} must have {nk = max(1,nm)} elements.  Note also that tic marks and 
    other elements may extend several times {hw} below {track_y[0]} and */

void neuromat_eeg_image_paint_marker_dots
  ( int32_t nc, 
    double valt[],
    int32_t nm,
    int32_t ic_mark[], 
    double vscale[], 
    frgb_t fpos[],
    frgb_t fneg[],
    r2_t ctr_mark[],
    double rad_mark,
    float_image_t *img
  );
  /* Paints onto the RGBA image {img} a set of {nm} dots showing
    the state of marker channels in a data frame {valt}.
  
    Assumes that {valt[0..nc-1]} are the sample values of the
    corresponding data frame; that there are {nm} markers to be painted
    whose channel indices are {ic_mark[0..nm-1]}.
    
    The value {valt[ic]} of each marker channel {ic = ic_mark[im]}, for
    {im} in {0..nm}, is painted as a dot with center {ctr_mark[k]} and
    radius {rad_mark}. The color of the dot is obtained by dividing
    {valt[ic]} by {vscale[ic]}, then converting the result to an RGB
    color according to the color scale of the given {style}.  The 
    opacity channel of {img} is set approriately. */

#endif
