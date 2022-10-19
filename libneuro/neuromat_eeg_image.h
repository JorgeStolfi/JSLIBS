#ifndef neuromat_eeg_image_H
#define neuromat_eeg_image_H

/* NeuroMat EEG-specific image tools. */
/* Last edited on 2021-08-28 21:10:44 by stolfi */

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
  );
  /* Returns an RGBA image {img} with a dot at each electrode position {pos[0..ne-1]}.
    The positions are assumed to be in image coordinates (in pixels, from the lower
    left corner).
    
    The image {img} will have 4 channels (RGB + opacity). Each dot is an
    antialiased disk with radius {drad} pixels, filled with {fcfill} (if
    not {NULL}) with outline drawn with a pen of width {2*hwd} with
    value {fcdraw} (if not {NULL}). The opacity channel will be set
    accordingly.
    
    If {ie_spec} is in {0..ne-1}, the color schema is reversed for the dot
    with index {ie_spec}. */


typedef struct neuromat_eeg_marker_spec_t 
  { int32_t ic;   /* Channel index. */
    double vref;   /* Reference value. */
    frgb_t color;  /* Associated color. */
  } neuromat_eeg_marker_spec_t;
  /* Specification of a marker channel, for {neuromat_eeg_image_make_time_tracks_image}. */


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
  );
  /* Returns an RGBA image {img} containing a set of {nm >= 1} /time tracks/ 
    for marker channels in {val}. Each track is a gray line  
    spanning most of the width of the image.
    
    The procedure assumes that {val} has {nt} data frames with {nc}
    channels. Specifically that {val[it][ic]} is the value of channel
    {ic} in data frame {it} for {ic} in {0..nc-1} and {it} in {0..nt-1}. 
    
    The procedure assumes that time track {im} refers to a marker
    channel with specs {marker[im]}, for {im} in {0..nm-1}. The track
    will have tics and highlights showing where channel
    {ic[im]=marker[im].ic} is nonzero. If {ic[im]} is not in
    {0..nc-1}, the track is just a gray bar without tics or highlights.
    
    The image will have {NX} pixel columns. The procedure assumes that a
    colored dot will be drawn to the left of each time track, showing
    the relative value of the corresponding marker for the current
    frame. The procedure returns in {*mkdots_xctrP} the X coordinate of
    the dots' centers, and in {*mkdots_radP} their radius.
    
    Each track will be {2*hw} pixels wide, with {hw} additional pixels
    all around reserved for tic marks, and another {hw} pixels between
    adjacent tracks. The highlights and tics of track {im}, if any, are
    painted with {neuromat_image_paint_marker_ranges} using color
    {marker[im].color}.
    
    The procedure also assumes that a slider will be drawn over the 
    time tracks, showing the position of the frame in the dataset. 
    It returns in  {*slider_hhP} the allowed extent of the slider
    above the top timeline and below the bottom timeline.
    
    The height of the image will be just enough to contain
    the tracks (with the tics) plus space above an below for the slider,
    and some safety margin.
    
    The procedure returns in {*track_xloP} and {*track_xszP} the 
    horizontal offset and total width of the tracks, in pixels.
    Note that tic marks may extend {hw} pixels to the left or to the
    right of this range.
    
    The array {track_y} must have {nm} elements. The procedure also
    returns in {track_y[im]} the vertical offset from the bottom image
    edge to the centerline of track number {im}. Note also tic marks and
    other elements may extend several times {hw} below {track_y[0]} and
    above {track_y[nm-1]}. */

void neuromat_eeg_image_paint_marker_dots
  ( float_image_t *img,
    int32_t nc, 
    double valt[],
    int32_t nm,
    neuromat_eeg_marker_spec_t marker[], 
    r2_t mkdot_ctr[],
    double mkdot_rad
  );
  /* Paints onto the RGBA image {img} a set of {nm} dots showing
    the state of marker channels in a data frame {valt}.
  
    Assumes that {valt[0..nc-1]} are the sample values of the
    corresponding data frame; that there are {nm} markers to be painted
    specified by {marker[0..nm-1]}.
    
    The value {vm=valt[ic]} of each marker channel {ic = marker[im].ic},
    for {im} in {0..nm}, is painted as a dot with center {mkdot_ctr[im]}
    and radius {mkdot_rad}. Let {vr} be the ratio
    {vm/marker[im].vref}, clipped to {[-1__+1]}. If {vr} is
    positive, the dot is painted with the color {vr*fcp} where
    {fcp=marker[im].color}. If {vr} is negative, it is painted with color
    {-vr*fcn} where {fcn} is the complementary color to {fcp}. The
    opacity channel of {img} is set approriately.
    
    As a special case, if {marker[im].ic} is not in {0..nc-1}, then
    the marker value {vm} is assumed to be zero, so that the dot will be always
    black. */

#endif
