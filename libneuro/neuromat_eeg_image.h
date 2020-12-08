#ifndef neuromat_eeg_image_H
#define neuromat_eeg_image_H

/* NeuroMat EEG-specific image tools. */
/* Last edited on 2017-09-29 01:25:55 by jstolfi */

#define _GNU_SOURCE
#include <r2.h>
#include <float_image.h>

#include <neuromat_eeg.h>

/* ANIMATION */

float_image_t *neuromat_eeg_image_schematic_head_mask(int NX, int NY, r2_t *ctr, r2_t *rad);
  /* Returns a monochromatic mask for a schematic head, namely 1 inside
    the ellipse with center {ctr} and radius {rad},0 outside, with
    antialiased boders. */

void neuromat_eeg_image_draw_electrodes
  ( float_image_t *img,
    int ne, 
    r2_t pos[], 
    double drad, 
    int ilo, 
    float vlo[], 
    float vhi[]
  );
  /* Draws dots on {img} at the electrode positions {pos[0..ne-1]}. Each
    dot is an antialiased disk with radius {drad} pixels, without
    outline. The dot with index {ilo} is painted with color
    {vlo[0..NC-1]}, all the other dots are painted with color
    {vhi[0..NC-1]}. */

void neuromat_eeg_image_paint_potentials
  ( int ne, 
    double val[], 
    float_image_t *bas[],
    float_image_t *msk, 
    int c, float_image_t *img
  );
  /* Paints the values {val[0..ne-1]} over channel {c} of the image {img}.
    In each pixel, computes the linear combination of {bas[0..ne-1]} with 
    coefficients {val[0..ne-1]}, then paints the result at the corresponding pixel
    of {img}, with opacity given by {msk}. */

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
  );
  /* Creates an RGB image with {NX} pixel columns and {NY} pixel rows,
    containing a set of timeline bars for {nt} data frames with {nc}
    channels, {nm} of them being marker channels.
    
    The procedure assumes that {val[it][ic]} is the value of channel
    {ic} in data frame {it} for {ic} in {0..nc-1} and {it} in {0..nt-1}.

    If {nm} is zero, the image will contain a single gray line (/track/)
    spanning most of the width of the image.
    
    If {nm} is positive, the image it will contain {nm} tracks, one for
    each marker channel. In the latter case, each track will have tics
    and highlights showing where the corresponding marker channel is
    nonzero. The procedure assumes that the channel index of marker {im}
    is {ic_mark[im]}, for {im} in {0..nm-1}.
    
    In any case, each track will be {2*hw} pixels wide, with {hw}
    additional pixels all around reserved for tic marks, and another
    {hw} pixels between adjacent tracks.  
    
    The highlights and tics on each marker channel 
    {ic=ic_mark[im]} are painted with {neuromat_image_paint_marker_ranges}. 
    
    The procedure returns in {*track_xloP} and {*track_xszP} the 
    horizontal offset and total width of the tracks, in pixels.
    Note that tic marks may extend {hw} pixels to the left or to the
    right of this range.
    
    The procedure also returns in {track_y[k]} the vertical offset from image
    edge to the centerline of track number {k}.  Note that 
    {track_y} must have {max(1,nm)} elements. */

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
  );
  /* Paints the marker channels into the given {frame}. 
  
    Assumes that
    {valt[0..nc-1]} are the sample values of the corresponding data
    frame; that there are {nm} markers to be painted whose channel
    indices are {ic_mark[0..nm-1]}.   
    
    The value {valt[ic]} of each marker channel {ic = ic_mark[im]}, for
    {im} in {0..nm}, is painted as a dot with center {ctr_mark[k]} and
    radius {rad_mark}. The color of the dot is obtained by dividing
    {valt[ic]} by {vscale[ic]}, then converting the result to an RGB
    color according to the color scale of the given {style}. */

#endif
