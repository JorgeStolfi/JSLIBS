#ifndef neuromat_eeg_image_H
#define neuromat_eeg_image_H

/* NeuroMat EEG-specific image tools. */
/* Last edited on 2013-12-06 02:48:25 by stolfilocal */

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

void neuromat_eeg_image_paint_marker_tics
  ( float_image_t *img, 
    int bxmin, int bxmax, 
    int mymin, int mymax, 
    int nt,
    int nc, 
    double **val, 
    int ic_mark
  );
  /* Assumes that that
   {val[0..nt-1][0..nc-1]} is a measurement set with {nc} channels
   sampled at {nt} equally spaced sampling times.
   
   Paints vertical tic marks into {img} at the positions that
   correesponds to the indices {it} when the marker channel
   {val[it][ic_mark]} is nonzero.
   
   The tics will span pixel rows {mymin..mymax}
   inclusive, will lie within columns {bxmin..bxmax} inclusive.
   Specifically, each tic mark is centered at the fractional
   abscissa {neuromat_image_slider_bar_position(bxmin, bxmax, it+0.5,
   nt)}. Tic marks are at least one pixel wide, and closely spaced tic
   marks (less than 1 pixel apart) may be coalesced. */

#endif
