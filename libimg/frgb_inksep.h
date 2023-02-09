#ifndef frgb_inksep_H
#define frgb_inksep_H

/* frgb_inksep.h - decomposing a color into a mixture of inks
** Last edited on 2007-01-03 10:37:23 by stolfi
**
** Copyright (C) 2006 by Jorge Stolfi, the University of Campinas, Brazil.
** See the rights and conditions notice at the end of this file.
*/

#include <r4x4.h>
#include <bool.h>
#include <frgb.h>

#define frgb_inksep_CHNS 3
  /* Number of input channels. */
  
#define frgb_inksep_LAYS (frgb_inksep_CHNS+1)
  /* Number of color layers is always 1 + number of input channels. */

void frgb_inksep_shrink_color_set(r4x4_t *mix, double shrink);
  /* Each row of the matrix {mix} should be the RGB coordinates of the
    color of one ink layer, with an homogenizing weight 1 appended.
    The procedure displaces each color {shrink} of the way
    towards the average of the four colors.  Thus {shrink=0}
    has no effect, and {srink=1} replaces all rows by the 
    average color. */

void frgb_inksep_separate_layers(frgb_t fv, r4x4_t *sep, r4_t *m);
  /* Analyzes the color {fv} into the affine mixing of {LAYS} ink layers,
    unsing the analysis matrix {sep} (the inverse of the mixing matrix);
    where {LAYS} is {frgb_inksep_LAYS}. Returns in {mix[0..LAYS-1]} 
    the mixing coefficients. */

void frgb_inksep_clip_to_simplex(r4_t *m, int *nbadP);
  /* Given an array {m[0..LAYS-1]} of ink mixing coefficients, adjusts
    it so that all elements are non-negative and add to 1.
    Also returns in {*nbadP} the number of mixing coefficients
    that were originally negative. */

frgb_t frgb_inksep_remap_color(r4_t *m, r4x4_t *syn);
  /* Re-synthesizes a color by mixing four new inks with
    mixing coefficients {m[0..3]}.  The first three columns of 
    {syn} must be the RGB coordinates of the new ink colors;
    the fourth column must be all 1's.  Beware that the
    result may lie outside the unit cube. */
    
void frgb_inksep_reveal_colors(r4_t *m);
  /* The input {m[0..LAYS-1]} should be a set of ink mixing
    coefficients (adding to 1.0).  
    
    The procedure assumes that the inks are semi-transparent and
    deposited in layers, with layer 0 in the background and layer
    {LAYS-1} at the foreground. It magnifies each mixing coefficient
    {m[i]} to compensate for the combined coverage of overlaid layers
    {m[i+1..LAYS-1]}; so that each {m[i]} becomes the mixing
    coefficient of layer {i} as it was just applied, rather than in
    the final color. In particular, {m[0]} is always set to 1.0.
    
    E.g., if {m} is originally {0.125,0.125,0.250,0.500},
    the procedure will adjust it to {1.000,0.500,0.500,0.500}. */

void frgb_inksep_compute_separations(frgb_t bgColor, r4_t *m, r4x4_t *mix, frgb_t gv[]);
  /* Computes color separations {gv[0..LAYS-1]}, where each {gv[i]} 
    is the result of painting row {i} of the {mix} array
    over the background color {bgColor}, with coverage fraction {m[i]}. */

#endif
