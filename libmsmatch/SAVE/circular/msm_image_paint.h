#ifndef msm_image_paint_H
#define msm_image_paint_H

/* Procedures that paint various things on images. */
/* Last edited on 2008-01-29 17:21:51 by hcgl */

#define msm_image_paint_H_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <msm_basic.h>
#include <msm_rung.h>
#include <msm_seq_desc.h>
#include <msm_pairing.h>
#include <msm_cand.h>
#include <msm_cand_vec.h>
#include <msm_image.h>
#include <msm_dyn.h>

#include <vec.h>
#include <float_image.h>

#include <stdint.h>

/* TOOLS FOR PAINTING RUNGS INTO IMAGES

  The procedures in this interface are used to paint into an
  image-like object {img} whose the X and Y pixel coordinates
  correspond to index positions in two sequences, {xs} and {ys}.
  The visrtual image dimensions are {img->NVX} and {img->NVY}.
  
  Virtual image coordinates are always implicitly reduced 
  modulo these two dimensions. Thus,
  for example, one can paint a pairing between two circular sequences
  in either the `folded' or `unfolded' version, by properly choosing the
  virtual image dimensions.
  
  For memory size reasons, the virtual image {img} is represented by a
  float image {img->fim}, whose size is that of {img} reduced by a
  positive integer scale factor {img->scale}; so that each image pixel
  corresponds to {scale^2} distinct index pairs (rungs). Also, the
  image may span multiple periods of a circular sequence. More
  precisely, a generic rung {(x,y)} maps to the pixel {(xp,yp) =
  ((x-xMin)/img->scale,(y-yMin)/img->scale)} of {img->fim}, always
  rounded downwards. */

/* RUNG SCORES */

void msm_image_seq_seq_score_paint
  ( msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp, 
    msm_rung_score_proc_t *score, 
    msm_image_t *img,
    int xMin, 
    int yMin
  );
  /* Sets each pixel in column {ix/scale} and row {iy/scale} of
    {img} to the value of {score(xp[ix],yp[iy])}, where
    {score} is a user-defined function.
    
    The {scale} must be positive. If {scale} is
    greaer than 1, all scores that map to the same pixel are averaged
    together. */

/* PAIRINGS */

void msm_image_pairing_paint
  ( msm_pairing_t *pr,
    int nx,
    int ny,
    float v[],
    bool_t add,
    msm_image_t *img,
    int xMin, 
    int yMin
  );
  /* Paints the pairing {pr}, reduced by {1/scale}, on image {img}.
    The painting uses color {v[0..NC]} where {NC} is the number of
    channels in {img}.
    
    More precisely, reduces each rung {ix,iy} of {pr} modulo the
    parameters {nx,ny}, to the rectangle {{0..nx-1}×{0..ny-1}}; and
    then assigns {v} (if {add} is FALSE) or adds {v/scale} (if {add}
    is TRUE) to the pixel in column {ix/scale} and row {iy/scale} of
    {img}. */

/* CANDIDATES */

void msm_image_cand_paint
  ( msm_cand_t *cd, 
    float v[], 
    bool_t add, 
    msm_image_t *img, 
    int xMin, 
    int yMin
  );
  /* Paints the pairing of candidate {cd}, reduced by {1/scale}, on
    image {img}. The pairing is painted with pixel value {v[0..NC-1]}
    where {NC} is the number of chanels in {img}.  
    
    Equivalent to {msm_image_pairing_paint(cd->pr,nx,ny,v,add,img,scale)}
    where {nx} and {ny} are the lengths of the two sequences paired by {cd}. */

/* AND CANDIDATE LISTS */

void msm_image_cand_vec_paint
  ( msm_cand_vec_t *cdv,
    msm_seq_desc_t *xp, 
    msm_seq_desc_t *yp,
    msm_image_t *img,
    int xMin, 
    int yMin
  );
  /* Fills {img} with an image showing the candidates of {cdv},
    against a background of {-INF} values. All candidates must refer
    to the same two sequences {xp,yp}. The color {v} of each candidate
    {cd} is its score {cd->score} divided by the number of fundamental
    steps. The image is similar to that produced by {msm_image_cand_paint},
    summed over all candidates {cd} in {cdv}. */

/* DYNAMIC PROGRAMMING TABLEAUS */

void msm_image_dyn_tableau_scores_paint
  ( int nx,                /* Length of first sequence. */
    int ny,                /* Length of second sequence. */
    msm_dyn_tableau_t *tb, /* Tableau. */
    msm_image_t *img,
    int xMin, 
    int yMin
  );
  /* Fills {img} with an image that shows the scores associated with
    each entry of the dynamic programming tableau {tb}. 
    The image should be monochromatic (single-channel).
    
    Assumes that {tb} refers to two sequences with {nx} and {ny}
    samples, respectively. Assumes also that the pixel in column {ix}
    and row {iy} refer to samples {ix/scale} and {iy/scale} of the two
    sequences. If {(ix,iy)} lies outside the tableau, the pixel is set
    to {-INF}. The indices are taken modulo {nx} and {ny} (i.e., the
    sequences are assumed to be circular.) */

void msm_image_dyn_tableau_pairing_paint
  ( int nx,                /* Length of first sequence. */
    int ny,                /* Length of second sequence. */
    msm_dyn_tableau_t *tb, /* Tableau. */
    msm_rung_t gopt,       /* Optimal end-rung or {msm_rung_none} */ 
    float clr[],
    msm_image_t *img,
    int xMin, 
    int yMin
  );
  /* Stores {clr[0..3]} into the pixels of {img} that correspond 
    to the rungs of the optimal path in tableau {tb} that ends with
    rung {gopt}. 
    
    Assumes that the tableau refers to two sequences with {nx} and
    {ny} samples, respectively. Assumes also that the pixel in column
    {ix} and row {iy} refer to samples {ix/scale} and {iy/scale} of
    the two sequences. The indices are taken modulo {nx} and {ny}
    (i.e., the sequences are assumed to be circular.) */

void msm_image_seq_periods_paint
  ( int nx,                /* Length of first sequence. */
    int ny,                /* Length of second sequence. */
    float clr[],
    msm_image_t *img,
    int xMin, 
    int yMin
  );
  /* Draws a single-pixel vertical line with color {clr} at every X
    coordinate that is a multiple of {nx}. Draws a single-pixel
    horizontal line with color {clr} at every Y coordinate that is a
    multiple of {ny}. The image dimensions must be multiples of {nx}
    and {ny}. */ 

#endif

