#ifndef msm_image_paint_H
#define msm_image_paint_H

/* Procedures that paint various things on images. */
/* Last edited on 2013-10-22 07:01:08 by stolfilocal */

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
  correspond to index positions in two sequences, {s0} and {s1}.
  The virtual image dimensions are {img->NV[0]} and {img->NV[1]}.
  
  Virtual image coordinates are always implicitly reduced 
  modulo these two dimensions.
  
  For memory size reasons, the virtual image {img} is represented by a
  float image {img->fim}, whose size is that of {img} reduced by a
  positive integer scale factor {img->scale}; so that each image pixel
  corresponds to {scale^2} distinct index pairs (rungs).  More
  precisely, a generic rung {(i0,i1)} maps to the pixel {(xp,yp) =
  ((i0-i0Min)/img->scale,(i1-i1Min)/img->scale)} of {img->fim}, always
  rounded downwards.

  !!! Remove the need for {i0Min,i1Min}, by trimming the seqs or including in {msm_image_t} !!!

  */

/* RUNG SCORES */

void msm_image_seq_seq_score_paint
  ( msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1, 
    msm_rung_score_proc_t *score, 
    msm_image_t *img,
    int i0Min, 
    int i1Min
  );
  /* Sets each pixel in column {i0/scale} and row {i1/scale} of
    {img} to the value of {score(seq0[i0],seq1[i1])}, where
    {score} is a user-defined function.
    
    The {scale} must be positive. If {scale} is
    greaer than 1, all scores that map to the same pixel are averaged
    together. */

/* PAIRINGS */

void msm_image_pairing_paint
  ( msm_pairing_t *pr,
    int n0,
    int n1,
    float v[],
    bool_t add,
    msm_image_t *img,
    int i0Min, 
    int i1Min
  );
  /* Paints the pairing {pr}, reduced by {1/scale}, on image {img}.
    The painting uses color {v[0..NC]} where {NC} is the number of
    channels in {img}.
    
    More precisely, reduces each rung {i0,i1} of {pr} modulo the
    parameters {n0,n1}, to the rectangle {{0..n0-1}×{0..n1-1}}; and
    then assigns {v} (if {add} is FALSE) or adds {v/scale} (if {add}
    is TRUE) to the pixel in column {i0/scale} and row {i1/scale} of
    {img}. */

/* CANDIDATES */

void msm_image_cand_paint
  ( msm_cand_t *cd, 
    float v[], 
    bool_t add, 
    msm_image_t *img, 
    int i0Min, 
    int i1Min
  );
  /* Paints the pairing of candidate {cd}, reduced by {1/scale}, on
    image {img}. The pairing is painted with pixel value {v[0..NC-1]}
    where {NC} is the number of chanels in {img}.
    
    !!! Explain/check the use of {i0Min,i1Min} !!!  
    
    Equivalent to {msm_image_pairing_paint(cd->pr,n0,n1,v,add,img,scale)}
    where {n0} and {n1} are the lengths of the two sequences paired by {cd}. */

/* AND CANDIDATE LISTS */

void msm_image_cand_vec_paint
  ( msm_cand_vec_t *cdv,
    msm_seq_desc_t *seq0, 
    msm_seq_desc_t *seq1,
    msm_image_t *img,
    int i0Min, 
    int i1Min
  );
  /* Fills {img} with an image showing the candidates of {cdv},
    against a background of {-INF} values.
    
    !!! Explain/check the use of {i0Min,i1Min} !!!
    
    All candidates must refer to the same two sequences {seq0,seq1}. The
    color {v} of each candidate {cd} is its score {cd->score} divided by
    the number of steps. The image is similar to that
    produced by {msm_image_cand_paint}, summed over all candidates {cd}
    in {cdv}. */

/* DYNAMIC PROGRAMMING TABLEAUS */

void msm_image_dyn_tableau_scores_paint
  ( int n0,                /* Length of first sequence. */
    int n1,                /* Length of second sequence. */
    msm_dyn_tableau_t *tb, /* Tableau. */
    msm_image_t *img,
    int i0Min, 
    int i1Min
  );
  /* Fills {img} with an image that shows the scores associated with
    each entry of the dynamic programming tableau {tb}. 
    The image should be monochromatic (single-channel).
    
    !!! Explain/check the use of {i0Min,i1Min} !!!
    
    Assumes that {tb} refers to two sequences with {n0} and {n1}
    samples, respectively. Assumes also that the pixel in column {i0}
    and row {i1} refer to samples {i0/scale} and {i1/scale} of the two
    sequences. If {(i0,i1)} lies outside the tableau, the pixel is set
    to {-INF}. The indices are taken modulo {n0} and {n1} (i.e., the
    sequences are assumed to be circular.) */

void msm_image_dyn_tableau_pairing_paint
  ( int n0,                /* Length of first sequence. */
    int n1,                /* Length of second sequence. */
    msm_dyn_tableau_t *tb, /* Tableau. */
    msm_rung_t gopt,       /* Optimal end-rung or {msm_rung_none} */ 
    float clr[],
    msm_image_t *img,
    int i0Min, 
    int i1Min
  );
  /* Stores {clr[0..3]} into the pixels of {img} that correspond 
    to the rungs of the optimal path in tableau {tb} that ends with
    rung {gopt}. 
    
    !!! Explain/check the use of {i0Min,i1Min} !!!
    
    Assumes that the tableau refers to two sequences with {n0} and
    {n1} samples, respectively. Assumes also that the pixel in column
    {i0} and row {i1} refer to samples {i0/scale} and {i1/scale} of
    the two sequences. The indices are taken modulo {n0} and {n1}
    (i.e., the sequences are assumed to be circular.) */

void msm_image_grid_paint
  ( int step0, 
    int step1, 
    float clr[],
    msm_image_t *img,
    int i0Min, 
    int i1Min
  );
  /* Draws a vertical line at every X-coordinate that is a multiple of
    {n0} and a horizontal line at every Y-coordinate that is a multiple
    of {n1}.
    
    !!! was {msm_image_seq_periods_paint} !!!
    
    !!! Explain/check the use of {i0Min,i1Min} !!!
    
    The lines are one pixel wide, with color {clr}. 
    The image dimensions must be multiples of {n0}
    and {n1}. */ 

#endif

