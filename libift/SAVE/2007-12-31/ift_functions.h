/* ift_functions.h - some path cost functions for segmentation */
/* Last edited on 2007-12-26 20:16:39 by stolfi */

#ifndef ift_functions_H
#define ift_functions_H

#include "ift.h"
#include <stdio.h>

/* FUNCTION SELECTION BY NAME */

PathCostFn *get_cost_function(char *name);
  /*
    Returns the function `pf_XXX' from the list below,
    where `XXX' is the given `name'. */
  
#define ift_functions_INFO \
  "        maxval\n" \
  "          Maximum intensity (L1 modulus) of pixels along the" \
  " path.  Useful for watershed transform.\n" \
  "\n" \
  "        sumval\n" \
  "          Sum of the intensities (L1 moduli) of pixels along" \
  " the path.  Useful for variable-speed Voronoi.\n" \
  "\n" \
  "        maxdiff\n" \
  "          Maximum absolute difference (L1 distance) of consecutive" \
  " pixels in the path.  Useful for gradient-based segmentation.\n" \
  "\n" \
  "        sumdiff\n" \
  "          The sum of the absolute differences (L1 distances) of" \
  " consecutive pixels along the path.  Useful for weighted" \
  " Euclidean transformation.\n" \
  "\n" \
  "        monoroot\n" \
  "          The intensity (L1 modulus) of the first pixel if the" \
  " intensities are non-decreasing along the path; else +oo.  Useful" \
  " to find local intensity minima in the image."
  
/* PATH COST FUNCTIONS */

PathCost pf_maxval(PixelNode *s, RelArc *a, PixelNode *t, ImageGraph *G, int stage);
  /*
    For grayscale, the path cost is the maximum intensity among all pixels
    the path, including `t'. In particular, when `s == a == NULL'
    (trivial path) returns the intensity of `t' as the handicap
    cost. For an `N'-channel image, each pixel is first reduced to its
    L_1 modulus, then scaled by `1/N'.
    
    When the seeds are the local minima, the IFT of this path-cost 
    function is the watershed transform. */

PathCost pf_sumval(PixelNode *s, RelArc *a, PixelNode *t, ImageGraph *G, int stage);
  /*
    For grayscale, the path cost is the sum of the intensities among
    all pixels the path, including `t'. In particular, when `s == a ==
    NULL' (trivial path) returns the intensity of `t' as the handicap
    cost. For an `N'-channel image, each pixel is first reduced to its
    L_1 modulus, then scaled by `1/N'.
    
    When the intensity is inversely proportional to the (isotropic)
    local travel speed, the IFT is the minimum-time Voronoi diagram. */

PathCost pf_maxdiff(PixelNode *s, RelArc *a, PixelNode *t, ImageGraph *G, int stage);
  /*
    For an `N'-channel image, the arc cost is the L_1 distance between
    pixel values, scaled by `1/(N*len)' where `len' is the Euclidean
    length of the arc `a'. (For grayscale, this reduces to absolute
    difference between pixel values, divided by `len'.) In either
    case, the path cost is the maximum of all arc costs. 
    
    When `s == a == NULL' returns an handicap cost of zero. */

PathCost pf_sumdiff(PixelNode *s, RelArc *a, PixelNode *t, ImageGraph *G, int stage);
  /*
    For an N-channel image, the arc cost is the L_1 distance between
    pixel values, scaled by 1/N. (For grayscale, this reduces to
    absolute difference between pixel values.) In either case, the
    path cost is the sum of all arc costs. */

PathCost pf_monoroot(PixelNode *s, RelArc *a, PixelNode *t, ImageGraph *G, int stage);
  /*
    For grayscale, the path cost is the intensity of `t' if the pixel
    intensities are non-decreasing along the path; else it is +oo.
    In particular, when `s == a == NULL' (trivial path) returns the pixel
    value at `t' as the handicap cost. For an `N'-channel image, 
    each pixel is first reduced to its L_1 modulus, then scaled by `1/N'.  
    
    When the seeds are all pixels, the roots of the IFT transform are
    all the local minima (FIFO queue policy) or one representative
    pixel from each local minimum (LIFO policy). */

PathCost pf_fuzzconn(PixelNode *s, RelArc *a, PixelNode *t, ImageGraph *G, int stage);
  /* 
    To be implemented. */

/* AUXILIARY FUNCTIONS */

PathCost arc_cost_quantize(double x);
  /* 
    Scales, discretizes, and clips `x' from `[0 _ 1]' to `[0.. MAX_FINITE_ARC_COST]'.
    The result is always an exact integer. Useful for many path-cost functions. */

double abs_pixel_value(PixelValue *p, int channels);
  /* 
    The L_1 modulus (sum of absolute intensities of all channels) of 
    pixel value `*p'. */

double abs_pixel_diff(PixelValue *p, PixelValue *q, int channels);
  /* 
    The L_1 distance (sum of absolute differences over all channels)
    between the pixel values `*p' and `*q'. */

#endif
