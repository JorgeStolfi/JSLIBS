/* ift_image.h - image handling for the IFT routines. */
/* Last edited on 2006-11-21 19:32:18 by stolfi */

#ifndef ift_image_H
#define ift_image_H

#define _GNU_SOURCE
#include <jspnm.h>
#include <jspnm_image.h>
#include "ift.h"

void ift_set_values_from_image(pnm_image_t *img, ImageGraph *G);
  /* 
    Sets the `y' field of all nodes from the image pixels. */

void ift_set_seeds_from_image(pnm_image_t *seed_img, ImageGraph *G);
  /* 
    The parameter `seed_img' should be a PGM (grayscale) image that
    specifies the seed pixels and their labels.  If `p' is the value of a
    pixel, the label field `L' of the corresponding node is set to `p'.
    Note that if `p == NON_SEED_LABEL' the node will *not* be a seed. */

pnm_image_t *ift_get_cost_image(ImageGraph *G, PathCost maxcost, pnm_sample_t maxval);
  /* 
    Returns a PGM image, with pixels in `[0..maxval]', containing the
    cost (`C') fields of all nodes. If `C = INFINITE_COST', the pixel
    is set to `maxval'; otherwise `C' is scaled by `(maxval-1)/maxcost',
    rounded, and clipped to the range `[0..maxval-1]'. */

pnm_image_t *ift_get_pred_image(ImageGraph *G);
  /* 
    Returns a PGM image containing the predecessor (`P') fields of all nodes.
    Each P value is encoded as two signed index displacements `(dh,dv)' 
    in [-127..+127], encoded as a 16-bit value by `(dv+128)*256 + (dh+128)'. */

pnm_image_t *ift_get_label_image(ImageGraph *G, pnm_sample_t maxval);
  /* 
    Returns a PGM image containing the label (`L') fields of all nodes,
    which must lie in the raneg [0..maxval]'. */

pnm_image_t *ift_get_root_image(ImageGraph *G);
  /* 
    Returns a PGM image with `maxval = 1', where each pixel has value 
    1 if it is a root of the path tree, and 0 otherwise. */

pnm_image_t *ift_get_spread_image(ImageGraph *G, pnm_sample_t maxval);
  /* 
    Returns a copy of the original image, with the same type as the
    original (PNM or PGM) and the specified `maxval', where each pixel
    has been replaced by the value of its corresponding root pixel. */

pnm_image_t *ift_get_single_label_image
  ( ImageGraph *G, 
    SeedLabel label, 
    pnm_sample_t bg[], 
    pnm_sample_t maxval
  );
  /* 
    Returns an image with the same size and type (PPM or PGM) as the original image,
    and the specified `maxval', where pixels with the given `label' are copied 
    from the original, while pixels wil all other labels are set to the value
    `bg[0..NC-1]'. */

void ift_write_boxes(FILE *wr, ImageGraph *G, pnm_sample_t maxval, int margin);
  /* 
    Writes to `wr' one line for each label in `0..maxval',
    containing the bounding box of the corresponding region, in
    the format `LABEL HMIN VMIN HSIZE VSIZE'.
    
    The bounding box of a non-empty region is the tightest box that
    contains all its pixels, expanded by `margin' pixels in all four
    directions, then clipped against the image bounds. A negative
    value for `margin' causes the box to be shrunk instead. 
    
    If the region is empty, or its box shrinks to nothing, the
    fields `HMIN VMIN HSIZE VSIZE' are all zero, irrespective
    of the `margin'. */

#endif

