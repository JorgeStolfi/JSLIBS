/* pnmift_image.h - image handling for the IFT routines. */
/* Last edited on 2024-12-05 10:29:15 by stolfi */

#ifndef pnmift_image_H
#define pnmift_image_H

#include <jspnm.h>
#include <uint16_image.h>
#include <frgb.h>

#include <ift.h>

#define pnmift_MAX_ROWS ift_MAX_ROWS
#define pnmift_MAX_COLS ift_MAX_COLS

#define pnmift_MAX_LABEL PNM_FILE_MAX_MAXVAL
  /* Maximum seed label that may appear. */

#define pnmift_NON_SEED_LABEL 0
  /* A label value that marks the pixel as NOT a seed. */
  
#define pnmift_DEFAULT_SEED_LABEL 1
/* The default label for seeds. */

void pnmift_set_values_from_image(uint16_image_t *img, ift_graph_t *G, frgb_t rgb[]);
  /* Sets the {y} field of all nodes from the image pixels. */

void pnmift_set_labels_from_image(uint16_image_t *seed_img, ift_graph_t *G, uint16_t label[], bool_t verbose);
  /* The parameter {seed_img} should be a PGM (grayscale) image that
    specifies the seed pixels and their labels.  For each pixel index {i} in {G},
    the procedure set {label[i]} to the raw value {p} of the corresponding pixel
    in {seed_img}.  Does not modify the costs {G->node[i].C}. */

uint16_image_t *pnmift_get_cost_image(ift_graph_t *G, ift_path_cost_t maxcost, uint16_t maxval);
  /* 
    Returns a PGM image, with pixels in {[0..maxval]}, containing the
    cost ({C}) fields of all nodes. If {C = INFINITE_COST}, the pixel
    is set to {maxval}; otherwise {C} is scaled by {(maxval-1)/maxcost},
    rounded, and clipped to the range {[0..maxval-1]}. */

uint16_image_t *pnmift_get_pred_image(ift_graph_t *G);
  /* 
    Returns a PGM image containing the predecessor ({P}) fields of all nodes.
    Each P value is encoded as two signed index displacements {(dh,dv)} 
    in [-127..+127], encoded as a 16-bit value by {(dv+128)*256 + (dh+128)}. */

uint16_image_t *pnmift_get_label_image(ift_graph_t *G, uint16_t label[], uint16_t maxval);
  /* Returns a PGM image where each pixel {p} is set to {label[k]} where {k} is the
    index of the root of the forest that contains {p}.
    The label must lie in the range [0..maxval]. */

uint16_image_t *pnmift_get_root_image(ift_graph_t *G);
  /* 
    Returns a PGM image with {maxval = 2}, where each pixel has value 
    2 if it is a root of the path tree, and 0 otherwise. */

uint16_image_t *pnmift_get_spread_image(ift_graph_t *G, frgb_t rgb[], int chns, uint16_t maxval);
  /* Returns a copy of the original image, with the same type as the
    original (PNM or PGM) and the specified {maxval}, where each pixel
    has been replaced by the value of its corresponding root pixel. */

uint16_image_t *pnmift_get_single_label_image
  ( ift_graph_t *G, 
    uint16_t label[], 
    uint16_t lab, 
    frgb_t rgb[],
    uint16_t bg[], 
    int chns,
    uint16_t maxval
  );
  /* Returns an image with the same size and type (PPM or PGM) as the original image,
    and the specified {maxval}, where every pixel {p} from selected trees is copied 
    from the original color {rgb[p]}, while pixels with all other labels are set to the value
    {bg[0..chns-1]}. A tree is selecetd if its root pixel {r} has {label[r] = lab}. */

void pnmift_write_boxes(FILE *wr, ift_graph_t *G, uint16_t label[], uint16_t maxval, int margin);
  /* Writes to {wr} one line for each label in {0..maxval},
    containing the bounding box of the corresponding region, in
    the format {LABEL HMIN VMIN HSIZE VSIZE}.
    
    The bounding box of a non-empty region is the tightest box that
    contains all its pixels, expanded by {margin} pixels in all four
    directions, then clipped against the image bounds. A negative
    value for {margin} causes the box to be shrunk instead. 
    
    If the region is empty, or its box shrinks to nothing due to a
    negative {margin}, the fields {HMIN VMIN HSIZE VSIZE} are all
    zero, irrespective of the {margin}. */

#endif

