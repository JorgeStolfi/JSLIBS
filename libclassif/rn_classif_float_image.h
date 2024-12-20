/* rn_classif_float_image.h --- paints float images of vector classifiers. */
/* Last edited on 2024-12-05 10:24:22 by stolfi */

#ifndef rn_classif_float_image_H
#define rn_classif_float_image_H

#include <stdio.h>
#include <math.h>

#include <bool.h>
#include <r2.h>
#include <uint16_image.h>
#include <float_image.h>
#include <frgb.h>

#include <rn_classif.h>

void rn_classif_float_image_paint_labeler
  ( float_image_t *fim,
    int NA, 
    int NC1, 
    rn_classif_labeler_t *lab1,
    int NC2, 
    rn_classif_labeler_t *lab2, 
    r2_t *ictr, 
    double HV, 
    int NSUB, 
    double sigma,
    frgb_t cmap[]
  );
  /* Like {rn_classif_uint16_image_from_labeler} but paint on the given 
    float image {fim} with the given colormap {cmap}.  Currently requires
    {fim} to be square and have three channels. */

void rn_classif_float_image_paint_dataset
  ( float_image_t *fim,
    int NA, 
    int NCD, 
    rn_classif_dataset_t *D,
    double HD[],
    int classD[], 
    r2_t *ictr, 
    double HV,
    bool_t fillH,
    double penH,
    frgb_t cmap[],
    bool_t fillP,
    double radP,
    double penP, 
    int NSUB
  );
  /* Like {rn_classif_uint16_image_from_daatset} but paint on the given 
    float image {fim} with the given colormap {cmap}.  Currently requires
    {fim} to be square and have three channels.
    
    Each sample {i} has a position {pos[i]}, whose coordinates are the
    first two attributes of {D.smp[i]}; and a color {cmap[i]}. If
    {classD} and {cmap} are not NULL, and {classD[i]} is not zero,
    then {cmap[i]} will be {cmap[classD[i]]}, otherwise it will be
    black.
    
    First, if {HD} is not null, a disk with center {posi]} 
    and radius {HD[i]} is drawn for each sample {i}.  If {fillH}
    is true, the disk is painted with gray, else it is left unfilled.  In any case
    if {penH} is positive, the disk is outlined with a black border of width {penH}.
    
    Then, sample {i} is painted as a dot at position {pos[i]} and fixed
    radius {radP} (pixels). If {fillP} is true, the dot is filled with
    {cmap[i]}, otherwise it is left unfilled. In any
    case, if {penP} is positive, the (painted or unpainted) dot is
    outlined with a thin black border of width {penwd}.
    
    All painting uses antialiasing with {NSUB} samples per pixel. */

frgb_t rn_classif_float_image_compute_pixel
  ( int NA, 
    int NC1, 
    rn_classif_labeler_t *lab1, 
    int NC2, 
    rn_classif_labeler_t *lab2, 
    r2_t *pctr, 
    double HW, 
    int NSUB, 
    double sigma, 
    frgb_t cmap[]
  );
  /*  The {NA,NC1,lab1,NC2,lab2}, and {sigma} parameters have the same
    meaning as in {rn_classif_float_image_from_labeler}. Computes the color of
    the image pixel which corresponds to the square {pctr + W×W} in
    the domain of {P}, where {W = [-HW _ +HW]}. Namely classifies a
    grid of {NSUB×NSUB} vector samples spanning that square, maps the
    class of each sample to a color with the color map {cmap}, and
    averages those colors. */

#endif
