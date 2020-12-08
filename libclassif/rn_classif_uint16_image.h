/* rn_classif_uint16_image.h --- builds an image from a vector classifier. */
/* Last edited on 2017-06-22 18:06:28 by stolfilocal */

#ifndef rn_classif_uint16_image_H
#define rn_classif_uint16_image_H

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>

#include <bool.h>
#include <r2.h>
#include <uint16_image.h>
#include <float_image.h>
#include <frgb.h>

#include <rn_classif.h>

uint16_image_t *rn_classif_uint16_image
  ( int NA, 
    int NC1, 
    rn_classif_labeler_t *lab1,
    int NC2, 
    rn_classif_labeler_t *lab2,
    int NCD, 
    rn_classif_dataset_t *D,
    double HD[],
    int classD[], 
    r2_t *ictr, 
    double HV, 
    int NXY, 
    int NSUB, 
    double sigma
  );
  /* Creates a PPM image for classification functions {NA,NC1,lab1}
    and {NA,NC2,lab2}, and the dataset {D}  with classification {classD}.
    Currently demands {NA == 2}.
    
    The image has {NXY} columns and rows of pixels and spans the
    square {ictr + V×V} in {\Reals^2}, where {V = [-HV _ +HV]}. The
    classifiers {NA,NC1,lab1} and {NA,NC2,lab2} are sampled
    systematically on every image pixel. The color of each pixel is
    computed by averaging the colors of {NSUB × NSUB} samples of {P}.
    
    For each of those samples, if {lab1} is not NULL, the class {cl1}
    assigned to the sample by {NA,NC1,lab1} should be in {0..NC1}, and
    is mapped to a color by an internal algorithm. Class 0 is white,
    other classes are strong colors of various hues. If {lab1} is
    null, it is tretaed as a labeler with no domains, and all samples
    are assigned the color white.
    
    Then, if {lab2} is not NULL, the saturation and brightness of the color
    derived from {NA,NC1,lab1} is modified according to the class {cl2}
    assigned by {NA,NC2,lab2}. Namely if {cl2 == 0} the color is
    preserved, and if {cl2} is 1 or higher the color is made
    proportionally darker and duller (less saturated). (Usually {lab1}
    is a classifier under test and {lab2} is the ideal classifier.)
    
    Before being classified, each sample is displaced by addition of a
    Gaussian noise with deviation {sigma} in each attribute, truncated
    to {4*sigma}. The client must seed the random generator properly.
    The classes {cl1,cl2} are computed for the perturbed sample but
    the color is painted in the pixel of the original (unperturbed) sample.
    
    Once the whole image has been painted as above, the samples in the
    dataset {D} (if not null) are painted as dots over it. These
    samples are painted with colors taken from {classD} (or black if
    {classD} is null) and possibly outlined in black.
    
    Furthermore, if {HD} is not null, each sample {i} is surrounded by
    a gray circle with radius {HD[i]} outlined in black. The dot
    positions are not perturbed by the {sigma}-noise. The dots and
    circles are antialiased. */

#endif
