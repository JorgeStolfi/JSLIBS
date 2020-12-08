#ifndef float_image_geostereo_multiscale_H
#define float_image_geostereo_multiscale_H

/* Tools for multiscale geometric stereo reconstruction from image pairs. */
/* Last edited on 2017-06-25 03:11:17 by stolfilocal */ 

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>
#include <float_image.h>

void float_image_geostereo_multiscale
  ( float_image_t *f1,  /* Image 1. */
    float_image_t *f2,  /* Image 2. */
    int nscales,        /* Number of scales to consider (0 = uniscale). */
    int ncands,         /* Number of candidates to keep. */
    int rx,             /* Window half-width. */
    int ry,             /* Window half-height. */
    double dmin,        /* Minimum signed displacement (pixels). */
    double dmax,        /* Maximum signed displacement (pixels). */
    float_image_t **fd, /* (OUT) Dispmap image. */
    float_image_t **fs  /* (OUT) Scoremap image. */
  );
  /* Same as {float_image_geostereo_uniscale_displacement_map}, but uses a
    multiscale search algorithm with {nscales} levels deep. (If
    {nscales = 0}, uses uniscale search.)
    
    At each scale {k}, the procedure uses versions of {img1}
    and {img2} shrunk by a factor of {1/2^k} in both
    directions.  !!! FINISH EXPLANATION !!! */
  
void float_image_geostereo_refine_and_prune_displacement_map
  ( float_image_t *gd,  /* Displacement map for halfsize images. */
    float_image_t *gs,  /* Score map for halfsize images. */
    float_image_t *f1,  /* Full-size image 1. */
    float_image_t *f2,  /* Full-size image 2. */
    int rx,             /* Window half-width. */
    int ry,             /* Window half-height. */
    double dmin,        /* Minimum signed displacement (pixels). */
    double dmax,        /* Maximum signed displacement (pixels). */
    float_image_t *fd,  /* (OUT) Dispmap for full-size images. */
    float_image_t *fs   /* (OUT) Scoremap for full-size images. */
  );
  /* Arguments {gd,gs} should be the displacements and score maps for
    half-size versions of images {f1,f2}. For each pixel of {f1,f2},
    takes the best {NC} displacements found in {gd} for the
    corresponding point, doubles those displacements (to match the
    scale of {f1,f2}), tries some small adjustments on them (by 1/3
    pixel), and finally stores the results in {fd,fs}.

    Assumes that {fd,fs} have already been allocated. The number {NC}
    is the number of channels of {fd}, which must not be greater than
    that of {gd}. The displacement scores in {fs} will be a
    combination of the rough scores in {gs} and the newly computed
    scores for the adjusted displacements. */

#endif
