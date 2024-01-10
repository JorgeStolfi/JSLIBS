#ifndef float_image_geostereo_H
#define float_image_geostereo_H

/* Tools for geometric stereo reconstruction from image pairs. */
/* Last edited on 2017-06-25 02:50:16 by stolfilocal */ 

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>
#include <float_image.h>

void float_image_geostereo_uniscale
  ( float_image_t *f1,   /* Image 1. */
    float_image_t *f2,   /* Image 2. */
    int ncands,          /* Number of candidates to keep. */
    int rx,              /* Window half-width. */
    int ry,              /* Window half-height. */
    double dmin,         /* Minimum signed displacement (pixels). */
    double dmax,         /* Maximum signed displacement (pixels). */
    float_image_t **fdP, /* (OUT) Dispmap image. */
    float_image_t **fsP  /* (OUT) Scoremap image. */
  );
  /* Computes the disparity map for two images {f1,f2}.
  
    The two input images must have the same dimensions and channel
    counts, and must have been geometrically corrected so that rows with
    the same index are corresponding epipolar lines.
  
    The returned displacement map {fd} is an image with same dimensions
    as {img1} and {img2}, but with {ncands} channels. Each pixel of {fd}
    in column {x} and row {y} contains the best {ncands} displacement
    values found for the corresponding pixels of the two images.
    
    More precisley, a displacement value {d} for that pixel means that
    the local values around row {y} and columns {x+d} of {img1} seem to
    match those around row {y} and column {x-d} of {img2}.
    
    The displacement {d} can be fractional, with granularity 0.3 pixel.
    The input images are interpolated when {d} is not an integer. The
    procedure considers only values {d} in the range {[dmin..dmax]}.
    
    The comparison is made using a sampling window with {2*rx+1} columns
    and {2*ry+1} rows centered at the points in question.
    
    The procedure also returns a score map {fs}, with same dimensions
    and with {ncands} channels. Each sample of {fs} is the discrepancy
    {s} between the two image neighborhoods implied by the corresponding
    sample {d} from the displacement map. The discrepancy is defined as
    the average squared difference of the coresponding samples of {img1}
    and {img2} in the respective sampling windows.
    
    For each pixel of the image domain, the procedure enumerates all displacements
    of the form {i/3}, where {i} is an integer, in the range {[dmin _ dmax]};
    and computes their scores {sc[i]}.  Then it looks for local minima 
    in the list {sc[]}, and stores the {ncands} smallest minima in {fd,fs}.
    
    The images {fd,fs} are allocated by the procedure, and retuned in
    {*fdP} and {*fsP}, if these are not {NULL}. */

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

void float_image_geostereo_local_match
  ( float_image_t *f1,  /* Image 1. */
    float_image_t *f2,  /* Image 2. */
    double x,           /* Central horiz position (pixels, fractional). */
    int32_t y,          /* Current row index in image. */
    double dmin,        /* Minimum signed displacement (pixels). */
    double dmax,        /* Maximum signed displacement (pixels). */
    int rx,             /* Window half-width. */
    int ry,             /* Window half-height. */
    double *dbest,      /* Adjusted displacement. */
    double *sbest,      /* Score (squared mismatch) for {dbest}. */
    float *w1,          /* Buffer for image 1 window samples. */
    float *w2           /* Buffer for image 2 window samples. */
  );
  /* Returns in {*dbest} the displacement in the range {[dmin..dmax]}
    that has minimum score, and in {*sbest} its score. All displacements
    are in multiples of 1/3 of a pixel. The score is the total squared
    sample difference in windows of {f1} and {f2}
    with dimensions {[-rx..+rx]x[-ry..+ry]} centered on row {y} columns
    {x-d/3} and {x+d/3}, respectively; after normalizing each window
    pixel to mean 0 and variance 1. The buffers {w1,w2} must have
    enough space for all window samples. */

#endif
