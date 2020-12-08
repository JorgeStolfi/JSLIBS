#ifndef float_image_geostereo_H
#define float_image_geostereo_H

/* Tools for geometric stereo reconstruction from image pairs. */
/* Last edited on 2017-06-26 00:08:32 by stolfilocal */ 

#define _GNU_SOURCE
#include <stdint.h>

#include <bool.h>
#include <float_image.h>

void float_image_geostereo_single_pixel_best
  ( float_image_t *f1,  /* Image 1. */
    float_image_t *f2,  /* Image 2. */
    int32_t x,          /* Central horiz position. */
    int32_t y,          /* Current row index in image. */
    int32_t nwx,        /* Window width. */
    int32_t nwy,        /* Window height. */
    double wt[],        /* Window weights. */
    double dmin,        /* Minimum signed displacement (pixels). */
    double dmax,        /* Maximum signed displacement (pixels). */
    int32_t ncands,     /* Number of best displacements to keep. */
    double dbest[],     /* Best {ncands} displacements. */
    double sbest[],     /* Score (squared mismatch) for {dbest[0..ncands-1]}. */
    bool_t debug,       /* If true, debugs computations. */
    double smp1[],      /* (WORK) Buffer for image 1 window samples. */
    double smp2[]       /* (WORK) Buffer for image 2 window samples. */
  );
  /* Computes the {ncands} most promising displacement candidates for a stereo
    image pair {f1,f2} at column {x} and row {y} of their domain.
    
    The procedure enumerates all displacements {d} that are multiples of
    {1/3} of pixel and are in the range {[dmin _ dmax]}. The score for
    each {d} is computed by
    {float_image_geostereo_single_score(f1,f2,x,y,d,nwx,nwy,wt,smp1,smp2)}.
    See that procedure for a description of the parameters.

    A displacement /candidate/ is a diplacement {d} in that sequence
    that is a local minimum of the discrepancy score. The best {ncands}
    displacement candidates are returned in {dbest[0..ncands-1]}, and
    their scores in {sbest[0..ncands-1]}, sorted by increasing score. If
    there are not enough candidates for that pixel, the lists {dbest}
    and {sbest} are completed with {-INF} and {+INF}, respectively.
    
    The working buffers {smp1,smp2} must have enough space for all
    window samples, namely {NC*nwx*nwy} where {NC} is the number of
    channels of {f1} and {f2}. */

double float_image_geostereo_single_disp_score
  ( float_image_t *f1,  /* Image 1. */
    float_image_t *f2,  /* Image 2. */
    uint32_t x,         /* Column index in the map domain. */
    uint32_t y,         /* Row index in the map domain. */
    double d,           /* Parallax displacement (pixels). */
    uint32_t nwx,       /* Window width. */
    uint32_t nwy,       /* Window height. */
    double wt[],        /* Window weights. */
    bool_t debug,       /* TRUE to debug the computations. */
    double smp1[],      /* (WORK) Buffer for window samples of image 1. */
    double smp2[]       /* (WORK) Buffer for window samples of image 2. */
  );
  /* Returns the discrepancy score between images {f1} and {f2} for the
    parallax displacement {d} at the pixel with column index {x} and row
    index {y}.  The two images must have the same dimensions and the
    same number of channels.
    
    Namely, compares the neighborhood around row {y} and column {x+d} of
    {img1} with the neighborhood around row {y} and column {x-d} of
    {img2}. Each neigborhood consists of a sampling window with {nwx}
    columns and {nwy} rows, centered at the specified indices; both
    window dimensions must be odd.
    
    The displacement {d} can be fractional, in which case the window
    samples are obtained by smooth interpolation of the image samples.
    See {float_image_geostereo_get_samples} below.
    
    If {normalize} is true, the samples in
    each channel of each window are normalized to have mean 0 and unit
    variance before the comparison.
    
    Either way, the discrepancy score is computed from the window
    samples by {float_image_geostereo_compute_score} below, using the
    window pixel weights {wt[0..nwx*nwy-1]}.
    
    The buffers {smp1} and {smp2} must have enough space for all window
    samples; that is, {nwx*nwy*NC} elements where {NC} is the number of
    channels in the input images. */

double float_image_geostereo_compute_score
  ( double smp1[], 
    double smp2[], 
    uint32_t nwx, 
    uint32_t nwy, 
    double wt[],
    uint32_t NC
  );
  /* Compute the discrepancy score from the sample arrays {smp1,smp2}.
    
    The discrepancy score is defined as the average squared difference
    of the coresponding samples of {img1} and {img2} in the respective
    sampling windows. The discrepancy is zero if the two images are
    identical inside the respective windows.
    
    Assumes that these samples are taken from two images, both with
    {NC} channels, within windows of {nwx} columns and {nwy} rows,
    as described in {float_image_geostereo_get_samples} below.
    
    The window pixel weight table {wt} must have {nwx*nwy} elements. For
    each channel {ic}, the two samples in column {ix} and row {iy} of
    the windows contribute to the score with weight {wt[ix + nwy*iy]}. */
    

void float_image_geostereo_get_samples
  ( float_image_t *f,   /* Pixel row buffer for image 1. */
    uint32_t x,         /* Reference column index. */
    uint32_t y,         /* Reference row index. */
    double d,           /* Horizontal displacement (pixels, fractional). */
    uint32_t nwx,       /* Window width. */
    uint32_t nwy,       /* Window height. */
    double smp[]         /* (OUT) Window sample buffer. */
  );
  /* Extracts samples from image {f} contained in a window centered on
    column {x+d} and row {y}.
    
    The sampling window will have {nwx} columns and {nwy} rows, centered
    at the specified indices. The window dimensions {nwx} and {nwy} must
    be odd. The displacement {d} can be fractional, in which case the
    window samples are obtained by smooth interpolation of the image
    samples.
    
    The output array {smp} must have enough space for all window
    samples; that is, {nwx*nw*NC} elements where {NC} is the number of
    channels in the input image. 
    
    The sample on column {ix}, row {iy}, and channel {ic} of the window
    is stored in {smp[ic + NC*(ix + nwx*iy)]}. */

/* INTERNAL TOOLS */

void float_image_geostereo_debug_window
  ( double smp[],
    uint32_t nwx,       /* Window width. */
    uint32_t nwy,       /* Window height. */
    uint32_t NC         /* Channel count. */
  );
  /* Prints to {stderr} the window samples {smp[0..NC*nwx*nwy-1]}. See
    {float_image_geostereo_get_samples} for the meaning of parameters
    and the order of samples in {smp}. */

int32_t float_image_geostereo_queue_insert(double d, double s, int32_t nq, double dq[], double sq[]);
  /* Inserts the displacement candidate {d}
    and its score {s} in the queue {dq[0..nq-1]} and {sq[0..nq-1]}.
    The queue must be sorted in order of non-decreasing {sq}.
    
    If {s} is strictly less than the last entry {sq[nq-1]},
    the procedure inserts {d} and {s} in the proper order, shifting 
    the entries with larger score and discarding the last one.
    In case of ties, {d,s} is inserted after the last tied entry.
    In that case the procedure returns the final rank (index)
    of {d,s} in the queue.
    
    If {s} is greater than or equal to the last entry {sq[nq-1]},
    the procedure does nothing and returns {nq} as the rank.*/
    
void float_image_geostereo_queue_dump(int32_t nq, double dq[], double sq[]); 

#endif
