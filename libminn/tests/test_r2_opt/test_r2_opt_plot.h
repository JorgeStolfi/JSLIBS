/* Test plot tools for {r2_opt.h}. */
/* Last edited on 2017-06-06 23:20:09 by stolfilocal */

#ifndef test_r2_opt_plot_H
#define test_r2_opt_plot_H

#define _GNU_SOURCE
#include <stdint.h>

#include <r2.h>
#include <i2.h>
#include <bool.h> 
#include <r2_opt.h>

#include <test_r2_opt_basic.h>

void tr2o_plot_goal
  ( r2_opt_goal_func_t *f2, 
    int NI, 
    i2_t iscale, 
    r2_t ctr[],
    char *ctrtag,
    r2_t rad[]
  );
  /* Writes a  file "out/f2-{SX}-{SY}-{ctrtag}.dat" with a random 2D slice of the goal function
    {f2(NI,q[ku,kv],iscale)}, where {q[ku,kv]} ranges over a 2D grid of sampling arguments 
    in the neighborhood of {ctr}. The max displacement of {q[ku,kv][i]} from {ctr[i]}
    is a random fraction of the norm of {rad[i]}.  The grid is randomly oriented
    and it has {2*NS-1} points along each direction. */

void tr2o_write_test_image
  ( tr2o_image_eval_proc_t *eval, 
    int i,        /* Image index. */
    i2_t iscale,  /* Image shrink scale in each axis. */
    i2_t wsize,   /* Comparison widow size along each axis. */
    r2_t ctr,     /* Center of image. */
    char *ctrtag, /* Type of center, e.g. "ini" or "opt". */
    r2_t rad      /* Search radius at this scale. */
  );
  /* Writes to "out/image-{SX}-{SY}-{III}-{ctrtag}.pgm" the test image
    number {i} as a PGM grayscale image; where {SX,SY} are the
    components of {iscale} in 2 digits and {III} is the value of {i} in
    three digits, all in decimal zero-padded.
    
    The image is scaled by {1/2^iscale.c[j]} along each axis {j}. The
    pixel spacing of the image along each axis is 1 unit in the domain
    of the scaled image. The image is centered at {ctr} also scaled by
    the same factors, and spans a rectangle with side approximately
    {2*rad.c[j]+wsize.c[j]} centered on that point. Note that {rad} is
    not scaled but used "as is". */

void tr2o_plot_grid_normalize_to_span(int NI, r2_t u[], r2_t rad[]);
  /* Scales all vectors {u[0..N-1]} by the same scale factor, so 
    that at least one {u[i]} lies on the axis-aligned 
    ellipse with principal radii {rad[i]}, and all the others are inside
    their respective ellipses. */
       
double tr2o_compute_rel_span(int NI, r2_t u[], r2_t rad[]);
  /* Returns the max ratio of the norm of each {u[i]} to the radius of 
    the axis-aligned ellipse with principal radii {rad[i]} along the
    same direction as {u[i]}. */
    
#endif
