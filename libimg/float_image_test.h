#ifndef float_image_test_H
#define float_image_test_H

/* Tools for testing programs that deal with {float_image_t}. */
/* Last edited on 2023-01-10 15:32:19 by stolfi */ 

#define _GNU_SOURCE_
#include <stdint.h>

#include <float_image.h>
#include <r2.h>

typedef void float_image_test_generator_t (r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[]);
  /* Type of a procedure that defines a procedural image with {NC} channels, 
    {NX} columns, and {NY} rows.
    
    More precisely, a procedure of this type must return in
    {fs[0..chns-1]} the value of the pixel whose center is the point
    {p}. The procedure should remove frequencies above the Nyquist
    limit, assuming that the pixel spacing is 1 along each axis. Note
    that pixel centers have half-integer coordinates. */

void float_image_test_paint
  ( float_image_t *img, 
    float_image_test_generator_t *proc,
    int32_t ns
  );
  /* Paints into {img} the procedural image defined by {proc}.
    Uses the simple average of {ns × ns} subsamples inside each pixel. */

/* IMAGE GENERATORS: */

void float_image_test_gen_ripples(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[]);
  /* A set of circular ripples centered at random points. */
  
void float_image_test_gen_stripes(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[]);
  /* Horizontal green stripes plus vertical magenta stripes. */
    
void float_image_test_gen_checker(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[]);
  /* A checkerboard pattern, equal to the product of two waves. */

void float_image_test_gen_chopsea(r2_t *p, int32_t NC, int32_t NX, int32_t NY, float fs[]);
  /* A sum of several waves of various frequencies and directions. */

#endif
