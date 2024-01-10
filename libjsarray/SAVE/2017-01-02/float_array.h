#ifndef float_array_H
#define float_array_H

/* Multi-dimensional arrays with {float} elements. */
/* Last edited on 2010-04-23 10:02:40 by stolfi */ 

#include <stdio.h>
#include <bool.h>

#include <array.h>
#include <array_io.h>
/* #include <array_linalg.h> */
#include <indexing.h>
#include <indexing_descr.h>

#define float_array_MAX_AXES (array_MAX_AXES)
  /* Number of indices (dimensions, axes) in any array. Lower-dimensional
    arrays are obtained by setting the size to 1 along the unwanted axes. */

array_typedef(float_array_t, float_array, float);
array_io_def(float_array_t, float_array, float);
/* array_linalg_def(float_array_t, float_array, float); */

#endif
