#ifndef float_array_H
#define float_array_H

/* Multi-dimensional arrays with {float} elements. */
/* Last edited on 2017-06-21 22:40:03 by stolfilocal */ 

#include <stdio.h>
#include <bool.h>

#include <array.h>
#include <array_io.h>
/* #include <array_linalg.h> */
#include <ix.h>
#include <ix_descr.h>

/* Multi-dimensional arrays of {float} elements.

  The representation uses the descriptors defined in {indexing.h} to
  allow constant-time index manipulation such as slicing, index reersal,
  index swapping (transposition), etc. The of max number of indices
  (dimensions, axes) is {float_array_NAXES}. Lower-dimensional arrays
  can be handled by setting the array size to 1 along the unwanted
  axes. */

array_typedef(float_array_t, float_array, float);
array_io_def(float_array_t, float_array, float);
/* array_linalg_def(float_array_t, float_array, float); */

#endif
