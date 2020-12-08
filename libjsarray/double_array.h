#ifndef double_array_H
#define double_array_H

/* Multi-dimensional arrays with {double} elements. */
/* Last edited on 2019-12-07 18:52:08 by jstolfi */ 

#include <stdio.h>
#include <bool.h>
#include <indexing.h>
#include <indexing_descr.h>

#include <array.h>
#include <array_io.h>
/* #include <array_linalg.h> */

#define double_array_MAX_AXES (array_MAX_AXES)
  /* Number of indices (dimensions, axes) in any array. Lower-dimensional
    arrays are obtained by setting the size to 1 along the unwanted axes. */

array_typedef(double_array_t, double_array, double);
array_io_def(double_array_t, double_array, double);
/* array_linalg_def(double_array_t, double_array, double); */

#endif
