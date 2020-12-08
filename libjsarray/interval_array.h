#ifndef interval_array_H
#define interval_array_H

/* Multi-dimensional arrays with {interval_t} elements. */
/* Last edited on 2010-04-23 10:22:02 by stolfi */ 

#include <stdio.h>
#include <bool.h>
#include <interval.h>

#include <array.h>
#include <array_io.h>
/* #include <array_linalg.h> */
#include <indexing.h>
#include <indexing_descr.h>

#define interval_t_array_MAX_AXES (array_MAX_AXES)
  /* Number of indices (dimensions, axes) in any array. Lower-dimensional
    arrays are obtained by setting the size to 1 along the unwanted axes. */

array_typedef(interval_array_t, interval_array, interval_t);
array_io_def(interval_array_t, interval_array, interval_t);
/* array_linalg_def(interval_array_t, interval_array, interval_t); */

#endif
