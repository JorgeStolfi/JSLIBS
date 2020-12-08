/* Last edited on 2009-08-31 23:44:38 by stolfi */

/* DESCRIPTOR ALLOCATION */

ixd_t *ixd_new ( );
  /* Allocates a new descriptor record {D}, and initializes it with 
  {ixd_from_size(sz)}. */

ixd_t *ixd_copy ( ixd_t *D );
  /* Allocates a new descrptor record, and copies into it the descriptor {*D}. */

void ixd_free ( ixd_t *D );
  /* If {D} is not NULL, reclaims the space used by the descriptor
    record {*D}. {D} must have been created with {malloc}, {ixd_new},
    or {ixd_copy}. If {D} is NULL, does nothing. */

/* FORMATS 

  These are formats suitable to print values of the above types,
  at least in GNU's gcc 3.4.2. */

#define float_array_size_t_FMT "%llu"
#define float_array_index_t_FMP "%lld"
#define float_array_step_t_FMT "%lld"
#define float_array_pos_t_FMT "%llu"
#define float_array_element_count_t_FMT "%llu"

#define float_array_element_t_FMT "%13.7e"

/* ELEMENT RANGE */

void float_array_update_element_range(float_array_t *A, float *vmin, float *vmax);
  /* Expands the range {[*vmin _ *vmax]} so that it contains the values of
    all elements of array {A}.   If {A} is empty, the range is not changed. */

/* ACESSING ROWS (ONE-DIMENSIONAL SUB-ARRAYS)

  A \row\ of an array {A} is a one-dimensional sub-array of {A}, that
  spans the whole extent of {A} along some axis, and is trivial (has
  size 1) along the remaining axes.
  
  For the routines in this section, a row is specified by the number
  {i} of the non-trivial axis, and the index tuple {ix[0..ND-1]} of
  any of its elements. (Actually, the value of {ix[i]} is
  irrelevant.) 
  
  The routines get and set the elements in increasing order of the
  index {ix[i]}. The array {v} is scanned in the same order. */

void float_array_get_row(float_array_t *A, float_array_axis_t i, float_array_index_t ix[], float v[]);
  /* Copies into {v[0..m-1]} the row of {A} specified by {i} and {ix[0..ND-1]};
    where {m} is the size of {A} along axis {i} (irrespective of any aliasing). */

void float_array_set_row(float_array_t *A, float_array_axis_t i, float_array_index_t ix[], float v[]);
  /* Stores {v[0..m-1]} into the row of {A} specified by {i} and {ix[0..ND-1]};
    where {m} is the size of {A} along axis {i} (irrespective of any aliasing). */
    
void float_array_fill_row(float_array_t *A, float_array_axis_t i, float_array_index_t ix[], float v);
  /* Stores the value {v} into every element of the row of {A} specified by {i} and
    {ix[0..ND-1]}. */

/* ACESSING PLANES (TWO-DIMENSIONAL SUB-ARRAYS)

  A \plane\ of an array {A} is a two-dimensional sub-array of {A},
  that spans the whole extent of {A} along two distinct axes, and is
  trivial (has size 1) along the remaining axes.
  
  For the routines in this section, a plane is specified by the
  numbers {i}, {j} of the of the non-trivial axis, and the index tuple
  {ix[0..ND-1]} of any of its elements. (Actually, the values of
  {ix[i]} and {ix[j]} are irrelevant.)
  
  The routines get or set the elements in increasing order of index
  {ix[i]}, and, for the same {ix[i]}, in increasing order of {ix[j]}.
  The array {v} is scanned in the same order. */

void float_array_get_plane(float_array_t *A, float_array_axis_t i, float_array_axis_t j, float_array_index_t ix[], float v[]);
  /* Copies into {v[0..m*n-1]} the plane of {A} that is parallel to axes
    {i} and {j}, and contains the element with indices {ix[0..ND-1]}; where
    {m} and {n} are the sizes of {A} along axes {i} and {j}, respectively
    (irrespective of any aliasing). */

void float_array_set_plane(float_array_t *A, float_array_axis_t i, float_array_axis_t j, float_array_index_t ix[], float v[]);
  /* Stores {v[0..m*n-1]} into the plane of {A} that is parallel to axes
    {i} and {j}, and contains the element with indices {ix[0..ND-1]};
    where {m} and {n} are the sizes of {A} along axes {i} and {j}, respectively
    (irrespective of any aliasing). */
    
void float_array_fill_plane(float_array_t *A, float_array_axis_t i, float_array_axis_t j, float_array_index_t ix[], float v);
  /* Stores the value {v} into every element of the plane of {A} that
    is parallel to axis {i} and contains the element with indices
    {ix[0..ND-1]}. */

/* ACCESSING THE WHOLE ARRAY
  
  The routines in this section get or set the elements in increasing
  "C-like" order: namely, indices {ix[0..i-1]} are held constant while
  indices {ix[i..ND-1]} vary over all valid combinations, for all {i}.
  The array {v} is scanned in the same order. */

void float_array_get_array(float_array_t *A, float v[]);
  /* Copies into {v[0..ne-1]} all the elements of {A}; where {ne}
    is the number of valid index tuples of {A} (irrespective of any aliasing). */

void float_array_set_array(float_array_t *A, float v[]);
  /* Stores {v[0..ne-1]} into the elements of {A}; where {ne}
    is the number of valid index tuples of {A} (irrespective of any aliasing). */

void float_array_fill(float_array_t *A, float v);
  /* Stores {v} into every element of {A}. */

/* ARRAY INTERPOLATION

  For the procedures in this section, the pixel on column {ix} and 
  row {iy} is assumed to be centered at point {(ix+0.5,iy+0.5)}. */

float float_array_interpolate_elements(float_array_t *fim, int c, double x, double y);
  /* Returns the value of channel {c} of array {fim} at the point
    {(x,y)}, computed by bilinear interpolation from nearby
    elements. 
    
    More precisely, if {x} and {y} are of the form {ix + fx} and {iy +
    fy}, with integer {ix,iy} and {fx,fy} in [-0.5 _ +0.5], then
    returns the weighted average of the four pixels in columns {ix-1}
    and {ix}, with weights {0.5-fx} and {0.5+fx}, and rows {iy-1} and
    {iy}, with weights {0.5-fy} and {0.5+fy}.
    
    For this formula, pixels outside the array bounds are assumed to
    be filled by duplicating the extremal rows and columns of
    {fim}. */

void float_array_interpolate_pixels(float_array_t *fim, double x, double y, float v[]);
  /* For {c=0..NC-1}, stores into {v[c]} the value of chanel {c} of
    array {fim} at the point {(x,y)}, computed by bilinear
    interpolation from nearby elements; where {NC} is the number of
    channels of {fim}. */

/* LIMITS */

#define float_array_MAX_ABS_STEP ix_MAX_ABS_STEP
  /* Max absolute value of a {step} field. */

#define float_array_MAX_ELEMENTS ix_MAX_POSITIONS
  /* Max number of *distinct* element positions, not the
    the apparent domain size (the number of valid index tuples). The
    latter can exceed this limit if the array has replicated elements
    (null {step}s). */

#define float_array_MAX_SIZE (float_array_MAX_ELEMENTS)
  /* In theory a memoryless array could have a {size} greater than
    {float_array_MAX_ELEMENTS}, but that freedom seems to be of little use, and
    makes the code more complicated. Hence we enforce this limit even
    for memoryless arrays. */

#define float_array_MAX_INDEX ix_MAX
  /* We assume that {float_array_MAX_INDEX + 1} does not overflow an {ix_index_t}.
    The type is signed in order to avoid underflow, e.g. in loops like
    {for (index = N; index >=0; index--).} */

#define float_array_MAX_POS ix_MAX_POS
  /* Valid positions run from 0 to the total number of elements minus 1. */

#define float_array_MAX_BYTES SIZE_MAX
  /* The size argument for {malloc} is a {size_t}. */

