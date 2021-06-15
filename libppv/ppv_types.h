/* Basic types and limits for portable 6-dimensional sample arrays. */
/* Last edited on 2021-06-13 12:25:06 by jstolfi */

#ifndef ppv_types_H
#define ppv_types_H

#include <bool.h>
#include <sign.h>
#include <indexing.h>

/* SAMPLES */

typedef uint8_t ppv_nbits_t;  /* Number of bits per sample or per word. */
  /* Number of bits in a sample or storage word. */

typedef uint32_t ppv_sample_t;  
  /* The value of a sample. */

/* STORAGE WORDS */

typedef uint8_t  ppv_word_08_t;
typedef uint16_t ppv_word_16_t;
typedef uint32_t ppv_word_32_t;
  /* The legal word sizes for the current implementation are 
    8, 16, and 32 bits.  Array storage areas are vectors
    of one of these types. */
   
typedef uint32_t ppv_word_t; 
  /* The largest storage word. */

/* AXES */

typedef ix_dim_t ppv_dim_t;
   /* Order (number of indices) of an array or sub-array.  
     Legitimate values are {0..ppv_NAXES}. */
   
typedef ix_axis_t ppv_axis_t;
  /* Specifies an axis (index). Legitimate values are {0..ppv_NAXES-1},
     but {ppv_NAXES} may be used as a `null' value. */

typedef ix_index_op_t ppv_index_op_t;
typedef ix_index_pos_op_t ppv_index_pos_op_t;
typedef ix_index_pos2_op_t ppv_index_pos2_op_t;
typedef ix_index_pos3_op_t ppv_index_pos3_op_t;
  /* Type of a procedure that operates on an index tuple. */

/* SIZES, STEPS, POSITIONS */

typedef ix_step_t ppv_step_t; 
   /* Sample position increment along some axis. */
   
typedef ix_count_t ppv_sample_count_t; 
  /* Total number of samples. */
  
typedef ix_size_t ppv_size_t; 
  /* Number of samples along some axis. */

typedef ix_index_t ppv_index_t; 
  /* Index of a sample. Note that it is signed, even though 
    legitimate values are always non-negative. */

typedef ix_pos_t ppv_pos_t; 
  /* Linearized index of a sample. */

/* LIMITS */

#define ppv_MAX_DIM (ix_MAX_DIM)
  /* Max number of axes (indices) in an array. */

#define ppv_MAX_AXIS (ix_MAX_AXIS)
  /* Should be {ppv_MAX_DIM-1}. Axes are numbered {0..ppv_MAX_AXIS}. */

#define ppv_MAX_BPS (32)
  /* Maximum number of bits per sample. */

#define ppv_MAX_SAMPLE_VAL UINT32_MAX
  /* I.e. {2^ppv_MAX_BPS - 1}. */

#define ppv_MAX_BPW (32)
  /* This must be large enough to hold the largest sample and the 
    largest storage word.  Every decent machine has 32-bit ints, but not every compiler can handle
    64-bit words, so we limit words and samples to 32 bits. This may change 
    in the future.  */

#define ppv_MAX_ABS_STEP ix_MAX_ABS_STEP
  /* Max absolute value of a {step} field. */

#define ppv_MAX_SAMPLES ix_MAX_POSITIONS
  /* Max number of *distinct* sample positions, not the
    the apparent domain size (the number of valid index tuples). The
    latter can exceed this limit if the array has replicated samples
    (null {step}s). */

#define ppv_MAX_SIZE (ppv_MAX_SAMPLES)
  /* In theory a memoryless array could have a {size} greater than
    {ppv_MAX_SAMPLES}, but that freedom seems to be of little use, and
    makes the code more complicated. Hence we enforce this limit even
    for memoryless arrays. */

#define ppv_MAX_INDEX (ppv_MAX_SAMPLES - 1)
  /* We assume that {ppv_MAX_INDEX + 1} does not overflow an {ix_index_t}.
    The type is signed in order to avoid underflow, e.g. in loops like
    {for (index = N; index >=0; index--).} */

#define ppv_MAX_POS (ppv_MAX_SAMPLES - 1)
  /* Valid positions run from 0 to the total number of samples minus 1. */

#define ppv_MAX_BYTES SIZE_MAX
  /* The size argument for {malloc} is a {size_t}. */

/* FORMATS 

  These are formats suitable to print values of the above types,
  at least in GNU's gcc 3.4.2. */

#define ppv_size_t_FMT "%lu"
#define ppv_index_t_FMP "%ld"
#define ppv_step_t_FMT "%ld"
#define ppv_pos_t_FMT "%lu"
#define ppv_sample_count_t_FMT "%lu"

#define ppv_sample_t_FMT "%u"

#endif
