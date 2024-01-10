/* Basic types and limits for portable 6-dimensional sample arrays. */
/* Last edited on 2005-05-29 02:12:08 by stolfi */

#ifndef ppv_types_H
#define ppv_types_H

/* TYPES AND LIMITS */

typedef unsigned char ppv_axis_t;
#define ppv_NAX 6
  /* Designates one of the indices of a {ppv_array_t}. Valid values
    are {0..ppv_NAX-1}. */

typedef unsigned char ppv_nbits_t;  /* Number of bits per sample or per word. */
#define ppv_MAX_BPW (32)
#define ppv_MAX_BPS (32)
  /* Every decent machine has 32-bit ints, but not every compiler can handle
    64-bit words, so we limit words and samples to 32 bits. This may change 
    in the future.  */

typedef unsigned int ppv_word_t; /* The largest storage word. */
  /* This must be large enough to hold the largest sample and the 
    largest storage word. */

typedef unsigned int ppv_sample_t;  /* The value of a sample. */
#define ppv_MAX_SAMPLE_VAL (4294967295u)
  /* I.e. {2^ppv_MAX_BPS - 1}. */
   
typedef signed int ppv_step_t; /* Sample position increment along some axis. */
#define ppv_MAX_ABS_STEP (2147483647u)
  /* Array flipping requires negating a {step}, hence {ppv_step_t}
    must be signed. Assuming we don't have 64-bit ints,
    we settle for {ppv_MAX_ABS_STEP = 2^31-1}. */

typedef unsigned int ppv_tot_size_t; /* Total number of samples. */
#define ppv_MAX_SAMPLES (ppv_MAX_ABS_STEP + 1)
  /* By combining subsampling with diagonal extraction a client could
    get the {step} of an array to be as high as {ppv_MAX_SAMPLES-1},
    which leads to this limit on ppv_MAX_SAMPLES (2^31 = 2 Giga
    samples). 
    
    Note that this is a limit on *distinct* sample locations, not on
    the apparent domain size (the number of valid index tuples). The
    latter can exceed this limit if the array has replicated samples
    (null {step}s). */

typedef unsigned int ppv_size_t; /* Number of samples along some axis. */
#define ppv_MAX_SIZE (ppv_MAX_SAMPLES)
  /* There is little advantage in allowing a {size} along any axis to
    be greater than {ppv_MAX_SAMPLES}. To keep the code simple, we
    enforce this limit even for memoryless arrays. */

typedef unsigned int ppv_index_t; /* Index of a sample */
#define ppv_MAX_INDEX (ppv_MAX_SAMPLES - 1)
  /* Sample indices cannot be negative, so we use unsigned ints. Note
    that {ppv_MAX_INDEX + 1} does not overflow, but watch out for 
    underflow in loops like {for (index = N; index >=0; index--).} */

typedef unsigned int ppv_pos_t; /* Linearized position of a sample. */
#define ppv_MAX_POS (ppv_MAX_SAMPLES - 1)
  /* Valid positions run from 0 to the total number of samples minus 1. */

#define ppv_MAX_BYTES (4294967295u)
  /* The size argument for {malloc} is a {size_t}, which is an
    unsigned int on most machines. */

#ifndef HAVE_BOOL
typedef enum { FALSE=0, TRUE=1 } bool_t;
#define FALSE FALSE
#define TRUE TRUE
#define HAVE_BOOL 1
#endif

#ifndef HAVE_SIGN
typedef enum { NEG=-1, ZER=0, POS=+1 } sign_t;
#define NEG NEG
#define ZER ZER
#define POS POS
#define HAVE_SIGN 1
#endif

typedef unsigned char  ppv_word_08_t;
typedef unsigned short ppv_word_16_t;
typedef unsigned int   ppv_word_32_t;
  /* The legal word sizes for the current implementation are 
    8, 16, and 32 bits.  Array storage areas are vectors
    of one of these types. */

#endif
