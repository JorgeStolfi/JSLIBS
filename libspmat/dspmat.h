#ifndef dspmat_H
#define dspmat_H
/* Sparse matrices with {double} entries, and operations thereon */

#define dspmat_H_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2009-08-31 21:49:04 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <bool.h>

#include <spmat.h>
#include <spmat_io.h>
#include <spmat_linalg.h>

spmat_typedef(dspmat_t, dspmat, double);
spmat_io_def(dspmat_t, dspmat, double);
spmat_linalg_def(dspmat_t, dspmat, double);

/* LIMITS */
/* They should be defined automatically by {spmat_typedef}, but CPP can't do that. */

#define dspmat_MAX_DIM (spmat_MAX_DIM)
#define dspmat_MAX_ROWS (spmat_MAX_ROWS)
#define dspmat_MAX_COLS (spmat_MAX_COLS)
  /* Max safe value for a row or column count ({dspmat_size_t}). */

#define dspmat_MAX_INDEX (spmat_MAX_INDEX)
  /* Max safe value for a valid row or column index ({dspmat_index_t}). */

#define dspmat_NO_INDEX (spmat_NO_INDEX)
  /* A {dspmat_index_t} value used to mean `no such column' or `no such row'. */

#define dspmat_MAX_COUNT (spmat_MAX_COUNT)
#define dspmat_MAX_ENTS (spmat_MAX_ENTS)
  /* Max safe value for a valid stored entry count ({dspmat_count_t}). */

#define dspmat_MAX_POS (spmat_MAX_POS)
  /* Max safe value for a valid stored entry index ({dspmat_pos_t}). */

#define dspmat_NO_POS (spmat_NO_POS)
  /* A {dspmat_pos_t} value used to mean `no such entry'. */

#endif
