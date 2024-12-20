#ifndef spectrum_table_convert_H
#define spectrum_table_convert_H

/* Tools for converting between exact and binned power spectra of images. */
/* Last edited on 2024-12-05 10:30:56 by stolfi */ 

#include <stdio.h>
#include <bool.h>
#include <vec.h>
#include <spectrum_table_binned.h>
#include <spectrum_table_exact.h>

spectrum_table_binned_t spectrum_table_convert_exact_to_binned
  ( spectrum_table_exact_t *tx,
    uint32_t cols, 
    uint32_t rows );
  /* Converts an exact spectrum table to a binned table, 
     without blurring the entries.
     
     Specifically, the output table has one entry for each distinct
     pixel frequency {.fmid = sqrt((fx/nx)^2+(fy/ny)^2)} appearing in
     the Hartley transform; except that any two frequencies that yield
     the same value when converted to {float} are treated as equal and
     merged into the same entry. These `proper' entries have {.fmin ==
     .fmax == .fmid}. The first entry, for the constant term, has
     {.fmin == .fmid == .fmax = 0}. There are also `filler' entries
     with {.nTerms == .power == 0}, inserted between and around the
     proper entries, so that the ranges {[.fmin _ .fmax]} are
     contiguous and the last entry has {.fmax = .fmid = sqrt(0.5)}.
     The {.fmid} values in the filler lines are approximately halfway
     between {.fmin} and {.fmax}. */

#endif
