#ifndef spectrum_table_binned_H
#define spectrum_table_binned_H

/* Tools for gathering binned (histogram-like) power spectra of images. */
/* Last edited on 2024-12-05 10:30:54 by stolfi */ 

#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <vec.h>
#include <float_image.h>

typedef struct spectrum_table_binned_entry_t 
  { double fmin;   /* Low end of frequency range. */
    double fmax;   /* High end of frequency range. */
    double fmid;   /* Nominal mean frequency. */
    double nTerms; /* Number of Hartley terms in range. */
    double power;  /* Total power in range. */
  } spectrum_table_binned_entry_t;
  /* An entry of a binned power spectrum table. The entry covers the
    frequency range {[fmin _ fmax]}. The nominal mean frequency {freq}
    should be in that range. The power and count of a Hartley term may
    be split between two or more consecutive entries, so the field
    {nTerms} may be fractional. */

vec_typedef(spectrum_table_binned_t,spectrum_table_binned,spectrum_table_binned_entry_t); 

int32_t spectrum_table_binned_locate_entry(spectrum_table_binned_t *tb, double f);
  /* Returns the index of entry in {tb->e[0..n-1]} whose range
    {fmin,fmax} contains the frequency {f}; where {n = tb.ne}.

    Assumes that the table entries have consecutive ranges (meaning
    that the {fmax} of one entry equal to the {fmin} of the next one)
    and that the ranges have positive width (meaning that {fmin} is
    less than {fmax}). The result is always in {0..n-1}. If {f} is
    less than the {fmin} of entry 0, returns 0; if {f} is greater than
    the {fmax} of entry {n-1}, returns {n-1}. */

spectrum_table_binned_t spectrum_table_binned_make(uint32_t nRanges);
  /* Builds a `binned' power spectrum table {tb} with exactly {nRanges}
    consecutive frequency ranges, spanning the interval from frequency 0
    to frequency {sqrt(0.5)} (waves per pixel).

    The ranges {tb.e[i].fmin,tb.e[i].fmax} are dimensioned so that each
    will contain approximately the same number of Hartley frequency
    terms, in the limit of a very large image. The mean frequency
    {tb.e[i].fmid} is the approximate median of those terms. The fields
    {nTerms} and {power} of all entries are initialized with zero. */

void spectrum_table_binned_add_all
  ( float_image_t *P,
    int32_t c,
    bool_t center,
    spectrum_table_binned_t *tb,
    bool_t verbose
  );
  /* Accumulates the terms of channel {c} of the Hartley power
    spectrum {P} onto a binned table {tb}. If {verbose} is TRUE,
    prints diagnostics to {stderr}.
    
    The procedure assumes that {P} was obtained from the Hartley transform
    {H} of an image as with {float_image_hartley_spectrum(H,P,center)}.

    Each entry of {P} with natural frequency vector {(fX,fY)} is
    splatted onto {tb}, assuming that it is smoothly spread over the
    range of real frequency vectors {(fX±1,fY±1)}, with fold-over at
    the ends. */
  
void spectrum_table_binned_add_term
  ( spectrum_table_binned_t *tb,
    int32_t fn[], 
    int32_t fd[],
    double nTerms,
    double power,
    bool_t verbose
  );
  /* Accumulates onto the table {tb} a set of {nTerms} terms of an
    Hartley power spectrum, whose total power is {power}. If {verbose}
    is TRUE, prints diagnostics to {stderr}.

    Assumes that all terms have natural frequency vector
    {(fn[0]/fd[0],fn[1]/fd[1])}, in waves per pixel, apart from
    signs. The values of {nTerms} and {power} are splatted onto {tb},
    assuming that they are smoothly spread over the range of real
    frequency vectors {((fn[0]±1)/fd[0],(fn[1]±1)/fd[1])}, with
    fold-over at the ends. */
 
#endif
