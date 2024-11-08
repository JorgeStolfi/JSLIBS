#ifndef spectrum_table_exact_H
#define spectrum_table_exact_H

/* Tools for computing unsmoothed radial power spectra of images. */
/* Last edited on 2024-11-06 10:04:51 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <urat64.h>
#include <vec.h>
#include <float_image.h>

typedef struct spectrum_table_exact_entry_t 
  { urat64_t freq2;  /* Square of absolute natural frequency. */
    double nTerms;   /* Count of Hartley terms with this frequency. */
    double power;    /* Total power of those terms, per channel. */
  } spectrum_table_exact_entry_t;
  /* An entry of an `exact' power spectrum table.
  
    Each entry combines {nTerms} terms of the Hartley spectrum whose
    absolute natural frequency is {freq = sqrt(freq2)}. The field
    {nTerms} is usually a positive integer.
    
    We store the squared frequency {freq2} as a {2×64}-bit rational
    number, rather than {freq} a double, to avoid the complications of
    rounded math. The value of {freq2.num} lies in the range
    {0..floor(freq2.den/2)}; so the frequency {freq} lies between 0 and
    {sqrt(0.5)}, the latter value being reached only for images with
    even width and height. */

vec_typedef(spectrum_table_exact_t,spectrum_table_exact,spectrum_table_exact_entry_t); 

void spectrum_table_exact_append_all
  ( float_image_t *P,
    int32_t c,
    bool_t center,
    spectrum_table_exact_t *tx, 
    bool_t verbose
  );
  /* Appends to the table {tx} the terms of channel {c} 
    of the Hartley power spectrum {P}.
    
    The procedure assumes that {P} was obtained from the Hartley transform
    {H} of an image as with {float_image_hartley_spectrum(H,P,center)}.
    
    Terms which have exactly the same absolute natural frequency may
    be collapsed together. The table is expanded and trimmed as needed.
    If {verbose} is TRUE, the function prints
    diagnostics to {stderr}. */

void spectrum_table_exact_append_term
  ( spectrum_table_exact_t *tx, 
    int32_t *ntxp,
    int32_t fn[], 
    int32_t fd[],
    double nTerms,
    double power
  );
  /* Appends to the table {tx} an entry that represents {nTerms} terms
    of an Hartley power spectrum, whose total power is {power}.
    Assumes that the table entries in use are {tx->e[0..(*ntxp)-1]};
    increments {*ntxp} and expands {tx->e} if necessary. If {verbose}
    is TRUE, prints diagnostics to {stderr}. */

void spectrum_table_exact_sort(spectrum_table_exact_t *tx, bool_t verbose);
  /* Sorts the entries of {tx} by increasing absolute frequency
    {freq2} and condenses entries with exactly the same frequency.
    Each output entry will contain the sum of the {nTerms} and {power}
    fields of the entries that were combined into it. If {verbose} is
    TRUE, the function prints diagnostics to {stderr}. */
  
#undef MAX_CHNS
 
#endif
