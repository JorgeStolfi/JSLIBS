/* See {spectrum_table_convert.h} */
/* Last edited on 2024-12-05 07:14:16 by stolfi */ 

#include <limits.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <jsmath.h>
#include <affirm.h>
#include <bool.h>
#include <urat64.h>
#include <vec.h>

#include <spectrum_table_binned.h>
#include <spectrum_table_exact.h>
#include <spectrum_table_convert.h>

/* INTERNAL PROTOTYPES */

spectrum_table_binned_t spectrum_table_convert_exact_to_binned
  ( spectrum_table_exact_t *tx,
    uint32_t cols, 
    uint32_t rows )
  {
    /* Create a binned table with roughly correct size: */
    uint32_t size_guess = cols*rows/2;
    spectrum_table_binned_t tb = spectrum_table_binned_new(size_guess);
    
    /* Scan the exact table and join/fill entries: */
    uint32_t nb = 0; /* Valid entries of {tb} are {tb.e[0..ib-1]}. */
    
    auto void append_proper_entry(float freq, double nTerms, double power);
      /* Appends to {tb} a proper entry {(freq,nTerms,power)},
        provided that {nTerms} is not zero. */
    
    void append_proper_entry(float freq, double nTerms, double power)
      { if (nTerms != 0.0)
          { spectrum_table_binned_expand(&tb, nb);
            spectrum_table_binned_entry_t *eb = &(tb.e[nb]);
            eb->fmin = eb->fmax = eb->fmid = (double)freq;
            eb->nTerms = nTerms;
            eb->power = power;
            nb++; 
          }
      }
    
    auto void append_filler_entry(float flo, float fhi);
    /* Appends to {tb} a filler entry spanning from {flo} to {fhi},
      if they are different. */
    
    void append_filler_entry(float flo, float fhi)
      { assert(flo <= fhi);
        if (flo < fhi)
          { spectrum_table_binned_expand(&tb, nb);
            spectrum_table_binned_entry_t *eb = &(tb.e[nb]);
            eb->fmin = (double)flo;
            eb->fmax = (double)fhi;
            float fmd = (float)fmax(flo, fmin(fhi, (flo + fhi)/2));
            eb->fmid = (double)fmd;
            eb->nTerms = 0;
            eb->power = 0.0;
            nb++;
          }
      }
    
    /* During processing we have an incomplete proper entry for {tb}, not yet saved: */
    urat64_t f2prev = urat64_ZERO; /* Exact {freq2} of previous {tx} entry. */
    float fprev = 0.0; /* Approx frequency of incomplete entry of {tb}. */
    double nprev = 0.0;  /* Hartley term count of the incomplete entry. */
    double pprev = 0.0;  /* Total {power} of incomplete entry. */
    float flimit = (float)M_SQRT1_2; /* Max float frequency, {sqrt(1/2)}. */
    /* Scan entries of {tx} (including a dummy one at the end): */
    int32_t ix; 
    for (ix = 0; ix <= tx->ne; ix++)
      { /* Get the data {fthis,nthis,pthis} of the next entry of {tx}: */
        float fthis; 
        double nthis, pthis;
        if (ix < tx->ne)
          { /* Get the next real entry of {tx}: */
            spectrum_table_exact_entry_t *ex = &(tx->e[ix]);
            /* The table {tx} must be sorted by {.freq2}: */
            assert(urat64_compare(&f2prev, &(ex->freq2)) <= 0);
            /* Compute the frequency of {ex} as a {float}: */
            double f2num = (double)(ex->freq2.num);
            double f2den = (double)(ex->freq2.den);
            fthis = (float)sqrt(f2num/f2den); /* !!! Is this adequare, roundingwise? !!! */
            /* Make sure that {fthis} does not exceed the maximum possible frequency: */
            if (fthis > flimit) { fthis = flimit; }
            /* Grab the entry's term count and total power: */
            nthis = ex->nTerms;
            pthis = ex->power;
          }
        else
          { /* Fake a dummy entry at the end of {tx}: */
            fthis = flimit; nthis = 0; pthis = 0.0;
          }
        /* Should we lump this entry with the incomplete one? */
        if (fthis <= fprev)
          { /* They have the same frequency within {float} accuracy, lump it: */
            nprev += nthis; pprev += pthis;
          }
        else 
          { /* The new entry of {tx} has a higher frequency, we cannot lump it. */
            /* The incomplete entry is now complete, append it to {tb}: */
            append_proper_entry(fprev, nprev, pprev);
            /* Append to {tb} a filler entry to span the gap: */
            append_filler_entry(fprev, fthis);
            /* Start a new incomplete entry: */
            fprev = fthis; nprev = nthis; pprev = pthis;
          }
      }
    /* The incomplete entry now must have {fprev = flimit}: */
    assert(fprev == flimit);
    /* Append it to {tb}: */
    append_proper_entry(fprev, nprev, pprev);

    /* Now trim{tb} to size and return it: */
    spectrum_table_binned_trim(&tb, nb);
    return tb;
  }
