/* See spectrum_table_exact.h */
/* Last edited on 2023-03-19 13:43:13 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <jsmath.h>
#include <affirm.h>
#include <bool.h>
#include <vec.h>
#include <float_image.h>

#include <spectrum_table_exact.h>

/* INTERNAL PROTOTYPES */

urat64_t spectrum_table_exact_compute_freq2(int32_t fn[], int32_t fd[]);
  /* Computes the absolute natural frequency squared {freq2} given the
    natural frequencies {fn[0]/fd[0], fn[1]/fd[1]} per axis (all in
    waves per pixel). */

/* IMPLEMENTATIONS */

vec_typeimpl(spectrum_table_exact_t,spectrum_table_exact,spectrum_table_exact_entry_t); 

void spectrum_table_exact_append_all
  ( float_image_t *P, 
    bool_t center,
    int32_t c,
    spectrum_table_exact_t *tx, 
    bool_t verbose
  )
  { /* To reduce the cost of {spectrum_table_exact_sort}, we try to 
      generate the entries in approx sorting order, and pre-combine
      entries with the same absolute natural frequency. 
      
      Specifically, we enumerate the non-negative natural frequency
      vectors {fn[0],fn[1]} in order of increasing sum, and, for each
      vector, we collet and combine all terms that differ from those
      only in sign. This strategy enumerates the terms in order of
      increasing L1 norm of their integer frequency vectors {fn},
      where {fn[0]} is the number of wave cycles per image row, and
      {fn[1]} is the number of wave cycles per image column.
      
      Since {spectrum_table_exact_sort} uses the L2 metric on the
      fractional frequency vector {fn[0]/fd[0],fn[1]/fd[1]} (waves
      per pixel), this strategy is only an heuristic, even for square
      images, and becomes less effective as the aspect ratio moves
      away from 1:1. In any case, it usually saves 75% of the entries
      that would be used by the naive approach (one entry for each
      term of the spectrum). */

    int32_t cols = (int32_t)P->sz[1];
    int32_t rows = (int32_t)P->sz[2];
    int32_t fd[2] = { cols, rows };  /* Denominators of int32_t freq vectors. */
    int32_t fnMax[2] = { fd[0]/2, fd[1]/2 }; /* Max value of numerators {fn[0],fn[1]}. */
    int32_t s; /* Sum of numerators {fn[0]+fn[1]}. */
    int32_t sMax = fnMax[0] + fnMax[1];  /* Max value of {s}. */
    int32_t ntx = tx->ne; /* Table entries in use are {tx.e[0..ntx-1]}. */
    for (s = 0; s <= sMax; s++)
      { int32_t d; /* Difference of numerators {fn[0] - fn[1]}. */
        int32_t dMin = (int32_t)imax(s - 2*fnMax[1], -s);
        int32_t dMax = (int32_t)imin(2*fnMax[0] - s, +s);
        for (d = dMin; d <= dMax; d += 2)
          { /* Compute the numerators {fxp,fyp} from {s,d}: */
            assert((s + d) % 2 == 0);
            
            /* Main frequency vector corodinates: */
            int32_t fxp = (s + d)/2; assert((fxp >= 0) && (fxp <= fnMax[0]));
            int32_t fyp = (s - d)/2; assert((fyp >= 0) && (fyp <= fnMax[1]));
            
            /* Mirrored frequency vector corodinates: */
            int32_t fxn = (fd[0] - fxp) % fd[0]; /* {-fxp mod fd[0]} */
            int32_t fyn = (fd[1] - fyp) % fd[1]; /* {-fyp mod fd[1]} */ 
            
            /* Add the terms {fxn,fxp}�{fyn,fyp} but avoiding duplicates: */
            double nTerms = 0;
            double power = 0; 
            
            auto void add_term(int32_t fx, int32_t fy);
              /* Adds the Hartley spectrum term with indices {fx,fy} to {nTerms,power}. */
              
            void add_term(int32_t fx, int32_t fy)
              { int32_t rx = (center ? (fx + cols/2) % cols : fx);
                int32_t ry = (center ? (fy + rows/2) % rows : fy);

                double val = float_image_get_sample(P, c, rx, ry);
                if (verbose) { fprintf(stderr, "    adding %2d %4d %4d = %18.10f\n", c, fx, fy, val); }
                nTerms += 1;
                power += val;
              }
              
            /* Enumerate the four (or less) terms with freqs {�fn[0],�fn[1]}: */
            
            add_term(fxp, fyp);
            if (fxn != fxp) { add_term(fxn, fyp); }
            if (fyn != fyp) { add_term(fxp, fyn); }
            if ((fxn != fxp) && (fyn != fyp)) { add_term(fxn, fyn); }
            
            /* Append term to spectrum table: */
            int32_t fn[2] = { fxp, fyp };
            spectrum_table_exact_append_term(tx, &ntx, fn, fd, nTerms, power);
          }
      }

    spectrum_table_exact_trim(tx, ntx);
  }
  
void spectrum_table_exact_append_term
  ( spectrum_table_exact_t *tx, 
    int32_t *ntxp,
    int32_t fn[], 
    int32_t fd[],
    double nTerms,
    double power
  )
  {
    int32_t ntx = (*ntxp);
    spectrum_table_exact_expand(tx, ntx);
    spectrum_table_exact_entry_t *txn = &(tx->e[ntx]);
    txn->freq2 = spectrum_table_exact_compute_freq2(fn, fd);
    txn->nTerms = nTerms;
    txn->power = power;
    ntx++;
    (*ntxp) = ntx;
  }

urat64_t spectrum_table_exact_compute_freq2(int32_t fn[], int32_t fd[])
  { 
    urat64_t F2[2]; /* {(fn[i]/fd[i])^2}, after reductions. */
    int32_t i;
    for (i = 0; i < 2; i++)
      { /* Convert axial frequency {fn[i]/fd[i]} to 64 bits: */
        int64_t Fn = fn[i]; 
        int64_t Fd = fd[i];
        demand(Fd > 0, "zero denominator in frequency vector");
        /* Reduce the numerator {Fn} to the range {0..Fd-1}: */
        Fn = ((Fn % Fd) + Fd) % Fd;
        assert((Fn >= 0) && (Fn < Fd));
        /* Reduce the numerator {Fn} to the natural range, and take the abs value: */
        int64_t FnMax = Fd/2;
        if (Fn > FnMax) { Fd = Fd - Fn; }
        assert((Fn >= 0) && (Fn <= FnMax));
        /* Pack as a {2�64} bit rational number: */
        urat64_t F = (urat64_t){ .num = Fn, .den = Fd };
        /* Remove common factors in fraction {Fn/Fd}: */
        urat64_reduce(&F);
        /* Set {F2[i]} to {F} squared (safe, since {fn,fd} were 32 bit signed ints). */
        urat64_sqr(&F, &(F2[i]));
        /* Check for possible overflow in following computations: */
        demand(F2[i].num < (1LLU<<31), "freq2 numerator is too large");
        demand(F2[i].den < (1LLU<<31), "freq2 denominator is too large");
      }
    /* Compute the sum {freq2 = F2[0] + F2[1]}: */
    urat64_t freq2;
    freq2.den = lcm(F2[0].den, F2[1].den); /* Less than {1<<62}. */
    freq2.num = 
      (freq2.den/F2[0].den)*F2[0].num + 
      (freq2.den/F2[1].den)*F2[1].num; /* Less than {1<<63}. */
    return freq2;
  }

void spectrum_table_exact_sort(spectrum_table_exact_t *tx, bool_t verbose)
  { 
  
    /* Sort the table by increasing term frequency {freq}: */
    auto int32_t cmp_freq(const void *a, const void *b);
    int32_t cmp_freq(const void *a, const void *b) 
      { /* Get hold the two entries: */
        spectrum_table_exact_entry_t *tba = (spectrum_table_exact_entry_t *)a;
        spectrum_table_exact_entry_t *tbb = (spectrum_table_exact_entry_t *)b;
        return urat64_compare(&(tba->freq2), &(tbb->freq2));
      }
    qsort(tx->e, tx->ne, sizeof(spectrum_table_exact_entry_t), &cmp_freq);
    
    /* Merge entries with identical frequencies: */
    int32_t nc = 0; /* The condensed entries are {tb.e[0..nc-1]}. */
    int32_t k;
    for (k = 0; k < tx->ne; k++)
      { if (nc == 0)
          { /* First entry -- leave it there: */
            assert(nc == k); nc++;
          }
        else 
          { /* Grab the next uncondensed entry {*ek}: */
            spectrum_table_exact_entry_t *ek = &(tx->e[k]);
            /* Grab the previous condensed term {*ec}: */
            spectrum_table_exact_entry_t *ec = &(tx->e[nc-1]);
            int32_t ord = cmp_freq(ek, ec);
            if (ord < 0)
              { /* Sorting bug? */ assert(FALSE); }
            else if (ord > 0)
              { /* A new frequency, append to table: */
                tx->e[nc] = (*ek);
                nc++;
              }
            else
              { /* Seems to be the same frequency as previous entry, condense: */
                ec->nTerms += ek->nTerms;
                ec->power += ek->power;
              }
          }
      }
      
    /* Trim the table to the condensed entries only: */
    spectrum_table_exact_trim(tx, nc);
  }
  
