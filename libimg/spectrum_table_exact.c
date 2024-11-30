/* See spectrum_table_exact.h */
/* Last edited on 2024-11-06 10:04:37 by stolfi */ 

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
    int32_t c,
    bool_t center,
    spectrum_table_exact_t *tx, 
    bool_t verbose
  )
  {
    int32_t cols = (int32_t)P->sz[1];
    int32_t rows = (int32_t)P->sz[2];
    int32_t fd[2] = { cols, rows };  /* Num of freqs in each axis. */
    int32_t fMin[2] = { -(cols-1)/2, -(rows-1)/2 }; /* min integer freqs. */
    int32_t fMax[2] = { cols/2, rows/2 }; /* Max integer freqs. */
    int32_t ntx = tx->ne; /* Table entries in use are {tx.e[0..ntx-1]}. */
    for (int32_t fx = fMin[0]; fx <= fMax[0]; fx++)
      { for (int32_t fy = fMin[1]; fy <= fMax[1]; fy++)
          { int32_t fn[2] = { fx, fy };
            /* Compute negated freq vector {fn}: */
            int32_t fc[2];  /* Denominators of int32_t freq vectors. */
            for (uint32_t j = 0;  j <= 1; j++) 
              { fc[j] = -fn[j];  if (fc[j] < fMin[j]) { fc[j] += fd[j]; } }
            /* We must consider {fn} only if it is leq {fc} in lex order: */
            if ((fn[0] < fc[0]) || ((fn[0] == fc[0]) && (fn[1] <= fc[1])))
              { double power = 0;
                int32_t nTerms = 0;

                auto void add_term(int32_t f[]);
                  /* Adds the Hartley spectrum term with indices {fx,fy} to {nTerms,power}. */
                add_term(fn);
                
                /* If {fn} is different from {fc}, also add {fc}: */
                if ((fc[0] != fn[0]) || (fc[1] != fn[1])) { add_term(fc); }

                void add_term(int32_t f[])
                  { int32_t rx = (center ? f[0] - fMin[0] : (f[0] + cols) % cols);
                    int32_t ry = (center ? f[1] - fMin[1] : (f[1] + rows) % rows);
                    if (verbose) { fprintf(stderr, "    adding %2d %4d %4d  fn = %+5d %+5d\n", c, rx, ry, f[0], f[1]); }

                    double val = float_image_get_sample(P, c, rx, ry);
                    if (verbose) { fprintf(stderr, "  val = %18.10f\n", val); }
                    nTerms += 1;
                    power += val;
                  }
            
                /* Append term to spectrum table: */
                spectrum_table_exact_append_term(tx, &ntx, fn, fd, nTerms, power);
              }
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
    bool_t debug = FALSE;
    int32_t ntx = (*ntxp);
    spectrum_table_exact_expand(tx, ntx);
    spectrum_table_exact_entry_t *txn = &(tx->e[ntx]);
    urat64_t freq2 = spectrum_table_exact_compute_freq2(fn, fd);
    txn->freq2 = freq2;
    txn->nTerms = nTerms;
    txn->power = power;
    if (debug)
      { fprintf(stderr, "      appended term with f = (%+5d/%-4d,%+5d/%-4d)", fn[0], fd[0], fn[1],fd[1]);
        fprintf(stderr, " |f|^2 = %lu/%lu", freq2.num, freq2.den);
        fprintf(stderr, " n = %.3f pwr = %20.12f\n", nTerms, power);
      }
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
        if (Fn > FnMax) { Fn = Fd - Fn; }
        assert((Fn >= 0) && (Fn <= FnMax));
        /* Pack as a {2×64} bit rational number: */
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
  
