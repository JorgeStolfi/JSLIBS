/* See ppv_array_write.h */
/* Last edited on 2021-07-10 04:28:01 by jstolfi */
/* Copyright © 2003 by Jorge Stolfi, from University of Campinas, Brazil. */
/* See the rights and conditions notice at the end of this file. */

#include <ppv_array.h>
#include <ppv_array_write.h>

#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <bool.h>
#include <affirm.h>

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* INTERNAL PROTOTPES */

void ppv_array_write_samples_plain ( FILE *wr, ppv_array_t *A );
void ppv_array_write_samples_raw_small ( FILE *wr, ppv_array_t *A );
void ppv_array_write_samples_raw_big ( FILE *wr, ppv_array_t *A );
void ppv_array_write_samples_raw_bytes ( FILE *wr, ppv_array_t *A );

#define ppv_FILE_TYPE "ppv_array_t"
 
#define smpFMT ppv_sample_t_FMT

void ppv_array_write_file(FILE *wr, ppv_array_t *A, bool_t plain)
  {
    /* Find effective number of axes: */
    ppv_dim_t d = A->d;
    
    /* Write header: */
    ppv_axis_t i;
    filefmt_write_header(wr, ppv_FILE_TYPE, ppv_FILE_VERSION_2021);
    
    /* Write number of effective indices: */
    fprintf(wr, "dim = %d\n", d);
    
    /* Write nominal sizes for effective indices, determine if empty: */
    fprintf(wr, "size =");
    bool_t empty = FALSE;
    for (i = 0; i < d; i++) 
      { fputc(' ', wr);
        fprintf(wr, ppv_size_t_FMT, A->size[i]);
        if (A->size[i] == 0) { empty = TRUE; }
      }
    fprintf(wr, "\n");
    
    /* Write virtual sizes, temporarily un-virtualize {A}, save original size: */
    ppv_size_t osz[d];
    fprintf(wr, "asize =");
    for (i = 0; i < d; i++) 
      { ppv_size_t szi = A->size[i];
        osz[i] = szi;
        if ((! empty) && (szi >= 2) && (A->step[i] == 0)) { ppv_crop(A, i, 0, 1); }
        fputc(' ', wr);
        fprintf(wr, ppv_size_t_FMT, A->size[i]);
      }
    fprintf(wr, "\n");
    
    /* Write max sample value: */
    fprintf(wr, "maxsmp = " ppv_sample_t_FMT "\n", A->maxsmp);
    
    /* Write plain/raw flag: */
    fprintf(wr, "plain = %d\n", (plain ? 1 : 0));

    /* Now write samples: */
    if ((A->maxsmp > 0) && (! empty))
      { if (plain)
          { ppv_array_write_samples_plain(wr, A); }
        else if (A->bps < 8)
          { ppv_array_write_samples_raw_small(wr, A); }
        else if (A->bps > 8)
          { ppv_array_write_samples_raw_big(wr, A); }
      else /* if (A->bps == 8) */
          { ppv_array_write_samples_raw_bytes(wr, A); }     
        fputc('\n', wr);
      }
    
    /* Restore virtual replications: */
    for (i = 0; i < d; i++) 
      { ppv_size_t szi = osz[i];
        if (szi > A->size[i])
          { assert(A->size[i] == 1);
            ppv_replicate(A, i, szi);
          }
      }

    /* Write footer: */
    filefmt_write_footer(wr, ppv_FILE_TYPE);

    fflush(wr);
  }

#define ppv_MAX_SAMPLES_PER_LINE 20

void ppv_array_write_samples_plain ( FILE *wr, ppv_array_t *A )
  {
    ppv_dim_t d = A->d;
    ppv_index_t ix[d];
    assert(A->maxsmp > 0);
    if (ppv_index_first(ix, A)) 
      { int32_t n = 0; /* Number of samples in current line. */
        ppv_pos_t pos = ppv_sample_pos(A, ix); /* Position of sample {A[ix]} */
        do 
          { ppv_sample_t smp = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pos);
            /* Space/newline control: */
            if (n >= ppv_MAX_SAMPLES_PER_LINE)
              { /* Line too long: */
                fputc('\n', wr); n = 0;
              }
            else if ((n > 0) && (d >= 2) && (ix[0] == 0))
              { /* Start of a new row: */
                fputc('\n', wr); n = 0;
              }
            else if (n > 0)
              { /* Insert a separating space: */
                fputc(' ', wr);
              }
            fprintf(wr, ppv_sample_t_FMT, smp); n++;
          } 
        while (! ppv_index_next(ix, A, d, &pos));
      }
  }

void ppv_array_write_samples_raw_small ( FILE *wr, ppv_array_t *A )
  {
    ppv_dim_t d = A->d;
    ppv_index_t ix[d];
    assert(A->maxsmp > 0);
    if (ppv_index_first(ix, A)) 
      { int32_t spc = 8/A->bps; /* Samples per byte. */
        ppv_word_08_t buf = 0;
        ppv_nbits_t shift_ini = (ppv_nbits_t)((spc - 1)*A->bps); /* Shift to apply to sample [0]. */
        ppv_nbits_t shift = shift_ini; /* Shift to apply to next sample. */
        ppv_pos_t pos = ppv_sample_pos(A, ix);
        do 
          { ppv_sample_t smp = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pos);
            buf = (ppv_word_08_t)(buf | (smp << shift));
            if (shift == 0) 
              { fputc((char)buf, wr); shift = shift_ini; buf = 0; }
            else
              { shift = (ppv_nbits_t)(shift - A->bps); }
          } 
        while (! ppv_index_next(ix, A, d, &pos));
        if (shift < shift_ini) { fputc((char)buf, wr); }
      }
  }

void ppv_array_write_samples_raw_big ( FILE *wr, ppv_array_t *A )
  {
    ppv_dim_t d = A->d;
    ppv_index_t ix[d];
    assert(A->maxsmp > 0);
    if (ppv_index_first(ix, A)) 
      { ppv_nbits_t cps = (ppv_nbits_t)((A->bps + 7)/8); /* Bytes per sample. */
        ppv_word_t maskw = 255; /* Mask to chop a byteful from the sample. */
        ppv_nbits_t shift_ini = (ppv_nbits_t)((cps-1)*8); /* Shift to apply to sample to get first byte. */
        ppv_pos_t pos = ppv_sample_pos(A, ix);
        do 
          { ppv_sample_t smp = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pos);
            assert(smp <= A->maxsmp);
            /* bool_t debug = ((A->bps == 10) && (pos == 0)); */
            /* if (debug) { fprintf(stderr,  "smp = " smpFMT "\n", smp); } */
            ppv_nbits_t shift = shift_ini; /* Shift to get next byte. */
            for (int32_t k = 0; k < cps; k++) 
              { ppv_word_08_t ch = (ppv_word_08_t)((smp >> shift) & maskw); 
                /* if (debug) { fprintf(stderr,  "  ch = %u\n", ch); } */
                fputc((char)ch, wr); 
                shift = (ppv_nbits_t)(shift - 8);
              }
          } 
        while (! ppv_index_next(ix, A, d, &pos));
      }
  }
                  
void ppv_array_write_samples_raw_bytes ( FILE *wr, ppv_array_t *A )
  {
    ppv_dim_t d = A->d;
    ppv_index_t ix[d];
    assert(A->maxsmp > 0);
    if (ppv_index_first(ix, A)) 
      { ppv_pos_t pos = ppv_sample_pos(A, ix);
        do 
          { ppv_sample_t smp = ppv_get_sample_at_pos(A->el, A->bps, A->bpw, pos);
            fputc((char)smp, wr);
          }                  
        while (! ppv_index_next(ix, A, d, &pos));
      }
  }

