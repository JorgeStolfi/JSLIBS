/* See ppv_io.h */
/* Last edited on 2020-12-07 23:20:37 by jstolfi */
/* Copyright © 2003 by Jorge Stolfi, from University of Campinas, Brazil. */
/* See the rights and conditions notice at the end of this file. */

#include <ppv_array.h>
#include <ppv_io.h>

#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <bool.h>
#include <affirm.h>

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* INTERNAL PROTOTPES */

void ppv_write_pixels_plain ( FILE *wr, ppv_array_t *A, ppv_dim_t N );
void ppv_write_pixels_raw_small ( FILE *wr, ppv_array_t *A, ppv_dim_t N );
void ppv_write_pixels_raw_big ( FILE *wr, ppv_array_t *A, ppv_dim_t N );
void ppv_write_pixels_raw_bytes ( FILE *wr, ppv_array_t *A, ppv_dim_t N );

void ppv_read_pixels_plain ( FILE *rd, ppv_array_t *A, ppv_dim_t N );
void ppv_read_pixels_raw_small ( FILE *rd, ppv_array_t *A, ppv_dim_t N );
void ppv_read_pixels_raw_big ( FILE *rd, ppv_array_t *A, ppv_dim_t N );
void ppv_read_pixels_raw_bytes ( FILE *rd, ppv_array_t *A, ppv_dim_t N );

#define ppv_FILE_TYPE "ppv_array_t"
#define ppv_FILE_VERSION "2005-05-28"
    
void ppv_write_array(FILE *wr, ppv_array_t *A, bool_t plain)
  {
    /* Find effective number of axes: */
    ppv_dim_t N = ppv_array_NAXES;
    while ((N > 0) && (A->size[N-1] == 1)) { N--; }
    
    /* Write header: */
    ppv_axis_t i;
    filefmt_write_header(wr, ppv_FILE_TYPE, ppv_FILE_VERSION);
    
    /* Write number of effective indices: */
    fprintf(wr, "dim = %d\n", N);
    
    /* Write nominal sizes for effective indices, determine if empty: */
    fprintf(wr, "size =");
    bool_t empty = FALSE;
    for (i = 0; i < N; i++) 
      { fputc(' ', wr);
        fprintf(wr, ppv_size_t_FMT, A->size[i]);
        if (A->size[i] == 0) { empty = TRUE; }
      }
    fprintf(wr, "\n");
    
    /* Write virtual sizes, temporarily un-virtualize {A}, save original size: */
    ppv_size_t osz[ppv_array_NAXES];
    fprintf(wr, "asize =");
    for (i = 0; i < N; i++) 
      { ppv_size_t szi = A->size[i];
        osz[i] = szi;
        if ((! empty) && (szi >= 2) && (A->step[i] == 0)) { ppv_crop(A, i, 0, 1); }
        fputc(' ', wr);
        fprintf(wr, ppv_size_t_FMT, A->size[i]);
      }
    fprintf(wr, "\n");
    
    /* Write bits per sample: */
    fprintf(wr, "bps = %d\n", A->bps);
    
    /* Write plain/raw flag: */
    fprintf(wr, "plain = %d\n", (plain ? 1 : 0));

    /* Now write samples: */
    if (! empty)
      { if (plain)
          { ppv_write_pixels_plain(wr, A, N); }
        else if (A->bps < 8)
          { ppv_write_pixels_raw_small(wr, A, N); }
        else if (A->bps > 8)
          { ppv_write_pixels_raw_big(wr, A, N); }
      else /* if (A->bps == 8) */
          { ppv_write_pixels_raw_bytes(wr, A, N); }     
        fputc('\n', wr);
      }
    
    /* Restore virtual replications: */
    for (i = 0; i < N; i++) 
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

void ppv_write_pixels_plain ( FILE *wr, ppv_array_t *A, ppv_dim_t N )
  {
    ppv_index_t ix[ppv_array_NAXES];
    if (ppv_index_first(ix, A)) 
      { int n = 0; /* Number of samples in current line. */
        ppv_pos_t pos = ppv_sample_pos(A, ix); /* Position of sample {A[ix]} */
        do 
          { ppv_sample_t qv = ppv_get_sample(A->el, A->bps, A->bpw, pos);
            /* Space/newline control: */
            if (n >= ppv_MAX_SAMPLES_PER_LINE)
              { /* Line too long: */
                fputc('\n', wr); n = 0;
              }
            else if ((n > 0) && (N >= 2) && (ix[0] == 0))
              { /* Start of a new row: */
                fputc('\n', wr); n = 0;
              }
            else if (n > 0)
              { /* Insert a separating space: */
                fputc(' ', wr);
              }
            fprintf(wr, ppv_sample_t_FMT, qv); n++;
          } 
        while (! ppv_index_next(ix, A, N, &pos));
      }
  }

void ppv_write_pixels_raw_small ( FILE *wr, ppv_array_t *A, ppv_dim_t N )
  {
    ppv_index_t ix[ppv_array_NAXES];
    if (ppv_index_first(ix, A)) 
      { int spc = 8/A->bps; /* Samples per byte. */
        ppv_word_08_t buf = 0;
        ppv_nbits_t shift_ini = (ppv_nbits_t)((spc - 1)*A->bps); /* Shift to apply to sample [0]. */
        ppv_nbits_t shift = shift_ini; /* Shift to apply to next sample. */
        ppv_pos_t pos = ppv_sample_pos(A, ix);
        do 
          { ppv_sample_t qv = ppv_get_sample(A->el, A->bps, A->bpw, pos);
            buf = (ppv_word_08_t)(buf | (qv << shift));
            if (shift == 0) 
              { fputc((char)buf, wr); shift = shift_ini; buf = 0; }
            else
              { shift = (ppv_nbits_t)(shift - A->bps); }
          } 
        while (! ppv_index_next(ix, A, N, &pos));
        if (shift < shift_ini) { fputc((char)buf, wr); }
      }
  }

void ppv_write_pixels_raw_big ( FILE *wr, ppv_array_t *A, ppv_dim_t N )
  {
    ppv_index_t ix[ppv_array_NAXES];
    if (ppv_index_first(ix, A)) 
      { ppv_nbits_t cps = (ppv_nbits_t)((A->bps + 7)/8); /* Bytes per sample. */
        ppv_word_t maskw = 255; /* Mask to chop a byteful from the sample. */
        ppv_nbits_t shift_ini = (ppv_nbits_t)((cps-1)*8); /* Shift to apply to sample to get first byte. */
        ppv_word_08_t buf;
        ppv_nbits_t shift = shift_ini; /* Shift to get next byte. */
        ppv_pos_t pos = ppv_sample_pos(A, ix);
        do 
          { ppv_sample_t qv = ppv_get_sample(A->el, A->bps, A->bpw, pos);
            int k;
            for (k = 1; k < cps; k++) 
              { buf = (ppv_word_08_t)((qv >> shift) & maskw); 
                fputc((char)buf, wr); 
                shift = (ppv_nbits_t)(shift - A->bps);
              }
          } 
        while (! ppv_index_next(ix, A, N, &pos));
      }
  }
                  
void ppv_write_pixels_raw_bytes ( FILE *wr, ppv_array_t *A, ppv_dim_t N )
  {
    ppv_index_t ix[ppv_array_NAXES];
    if (ppv_index_first(ix, A)) 
      { ppv_pos_t pos = ppv_sample_pos(A, ix);
        do 
          { ppv_sample_t qv = ppv_get_sample(A->el, A->bps, A->bpw, pos);
            fputc((char)qv, wr);
          }                  
        while (! ppv_index_next(ix, A, N, &pos));
      }
  }

ppv_array_t ppv_read_array(FILE *rd, ppv_nbits_t bpw)
  {
    /* Read header: */
    filefmt_read_header(rd, ppv_FILE_TYPE, ppv_FILE_VERSION);
    
    /* Read effective dimensions: */
    int N = nget_int(rd, "dim");
    demand(N <= ppv_array_NAXES, "too many dimensions");
    fget_eol(rd);
    
    /* Read nominal array shape {sz}, complete with 1's, check for emptiness: */
    ppv_size_t sz[ppv_array_NAXES]; 
    bool_t empty = FALSE;
    nget_name_eq(rd, "size");
    ppv_axis_t i;
    for (i = 0; i < ppv_array_NAXES; i++) 
      { if (i < N)
          { uint64_t szi = fget_uint64(rd, 10);
            demand(szi <= ppv_MAX_SIZE, "size too big");
            sz[i] = szi;
            if (szi == 0) { empty = TRUE; }
          }
        else
          { sz[i] = 1; }
      }
    fget_eol(rd); 
    
    /* Read actual sizes {asz}: */
    ppv_size_t asz[ppv_array_NAXES]; 
    nget_name_eq(rd, "asize");
    for (i = 0; i < ppv_array_NAXES; i++) 
      { if (i < N)
          { uint64_t aszi = fget_uint64(rd, 10);
            if (aszi != sz[i])
              { demand((aszi == 1) && (sz[i] >= 2), "invalid actual size"); }
            asz[i] = aszi;
          }
        else
          { asz[i] = 1; }
      }
    fget_eol(rd); 
    
    /* Read bits-per-sample: */
    uint64_t bps = nget_uint64(rd, "bps", 10);
    demand(bps <= ppv_MAX_BPS, "too many bits per sample");
    fget_eol(rd); 
    
    /* Read plain/raw flag: */
    bool_t plain = nget_bool(rd, "plain");
    fget_eol(rd); 
    
    /* Allocate array using {asz}: */
    ppv_array_t A = ppv_new_array(asz, (ppv_nbits_t)bps, bpw);
    
    if (! empty)
      { 
        /* Now read samples: */
        if (plain)
          { ppv_read_pixels_plain(rd, &A, (ppv_dim_t)N); }
        else if (A.bps < 8)
          { ppv_read_pixels_raw_small(rd, &A, (ppv_dim_t)N); }
        else if (A.bps > 8)
          { ppv_read_pixels_raw_big(rd, &A, (ppv_dim_t)N); }
        else /* A.bps == 8 */
          { ppv_read_pixels_raw_bytes(rd, &A, (ppv_dim_t)N); }     
        fget_skip_formatting_chars(rd);
      }

    /* Apply virtual replication: */
    for (i = 0; i < N; i++) 
      { if (sz[i] > asz[i]) 
          { assert(A.size[i] == 1); ppv_replicate(&A, i, sz[i]); }
      }
    
    /* Read footer: */
    filefmt_read_footer(rd, ppv_FILE_TYPE);
    return A;
  }

void ppv_read_pixels_plain ( FILE *rd, ppv_array_t *A, ppv_dim_t N )
  {
    bool_t debug = FALSE;
    ppv_index_t ix[ppv_array_NAXES];
    ppv_sample_t max_sample = (ppv_sample_t)((1LLU << A->bps) - 1);
    if (ppv_index_first(ix, A)) 
      { ppv_pos_t pos = ppv_sample_pos(A, ix);
        do 
          { fget_skip_formatting_chars(rd);
            uint64_t qv = fget_uint64(rd, 10);
            if (debug) 
              { fprintf(stderr, "    ix =");
                int j;
                for (j = 0; j < N; j++) { fprintf(stderr, " %ld", ix[j]); }
                fprintf(stderr, "  smp = %lu\n", (uint64_t)qv);
              }
            demand(qv <= max_sample, "sample too big");
            ppv_set_sample(A->el, A->bps, A->bpw, pos, (ppv_sample_t)qv);
         } 
        while (! ppv_index_next(ix, A, N, &pos));
      }
  }

void ppv_read_pixels_raw_small ( FILE *rd, ppv_array_t *A, ppv_dim_t N )
  {
    
    
    ppv_index_t ix[ppv_array_NAXES];
    ppv_sample_t smask = (ppv_sample_t)((1LLU << A->bps) - 1);
    if (ppv_index_first(ix, A)) 
      { int spc = 8/A->bps; /* Samples per byte. */
        ppv_nbits_t shift_ini = (ppv_nbits_t)((spc - 1)*A->bps); /* Shift to apply to sample [0]. */
        ppv_nbits_t shift = 0; /* Shift applied to last sample. */
        ppv_pos_t pos = ppv_sample_pos(A, ix);
        int buf = 0; /* Character being torn apart. */
        do 
          { if (shift == 0) 
              { buf = fgetc(rd);
                demand(buf != EOF, "unexpected end-of-file");
                shift = shift_ini;
              }
            else
              { shift = (ppv_nbits_t)(shift - A->bps); }
            ppv_set_sample(A->el, A->bps, A->bpw, pos, (buf >> shift) & smask);
          }
        while (! ppv_index_next(ix, A, N, &pos));
      }
  }

void ppv_read_pixels_raw_big ( FILE *rd, ppv_array_t *A, ppv_dim_t N )
  {
    affirm(FALSE, "not implemented yet");
  }

void ppv_read_pixels_raw_bytes ( FILE *rd, ppv_array_t *A, ppv_dim_t N )
  {
    ppv_index_t ix[ppv_array_NAXES];
    ppv_sample_t smask = (ppv_sample_t)((1LLU << A->bps) - 1);
    if (ppv_index_first(ix, A)) 
      { ppv_pos_t pos = ppv_sample_pos(A, ix);
        do 
          { int qv = fgetc(rd);
            demand(qv != EOF, "unexpected end-of-file");
            ppv_set_sample(A->el, A->bps, A->bpw, pos, qv & smask);
         } 
        while (! ppv_index_next(ix, A, N, &pos));
      }
  }
