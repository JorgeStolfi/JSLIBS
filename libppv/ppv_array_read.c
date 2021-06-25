/* See ppv_array_read.h */
/* Last edited on 2021-06-22 13:44:34 by jstolfi */
/* Copyright © 2003 by Jorge Stolfi, from University of Campinas, Brazil. */
/* See the rights and conditions notice at the end of this file. */

#include <ppv_array.h>
#include <ppv_array_read.h>

#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <bool.h>
#include <affirm.h>

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

/* INTERNAL PROTOTPES */

void ppv_array_read_samples_plain ( FILE *rd, ppv_array_t *A );
void ppv_array_read_samples_raw_small ( FILE *rd, ppv_array_t *A );
void ppv_array_read_samples_raw_big ( FILE *rd, ppv_array_t *A );
void ppv_array_read_samples_raw_bytes ( FILE *rd, ppv_array_t *A );

#define ppv_FILE_TYPE "ppv_array_t"
#define ppv_FILE_VERSION "2005-05-28"
    
ppv_array_t *ppv_array_read_file(FILE *rd, ppv_nbits_t bpw)
  {
    /* Read header: */
    filefmt_read_header(rd, ppv_FILE_TYPE, ppv_FILE_VERSION);
    
    /* Read effective dimensions: */
    uint64_t d64 = nget_uint64(rd, "dim", 10);
    demand(d64 <= ppv_MAX_DIM, "too many dimensions");
    fget_eol(rd);
    ppv_dim_t d = (ppv_dim_t)d64;
    
    /* Read nominal array shape {sz}, check for emptiness: */
    ppv_size_t sz[d]; 
    bool_t empty = FALSE;
    nget_name_eq(rd, "size");
    ppv_axis_t i;
    for (i = 0; i < d; i++) 
      { uint64_t szi = fget_uint64(rd, 10);
        demand(szi <= ppv_MAX_SIZE, "size too big");
        sz[i] = szi;
        if (szi == 0) { empty = TRUE; }
      }
    fget_eol(rd); 
    
    /* Read actual sizes {asz}: */
    ppv_size_t asz[d]; 
    nget_name_eq(rd, "asize");
    for (i = 0; i < d; i++) 
      { uint64_t aszi = fget_uint64(rd, 10);
        if (aszi != sz[i])
          { demand((aszi == 1) && (sz[i] >= 2), "invalid actual size"); }
        asz[i] = aszi;
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
    ppv_array_t *A = ppv_array_new(d, asz, (ppv_nbits_t)bps, bpw);
    
    if (! empty)
      { 
        /* Now read samples: */
        if (plain)
          { ppv_array_read_samples_plain(rd, A); }
        else if (A->bps < 8)
          { ppv_array_read_samples_raw_small(rd, A); }
        else if (A->bps > 8)
          { ppv_array_read_samples_raw_big(rd, A); }
        else /* A->bps == 8 */
          { ppv_array_read_samples_raw_bytes(rd, A); }     
        fget_skip_formatting_chars(rd);
      }

    /* Apply virtual replication: */
    for (i = 0; i < d; i++) 
      { if (sz[i] > asz[i]) 
          { assert(A->size[i] == 1); ppv_replicate(A, i, sz[i]); }
      }
    
    /* Read footer: */
    filefmt_read_footer(rd, ppv_FILE_TYPE);
    return A;
  }

void ppv_array_read_samples_plain ( FILE *rd, ppv_array_t *A )
  {
    ppv_dim_t d = A->d;
    bool_t debug = FALSE;
    ppv_index_t ix[d];
    ppv_sample_t max_sample = (ppv_sample_t)((1LLU << A->bps) - 1);
    if (ppv_index_first(ix, A)) 
      { ppv_pos_t pos = ppv_sample_pos(A, ix);
        do 
          { fget_skip_formatting_chars(rd);
            uint64_t qv = fget_uint64(rd, 10);
            if (debug) 
              { fprintf(stderr, "    ix =");
                int32_t j;
                for (j = 0; j < d; j++) { fprintf(stderr, " %ld", ix[j]); }
                fprintf(stderr, "  smp = %lu\n", (uint64_t)qv);
              }
            demand(qv <= max_sample, "sample too big");
            ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, (ppv_sample_t)qv);
         } 
        while (! ppv_index_next(ix, A, d, &pos));
      }
  }

void ppv_array_read_samples_raw_small ( FILE *rd, ppv_array_t *A )
  {
    ppv_dim_t d = A->d;
    ppv_index_t ix[d];
    ppv_sample_t smask = (ppv_sample_t)((1LLU << A->bps) - 1);
    if (ppv_index_first(ix, A)) 
      { int32_t spc = 8/A->bps; /* Samples per byte. */
        ppv_nbits_t shift_ini = (ppv_nbits_t)((spc - 1)*A->bps); /* Shift to apply to sample [0]. */
        ppv_nbits_t shift = 0; /* Shift applied to last sample. */
        ppv_pos_t pos = ppv_sample_pos(A, ix);
        int32_t buf = 0; /* Character being torn apart. */
        do 
          { if (shift == 0) 
              { buf = fgetc(rd);
                demand(buf != EOF, "unexpected end-of-file");
                shift = shift_ini;
              }
            else
              { shift = (ppv_nbits_t)(shift - A->bps); }
            ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, (buf >> shift) & smask);
          }
        while (! ppv_index_next(ix, A, d, &pos));
      }
  }

void ppv_array_read_samples_raw_big ( FILE *rd, ppv_array_t *A )
  {
    affirm(FALSE, "not implemented yet");
  }

void ppv_array_read_samples_raw_bytes ( FILE *rd, ppv_array_t *A )
  {
    ppv_dim_t d = A->d;
    ppv_index_t ix[d];
    ppv_sample_t smask = (ppv_sample_t)((1LLU << A->bps) - 1);
    if (ppv_index_first(ix, A)) 
      { ppv_pos_t pos = ppv_sample_pos(A, ix);
        do 
          { int32_t qv = fgetc(rd);
            demand(qv != EOF, "unexpected end-of-file");
            ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, qv & smask);
         } 
        while (! ppv_index_next(ix, A, d, &pos));
      }
  }
