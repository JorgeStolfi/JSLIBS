/* See ppv_array_read.h */
/* Last edited on 2021-07-10 04:27:33 by jstolfi */
/* Copyright © 2003 by Jorge Stolfi, from University of Campinas, Brazil. */
/* See the rights and conditions notice at the end of this file. */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>
#include <nget.h>
#include <fget.h>
#include <filefmt.h>

#include <ppv_array.h>
#include <ppv_array_read.h>

/* INTERNAL PROTOTPES */

void ppv_array_read_samples_plain ( FILE *rd, ppv_array_t *A );
void ppv_array_read_samples_raw_small ( FILE *rd, ppv_array_t *A );
void ppv_array_read_samples_raw_big ( FILE *rd, ppv_array_t *A );
void ppv_array_read_samples_raw_bytes ( FILE *rd, ppv_array_t *A );

#define ppv_FILE_TYPE "ppv_array_t"
 
#define smpFMT ppv_sample_t_FMT
    
ppv_array_t *ppv_array_read_file ( FILE *rd )
  {
    /* Read header: */
    char *fileType = NULL;
    char *fileVersion = NULL;
    filefmt_read_gen_header(rd, &fileType, &fileVersion);
    demand(strcmp(fileType, ppv_FILE_TYPE) == 0, "invalid file type");
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
    
    ppv_sample_t maxsmp; /* Max sample value. */
    if (strcmp(fileVersion, ppv_FILE_VERSION_2005) == 0)
      { /* Read bits-per-sample: */
        uint64_t bps = nget_uint64(rd, "bps", 10);
        demand(bps <= ppv_MAX_BPS, "too many bits per sample");
        maxsmp = ppv_max_sample((ppv_nbits_t)bps);
      }
    else if (strcmp(fileVersion, ppv_FILE_VERSION_2021) == 0)
      { /* Read max sample value: */
        uint64_t mu = nget_uint64(rd, "maxsmp", 10);
        demand(mu <= ppv_MAX_SAMPLE_VAL, "invalid max sample value");
        maxsmp = (ppv_sample_t)mu;
      }
    else
      { demand(FALSE, "invalid file version"); }
    fget_eol(rd); 

    /* Read plain/raw flag: */
    bool_t plain = nget_bool(rd, "plain");
    fget_eol(rd); 
    
    /* Allocate array using {asz}: */
    ppv_array_t *A = ppv_array_new(d, asz, maxsmp);
    
    if ((A->maxsmp > 0) && (! empty))
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
    ppv_sample_t maxsmp = A->maxsmp;
    if (ppv_index_first(ix, A)) 
      { ppv_pos_t pos = ppv_sample_pos(A, ix);
        do 
          { fget_skip_formatting_chars(rd);
            uint64_t usmp = fget_uint64(rd, 10);
            if (debug) 
              { fprintf(stderr, "    ix =");
                for (int32_t j = 0; j < d; j++) { fprintf(stderr, " " ppv_index_t_FMT, ix[j]); }
                fprintf(stderr, "  usmp = %lu\n", (uint64_t)usmp);
              }
            demand(usmp <= maxsmp, "sample too big");
            ppv_sample_t smp = (ppv_sample_t)usmp;
            ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, smp);
         } 
        while (! ppv_index_next(ix, A, d, &pos));
      }
    return;
  }

void ppv_array_read_samples_raw_small ( FILE *rd, ppv_array_t *A )
  {
    ppv_dim_t d = A->d;
    ppv_index_t ix[d];
    ppv_sample_t smask = (ppv_sample_t)((((uint64_t)1) << A->bps) - 1);
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
            ppv_sample_t smp = (buf >> shift) & smask;
            demand(smp <= A->maxsmp, "invalid sample value");
            ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, smp);
          }
        while (! ppv_index_next(ix, A, d, &pos));
      }
    return;
  }

void ppv_array_read_samples_raw_big ( FILE *rd, ppv_array_t *A )
  {
    ppv_dim_t d = A->d;
    ppv_index_t ix[d];
    ppv_sample_t smask = (ppv_sample_t)((((uint64_t)1) << A->bps) - 1); /* Mask to eat out a sample */
    int32_t cps = (A->bps + 7)/8; /* Bytes per sample. */
    assert(sizeof(ppv_word_t) >= cps);
    if (ppv_index_first(ix, A)) 
      { ppv_pos_t pos = ppv_sample_pos(A, ix);
        do
          { ppv_word_t w = 0;
            /* bool_t debug = ((A->bps == 10) && (pos == 0)); */
            for (int32_t i = 0; i < cps; i++)
              { int32_t ch = fgetc(rd);
                demand(ch != EOF, "unexpected end-of-file");
                /* if (debug) { fprintf(stderr,  "  ch = %u\n", (unsigned char)ch); } */
                w = (w << 8) | ((unsigned char)ch);
              }
            /* if (debug) { fprintf(stderr,  "w = " smpFMT "\n", w); } */
            ppv_sample_t smp = (ppv_sample_t)(w & smask);
            /* if (debug) { fprintf(stderr,  "smp = " smpFMT "\n", smp); } */
            demand(smp <= A->maxsmp, "invalid sample value");
            ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, smp);
          }
        while (! ppv_index_next(ix, A, d, &pos));
      }
    return;
  }

void ppv_array_read_samples_raw_bytes ( FILE *rd, ppv_array_t *A )
  {
    ppv_dim_t d = A->d;
    ppv_index_t ix[d];
    assert(A->bps == 8);
    if (ppv_index_first(ix, A)) 
      { ppv_pos_t pos = ppv_sample_pos(A, ix);
        do 
          { int32_t ch = fgetc(rd);
            demand(ch != EOF, "unexpected end-of-file");
            ppv_sample_t smp = (ppv_sample_t)ch;
            demand(smp <= A->maxsmp, "invalid sample value");
            ppv_set_sample_at_pos(A->el, A->bps, A->bpw, pos, smp);
         } 
        while (! ppv_index_next(ix, A, d, &pos));
      }
    return;
  }
