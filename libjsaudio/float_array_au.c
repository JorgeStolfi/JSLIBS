/* See float_array_au.h */
/* Last edited on 2024-12-05 10:32:00 by stolfi */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <ix.h>
#include <affirm.h>
#include <bool.h>
#include <jsaudio.h>
#include <jsaudio_au.h>
#include <jsaudio_io.h>

#include <float_array.h>
#include <float_array_au.h>

#define NA  float_array_NAXES
  /* Number of dimensions. */

/* IMPLEMENTATIONS */
 
float_array_t *float_array_from_sound(sound_t *snd, bool_t verbose)
  { 
    /* Get sound channel and sample counts: */
    int32_t NC = snd->nc; /* Num of channels. */
    int32_t NS = snd->ns; /* Num of samples per channel. */
    
    /* Allocate float image: */
    ix_dim_t na = 2;
    ix_size_t sz[2];
    sz[0] = NC; sz[1] = NS;
    float_array_t AN = float_array_new(na, sz);
    float_array_t *A = &AN;
    
    /* Input and output range registers: */
    float vmin[NC], vmax[NC];  /* Output range registers. */ 
    int32_t c; /* Channel index. */
    for (c = 0; c < NC; c++) { vmin[c] = +INFINITY; vmax[c] = -INFINITY; }
    
    /* Convert pixels, keep statistics: */
    ix_index_t ix[na];
    int32_t s;
    for(s = 0; s < NS; s++)
      { double *srow = snd->sv[s];
        for (c = 0; c < NC; c++)
          { /* Convert int32_t sample {*srow} to float {v}, store, keep stats: */
            double ismp = (*srow);
            float osmp = (float)ismp;
            if (osmp < vmin[c]) { vmin[c] = osmp; }
            if (osmp > vmax[c]) { vmax[c] = osmp; }
            ix[0] = c; ix[1] = s; ix[2] = ix[3] = ix[4] = ix[5] = 0;
            float_array_set_elem(A, ix, osmp);
            srow++;
          }
      }
    
    if (verbose) 
      { /* Print statistics: */
        int64_t NCS = ((int64_t)NC)*((int64_t)NS);
        fprintf(stderr, "  %d channels, %d samples per channel, %ld tot samples\n", NC, NS, NCS);
        if (NCS > 0)
          { for (c = 0; c < NC; c++)
              { int32_t iChan = c;  /* Channel index in input sound array. */
                int32_t oChan = c;  /* Channel index in output float array. */
                fprintf(stderr, "  converted au channel %d to array channel %d:\n", iChan, oChan);
                fprintf(stderr, "    actual output range = [ %14.7e _ %14.7e]\n", vmin[c], vmax[c]);
              }
          }
      }
    return A;
  }

sound_t *float_array_to_sound(float_array_t *A, int32_t chns, int32_t ch[], double freq, bool_t verbose)
  { 
    /* Get indexing descriptor {DA}: */
    ix_descr_t *DA = &(A->ds);
    
    /* Get integer image dimensions: */
    int32_t na = DA->na;
    assert(na == 2);
    int32_t NS = (int32_t)DA->sz[1];  /* Num of samples per channel. */
    int32_t sNC = chns;           /* Num channels in sound array. */
    int32_t fNC = (int32_t)DA->sz[0]; /* Num channels in float image. */
    
    /* Allocate sound array: */
    sound_t *snd = (sound_t *)notnull(malloc(sizeof(sound_t)), "no mem");
    (*snd) = jsa_allocate_sound(sNC, NS);
    
    /* Set max sample value in integer image: */
    snd->fsmp = freq;
    
    /* Channel indexing variables: */
    int32_t k; /* Channel of sound array. */
    int32_t c; /* Channel of float array. */
    
    /* Input and output range registers: */
    double vmin[sNC], vmax[sNC];   /* Sound pixel range. */
    for (k = 0; k < sNC; k++) { vmin[k] = +INFINITY; vmax[k] = -INFINITY; }
    
    /* Convert pixels, store in {snd}, keep statistics: */
    ix_index_t ix[na];
    int32_t s;
    for(s = 0; s < NS; s++)
      { double *srow = snd->sv[s];
        /* Convert float pixel {fpxy[c..c+2]} to integer pixel {ipxy[0..2]}, keep stats: */
        for (k = 0; k < sNC; k++)
          { c = (ch == NULL ? k : ch[k]);
            ix[0] = c; ix[1] = s; ix[2] = ix[3] = ix[4] = ix[5] = 0;
            float ismp = ((c < 0) || (c >= fNC) ? 0.0f : float_array_get_elem(A, ix));
            double osmp = (double)ismp;
            if (osmp < vmin[k]) { vmin[k] = osmp; }
            if (osmp > vmax[k]) { vmax[k] = osmp; }
            (*srow) = osmp;
            srow++;
          }
      }
    
    if (verbose)
      { /* Print statistics: */
        int64_t NCS = ((int64_t)sNC)*((int64_t)NS);
        fprintf(stderr, "  %d channels, %d samples per channel, %ld tot samples\n", sNC, NS, NCS);
        if (NCS > 0)
          { for (k = 0; k < chns; k++)
              { c = (ch == NULL ? k : ch[k]);
                int32_t iChan = c;  /* Channel index in input float array. */
                int32_t oChan = k;  /* Channel index in output sound array. */
                fprintf(stderr, "  converted array channel %d to sound channel %d:\n", iChan, oChan);
                fprintf(stderr, "    actual output range = [ %24.16e .. %24.16e ]\n", vmin[k], vmax[k]);
              }
          }
      }
    
    return snd;
  }

