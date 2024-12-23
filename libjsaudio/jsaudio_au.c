/* See jsaudio_au.h */
/* Last edited on 2024-12-21 03:45:16 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <jsaudio.h>
#include <jsaudio_au.h>
#include <jsaudio_io.h>

/* IMPLEMENTATIONS */

uint32_t jsaudio_au_file_bytes_per_sample(uint32_t enc)
  {
    switch(enc)
      {
        case  1: /*  8 bit ISDN u-law*/
        case  2: /*  8 bit linear PCM*/
        case 27: /*  8 bit ISDN a-law*/
          return 1;
          break;

        case  3: /*  16 bit linear PCM*/
          return 2;
          break;

        case  4: /*  24 bit linear PCM*/
          return 3;
          break;

        case  5: /*  32 bit linear PCM*/
        case  6: /*  32 bit IEEE floating point*/
          return 4;
          break;

        case  7: /*  64 bit IEEE floating point*/
          return 8;
          break;

        case 23: /*  4 bit CCITT G721 ADPCM*/
        case 24: /*  CCITT G722 ADPCM*/
        case 25: /*  CCITT G723 ADPCM*/
        case 26: /*  5 bit CCITT G723 ADPCM*/
          fprintf(stderr, "reading of encoding type %d is not supported\n", enc);
          exit(1);
          return 0;  /* To prevent the compiler from complaining. */
          break;

        default:
          fprintf(stderr, "encoding type %d is undefined\n", enc);
          exit(1);
          return 0;  /* To prevent the compiler from complaining. */
      }
  }

jsaudio_t jsaudio_au_read_file(FILE *rd)
  {
    /* Reads the file header: */
    jsaudio_au_file_header_t h = jsaudio_au_read_file_header(rd);
    
    /* Compute number of channels {nc} and samples per channel {ns}: */
    uint32_t nc = h.channels;
    uint32_t bps = jsaudio_au_file_bytes_per_sample(h.encoding);
    uint32_t ns = h.data_size / nc / bps;
    assert(h.data_size == ns * nc * bps); 
    
    /* Allocates the {jsaudio_t} structure and fills header data: */
    jsaudio_t s = jsaudio_allocate_sound(nc, ns);
    s.fsmp = (double)h.sample_rate;
    s.ns = ns;
    s.nc = nc;
    
    /* Reads the samples: */
    jsaudio_au_read_file_samples(rd, &h, &s, 0, ns);

    return s;
  }
     
jsaudio_au_file_header_t jsaudio_au_read_file_header(FILE *rd)
  {
    jsaudio_au_file_header_t h;
    h.magic = jsaudio_read_uint32_be(rd);        /* Magic ID number of AU file. */  
    h.hdr_size = jsaudio_read_uint32_be(rd);     /* Byte offset to first sample. */
    h.data_size = jsaudio_read_uint32_be(rd);    /* Total bytes in sample data. */
    h.encoding = jsaudio_read_uint32_be(rd);     /* Encoding format. */
    h.sample_rate = jsaudio_read_uint32_be(rd);  /* Sampling frequency. */
    h.channels = jsaudio_read_uint32_be(rd);     /* Number of channels. */
    
    /* Header consistency checks: */
    assert(h.magic == 0x2e736e64); /* ".snd" */
    assert
      ( ((h.encoding >= 1) && (h.encoding <= 7)) || 
        ((h.encoding >= 23) && (h.encoding <= 27))
      ); 
    assert(h.sample_rate > 0);
    assert(h.data_size > 0);           /* Requires explicit {data_size}. */
    assert(h.data_size != 0xffffffff); /* Requires explicit {data_size}. */
    assert((h.channels >= 1) && (h.channels <= 16));
    
    /* Consume stream until first data byte: */
    int32_t nread = 24;
    while(nread < h.hdr_size) { int32_t chr = fgetc(rd); assert(chr != EOF); nread++; }
      
    return h;
  }

#define SCALE_INT32_T (2.0 * 1024.0 * 1024.0 * 1024.0) /* {2^31} */
#define SCALE_INT16_T (32.0 * 1024.0) /* {2^15} */

void jsaudio_au_skip_file_samples(FILE *rd, jsaudio_au_file_header_t *h, uint32_t ns)
  {
    uint32_t bps = jsaudio_au_file_bytes_per_sample(h->encoding); /* Bytes per sample. */
    uint32_t bskip = ns * h->channels * bps; /* Bytes to skip. */
    for (int32_t i = 0; i < bskip; i++)
      { int32_t chr = fgetc(rd); assert(chr != EOF); }
  }

void jsaudio_au_read_file_samples(FILE *rd, jsaudio_au_file_header_t *h, jsaudio_t *s, uint32_t skip, uint32_t ns)
  {
    /* Get number of channels: */
    uint32_t nc = h->channels;
    assert(nc == s->nc);
    
    /* Check if elements do exist in {s}: */
    assert(skip + ns <= s->ns);

    /* Loop on samples and channels: */
    uint32_t enc = h->encoding; /* A shorter name for the encoding tag. */
    for (uint32_t i = 0; i < ns; i++) 
      { for (uint32_t ic = 0; ic < nc; ic++)
          { double dv;
            switch(enc)
              { case 3:
                  { /* Reads 16 bit integer, converts to {double}. */
                    int16_t rv = jsaudio_read_int16_be(rd);
                    dv = ((double)rv)/SCALE_INT16_T;
                  }
                  break;
                case 5:
                  { /* Reads 32 bit integer, converts to {double}. */
                    int32_t rv = jsaudio_read_int32_be(rd);
                    dv = ((double)rv)/SCALE_INT32_T;
                  }
                  break;
                case 6:
                  { /* Reads 32 bit float, converts to {double}. */
                    float rv = jsaudio_read_float_be(rd);
                    dv = (double)rv;
                  }
                  break;
                default:
                  fprintf(stderr, "decoding of encoding type %d is not supported\n", enc);
                  exit(1);
              }
            s->sv[ic][i] = dv;
          }
      }
  }
  
void jsaudio_au_write_file(FILE *wr, jsaudio_t *s)
  {
    /* Chooses the encoding: */
    uint32_t enc = 6; /* For now, use 32-bit float encoding. */
    uint32_t bps = jsaudio_au_file_bytes_per_sample(enc); /* Bytes per sample. */

    /* Allocate a file header: */
    jsaudio_au_file_header_t h;

    /* Build the file header: */
    h.magic = 0x2e736e64; /* ".snd" */
    h.hdr_size = 24;    
    h.data_size = s->nc * s->ns * bps; 
    h.encoding = enc;   
    h.sample_rate = (uint32_t)(s->fsmp + 0.5);
    if (fabs(s->fsmp - (double)h.sample_rate)/s->fsmp > 1.0e-4)
      { fprintf(stderr, "warning: sample rate %14.8e rounded to %d\n", s->fsmp, h.sample_rate); }
    h.channels = s->nc;
    
    /* Write the file header: */
    jsaudio_au_write_file_header(wr, &h);
    
    /* Writes the sample data: */
    jsaudio_au_write_file_samples(wr, &h, s, 0, s->ns);
  }

void jsaudio_au_write_file_header(FILE *wr, jsaudio_au_file_header_t *h)
  { jsaudio_write_uint32_be(wr, &(h->magic));         /* Magic ID number of AU file. */  
    jsaudio_write_uint32_be(wr, &(h->hdr_size));      /* Byte offset to first sample. */ 
    jsaudio_write_uint32_be(wr, &(h->data_size));     /* Total bytes in sample data. */  
    jsaudio_write_uint32_be(wr, &(h->encoding));      /* Encoding format. */             
    jsaudio_write_uint32_be(wr, &(h->sample_rate));   /* Sampling frequency. */          
    jsaudio_write_uint32_be(wr, &(h->channels));      /* Number of channels. */ 
    fflush(wr);
  }

void jsaudio_au_write_file_samples(FILE *wr, jsaudio_au_file_header_t *h, jsaudio_t *s, uint32_t skip, uint32_t ns)
  {
    assert(skip >= 0);
    assert(skip + ns <= s->ns);
    assert(h->encoding == 6); /* For now. */
    
    for (uint32_t i = 0;  i < ns; i++) 
      { for (uint32_t ic = 0;  ic < s->nc; ic++)
          { float fv = (float)(s->sv[ic][skip + i]);
            jsaudio_write_float_be(wr, &fv);
          }
      }
    fflush(wr);
  }

