/* See jsaudio_au.h */
/* Last edited on 2024-12-05 10:32:11 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <jsaudio.h>
#include <jsaudio_au.h>
#include <jsaudio_io.h>

/* IMPLEMENTATIONS */

int32_t jsa_au_file_bytes_per_sample(int32_t enc)
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

sound_t jsa_read_au_file(FILE *rd)
  {
    /* Reads the file header: */
    au_file_header_t h = jsa_read_au_file_header(rd);
    
    /* Compute number of channels {nc} and samples per channel {ns}: */
    int32_t nc = h.channels;
    int32_t bps = jsa_au_file_bytes_per_sample(h.encoding);
    int32_t ns = h.data_size / nc / bps;
    assert(h.data_size == ns * nc * bps); 
    
    /* Allocates the {sound_t} structure and fills header data: */
    sound_t s = jsa_allocate_sound(nc, ns);
    s.fsmp = (double)h.sample_rate;
    s.ns = ns;
    s.nc = nc;
    
    /* Reads the samples: */
    jsa_read_au_file_samples(rd, &h, &s, 0, ns);

    return s;
  }
     
au_file_header_t jsa_read_au_file_header(FILE *rd)
  {
    au_file_header_t h;
    h.magic = jsa_read_uint32_be(rd);        /* Magic ID number of AU file. */  
    h.hdr_size = jsa_read_uint32_be(rd);     /* Byte offset to first sample. */
    h.data_size = jsa_read_uint32_be(rd);    /* Total bytes in sample data. */
    h.encoding = jsa_read_uint32_be(rd);     /* Encoding format. */
    h.sample_rate = jsa_read_uint32_be(rd);  /* Sampling frequency. */
    h.channels = jsa_read_uint32_be(rd);     /* Number of channels. */
    
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

void jsa_skip_au_file_samples(FILE *rd, au_file_header_t *h, int32_t ns)
  {
    int32_t bps = jsa_au_file_bytes_per_sample(h->encoding); /* Bytes per sample. */
    int32_t bskip = ns * bps; /* Bytes to skip. */
    for (uint32_t i = 0;  i < bskip; i++)
      { int32_t chr = fgetc(rd); assert(chr != EOF); }
  }

void jsa_read_au_file_samples(FILE *rd, au_file_header_t *h, sound_t *s, int32_t skip, int32_t ns)
  {
    /* Get number of channels: */
    int32_t nc = h->channels;
    assert(nc == s->nc);
    
    /* Check if elements do exist in {s}: */
    assert(skip >= 0);
    assert(skip + ns <= s->ns);

    /* Loop on samples and channels: */
    int32_t enc = h->encoding; /* A shorter name for the encoding tag. */
    for (uint32_t i = 0;  i < ns; i++) 
      { for (uint32_t ic = 0;  ic < nc; ic++)
          { double dv;
            switch(enc)
              {
                case 3:
                  { /* Reads 16 bit integer, converts to {double}. */
                    int16_t rv = jsa_read_int16_be(rd);
                    dv = ((double)rv)/SCALE_INT16_T;
                  }
                  break;
                case 5:
                  { /* Reads 32 bit integer, converts to {double}. */
                    int32_t rv = jsa_read_int32_be(rd);
                    dv = ((double)rv)/SCALE_INT32_T;
                  }
                  break;
                case 6:
                  { /* Reads 32 bit float, converts to {double}. */
                    float rv = jsa_read_float_be(rd);
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
  
void jsa_write_au_file(FILE *wr, sound_t *s)
  {
    /* Chooses the encoding: */
    int32_t enc = 6; /* For now, use 32-bit float encoding. */
    int32_t bps = jsa_au_file_bytes_per_sample(enc); /* Bytes per sample. */

    /* Allocate a file header: */
    au_file_header_t h;

    /* Build the file header: */
    h.magic = 0x2e736e64; /* ".snd" */
    h.hdr_size = 24;    
    h.data_size = s->nc * s->ns * bps; 
    h.encoding = enc;   
    h.sample_rate = (int32_t)(s->fsmp + 0.5);
    if (fabs(s->fsmp - (double)h.sample_rate)/s->fsmp > 1.0e-4)
      { fprintf(stderr, "warning: sample rate %14.8e rounded to %d\n", s->fsmp, h.sample_rate); }
    h.channels = s->nc;
    
    /* Write the file header: */
    jsa_write_au_file_header(wr, &h);
    
    /* Writes the sample data: */
    jsa_write_au_file_samples(wr, &h, s, 0, s->ns);
  }

void jsa_write_au_file_header(FILE *wr, au_file_header_t *h)
  { jsa_write_uint32_be(wr, &(h->magic));         /* Magic ID number of AU file. */  
    jsa_write_uint32_be(wr, &(h->hdr_size));      /* Byte offset to first sample. */ 
    jsa_write_uint32_be(wr, &(h->data_size));     /* Total bytes in sample data. */  
    jsa_write_uint32_be(wr, &(h->encoding));      /* Encoding format. */             
    jsa_write_uint32_be(wr, &(h->sample_rate));   /* Sampling frequency. */          
    jsa_write_uint32_be(wr, &(h->channels));      /* Number of channels. */ 
    fflush(wr);
  }

void jsa_write_au_file_samples(FILE *wr, au_file_header_t *h, sound_t *s, int32_t skip, int32_t ns)
  {
    assert(skip >= 0);
    assert(skip + ns <= s->ns);
    assert(h->encoding == 6); /* For now. */
    
    for (uint32_t i = 0;  i < ns; i++) 
      { for (uint32_t ic = 0;  ic < s->nc; ic++)
          { float fv = (float)(s->sv[ic][skip + i]);
            jsa_write_float_be(wr, &fv);
          }
      }
    fflush(wr);
  }

