/* See jsaudio_io.h */
/* Last edited on 2024-12-21 03:19:55 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <jsaudio_io.h>

/* IMPLEMENTATIONS */

uint8_t jsaudio_read_uint8(FILE *rd)
  { uint8_t u;
    int32_t chr;
    chr = fgetc(rd); assert(chr != EOF); u = (uint8_t)(chr & 0xff);
    return u;
  }

uint16_t jsaudio_read_uint16_be(FILE *rd)
  { uint16_t u;
    int32_t chr;
    chr = fgetc(rd); assert(chr != EOF); u = (uint16_t)(chr & 0xff);
    chr = fgetc(rd); assert(chr != EOF); u = (uint16_t)((u << 8) | (chr & 0xff));
    return u;
  }

short jsaudio_read_int16_be(FILE *rd)
  { uint16_t u = jsaudio_read_uint16_be(rd);
    short r;
    (*(uint16_t*)(&r)) = u;
    return r;
  }

uint32_t jsaudio_read_uint32_be(FILE *rd)
  { uint32_t u;
    int32_t chr;
    chr = fgetc(rd); assert(chr != EOF); u = (chr & 0xff);
    chr = fgetc(rd); assert(chr != EOF); u = (u << 8) | (chr & 0xff);
    chr = fgetc(rd); assert(chr != EOF); u = (u << 8) | (chr & 0xff);
    chr = fgetc(rd); assert(chr != EOF); u = (u << 8) | (chr & 0xff);
    return u;
  }

uint64_t jsaudio_read_uint64_be(FILE *rd)
  { uint64_t u;
    int32_t chr;
    chr = fgetc(rd); assert(chr != EOF); u = (chr & 0xff);
    chr = fgetc(rd); assert(chr != EOF); u = (u << 8) | (chr & 0xff);
    chr = fgetc(rd); assert(chr != EOF); u = (u << 8) | (chr & 0xff);
    chr = fgetc(rd); assert(chr != EOF); u = (u << 8) | (chr & 0xff);
    chr = fgetc(rd); assert(chr != EOF); u = (u << 8) | (chr & 0xff);
    chr = fgetc(rd); assert(chr != EOF); u = (u << 8) | (chr & 0xff);
    chr = fgetc(rd); assert(chr != EOF); u = (u << 8) | (chr & 0xff);
    chr = fgetc(rd); assert(chr != EOF); u = (u << 8) | (chr & 0xff);
    return u;
  }

int32_t jsaudio_read_int32_be(FILE *rd)
  { uint32_t u = jsaudio_read_uint32_be(rd);
    int32_t r;
    (*(uint32_t*)(&r)) = u;
    return r;
  }

float jsaudio_read_float_be(FILE *rd)
  { uint32_t u = jsaudio_read_uint32_be(rd);
    float r;
    (*(uint32_t*)(&r)) = u;
    return r;
  }

void jsaudio_write_uint32_be(FILE *wr, uint32_t *p)
  { uint32_t u = (*p);
    fputc((int32_t)((u & 0xff000000) >> 24), wr);
    fputc((int32_t)((u & 0x00ff0000) >> 16), wr);
    fputc((int32_t)((u & 0x0000ff00) >>  8), wr);
    fputc((int32_t)((u & 0x000000ff) >>  0), wr);
  }

void jsaudio_write_float_be(FILE *wr, float *p)
  { /* Cast the 32-bit float as an {uint32_t}: */
    uint32_t u = (*(uint32_t*)p);
    jsaudio_write_uint32_be(wr, &u); 
  }
