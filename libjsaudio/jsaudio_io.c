/* See jsaudio_io.h */
/* Last edited on 2017-01-02 21:31:00 by jstolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <jsaudio_io.h>

/* IMPLEMENTATIONS */

uint8_t jsa_read_uint8(FILE *rd)
  { uint8_t u;
    int c;
    c = fgetc(rd); assert(c != EOF); u = (uint8_t)(c & 0xff);
    return u;
  }

uint16_t jsa_read_uint16_be(FILE *rd)
  { uint16_t u;
    int c;
    c = fgetc(rd); assert(c != EOF); u = (uint16_t)(c & 0xff);
    c = fgetc(rd); assert(c != EOF); u = (uint16_t)((u << 8) | (c & 0xff));
    return u;
  }

short jsa_read_short_be(FILE *rd)
  { uint16_t u = jsa_read_uint16_be(rd);
    short r;
    (*(uint16_t*)(&r)) = u;
    return r;
  }

uint32_t jsa_read_uint32_be(FILE *rd)
  { uint32_t u;
    int c;
    c = fgetc(rd); assert(c != EOF); u = (c & 0xff);
    c = fgetc(rd); assert(c != EOF); u = (u << 8) | (c & 0xff);
    c = fgetc(rd); assert(c != EOF); u = (u << 8) | (c & 0xff);
    c = fgetc(rd); assert(c != EOF); u = (u << 8) | (c & 0xff);
    return u;
  }

uint64_t jsa_read_uint64_be(FILE *rd)
  { uint64_t u;
    int c;
    c = fgetc(rd); assert(c != EOF); u = (c & 0xff);
    c = fgetc(rd); assert(c != EOF); u = (u << 8) | (c & 0xff);
    c = fgetc(rd); assert(c != EOF); u = (u << 8) | (c & 0xff);
    c = fgetc(rd); assert(c != EOF); u = (u << 8) | (c & 0xff);
    c = fgetc(rd); assert(c != EOF); u = (u << 8) | (c & 0xff);
    c = fgetc(rd); assert(c != EOF); u = (u << 8) | (c & 0xff);
    c = fgetc(rd); assert(c != EOF); u = (u << 8) | (c & 0xff);
    c = fgetc(rd); assert(c != EOF); u = (u << 8) | (c & 0xff);
    return u;
  }

int jsa_read_int_be(FILE *rd)
  { uint32_t u = jsa_read_uint32_be(rd);
    int r;
    (*(uint32_t*)(&r)) = u;
    return r;
  }

float jsa_read_float_be(FILE *rd)
  { uint32_t u = jsa_read_uint32_be(rd);
    float r;
    (*(uint32_t*)(&r)) = u;
    return r;
  }

void jsa_write_uint32_be(FILE *wr, uint32_t *p)
  { uint32_t u = (*p);
    fputc((u & 0xff000000) >> 24,wr);
    fputc((u & 0x00ff0000) >> 16,wr);
    fputc((u & 0x0000ff00) >>  8,wr);
    fputc((u & 0x000000ff) >>  0,wr);
  }

void jsa_write_float_be(FILE *wr, float *p)
  { /* Cast the 32-bit float as an {uint32_t}: */
    uint32_t u = (*(uint32_t*)p);
    jsa_write_uint32_be(wr, &u); 
  }
