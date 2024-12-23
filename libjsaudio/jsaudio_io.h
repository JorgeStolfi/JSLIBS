/* jsaudio_io.h - basic I/O tools for audio files */
/* Last edited on 2024-12-21 03:19:42 by stolfi */

/* Created by Jorge Stolfi on sep/2006. */

#ifndef jsaudio_io_H
#define jsaudio_io_H

#include <stdio.h>
#include <stdint.h>

/* The multibyte procedures with {_be} ending assume big-endian format. */

uint8_t jsaudio_read_uint8(FILE *rd);
  /* Reads an 8 bit unsigned integer from {rd}. */

uint16_t jsaudio_read_uint16_be(FILE *rd);
  /* Reads a 16 bit unsigned integer from {rd}. */

short jsaudio_read_int16_be(FILE *rd);
  /* Reads a 16-bit signed integer from {rd}. */

uint32_t jsaudio_read_uint32_be(FILE *rd);
  /* Reads a 32-bit unsigned integer from {rd}. */

int32_t jsaudio_read_int32_be(FILE *rd);
  /* Reads a 32-bit signed integer from {rd}. */

uint64_t jsaudio_read_uint64_be(FILE *rd);
  /* Reads a 64-bit unsigned integer from {rd}. */

float jsaudio_read_float_be(FILE *rd);
  /* Reads a 32-bit {float} from {rd}. */

void jsaudio_read_chars(FILE *rd, int32_t n, char s[]);
  /* Reads the next {n} characters in {rd}, stores them in {s[0..n-1]}. */

/* WRITING DATA */

void jsaudio_write_uint32_be(FILE *wr, uint32_t *p);
  /* Writes a 32 bit unsigned integer {*p} to file {wr}. */

void jsaudio_write_float_be(FILE *wr, float *p);
  /* Writes a 32-bit float {*p} to file {wr}. */

#endif
