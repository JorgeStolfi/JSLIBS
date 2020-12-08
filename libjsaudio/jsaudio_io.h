/* jsaudio_io.h - basic I/O tools for audio files */
/* Last edited on 2006-10-29 20:24:16 by stolfi */

/* Created by Jorge Stolfi on sep/2006. */

#ifndef jsaudio_io_H
#define jsaudio_io_H

#include <stdio.h>
#include <stdint.h>

/* The multibyte procedures with {_be} ending assume big-endian format. */

uint8_t jsa_read_uint8(FILE *rd);
  /* Reads an 8 bit unsigned integer from {rd}. */

uint16_t jsa_read_uint16_be(FILE *rd);
  /* Reads a 16 bit unsigned integer from {rd}. */

short jsa_read_short_be(FILE *rd);
  /* Reads a 16-bit signed integer from {rd}. */

uint32_t jsa_read_uint32_be(FILE *rd);
  /* Reads a 32-bit unsigned integer from {rd}. */

int jsa_read_int_be(FILE *rd);
  /* Reads a 32-bit signed integer from {rd}. */

uint64_t jsa_read_uint64_be(FILE *rd);
  /* Reads a 64-bit unsigned integer from {rd}. */

float jsa_read_float_be(FILE *rd);
  /* Reads a 32-bit {float} from {rd}. */

void jsa_read_chars(FILE *rd, int n, char s[]);
  /* Reads the next {n} characters in {rd}, stores them in {s[0..n-1]}. */

/* WRITING DATA */

void jsa_write_uint32_be(FILE *wr, uint32_t *p);
  /* Writes a 32 bit unsigned integer {*p} to file {wr}. */

void jsa_write_float_be(FILE *wr, float *p);
  /* Writes a 32-bit float {*p} to file {wr}. */

#endif
