/* uint16_image_read.h - reading {uint16_image_t} from a PNM (PBM/PGM/PPM) image file. */
/* Last edited on 2017-06-22 02:29:37 by stolfilocal */

#ifndef uint16_image_read_pnm_H
#define uint16_image_read_pnm_H

#define _GNU_SOURCE
#include <stdio.h>

#include <bool.h>
#include <uint16_image.h>

uint16_image_t *uint16_image_read_pnm_named(char *name, bool_t verbose);
  /* Reads a PBM, PGM pr PPM image from the named file. The file format is
    determined by the file's `magic number': "P1" or "P4"
    for the PBM format (1 channel,0 or 1), "P2" or "P5" for PGM (1
    channel --- luminance), "P3" or "P6" for PPM (3 channels --- Red,
    Green, and Blue). If the {name} is "-", reads from {stdin}.
    The field {img.maxval} is set as specified in the file.
    If {verbose} is TRUE, prints a notice to {stderr}. */

uint16_image_t *uint16_image_read_pnm_file(FILE *rd);
  /* Same as {uint16_image_read_pnm_named}, 
    but from a previously opened file handle {rd}. */

#endif
