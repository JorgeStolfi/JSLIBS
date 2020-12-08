/* uint16_image_write.h - reading {uint16_image_t} from a PNM (PBM/PGM/PPM) image file. */
/* Last edited on 2017-06-30 00:55:59 by stolfilocal */

#ifndef uint16_image_write_pnm_H
#define uint16_image_write_pnm_H

#define _GNU_SOURCE
#include <bool.h>
#include <stdio.h>

#include <uint16_image.h>

void uint16_image_write_pnm_named(char *name, uint16_image_t *img, bool_t forceplain, bool_t verbose);
  /* Writes a PGM or PPM image to the specified file, in a format
    compatibe with {img.chns}: "P2" or "P5" if {img.chns == 1}, "P3" or
    "P6" if {img.chns == 3}. Fails if {img.chns} is anything else.
    The field {img.maxval} must have been set properly.
    
    The plain variant of the format ("P2" or "P3") is chosen if
    {forceplain} is true or {img.maxval} is too big for the raw
    variant ("P5" or "P6"). If {name} is "-", writes to {stdout}.
    If {verbose} is TRUE, prints a notice to {stderr}. */

void uint16_image_write_pnm_file(FILE *wr, uint16_image_t *img, bool_t forceplain, bool_t verbose);
  /* Same as {uint16_image_write_pnm_named}, 
    but to a previously opened file handle {wr}. */

#endif
