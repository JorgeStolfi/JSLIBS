/* See {jspng.h} */
/* Last edited on 2017-06-23 00:24:16 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#include <png.h>

#include <bool.h>
#include <affirm.h>

#include <jspng.h>

void jspng_dump_info(FILE *wr, const char *func, char *label, png_structp pr, png_infop pi)
  { 
    fprintf(wr, "------------------------------------------------------------\n");
    fprintf(wr, "PNG file info - %s - %s\n", func, label);
    fprintf(wr, "image_width = %u\n", (uint32_t)png_get_image_width(pr, pi));
    fprintf(wr, "image_height = %u\n", (uint32_t)png_get_image_height(pr, pi));
    fprintf(wr, "channels = %u\n", (uint32_t)png_get_channels(pr, pi));
    fprintf(wr, "bit_depth = %d\n", png_get_bit_depth (pr, pi));
    /* fprintf(wr, "usr_bit_depth = %d\n", (int32_t)pr->usr_bit_depth; */
    fprintf(wr, "color_type = %d\n", png_get_color_type(pr, pi));
    fprintf(wr, "interlace_type = %d\n", png_get_interlace_type(pr, pi));
    fprintf(wr, "compression_type = %d\n", png_get_compression_type(pr, pi));
    fprintf(wr, "filter_type = %d\n", png_get_filter_type(pr, pi));
    
    png_color_8p sBIT; 
    if (png_get_sBIT(pr, pi, &sBIT)) 
      { fprintf(wr, "sBIT = %16p = ( %u %u %u %u %u )\n", (void *)sBIT, sBIT->gray, sBIT->red, sBIT->green, sBIT->blue, sBIT->alpha); }
    else
      { fprintf(wr, "sBIT not specified\n"); }
    
    int has_tRNS = png_get_valid(pr, pi, PNG_INFO_tRNS);
    fprintf(wr, "has tRNS = %d\n", has_tRNS);
    
    double gamma;
    if (png_get_gAMA(pr, pi, &gamma))
      { fprintf(wr, "gAMA = %25.16e\n", gamma); }
    else
      { fprintf(wr, "gAMA not specified\n"); }
    
    fprintf(wr, "channels = %d\n", png_get_channels(pr, pi));
    fprintf(wr, "rowbytes = %u\n", (uint32_t)png_get_rowbytes(pr, pi));
    png_bytep sg = (png_bytep)png_get_signature(pr, pi);
    if (sg != NULL)
      { fprintf(wr, "signature = %16p = %02x %02x %02x %02x %02x %02x %02x %02x\n", (void *)sg, sg[0], sg[1], sg[2], sg[3], sg[4], sg[5], sg[6], sg[7]); }
    else
      { fprintf(wr, "png_get_signature failed\n"); }
    /* fprintf(wr, " = %d\n", png_get_(pr, pi)); */
    fprintf(wr, "------------------------------------------------------------\n");

  }
