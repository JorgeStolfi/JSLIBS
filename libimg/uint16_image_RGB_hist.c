/* uint16_image_RGB_hist.c - RGB histogram tools
** Last edited on 2009-01-07 01:39:28 by stolfi
**
** Copied from Jef Poskanzer's libppm3.c - ppm utility library part 3
**
** Copyright (C) 1989, 1991 by Jef Poskanzer.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty.
*/

#include <jspnm.h>

#include <uint16_image_RGB_hist.h>
#include <uint16_image_RGB_table.h>

int ppm_equal(ppm_pixel_t *a, ppm_pixel_t *b)
  { return ((a->c[0] == b->c[0]) && (a->c[1] == b->c[1]) && (a->c[2] == b->c[2])); }

uint16_image_RGB_hist_vector uint16_image_RGB_hist_build
  ( uint16_t** samples, 
    int chns,
    int cols,
    int rows,
    int maxcolors, 
    int* colorsP
  )
  { uint16_image_RGB_table cht = uint16_image_RGB_table_build(samples, chns, cols, rows, maxcolors, colorsP);
    if (cht == NULL) { return NULL; }
    uint16_image_RGB_hist_vector chv = uint16_image_RGB_table_to_hist(cht, maxcolors);
    uint16_image_RGB_table_free(cht);
    return chv;
  }

void uint16_image_RGB_hist_add
  ( uint16_image_RGB_hist_vector chv, 
    int* colorsP, 
    int maxcolors, 
    ppm_pixel_t* colorP,
    int value,
    int position
  )
  {
    int i, j;
    /* Search colorhist for the color. */
    for (i = 0; i < *colorsP; ++i)
      { if (ppm_equal(&(chv[i].color), colorP))
        { /* Found it - move to new slot. */
          if (position > i)
            { for (j = i; j < position; ++j) { chv[j] = chv[j + 1]; } }
          else if (position < i)
            { for (j = i; j > position; --j) { chv[j] = chv[j - 1]; } }
          chv[position].color = *colorP;
          chv[position].value = value;
          return;
        }
      }
    if (*colorsP < maxcolors)
      {
        /* Didn't find it, but there's room to add it; so do so. */
        for (i = *colorsP; i > position; --i) { chv[i] = chv[i - 1];}
        chv[position].color = *colorP;
        chv[position].value = value;
        ++(*colorsP);
      }
    }

void uint16_image_RGB_hist_free(uint16_image_RGB_hist_vector chv)
  {
    free((char*) chv);
  }

