/* See uint16_image_RGB_table.h
** Last edited on 2009-01-07 03:13:15 by stolfi
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

#include <affirm.h>
#include <stdint.h>

#include <jspnm.h>

#include <uint16_image_RGB_hist.h>
#include <uint16_image_RGB_table.h>

#define HASH_SIZE 20023

uint32_t uint16_image_RGB_table_hash_pixel(ppm_pixel_t *p);
  /* Computes a hash index for an RGB pixel. */

uint32_t uint16_image_RGB_table_hash_pixel(ppm_pixel_t *p)
  { uint32_t r = p->c[0];
    uint32_t g = p->c[1];
    uint32_t b = p->c[2];
    return ((r * 33023 + g * 30013 + b * 27011) & 0x7fffffff) % HASH_SIZE;
  }

uint16_image_RGB_table uint16_image_RGB_table_build
  ( uint16_t **samples,
    int chns, 
    int cols, 
    int rows, 
    int maxcolors, 
    int *colorsP )
  {
    uint16_image_RGB_table cht;
    uint16_image_RGB_bucket chl;
    int col, row, hash;

    cht = uint16_image_RGB_table_alloc( );
    *colorsP = 0;

    /* Go through the entire image, building a hash table of colors. */
    for (row = 0; row < rows; ++row)
      { uint16_t *sP = samples[row];
        for (col = 0; col < cols; ++col)
          { /* Grab a pixel {p}, expanding/truncating to 3 samples: */ 
            ppm_pixel_t p;
            if (chns > 0) { p.c[0] = *sP; sP++; } else {p.c[0] = 0; }
            if (chns > 1) { p.c[1] = *sP; sP++; } else {p.c[1] = p.c[0]; }
            if (chns > 2) { p.c[2] = *sP; sP++; } else {p.c[2] = p.c[0]; }
            if (chns > 3) { sP += (chns-3); }
            hash = uint16_image_RGB_table_hash_pixel(&p);
            for (chl = cht[hash]; chl != (uint16_image_RGB_bucket)NULL; chl = chl->next)
              { if (ppm_equal(&(chl->ch.color), &p)) break; }
            if (chl != (uint16_image_RGB_bucket)NULL) 
              { ++(chl->ch.value); }
            else
              { if (++(*colorsP) > maxcolors)
                  { uint16_image_RGB_table_free(cht);
                    return (uint16_image_RGB_table)NULL;
                  }
                chl = (uint16_image_RGB_bucket)pnm_malloc(sizeof(struct uint16_image_RGB_bucket_item));
                chl->ch.color = p;
                chl->ch.value = 1;
                chl->next = cht[hash];
                cht[hash] = chl;
              }
          }
      }
    return cht;
  }

uint16_image_RGB_table uint16_image_RGB_table_alloc(void)
  {
    uint16_image_RGB_table cht;
    int i;

    cht = (uint16_image_RGB_table)pnm_malloc(HASH_SIZE * sizeof(uint16_image_RGB_bucket));

    for (i = 0; i < HASH_SIZE; ++i) { cht[i] = (uint16_image_RGB_bucket)NULL; }

    return cht;
  }

int uint16_image_RGB_table_add(uint16_image_RGB_table cht, ppm_pixel_t* colorP, int value)
  {
    register int hash;
    register uint16_image_RGB_bucket chl;

    chl = (uint16_image_RGB_bucket)pnm_malloc(sizeof(struct uint16_image_RGB_bucket_item));
    hash = uint16_image_RGB_table_hash_pixel(colorP);
    chl->ch.color = *colorP;
    chl->ch.value = value;
    chl->next = cht[hash];
    cht[hash] = chl;
    return 0;
  }

uint16_image_RGB_hist_vector uint16_image_RGB_table_to_hist(uint16_image_RGB_table cht, int maxcolors)
  {
    uint16_image_RGB_hist_vector chv;
    uint16_image_RGB_bucket chl;
    int i, j;

    /* Now collate the hash table into a simple colorhist array. */
    chv = (uint16_image_RGB_hist_vector)pnm_malloc(maxcolors * sizeof(struct uint16_image_RGB_hist_item));

    /* Loop through the hash table. */
    j = 0;
    for (i = 0; i < HASH_SIZE; ++i)
      { for (chl = cht[i]; chl != (uint16_image_RGB_bucket)NULL; chl = chl->next)
          { /* Add the new entry. */
            chv[j] = chl->ch;
            ++j;
          }
      }
    /* All done. */
    return chv;
  }

uint16_image_RGB_table uint16_image_RGB_hist_to_table(uint16_image_RGB_hist_vector chv, int colors)
  {
    uint16_image_RGB_table cht;
    int i, hash;
    ppm_pixel_t color;
    uint16_image_RGB_bucket chl;

    cht = uint16_image_RGB_table_alloc();

    for (i = 0; i < colors; ++i)
      { color = chv[i].color;
        hash = uint16_image_RGB_table_hash_pixel(&color);
        for (chl = cht[hash]; chl != (uint16_image_RGB_bucket)NULL; chl = chl->next)
          { if (ppm_equal(&(chl->ch.color), &color))
              { pnm_error
                  ( "same color found twice - %u %u %u", 
                    color.c[0], color.c[1], color.c[2]
                  );
              }
          }
        chl = (uint16_image_RGB_bucket)pnm_malloc( sizeof(struct uint16_image_RGB_bucket_item));
        chl->ch.color = color;
        chl->ch.value = i;
        chl->next = cht[hash];
        cht[hash] = chl;
      }

    return cht;
  }

int uint16_image_RGB_table_lookup(uint16_image_RGB_table cht, ppm_pixel_t* colorP)
  {
    int hash;
    uint16_image_RGB_bucket chl;
    hash = uint16_image_RGB_table_hash_pixel(colorP);
    for (chl = cht[hash]; chl != (uint16_image_RGB_bucket)NULL; chl = chl->next)
      { if (ppm_equal(&(chl->ch.color), colorP))
          { return chl->ch.value; }
      }
    return -1;
  }

void uint16_image_RGB_table_free(uint16_image_RGB_table cht)
  {
    int i;
    uint16_image_RGB_bucket chl, chlnext;
    for (i = 0; i < HASH_SIZE; ++i)
      { for (chl = cht[i]; chl != (uint16_image_RGB_bucket)NULL; chl = chlnext)
          { chlnext = chl->next;
            free(chl);
          }
      }
    free(cht);
  }
