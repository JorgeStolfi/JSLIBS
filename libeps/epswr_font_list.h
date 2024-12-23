/* Tools for handling an expandable font list. */
/* Last edited on 2024-12-05 10:13:55 by stolfi */

#ifndef epswr_font_list_H
#define epswr_font_list_H
  
#include <stdio.h>
#include <stdint.h>

#include <bool.h>

/* 
  The font list is an array {fonts[]} of strings, where {fonts =
  (*fontsP)}. 
  
  The entries in use are {fonts[0..nFonts-1]}, where {nFonts
  =(*nFontsP)}. */
  
void epswr_font_list_initialize(int32_t *nFontsP, char ***fontsP);
  /* Intialize {nFontsP} to zero and {*fontsP} to the address 
    of a newly allocated array of suitable size. */

void epswr_font_list_add(const char *font, int32_t *nFontsP, char ***fontsP);
  /* If the string {font} is not in the font list {fonts[0..nFonts-1]}, appends a copy of 
    it to the list, incementing {*nFontsP}. Reallocates {*fontsP}
    if necessary. */

void epswr_font_list_free(int32_t *nFontsP, char ***fontsP);
  /* Reclaims the storage of all strings in {fonts[0.. */
 
#endif
