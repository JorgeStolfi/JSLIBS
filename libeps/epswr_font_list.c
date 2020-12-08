/* See epswr.h */
/* Last edited on 2020-01-11 01:37:22 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#include <bool.h>
#include <affirm.h>
#include <jsstring.h>

#include <epswr_font_list.h>

#define epswr_font_list_MINSZ 8
  /* The allocated size of the {fonts[]} array is the number {nFonts} of
    elements in use rounded up to the next power of two, or
    {epswr_font_list_MINSZ}, whichever is larger. */

void epswr_font_list_initialize(int *nFontsP, char ***fontsP)
  {
    demand((*fontsP) == NULL, "non-NULL font list");
    int nFonts = 0;
    int nAlloc = epswr_font_list_MINSZ;
    char **fonts = (char **)malloc(nAlloc*sizeof(char *));
    /* Update caller's variables: */
    (*nFontsP) = nFonts;
    (*fontsP) = fonts;
  }

void epswr_font_list_add(const char *font, int *nFontsP, char ***fontsP)
  {
    int nFonts = (*nFontsP);
    char **fonts = (*fontsP);
    demand(fonts != NULL, "font list not initialized");
    /* Check whether font already is in the table: */
    for (int i = 0; i < nFonts; i++)
      { if (strcmp(fonts[i], font) == 0) { return; } }
    /* Ensure that there is space for a new entry: */
    if ((nFonts >= epswr_font_list_MINSZ) && (((nFonts - 1) & nFonts) == 0))
      { /* Table is full ({nFonts} is a power of 2), reallocate it. */
        fonts = realloc(fonts, 2*nFonts*sizeof(char*)); 
        affirm(fonts != NULL, "out of mem"); 
      }
    /* Store a copy of the name: */
    fonts[nFonts] = txtcat(font, "");
    nFonts++;
    /* Update caller's variables: */
    (*nFontsP) = nFonts;
    (*fontsP) = fonts;
  }

void epswr_font_list_free(int *nFontsP, char ***fontsP)
  {
    int nFonts = (*nFontsP);
    char **fonts = (*fontsP);
    demand(fonts != NULL, "font list not initialized");
    for (int i = 0; i < nFonts; i++) { free(fonts[i]); }
    free(fonts);
    (*fontsP) = NULL;
    (*nFontsP) = 0;
  }
 
