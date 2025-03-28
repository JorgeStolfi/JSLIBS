/* See epswr.h */
/* Last edited on 2024-12-05 10:13:52 by stolfi */

#include <stdio.h>
#include <stdint.h>
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

void epswr_font_list_initialize(int32_t *nFontsP, char ***fontsP)
  {
    demand((*fontsP) == NULL, "non-NULL font list");
    int32_t nFonts = 0;
    int32_t nAlloc = epswr_font_list_MINSZ;
    char **fonts = talloc(nAlloc, char*);
    /* Update caller's variables: */
    (*nFontsP) = nFonts;
    (*fontsP) = fonts;
  }

void epswr_font_list_add(const char *font, int32_t *nFontsP, char ***fontsP)
  {
    int32_t nFonts = (*nFontsP);
    char **fonts = (*fontsP);
    demand(fonts != NULL, "font list not initialized");
    /* Check whether font already is in the table: */
    for (int32_t i = 0;  i < nFonts; i++)
      { if (strcmp(fonts[i], font) == 0) { return; } }
    /* Ensure that there is space for a new entry: */
    if ((nFonts >= epswr_font_list_MINSZ) && (((nFonts - 1) & nFonts) == 0))
      { /* Table is full ({nFonts} is a power of 2), reallocate it. */
        fonts = retalloc(fonts, 2*nFonts, char*); 
        affirm(fonts != NULL, "out of mem"); 
      }
    /* Store a copy of the name: */
    fonts[nFonts] = txtcat((char*)font, "");
    nFonts++;
    /* Update caller's variables: */
    (*nFontsP) = nFonts;
    (*fontsP) = fonts;
  }

void epswr_font_list_free(int32_t *nFontsP, char ***fontsP)
  {
    int32_t nFonts = (*nFontsP);
    char **fonts = (*fontsP);
    demand(fonts != NULL, "font list not initialized");
    for (int32_t i = 0;  i < nFonts; i++) { free(fonts[i]); }
    free(fonts);
    (*fontsP) = NULL;
    (*nFontsP) = 0;
  }
 
