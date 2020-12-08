/* See paper_size.h */
/* Last edited on 2020-01-11 05:46:11 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>

#include <affirm.h>
#include <bool.h>

#include <paper_size.h>

#define mm (72.0/25.4)
  /* One millimeter in Postscript points. */

void paper_size_get_dimensions
  ( const char *paperSize, 
    bool_t landscape, 
    double *xpt, 
    double *ypt
  )
  { /* US sizes are defined in inches, which are 72 "pt" exactly: */
    /* ISO paper sizes are defined in "mm", and need conversion: */
    if (
        (strcmp(paperSize, "letter") == 0) ||
        (strcmp(paperSize, "Letter") == 0)
      )
      { /* 8.5 x 11 in, 216 x 279 mm */
        (*xpt) = 612.0; (*ypt) = 792.0;
      }
    else if (
       (strcmp(paperSize, "a4") == 0) ||
        (strcmp(paperSize, "A4") == 0)
      )
      { (*xpt) = 210 * mm; (*ypt) = 297 * mm; }
    else if (
        (strcmp(paperSize, "a3") == 0) ||
        (strcmp(paperSize, "A3") == 0)
      )
      { (*xpt) = 297 * mm; (*ypt) = 420 * mm; }
    else if (
        (strcmp(paperSize, "legal") == 0) ||
        (strcmp(paperSize, "Legal") == 0)
      )
      { /* 8.5 x 14 in, 216 x 356 mm */
         (*xpt) = 612.0; (*ypt) = 1008.0; 
      }
    else if (
        (strcmp(paperSize, "executive") == 0) ||
        (strcmp(paperSize, "Executive") == 0)
      )
      { /* 7.5 x 10 in, 190 x 254 mm */
        (*xpt) = 540.0; (*ypt) = 720.0; 
      }
    else if (
        (strcmp(paperSize, "ledger") == 0) ||
        (strcmp(paperSize, "Ledger") == 0) ||
        (strcmp(paperSize, "tabloid") == 0) ||
        (strcmp(paperSize, "Tabloid") == 0))
      { /* 11 x 17 in, 279 x 432 mm */
        (*xpt) = 792.0; (*ypt) = 1224.0; 
      }
    else if (
        (strcmp(paperSize, "a10") == 0) ||
        (strcmp(paperSize, "A10") == 0)
      )
      { (*xpt) = 26 * mm; (*ypt) = 37 * mm; }
    else if (
        (strcmp(paperSize, "a9") == 0) ||
        (strcmp(paperSize, "A9") == 0)
      )
      { (*xpt) = 37 * mm; (*ypt) = 52 * mm; }
    else if (
      (strcmp(paperSize, "a8") == 0) ||
      (strcmp(paperSize, "A8") == 0))
      { (*xpt) = 52 * mm; (*ypt) = 74 * mm; }
    else if (
        (strcmp(paperSize, "a7") == 0) ||
        (strcmp(paperSize, "A7") == 0)
      )
      { (*xpt) = 74 * mm; (*ypt) = 105 * mm; }
    else if (
        (strcmp(paperSize, "a6") == 0) ||
        (strcmp(paperSize, "A6") == 0)
      )
      { (*xpt) = 105 * mm; (*ypt) = 148 * mm; }
    else if (
        (strcmp(paperSize, "a5") == 0) ||
        (strcmp(paperSize, "A5") == 0)
      )
      { (*xpt) = 148 * mm; (*ypt) = 210 * mm; }
    else if (
        (strcmp(paperSize, "a2") == 0) ||
        (strcmp(paperSize, "A2") == 0)
      )
      { (*xpt) = 420 * mm; (*ypt) = 594 * mm; }
    else if (
        (strcmp(paperSize, "a1") == 0) ||
        (strcmp(paperSize, "A1") == 0)
      )
      { (*xpt) = 594 * mm; (*ypt) = 841 * mm; }
    else if (
        (strcmp(paperSize, "a0") == 0) ||
        (strcmp(paperSize, "A0") == 0)
      )
      { (*xpt) = 841 * mm; (*ypt) = 1189 * mm; }
    else if (
        (strcmp(paperSize, "2a0") == 0) ||
        (strcmp(paperSize, "2A0") == 0)
      )
      { (*xpt) = 1189 * mm; (*ypt) = 1682 * mm; }
    else if (
        (strcmp(paperSize, "4a0") == 0) ||
        (strcmp(paperSize, "4A0") == 0)
      )
      { (*xpt) = 1682 * mm; (*ypt) = 2378 * mm; }
    else
      { demand(FALSE, "unkown paper size"); }
    if (landscape)
      { /* Swap the dimensons: */
        double tmp = (*xpt); (*xpt) = (*ypt); (*ypt) = tmp;
      }
  }
