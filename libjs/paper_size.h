/* paper_size.h - dimensions of standard paper sizes. */

#ifndef paper_size_H
#define paper_size_H

#define _GNU_SOURCE

#include <bool.h>

void paper_size_get_dimensions
 ( const char *paperSize, 
   bool_t landscape, 
   double *xpt, 
   double *ypt
 );
  /* Sets *xpt and *ypt to the dimensions of the specified paper
    type, in points. Knows about US sizes "letter", "ledger",
    "tabloid", "legal", "executive", and the ISO "A" sizes (from
    "4A0" to "A10").
    
    If {landscape} is false, assumes that the paper is in "portrait"
    orientation, witht {*xpt <= *ypt}. If {landscape} is true,
    assumes the "landscape" orientation, thus swapping {*xpt} and
    {*ypt} . */

#endif
