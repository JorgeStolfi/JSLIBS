#ifndef jsdebug_H
#define jsdebug_H

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

/* J. Stolfi's miscellaneous debugging utilities. */
/* Last edited on 2024-06-28 02:08:08 by stolfi */

    
void jsdebug_addr_span(char *name, void *ini, void *fin, int32_t n);
  /* Print addresses and byte count between {ini} and {fin}, and the
    element byte size assuming that the area consists of {n}
    elements.  The {name} is just printed at the start. */

#endif
