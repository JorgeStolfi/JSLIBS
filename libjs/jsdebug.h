#ifndef jsdebug_H
#define jsdebug_H

#include <stdint.h>
#include <bool.h>

/* J. Stolfi's miscellaneous debugging utilities. */
/* Last edited on 2024-11-16 06:33:14 by stolfi */

    
void jsdebug_addr_span(char *name, void *ini, void *fin, uint32_t n);
  /* Print addresses and byte count between {ini} and {fin}, and the
    element byte size assuming that the area consists of {n} elements.
    The {name} is just printed at the start. The addresses {ini} and
    {fin} must be different, and {n} must be positive. */

#endif
