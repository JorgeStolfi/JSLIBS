/* Private clone of {asprintf} from {libiberty.h} */
/* Last edited on 2024-11-20 03:59:26 by stolfi */

#ifndef jsprintf_H
#define jsprintf_H

#include <stdint.h>

char* jsprintf (const char *fmt, ...);
  /* Converts zero or more arguments of arbitrary types (the {...}) to a
    string, according to the format string {fmt}. Like {sprintf} (q.v.)
    but returns a pointer to a new storage area internally allocated
    from the heap.
    
    This function will compute the size of the buffer needed, allocate
    it {malloc}, write into it the {...} arguments formatted according
    to {fmt}, and return the address of that buffer as a result. This
    pointer is always not {NULL} and safe to reclaim with {free}.
    
    Unlike GNU's {asprintf} from {libiberty.h}, this version bombs out
    if memory cannot be allocated or an error occurs in the
    conversion. */

#endif
