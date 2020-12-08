#ifndef jsstring_H
#define jsstring_H

#include <stdint.h>

/* J. Stolfi's miscellaneous string utilities. */
/* Last edited on 2019-12-14 16:55:19 by jstolfi */

int32_t isprefix(const char *s, const char *t);
  /* Returns TRUE iff the string {*s} is a prefix of {*t},
    including the cases when {*s} is empty or is equal to {*t}. */

char *txtcat (const char *a, const char *b);
  /* Returns a string that is the concatenation of {*a} and {*b}. 
    The result is always a newly allocated string. */
  
char *txtcat3 (const char *a, const char *b, const char *c);
  /* Returns a string that is the concatenation of {*a,*b,*c}. 
    The result is always a newly allocated string. */
  
char *addext(const char *name, const char *ext);
  /* Appends the string {*ext} (which should normally start with a
    period) to the given file {*name}. However, if {*name} is empty or
    "-", ignores {*ext} and returns a copy of {*name} itself. In any
    case, the result is a newly allocated string.*/

char *txtrep(const char* x, uint32_t n);
  /* Return the string {x} concatenated with itself {n} times. */

char *fmtint(int64_t x, uint32_t wid);
  /* Returns a newly allocated string containing the integer {x}
    converted to decimal and zero-padded to {wid} characters. */

#endif
