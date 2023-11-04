#ifndef jsstring_H
#define jsstring_H

#define _GNU_SOURCE
#include <stdint.h>
#include <bool.h>

/* J. Stolfi's miscellaneous string utilities. */
/* Last edited on 2023-10-15 03:26:18 by stolfi */

int32_t isprefix(const char *s, const char *t);
  /* Returns TRUE iff the string {*s} is a prefix of {*t},
    including the cases when {*s} is empty or is equal to {*t}. */

char *txtcat (const char *a, const char *b);
  /* Returns a string that is the concatenation of {*a} and {*b}. 
    The result is always a newly allocated string. */
  
char *txtcat3 (const char *a, const char *b, const char *c);
  /* Returns a string that is the concatenation of {*a,*b,*c}. 
    The result is always a newly allocated string. */
  
char *txtcat4 (const char *a, const char *b, const char *c, const char *d);
  /* Returns a string that is the concatenation of {*a,*b,*c,*d}. 
    The result is always a newly allocated string. */

char *txtrep(const char* x, uint32_t n);
  /* Return the string {x} concatenated with itself {n} times. */
  
char *add_ext(const char *name, const char *ext);
  /* Appends the string {*ext} (which should normally start with a
    period) to the given file {*name}. However, if {*name} is empty or
    "-", ignores {*ext} and returns a copy of {*name} itself. In any
    case, the result is a newly allocated string.*/

char *trim_spaces(char *x, bool_t at_beg, bool_t at_end);
  /* Returns a newly allocated string {y} that is a copy of {x},
    minus any space charaters at the beginning (if {at_beg} is true)
    and/or end (if {at_end} is true.  */

char *fmt_int(int64_t x, uint32_t wid);
  /* Returns a newly allocated string containing the integer {x}
    converted to decimal and zero-padded to {wid} characters. */

char *escapify(char *x);
  /* Returns a newly allocated copy of {x}, with any 
    special characters replaced by an escape code starting with a backslash.
    Specifically, 
    
      ASCII LF (end-of-line, '\012') is replaced by backslash-'n';

      ASCII TAB ('\011') is replaced by backslash-'t';
      
      backslash ('\134') is replaced by double backslash; 
      
      any character between ' ' (ASCII space, '\040') and '~' ('\126'), inclusive, is
      copied without change;
      
      any other character (from '\001' to '\277') is replaced by backslash and 
      its three-digit octal code.
      
   */ 

#endif
