#ifndef jsstring_H
#define jsstring_H

#include <stdint.h>
#include <bool.h>

/* J. Stolfi's miscellaneous string utilities. */
/* Last edited on 2025-03-13 08:54:50 by stolfi */

typedef char* string_t;
  /* A more logical name for strings. */

int32_t isprefix(const string_t s, const string_t t);
  /* Returns TRUE iff the string {*s} is a prefix of {*t},
    including the cases when {*s} is empty or is equal to {*t}. */

string_t prefix(const string_t s, int32_t len);
  /* Returns a newly allocated string which is a copy
    of the first {len} characters of {s}.  The {len}
    must be non-negative and no more than {strlen(s)}. */

string_t txtcat (const string_t a, const string_t b);
  /* Returns a string that is the concatenation of {*a} and {*b}. 
    The result is always a newly allocated string, even if {a} and/or {b}
    are empty or {NULL}. */
  
string_t txtcat3 (const string_t a, const string_t b, const string_t c);
  /* Returns a string that is the concatenation of {*a,*b,*c}. The
    result is always a newly allocated string, even if any or all
    arguments are empty or {NULL}. */
  
string_t txtcat4 (const string_t a, const string_t b, const string_t c, const string_t d);
  /* Returns a string that is the concatenation of {*a,*b,*c,*d}. 
    The result is always a newly allocated string, even if any or all
    arguments are empty or {NULL}. */

string_t txtrep(const string_t  x, uint32_t n);
  /* Return the string {x} concatenated with itself {n} times. */
  
string_t add_ext(const string_t name, const string_t ext);
  /* Appends the string {*ext} (which should normally start with a
    period) to the given file {*name}. However, if {*name} is empty or
    "-", ignores {*ext} and returns a copy of {*name} itself. In any
    case, the result is a newly allocated string.*/

string_t trim_spaces(string_t x, bool_t at_beg, bool_t at_end);
  /* Returns a newly allocated string {y} that is a copy of {x},
    minus any space charaters at the beginning (if {at_beg} is true)
    and/or end (if {at_end} is true.  */

string_t fmt_int(int64_t x, uint32_t wid);
  /* Returns a newly allocated string containing the integer {x}
    converted to decimal and zero-padded to {wid} characters. */

string_t escapify(string_t x);
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
