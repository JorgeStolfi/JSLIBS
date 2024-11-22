/* See argparser_get_int_list.h. */
/* Last edited on 2024-11-22 03:15:34 by stolfi */

/* Copyright © 2003 Jorge Stolfi, Unicamp. See note at end of file. */
/* Based on Params.m3 by J.Stolfi, DEC-SRC, 1988.  */

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include <vec.h>
#include <affirm.h>
#include <argparser.h>
#include <argparser_extra.h>

#include <argparser_get_int_list.h>

/* STRING PARSING

  The following procedures parse the next command line args
  or a given string {arg} as the specification of a set of integers.

  Each integer thus specified must be in the range {min..max}. All those
  integers are appended to the vector {num} (which will ne expanded as
  needed), starting at position {*n_P}.

  On exit, {*n_P} will be incremented by the number of elements
  appendded to {num}. The procedure fails if {*n_P} would exceed {n_max}.

  In procedures that take a string {arg}, eventual error messages assume
  that {arg} is {pp.arg[index]} or part thereof. */

void argparser_get_next_int_groups
  ( argparser_t *pp, 
    uint32_t n_max, 
    uint32_t *n_P, 
    int64_vec_t *num, 
    int64_t min, 
    int64_t max
  );
  /* Parses zero or more command-line arguments,
    starting with {pp.next}, that are integers,
    integer ranges, or lists of integer ranges
    separared by ",". */

void argparser_parse_int_group_string
  ( argparser_t *pp, 
    int32_t index,
    char *arg, 
    uint32_t n_max, 
    uint32_t *n_P, 
    int64_vec_t *num, 
    int64_t min, 
    int64_t max
  );
  /* Parses the given string {arg} as one or more integers or integer
    ranges separated by ",". See {argparser_parse_int_range_string} for
    the syntax of ranges. */
    
void argparser_parse_int_range_string
  ( argparser_t *pp, 
    int32_t index,
    char *arg, 
    uint32_t n_max, 
    uint32_t *n_P, 
    int64_vec_t *num, 
    int64_t min, 
    int64_t max
  );
  /* Parses the given string {arg} as an integer, or as an integer 
    range "{a}..{b}" or "{a}-{b}".  If {a} is greater than {b},
    the range denotes an empty set. Otherwise it denotes
    all integers between {a} and {b} inclusive. */
    
/* IMPLEMENTATIONS */    

int64_vec_t argparser_get_int_list
  ( argparser_t *pp,
    uint32_t n_max,
    char *key,
    int64_t min,
    int64_t max
  )
  { int64_vec_t num = int64_vec_new(10);
    uint32_t n = 0;
    while (argparser_keyword_present(pp, key))
      { argparser_get_next_int_groups(pp, n_max, &n, &num, min, max); }
    int64_vec_trim(&num, (uint32_t)n);
    return num;
  }

void argparser_get_next_int_groups
  ( argparser_t *pp, 
    uint32_t n_max, 
    uint32_t *n_P, 
    int64_vec_t *num, 
    int64_t min, 
    int64_t max
  )
  { while (argparser_next_is_number(pp))
      { char *arg = argparser_get_next_non_keyword(pp);
        char *pbeg = arg;
        assert((*pbeg) != '\000');
        while(TRUE)
          { char c = (*pbeg);
            if (c == '\000') 
              { /* If {arg} was, e.g., "2..5,": */
                argparser_error_at(pp, "missing integer or range", "at", pp->next-1);
              }
            if ((c != '+') && (c != '-') && ((c < '0') || (c > '9')))
              { /* Not an integer or integer range: */
                argparser_error_at(pp, "invalid integer", "at", pp->next-1);
              }
            /* Look for ",", if any: */
            char *pend = strchrnul(pbeg, ',');
            /* Isolate the next integer or range: */
            char c2 = (*pend);
            (*pend) = '\000';
            argparser_parse_int_range_string(pp, pp->next-1, pbeg, n_max, n_P, num, min, max);
            /* Are we done parsing this argument? */
            if (c2 == '\000') { break; }
            /* Get the next int or range: */
            assert(c2 == ',');
            pbeg = pend+1;
          }
      }
  }
  
void argparser_parse_int_group_string
  ( argparser_t *pp,
    int32_t index,
    char *arg, 
    uint32_t n_max, 
    uint32_t *n_P, 
    int64_vec_t *num, 
    int64_t min, 
    int64_t max
  )
  { char *pbeg = arg;
    assert((*pbeg) != '\000');
    while(TRUE)
      { char c = (*pbeg);
        if (c == '\000') 
          { /* If {arg} was, e.g., "2..5,": */
            argparser_error_at(pp, "missing integer or range", "at", index);
          }
        if ((c != '+') && (c != '-') && ((c < '0') || (c > '9')))
          { /* Not an integer or integer range: */
            argparser_error_at(pp, "invalid integer", "at", index);
          }
        /* Look for ",", if any: */
        char *pend = strchrnul(pbeg, ',');
        /* Isolate the next integer or range: */
        char c2 = (*pend);
        (*pend) = '\000';
        argparser_parse_int_range_string(pp, index, pbeg, n_max, n_P, num, min, max);
        /* Are we done parsing this argument? */
        if (c2 == '\000') { break; }
        /* Get the next int or range: */
        assert(c2 == ',');
        pbeg = pend+1;
      }
  }
  
void argparser_parse_int_range_string
  ( argparser_t *pp, 
    int32_t index,
    char *arg, 
    uint32_t n_max, 
    uint32_t *n_P, 
    int64_vec_t *num, 
    int64_t min, 
    int64_t max
  )
  {
    char *pbeg = arg;
    int64_t a, b; /* Range endpoints. */

    /* The next thing, if any must start like an integer: */
    char c = (*pbeg);
    assert((c == '+') || (c == '-') | ((c >= '0') && (c <= '9')));
    
    /* Parse the first part of the argument: */
    char *pend = NULL;
    a = argparser_parse_int_string(pp, index, pbeg, min, max, &pend);
    assert((pend != NULL) && (pend > pbeg));
    
    /* Check for "-" or "..": */
    char c2 = (*pend);
    if (c2 == '\000')
      { /* Single int, not range: */
        b = a;
      }
    else
      { /* Skip the "-" or "..": */
        pbeg = pend+1;
        if (c2 != '-')
          { /* Must be a "..": */
            if ((c2 != '.') || ((*pbeg) != '.'))
              { argparser_error_at(pp, "invalid range separator", "at", index); }
            pbeg++;
          }
        b = argparser_parse_int_string(pp, index, pbeg, min, max, &pend);
        assert((pend != NULL) && (pend > pbeg));
      }
    uint32_t n = (*n_P);
    for (int64_t v = a; v <= b; v++) 
      { if (n >= n_max) { argparser_error_at(pp, "too many integers", "at", index); }
        int64_vec_expand(num, n);
        num->e[n] = v;
        n++;
      }
    (*n_P) = n;
  }
