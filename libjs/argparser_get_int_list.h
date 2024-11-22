/* argparser_get_int_list.h -- Parsing lists of integers from command line. */
/* Last edited on 2024-11-16 00:44:26 by stolfi */

#ifndef argparser_get_int_list_H
#define argparser_get_int_list_H

/* Copyright © 2023 by the State University of Campinas (UNICAMP).*/
/* Last edited on 2006-04-01 10:31:38 by stolfi*/

#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <affirm.h>
#include <argparser.h>

int64_vec_t argparser_get_int_list
  ( argparser_t *pp,
    uint32_t n_max,
    char *key,
    int64_t min,
    int64_t max
  );
  /* Parses all (zero or more) unparsed occurrences of the keyword
    {key}, not necessarily in consecutive positions, in order. Each keyword must
    be followed by one or more non-keyword arguments, and each of these should be an
    integer, an integer range "{a}-{b}" or "{a}..{b}", or one or more
    integers or ranges separated by commas. The integers thus specified,
    in the order given, are returne as an {it64_vec_t}.
    
    If {a} is greater than {b}, the range "{a}-{b}" or "{a}..{b}"
    specifies no integers; otherwise it stands for all the integers from
    {a} to {b}, inclusive. All the integers appearing in those arguments
    must be in the range {min..max}, and there must be no more than
    {n_max} integers in total.
    
    Integers may be repeated and ranges may overlap. Thus, for example, 
    if the command line has 
    
      "... -lines 7 2,5..8,3 6-9 ... -lines 4,6..5 7,1 ..."
      
    then {argparser_get_int_list(pp, "-lines", ...)} will
    return the list {{7,2,5,6,7,8,3,6,7,8,9,4,7,1}}.
 
    If {key} begins with "-", the procedure 
    also accepts the version with "--",  and vice-versa. */

#endif
