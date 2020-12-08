#ifndef pst_argparser_H
#define pst_argparser_H

/* pst_argparser.h -- parsing photostereo-related command line arguments. */
/* Last edited on 2006-05-02 12:06:06 by stolfi */

#include <argparser.h>
#include <bool.h>
#include <vec.h>
#include <r2.h>
    
#include <pst_basic.h>

bool_t pst_keyword_present(argparser_t *pp, char *key, bool_t next);
  /* If {next} is true, same as {argparser_keyword_present_next}
    (checks for unparsed {key} only at the next command line argument).
    If {next} is FALSE, same as {argparser_keyword_present}
    (check for unparsed {key} anywhere in the command line). */

char *pst_parse_next_file_name(argparser_t *pp);
  /* Tries to parse the next command line argument as a file name.
    If the next argument exists, is still unparsed, and looks like 
    a file name, marks it parsed, increments {pp->next}, and returns
    that argument.  Otherwise does none of those things, and returns NULL.
    
    An argument looks like a file name if it is just "-", or begins
    with letter, digit, "@", "/" or ".". */

name_vec_t pst_parse_file_name_list(argparser_t *pp, int *NNP);
  /* Parses the next zero or more consecutive command line arguments
    as file names. Accepts only arguments that exist, have not been
    parsed, and look like file names, as defined by
    {pst_parse_next_file_name}.
    
    If {*NNP} is zero or positive, parses only the next {*NNP}
    arguments as names, and fails if there aren't that many names. If
    {*NNP} is negative, parses as many arguments as possible, until
    the first argument that cannot be parsed, and sets {*NNP} to the
    number of arguments actually parsed. */

#endif
