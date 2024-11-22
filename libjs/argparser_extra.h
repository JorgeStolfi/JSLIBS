/* argparser_extra.h -- internal functions for {argparser.h}. */
/* Last edited on 2024-11-15 20:38:50 by stolfi */

#ifndef argparser_extra_H
#define argparser_extra_H

/* Copyright © 2005 by the State University of Campinas (UNICAMP).*/
/* Last edited on 2006-04-01 10:31:38 by stolfi*/

#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <affirm.h>

#include <argparser.h>

bool_t argparser_seems_keyword(char *n);
  /* Returns true iff the string {n} looks superficially like a UNIX keyword argument,
    namely begins with "-{X}" where {X} is an ASCII letter, or with "--". */
    
bool_t argparser_key_matches(char *a, char *b);
  /* Returns {TRUE} if the strings {a} and {b} are equal,
    or differ only by replacement of an initial "--" by "-". */ 

int64_t argparser_parse_int_string(argparser_t *pp, uint32_t index, char *arg, int64_t min, int64_t max, char **rest_P);
  /* Parses the beginning of the given string {arg}
    as a 64-bit signed integer, which must be in the range {min..max}.
    The string must be an optional '+' or '-' followed by at least one
    decimal digit. Otherwise the procedure fails with an error
    message, which assunes that {arg} is {pp.arg[index]} or part 
    thereof.
    
    If {rest_P} is {NULL}, the entire string must be a valid number.
    Otherwise, the procedure returns in {*rest_P} a pointer to the
    first character of {arg} that was not parsed as part of the integer;
    which may be a pointer to the final '\000' character of {arg}. */

void argparser_arg_msg(argparser_t *pp, char *msg1, uint32_t index, char *msg2, char *val);
  /* Writes to {pp->wr} the program name as "{name}: {msg1}".
    Then, if argument {index} is between 1 and {pp.arg.ne}, writes "argument {index} = "{argv[index]}".
    Then writes {val} using {msg2} as the format. */

void argparser_error_at(argparser_t *pp, char *msg, char* pos, uint32_t index);
  /* Prints the error message {msg} with a "**" prefix.
    Then, if {index} is between 1 and {pp.arg.ne} (i.e. main's {argc})
    prints also the string {pos}, which should be "at" or "after",
    and the argument {pp.arg.e[index]}. Then aborts the program 
    with status 1. */

void argparser_print_text(FILE *wr, char *info, int32_t wd);
  /* Prints a string {info} (terminated by '\000') to {wr}. The text
    may have multiple lines, separated by '\n'. Each line is
    formatted {argparser_print_info_line} to the requested width {wd}. */

void argparser_print_info(argparser_t *pp, int32_t wd);
  /* Prints the info texts stored in {pp->help}. Each text is
    formatted with {argparser_print_info_text} to the requested
    width {wd}, and terminated by a newline. */

void argparser_print_help_and_halt(argparser_t *pp, int32_t status);
  /* Prints the help texts stored in {pp->help}, and exits the program
    with the given {status} code. Each text is terminated by a
    newline. */
 
#endif
