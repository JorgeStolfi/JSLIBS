/* argparser.h -- facilities for parsing command line arguments. */
/* Last edited on 2024-11-22 02:18:03 by stolfi */

#ifndef argparser_H
#define argparser_H

/* Copyright © 2005 by the State University of Campinas (UNICAMP).*/
/* See the copyright, authorship, and warranty notice at end of file.*/
/* Last edited on 2006-04-01 10:31:38 by stolfi*/

/* This interface provides simple and robust tools for parsing the
  command line arguments given to a process when it is started. Check
  the usage example at the end of this interface. */

#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <affirm.h>

typedef struct argparser_t /* A parser for command line arguments. */
  { string_vec_t arg;      /* Command line arguments; {arg[0]} is prog name. */
    bool_vec_t parsed;     /* parsed[i] is {TRUE} if {arg[i]} has been parsed. */
    int32_t next;          /* The next argument to parse is {arg[next]} */
    FILE *wr;              /* File for errors */
    uint32_t nhelp;        /* Number of used lines of error help text. */
    string_vec_t help;     /* {help.e[0..nhelp-1]}  is the help text for errors. */
    uint32_t ninfo;        /* Number of used lines of program info text. */
    string_vec_t info;     /* {info.e[0..ninfo-1]}  is the program info text. */
  } argparser_t;
  
argparser_t *argparser_new(FILE *wr, int32_t argc, char **argv);
  /* Saves pointers to the given command line arguments. Marks the
   command name {arg[0]} as parsed, all other arguments as unparsed.
   The next argument to be parsed will be {arg[1]}. Any parsing
   errors will be printed to {wr}. */
 
void argparser_error(argparser_t *pp, char *msg);
  /* Prints the given message, the help text, and terminates the program
    with status 1. */

void argparser_finish(argparser_t *pp);
  /* Checks whether all parameters have been parsed; if not, prints a
    message and raises error. Also reclaims {*pp} and all its internal
    storage. */

/* PROGRAM DOCUMENTATION OPTIONS */

void argparser_set_help(argparser_t *pp, char *help);
  /* Appends {help} to the help text stored in {pp}.
  
    The help text will be used by
    {argparser_process_help_info_options} below. Also, if a syntax
    error is found later during argument parsing, the stored help text
    will be written to {pp->wr}, just before halting the program. */

void argparser_set_info(argparser_t *pp, char *info);
  /* Appends {info} to the info text stored in {pp}.
  
    The info text will be used by {argparser_process_help_info_options}
    below.  If and when the info text is printed, it is reformatted
    by breaking every line that is longer than 72 characters at
    internal blanks. Any blanks around the break points (or at the end
    of each line) are discarded, but the initial indentation of the
    line is preserved on each piece. */

void argparser_process_help_info_options(argparser_t *pp);
  /*  If a documentation request keyword ("-info" or "--info") is
    present, the function prints the info text stored in {pp}. Else,
    if a help request keyword ("-help" or "--help") is present, it
    prints the stored help text. In either case, exits the program
    with status 0. */

#define argparser_help_info_HELP \
  "[ -help | --help ] [ -info | --info ]"

#define argparser_help_info_HELP_INFO \
  "  -help\n" \
  "  --help\n" \
  "    Prints an options summary and exits.\n" \
  "\n" \
  "  -info\n" \
  "  --info\n" \
  "    Prints this manpage and exits.\n" \
  "\n" \
  "  For compatibility with GNU/Linux tradition, any keyword that" \
  " starts with \"-\" may also start with \"--\", and vice-versa.  Thus \"-size\" is the" \
  " same as \"--size\", and \"-v\" is the same as \"--v\".  However the legacy Unix" \
  " practice of condensing single-character keywords is NOT" \
  " supported: \"-rs\" (and \"--rs\") is always considered a single keyword," \
  " not equivalent to \"-r -s\". "

#define argparser_help_info_NO_WARRANTY \
  "  This software is provided \"as is\", WITHOUT ANY EXPLICIT OR" \
  " IMPLICIT WARRANTY, not even the implied warranties of merchantibility" \
  " and fitness for a particular purpose."

#define argparser_help_info_STANDARD_RIGHTS \
  "  Permission to use, copy, modify, and redistribute this software and" \
  " its documentation for any purpose is hereby granted, royalty-free," \
  " provided that: (1) the copyright, AUTHOR, WARRANTY and RIGHTS notices" \
  " in the source files are retained or replaced by completely free" \
  " alternatives, such as the GNU Public  License (GPL);" \
  " (2) no executable code derived from this file is published" \
  " or distributed without the corresponding source code; and (3) these" \
  " same rights are granted to any recipient of such code, under" \
  " the same conditions."

/* PARSING FREE-ORDER NAMED OPTIONS AND ARGUMENTS */

bool_t argparser_keyword_present(argparser_t *pp, char *key);
  /*  Looks for the first unparsed argument {arg[i]} that is equal to
    {key}. If found, marks it as parsed, sets {pp->next} to {i+1}, and
    returns {TRUE}. Otherwise returns {FALSE} and leaves {pp->next}
    unchanged.  If {key} begins with "-", also accepts the version with "--",
    and vice-versa. */

void argparser_get_keyword(argparser_t *pp, char *key);
  /* Same as {argparser_keyword_present}, but raises error if the 
    keyword is not found. */

/* PARSING PARAMETERS AFTER KEYWORDS */

char *argparser_get_next(argparser_t *pp);
  /* Returns {arg[pp->next]}, marks it as parsed and increments {pp->next}.  
    Raises error if {arg[pp->next]} does not exist or has already 
    been parsed. */

char *argparser_get_next_keyword(argparser_t *pp);
  /* Returns {arg[pp->next]}, marks it as parsed and increments {pp->next}.  
    Raises error if {arg[pp->next]} does not exist, has already 
    been parsed, or does not look like a keyword (as
    per {argparser_next_is_keyword}) */

char *argparser_get_next_non_keyword(argparser_t *pp);
  /* Returns {arg[pp->next]}, marks it as parsed and increments {pp->next}.  
    Raises error if {arg[pp->next]} does not exist, has already 
    been parsed, or looks superficially like a keyword (as
    per {argparser_next_is_keyword}). */

int64_t argparser_get_next_int(argparser_t *pp, int64_t min, int64_t max);
uint64_t argparser_get_next_uint(argparser_t *pp, uint64_t min, uint64_t max);
double argparser_get_next_double(argparser_t *pp, double min, double max);
  /* Same as {argparser_get_next}, but converts the result to the
    approriate type (using {strtol} and {strtod}, respectively).
    Raises error if the parameter is not a valid literal, or lies
    outside of the range {[min..max]}.  */

bool_t argparser_get_next_bool(argparser_t *pp);
  /* Same as {argparser_get_next}, but converts the argument string to a
    boolean. The strings "t", "T", "TRUE", "true", "y", "Y", "yes", and "YES", and the numeric value
    1 are converted to {TRUE}; the strings "f", "F", "false", "FALSE", "n", "N","no", and "NO",
    and the numeric value 0 are converted to {FALSE}. Any other value is
    an error. */

/* PARSING SPECIAL SYNTAX */
  
char *argparser_next(argparser_t *pp);
  /* Returns {arg[pp->next]} if it exists and is still 
    unparsed; else return NULL. Does not change {pp->next}
    and does not mark that argument as parsed. */

bool_t argparser_next_is_keyword(argparser_t *pp);
  /* Returns TRUE if and only if {arg[pp->next]} exists, is still
    unparsed, and looks superficially like a keyword --- that is,
    begins with "-{X}" where {X} is an ASCII letter, or with "--".
    Does not change {pp->next} and does not mark that argument as parsed. */

bool_t argparser_next_is_non_keyword(argparser_t *pp);
  /* Returns TRUE if and only if {arg[pp->next]} exists, is still
    unparsed, and does NOT look like a keyword as per above. Does not
    change {pp->next} and does not mark that argument as parsed. Note
    that it is not the same as {(! argparser_next_is_keyword(pp))}. */

bool_t argparser_next_is_number(argparser_t *pp);
  /* Returns TRUE if and only if {arg[pp->next]} exists, is still
    unparsed, and looks superficially like a number --- namely, begins
    with a digit, or "+" or "-" followed by a digit. Does not change
    {pp->next} and does not mark that argument as parsed. */

bool_t argparser_is_next(argparser_t *pp, char *key);
  /* Returns TRUE if and only if {arg[pp->next]} exists, is still 
    unparsed, and is equal to {key}. Does not change {pp->next}
    and does not mark that argument as parsed.  If {key} 
    begins with "-", also accepts the version with "--",
    and vice-versa. */

void argparser_skip_parsed(argparser_t *pp);
  /* Points {pp->next} at the first unparsed argument. If there are
    any parsed arguments beyond that one, prints a message and raises
    error. */

/* PARSING FIXED-POSITION KEYWORDS AND DELIMITERS */

bool_t argparser_keyword_present_next(argparser_t *pp, char *key);
  /* If {argparser_is_next(pp, key)} is true, marks the next argument
    as parsed, increments {pp->next} and returns TRUE. Otherwise does
    none of these things and returns {FALSE}. */

void argparser_get_keyword_next(argparser_t *pp, char *key);
  /* If {argparser_is_next(pp, key)} is true, marks the next argument
    as parsed and increments {pp->next}. Otherwise raises an error. */

/*
  Most Unix programs expect their command-line arguments to consist
  of a string of keywords and keyword-labeled arguments (`options',
  `switches', etc.), followed by a list of positional arguments.

  For the user's convenience, programs generally allow the switches
  and keyword-labeled arguments to be given in any order. Some of
  those parameters may be optional and/or repeatable, some may be
  mandatory; some may be required or forbidden depending on the values
  of the other parameters. Furthermore, the value of an argument may
  be just a number or a text string, or may be a cluster of two or
  more values with their own little syntax.

  This module simplifies the parsing of such command-line parameters,
  by allowing the program to scan the arguments in the order which
  most suits the program. This module also detects automatically many
  kinds of common mistakes --- such as arguments that are missing, repeated,
  extraneous, malformed, or out of range --- and prints the appropriate
  error messages.

  For example, here is how this module could be used by an
  hypothetical program {prt} that concatenates a bunch of files and
  prints selected line ranges of the result, possibly in reverse
  order, with several formatting options.

    #define MaxLines MAX_INT
    #define MaxRanges 100
    #define MaxFiles 100
    #define MinFontSize 1
    #define MaxFontSize 100
    
    / * Arguments from command line: * /
    int32_t fontSize;
    bool_t landscape;
    int32_t nRanges = 0;
    int32_t ini[MaxRanges], fin[MaxRanges];
    bool_t reverse[MaxRanges];
    int32_t nFiles = 0;
    char *files[MaxFiles];
    
    void parse_args(int32_t argc, char **argv)
      {
        static char *help = 
          "prt \\\n"
          "  -fontSize NUM \\\n"
          "  [ -landscape | -portrait ] \\\n"
          "  [ -lines NUM NUM [ -reverse ] ]..."
          "  FNAME...";
    
        / * Initialize the argument parser: * /
        argparser_t *pp = argparser_new(stderr, argc, argv);
        argparser_set_help(pp, help);
    
        / * The "-fontSize" parameter is mandatory: * /
        argparser_get_keyword(pp, "-fontSize");
        fontSize = (int32_t)argparser_get_next_int(pp, MinFontsize, MaxFontSize);
    
        / * Either "-landscape" or "-portrait", but not both: * /
        if (argparser_keyword_present(pp, "-landscape"))
          { landscape = TRUE; } 
        else if (argparser_keyword_present(pp, "-portrait"))
          { landscape = FALSE; }
        else
          { / * Default is "-portrait" unless font is too big. * /
            landscape = (fontSize > 8);
          }
    
        / * Parse the line ranges: * /
        nRanges = 0;
        while (argparser_keyword_present(pp, "-lines"))
          { if (nRanges >= MaxRanges) 
              { argparser_error(pp, "Too many page ranges"); }
            ini[nRanges] = (int32_t)argparser_get_next_int(pp, 1,MaxLines);
            fin[nRanges] = (int32_t)argparser_get_next_int(pp, ini[nRanges],MaxLines);
            rev[nRanges] = argparser_keyword_present_next(pp, "-reverse");
            nRanges = nRanges+1;
          }
    
        / * By default, print all lines: * /
        if (nRanges == 0)
          { ini[0] = 1; fin[0] = MaxLines; rev[0] = FALSE;
            nRanges = 1;
          }
    
        / * Parse the file list (after all keyword options): * /
        argparser_skip_parsed(pp);
        nFiles = argc - pp->next;
        if (nFiles == 0)
          { argparser_error(pp, "no files specified"); }
        for (i = 0; i < nFiles; i++)
          { files[i] = argparser_get_next(pp); }
    
        / * Check for any unparsed parameters: * /
        argparser_finish(pp);
      }

  Note that this code allows the user to give the options
  "-fontSize" and "-landscape"/"-portrait" in any order, even
  anywhere among or after the "-range" arguments. However, each
  "-range" flag must be immediately followed by two numbers; and the
  "-reverse" flag, if given, must immediately follow the second
  number.  Also, all file names must follow all the other options
  and arguments.
*/  

/* Copyright © 2003 by State University of Campinas (UNICAMP).
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appears in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty of any kind.
*/

#endif
