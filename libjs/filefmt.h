/* filefmt.h -- version-checked file headers, footers, and comments. */
/* Last edited on 2019-04-09 12:49:10 by jstolfi */

#ifndef filefmt_H
#define filefmt_H

#include <stdio.h>
#include <affirm.h>

/* FILE HEADERS AND FOOTERS */

/* These procedures help maintaining structured files with typed 
  and versioned headers and footers, like this:
  
     begin myimage (format of 1996-11-29)
     [...stuff...]
     end myimage
     
  The strings "myimage" and "1996-11-29" are client-specified,
  and checked on input.  These begin/end pairs may be nested.  
  Some file types may not have the part in parentheses,
  especially if they are nested inside other header/footer pairs. */

void filefmt_write_header(FILE *wr, char *type, char *version);
void filefmt_write_footer(FILE *wr, char *type);
  /* These procedures write a header or footer line to {wr}, 
    complete with the final newline.  The part "(format of {VERSION})" is 
    omitted if {version} is {NULL}. */

void filefmt_read_gen_header(FILE *rd, char **typeP, char **versionP);
void filefmt_read_gen_footer(FILE *rd, char **typeP);
  /* These procedures parse a header or footer from the given file, with
    any type and version, and returns the values of those fields in
    {*typeP} and {*versionP}. Those strings are newly allocated by the
    prcedure. If the part "(format of {VERSION})" is omitted 
    in the file, {*versionP} is set to {NULL}. 
    The procedures abort the program on any error.
    
    These procedures skip any formatting characters (SPACE, TAB, NUL,
    CR, LF, and page breaks), before the first token, and any number of
    spaces (SPACE, TAB, NUL) after each token; but require the entire
    header or footer to be contained in one line, and terminated by a
    newline (which is consumed too). The type and version fields must not contain
    embedded formatting characters. */

void filefmt_read_header(FILE *rd, char *type, char *version);
void filefmt_read_footer(FILE *rd, char *type);
  /* These procedures are similar to {filefmt_read_gen_header} and
    {filefmt_read_gen_footer}, but require the type and version to be
    those specified. If the {version} is {NULL}, the part "(format of
    {VERSION})" must be omitted in the file, otherwise it must be
    present. They abort the program on any error; in particular, if the
    header fields are of a different {type} or {version}. */

char *filefmt_make_header(char *type, char *version);
char *filefmt_make_footer(char *type);
  /* These procedures construct a header or footer line as a string, 
    including the terminating newline.  The part "(format of {VERSION})" is 
    omitted if {version} is {NULL}. */
      
/* COMMENTS: */

/* These routines write and parse a comment text, consisting of zero or
  more lines marked by a given {prefix} character:
  
    | Blah blah blah
    |   blah blah
    | and more blah.
  
  The {prefix} must be the first non-blank character of the line, and is
  normally followed by a blank. */

void filefmt_write_comment(FILE *wr, char *cmt, int ind, char prefix);
  /* Writes the given {cmt} text to {wr}, with {ind} blanks, a {prefix} character
    and a blank in front of every line.  If {cmt} is {NULL}
    or an empty string, writes nothing. Supplies a final '\n' if 
    the text is non-empty but does not end with newline. */

char *filefmt_read_comment(FILE *rd, char prefix);
  /* Reads zero or more lines from {rd} that begin with blanks followed by the {prefix}
    character; strips the leading {prefix} and the following blank (if
    present) from each line; and returns all those lines as a single
    newly allocated string, where each line is terminated by a
    newline.  If there are no comments, returns a newly allocated empty string. 
    
    In any case, skips the leading blanks in the first line after the lines
    that were parsed as comments. */

#endif
