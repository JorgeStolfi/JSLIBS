/* fget.h -- reading items from text files. */
/* Last edited on 2024-11-15 19:12:22 by stolfi */

#ifndef fget_H
#define fget_H

/* This module provides various functions to read data items from text
  files, skipping spaces and '#'-comments, etc. They are alternatives to
  {fscanf} that provide better control on line breaks, and normally
  abort the program in case of syntax error or unexpected end-of-file.
  
  These procedures assume that every line in a file ends with a single
  ASCII LF character (line feed, '\n', '\012', ^J), referred to as
  "end-of-line " in the comments below. In particular, every non-empty
  file ends with end-of-line.
  
  These procedures will NOT recognize ASCII CR (carriage-return, '\015',
  ^M) as a line terminator or separator, but will instead treat it as a
  non-space character (one of the "formatting characters"), even if it
  is immediately followed by LF.
  
  See the end of this file for changes made on 2023-10-15.    
  
  Created on Dec/2002 by J. Stolfi, from Modula-3 version ca. 1995,
  inspired on Lex.m3 by L. Cardelli. */

#include <stdio.h>
#include <stdint.h>

#include <bool.h>

/* SPACES */
    
#define fget_spaces " \240\011"
  /* The space characters: ASCII SPACE='\040',
    TAB='\011', and ISO-LATIN-1 NBSP='\240'. */

bool_t fget_is_space(char c);
  /* Returns TRUE iff {c} is a space char. Note
    that {c} cannot be {EOF}, since {(char)EOF} coincides with '\240',
    non-breaking space.. */

void fget_skip_spaces(FILE *rd); 
  /* Skips space characters until the first non-space character, or end-of-file.
    Will NOT skip over any other characters, including end-of-line and 
    the other non-space formatting characters. */
    
/* FORMATTING CHARACTERS */
    
#define fget_formatting_chars " \240\011\012\013\014\015"
  /* The formatting characters: the spaces, plus end-of-line (ASCII LF='\012')
    vertical tab (VT='\013'), carriage-return ('CR='\015'),
    and page break (FF='\014'). */

bool_t fget_is_formatting_char(char c);
  /* Returns TRUE iff {c} is a formatting character. Note that {c} cannot be {EOF},
    since {EOF} is not in the {char} range {-128..+127}. */

void fget_skip_formatting_chars(FILE *rd); 
  /* Skips all blank `formatting' characters until the first character that is not
    in that list, or end-of-file. */

/* DAAT ITEMS */

/* The following procedures will skip spaces (SPACEs, NBSPs, TABs and NULs)
  before the desired input.  They will NOT skip line-breaks or page-breaks. */

char fget_char(FILE *rd);
  /* Skips spaces, then reads a single character --- which must
    not be a space or formatting char (such as end-of-line , '\n'). Fails
    if a formatting char or end-of-file occurs before the character.
    
    Note that a result of {-1 = '\377'} means the byte 0xff
    not {EOF}. */

bool_t fget_bool(FILE *rd);
  /* The procedure {fget_bool_t} will skip spaces, then read a single
    character: 'T','t','Y','y','1' for TRUE, 'F','f','N','n','0' for FALSE. Anything
    else is an error.  Note that, if the file contains "TRUE" or
    "FALSE", only the first letter will be consumed. */

char *fget_string(FILE *rd);
  /* Skips spaces. Fails if hits end-fo-file or a formatting character.
    Otherwise reads one or more characters, until end-of-file or the
    a formatting character (space, line break, or page break), which
    is not consumed.  The result is returned as a newly allocated,
    zero-terminated string (which has at least one character). */
    
int32_t fget_int32(FILE *rd);
int64_t fget_int64(FILE *rd);
  /* Skips spaces, reads the longest possible non-empty string that
    looks like a decimal integer (with optional '+' or '-' sign), and
    tries to convert it to a signed integer value of the specified
    size. A failure in any of these steps is a fatal error. */
    
uint32_t fget_uint32(FILE *rd, uint32_t base);
uint64_t fget_uint64(FILE *rd, uint32_t base);
  /* Skips spaces, reads the longest possible non-empty string that
    looks like an unsigned integer, and tries to convert it to an 
    unsigned integer value of the specified size. A failure in any 
    of these steps is a fatal error.
    
    The digit string found is interpreted relative to the given
    {base}, which should be in {2..36}. If {base} is 10 or less, the
    function consumes any string of digits in '0'..'9'; if {base} is
    11 or more, it consumes any string of digits or letters 'a'..'z'
    (or 'A'..'Z') which represent values in 10..35. In either case,
    the procedure fails if any digit has a value {base} or more. */

uint64_t fget_bin64(FILE *rd);
  /* Skips spaces, reads the longest possible non-empty string of
    digits, ans tries to parseit as the binary representation of an
    {uint64_t} value. Fails if there are no digits, if it runs into
    a digit greater than 1, or if the value is greater than {2^64-1}. */

double fget_double(FILE *rd);
  /* Skips spaces, reads the longest possible non-empty string that
    looks like a floating-point number (including optional exponent),
    and tries to convert it to a {double} value. A failure in any of
    these steps is a fatal error. */

/* STRING AND CHAR PARSING */

bool_t fget_test_eof(FILE *rd);
  /* Returns {TRUE} iff there are no more characters to be read from {rd}.
    Otherwise leaves {rd} effectively unchanged.  */
 
void fget_eol(FILE *rd);
  /* Skips any spaces and requires the next non-space character to be a
    end-of-line  (ASCII LF, '\n', '\012', ^J), which it consumes. It is a
    fatal error if the next non-space is not a end-of-line, or end-of-file
    is found. Equivalent to {fget_skip_spaces(rd);
    fget_match(rd,"\n")}. */

void fget_skip_to_eol(FILE *rd);
  /* Reads and discards arbitrary ] characters from {rd} until an
    end-of-line  or end-of-file. The end-of-line  character is
    consumed. */
   
char *fget_line(FILE *rd);
  /* Collects all characters from the current point up to but 
    not including the end of the line or end of file into as a newly allocated,
    zero-terminated string, which will be the returned result.  Then,
    if the next character is end-of-line , consumes it. The end-of-file,
    if any, is not consumed. */

bool_t fget_test_comment_or_eol(FILE *rd, char cmtc, char **text_P);
  /* Skips spaces, and then checks whether the next character is either
    end-of-line, or a comment that starts with the character {cmtc} and
    extends to the end of the line. In either case, returns {TRUE} and
    consumes everything up to and including the end-of-line. If the
    next character after skipping the spaces is anything other than {cmtc} or
    end-of-line, leaves that character there and returns {FALSE}. Fails with
    error if {cmtc} occurs but end-of-file is reached before finding the 
    end-of-line.
    
    If the procedure returns {TRUE} (i.e., finds end-of-line  or
    {cmtc}), and {text_P} is not {NULL}, the procedure sets {*text_P} to the
    comment text parsed, if any. More precisely, if a {cmtc} character
    was found, {*text_P} is set to a newly-allocated, zero-terminated
    string containing the characters between {cmtc} and the end-of-line ,
    excluding both. This string may be empty but will be not {NULL}, and
    will always contain the terminating '\000' byte. If end-of-line was
    found instead of {cmtc}, the procedure sets {*text_P} to {NULL}. */

void fget_comment_or_eol(FILE *rd, char cmtc, char **text_P);
  /* Like {fget_eol}, but allows for a comment that begins 
    with {cmtc} and extends to the end of the line.  
    Namely, the same as {fget_test_comment_or_eol(rd,cmtc,text_P)},
    but fails with an error instead of returning {FALSE}. */

char *fget_to_delims(FILE *rd, char delim, char *delims);
  /* Reads zero or more characters, until a terminating character -- which
    is either end-of-line , or the character {delim}, or one of the characters in {delims},
    whichever comes first. 
    
    If {delims} is {NULL}, uses only {delim}. In particular, if 
    {delim} is '\n' and {delims} is NULL, it just reads characters until
    end-of-line .
    
    The terminating character is NOT consumed in any case.
    The characters read, except the terminating one, are 
    returned as a newly allocated zero-terminated string,
    (which may be empty). 
    
    The procedure fails with error if end-of-file occurs before the
    terminating character is found. */

bool_t fget_test_char(FILE *rd, char c);
  /* If {rd} is exhausted, returns {FALSE}.
    Otherwise, checks whether the next character is {c} (which may be anything,
    including space or formatting char). If it is, consumes that
    character and returns TRUE.  If it is some other character, (even if it is a
    space or formatting char), returns FALSE and leaves the character
    there. 
    
    Note that {c = -1 = '\377'} will match only the byte 0xff, not {EOF}. */

bool_t fget_skip_and_test_char(FILE *rd, char c);
  /* Equivalent to {fget_skip_spaces(rd); fget_test_char(rd, c)}. */

void fget_match(FILE *rd, char *t);
  /* Requires the string {t} to be the next thing on {rd}, and
    consumes that string. */
    
void fget_skip_spaces_and_match(FILE *rd, char *t);
  /* Equivalent to {fget_skip_spaces(rd); fget_match(rd, t)}. */
 
/* CHANGED on 2023-10-15:

    * The ASCII NUL='\000' was removed from  {fget_formatting_chars} 
      and {fget_spaces}.
      
    * {fget_delims} no longer skips spaces on entry. It now takes both a
      single character {delim} and character string {delims}; the
      terminating char is any of those characters OR end-of-line
      (LF='\012'). It no longer stops at formatting characters other
      than end-of-line, unless those are included in {delim} and
      {delims}. In any case the terminating character will NOT consumed.
      It fails with error if runs into EOF before the terminating
      character, even if right away on entry.
      
    * {fget_test_comment_or_eol} will return FALSE if, after skipping
      spaces, runs into EOF or any character other than {cmtc} and
      end-of-line. It fails with error if, after skipping spaces, finds
      {cmtc} but then runs into EOF before any end-of-line. When it
      returns {TRUE}, it optionally returns also the comment text.
      
    * {fget_comment_or_eol} is like {fget_test_comment_or_eol} but fails
      with error message instead of returning {FALSE}. Namely, if fails
      if, after skipping spaces, does not find either end-of-line or a
      complete comment ({cmtc}, zero or more chars, and end-of-line). */   

/* DEBUGGING */

void fget_show_next(FILE *wr, char *pref, FILE *rd, char *suff);
  /* Prints to {wr} the next character from {rd}, legibly, without consuming it,
    surrounded by {pref} and {suff}. */ 

#endif
