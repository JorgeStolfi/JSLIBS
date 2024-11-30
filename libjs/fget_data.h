/* fget_data.h -- flexible parsing of multicolum data files. */
/* Last edited on 2024-11-23 07:08:34 by stolfi */

#ifndef fget_data_H
#define fget_data_H

/* This module lets parse selected columns of data files with multiple 
  space-separated columns, some numeric, some non-numeric, in flexible order. */

#include <stdio.h>
#include <stdint.h>

#include <bool.h>

typedef enum 
  { fget_data_type_NOT,  /* Skip this data field. */
    fget_data_type_ALF,  /* Parse this field as a string. */
    fget_data_type_NUM   /* Parse this field as a number. */
  } fget_data_type_t;
  /* Data type for {fget_data_fields} (quod videt). */
#define fget_data_type_FIRST fget_data_type_NOT
#define fget_data_type_LAST fget_data_type_NUM

bool_t fget_data_fields
  ( FILE *rd,
    char cmtc,
    uint32_t nf,
    fget_data_type_t type[],
    char* alf[],
    double num[]
  );
  /* Parses or skips the next {nf} numeric or string data fields 
    on the same line.
     
    Those fields are implicitly indexed {0.nf-1}, and their parsing and
    handling is determined by the table {type[0..nf-1]}. Accordin to {type[kf]},
    field {kf} may be skipped, parsed as a string and stored into {alf[kf]}, or
    parsed as a {double} floating-point number and stored in {num[kf]}.
 
    Specifically, this procedure skips any formatting chars (spaces,
    line breaks, page breaks) and comments. If it hits end-of-file at
    this stage, it returns {FALSE}. Otherwise it tries to parse {nf}
    fields, all on the current line. If it succeeds, it returns {TRUE}.
    
    Comments are assumed to start with the character {cmtc} (e.g. '#') 
    and extend to the end of the same line.
    
    The procedure skips spaces, but not line breaks or comments, before
    parsing each field after the first one. If {type[kf]} is
    {fget_data_type_NUM}, field {kf} is parsed with {fget_double} and
    stored in {num[kf]}. Otherwise field {kf} is parsed as a sequence of
    one or more arbitrary characters up to the next formatting
    character, end-of-file, or comment -- like {fget_string}, but
    excluding the {cmtc} character. Then, if {type[kf]} is
    {fget_data_type_ALF}, these character are turned into a
    newly-allocated {char*} string and stored in {alf[kf]}. If
    {type[kf]} is {fget_data_type_NOT}, those characters are discarded.
    Any other value of {type[kf]} is a fatal error.
    
    If {alf[kf]} is not set as described above, it is set to {NULL}. If
    {num[kf]} is not set as described above, it is set to {NAN}.
    
    The procedure fails with error if the data line is incomplete;
    namely, if it hits end-of-line, end-of-file, or one of the {cmtc}
    character before parsing any field after the first one.
    
    For good measure, any two consecutive fields should be separated by
    one or more spaces, and the last field should be followed by a
    space, a {cmtc} character, or end-of-line. (However, the spaces are
    not required after a numeric field if the next character cannot be
    parsed as part of the number.)
    
    The {cmtc} character had better not be any space character, EOF, or
    any character that may be part of a number ('+', '-', '.', 'e', or
    digit '0'..'9'). If {cmtc} is '\n', comments will not be allowed,
    whether in blank lines or at the end of data lines. */

void fget_data_set_field_type
  ( uint32_t kf,
    fget_data_type_t tkf,
    bool_t rep_ok,
    uint32_t nf,
    fget_data_type_t type[]
  );
  /* Sets {type[kf]} to {tkf}, thus detemining the type of field {kf}
    in subsequent calls of {fget_data_fields}.
    
    The  index {kf} should be in {0..nf-1}.
    
    If {rep_ok} is false, the procedure fails with error if {type[kf]}
    has already been set to a non-zero value. If {rep_ok} is true, such
    a repeated assignment is accepted, provided {tkf} is the same as
    {type[kf]}. These rules help detecting improper field assignments by
    the caller, such as trying to use the same column twice. */

#endif
