/* nget.h -- read named parameters, in the format "NAME = VALUE". */
/* Last edited on 2024-11-16 09:17:31 by stolfi */

#ifndef nget_H
#define nget_H

/* These procedure are alternatives to {fscanf}, specialized for
  parsing named parameters, in the format "{NAME} = {VALUE}", where {NAME}
  is a given string. 
  
  All procedures below will skip spaces (SPACE, TAB, NUL), but not
  line or page breaks, before each token. They will abort with
  assertion failure if the input does not conform to the expected
  format.
  
  Created by J. Stolfi, Dec/2002, from Modula-3 version ca. 1995. */

#include <stdio.h>
#include <stdint.h>
#include <bool.h>

void nget_name_eq(FILE *f, char *name);
  /* Matches the {name} and an "=". The {name} must not contain any blanks. */

char nget_char(FILE *f, char *name);
bool_t nget_bool(FILE *f, char *name);
int32_t nget_int32(FILE *f, char *name);
int64_t nget_int64(FILE *f, char *name);
uint32_t nget_uint32(FILE *f, char *name, uint32_t base);
uint64_t nget_uint64(FILE *f, char *name, uint32_t base);
double nget_double(FILE *f, char *name);
char *nget_string(FILE *f, char *name);
  /* These procedures parse a {NAME = VALUE} pair. The syntax of the
    {VALUE} token is the same as in the corresponding {fget}
    procedures. They fail if end-of-line or end-of-file occurs before the
    pair is complete. */

#endif
