/* nget.h -- read named parameters, in the format "NAME = VALUE". */
/* Last edited on 2013-05-28 02:12:12 by stolfilocal */

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
int nget_int(FILE *f, char *name);
int32_t nget_int32(FILE *f, char *name);
int64_t nget_int64(FILE *f, char *name);
unsigned int nget_uint(FILE *f, char *name, int base);
uint32_t nget_uint32(FILE *f, char *name, int base);
uint64_t nget_uint64(FILE *f, char *name, int base);
double nget_double(FILE *f, char *name);
char *nget_string(FILE *f, char *name);
  /* These procedures parse a {NAME = VALUE} pair. The syntax of the
    {VALUE} token is the same as in the corresponding {fget}
    procedures. */

#endif
