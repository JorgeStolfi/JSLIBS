/* Reading {psgr_t} structures. */ 
/* Last edited on 2025-04-27 02:58:28 by stolfi */

/* Created by Rafael F. V. Saracchini */

#ifndef psgr_read_H
#define psgr_read_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>

#include <psgr_types.h>
#include <psgr.h>

psgr_t* psgr_read_file(FILE* rd, bool_t verbose);
  /* Reads from {rd} a descripton of a graph {gr}. */

psgr_t* psgr_read_named(char *fname, bool_t verbose);
  /* Reads from file "{fname}" a descripton of a graph {gr}. */

#endif


