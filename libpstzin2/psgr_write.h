/* Writing {psgr_t} structures. */ 
/* Last edited on 2025-04-27 02:58:17 by stolfi */

/* Created by Rafael F. V. Saracchini */

#ifndef psgr_write_H
#define psgr_write_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>

#include <psgr_types.h>
#include <psgr.h>

void psgr_write_file(FILE* wr, psgr_t* gr, bool_t verbose);
  /* Writes the graph {gr} to {wr}, in a format compatible 
     with {psgr_read}. */

void psgr_write_named(char *fname, psgr_t* gr, bool_t verbose);
  /* Writes the graph {gr} to file "{fname}", in a format compatible 
     with {psgr_read}. */

#endif


