/* Reading {pst_gr_t} structures. */ 
/* Last edited on 2025-03-13 06:57:45 by stolfi */

/* Created by Rafael F. V. Saracchini */

#ifndef pst_gr_read_H
#define pst_gr_read_H

#include <stdio.h>
#include <stdint.h>

#include <bool.h>

#include <pst_gr.h>

pst_gr_t* pst_gr_read_file(FILE* rd, bool_t verbose);
  /* Reads from {rd} a descripton of a graph {gr}. */

pst_gr_t* pst_gr_read_named(char *fname, bool_t verbose);
  /* Reads from file "{fname}" a descripton of a graph {gr}. */

#endif


