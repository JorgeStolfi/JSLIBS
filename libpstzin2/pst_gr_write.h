/* Writing {pst_gr_t} structures. */ 
/* Last edited on 2025-03-13 06:58:11 by stolfi */

/* Created by Rafael F. V. Saracchini */

#ifndef pst_gr_write_H
#define pst_gr_write_H

#include <stdio.h>
#include <stdint.h>

#include <pst_gr.h>
#include <bool.h>

void pst_gr_write_file(FILE* wr, pst_gr_t* gr, bool_t verbose);
  /* Writes the graph {gr} to {wr}, in a format compatible 
     with {pst_gr_read}. */

void pst_gr_write_named(char *fname, pst_gr_t* gr, bool_t verbose);
  /* Writes the graph {gr} to file "{fname}", in a format compatible 
     with {pst_gr_read}. */

#endif


