#ifndef obj_file_write_H
#define obj_file_write_H

/* Basic writing of OBJ format files. */ 
/* Last edited on 2024-12-22 10:49:29 by stolfi */

#define half_read_obj_H_copyright \
  "Copyright (C) 2024 Jorge Stolfi, UNICAMP.\n\n" jslibs_copyright

#include <stdio.h>
#include <stdint.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>
#include <r3.h>

#include <obj_file.h>

void obj_file_write(FILE *wr, obj_file_data_t *D, uint32_t prec);
  /* Writes into {wr} a solid model in the Wavefront OBJ format,
    with the data specified in the {D} tables.
      
    Each vertex or textpoint coordinate is written out as a decimal number
    with {prec} fraction digits. Normal vector coordinates are always written with
    {obj_file_data_prec_normal} decimal digits. */

#endif
