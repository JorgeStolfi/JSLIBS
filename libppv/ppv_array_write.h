/* Writing portable multidimensional arrays of samples. */
/* Last edited on 2021-06-12 00:57:54 by jstolfi */
/* Copyright © 2005 by Jorge Stolfi, from University of Campinas, Brazil. */

#ifndef ppv_array_write_H
#define ppv_array_write_H

#include <ppv_types.h>
#include <ppv_array.h>
#include <stdio.h>

void ppv_array_write_file ( FILE *wr, ppv_array_desc_t *A, bool_t plain );
  /* Outputs the array {A} to file {wr}, in a format compatible with
    {ppv_array_read_file}. 
    
    Samples are written in the plain ASCII format if {plain} is true, or
    in the binary format if {plain} is false. The memory layout
    parameters {A.base}, {A.step} and {A.bpw} are not written. See
    {ppv_array_read_FORMAT_INFO} for details. */

#endif
