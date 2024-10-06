/* Writing a {tosl.h} mesh as an OBJ file.  */
/* Last edited on 2024-10-06 16:46:54 by stolfi  */

#ifndef tosl_mesh_obj_write_H
#define tosl_mesh_obj_write_H

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <tosl.h>
#include <tosl_mesh.h>

void tosl_mesh_obj_write(FILE *wr, tosl_mesh_t *mesh);
  /* Writes to {wr} the given {mesh} in OBJ format.  In the file, 
  vertices are numbered from {1..mesh.NV} instead of  {0..mesh.NV-1}. */

#endif
