#ifndef gem_print_face_H
#define gem_print_face_H
/* Last edited on 2024-12-05 10:26:03 by stolfi */

#include <stdio.h>
#include <gem.h>

void gem_print_face_2d_vertices(FILE *wr, gem_ref_t root, int d);

void gem_print_face_2d_vertices_named(char *filename, gem_ref_t root, int d);
  /* Same as {gem_print_face_2d_vertices} but writes to a file named "{filename}".
    If {filename} is "-", writes to standard output. */

#endif
