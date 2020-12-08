#ifndef gem_print_graph_H
#define gem_print_graph_H
/* Last edited on 2014-07-06 20:57:08 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <gem.h>

void gem_print_graph(FILE *wr, gem_ref_t root, int d);
  /* Write a {d}-dimensional gem rooted at {root} to the file {wr}.
    The file must be open for writing, and is not closed. */

void gem_print_graph_named(char *filename, gem_ref_t p, int d);
  /* Write a {d}-dimensional gem rooted at {root} to a file named {filename}.
    If {filename} is "-", writes to the standard output. */

#endif
