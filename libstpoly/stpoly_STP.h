/* Types and tools for STP files. */
/* Last edited on 2022-10-20 05:59:46 by stolfi */

#ifndef stpoly_STP_H
#define stpoly_STP_H

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <vec.h>
#include <bool.h>
#include <r2.h>
#include <i2.h>

typedef struct stpoly_STP_edge_t { r2_t v[2]; } stpoly_STP_edge_t;
  /* A straight line segment (an /edge/) of a polygonal figure,
    as represented in an STP file. */  

typedef void stpoly_STP_edge_proc_t(int32_t line, stpoly_STP_edge_t *edge);
  /* Type of a procedure that is called by {stpoly_STP_read} to process
    each segment {*edge} read from the STP file. The {line} argument is
    a line number to be used in error messages. Note: the edge record
    {*edge} is reused. */

void stpoly_STP_read(char *fileName, bool_t binary, stpoly_STP_edge_proc_t *process_edge);
  /* Reads the STP file {fileName} and calls {process_edge(line,edge)} for each
    segment {edge} read from it.  The file formatis described by 
    {stpoly_STP_format_ascii_INFO} below.
    
    If {binary} is false, assumes ASCII STP format; in that case {line}
    is the line number in the file (counting from 1). If {binary} is
    true, assumes binary STP format; in that case, {line} is 1 for the
    header, 2 for the number of edges, and is incremented by 1 for
    each edge read. */

#define stpoly_STP_format_ascii_INFO \
  "An STP file in ASCII format begins with the keyword \"poly\" and" \
  " closes with the keyword \"endpoly\".  In between are the straight" \
  " line segments (edges) that make up the boundary of the polygonal figure, in" \
  " arbitrary order and orientation.\n" \
  "\n" \
  "  Each edge starts with the keyword \"side\" and closes" \
  " with \"endside\".  Between them are two vertices, in any order.\n" \
  "\n" \
  "  Each vertex begins with the keyword \"vertex\" followed by two" \
  " floating point coordinates, {X} and {Y}, in millimeters.  (There" \
  " is no \"endvertex\" keyword).  Each coordinate may be negative, and" \
  " may be given with up to 6 digits after the decimal point.  Only the" \
  " first 6 decimal digits are significant.  (Therefore, the max absolute" \
  " precision is 0.000001 mm, or 1 nanometer; an the max relative" \
  " precision is 0.0001%.)  The scientific format \"e{NNN}\" is not allowed.\n" \
  "\n" \
  "  Tokens (numbers and keywords) must be separated by at least one" \
  " blank or line break.  Any number of extra blanks and line breaks" \
  " are allowed before, between, and after the tokens."

#define stpoly_STP_format_binary_INFO \
  "The binary format consists of a 80 byte title, then a 4-byte" \
  " edge count {ne}, then {ne} edges, each with 16 bytes.  These" \
  " bytes consist of two coordinates ({X} and {Y}) of one endpoint, and two coordinates" \
  " of the other endpoint; each as an IEEE binary single-precision float (4 bytes)." 

/* SIZE LIMITS */

#define stpoly_n_MAX ((uint32_t)(1u << 30))
  /* Maximum number of vertices or ORIENTED edgess in a polygonal region. An index of any
    of those things will safely fit in an {int32_t} as well as in an
    {uint32_t}. */

#define stpoly_nv_MAX (stpoly_n_MAX)
  /* Max number of vertices in a polygonal region. */

#define stpoly_ne_MAX (stpoly_n_MAX / 2)
  /* Max number of UNORIENTED edges in a polygonal region. */

/* DEBUGGING */

void stpoly_STP_print_edge(FILE *wr, stpoly_STP_edge_t *edge);
  /* Prints the STP segment {edge} to {wr} for debugging purposes. */

/* VERTEX ROUNDING */

i2_t stpoly_STP_round_point(r2_t *p, double eps, bool_t even);
  /* Converts a float-valued point {p} of {\RR^2}, as read from an STP
    file, to an integer point {q} of {\RZ^2}. Namely, rounds each coordinate
    {p.c[k]} of {p} to an integer multiple {eps*q.c[k]} of {eps},
    and returns the integer vector {q}.  If {even} is true,
    the coordinates will be rounded to even integers. */

#endif
