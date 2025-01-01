#ifndef cpk_debug_H
#define cpk_debug_H

/* Debugging routines for circle packing. */
/* Last edited on 2024-12-31 15:53:16 by stolfi */ 

#include <stdint.h>
#include <stdio.h>

#include <r2.h>
#include <pqueue.h>

#include <cpk_basic.h>

/* 
  DEBUGGING PRINTOUTS
  
  In all these procedures, if {nMax} is not negative, 
  prints only the first and last {nMax} vertices. */

void cpk_print_vertices(FILE *wr, char *title, r2_vec_t *V, double_vec_t *W, uint32_t nMax);
  /* Prints vertices {V.e[0..]} and (if {W} is not NULL) their
    weights {W[0..]} to file {wr}. */

void cpk_print_edges(FILE *wr, char *title, ui2_vec_t *E, r2_vec_t *V, uint32_t nMax);
  /* Prints edges {E.e[0..]} to file {wr}, showing the indices and
    (if {V} is not null) the coordinates of the two endpoints. */

void cpk_print_vertex_set(FILE *wr, char *title, uint32_t nJ, uint32_t *J, r2_vec_t *V, uint32_t nMax);
  /* Prints the vertex list {J[0..nJ-1]} to {wr}.  Iv {V} is not null, prints also 
   their coordinates. */

/* OPTIMIZATION TRACING */

void TRACE_SOL(char *title, uint32_t nS, double WS);

void TRACE_Q(char *title, pqueue_t *Q);

/* PROGRAM TRACING */

#define cpk_trace(ARGS) \
  do { \
    fprintf(stderr, "\n[%s at %s:%d]", __FUNCTION__, __FILE__, __LINE__); \
    fprintf(stderr, " " ARGS);  fprintf(stderr, "\n"); \
  } while (0)

#define cpk_trace_entry(ARGS) \
  do { \
    fprintf(stderr, "\n[++ enter %s at %s:%d]", __FUNCTION__, __FILE__, __LINE__); \
    fprintf(stderr, " " ARGS);  fprintf(stderr, "\n"); \
  } while (0)

#define cpk_trace_exit(ARGS) \
  do { \
    fprintf(stderr, "\n[-- exit  %s at %s:%d]", __FUNCTION__, __FILE__, __LINE__); \
    fprintf(stderr, " " ARGS);  fprintf(stderr, "\n"); \
  } while (0)

#endif
