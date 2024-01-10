/* pqueue - internal representation of a pqueue_t. */
/* Last edited on 2010-06-04 19:26:05 by stolfi */

#ifndef pqueue_rep_H
#define pqueue_rep_H

#include <pqueue.h>
    
typedef struct pqueue_rep_t  /* A priority queue. */
  { int order;               /* {-1} is max-first, {-1} is min-first. */
    pqueue_count_t n;        /* Number of items currently in queue */
    pqueue_item_t *itm;      /* {itm[0..n-1]} are the items in the queue. */
    double  *val;            /* {val[i]} is the value of item {itm[i]}. */
    pqueue_position_t *pos;  /* {pos[z]} is the current position of item {z} in {itm,val}. */
    pqueue_count_t nmax;     /* Alloc size of vectors {itm} and {val}. */
    pqueue_item_t zlim;      /* Alloc size of vector {pos}. */
  } pqueue_rep_t;
  /* 
    Table invariants: 
    
       (0) {itm[0..n-1]} are all the items inserted but not deleted, no dups.
       (1) {itm[i]} is in {0..zlim-1} for all {i} in {0..n-1}.
       (2) {val[i]} is the last value assigned to item {itm[i]}, times {order}.
       (3) {pos[itm[i]] == i} for all {i} in {0..n-1}.
       
    If {z} is not in the queue, then {pos[z]} is undefined.

    Heap invariants: 
    
       (4) {val[0]} is the minimum of {val[0..n-1]}.  
       (5) {val[i]} is no greater than {val[2*i+1]} and {val[2*i+2]} (when they exist).
  */

#endif
