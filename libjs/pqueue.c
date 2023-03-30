/* See pqueue.h */
/* Last edited on 2023-03-18 11:15:49 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>

#include <pqueue.h>
#include <pqueue_rep.h>

/* !!! Consider combining this module with a simpler list {L} of
  entries, sharing storage with the heap queue {Q} but growing at the
  opposite end. Each item would be either in {Q} or in {L}, but not in
  both. The list {L} would not have the heap property, so insert and
  delete can be {O(1)}, and append can be {O(1)} and
  position-preserving. The pair {Q,L} would be useful for heapsort and
  Dikstra-type algorithms. */

/* INTERNAL PROTOS */

void pqueue_expand(pqueue_t *Q, pqueue_item_t z);
  /* Reallocates the {Q} vectors as needed to allow the insertion
    of item {z}. The size at least doubles, so that the 
    total alloc time is O(1). */

void pqueue_bubble_up
  ( pqueue_item_t itm[],  
    pqueue_value_t val[], 
    pqueue_position_t pos[], 
    pqueue_position_t p
  );
  /* Restores the heap invariant after a decrease in {Q->val[p]}. */

void pqueue_bubble_down
  ( pqueue_item_t itm[], 
    pqueue_value_t val[], 
    pqueue_position_t pos[], 
    pqueue_position_t p, 
    pqueue_count_t n
  );
  /* Restores the heap invariant after an increase in {Q->val[p]}. */

void pqueue_check(pqueue_t *Q);
  /* Checks whether {Q} satisfies the invariants. */

/* IMPLS */

pqueue_t *pqueue_new(void)
  { pqueue_t *Q = (pqueue_t *)notnull(malloc(sizeof(pqueue_t)), "out of memory");
    Q->order = +1;
    Q->n = 0;
    Q->itm = NULL; Q->val = NULL; Q->nmax = 0;
    Q->pos = NULL; Q->zlim = 0;
    return Q;
  }
  
void pqueue_set_order(pqueue_t *Q, int32_t order)
  {
    demand(order != 0, "zero order");
    int32_t old = Q->order;
    Q->order = (order < 0 ? -1 : +1);
    if (old != Q->order)
     { /* Negate all values and re-bubble all elements: */
       int32_t n = Q->n;
       int32_t p;
       for (p = 0; p < n; p++)
         { Q->val[p] = - Q->val[p];
           pqueue_bubble_up(Q->itm, Q->val, Q->pos, p);
         }
     }
  }

void pqueue_realloc(pqueue_t *Q, pqueue_count_t nmax, pqueue_item_t zlim)
  { 
    demand(Q->n == 0, "queue not empty");
    /* If we are using too much or too little storage, free it: */
    if ((Q->nmax < nmax) || (Q->nmax/2 >= nmax))
      { /* Free {itm,val}: */
        free(Q->itm); Q->itm = NULL; 
        free(Q->val); Q->val = NULL;
        Q->nmax = 0;
      }
    if ((Q->zlim < zlim) || (Q->zlim/2 >= zlim))
      { /* Free {pos}: */
        free(Q->pos); Q->pos = NULL;
        Q->zlim = 0;
      }
    /* If we are using too little, expand it: */
    if (Q->nmax < nmax)
      { demand(nmax <= pqueue_NMAX, "to many items");
        assert(Q->nmax == 0);
        Q->itm = (pqueue_item_t *)notnull(malloc(nmax*sizeof(pqueue_item_t)), "no mem"); 
        Q->val = (pqueue_value_t *)notnull(malloc(nmax*sizeof(pqueue_value_t)), "no mem"); 
        Q->nmax = nmax;
      }
    if (Q->zlim < zlim)
      { demand(zlim <= pqueue_ITEM_MAX+1, "item too big");
        assert(Q->zlim == 0);
        Q->pos = (pqueue_position_t *)notnull(malloc(zlim*sizeof(pqueue_position_t)), "no mem");
        Q->zlim = zlim;
      }
  }
  
void pqueue_free(pqueue_t *Q)
  { pqueue_reset(Q); pqueue_realloc(Q,0,0); free(Q); }

/* QUERIES (WITHOUT SIDE EFFECTS) */

pqueue_count_t pqueue_count(pqueue_t *Q)
  { return Q->n; }

pqueue_item_t pqueue_head(pqueue_t *Q)
  { demand(Q->n > 0, "queue is empty");
    return Q->itm[0];
  }

pqueue_position_t pqueue_position(pqueue_t *Q, pqueue_item_t z)
  { /* if (z < 0) { return Q->n; } */
    if (z >= Q->zlim) { return Q->n; }
    pqueue_position_t p = Q->pos[z];
    if ((p > Q->n) || (Q->itm[p] != z)) { return Q->n; }
    return p;
  }

bool_t pqueue_has(pqueue_t *Q, pqueue_item_t z)
  { pqueue_position_t p = pqueue_position(Q, z);
    return (p < Q->n);
  }

pqueue_value_t pqueue_value(pqueue_t *Q, pqueue_item_t z)
  { pqueue_position_t p = pqueue_position(Q, z);
    demand(p < Q->n, "item not in queue");
    return Q->order * Q->val[p];
  }

pqueue_item_t pqueue_item(pqueue_t *Q, pqueue_position_t p)
  { demand(p < Q->n, "no such position in queue");
    return Q->itm[p];
  }

int32_t pqueue_order(pqueue_t *Q)
  { return Q->order; }

/* MODIFYING */

void pqueue_expand(pqueue_t *Q, pqueue_item_t z)
  { if ((Q->itm == NULL) || (Q->n >= Q->nmax))
      { /* Reallocate {Q->itm,Q->val}: */
        pqueue_count_t nmax_new = 2*Q->nmax + 15;
        Q->itm = realloc(Q->itm, nmax_new*sizeof(pqueue_item_t));
        affirm(Q->itm != NULL, "out of mem");
        Q->val = realloc(Q->val, nmax_new*sizeof(pqueue_value_t));
        affirm(Q->val != NULL, "out of mem");
        Q->nmax = nmax_new;
      }
    if ((Q->pos == NULL) || (z >= Q->zlim))
      { /* Reallocate {Q->pos}: */
        pqueue_item_t zlim_new = Q->zlim + z + 15;
        Q->pos = realloc(Q->pos, zlim_new*sizeof(pqueue_position_t));
        affirm(Q->pos != NULL, "out of mem");
        Q->zlim = zlim_new;
      }
  }

void pqueue_bubble_up
  ( pqueue_item_t itm[],  
    pqueue_value_t val[], 
    pqueue_position_t pos[], 
    pqueue_position_t p
  )
  { pqueue_item_t tel = itm[p];
    pqueue_value_t tva = val[p];
    pqueue_position_t k = p; /* Index of bubble. */
    pqueue_position_t j;     /* {j} is the parent of {k}. */
    /* fprintf(stderr, "bubble up %d:", p); */
    while ((k > 0) && (tva < val[(j = (k-1)/2)])) 
      { itm[k] = itm[j]; val[k] = val[j]; pos[itm[k]] = k; k = j; /* fprintf(stderr, " %d", k); */ }
    if (k != p) { itm[k] = tel; val[k] = tva; pos[tel] = k; }
    /* fprintf(stderr, "\n"); */
  }

void pqueue_bubble_down
  ( pqueue_item_t itm[], 
    pqueue_value_t val[], 
    pqueue_position_t pos[], 
    pqueue_position_t p, 
    pqueue_count_t n
  )
  { pqueue_item_t tel = itm[p];
    pqueue_value_t tva = val[p];
    pqueue_position_t k = p; /* Index of bubble. */
    pqueue_position_t ja;    /* {ja} is the first child of {k}. */
    /* fprintf(stderr, "bubble down %d:", p); */
    while ((ja = 2*k + 1) < n)
      { /* Find smallest child {itm[j]} of {itm[k]}: */
        pqueue_position_t jb = ja + 1; /*  {itm[jb]} is the second child of {itm[k]}. */
        pqueue_position_t j = ((jb >= n) || (val[ja] <= val[jb]) ? ja : jb);
        if (tva <= val[j]) { break; }
        /* Promote smallest child into hole: */
        itm[k] = itm[j]; val[k] = val[j]; pos[itm[k]] = k;
        k = j; /* fprintf(stderr, " %d", k); */
      }
    if (k != p) { itm[k] = tel; val[k] = tva; pos[tel] = k; }
    /* fprintf(stderr, "\n"); */
  }

void pqueue_insert(pqueue_t *Q, pqueue_item_t z, pqueue_value_t v)
  { pqueue_position_t p = pqueue_position(Q, z);
    demand(p == Q->n, "item already in queue");
    pqueue_expand(Q, z);
    Q->itm[p] = z;
    Q->val[p] = Q->order * v;
    Q->pos[z] = p;
    Q->n++;
    pqueue_bubble_up(Q->itm, Q->val, Q->pos, p);
    /* pqueue_check(Q); */
  }

void pqueue_delete(pqueue_t *Q, pqueue_item_t z)
  { pqueue_position_t p = pqueue_position(Q, z);
    demand(p < Q->n, "item not in queue");
    
    pqueue_count_t n = Q->n;
    pqueue_item_t *itm = Q->itm;
    pqueue_value_t *val = Q->val;
    pqueue_position_t *pos = Q->pos;
    
    /* Now slot {itm[p],val[p]} is vacant. */
    /* Promote children into vacancy {itm[p],val[p]} until {p} reaches the fringe: */
    pqueue_position_t ja; /* {itm[ja]} is the first child of {itm[p]}. */
    while ((ja = 2*p + 1) < n)
      { /* Find smallest child {itm[j]} of {itm[p]}: */
        pqueue_position_t jb = ja + 1; /*  {itm[jb]} is the second child of {itm[p]}. */
        pqueue_position_t j = ((jb >= n) || (val[ja] <= val[jb]) ? ja : jb);
        /* Promote smallest child into hole: */
        itm[p] = itm[j]; val[p] = val[j]; pos[itm[p]] = p;
        p = j;
      }
    /* One less element in queue: */
    n--;
    if (p < n)
      { /* The vacancy {itm[p]} did not end up at {itm[n]}, so fill it with {itm[n]}: */        
        pqueue_item_t tel = itm[n];
        pqueue_value_t tva = val[n];
        pqueue_position_t j; 
        /* Bubble it up to the proper place: */
        while ((p > 0) && (tva < val[(j = (p-1)/2)]))
          { itm[p] = itm[j]; val[p] = val[j]; pos[itm[p]] = p;
            p = j;
          }
        itm[p] = tel; val[p] = tva; pos[itm[p]] = p;
      }
    Q->n = n;
    /* pqueue_check(Q); */
  }

void pqueue_set_value(pqueue_t *Q, pqueue_item_t z, pqueue_value_t v)
  { pqueue_position_t p = pqueue_position(Q, z);
    demand(p < Q->n, "item not in queue");
    v *= Q->order;
    pqueue_value_t old_v = Q->val[p];
    Q->val[p] = v;
    if (v < old_v)
      { pqueue_bubble_up(Q->itm, Q->val, Q->pos, p); }
    else if (v > old_v)
      { pqueue_bubble_down(Q->itm, Q->val, Q->pos, p, Q->n); }
    /* pqueue_check(Q); */
  }

void pqueue_set_all_values(pqueue_t *Q, pqueue_value_function_t *f)
  { pqueue_position_t p;
    for (p = 0; p < Q->n; p++)
      { pqueue_item_t z = Q->itm[p];
        pqueue_value_t v = Q->order*Q->val[p]; 
        Q->val[p]= Q->order*f(z, v);
        pqueue_bubble_up(Q->itm, Q->val, Q->pos, p);
      }
  }

void pqueue_reset(pqueue_t *Q)
  { Q->n = 0; }

void pqueue_check(pqueue_t *Q)
  {
    /* Validate {Q->n}: */
    /* assert(Q->n >= 0); */
    assert(Q->n <= Q->nmax); 
    pqueue_position_t p;
    for (p = 0; p < Q->n; p++)
      { /* Get item {z} in slot {p}: */
        pqueue_item_t z = Q->itm[p];
        /* Items must be in range {0..Q->zlim-1}: */
        /* assert(z >= 0);  */
        assert(z < Q->zlim); 
        /* Get claimed position of item {z}: */
        pqueue_position_t k = Q->pos[z];
        /* Must be this position: */
        assert(k == p);
        /* Check heap invatiant: */
        if (p > 0)
          { int32_t p = (k-1)/2; /* Parent of slot {k}. */
            assert(Q->order*pqueue_dblcmp(Q->val[p], Q->val[p]) >= 0);
          }
      }
  }
  
int32_t pqueue_dblcmp(pqueue_value_t x, pqueue_value_t y) 
  { return (x < y ? -1 : ( x > y ? +1 : 0)); }
