/* Last edited on 2024-11-16 12:56:30 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <bool.h>
#include <affirm.h>
#include <jsrandom.h>

#include <pqueue.h>

#define ER(...) fprintf(stderr, __VA_ARGS__)

#define NMAX_MAX (1024*1024*1024)
  /* Max num of elements that can be requested. */

uint32_t throw_index(uint32_t N);
  /* Generates a random index in the range {0..N-1}}. */

pqueue_value_t throw_distinct_value(uint32_t N, pqueue_value_t ival[]);
  /* Generates a random value in the range [0 _ 1], distinct from {ival[0..N-1]}. */

pqueue_item_t throw_new_item(uint32_t N, pqueue_item_t zlim, pqueue_item_t item[]);
  /* Generates a random item in the range [0 _ zlim-1], distinct from {item[0..N-1]}. */

void check_queue
  ( pqueue_t *Q, 
    uint32_t N, 
    pqueue_item_t zlim,
    pqueue_item_t item[], 
    pqueue_value_t ival[], 
    int32_t order, 
    bool_t verbose
  );
  /* Checks the consistency of queue {Q} (as seen through its query 
    operations) with the brute-force queue {N,item,ival,order}
    for items in {0..zlim-1}. */

void perform_operation
  ( uint32_t iop,
    pqueue_t *Q, 
    uint32_t *NP, 
    uint32_t NMAX,
    pqueue_item_t zlim,
    pqueue_item_t item[], 
    pqueue_value_t ival[], 
    int32_t *orderP, 
    bool_t verbose
  );
  /* Performs one of the following operations

      {pqueue_set_order}
      {pqueue_insert}
      {pqueue_delete}
      {pqueue_set_value}
      {pqueue_set_all_values}

    on the queue {Q} and on the brute-force queue {N,item,ival,order}
    with items in {0..zlim-1}. 
  */

int32_t main(int32_t argc, char **argv);

int32_t main(int32_t argc, char **argv)
  { int64_t arg1 = atol(argv[1]);
    demand((arg1 > 0) && (arg1 < NMAX_MAX), "invalid {NMAX}");
    uint32_t NMAX = (uint32_t)arg1;
  
    pqueue_item_t zlim = 2*NMAX + 15;
    bool_t verbose = (NMAX <= 200);
    
    srandom(4615 + 32*NMAX);
  
    /* The brute-force queue, sorted by increasing value: */
    uint32_t N = 0; /* Number of items in queue. */
    pqueue_item_t *item = talloc(NMAX, pqueue_item_t); /* Items. */
    pqueue_value_t *ival = talloc(NMAX, pqueue_value_t); /* Item values. */
    int32_t order = +1;  /* The nominal queue order. */
    
    /* The library queue: */
    pqueue_t *Q = pqueue_new();
    check_queue(Q, N, zlim, item, ival, order, verbose);
    
    pqueue_realloc(Q, NMAX/3, zlim);
    check_queue(Q, N, zlim, item, ival, order, verbose);
    
    uint32_t NTEST = 3*NMAX;
    uint32_t iop;
    for (iop = 0; iop < NTEST; iop++)
      { perform_operation(iop, Q, &N, NMAX, zlim, item, ival, &order, verbose);
        check_queue(Q, N, zlim, item, ival, order, verbose);
      }
      
    pqueue_reset(Q);
    N = 0;
    check_queue(Q, N, zlim, item, ival, order, verbose);
    
    pqueue_free(Q);
    
    ER("done.\n");
      
    return 0;
  }

void perform_operation
  ( uint32_t iop,
    pqueue_t *Q, 
    uint32_t *NP, 
    uint32_t NMAX,
    pqueue_item_t zlim,
    pqueue_item_t item[], 
    pqueue_value_t ival[], 
    int32_t *orderP, 
    bool_t verbose
  )
  {
    uint32_t N = (*NP);
    int32_t order = (*orderP);
    
    double coin;
    if (iop <= 20)
      { coin = ((double)iop)/20.0; }
    else
      { coin = drandom(); }

    /* Probabilities of various operations (sum must be ~1): */
    double P_insert = 0.30;
    double P_delete = 0.30;
    double P_set_value = 0.30;
    double P_set_all_values = 0.05;
    double P_set_order = 0.05;
    
    /* Test {pqueue_insert}: */
    if ((N < NMAX) && (coin < P_insert))
      { pqueue_item_t z = throw_new_item(N, zlim, item);
        pqueue_value_t v = throw_distinct_value(N, ival);
        if (verbose) { ER("%06d inserting item %08u with value %18.16f\n", iop, z, v); }
        pqueue_insert(Q, z, v);
        uint32_t i = N;
        while ((i > 0) && (ival[i-1] > v)) { item[i] = item[i-1]; ival[i] = ival[i-1]; i--; }
        N++;
        item[i] = z;
        ival[i] = v;
        if (i > 0) { assert(ival[i-1] < ival[i]); }
        if (i < N-1) { assert(ival[i] < ival[i+1]); }
        coin = +INFINITY;
      }
    else
      { coin = fmax(0, coin - P_insert); }

    /* Test {pqueue_delete}: */
    if ((N > 0) && (coin < P_delete))
      { uint32_t i = throw_index(N);
        pqueue_item_t z = item[i];
        if (verbose) { ER("%06d deleting item %08u\n", iop, z); }
        pqueue_delete(Q, z);
        N--;
        while (i < N) { item[i] = item[i+1]; ival[i] = ival[i+1]; i++; }
        coin = +INFINITY;
      }
    else
      { coin = fmax(0, coin - P_delete); }

    /* Test {pqueue_set_value}: */
    if ((N > 0) && (coin < P_set_value))
      { uint32_t i = throw_index(N);
        pqueue_item_t z = item[i];
        pqueue_value_t old_v = ival[i];
        pqueue_value_t new_v = (coin < 0.1*P_set_value ? ival[i] : throw_distinct_value(N, ival));
        if (verbose) { ER("%06d changing value of %08u from %18.16f to %18.16f\n", iop, z, old_v, new_v); }
        pqueue_set_value(Q, z, new_v);
        while ((i > 0) && (new_v < ival[i-1])) { item[i] = item[i-1]; ival[i] = ival[i-1]; i--; }
        while ((i < N-1) && (new_v > ival[i+1])) { item[i] = item[i+1]; ival[i] = ival[i+1]; i++; }
        item[i] = z;
        ival[i] = new_v;
        if (i > 0) { assert(ival[i-1] < ival[i]); }
        if (i < N-1) { assert(ival[i] < ival[i+1]); }
        coin = +INFINITY;
      }
    else
      { coin = fmax(0, coin - P_set_value); }
        
    /* Test {pqueue_set_all_values}: */
    auto pqueue_value_t rbump(pqueue_item_t z, pqueue_value_t v);
      /* A one-to-one function from {pqueue_value_t} to {pqueue_value_t}. */
      
    if (coin < P_set_all_values)
      { if (verbose) { ER("%06d modifying all values through {rbump}\n", iop); }
        pqueue_set_all_values(Q, rbump);
        for (int32_t i = 0; i < N; i++) { ival[i] = rbump(item[i], ival[i]); }
        /* Insertion re-sort: */
        for (int32_t i = 0; i < N; i++) 
          { pqueue_item_t z = item[i];
            pqueue_value_t v = ival[i];
            uint32_t j = i;
            while ((j > 0) && (ival[j-1] > v))
              { item[j] = item[j-1]; ival[j] = ival[j-1]; j--; }
            item[j] = z;
            ival[j] = v;
          }
        coin = +INFINITY;
      }
    else
      { coin = fmax(0, coin - P_set_all_values); }

    /* Test {pqueue_set_order}: */
    if (coin < P_set_order)
      { int32_t old_order = order;
        int32_t new_order = 2*int32_abrandom(0,1)-1;
        if (verbose) { ER("%06d changing the order from %+2d to %+2d\n", iop, old_order, new_order); }
        pqueue_set_order(Q, new_order);
        order = new_order;
        coin = +INFINITY;
      }
    else
      { coin = fmax(0, coin - P_set_order); }
    
    (*NP) = N;
    (*orderP) = order;
    
    return;
    
    /* Local implementations: */
    
    pqueue_value_t rbump(pqueue_item_t z, pqueue_value_t v)
      { 
        v = v*256;
        int32_t iv0 = (int32_t)floor(v);
        v = v - (double)iv0;
        
        v = v*127;
        int32_t iv1 = (int32_t)floor(v);
        v = v - (double)iv1;
        
        v = (iv1 + (iv0 + v)/256.0)/127.0;
        return v;
      }
  }

void check_queue
  ( pqueue_t *Q, 
    uint32_t N, 
    pqueue_item_t zlim,
    pqueue_item_t item[], 
    pqueue_value_t ival[], 
    int32_t order, 
    bool_t verbose
  )
  {
    pqueue_position_t Q_N = pqueue_count(Q);
    if (verbose || (Q_N != N)) { ER("  pqueue_count(Q) = %u\n", Q_N); }
    if (Q_N != N) { ER("  ** should be %u\n", N); assert(FALSE); }
    
    int32_t Q_order = pqueue_order(Q);
    if (verbose || (Q_order != order)) { ER("  pqueue_order(Q) = %+2d\n", Q_order); }
    if (Q_order != order) { ER("  ** should be %+2d\n", order); assert(FALSE); }
    
    if (N > 0)
      { pqueue_item_t Q_head = pqueue_head(Q);
        pqueue_item_t head = (order > 0 ? item[0] : item[N-1]);
        if (verbose || (Q_head != head)) { ER("  pqueue_head(Q) = %08u\n", Q_head); }
        if (Q_head != head) { ER("  ** should be %08u\n", head); assert(FALSE); }
      }
     
    /* Check items in queue: */
    for (int32_t i = 0; i < N; i++)
      { pqueue_item_t z = item[i];
        pqueue_value_t v = ival[i];
        bool_t has = TRUE;
        
        bool_t Q_has = pqueue_has(Q, z);
        if (verbose || (Q_has != has)) { ER("  pqueue_has(Q, %08u) = %c\n", z, "FT"[Q_has]); }
        if (Q_has != has) { ER("  ** should be %c\n", "FT"[has]); assert(FALSE); }
        
        pqueue_value_t Q_v = pqueue_value(Q, z);
        if (verbose || (Q_v != v)) { ER("  pqueue_value(Q, %08u) = %18.16f\n", z, Q_v); }
        if (Q_v != v) { ER("  ** should be %18.16f\n", v); assert(FALSE); }
        
        pqueue_position_t Q_p = pqueue_position(Q, z);
        if (verbose || (Q_p >= N)) { ER("  pqueue_position(Q, %08u) = %u\n", z, Q_p); }
        if (Q_p >= N) { ER("  ** should be in {0 .. %d}\n", N); assert(FALSE); }
        
        bool_t Q_z = pqueue_item(Q, Q_p);
        if (verbose || (Q_z != z)) { ER("  pqueue_item(Q, %u) = %08u\n", z, Q_z); }
        if (Q_z != z) { ER("  ** should be %08u\n", z); assert(FALSE); }
      }
      
    /* Check an item outside the queue: */
    pqueue_item_t z_o = throw_new_item(N, zlim, item);
    bool_t has_o = FALSE;

    bool_t Q_has_o = pqueue_has(Q, z_o);
    if (verbose || (Q_has_o != has_o)) { ER("  pqueue_has(Q, %08u) = %c\n", z_o, "FT"[Q_has_o]); }
    if (Q_has_o != has_o) { ER("  ** should be %c\n", "FT"[has_o]); assert(FALSE); }

    pqueue_position_t Q_p_o = pqueue_position(Q, z_o);
    if (verbose || (Q_p_o != N)) { ER("  pqueue_position(Q, %08u) = %u\n", z_o, Q_p_o); }
    if (Q_p_o != N) { ER("  ** should be %d\n", N); assert(FALSE); }
  }
  
pqueue_item_t throw_new_item(uint32_t N, pqueue_item_t zlim, pqueue_item_t item[])
  {
    while(TRUE)
      { pqueue_item_t z = uint32_abrandom(0, zlim-1);
        pqueue_position_t i = 0;
        while((i < N) && (item[i] != z)) { i++; }
        if (i >= N) { return z; }
      }
  }
  
pqueue_value_t throw_distinct_value(uint32_t N, pqueue_value_t ival[])
  {
    while(TRUE)
      { pqueue_value_t v = drandom();
        pqueue_position_t i = 0;
        while((i < N) && (ival[i] != v)) { i++; }
        if (i >= N) { return v; }
      }
  }
  
uint32_t throw_index(uint32_t N)
  {
    return uint32_abrandom(0, N-1);
  }
  
