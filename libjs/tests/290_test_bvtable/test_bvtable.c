#define PROG_NAME "test_bvtable"
#define PROG_DESC "tests the {bvtable.h} procedures"
#define PROG_VERS "1.1"

/* Last edited on 2024-11-18 09:14:11 by stolfi */
/* Created on 2007-01-31 by J. Stolfi, UNICAMP */

#define PROG_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <limits.h>

#include <vec.h>
#include <bool.h>
#include <ref.h>
#include <affirm.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <jstime.h>
#include <bvhash.h>

#include <bvtable.h>

int32_t main (int32_t argc, char **argv);

typedef unsigned char ubyte;

void test_bvtable_suite(uint32_t nItems);
  /* Tests the {bvtable.h} operations with {nItems} items,
    both with generour preallocation and with automatic expansion
    from minimum size. The test is repeated for several 
    item sizes. */

void test_bvtable_static(uint32_t nItems, size_t szRel, size_t sz, uint32_t nTimes, bool_t preAlloc);
  /* Tests the {bvtable.h} operations with {nItems}, each {sz} bytes 
    long, with only the first {szRel} being relevant for 
    item equivalence. If {preAlloc} is true, pre-allocates the table
    with the {nItems} entries, else starts with a
    minimal table and lets it expand automatically.
    The timing loops execute each operation about {nTimes} times. */

/* CARACTER STRING EQUIVALENCE AND HASHING */

uint64_t ubytes_hash(ubyte *p, size_t szRel);
  /* Computes a hash function from the frist {szRel} bytes starting at {*p}.
    Fails if {p} is NULL or {n} is 0. Cost: {K*szRel}. */

int32_t ubytes_cmp(ubyte *x, ubyte *y, size_t szRel);
  /* Compares lexicographically the first {szRel} unsigned 
    bytes that start at {*x} and {*y} respectively. 
    Returns {-1,0,+1} depending on the {x} string being less than,
    equal to, or greater than the {y} string.
    
    Fails if either {x} or {y} is NULL. Cost: {K*szRel}. */

void create_items
  ( uint32_t nItems,
    size_t szRel,
    size_t sz, 
    ubyte *item[], 
    uint32_t eqix[]
  );
  /* Carves the {zone} array into {nItems} disjoint items with {sz}
    bytes each, fills them with random data, and saves their addresses
    in {item[0..nItems-1]}, in some random order.
    
    Two items are considered equivalent if they have the same contents
    in their first {szRel} bytes, irrespective of the contents of the
    remaining {sz-szRel} bytes. The table will contain a certain
    fraction of equivalent items, and maybe some identical ones. The
    procedure will set {eqix[i]}, for each {i} in {0..nItems-1}, to the
    smallest {j} so that item {j} is equivalent to item {i}.
    
    The reason to carve the items out of a static area, instead of
    using {malloc} for each item, is that the addresses returned by
    {malloc} are random and do not repeat between two runs of the
    program; which is terrible for debugging. */

void fill_item(ubyte *p, size_t szRel, uint32_t ic, size_t sz, uint32_t ii);
  /* Fills the frist {szRel} bytes of {*p} with the integer value 
    {ic}, scrambled.  The value of {ic} must be less than {256^szRel}.
    Then fills the remaining {sz-szRel} bytes with some garbage that
    usually depends on {ii}. */
    
void test_bvtable_correctness(uint32_t nItems, uint32_t nGuess, size_t szRel, size_t sz, ubyte *item[], uint32_t eqix[]);
// void test_bvtable_speed(uint32_t nItems, ubyte *item[], uint32_t eqix[], bvtable_t *S, int32_t nTimes, bool_t strings);
// void print_timing(char *func, double usec, int32_t nops);

#define MAX_ITEM_BYTES (64*1024*1024)
  /* Max total byte size of all items in test set. */

static ubyte zone[MAX_ITEM_BYTES];
  /* The area where items are allocated from. */

int32_t main (int32_t argc, char **argv)
  {
    srandom(4615);
    fprintf(stderr, "random() = %ld\n", random());
    srand(4615);
    fprintf(stderr, "rand() = %d\n", rand());
    
    /* Testing with various table sizes: */
    test_bvtable_suite(   100);
    test_bvtable_suite( 10000);

    return 0;
  }
  
void test_bvtable_suite(uint32_t nItems)
  {
    uint32_t nTimes = 100000; /* Number of timing calls per function. */
    for (int32_t pa = 0; pa < 2; pa++)
      { test_bvtable_static(nItems, 1,   2,   nTimes, (pa == 0));
        test_bvtable_static(nItems, 10,  20,  nTimes, (pa == 0));
        test_bvtable_static(nItems, 150, 200, nTimes, (pa == 0));
      }
  }

void test_bvtable_static(uint32_t nItems, size_t szRel, size_t sz, uint32_t nTimes, bool_t preAlloc)
  { 
    fprintf(stderr, "============================================================\n");
    fprintf(stderr, "testing with %u items of %lu bytes (first %lu relevant)\n", nItems, sz, szRel);
    fprintf(stderr, "preallocation = %c\n", "FT"[preAlloc]);
    
    /* Create the items {item[0..nItems-1]} to store in the set: */
    fprintf(stderr, "creating the items ...\n");
    ubyte *item[nItems]; /* The items to use for the test. */
    uint32_t eqix[nItems];   /* {eqix[i]} is the smallest {j} such that {cmp(item[i],item[j],sz) is 0}. */
    create_items(nItems, szRel, sz, item, eqix);
          
    /* Allocate the set for the specified occupancy ratio: */
    fprintf(stderr, "testing {bvtable_new} ...\n");
    uint32_t nGuess = (preAlloc ? nItems : 0);
    
    test_bvtable_correctness(nItems, nGuess, szRel, sz, item, eqix);
    
    // test_bvtable_speed(nItems, nGuess, szRel, sz, item, eqix, nTimes);
    
    fprintf(stderr, "============================================================\n");
    return;
  }

void create_items
  ( uint32_t nItems, 
    size_t szRel, 
    size_t sz, 
    ubyte *item[], 
    uint32_t eqix[]
  )
  {
    /* Sanity checks: */
    demand((sz > 0) && (nItems > 0), "duh?");
    demand((szRel > 0) && (szRel <= sz), "bad szRel");
    demand(nItems*sz <= MAX_ITEM_BYTES, "test requires for too much memory");
    
    /* Decide how many equivalence classes: */
    uint32_t nClasses = (3 * nItems /4);
    if (szRel < 4)
      { /* The relevant bytes are less than 32 bits. Make sure that classes can be stored: */
        uint32_t maxClasses = 1u << (8*szRel); /* Max num of classes identifiable with {szRel} byes. */
        if (nClasses > maxClasses) { nClasses = maxClasses; }
      }
    assert(nClasses <= nItems);

    /* Table of first item index in each equivalence class, or {UINT32_MAX}: */
    uint32_t *first = notnull(malloc(nClasses*sizeof(uint32_t)), "no mem");
    for (int32_t ic = 0; ic < nClasses; ic++) { first[ic] = UINT32_MAX; }
    
    /* Carve the items out of {zone}: */
    ubyte *next = zone;
    for (int32_t ii = 0; ii < nItems; ii++)
      { item[ii] = next;
        next += sz;
        assert(next <= zone + MAX_ITEM_BYTES); 
        /* Pick an equivalence class {ic} for this item: */
        uint32_t ic = uint32_abrandom(0, (uint32_t)(nClasses-1));
        /* Fill the item with contents belonging to class {ic}: */
        fill_item(item[ii], szRel, ic, sz, ii);
        if (first[ic] == UINT32_MAX)
          { /* First item of this class: */
            eqix[ii] = ii; 
            first[ic] = ii;
          }
        else
          { /* Class has been used before: */
            eqix[ii] = first[ic];
          }
      }
      
    /* Dump some items and check the equivalence: */
    for (int32_t ii = 0; ii < nItems; ii++)
      { uint32_t ie = eqix[ii];
        bool_t debug = (ii < 10);
        if (debug) 
          { fprintf(stderr, "  item[%d] = %16p", ii, item[ii]);
            if (ie != ii)  { fprintf(stderr, "  eq to item[%d]", ie); }
            fprintf(stderr, "\n");
          }
        assert(ubytes_cmp(item[ii], item[ie], szRel) == 0);
      }
    
    free(first);
  }
        
void fill_item(ubyte *p, size_t szRel, uint32_t ic, size_t sz, uint32_t ii)
  {
    uint32_t xc = ic;
    for (size_t k = 0; k < szRel; k++)
      { (*p) = (ubyte)((xc & 255u) ^ 417u ^ (uint32_t)k);
        p++; xc = xc >> 8;
        if (xc == 0) { xc = (ic ^ (uint32_t)(k << 2)); } 
      }
    /* Complete the string with random bits: */
    for (size_t k = szRel; k < sz; k++)
      { (*p) = (ubyte)(int32_abrandom(0, 255));
        p++; 
      }
  }

void test_bvtable_correctness(uint32_t nItems, uint32_t nGuess, size_t szRel, size_t sz, ubyte *item[], uint32_t eqix[])
  {
    fprintf(stderr, "TESTING CORRECTNESS\n");
    uint32_t n_new = 0; /* Number of non-equivalent items added to {S}. */

    bvtable_t *S = bvtable_new(sz, nGuess);

    auto int32_t cmp_proc(void *x, void *y, size_t szS);
    auto uint64_t hash_proc(void *p, size_t szS);
    
    /* Add the items to the set: */
    fprintf(stderr, "testing {bvtable_get_index}, {bvtable_add} ...\n");
    for (int32_t ii = 0; ii < nItems; ii++)
      { ubyte *pi = item[ii];
        assert(pi != NULL);
        /* Get the earliest item {pe} that is equiv to {pi} and is in {S}, or 0 if nonesuch. */
        uint32_t ie = eqix[ii]; /* Index of item equivalent to {pi}, or {-1}. */
        ubyte *pe = item[ie];
        /* Perform some operations while adding {pi} to {S}: */
        uint32_t ir = bvtable_get_index(S, pi, &hash_proc, &cmp_proc);
        if (ie == ii) 
          { /* There should be no item in {S} equivalent to {pi}. */
            affirm(ir == UINT32_MAX, "{bvtable_get_index} error 1 (found when shouldn't)");
            uint32_t r0 = n_new; /* Where the item should be added. */
            ir = bvtable_add(S, pi, &hash_proc, &cmp_proc);
            n_new++;
            affirm(ir == r0, "{bvtable_add} error 1 (returned the wrong index)");
            uint32_t r3 = bvtable_get_index(S, pi, &hash_proc, &cmp_proc);
            affirm(r3 < n_new, "{bvtable_get_index} error 2 (cannot find added item)");
            affirm(r3 == r0, "{bvtable_get_index} error 3 (found added item in wrong place)");
          }
        else
          { /* The item {ii} should have been inserted before as item {ie}: */
            assert(ie < ii);
            affirm(ir != UINT32_MAX, "{bvtable_get_index} error 5 (failed to find equivalent in S)");
            affirm(ir < n_new, "{bvtable_get_index} error 6 (returned invalid index)");
            uint32_t r0 = bvtable_get_index(S, pe, &hash_proc, &cmp_proc);
            affirm(ir == r0, "{bvtable_get_index} error 7 (inconsistent for equivalents)");
          }
        affirm(bvtable_item_count(S) == n_new, "{bvtable_item_count} error 8 (count mismatch)");
        ubyte *ps = bvtable_item_address(S, ir);
        demand(ps != NULL, "{bvtable_item_address} error 9 (returned NULL)");
        demand(ubytes_cmp(pe, ps, sz) == 0, "{bvtable_item_address} error 10 (item not there)");
      }
    
    fprintf(stderr, "closing the table ...\n");
    ubyte *item_fin = NULL;
    uint32_t n_fin = UINT32_MAX;
    bvtable_close(S, &n_fin, (void**)&(item_fin));
    fprintf(stderr, "closed (ne = %u)\n", n_fin);
    affirm(n_fin == n_new, "{bvtable_close} error 1 (count does not match insertions)");
    affirm((n_fin == 0) || (item_fin != NULL), "{bvtable_close} error 2 (no pointer returned)");
    
    /* Check that the items are the first items of each class: */
    ubyte *p_fin = item_fin;
    uint32_t i_fin = 0;
    for (int32_t ii = 0; ii < nItems; ii++)
      { uint32_t ie = eqix[ii]; /* Index of item equivalent to {pi}, or {-1}. */
        if (ie == ii)
          { /* First of its class. */
            affirm(i_fin < n_fin, "{bvtable_close} error 3 (insuff items stored)");
            ubyte *pe = item[ie];
            affirm(cmp_proc(p_fin, pe, sz) == 0, "{bvtable_close} error 3 (item stored not eq)");
            affirm(bcmp(p_fin, pe, sz) == 0, "{bvtable_close} error 4 (item stored not first eq)");
            p_fin += sz;
            i_fin++;
          }
      }
    fprintf(stderr, "\n");
    fprintf(stderr, "DONE TESTING CORRECTNESS\n");
    fprintf(stderr, "\n");
    
    free(item_fin);
    return;

    int32_t cmp_proc(void *x, void *y, size_t szS)
      { 
        assert(szS == sz);
        return ubytes_cmp((ubyte*)x, (ubyte*)y, szRel); 
      }
      
    uint64_t hash_proc(void *p, size_t szS)
      {
        assert(szS == sz);
        return ubytes_hash((ubyte*)p, szRel); 
      }

  }

// void test_bvtable_speed(int32_t nItems, ubyte *item[], int32_t eqix[], bvtable_t *S, int32_t nTimes, bool_t strings)
//   {
//     fprintf(stderr, "TESTING SPEED\n");
//     
//     int32_t i;
//     double start, stop; /* Clock readings. */
//     
//     /* Pick a number {step} that is relatively prime to {nItems}: */
//     int32_t step = (int32_t)(0.61803398874989484820 * nItems);
//     while(gcd(step, nItems) != 1) { step--; }
// 
//     /* Measure mean time of {bvtable_append_last} from empty to full: */
//     bvtable_stats_clear(S);
//     double tAdd = 0;
//     int32_t kAdd = 0; /* Next item to add. */
//     int32_t nAdd = 0;   /* Number of calls to {bvtable_insert_last}. */
//     while(nAdd < nTimes)
//       { /* Clear the list and insert all items: */
//         bvtable_clear(S);
//         start = user_cpu_time_usec();
//         for (int32_t i = 0; i < nItems; i++)
//           { bvtable_item_t pi = UITEM(item[kAdd]);
//             (void)bvtable_insert_last(S, pi);
//             nAdd++;
//             kAdd = (kAdd + step) % nItems; 
//           }
//         stop = user_cpu_time_usec();
//         tAdd += stop - start;
//       }
//     print_timing("bvtable_insert_last",  tAdd, nAdd);
//     bvtable_stats_print(S);
// 
//     /* Measure mean time of {bvtable_get_index} in full table: */
//     bvtable_stats_clear(S);
//     double tInd = 0;
//     int32_t kInd = 0;     /* Next item to look up. */
//     int32_t nInd = 0;   /* Number of calls to {bvtable_get_index}. */
//     while(nInd < nTimes)
//       { /* Clear the list and insert all items: */
//         bvtable_clear(S);
//         for (int32_t i = 0; i < nItems; i++)
//           { bvtable_item_t pi = UITEM(item[kAdd]);
//             (void)bvtable_insert_last(S, pi);
//             kAdd = (kAdd + step) % nItems;
//           }
//         /* Look up all items: */
//         start = user_cpu_time_usec();
//         for (int32_t i = 0; i < nItems; i++)
//           { bvtable_item_t pi = UITEM(item[kInd]);
//             (void)bvtable_get_index(S, pi);
//             nInd++;
//             kInd = (kInd + step) % nItems;
//           }
//         stop = user_cpu_time_usec();
//         tInd += stop - start;
//         nInd += nItems;
//       }
//     print_timing("bvtable_get_index",  tInd, nInd);
//     bvtable_stats_print(S);
//   }
// 
// void print_timing(char *func, double usec, int32_t nops)
//   {
//     fprintf(stderr, "%-25s  %13.0f usec / %10d ops = %13.6f usec/op\n", func, usec, nops, usec/nops);
//   }


uint64_t ubytes_hash(ubyte *p, size_t szRel)
  { 
    return bvhash_bytes(p, szRel);
  }

int32_t ubytes_cmp(ubyte *x, ubyte *y, size_t szRel)
  { if (x == y) { return 0; }
    while (szRel > 0)
      { ubyte xk = (*x); x++;
        ubyte yk = (*y); y++;
        if (xk < yk)
          { return -1; }
        if (xk > yk)
          { return +1; }
        szRel--;
      }
    return 0;
  }
