#define PROG_NAME "test_ulist"
#define PROG_DESC "tests the {ulist.h} procedures"
#define PROG_VERS "1.1"

/* Last edited on 2018-03-04 22:57:20 by stolfilocal */
/* Created on 2007-01-31 by J. Stolfi, UNICAMP */

#define PROG_COPYRIGHT \
  "Copyright © 2007  by the State University of Campinas (UNICAMP)"

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

#include <ulist.h>

int main (int argc, char **argv);

void test_ulist_suite(int nItems, int nSlots);
  /* Tests the {ulist.h} operations on a set vector with static size
    {nSlots} and containing at most {nItems} items. The test is
    repeated for several combinations of item sizes, sometimes with
    opaque addresses, sometimes with strings under {strcmp}
    equivalence. */

void test_ulist_static(int nItems, int nSlots, int nTimes, int szMin, int szMax, bool_t strings);
  /* Tests the {ulist.h} operations on a list {S} with static size
    {nSlots} and containing at most {nItems} items. The timing
    loops execute each operation about {nTimes} times. The test
    items are obtained by allocating {nItems} records with
    sizes in the range {szMin..szMax}. If {strings} is TRUE, uses the
    string hash functions, else uses the default address hash
    functions. */

void create_items
  ( int nItems, 
    int szMin, 
    int szMax, 
    bool_t strings, 
    ref_t item[], 
    int eqix[]
  );
  /* Carves the {zone} array into {nItems} disjoint items with sizes in the range
    {szMin..szMax}, and saves their addresses in {item[0..nItems-1]},
    in some random order. 
    
    If {strings} is FALSE, the contents of the item records are left
    undefined, and all entries {eqix[0..nItems-1]} are all set to {-1}.
    
    If {strings} is TRUE, the items records are filled with
    zero-terminated letter strings, such that at most two strings will
    have the same contents. For any pair of strings {item[i],item[j]}
    with {i != j} and same contents, the procedure sets {eqix[i] = j}
    and {eqix[j] = i}. For all the remaining items, it sets {eqix[i] =
    -1}. In this case, one must specify {szMin == szMax >= 7}.
    
    The reason to carve the items out of a static area, instead of
    using {malloc} for each item, is that the addresses returned by
    {malloc} are random and do not repeat between two runs of the
    program; which is terrible for debugging. */

void test_ulist_hash_eq(int nItems, ref_t item[], int eqix[], ulist_t *S);
void test_ulist_correctness(int nItems, ref_t item[], int eqix[], ulist_t *S, bool_t strings);
void test_ulist_speed(int nItems, ref_t item[], int eqix[], ulist_t *S, int nTimes, bool_t strings);
void print_timing(char *func, double usec, int nops);

unsigned int string_eq_class_index(unsigned int i);
  /* Maps an integer {i} to some integer {j} in {0..i}, in
    such a way that there are at most 2 distinct integers {i1,i2}
    such that {string_eq_class_index(i1) == string_eq_class_index(i2)}. */
    
unsigned int string_eq_class_mate(unsigned int i);
  /* If there is an integer {j != i} such that 
    {string_eq_class_index(j) = string_eq_class_index(i)},
    returns that {j}; otherwise returns -1. */

#define MAX_ITEM_BYTES (64*1024*1024)

static char zone[MAX_ITEM_BYTES];
  /* The area where items are allocated from. */

ulist_item_t UITEM(char *x);
  /* Maps the address {x} to an offset from {zone} so that it will
    fit in 32 bits.  Added 1 so that {UITEM(x)} is zero iff {x == NULL}. */

char *PCHAR(ulist_item_t x);
  /* Reverses the mapping {UITEM} and casts the resulting address
    to {char *}. */

/* CARACTER STRING EQUIVALENCE AND HASHING */

ulist_hash_val_t string_hash(ulist_item_t a, ulist_hash_size_t n);
  /* Assumes that {a} is actually the memory address (of type {char*})
    a zero-terminated character string. Returns a hash value in {0..n-1},
    computed solely from the characters in that string, independently
    of the string's address {a}. Fails if {a} is NULL or {n} is 0.
    Cost: {K*strlen(a)}. */

bool_t string_eq(ulist_item_t a, ulist_item_t b);
  /* Assumes that {a} and {b} are actually the memory addresses (of
    type {char*}) of two zero-terminated character strings. Returns
    TRUE if and only if those strings are byte-by-byte equal, as per
    {strcmp}, independently of whether their addresses {a} and {b} are
    the same. Fails if either {(char*)a} or {(char*)b} is NULL. Cost:
    {K*min(strlen((char*)a),strlen((char*)b))}. */

int main (int argc, char **argv)
  {
    srandom(4615);
    fprintf(stderr, "random() = %ld\n", random());
    srand(4615);
    fprintf(stderr, "rand() = %d\n", rand());
    
    /* Testing with third-full tables (too generous): */
    test_ulist_suite(   100,   300);
    test_ulist_suite( 10000, 30000);
    
    /* Testing with half-full tables (more typical): */
    test_ulist_suite(   100,   200);
    test_ulist_suite( 10000, 20000);
    
    /* Testing with full tables (mainly for correctness): */
    test_ulist_suite(   100,   100);
    test_ulist_suite(  1000,  1000);
    test_ulist_suite( 10000, 10000);

    return 0;
  }
  
void test_ulist_suite(int nItems, int nSlots)
  {
    int nTimes = 100000; /* Number of timing calls per function. */
    /* Tests with addresses: */
    test_ulist_static(nItems, nSlots, nTimes,    1,    1, FALSE);
    test_ulist_static(nItems, nSlots, nTimes,  100,  200, FALSE);
    test_ulist_static(nItems, nSlots, nTimes,    1, 1000, FALSE);
    /* Tests with strings: */
    test_ulist_static(nItems, nSlots, nTimes,    8,    8,  TRUE);
  }

void test_ulist_static(int nItems, int nSlots, int nTimes, int szMin, int szMax, bool_t strings)
  { 
    fprintf(stderr, "============================================================\n");
    fprintf(stderr, "testing with %d items in %d slots", nItems, nSlots);
    fprintf(stderr, " (occupancy ratio %6.4f)\n", ((double)nItems)/((double)nSlots));
    fprintf(stderr, "record size range = [%d .. %d]\n", szMin, szMax);
    fprintf(stderr, "equality/hash = %s\n", (strings ? "strings" : "addresses"));
    
    /* Create the items {item[0..nItems-1]} to store in the set: */
    fprintf(stderr, "creating the items ...\n");
    ref_t item[nItems]; /* The items to use for the test. */
    int eqix[nItems];   /* {eqix[i]} is the {j!=i} such that {eq(item[i],item[j])}, or {-1} */
    create_items(nItems, szMin, szMax, strings, item, eqix);
          
    /* Allocate the set for the specified occupancy ratio: */
    fprintf(stderr, "testing {ulist_new} ...\n");
    ulist_t *S = ulist_new(nSlots);
    affirm(ulist_capacity(S) == nSlots, "wrong refset size");
    if (strings)
      { ulist_install_eq(S, &string_eq);
        ulist_install_hash(S, &string_hash);
      }
    
    test_ulist_hash_eq(nItems, item, eqix, S);
    
    test_ulist_correctness(nItems, item, eqix, S, strings);
    
    if (! strings) { test_ulist_speed(nItems, item, eqix, S, nTimes, strings); }
    
    /* Reclaim test storage: */
    ulist_free(S);
    
    fprintf(stderr, "============================================================\n");
    return;
  }
    
void create_items
  ( int nItems, 
    int szMin, 
    int szMax, 
    bool_t strings, 
    ref_t item[], 
    int eqix[]
  )
  {
    /* Sanity check for memory size: */
    demand(nItems*szMax <= MAX_ITEM_BYTES, "test requires for too much memory");
    
    if (strings) 
      { demand(szMin >= 7, "{szMin} too small for strings");
        demand(szMin == szMax, "{szMin} must be equal to {szMax} for strings");
      }
    
    int i; 
    /* Carve the items out of {zone}: */
    char *next = zone;
    for (i = 0; i < nItems; i++)
      { size_t sz = int32_abrandom(szMin, szMax);
        item[i] = next;
        next += sz;
        assert(next <= zone + MAX_ITEM_BYTES); 
        eqix[i] = -1; /* By default. */
        if (strings)
          { /* Map {i} to some equivalence class index {vi}. */
            int vi = string_eq_class_index(i);
            char *p = item[i];
            /* Set {*(item[i])} to be {vi} in reverse base 26 with digits [a-z]: */
            int k;
            for (k = 0; k < sz-1; k++)
              { (*p) = (char)('a' + (vi % 26));
                p++; vi /= 26;
              }
            /* Terminate the string with a zero byte: */
            (*p) = 0;
            /* Compute {eqi != i} such that {eq(item[i],item[eqi])}: */
            int eqi = string_eq_class_mate(i);
            /* If the mate exists and is already in {item}, set {eqix} accordingly: */
            if ((eqi != -1) && (eqi < i)) 
              { assert(strcmp(item[i],item[eqi]) == 0);
                eqix[i] = eqi; eqix[eqi] = i;
              }
          }
      }

    /* Apply a random permutation to the items: */
    for (i = 1; i < nItems; i++)
      { int j = int32_abrandom(0,i);
        if (j < i) 
          { /* Grab their eq indices {eqi,eqj}: */
            int eqi = eqix[i], eqj = eqix[j]; 
            assert((eqi != i) && (eqj != j));
            /* Swap {item[i]} with {item[j]}: */
            { ref_t t = item[i]; item[i] = item[j]; item[j] = t; }
            /* Fix the equivalence indices: */
            if ((eqi == j) && (eqj == i))
              { /* Nothing to fix, just check: */
                assert(strcmp(item[i],item[j]) == 0);
              }
            else
              { /* Either they are eq or not: */
                assert((eqi != j) && (eqj != i));
                /* Swap the double-ended pointers: */
                eqix[j] = eqi; if (eqi != -1) { eqix[eqi] = j; }
                eqix[i] = eqj; if (eqj != -1) { eqix[eqj] = i; }
              }
          }
      }
      
    /* Dump some items and check the equivalence: */
    for (i = 0; i < nItems; i++)
      { bool_t debug = (i < 10);
        if (debug) 
          { fprintf(stderr, "  item[%d] = %16p", i, item[i]);
            if (strings) { fprintf(stderr, " = \"%s\"", (char*)(item[i])); }
            int eqi = eqix[i];
            if (eqi != -1) 
              { assert(strings);
                fprintf(stderr, "  eq to item[%d]", eqi);
                fprintf(stderr, " = \"%s\"", (char*)(item[eqi]));
              }
            fprintf(stderr, "\n");
          }
        int j = eqix[i];
        if (j != -1)
          { assert(strings);
            assert(string_eq(UITEM(item[i]), UITEM(item[j])));
          }
      }
  }
  
unsigned int string_eq_class_index(unsigned int i)
  {
    return i & ((i >> 1) | (~ 2u));
  }

unsigned int string_eq_class_mate(unsigned int i)
  {
    if ((i & 4u) != 0)
      { /* Index has no equivalents: */ return -1; }
    else
      { /* Index has one equivalent: */ return (i ^ 2u); }
  }

ulist_item_t UITEM(char *x)
  { return (x == NULL ? 0 : ((ulist_item_t)(1 + (uint32_t)(((char *)x) - &(zone[0]))))); }

char *PCHAR(ulist_item_t x)
  { return ((char*)(x == 0 ? NULL : (&(zone[0]) + (uint32_t)x - 1))); }

void test_ulist_hash_eq(int nItems, ref_t item[], int eqix[], ulist_t *S)
  {
    fprintf(stderr, "TESTING HASH AND EQUALITY\n");
    
    ulist_hash_func_t *hash = ulist_get_hash(S); 
    ulist_eq_func_t *eq = ulist_get_eq(S); 
    /* Histogram of hash results: */
    ulist_hash_size_t nh = ulist_hash_size(S);
    uint32_t hct[nh]; /* {hct[h]} is how many non-equivalent items are hashed to {h}. */
    ulist_hash_val_t h;
    for (h = 0; h < nh; h++) { hct[h] = 0; }
    /* Check all items: */
    int i;
    for (i = 0; i < nItems; i++) 
      { h = hash(UITEM(item[i]), nh);
        affirm(h < nh, "{S.hash} returns out-of-bounds result");
        int j = eqix[i];
        if (j >= 0)
          { /* Check whether {item[i]} and {item[j]} are equivalent by {eq}: */
            affirm(eq(UITEM(item[i]), UITEM(item[j])), "{S.eq} error");
            /* Check whether {item[i]} and {item[j]} hash to the same key: */
            ulist_hash_val_t g = hash(UITEM(item[j]), nh); 
            affirm(h == g, "{S.hash} is inconsistent with {S.eq}"); 
          }
        else
          { /* Should check {!eq(UITEM(item[i]), UITEM(item[j]))} for all {j < i}. */
            /* Update histogram: */
            hct[h]++;
          }
      }
    /* Print large entries from the histogram and compute a mean collision estimate. */
    /* The estimated total probes for a bucket of size {m} is {m*(m+1)/2}. */
    /* This estimate ignores the merging of buckets that occurs in linear hashing. */
    int64_t tm2 = 0; /* Sum of {hct[i]*(hct[i]+1)} for {i} in {0..nh-1}. */
    int szct[nItems+1]; /* {szct[m]} is the number of keys that are shared by {m} items. */
    int m;
    for (m = 0; m <= nItems; m++) { szct[m] = 0; }
    for (h = 0; h < nh; h++)
      { int m = hct[h];
        assert(m <= nItems);
        szct[m]++;
        tm2 += ((uint64_t)m)*((uint64_t)m+1);
        if (hct[h] > 2)
          { fprintf(stderr, "%7d items hash to %7d\n", hct[h], h); }
      }
    double eppo = 0.5*((double)tm2)/((double)nItems);
    fprintf(stderr, "estimated probes per operation = %8.2f\n", eppo);
    /* Print the histogram of the histogram: */
    for (m = 0; m <= nItems; m++)
      { if ((m < 4) || (szct[m] > 0))
          { fprintf(stderr, "there are %7d buckets with %7d entries\n", szct[m], m); }
      }
  }

void test_ulist_correctness(int nItems, ref_t item[], int eqix[], ulist_t *S, bool_t strings)
  {
    fprintf(stderr, "TESTING CORRECTNESS\n");
    int i;
    int nne = 0; /* Number of non-equivalent items added to {S}. */
    
    affirm(ulist_count(S) == 0, "{ulist_count} error 1 (not zero initially)");

    /* Add the items to the set: */
    fprintf(stderr, "testing {ulist_add}, {ulist_has} ...\n");
    for (i = 0; i < nItems; i++)
      { /* Conver the pointer {item[i]} to an {ulist_item_t}: */
        ulist_item_t a = UITEM(item[i]);
        assert(a != 0);
        /* Get the item {e} that is equiv to {a} and is in {S}, or 0 if nonesuch. */
        int eqi = eqix[i]; /* Index of item equivalent to {a}, or {-1}. */
        ulist_item_t e = ((eqi >= 0) && (eqi < i) ? UITEM(item[eqi]) : 0);
        /* Perform some operations while adding {a} to {S}: */
        ulist_index_t r0 = ulist_count(S);
        affirm(r0 == nne, "{ulist_count} error 0 (count does not match insertions)");
        ulist_index_t r1 = ulist_index_of(S, a);
        bool_t h1 = ulist_has(S, a);
        if (e == 0) 
          { /* There should be no item in {S} equivalent to {a}. */
            affirm(r1 >= nne, "{ulist_index_of} error 1 (found when shouldn't)");
            affirm(! h1, "{ulist_has} error 1 (found when shouldn't)");
            ulist_index_t r2 = ulist_insert_last(S, a);
            nne++;
            affirm(r2 == r0, "{ulist_insert_last} error 1 (returned the wrong index)");
            ulist_index_t r3 = ulist_index_of(S, a);
            affirm(r3 < nne, "{ulist_index_of} error 2 (cannot find added item)");
            affirm(r3 == r0, "{ulist_index_of} error 3 (found added item in wrong place)");
            bool_t h3 = ulist_has(S, a);
            affirm(h3, "{ulist_has} error 2 (says added item is not there)");
            ulist_item_t a4 = ulist_item_at(S, r3);
            affirm(a4 == a, "{ulist_item_at} error 2 (returned wrong item)");
          }
        else
          { if (strings) { assert(string_eq(a,e)); }
            affirm(r1 < nne, "{ulist_index_of} error 1 (failed to find equivalent in S)");
            affirm(r1 == ulist_index_of(S, e), "{ulist_index_of} error 2 (inconsistent for equivalents)");
            affirm(h1, "{ulist_has} error 1 (says there is no equivalent in S)");
            ulist_item_t a4 = ulist_item_at(S, r1);
            affirm(a4 == e, "{ulist_item_at} error 2 (returned wrong item)");
          }
        /* Check item count: */
        affirm(ulist_count(S) == nne, "{ulist_count} error 2 (count does not match insertions)");
      }
    
    /* Extract all items in some order: */
    fprintf(stderr, "extracting the items from the set ...\n");
    ulist_item_vec_t E = ulist_items(S);
    affirm(E.ne == nne, "{ulist_items} error 1 (count does not match insertions)");
    for (i = 0; i < E.ne; i++)
      { ulist_item_t a = E.e[i];
        ulist_index_t i = ulist_index_of(S, a);
        ulist_item_t b = ulist_item_at(S, i);
        affirm(a == b, "{ulist_index_of/ulist_item_at} error 7 (did not find item from list)");
        affirm(ulist_has(S, a), "{ulist_has} error 6 (says that item is not there)");
      }
    
    /* Prints operation statistics: */
    free(E.e);
  }

void test_ulist_speed(int nItems, ref_t item[], int eqix[], ulist_t *S, int nTimes, bool_t strings)
  {
    fprintf(stderr, "TESTING SPEED\n");
    
    int i;
    double start, stop; /* Clock readings. */
    
    /* Pick a number {step} that is relatively prime to {nItems}: */
    int step = (int)(0.61803398874989484820 * nItems);
    while(gcd(step, nItems) != 1) { step--; }

    /* Measure mean time of {ulist_append_last} from empty to full: */
    ulist_stats_clear(S);
    double tAdd = 0;
    int kAdd = 0; /* Next item to add. */
    int nAdd = 0;   /* Number of calls to {ulist_insert_last}. */
    while(nAdd < nTimes)
      { /* Clear the list and insert all items: */
        ulist_clear(S);
        start = user_cpu_time_usec();
        for (i = 0; i < nItems; i++)
          { ulist_item_t a = UITEM(item[kAdd]);
            (void)ulist_insert_last(S, a);
            nAdd++;
            kAdd = (kAdd + step) % nItems; 
          }
        stop = user_cpu_time_usec();
        tAdd += stop - start;
      }
    print_timing("ulist_insert_last",  tAdd, nAdd);
    ulist_stats_print(S);

    /* Measure mean time of {ulist_index_of} in full table: */
    ulist_stats_clear(S);
    double tInd = 0;
    int kInd = 0;     /* Next item to look up. */
    int nInd = 0;   /* Number of calls to {ulist_index_of}. */
    while(nInd < nTimes)
      { /* Clear the list and insert all items: */
        ulist_clear(S);
        for (i = 0; i < nItems; i++)
          { ulist_item_t a = UITEM(item[kAdd]);
            (void)ulist_insert_last(S, a);
            kAdd = (kAdd + step) % nItems;
          }
        /* Look up all items: */
        start = user_cpu_time_usec();
        for (i = 0; i < nItems; i++)
          { ulist_item_t a = UITEM(item[kInd]);
            (void)ulist_index_of(S, a);
            nInd++;
            kInd = (kInd + step) % nItems;
          }
        stop = user_cpu_time_usec();
        tInd += stop - start;
        nInd += nItems;
      }
    print_timing("ulist_index_of",  tInd, nInd);
    ulist_stats_print(S);
  }

void print_timing(char *func, double usec, int nops)
  {
    fprintf(stderr, "%-25s  %13.0f usec / %10d ops = %13.6f usec/op\n", func, usec, nops, usec/nops);
  }

unsigned int string_hash(ulist_item_t a, unsigned int n)
  { char *sa = PCHAR(a);
    demand(sa != NULL, "cannot hash NULL");
    demand(n > 0, "cannot hash to empty range");
    /* Let's hope that this works: */
    unsigned uh = 0u;
    unsigned char *p = (unsigned char *)sa;
    while ((*p) != 0)
      { unsigned int uc = (*p);
        unsigned int ur = 
          (4615u * ((uc & 208u) >> 4)) + \
          (417u * ((uc & 13u) << 13)) + \
          (471703u *((uc & 34u) << 3));
        uh = 27183u * uh + ur + 141421 *(uh >> 16);
        p++;
      }
    return (uh % n);
  }

bool_t string_eq(ulist_item_t a, ulist_item_t b)
  { if (a == b) { return TRUE; }
    char *sa = PCHAR(a);
    char *sb = PCHAR(b);
    if ((sa == NULL) || (sb == NULL))
      { return FALSE; }
    else
      { return (strcmp(sa,sb) == 0); }
  }
