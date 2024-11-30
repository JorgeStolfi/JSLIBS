#define PROG_NAME "test_indexing"
#define PROG_DESC "test of {indexing.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-11-23 05:32:59 by stolfi */ 
/* Created on 2005-02-14 (or earlier) by J. Stolfi, UNICAMP */

#define test_indexing_COPYRIGHT \
  "Copyright © 2005  by the State University of Campinas (UNICAMP)"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#include <ix.h>
#include <ix_reduce.h>
#include <affirm.h>
#include <jsrandom.h>
#include <jswsize.h>

#define bug(FMT_AND_ARGS...) \
  do { \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); \
    fprintf(stderr, FMT_AND_ARGS); \
    fprintf(stderr, "\n"); \
    exit(1); \
  } while(0)

#define bug_i_sz_sz(msg,i,taga,sza,tagb,szb)                            \
  bug(("%s  i = %d  %s %" uint64_u_fmt "  %s %" uint64_u_fmt ""), msg, i, taga, sza, tagb, szb)

#define bug_sz_sz(msg,taga,sza,tagb,szb)                            \
  bug(("%s  %s %" uint64_u_fmt "  %s %" uint64_u_fmt ""), msg, taga, sza, tagb, szb)

#define MAXDIM 6

typedef struct desc_t
  { ix_dim_t d;
    ix_size_t sz[MAXDIM];
    ix_pos_t bp;
    ix_step_t st[MAXDIM];
    uint32_t nel;   /* Number of elements allocated for  {el}. */
    uint8_t *el;    /* Element storage area. */
  } desc_t;

desc_t *make_desc(ix_dim_t d, ix_size_t sz[]);
  /* Creates a descriptor with axis count {d} and
    size vector {sz[0..d-1]}.  Some indices {i}
    may have virtual replication ({st[i] == 0}
    with {sz[i] > 1}). */

desc_t *make_packed_desc(ix_dim_t d, ix_pos_t b, ix_size_t sz[], ix_order_t ixor);
  /* Creates a descriptor with axis count {d},
    size vector {sz[0..d-1]}, and base position {b}, with steps consistent
    with a packed element arrangement with index
    order {L}, and no virtual replication. */

void do_test_suite(int32_t nt);
  /* Repeats the full test suite {nt} times. */

void extend_size_vector(ix_dim_t dA, ix_size_t szA[], ix_dim_t dB, ix_size_t szB[]);
  /* Copies {szA[0..dA-1]} to {szB[0..dA-1]} and fills { {szB[dA..dB-1]} with '1'. */

void test_size_ops(ix_dim_t d);
/* void test_new_array(desc_t *A, ix_size_t *sz); */
void test_addressing(desc_t *A);
void test_enum(desc_t *A);
void test_packed(desc_t *A);
/* void test_indices(desc_t *A); */
void test_crop(desc_t *A);
void test_slice(desc_t *A);
void test_subsample(desc_t *A);
void test_swap_indices(desc_t *A);
void test_flip_indices(desc_t *A);
void test_flip(desc_t *A);
void test_diagonal(desc_t *A);
void test_chop(desc_t *A);
void test_replicate(desc_t *A);
void test_reduce(void);

void print_desc(FILE *wr, char *pf, desc_t *A, char *sf);
void print_indices(FILE *wr, char *pf, ix_dim_t d, ix_index_t ix[], char *sf);

/* void dump_storage(FILE *wr, void *el, int32_t nw, ix_nbits_t bpw); */
void check_sizes(ix_dim_t d, ix_size_t sza[], ix_size_t szb[]);
void check_max_indices(ix_dim_t d, ix_index_t ixa[], ix_index_t ixb[]);
void check_max_indices_sz(ix_dim_t d, ix_index_t ixa[], ix_size_t szb[]);

int32_t main (int32_t argn, char **argv)
  {

    /* General limit and data type tests: */
    demand((ix_pos_NONE < 0) || (ix_pos_NONE > ix_MAX_POS), "invalid ix_pos_NONE");
    do_test_suite(30);
    return 0;
  }

void do_test_suite(int32_t nt)
  {
    ix_size_t sz[MAXDIM];
    for (uint32_t it = 0;  it < nt; it++)
      { 
        ix_dim_t d = (ix_dim_t)(it % MAXDIM);
        test_size_ops(d);

        /* An array for testing: */
        sz[0] = uint32_abrandom(0, 2) + uint32_abrandom(0, 2); 
        sz[1] = uint32_abrandom(0, 5) + uint32_abrandom(0, 5); 
        sz[2] = 1; 
        sz[3] = uint32_abrandom(0, 9) + uint32_abrandom(0, 9); 
        sz[4] = uint32_abrandom(0, 7) + uint32_abrandom(0, 7);
        sz[5] = uint32_abrandom(0, 3) + uint32_abrandom(0, 3);
        assert(MAXDIM == 6);

        desc_t *A = make_desc(d, sz);

        /* test_new_array(A, sz, bps, bpw); */
        test_packed(A);
        test_enum(A);
        test_addressing(A);
        test_crop(A);
        test_slice(A);
        test_subsample(A);
        test_swap_indices(A);
        test_flip_indices(A);
        test_flip(A);
        test_diagonal(A);
        test_chop(A);
        test_replicate(A);
        test_reduce();
      }
  }

desc_t *make_desc(ix_dim_t d, ix_size_t sz[])
  {
    fprintf(stderr, "Creating a descriptor (d = %d)...\n", d);
    desc_t *A = talloc(1, desc_t);
    ix_axis_t iz = (ix_axis_t)int32_abrandom(0,2*d+1);     /* The index that will get zero step. */
    A->d = d;
    uint64_t npos = 1;  /* Number of elements already allocated. */
    for (ix_axis_t i = 0; i < d; i++)
      { A->sz[i] = sz[i];
        if ((i == iz) || (sz[i] <= 1))
          { /* Use virtual replication on this index: */
            A->st[i] = 0;
          }
        else
          { /* Use distinct elements on this index: */
            A->st[i] = (sz[i] == 1 ? 0 : (ix_step_t)npos);
          }
        /* Update {npos} if axis is neither trivial nor replicated: */
        if (A->sz[i] == 0)
          { npos = 0; }
        else if ((A->sz[i] > 1) && (A->st[i] != 0))
          { npos *= sz[i]; }
      }
    /* If the array is empty, all steps and the base must be zero: */
    if (npos == 0)
      { for (ix_axis_t i = 0; i < d; i++) { A->st[i] = 0; } 
        A->bp = 0;
      }
    else
      { /* A random base position: */
        A->bp = 1000*uint32_abrandom(1,9);
      }
    if (! ix_parms_are_valid(A->d, A->sz, A->bp, A->st, FALSE))
      { print_desc(stderr, "A = { ", A, " }\n");
        bug("make_desc creates an invalid desc");
        assert(FALSE);
      }
    return A;
  }

desc_t *make_packed_desc(ix_dim_t d, ix_pos_t b, ix_size_t sz[], ix_order_t ixor)
  {
    desc_t *A = talloc(1, desc_t);
    A->d = d;
    ix_sizes_assign(d, A->sz, sz);
    ix_packed_steps(d, sz, ixor, A->st);
    ix_count_t npos = ix_num_tuples(d, sz);
    A->bp = (npos == 0 ? 0 : b);
    if (! ix_parms_are_valid(A->d, A->sz, A->bp, A->st, FALSE))
      { print_desc(stderr, "A = { ", A, " }\n");
        bug("make_desc creates an invalid desc");
        assert(FALSE);
      }
    return A;
  }

void test_size_ops(ix_dim_t d)
  {
    fprintf(stderr, "Checking ix_sizes_assign, ix_sizes_shrink, ix_sizes_expand (d = %d)...\n", d);
    bool_t ch_r, ch_g; 

    /* Size vectors for tests: */
    ix_size_t sza[MAXDIM];
    ix_size_t szb[MAXDIM];
    ix_size_t szc[MAXDIM];
    ix_size_t szd[MAXDIM];

    /* Initialize {sza,szb,szc} with random values: */
    auto ix_size_t rnd_sz(uint32_t a, uint32_t b, uint32_t c, uint32_t i); 
      /* A funct of {a,b,c}, and {i}. */

    ix_size_t rnd_sz(uint32_t a, uint32_t b, uint32_t c, uint32_t i) 
      { return a + ((b*i*i) % c); }

    /* Test {ix_size_assign}: */
    /* fprintf(stderr, "  checking {ix_size_assign}\n"); */
    for (ix_axis_t i = 0; i < MAXDIM; i++) { sza[i] = szc[i] = rnd_sz(17,  8, 13, i); }
    for (ix_axis_t i = 0; i < MAXDIM; i++) { szb[i] = rnd_sz(19, 10, 17, i); }
    ix_sizes_assign(d, sza, szb);
    for (ix_axis_t i = 0; i < MAXDIM; i++)
      { ix_size_t g = (i < d ? szb[i] : szc[i]); 
        if (sza[i] != g) 
          { bug_i_sz_sz("ix_size_assign mismatch", i, "sza =", sza[i], "should be", g); }
      }

    /* Test {ix_size_min}: */
    /* fprintf(stderr, "  checking {ix_size_min}\n"); */
    for (ix_axis_t i = 0; i < MAXDIM; i++) { sza[i] = rnd_sz(17,  8, 13, i); }
    ix_sizes_assign(MAXDIM, szb, sza); /* Save current value of {sza} in {szb}. */
    ch_r = ix_sizes_shrink(d, sza, szc);
    ch_g = FALSE;
    for (ix_axis_t i = 0; i < MAXDIM; i++)
      { ix_size_t g = (i < d ? (szb[i] < szc[i] ? szb[i] : szc[i]) : szb[i]); 
        if (sza[i] != g) 
          { bug_i_sz_sz("ix_size_min mismatch", i, "sza =", sza[i], "should be", g); }
        if ((i < d) && (sza[i] != szb[i])) { ch_g = TRUE; }
      }
    if (ch_r != ch_g)
      { bug("ix_size_min wrong return value"); }

    /* Test {ix_size_max}: */
    /* fprintf(stderr, "  checking {ix_size_max}\n"); */
    for (ix_axis_t i = 0; i < MAXDIM; i++) { sza[i] = rnd_sz(17,  8, 13, i); }
    ix_sizes_assign(MAXDIM, szb, sza); /* Save current value of {sza} in {szb}. */
    ch_r = ix_sizes_expand(d, sza, szd);
    ch_g = FALSE;
    for (ix_axis_t i = 0; i < MAXDIM; i++)
      { ix_size_t g = (i < d ? (szb[i] > szd[i] ? szb[i] : szd[i]) : szb[i]); 
        if (sza[i] != g) 
          { bug_i_sz_sz("ix_size_max mismatch", i, "sza =", sza[i], "should be", g); }
        if ((i < d) && (sza[i] != szb[i])) { ch_g = TRUE; }
      }
    if (ch_r != ch_g)
      { bug("ix_size_max wrong return value"); }
  }

void test_addressing(desc_t *A)
  {
    ix_count_t nt = ix_num_tuples(A->d, A->sz);
    fprintf(stderr, ("Checking ix_position, ix_position_safe, setting and checking all elems (d = %d nt = %" uint64_u_fmt ")...\n"), A->d, nt);
    /* Clear elements: */
    for (ix_axis_t i = 0; i < A->nel; i++) { A->el[i] = 0; }
    bool_t verbose = FALSE;
    if (verbose) { print_desc(stderr, "A = { ", A, " }\n"); } 
    ix_index_t ixA[MAXDIM];

    /* Allocate a large enough boolean array: */
    uint64_t npos = 1;
    for (ix_axis_t i = 0; i < A->d; i++) { npos *= A->sz[i]; }
    uint64_t nel = A->bp + npos + 1000;
    uint8_t *el = talloc(nel, uint8_t);

    /* Get the size vector, extend with 1's: */
    ix_size_t sz[MAXDIM];
    extend_size_vector(A->d, A->sz, MAXDIM, sz);

    /* Reduce the size to 1 whenever the step is 0: */
    for (ix_axis_t i = 0; i < A->d; i++) 
      { if ((A->sz[i] > 1) && (A->st[i] == 0)) { sz[i] = 1; } }

    assert(MAXDIM == 6);
    for (ixA[5] = 0; ixA[5] < sz[5]; ixA[5]++)
      for (ixA[4] = 0; ixA[4] < sz[4]; ixA[4]++)
        for (ixA[3] = 0; ixA[3] < sz[3]; ixA[3]++)
          for (ixA[2] = 0; ixA[2] < sz[2]; ixA[2]++)
            for (ixA[1] = 0; ixA[1] < sz[1]; ixA[1]++)
              for (ixA[0] = 0; ixA[0] < sz[0]; ixA[0]++)
                { 
                  /* Compute the expected address: */
                  ix_pos_t posExp = A->bp;
                  for (ix_axis_t i = 0; i < A->d; i++) { ix_shift_pos(&posExp, A->st[i]*ixA[i]); }

                  /* Compare with library procs: */
                  if (! ix_is_valid(A->d, ixA, A->sz)) { bug("not ix_is_valid"); } 
                  ix_pos_t pos = ix_position(A->d, ixA, A->bp, A->st);
                  if (pos != posExp) 
                    { print_desc(stderr, "A = { ", A, " }\n");
                      print_indices(stderr, "ixA = [ ", A->d, ixA, " ]\n");
                      bug_sz_sz("ix_position mismatch", "pos =", pos, "posExp =", posExp);
                    } 
                  ix_pos_t pos_safe = ix_position_safe(A->d, ixA, A->sz, A->bp, A->st);
                  if (pos != pos_safe) 
                    { print_desc(stderr, "A = { ", A, " }\n");
                      print_indices(stderr, "ixA = [ ", A->d, ixA, " ]\n");
                      bug_sz_sz("ix_position_safe mismatch", "pos =", pos, "pos_safe =", pos_safe);
                    } 
                  uint8_t smp_exp = (uint8_t)(1 + (pos % 253));
                  el[pos] = smp_exp;
                }

    fprintf(stderr, "Re-checking the value of all samples (d = %d)...\n", A->d);
    assert(MAXDIM == 6);
    for (ixA[5] = 0; ixA[5] < sz[5]; ixA[5]++)
      for (ixA[4] = 0; ixA[4] < sz[4]; ixA[4]++)
        for (ixA[3] = 0; ixA[3] < sz[3]; ixA[3]++)
          for (ixA[2] = 0; ixA[2] < sz[2]; ixA[2]++)
            for (ixA[1] = 0; ixA[1] < sz[1]; ixA[1]++)
              for (ixA[0] = 0; ixA[0] < sz[0]; ixA[0]++)
                { 
                  ix_pos_t pos = ix_position(A->d, ixA, A->bp, A->st);
                  uint8_t smp_exp = (uint8_t)(1 + (pos % 253));
                  uint8_t smp = el[pos];
                  if (smp != smp_exp) 
                    { bug("set/get sample mismatch smp = %u smp_exp = %u", smp, smp_exp); } 
              }
  }

void test_enum(desc_t *A)
  {
    /* !!! Test {ix_enum} with early abort when {op} returns {TRUE}. !!! */
    
    ix_count_t nt = ix_num_tuples(A->d, A->sz);
    fprintf(stderr, ("Checking ix_enum (d = %d nt = %" uint64_u_fmt ")...\n"), A->d, nt);

    ix_dim_t d = A->d;

    /* Get the size vector, extend with 1's: */
    ix_size_t sz[MAXDIM];
    extend_size_vector(A->d, A->sz, MAXDIM, sz);

    /* Fabricate the other two arrays: */
    desc_t *B = make_desc(d, sz);
    desc_t *C = make_desc(d, sz);
    bool_t verbose = FALSE;
    if (verbose) 
      { print_desc(stderr, "A = { ", A, " }\n");
        print_desc(stderr, "B = { ", B, " }\n");
        print_desc(stderr, "C = { ", C, " }\n");
      }

    ix_index_t ix[MAXDIM];

    /* Allocate a vector to store all positions of {A}: */
    ix_count_t npos = ix_num_tuples(d, sz);
    ix_pos_t *posA = talloc(npos, ix_pos_t);
    ix_pos_t *posB = talloc(npos, ix_pos_t);
    ix_pos_t *posC = talloc(npos, ix_pos_t);

    uint64_t kpos = 0;

    auto bool_t check_fw( const ix_index_t ix[], ix_pos_t pA, ix_pos_t pB, ix_pos_t pC );
    auto bool_t check_bw( const ix_index_t ix[], ix_pos_t pA, ix_pos_t pB, ix_pos_t pC );
      /* These are procedures that can be passed as the {op} argument of 
        {ix_enum}. They check whether the given positions {pA,pB,pC}
        match {posA[kpos],posB[kpos],posC[kpos]}.  The procedure 
        {check_fw} increments {kpos} after the test,
        while {check_bw} decrements it before the test.
        They always return {FALSE}. */

    /* Enumerate all positions by hand in "F" order: */
    assert(MAXDIM == 6);
    for (ix[5] = 0; ix[5] < sz[5]; ix[5]++)
      for (ix[4] = 0; ix[4] < sz[4]; ix[4]++)
        for (ix[3] = 0; ix[3] < sz[3]; ix[3]++)
          for (ix[2] = 0; ix[2] < sz[2]; ix[2]++)
            for (ix[1] = 0; ix[1] < sz[1]; ix[1]++)
              for (ix[0] = 0; ix[0] < sz[0]; ix[0]++)
                { posA[kpos] = ix_position(A->d, ix, A->bp, A->st);
                  posB[kpos] = ix_position(B->d, ix, B->bp, B->st);
                  posC[kpos] = ix_position(C->d, ix, C->bp, C->st);
                  kpos++;
                }
    assert(kpos == npos);
    if (verbose) { fprintf(stderr, ("  counted %" uint64_u_fmt " tuples\n"), npos); }

    /* Check forward enumeration, Fortran order: */
    kpos = 0;
    bool_t stop;
    stop = ix_enum(check_fw, d, sz, ix_order_F, FALSE, A->bp, A->st, B->bp, B->st, C->bp, C->st);
    assert(! stop);
    assert(kpos == npos);

    /* Check backward enumeration, Fortran order: */
    kpos = npos;
    stop = ix_enum(check_bw, d, sz, ix_order_F, TRUE,  A->bp, A->st, B->bp, B->st, C->bp, C->st);
    assert(! stop);
    assert(kpos == 0);

    /* Enumerate all positions by hand in "L" order: */
    kpos = 0;
    assert(MAXDIM == 6);
    for (ix[0] = 0; ix[0] < sz[0]; ix[0]++)
      for (ix[1] = 0; ix[1] < sz[1]; ix[1]++)
        for (ix[2] = 0; ix[2] < sz[2]; ix[2]++)
          for (ix[3] = 0; ix[3] < sz[3]; ix[3]++)
            for (ix[4] = 0; ix[4] < sz[4]; ix[4]++)
              for (ix[5] = 0; ix[5] < sz[5]; ix[5]++)
                { posA[kpos] = ix_position(A->d, ix, A->bp, A->st);
                  posB[kpos] = ix_position(B->d, ix, B->bp, B->st);
                  posC[kpos] = ix_position(C->d, ix, C->bp, C->st);
                  kpos++;
                }
    assert(kpos == npos);

    /* Check forward enumeration, C/Pascal order: */
    kpos = 0;
    stop = ix_enum(check_fw, d, sz, ix_order_L, FALSE, A->bp, A->st, B->bp, B->st, C->bp, C->st);
    assert(! stop);
    assert(kpos == npos);

    /* Check backward enumeration, C/Pascal order: */
    kpos = npos;
    stop = ix_enum(check_bw, d, sz, ix_order_L, TRUE,  A->bp, A->st, B->bp, B->st, C->bp, C->st);
    assert(! stop);
    assert(kpos == 0);

    /* IMPLEMENTATION OF LOCAL PROCS */

    bool_t check_fw( const ix_index_t ix[], ix_pos_t pA, ix_pos_t pB, ix_pos_t pC )
      { assert(kpos >= 0);
        assert(kpos < npos);
        assert(pA == posA[kpos]);
        assert(pB == posB[kpos]);
        assert(pC == posC[kpos]);
        kpos++;
        return FALSE;
      }

    bool_t check_bw( const ix_index_t ix[], ix_pos_t pA, ix_pos_t pB, ix_pos_t pC )
      { kpos--;
        assert(kpos >= 0);
        assert(kpos < npos);
        assert(pA == posA[kpos]);
        assert(pB == posB[kpos]);
        assert(pC == posC[kpos]);
        return FALSE;
      }
  }

void test_packed(desc_t *A)
  {
    ix_count_t nt = ix_num_tuples(A->d, A->sz);
    fprintf(stderr, ("Checking ix_packed_XXX ops (d = %d nt = %" uint64_u_fmt ")...\n"), A->d, nt);

    /* Get number of indices of {A}: */
    ix_dim_t d = A->d;

    /* Get the size vector, extend with 1's: */
    ix_size_t sz[MAXDIM];
    extend_size_vector(A->d, A->sz, MAXDIM, sz);

    /* Fabricate a *packed* array (can't use {A}, but use its size): */
    ix_pos_t b = 417; /* Arbitrary base postion. */
    desc_t *F = make_packed_desc(d, b, A->sz, ix_order_F);
    desc_t *L = make_packed_desc(d, b, A->sz, ix_order_L);
    bool_t verbose = FALSE;
    if (verbose) 
      { print_desc(stderr, "A = { ", A, " }\n");
        print_desc(stderr, "F = { ", F, " }\n");
        print_desc(stderr, "L = { ", L, " }\n");
      }

    ix_index_t ix[MAXDIM], jx[MAXDIM];

    ix_count_t npos = ix_num_tuples(d, sz); /* Total positions of descriptors. */
    /* Enumerate all tuples by hand in "F" order, check {F} positions: */
    ix_pos_t posExp = b;
    assert(MAXDIM == 6);
    for (ix[5] = 0; ix[5] < sz[5]; ix[5]++)
      for (ix[4] = 0; ix[4] < sz[4]; ix[4]++)
        for (ix[3] = 0; ix[3] < sz[3]; ix[3]++)
          for (ix[2] = 0; ix[2] < sz[2]; ix[2]++)
            for (ix[1] = 0; ix[1] < sz[1]; ix[1]++)
              for (ix[0] = 0; ix[0] < sz[0]; ix[0]++)
                { ix_pos_t posF1 = ix_packed_position(d, ix, F->bp, sz, ix_order_F);
                  ix_pos_t posF2 = ix_position(d, ix, F->bp, F->st);
                  if (verbose)
                    { fprintf(stderr, ("  posExp = %" uint64_u_fmt "\n"), posExp);  
                      fprintf(stderr, ("  posF1 =  %" uint64_u_fmt "\n"), posF1); 
                      fprintf(stderr, ("  posF2 =  %" uint64_u_fmt "\n"), posF2);
                    }
                  assert(posF1 == posExp);
                  assert(posF2 == posF1);
                  ix_packed_indices(d, posExp, F->bp, sz, ix_order_F, jx);
                  for (ix_axis_t i = 0; i < d; i++) { assert(ix[i] == jx[i]); }
                  posExp++;
                }
    assert(posExp == b + npos);

    /* Enumerate all tuples by hand in "L" order, check {L} positions: */
    posExp = b;
    assert(MAXDIM == 6);
    for (ix[0] = 0; ix[0] < sz[0]; ix[0]++)
      for (ix[1] = 0; ix[1] < sz[1]; ix[1]++)
        for (ix[2] = 0; ix[2] < sz[2]; ix[2]++)
          for (ix[3] = 0; ix[3] < sz[3]; ix[3]++)
            for (ix[4] = 0; ix[4] < sz[4]; ix[4]++)
              for (ix[5] = 0; ix[5] < sz[5]; ix[5]++)
                { ix_pos_t posL1 = ix_packed_position(d, ix, L->bp, sz, ix_order_L);
                  ix_pos_t posL2 = ix_position(d, ix, L->bp, L->st);
                  assert(posL1 == posExp);
                  assert(posL2 == posL1);
                  ix_packed_indices(d, posExp, L->bp, sz, ix_order_L, jx);
                  for (ix_axis_t i = 0; i < d; i++) { assert(ix[i] == jx[i]); }
                  posExp++;
                }
    assert(posExp == b + npos);

  }

void test_crop(desc_t *A)
  {
    if (A->d == 0) { return; }
    ix_count_t nt = ix_num_tuples(A->d, A->sz);
    fprintf(stderr, ("Checking ix_crop (d = %d nt = %" uint64_u_fmt ")...\n"), A->d, nt);
    bool_t verbose = FALSE;
    if (verbose) { print_desc(stderr, "A = { ", A, " }\n"); } 
    uint32_t ntrials = (uint32_t)(8*A->d + 1);
    ix_index_t sh[MAXDIM]; /* Shift of cropped array along each axis. */
    ix_size_t sz[MAXDIM]; /* Size of cropped array along each axis. */
    for (uint32_t trial = 0;  trial < ntrials; trial++)
      { /* Initialze {sh} and {sz} for full array: */
        for (ix_axis_t i = 0; i < MAXDIM; i++) { sh[i] = 0; }
        extend_size_vector(A->d, A->sz, MAXDIM, sz);

        /* Start with the standard array: */
        desc_t B = (*A);
        if (verbose) { print_desc(stderr, "B = { ", &B, " }\n"); } 

        /* Perform first cropping: */
        ix_axis_t ax1 = (ix_axis_t)(trial % A->d);
        ix_size_t skip1 = (ix_size_t)((trial & 3)*sz[ax1]/8);
        ix_size_t keep1 = (ix_size_t)((trial/4 & 3)*(sz[ax1] - skip1)/8);
        if (keep1 < 0) { keep1 = 0; }
        sh[ax1] += (int32_t)skip1; sz[ax1] = keep1;
        ix_crop(B.d, B.sz, &B.bp, B.st, ax1, skip1, keep1);
        if (verbose) 
          { fprintf(stderr, ("crop(%d, %" uint64_u_fmt ", %" uint64_u_fmt ")"), ax1, skip1, keep1);
            print_desc(stderr, " = { ", &B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ix_parms_are_valid(B.d, B.sz, B.bp, B.st, FALSE))
          { print_desc(stderr, "A = { ", A, " }\n");
            fprintf(stderr, ("crop(%d, %" uint64_u_fmt ", %" uint64_u_fmt ")"), ax1, skip1, keep1);
            print_desc(stderr, "B = { ", &B, " }\n");
            fprintf(stderr, ("crop(%d, %" uint64_u_fmt ", %" uint64_u_fmt ")"), ax1, skip1, keep1);
            bug("ix_crop returns an invalid desc");
          }
        /* Check whether the size of {B} is correct: */
        check_sizes(B.d, B.sz, sz);

        /* Perform second cropping (possibly on same axis): */ 
        ix_axis_t ax2 = (ix_axis_t)((trial / A->d) % A->d);
        ix_size_t skip2 = (ix_size_t)((trial/2 & 3)*sz[ax2]/8);
        ix_size_t keep2 = (ix_size_t)((trial/5 & 3)*(sz[ax2] - skip2)/8);
        if (keep2 < 0) { keep2 = 0; }
        sh[ax2] += (int32_t)skip2; sz[ax2] = keep2;
        ix_crop(B.d, B.sz, &B.bp, B.st, ax2, skip2, keep2);
        if (verbose) 
          { fprintf(stderr, ("crop(%d, %" uint64_u_fmt ", %" uint64_u_fmt ")"), ax2, skip2, keep2);
            print_desc(stderr, " = { ", &B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ix_parms_are_valid(B.d, B.sz, B.bp, B.st, FALSE))
          { bug("ix_crop returns an invalid desc"); }

        /* Check whether the size of {B} is correct: */
        check_sizes(B.d, B.sz, sz);

        /* Now check coincidence of the two arrays: */

        ix_index_t ixA[MAXDIM], ixB[MAXDIM];
        assert(MAXDIM == 6);
        for (ixB[5] = 0; ixB[5] < sz[5]; ixB[5]++)
          for (ixB[4] = 0; ixB[4] < sz[4]; ixB[4]++)
            for (ixB[3] = 0; ixB[3] < sz[3]; ixB[3]++)
              for (ixB[2] = 0; ixB[2] < sz[2]; ixB[2]++)
                for (ixB[1] = 0; ixB[1] < sz[1]; ixB[1]++)
                  for (ixB[0] = 0; ixB[0] < sz[0]; ixB[0]++)
                    { 
                      ix_assign(A->d, ixA, ixB);
                      ix_shift(A->d, ixA, sh);
                      ix_pos_t posB = ix_position(B.d, ixB, B.bp, B.st);
                      ix_pos_t posA = ix_position(A->d, ixA, A->bp, A->st);
                      if (posB != posA)
                        { bug_sz_sz("crop position error", "posB = ", posB, "posA =", posA); } 
                    }
      }
  }

void test_subsample(desc_t *A) 
  {
    if (A->d == 0) { return; }
    ix_count_t nt = ix_num_tuples(A->d, A->sz);
    fprintf(stderr, ("Checking ix_subsample (d = %d nt = %" uint64_u_fmt ")...\n"), A->d, nt);
    bool_t verbose = FALSE;
    if (verbose) { print_desc(stderr, "A = { ", A, " }\n"); } 
    uint32_t ntrials = (uint32_t)(8*A->d + 1);
    ix_step_t st[MAXDIM]; /* Subsampling step along each axis. */
    ix_size_t sz[MAXDIM]; /* Size of subsampled array along each axis. */
    for (uint32_t trial = 0;  trial < ntrials; trial++)
      { /* Initialze {st} and {sz} for full array: */
        for (ix_axis_t i = 0; i < MAXDIM; i++) { st[i] = 1; }
        extend_size_vector(A->d, A->sz, MAXDIM, sz);
        /* Start with the standard array: */
        desc_t B = (*A);

        /* Perform first subsampling: */
        ix_axis_t ax1 = (ix_axis_t)(trial % A->d);
        ix_size_t step1 = 1 + (trial & 3)*(sz[ax1] - 1)/8;
        ix_size_t keep1 = (sz[ax1] + step1 - 1)/step1;
        st[ax1] = (ix_step_t)(st[ax1] * (int64_t)step1); 
        sz[ax1] = keep1;
        ix_subsample(B.d, B.sz, &(B.bp), B.st, ax1, step1);
        if (verbose) 
          { fprintf(stderr, ("subsample(%d, %" uint64_u_fmt ")"), ax1, step1);
            print_desc(stderr, " = { ", &B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ix_parms_are_valid(B.d, B.sz, B.bp, B.st, FALSE)) 
          { bug("ix_flip returns an invalid desc"); }

        /* Perform second subsampling (possibly on same axis): */ 
        ix_axis_t ax2 = (ix_axis_t)((trial / A->d) % A->d);
        ix_size_t step2 = (ix_size_t)(1 + (trial/2 & 3)*(sz[ax2] - 1)/8);
        ix_size_t keep2 = (ix_size_t)((sz[ax2] + step2 - 1)/step2);
        st[ax2] = (ix_step_t)(st[ax2] * (int64_t)step2); 
        sz[ax2] = keep2;
        ix_subsample(B.d, B.sz, &(B.bp), B.st, ax2, step2);
        if (verbose) 
          { fprintf(stderr, ("subsample(%d, %" uint64_u_fmt ")"), ax2, step2);
            print_desc(stderr, " = { ", &B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ix_parms_are_valid(B.d, B.sz, B.bp, B.st, FALSE)) 
          { bug("ix_subsample returns an invalid desc"); }

        check_sizes(B.d, B.sz, sz); 
        /* Now check coincidence of the two arrays: */
        ix_index_t ixA[MAXDIM], ixB[MAXDIM];
        assert(MAXDIM == 6);
        for (ixB[5] = 0; ixB[5] < sz[5]; ixB[5]++)
          for (ixB[4] = 0; ixB[4] < sz[4]; ixB[4]++)
            for (ixB[3] = 0; ixB[3] < sz[3]; ixB[3]++)
              for (ixB[2] = 0; ixB[2] < sz[2]; ixB[2]++)
                for (ixB[1] = 0; ixB[1] < sz[1]; ixB[1]++)
                  for (ixB[0] = 0; ixB[0] < sz[0]; ixB[0]++)
                    { 
                      for (ix_axis_t ax = 0; ax < A->d; ax++) { ixA[ax] = ixB[ax]*st[ax]; }
                      ix_pos_t posB = ix_position(B.d, ixB, B.bp, B.st);
                      ix_pos_t posA = ix_position(A->d, ixA, A->bp, A->st);
                      if (posB != posA)
                        { bug_sz_sz("subsample position error", "posB =", posB, "posA =", posA); } 
                    }
      }
  }

void test_slice(desc_t *A)
  {
    if (A->d == 0) { return; }
    ix_count_t nt = ix_num_tuples(A->d, A->sz);
    if (nt == 0) { return; }
    fprintf(stderr, ("Checking ix_slice (d = %d nt = %" uint64_u_fmt ")...\n"), A->d, nt);
    bool_t verbose = TRUE;
    uint32_t ntrials = (uint32_t)(8*A->d + 1);
    /* Expected slice attributes and mapping between slice and original: */
    uint64_t toss = 417*A->d;
    for (uint32_t trial = 0; trial < ntrials; trial++)
      { /* Start current slice with the given array: */
        desc_t B = (*A);
        if (verbose) { print_desc(stderr, "\nnew trial B = A = { ", A, " }\n"); } 

        /* Mapping to current slice to original array: */
        ix_axis_t axO[MAXDIM];  /* Axis of original for each unset axis of slice, or {A->d} if slice axis is trivial. */
        ix_index_t ixO[MAXDIM]; /* Index value for each axis of original that was set, or {-1}. */
        for (ix_axis_t i = 0; i < MAXDIM; i++) { axO[i] = i; ixO[i] = -1; }

        for (uint32_t pass = 0; pass < 2; pass++)
          { 
            /* Here axis {k} of current slice is axis {axO[k]} of original, or trivial if {axO[k]==A->d}. */
            /* Here if {ixO[k]} is non-negative, then axis {k} of original was set to {ixO[k]}. */

            /* Choose some axes to fix, and their index values: */
            ix_dim_t nx1 = 0;            /* Number of axes to fix. */
            ix_axis_t ax1[MAXDIM];  /* Which axes to fix. */
            ix_index_t ix1[MAXDIM]; /* Index values of those axes. */
            ix_axis_t *axp = (((toss/9) % 3) == 0 ? NULL : ax1);
            uint32_t nxp = (uint32_t)((toss/11) % (B.d + 1)); /* Number of axes to fix if {axp = NULL}. */
            { ix_axis_t j = 0;      /* Axis of next slice. */
              for (ix_axis_t i = 0; i < B.d; i++)
                { /* Decide whether to fix axis {i}: */
                  bool_t fixit = (axp == NULL ? nx1 < nxp : (((toss/3) % 3) == 0));
                  if (fixit)
                    { /* Fix axis {i} of current slice: */
                      ax1[nx1] = i; 
                      if (axO[i] < A->d)
                        { /* Axis was not set before. */
                          ix_size_t szi = A->sz[axO[i]];
                          assert(szi > 0);
                          ix1[nx1] = (ix_index_t)((toss/5) % szi);
                          /* Remember that this index was set: */
                          ixO[axO[i]] = ix1[nx1];
                        }
                      else
                        { /* This axis of the original was already set: */
                          ix1[nx1] = 0;
                        }
                      nx1++;
                    }
                  else
                    { /* Keep axis {i} of current slice as axis {j} of next slice: */
                      axO[j] = axO[i];
                      j++;
                    }
                  toss = 27*toss + 3*pass + 7*(uint64_t)i;
                }
              assert(j == B.d - nx1);
              /* Mark axes {j..bp.d} as trivialized:*/
              while(j < B.d) { axO[j] = A->d; j++; }
            }
            if (axp == NULL) { assert(nx1 == nxp); }
            
            if (verbose)
              { fprintf(stderr, "  axis:index =");
                for (uint32_t i = 0;  i < nx1; i++)
                  { fprintf(stderr, " %d:%ld", (axp == NULL ? i : axp[i]), ix1[i]); }
                fprintf(stderr, "\n");
              }

            /* Perform the slice op on {B}: */
            ix_slice(B.d, B.sz, &(B.bp), B.st, nx1, axp, ix1);
            if (verbose) 
              { fprintf(stderr, "pass %d slice(%d, %s{", pass, nx1, (axp == NULL ? "[NULL]" : ""));
                for (ix_index_t k = 0; k < nx1; k++) 
                  { fprintf(stderr, (" [%d]=%" int64_d_fmt ""), ax1[k], ix1[k]); }
                fprintf(stderr, " })");
                print_desc(stderr, "  B = { ", &B, " }\n");
                fprintf(stderr, "expected correspondence: B[");
                char *del = "";
                for (ix_axis_t k = 0; k < B.d; k++)
                  { if (axO[k] >= A->d)
                      { fprintf(stderr, "%s0", del); } 
                    else
                      { fprintf(stderr, "%si%d", del, axO[k]); }
                    del = ",";
                  }
                fprintf(stderr, "] = A[");
                del = "";
                for (ix_axis_t i = 0; i < A->d; i++) 
                  { if (ixO[i] < 0)
                      { fprintf(stderr, "%si%d", del, i); }
                    else
                      { fprintf(stderr, ("%s%" int64_d_fmt ""), del, ixO[i]); }
                    del = ","; 
                  }
                fprintf(stderr, "]\n");
              } 
            /* Check validity of descriptor: */
            if (! ix_parms_are_valid(B.d, B.sz, B.bp, B.st, FALSE))
              { bug("ix_slice returns an invalid desc"); }
          }

        /* Compute expected {sz,st} of slice, extended with trivial to {MAXDIM}: */
        ix_step_t st[MAXDIM]; /* Step of slice along each axis. */
        ix_size_t sz[MAXDIM]; /* Size of slice along each axis. */
        for (ix_axis_t j = 0; j < MAXDIM; j++) 
          { if ((j < B.d) && (axO[j] < A->d))
              { ix_axis_t i = axO[j]; 
                sz[j] = A->sz[i]; 
                st[j] = A->st[i];
              }
            else
              { sz[j] = 1; st[j] = 0; }
          }

        /* Check whether the size and increments of {B} is correct: */
        { if (B.d != A->d) { bug("B.d = %d  A->d = %d", B.d, A->d); }
          for (ix_axis_t i = 0; i < A->d; i++)
            { if (B.sz[i] != sz[i])
                { bug(("B.sz[%d] = %" uint64_u_fmt "  sz[%d] = %" uint64_u_fmt ""), i, B.sz[i], i, sz[i]); }
              if (B.st[i] != st[i])
                { bug(("B.st[%d] = %" uint64_u_fmt "  A.st[%d] = %" uint64_u_fmt ""), i, B.st[i], i, st[i]); }
            }
        }

        /* Now check coincidence of the two arrays: */
        ix_index_t ixA[MAXDIM], ixB[MAXDIM];
        assert(MAXDIM == 6);

        for (ixB[5] = 0; ixB[5] < sz[5]; ixB[5]++)
          for (ixB[4] = 0; ixB[4] < sz[4]; ixB[4]++)
            for (ixB[3] = 0; ixB[3] < sz[3]; ixB[3]++)
              for (ixB[2] = 0; ixB[2] < sz[2]; ixB[2]++)
                for (ixB[1] = 0; ixB[1] < sz[1]; ixB[1]++)
                  for (ixB[0] = 0; ixB[0] < sz[0]; ixB[0]++)
                    { 
                      for (ix_axis_t i = 0; i < A->d; i++) 
                       { if (ixO[i] >= 0) { ixA[i] = ixO[i]; } }
                      for (ix_axis_t j = 0; j < B.d; j++) 
                        { ix_index_t i = axO[j];
                          if (i < A->d) { assert(ixO[i] < 0); ixA[i] = ixB[j]; }
                        }
                      ix_pos_t posB = ix_position(B.d, ixB, B.bp, B.st);
                      ix_pos_t posA = ix_position(A->d, ixA, A->bp, A->st);
                      if (posB != posA)
                        { bug_sz_sz("slice position error", "posB =", posB, "posA =", posA); } 
                    }

        if (verbose) { fprintf(stderr, "\n"); } 
      }
  }

void test_swap_indices(desc_t *A)
  {
    if (A->d == 0) { return; }
    ix_count_t nt = ix_num_tuples(A->d, A->sz);
    fprintf(stderr, ("Checking ix_swap_indices (d = %d nt = %" uint64_u_fmt ")...\n"), A->d, nt);
    bool_t verbose = FALSE;
    if (verbose) { print_desc(stderr, "A = { ", A, " }\n"); } 
    uint32_t ntrials = (uint32_t)(8*A->d + 1);
    ix_axis_t tr[MAXDIM]; /* Axis {i} of {B} is axis {tr[i]} of {A}. */
    for (uint32_t trial = 0;  trial < ntrials; trial++)
      { /* Initialze {tr} for full array: */
        for (ix_axis_t i = 0; i < MAXDIM; i++) { tr[i] = i; }
        /* Start with the standard array: */
        desc_t B = (*A);

        for (uint32_t pass = 0;  pass < 2; pass++)
          { /* Choose two axes {ax1,bx1} and the axis count {nx1}: */
            ix_axis_t ax1 = (ix_axis_t)(trial % A->d);
            ix_axis_t bx1 = (ix_axis_t)((trial/2) % A->d);
            ix_dim_t dx1 = (ix_dim_t)(ax1 == bx1 ? A->d - ax1 : (ax1 > bx1 ? ax1 - bx1 : bx1 - ax1));
            assert(dx1 > 0);
            ix_dim_t nx1 = (ix_dim_t)((trial/3) % dx1);
            assert(ax1 + nx1 <= A->d);
            assert(bx1 + nx1 <= A->d);
            assert((ax1 == bx1) || (ax1 + nx1 <= bx1) || (bx1 + nx1 <= ax1));
            if (ax1 != bx1)
              { /* Perform a tranpose on {tr}: */
                for (ix_axis_t k = 0; k < nx1; k++) 
                  { ix_axis_t tmp = tr[ax1+k]; tr[ax1+k] = tr[bx1+k]; tr[bx1+k] = tmp; }
              }
            /* Perform swap on {B}: */
            ix_swap_indices(B.d, B.sz, &(B.bp), B.st, ax1, bx1, nx1);
            if (verbose) 
              { fprintf(stderr, "swap_indices(%d, %d, %d)", ax1, bx1, nx1);
                print_desc(stderr, " = { ", &B, " }\n");
              } 
            /* Check validity of descriptor: */
            if (! ix_parms_are_valid(B.d, B.sz, B.bp, B.st, FALSE))
              { bug("ix_swap_indices returns an invalid desc"); }
          }

        /* Check whether the size and increments of {B} is correct: */
        for (ix_axis_t i = 0; i < A->d; i++)
          { if (B.sz[i] != A->sz[tr[i]])
              { bug(("B.sz[%d] = %" uint64_u_fmt "  A.sz[%d] = %" uint64_u_fmt ""), i, B.sz[i], tr[i], A->sz[tr[i]]); }
            if (B.st[0] != A->st[tr[0]])
              { bug(("B.st[%d] = %" uint64_u_fmt "  A.st[%d] = %" uint64_u_fmt ""), i, B.st[i], tr[i], A->st[tr[i]]); }
          }
        /* Get {A}'s size vector, extend with 1's: */
        ix_size_t sz[MAXDIM];
        extend_size_vector(A->d, A->sz, MAXDIM, sz);

        /* Now check coincidence of the two arrays: */
        ix_index_t ixA[MAXDIM], ixB[MAXDIM];
        assert(MAXDIM == 6);

        for (ixA[5] = 0; ixA[5] < sz[5]; ixA[5]++)
          for (ixA[4] = 0; ixA[4] < sz[4]; ixA[4]++)
            for (ixA[3] = 0; ixA[3] < sz[3]; ixA[3]++)
              for (ixA[2] = 0; ixA[2] < sz[2]; ixA[2]++)
                for (ixA[1] = 0; ixA[1] < sz[1]; ixA[1]++)
                  for (ixA[0] = 0; ixA[0] < sz[0]; ixA[0]++)
                    { 
                      for (ix_axis_t ax = 0; ax < A->d; ax++) { ixB[ax] = ixA[tr[ax]]; }
                      ix_pos_t posB = ix_position(B.d, ixB, B.bp, B.st);
                      ix_pos_t posA = ix_position(A->d, ixA, A->bp, A->st);
                      if (posB != posA)
                        { bug_sz_sz("swap_indices position error", "posB =", posB, "posA =", posA); } 
                    }
      }
  }

void test_flip_indices(desc_t *A)
  {
    ix_count_t nt = ix_num_tuples(A->d, A->sz);
    fprintf(stderr, ("NOT Checking ix_flip_indices (d = %d nt = %" uint64_u_fmt ")...\n"), A->d, nt);
  }

void test_flip(desc_t *A)
  {
    if (A->d == 0) { return; }
    ix_count_t nt = ix_num_tuples(A->d, A->sz);
    fprintf(stderr, ("Checking ix_flip (d = %d nt = %" uint64_u_fmt ")...\n"), A->d, nt);
    bool_t verbose = FALSE;
    if (verbose) { print_desc(stderr, "A = { ", A, " }\n"); } 
    uint32_t ntrials = (uint32_t)(8*A->d + 1);
    bool_t fp[MAXDIM]; /* Tells whether each axis was flipped or not. */
    for (uint32_t trial = 0;  trial < ntrials; trial++)
      { /* Initialze {fp} for unflipped array: */
        for (ix_axis_t i = 0; i < MAXDIM; i++) { fp[i] = FALSE; }
        /* Start with the standard array: */
        desc_t B = (*A);

        /* Perform first flip: */
        ix_axis_t ax1 = (ix_axis_t)(trial % A->d);
        fp[ax1] = ! fp[ax1];
        ix_flip(B.d, B.sz, &(B.bp), B.st, ax1);
        if (verbose) 
          { fprintf(stderr, "flip(%d)", ax1);
            print_desc(stderr, " = { ", &B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ix_parms_are_valid(B.d, B.sz, B.bp, B.st, FALSE))
          { bug("ix_flip returns an invalid desc"); }

        /* Perform second flip (possibly on same axis): */ 
        ix_axis_t ax2 = (ix_axis_t)((trial / A->d) % A->d);
        fp[ax2] = ! fp[ax2];
        ix_flip(B.d, B.sz, &(B.bp), B.st, ax2);
        if (verbose) 
          { fprintf(stderr, "flip(%d)", ax2);
            print_desc(stderr, " = { ", &B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ix_parms_are_valid(B.d, B.sz, B.bp, B.st, FALSE))
          { bug("ix_flip returns an invalid desc"); }

        /* Check whether the size of {B} is correct: */
        check_sizes(B.d, B.sz, A->sz);

        /* Check whether the increments of {B} are correct: */
        for (ix_axis_t i = 0; i < A->d; i++)
          { if (B.st[i] != A->st[i]*(fp[i]?-1:+1))
              { bug(("step[%d] = %" uint64_u_fmt ""), i, B.st[i]); }
          }

        /* Get {B}'s size vector, extend with 1's: */
        ix_size_t sz[MAXDIM];
        extend_size_vector(B.d, B.sz, MAXDIM, sz);

        /* Now check coincidence of the two arrays: */
        ix_index_t ixA[MAXDIM], ixB[MAXDIM];
        assert(MAXDIM == 6);
        for (ixB[5] = 0; ixB[5] < sz[5]; ixB[5]++)
          for (ixB[4] = 0; ixB[4] < sz[4]; ixB[4]++)
            for (ixB[3] = 0; ixB[3] < sz[3]; ixB[3]++)
              for (ixB[2] = 0; ixB[2] < sz[2]; ixB[2]++)
                for (ixB[1] = 0; ixB[1] < sz[1]; ixB[1]++)
                  for (ixB[0] = 0; ixB[0] < sz[0]; ixB[0]++)
                    { 
                      for (ix_axis_t ax = 0; ax < A->d; ax++) 
                        { ixA[ax] = (fp[ax] ? (ix_index_t)(A->sz[ax]) - 1 - ixB[ax] : ixB[ax]); }
                      ix_pos_t posB = ix_position(B.d, ixB, B.bp, B.st);
                      ix_pos_t posA = ix_position(A->d, ixA, A->bp, A->st);
                      if (posB != posA)
                        { bug_sz_sz("flip position error", "posB =", posB, "posA =", posA); } 
                    }
      }
  }

void test_diagonal(desc_t *A)
  {
    if (A->d < 2) { return; }
    ix_count_t nt = ix_num_tuples(A->d, A->sz);
    fprintf(stderr, ("Checking ix_diagonal (d = %d nt = %" uint64_u_fmt ")...\n"), A->d, nt);
    bool_t verbose = FALSE;
    if (verbose) { print_desc(stderr, "A = { ", A, " }\n"); } 
    uint32_t ntrials = (uint32_t)(8*A->d + 1);
    ix_size_t sz[MAXDIM]; /* Size of diagonalized array along each axis. */
    for (uint32_t trial = 0;  trial < ntrials; trial++)
      { /* Initialze {sz} for full array: */
        extend_size_vector(A->d, A->sz, MAXDIM, sz);
        /* Start with the standard array: */
        desc_t B = (*A);

        /* Choose any two axes: */
        ix_axis_t ax1 = (ix_axis_t)(trial % A->d);
        ix_axis_t bx1 = (ix_axis_t)((trial/A->d) % A->d);
        if (bx1 == ax1) { bx1 = (ix_axis_t)((ax1 + 1) % A->d); }
        /* Make it square: */
        ix_size_t szmin = (sz[ax1] < sz[bx1] ? sz[ax1] : sz[bx1]);
        ix_crop(B.d, B.sz, &(B.bp), B.st, ax1, 0, szmin); sz[ax1] = szmin;
        ix_crop(B.d, B.sz, &(B.bp), B.st, bx1, 0, szmin); sz[bx1] = szmin;
        /* Diagonalize array: */
        ix_diagonal(B.d, B.sz, &(B.bp), B.st, ax1, bx1);
        /* Compute expected sizes of {B}: */
        if (sz[bx1] < sz[ax1]) { sz[ax1] = sz[bx1]; }
        if (bx1 != ax1) { sz[bx1] = 1; }
        if (verbose) 
          { fprintf(stderr, "diagonal(%d, %d)", ax1, bx1);
            print_desc(stderr, " = { ", &B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ix_parms_are_valid(B.d, B.sz, B.bp, B.st, FALSE))
          { bug("ix_diagonal returns an invalid desc"); }

        /* Check whether the size of {B} is correct: */
        check_sizes(B.d, B.sz, sz);
        /* Now check coincidence of the two arrays: */
        ix_index_t ixA[MAXDIM], ixB[MAXDIM];
        assert(MAXDIM == 6);
        for (ixB[5] = 0; ixB[5] < sz[5]; ixB[5]++)
          for (ixB[4] = 0; ixB[4] < sz[4]; ixB[4]++)
            for (ixB[3] = 0; ixB[3] < sz[3]; ixB[3]++)
              for (ixB[2] = 0; ixB[2] < sz[2]; ixB[2]++)
                for (ixB[1] = 0; ixB[1] < sz[1]; ixB[1]++)
                  for (ixB[0] = 0; ixB[0] < sz[0]; ixB[0]++)
                    { 
                      ix_assign(A->d, ixA, ixB);
                      ixA[bx1] = ixA[ax1];
                      ix_pos_t posB = ix_position(B.d, ixB, B.bp, B.st);
                      ix_pos_t posA = ix_position(A->d, ixA, A->bp, A->st);
                      if (posB != posA)
                        { bug_sz_sz("diagonal position error", "posB =", posB, "posA =", posA); } 
                    }
      }
  }

void test_chop(desc_t *A)
  {
    if (A->d < 2) { return; }
    ix_count_t nt = ix_num_tuples(A->d, A->sz);
    fprintf(stderr, ("Checking ix_chop (d = %d nt = %" uint64_u_fmt ")...\n"), A->d, nt);
    bool_t verbose = FALSE;
    if (verbose) { print_desc(stderr, "A = { ", A, " }\n"); } 
    uint32_t ntrials = (uint32_t)(8*A->d + 1);
    ix_size_t sz[MAXDIM]; /* Size of chopped array along each axis. */
    for (uint32_t trial = 0;  trial < ntrials; trial++)
      { /* Initialze {sz} for full array: */
        extend_size_vector(A->d, A->sz, MAXDIM, sz);
        /* Start with the standard array: */
        desc_t B = (*A);

        /* Choose two DISTINCT axes: */
        ix_axis_t ax1 = (ix_axis_t)(trial % A->d);
        ix_axis_t bx1 = (ix_axis_t)((trial/A->d) % (A->d - 1));
        if (bx1 >= ax1) { bx1++; }
        /* Choose a nonzero chunk size: */
        ix_size_t chunksz1 = 1 + (trial/3 & 3)*((A->sz[ax1] <= 1 ? 1 : A->sz[ax1]) + 7)/8;
        /* Ensure integral number of chunks: */
        ix_size_t szround = (sz[ax1]/chunksz1)*chunksz1;
        ix_crop(B.d, B.sz, &(B.bp), B.st, ax1, 0, szround); sz[ax1] = szround;
        /* Ensure that the argument of {chop} is trivial along axis {bx1}: */
        if (sz[bx1] == 0) { /* Skip this trial: */ continue; }
        ix_crop(B.d, B.sz, &(B.bp), B.st, bx1, 0, 1); sz[bx1] = 1;
        /* Chop array: */
        ix_chop(B.d, B.sz, &(B.bp), B.st, ax1, chunksz1, bx1);
        /* Compute expected sizes of {B}: */
        sz[bx1] = sz[ax1]/chunksz1; sz[ax1] = chunksz1;
        if (verbose) 
          { fprintf(stderr, ("chop(%d, %" int64_d_fmt ", %d)"), ax1, chunksz1, bx1);
            print_desc(stderr, " = { ", &B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ix_parms_are_valid(B.d, B.sz, B.bp, B.st, FALSE))
          { bug("ix_chop returns an invalid desc"); }

        /* Check whether the size of {B} is correct: */
        check_sizes(B.d, B.sz, sz);
        /* Now check coincidence of the two arrays: */
        ix_index_t ixA[MAXDIM], ixB[MAXDIM];
        assert(MAXDIM == 6);
        for (ixB[5] = 0; ixB[5] < sz[5]; ixB[5]++)
          for (ixB[4] = 0; ixB[4] < sz[4]; ixB[4]++)
            for (ixB[3] = 0; ixB[3] < sz[3]; ixB[3]++)
              for (ixB[2] = 0; ixB[2] < sz[2]; ixB[2]++)
                for (ixB[1] = 0; ixB[1] < sz[1]; ixB[1]++)
                  for (ixB[0] = 0; ixB[0] < sz[0]; ixB[0]++)
                    { 
                      ix_assign(A->d, ixA, ixB);
                      ixA[ax1] = ixA[ax1] + (ix_index_t)chunksz1 * ixB[bx1];
                      ixA[bx1] = 0;
                      ix_pos_t posB = ix_position(B.d, ixB, B.bp, B.st);
                      ix_pos_t posA = ix_position(A->d, ixA, A->bp, A->st);
                      if (posB != posA)
                        { bug_sz_sz("chop position error", "posB =", posB, "posA =", posA); } 
                    }
      }
  }

void test_replicate(desc_t *A)
  {
    if (A->d == 0) { return; }
    ix_count_t nt = ix_num_tuples(A->d, A->sz);
    fprintf(stderr, ("Checking ix_replicate (d = %d nt = %" uint64_u_fmt ")...\n"), A->d, nt);
    bool_t verbose = FALSE;
    if (verbose) { print_desc(stderr, "A = { ", A, " }\n"); } 
    uint32_t ntrials = (uint32_t)(8*A->d + 1);
    ix_size_t sz[MAXDIM]; /* Size of replicated array along each axis. */
    for (uint32_t trial = 0; trial < ntrials; trial++)
      { /* Initialze {sz} for full array: */
        extend_size_vector(A->d, A->sz, MAXDIM, sz);
        /* Start with the standard array: */
        desc_t B = (*A);

        /* Choose axis for first replication: */
        ix_axis_t ax1 = (ix_axis_t)(trial % A->d);
        if (sz[ax1] == 0) { continue; }
        /* Choose a positive replication factor: */
        ix_size_t rep1 = 1 + (trial/3 & 3)*A->d;
        /* Make trivial along axis {ax1}: */
        ix_size_t skip1 = (B.sz[ax1]-1)/2;
        ix_crop(B.d, B.sz, &(B.bp), B.st, ax1, skip1, 1); sz[ax1] = 1;
        /* Replicate along axis {ax1}: */
        ix_replicate(B.d, B.sz, &(B.bp), B.st, ax1, rep1);
        sz[ax1] = rep1;
        if (verbose) 
          { fprintf(stderr, ("chop(%d,%" uint64_u_fmt ",1)+replicate(%d,%" int64_d_fmt ")"), ax1, skip1, ax1, rep1);
            print_desc(stderr, " = { ", &B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ix_parms_are_valid(B.d, B.sz, B.bp, B.st, FALSE))
          { bug("ix_replicate returns an invalid desc"); }

        /* Choose axis for second replication (may be {ax1}): */
        ix_axis_t ax2 = (ix_axis_t)((trial/A->d) % A->d);
        if (sz[ax2] == 0) { continue; }
        /* Choose a replication factor: */
        ix_size_t rep2 = 1 + (trial/8 & 3)*A->d;
        /* Make trivial along axis {ax2}: */
        ix_size_t skip2 = (B.sz[ax2]-1)/2;
        ix_crop(B.d, B.sz, &(B.bp), B.st, ax2, skip2, 1); sz[ax2] = 1;
        /* Replicate along axis {ax2}: */
        ix_replicate(B.d, B.sz, &(B.bp), B.st, ax2, rep2);
        sz[ax2] = rep2;
        if (verbose) 
          { fprintf(stderr, ("chop(%d,%" uint64_u_fmt ",1)+replicate(%d,%" int64_d_fmt ")"), ax2, skip2, ax2, rep2);
            print_desc(stderr, " = { ", &B, " }\n");
          } 
        /* Check validity of descriptor: */
        if (! ix_parms_are_valid(B.d, B.sz, B.bp, B.st, FALSE))
          { bug("ix_replicate returns an invalid desc"); }

        /* Check whether the size of {B} is correct: */
        check_sizes(B.d, B.sz, sz);
        /* Now check coincidence of the two arrays: */
        ix_index_t ixA[MAXDIM], ixB[MAXDIM];
        assert(MAXDIM == 6);
        for (ixB[5] = 0; ixB[5] < sz[5]; ixB[5]++)
          for (ixB[4] = 0; ixB[4] < sz[4]; ixB[4]++)
            for (ixB[3] = 0; ixB[3] < sz[3]; ixB[3]++)
              for (ixB[2] = 0; ixB[2] < sz[2]; ixB[2]++)
                for (ixB[1] = 0; ixB[1] < sz[1]; ixB[1]++)
                  for (ixB[0] = 0; ixB[0] < sz[0]; ixB[0]++)
                    { 
                      ix_assign(A->d, ixA, ixB);
                      ixA[ax1] = (ix_index_t)skip1;
                      /* The second crop changes {b} only if {ax1 != ax2}: */
                      if (ax2 != ax1) { ixA[ax2] = (ix_index_t)skip2; }
                      ix_pos_t posB = ix_position(B.d, ixB, B.bp, B.st);
                      ix_pos_t posA = ix_position(A->d, ixA, A->bp, A->st);
                      if (posB != posA)
                        { bug_sz_sz("replicate position error", "posB =", posB, "posA =", posA); } 
                    }
      }
  }

void test_reduce(void)
  {
    fprintf(stderr, "Checking ix_reduce...\n");
    bool_t verbose = FALSE;
    /* Precomputed examples: */
    ix_size_t N = 6;
    ix_index_t iMin = -3;
    ix_index_t iMax = 15;
    char *tseq[ix_reduce_mode_LAST + 1];
    tseq[ix_reduce_mode_SINGLE] = "*,*,*,0,1,2,3,4,5,*,*,*,*,*,*,*,*,*,*,";
    tseq[ix_reduce_mode_EXTEND] = "0,0,0,0,1,2,3,4,5,5,5,5,5,5,5,5,5,5,5,";
    tseq[ix_reduce_mode_REPEAT] = "3,4,5,0,1,2,3,4,5,0,1,2,3,4,5,0,1,2,3,";
    tseq[ix_reduce_mode_MIRROR] = "2,1,0,0,1,2,3,4,5,5,4,3,2,1,0,0,1,2,3,";
    tseq[ix_reduce_mode_PXMIRR] = "3,2,1,0,1,2,3,4,5,4,3,2,1,0,1,2,3,4,5,";

    if (verbose) { fprintf(stderr, ("N = %" int64_d_fmt "\n"), N); } 
    for (ix_reduce_mode_t red = ix_reduce_mode_FIRST; red <= ix_reduce_mode_LAST; red++)
      { if (verbose) { fprintf(stderr, "red = %d (%c)\n", red, "SERMP"[red]); } 
        for (ix_index_t i = iMin; i <= iMax; i++)
          { char cr = tseq[red][2*(i-iMin)];
            ix_index_t jExp = (cr == '*' ? -1 : cr - '0');
            ix_index_t j = ix_reduce(i, N, red);
            if (verbose) { fprintf(stderr, ("%3" int64_d_fmt " --> %3" int64_d_fmt "  (should be %3" int64_d_fmt ")\n"), i, j, jExp); } 
            demand(((red == ix_reduce_mode_SINGLE) && (j == -1)) || ((j >= 0) && (j < N)), "invalid {ix_reduce} result");
            demand(j == jExp, "error in {ix_reduce}");
          }
      }

  }

void extend_size_vector(ix_dim_t dA, ix_size_t szA[], ix_dim_t dB, ix_size_t szB[])
  {
    for (ix_axis_t i = 0; i < dB; i++) { szB[i] = (i < dA ? szA[i] : 1); }
  }

void check_max_indices(ix_dim_t d, ix_index_t ixa[], ix_index_t ixb[])
  {
    for (ix_axis_t i = 0; i < d; i++)
      { if (ixa[i] != ixb[i]) 
          { bug(("(ixa[%d] = %" int64_d_fmt " ixb[%d] = %" int64_d_fmt ""), i, ixa[i], i, ixb[i]); }
      }
  }

void check_max_indices_sz(ix_dim_t d, ix_index_t ixa[], ix_size_t szb[])
  {
    for (ix_axis_t i = 0; i < d; i++)
      { ix_index_t ixbi = (ix_index_t)(szb[i]-1);
        if (ixa[i] != ixbi) 
          { bug(("(ixa[%d]-1) = %" int64_d_fmt " ixb[%d] = %" int64_d_fmt ""), i, ixa[i], i, ixbi); }
      }
  }

void check_sizes(ix_dim_t d, ix_size_t sza[], ix_size_t szb[])
  {
    for (ix_axis_t i = 0; i < d; i++)
      { if (sza[i] != szb[i]) 
          { bug(("sza[%d] = %" uint64_u_fmt " szb[%d] = %" uint64_u_fmt ""), i, sza[i], i, szb[i]); }
      }
  }

void print_desc(FILE *wr, char *pf, desc_t *A, char *sf)
  { fprintf(wr, "%s", pf); 
    for (ix_axis_t i = 0; i < A->d; i++) 
      { fprintf(wr, ("%s%" uint64_u_fmt "(×%" int64_d_fmt ")"), " "+(i==0), A->sz[i], A->st[i]); }
    fprintf(wr, 
      ("  b = %" uint64_u_fmt "  tuples = %" uint64_u_fmt "  positions = %" uint64_u_fmt ""), 
      A->bp, 
      ix_num_tuples(A->d,A->sz), 
      ix_num_positions(A->d,A->sz,A->st)
    );
    fprintf(wr, "%s", sf);
  }

void print_indices(FILE *wr, char *pf, ix_dim_t d, ix_index_t ix[], char *sf)
  { fprintf(wr, "%s", pf); 
    for (ix_axis_t i = 0; i < d; i++)
      { fprintf(wr, (" %2" int64_d_fmt ""), ix[i]); }
    fprintf(wr, "%s", sf);
  }
