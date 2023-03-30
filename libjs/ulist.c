/* See ulist.h */
/* Last edited on 2023-03-18 11:31:38 by stolfi */

/* !!! Test and debug throughly !!! */
/* !!! There may be performance problems (unnecessary collisions?) !!! */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <vec.h>
#include <bool.h>
#include <ref.h>
#include <jsrandom.h>
#include <jsmath.h>
#include <affirm.h>
#include <jswsize.h>

#include <ulist.h>
#include <ulist_impl.h>
  
ulist_hash_size_t ulist_choose_hash_size(ulist_count_t n);
/* Chooses a suitable hash table size for a list with capacity {n}. */

/* PROCEDURES THAT DO NOT DEPEND ON HASHING */

vec_typeimpl(ulist_item_vec_t, ulist_item_vec, ulist_item_t);
vec_typeimpl(ulist_hash_val_vec_t, ulist_hash_val_vec, ulist_hash_val_t);
vec_typeimpl(ulist_index_vec_t, ulist_index_vec, ulist_index_t);

ulist_t *ulist_new(ulist_count_t n)
  { demand(n <= ulist_item_count_MAX, "requested capacity is too large");
    ulist_t *S = (ulist_t *)notnull(malloc(sizeof(ulist_t)), "out of mem");
    S->it = ulist_item_vec_new(n);
    S->ct = 0;
    S->hx = ulist_hash_val_vec_new(n);
    ulist_hash_size_t nh = ulist_choose_hash_size(n); 
    S->ix = ulist_index_vec_new(nh);
    S->eq = &ulist_uint64_eq;
    S->hash = &ulist_uint64_hash;
    #if (ulist_DEBUG)
    ulist_stats_clear(S);
    #endif /* (ulist_DEBUG) */
    return S;
  }
  
ulist_hash_size_t ulist_choose_hash_size(ulist_count_t n)  
  { if (n == 0) { /* No need for a hash table: */ return 0; }
    double dnh = 2.0*((double)n); /* Ideal size of hash table. */
    demand(dnh <= (double)ulist_hash_size_MAX, "hash table would be too big");
    ulist_hash_size_t nh = (ulist_hash_size_t)floor(dnh);
    /* If {n} is large enough, we had better round it up to avoid small divisors: */
    uint32_t pmax = 13; /* Largest prime we plan to avoid. */
    uint32_t m = 1; /* Product of all small primes we intend to avoid. */
    uint32_t r = 0; /* Desired remainder modulo {m}. */
    uint32_t p = 1; /* Largest prime in {m}, or 1: */
    if (ulist_DEBUG) { fprintf(stderr, "  n = %u  ideal nh = %d", n, nh); }
    while (TRUE)
      { /* Find the next prime {p}: */
        do { p++; } while (gcd(m, p) != 1); 
        if (p > pmax) { /* No more primes to avoid: */ break; }
        if (m > 3*nh/2/p) { /* Avoiding {p} would require too many entries: */ break; }
        /* Adjust {r} so that it is about {m/2} and {gcd(r,m*p) = 1}: */
        uint32_t q = (uint32_t)(floor(0.61803398874989484820*p));
        r = q*m + r;
        /* Include {p} in {m}: */
        m = p*m;
        if (ulist_DEBUG) { fprintf(stderr, " [k*%d + %d]", m, r); }
      }
    assert(m < ulist_hash_size_MAX); /* Paranoia... */
    /* Round {nh} down to a multiple {nhs} of {m}: */
    ulist_hash_size_t nhs = (nh / m) * m;
    /* Ensure that {nhs + r} is at least {nh}: */
    if (r < nh - nhs) { r += m; }
    /* Ensure that {nhs + r} is not too big: */
    while (r > ulist_hash_size_MAX - nhs)
      { assert(nhs >= m);
        nhs -= m;
      }
    /* Phew: */
    nh = nhs + r;
    if (ulist_DEBUG) { fprintf(stderr, "  nh = %u\n", nh); }
    return nh;
  }

void ulist_free(ulist_t *S)
  { free(S->it.e); 
    free(S->ix.e); 
    free(S->hx.e); 
    free(S); 
  }

void ulist_clear(ulist_t *S)
  { S->ct = 0; }

ulist_eq_func_t *ulist_get_eq(ulist_t *S)
  { return S->eq; }

void ulist_install_eq(ulist_t *S, ulist_eq_func_t *eq)
  { if (eq == NULL) { /* Provide default: */ eq = &ulist_uint64_eq; }
    if (eq != S->eq)
      { demand((S->ct <= 1), "cannot change {eq} of non-trivial list"); 
        S->eq = eq;
      }
  }

ulist_hash_func_t *ulist_get_hash(ulist_t *S)
  { return S->hash; }

void ulist_install_hash(ulist_t *S, ulist_hash_func_t *hash)
  { if (hash == NULL) { /* Provide default: */ hash = &ulist_uint64_hash; }
    if (hash != S->hash)
      { demand((S->ct == 0), "cannot change {hash} of non-empty list"); 
        S->hash = hash;
      }
  }

ulist_count_t ulist_count(ulist_t *S)
  { return S->ct; }

ulist_count_t ulist_capacity(ulist_t *S)
  { return S->it.ne; }

ulist_hash_size_t ulist_hash_size(ulist_t *S)
  { return S->ix.ne; }

ulist_item_vec_t ulist_items(ulist_t *S)
  { ulist_count_t ct = S->ct;
    /* Allocate the result vector: */
    ulist_item_vec_t v = ulist_item_vec_new(ct);
    /* Copy the items to {v}: */
    ulist_count_t i;
    ulist_item_t *p; /* Pointer that scans {S->it.e}. */
    ulist_item_t *q; /* Pointer that scans {v.e}. */
    for (i = 0, p = S->it.e, q = v.e; i < ct; i++, p++, q++) { (*q) = (*p); }
    return v;
  }

ulist_item_t ulist_item_at(ulist_t *S, ulist_index_t i)
  { ulist_count_t ct = S->ct;
    /* demand(i >= 0, "invalid item index"); */
    demand(i < ct, "invalid item index");
    return S->it.e[i];
  }

/* PROCEDURES THAT DEPEND ON HASHING */

void ulist_locate(ulist_t *S, ulist_item_t a, ulist_index_t *ip, ulist_hash_val_t *kp)
  { /* Get the hash table size {nh} and the item count {ct}: */
    ulist_hash_size_t nh = S->ix.ne;
    ulist_count_t ct = S->ct;
    demand(nh > 0, "should not have called locate");
    assert(ct < nh); /* By the invariants. */
    /* Get the hash value {h} of {a}: */
    ulist_hash_val_t h = S->hash(a, nh);
    demand(h < nh, "bad hash value");
    /* Sequential serch in {ix} starting at {ix[h]}: */
    ulist_hash_val_t k = h;
    while (TRUE) { 
      /* Get the presumed index {i}: */
      #if (ulist_DEBUG)
      S->ct_prb++;
      #endif /* (ulist_DEBUG) */
      ulist_index_t i = S->ix.e[k]; 
      /* fprintf (stderr, "  %d %d\n", k, i); */
      /* Check its validity: */
      if ((i >= ct) || (S->hx.e[i] != k)) 
        { /* Entry {ix[k]} is vacant, so the search is over (not found): */
          /* assert(k >= 0); */
          assert(k < nh); 
          (*kp) = k; (*ip) = ct; 
          return;
        }
      /* Entry {ix[k]} is occupied: */
      ulist_item_t b = S->it.e[i]; 
      if (S->eq(a,b))
        { /* The search is over (found): */
          (*kp) = k; (*ip) = i;
          return;
        }
      /* Advance {k} to the next entry: */
      k++; if (k == nh) { k = 0; }
      /* Since {nh > ct}, we will stop before {k} returns to {h}. */ 
    }
  }

ulist_index_t ulist_index_of(ulist_t *S, ulist_item_t a)
  { 
    #if (ulist_DEBUG)
    S->ct_ind++;
    #endif /* (ulist_DEBUG) */
    /* Get the item count {ct}: */
    ulist_count_t ct = S->ct;
    if (ct == 0) { /* List is empty: */ return 0; } 
    /* Locate the index {i} and hash table slot {k} of {~a}: */
    ulist_index_t i;
    ulist_hash_val_t k;
    ulist_locate(S, a, &i, &k);
    return i;
  }

bool_t ulist_has(ulist_t *S, ulist_item_t a)
  { return (ulist_index_of(S,a) != S->ct); }

ulist_index_t ulist_insert_last(ulist_t *S, ulist_item_t a)
  { 
    #if (ulist_DEBUG)
    S->ct_add++;
    #endif /* (ulist_DEBUG) */
    /* Get the list count {ct} and capacity {sz}: */
    ulist_count_t ct = S->ct;
    ulist_count_t sz = S->it.ne;
    /* Can we add another element? */
    demand(ct < sz, "list is full");
    /* !!! implement auto-resize if full. !!! */
    /* Locate the index {i} and hash table slot {k} of {~a}: */
    ulist_index_t i;
    ulist_hash_val_t k;
    ulist_locate(S, a, &i, &k);
    if ((i == 749) && (k == 2180)) { (void)ulist_verify(S, TRUE); }
    demand(i == ct, "element is already in table");
    /* Insert item: */
    S->it.e[i] = a;
    S->hx.e[i] = k;
    S->ix.e[k] = i;
    S->ct++;
    #if (ulist_PARANOIA)
    (void)ulist_verify(S, TRUE);
    #endif /* (ulist_PARANOIA) */
    return i;
  }

void ulist_swap(ulist_t *S, ulist_index_t i, ulist_index_t j)
  { /* Get the list count {ct}: */
    ulist_count_t ct = S->ct;
    /* Validate the indices: */
    demand(i < ct, "bad index i");
    demand(j < ct, "bad index j");
    /* Trivial case: */
    if (i == j) { return; }
    /* Swap the entries of {it[]}: */
    { ulist_item_t t = S->it.e[i]; S->it.e[i] = S->it.e[j]; S->it.e[j] = t; }
    /* Swap the entries of {hx[]}: */
    { ulist_hash_val_t t = S->hx.e[i]; S->hx.e[i] = S->hx.e[j]; S->hx.e[j] = t; }
    /* Fix the links in {ix[]}: */
    S->ix.e[S->hx.e[i]] = i;
    S->ix.e[S->hx.e[j]] = j;
  }

void ulist_resize(ulist_t *S, int32_t n, ulist_hash_func_t *hash)
  { 
    #if (ulist_DEBUG)
    S->ct_rsz++;
    #endif /* (ulist_DEBUG) */
    if (hash == NULL) { /* Provide default: */ hash = &ulist_uint64_hash; }
    ulist_count_t sz = S->it.ne;
    if ((sz != n) || (S->hash != hash))
      { /* Must actually do a resize: */
        /* Grab the item vector {it} of {S}: */
        ulist_item_vec_t old_it = S->it;
        /* Replace {S->it} by a new item vector with {n} slots: */
        S->it = ulist_item_vec_new(n);
        /* Make {S} empty and change its hash function: */
        S->ct = 0;
        S->hash = hash;
        /* Items are non-eq already, so can use {==} instead of {eq} for reinsertion: */
        ulist_eq_func_t *eq = S->eq;
        S->eq = ulist_uint64_eq;
        /* Insert all items of {old_it} back into {S}. */
        ulist_index_t i;
        ulist_item_t *p; /* Pointer that scans the item slots. */
        for (i = 0, p = old_it.e; i < sz; i++, p++)
          { ulist_item_t a = (*p);
            ulist_insert_last(S, a);
          }
        /* Restore the original equivalence: */
        S->eq = eq;
        /* Reclaim the old vector: */
        free(old_it.e);
      }
  }

/* DEBUGGING */

bool_t ulist_verify(ulist_t *S, bool_t die)
  { ulist_count_t ct = S->ct;         /* Current item count. */
    ulist_count_t sz = S->it.ne;     /* Capacity */
    if (sz != S->hx.ne) { fail_test(die, "inconsistent sizes of it,hx"); }
    if (ct > sz) { fail_test(die, "size greater than capacity"); }
    ulist_hash_size_t nh = S->ix.ne; /* Hash table size. */
    if (sz == 0)
      { if (nh != 0) { fail_test(die, "zero capacity list with hash table"); } }
    else
      { if (nh <= sz) { fail_test(die, "hash table too small"); } }
    /* Check the structural invariants: */
    int32_t i;
    for (i = 0; i < ct; i++)
      { ulist_item_t a = S->it.e[i];
        ulist_hash_val_t k = S->hx.e[i];
        if (k > nh) { fail_test(die, "invalid hash position hx[i]"); }
        if (S->ix.e[k] != i) { fail_test(die, "ix[hx[i]] is not i"); }
        /* Check the seqsearch invariant: */
        ulist_index_t j;
        ulist_hash_val_t kj;
        ulist_locate(S, a, &j, &kj);
        if (j != i) 
          { fprintf(stderr, "%s:%d: ** seq search invariant violated for it[%u]\n", __FILE__, __LINE__, i); 
            fprintf(stderr, "  sz = %u ct = %u nh = %u\n", sz, ct, nh); 
            ulist_hash_val_t ha = S->hash(a,nh);
            fprintf(stderr, "  hash(it[%u]) = %u  locate = it[%u],ix[%u]", i, ha, j,kj); 
            if (j < ct) 
              { ulist_item_t b = S->it.e[j];
                ulist_hash_val_t hb = S->hash(b,nh);
                fprintf(stderr, "  hash(it[%u]) = %u\n", j, hb);
              }
            else
              { fprintf(stderr, "  (not in list)\n"); }
            fail_test(die, "aborting");
          }
      }
    return TRUE;
  }

void ulist_stats_print(ulist_t *S)
  {
    #if (ulist_DEBUG)
    fprintf(stderr, "........................................................................\n");
    fprintf(stderr, "operation counts:\n");
    /* Basic counts: */
    fprintf(stderr, ("  {ulist_insert_last} = %12" int64_d_fmt "\n"), S->ct_add);
    fprintf(stderr, ("  {ulist_delete_last} = %12" int64_d_fmt "\n"), S->ct_del);
    fprintf(stderr, ("  {ulist_index_of} = %12" int64_d_fmt "\n"), S->ct_ind);
    fprintf(stderr, ("  {ulist_resize} = %12" int64_d_fmt "\n"), S->ct_rsz);
    fprintf(stderr, ("  probes in hash table = %12" int64_d_fmt "\n"), S->ct_prb);
    /* Probes per operation (including the resize overhead): */
    double ppo = ((double)S->ct_prb)/((double)(S->ct_add + S->ct_del + S->ct_ind)); 
    fprintf(stderr, "  probes per operation = %8.2f\n", ppo);
    fprintf(stderr, "........................................................................\n");
    #else 
    fprintf(stderr, "%s compiled without statistics support\n", __FILE__);
    #endif /* (ulist_DEBUG) */
  }

void ulist_stats_clear(ulist_t *S)
  {
    #if (ulist_DEBUG)
    S->ct_add = 0;
    S->ct_del = 0;
    S->ct_ind = 0;
    S->ct_prb = 0;
    S->ct_rsz = 0;
    #endif /* (ulist_DEBUG) */
  }

/* DEFAULT EQ/HASH FUNCTIONS */

ulist_hash_val_t ulist_uint64_hash(ulist_item_t a, ulist_hash_size_t n)
  { demand(n > 0, "cannot hash to empty range");
    /* Let's hope that this works: */
    uint64_t ua = (uint64_t) a;
    uint64_t ub = ua & 8191ULL;
    uint64_t ur = (4615ULL * (ua >> 14)) + (417ULL * (ua >> 9)) + (ub * ub * 471703ULL);
    ulist_hash_val_t u = (ulist_hash_val_t)(ur % n);
    return u;
  }

bool_t ulist_uint64_eq(ulist_item_t a, ulist_item_t b)
  { 
    return (a == b);
  }
  
