/* See intsort_extra.h */
/* Last edited on 2024-11-17 15:56:59 by stolfi */

#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <affirm.h>
#include <jsmath.h>
#include <intsort.h>
#include <intmerge_extra.h>

#include <intsort_extra.h>

void isrt_heapsort_classic(int32_t *h, uint32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  { 
    if (n <= 1) { return; }
    
    /* Arrange elements into a sorted heap, in reverse order: */
    uint32_t m = 0;
    while (m < n) 
      { /* Grab next element {v}: */
        int32_t v = h[m];
        /* Create a vacancy {h[i]} at the end of the heap: */
        uint32_t i = m, j;
        m++;
        /* Bubble it rootwards to its proper place {h[i]}: */
        while ((i > 0) && (sgn*cmp(v,h[j=(i-1)/2]) > 0)) { h[i] = h[j]; i = j; }
        h[i] = v;
      }
    /* Pop them out in reverse order and store them at the end: */
    while (m > 1) 
      { /* Swap root with last element, decrement heap: */
        m--; int32_t w = h[m]; h[m] = h[0]; h[0] = w;
        /* Bubble the new root down to its proper place: */
        uint32_t j = 0;
        uint32_t ia = 1; /* {h[ia]} is the first child of {h[j]}. */
        while (ia < m)
          { /* Find largest child {h[i]} of {h[j]}: */
            uint32_t ib = ia + 1; /*  {h[ib]} is the SECOND child of {h[j]}. */
            uint32_t i = ((ib < m) && (sgn*cmp(h[ia],h[ib]) < 0) ? ib : ia);
            /* Ensure that {w=h[j]} precedes the smallest child: */
            if (sgn*cmp(w,h[i]) < 0)
              { /* Swap {h[j]} with child {h[i]}: */
                h[j] = h[i]; h[i] = w; 
                j = i;
              }
            else
              { /* Item {h[j]} is correctly placed, stop: */
                j = m;
              }
            ia = 2*j + 1;
          }
      }
  }

void isrt_heapsort_vacsink(int32_t *h, uint32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  { 
    if (n <= 1) { return; }
    
    /* Arrange elements into a sorted heap, in reverse order: */
    uint32_t m = 0;
    while (m < n) 
      { /* Insert at bottom of heap: */
        uint32_t i = m;
        m++;
        /* Grab dubious element {v}: */
        int32_t v = h[i];
        uint32_t j;
        /* Bubble it rootwards to its proper place {h[i]}: */
        while ((i > 0) && (sgn*cmp(v,h[j=(i-1)/2]) > 0)) { h[i] = h[j]; i = j; }
        h[i] = v;
      }
    /* Pop them out in reverse order and store them at the end: */
    while (m > 1) 
      { /* Remove the root: */
        int32_t w = h[0]; 
        uint32_t i = 0;  /* {h[i]} is a vacant slot. */
        /* Promote children into vacancy {h[i]} until it reaches the fringe: */
        uint32_t ja = 1; /* {h[ja]} is the first child of {h[i]}. */
        while (ja < m)
          { /* Find largest child {h[j]} of {h[i]}: */
            uint32_t jb = ja + 1; /*  {h[jb]} is the SECOND child of {h[i]}. */
            uint32_t j = ((jb < m) && (sgn*cmp(h[ja],h[jb]) < 0) ? jb : ja);
            /* Promote largest child into hole: */
            h[i] = h[j]; i = j; ja = 2*i + 1;
          }
        m--;
        if (i < m)
          { /* The vacancy did not end up at {h[m]}, so fill it with {h[m]}: */        
            int32_t v = h[m]; 
            /* Bubble it up to the proper place: */
            while (i > 0)
              { uint32_t j = (i-1)/2;
                if (sgn*cmp(v,h[j]) <= 0) { break; }
                h[i] = h[j];
                i = j; 
              }
            h[i] = v;
          }
        /* Add the removed element to the sorted region: */
        h[m] = w;
      }
  }

void isrt_quicksort_middle(int32_t *h, uint32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  {
    auto void sort(uint32_t i, uint32_t j); 
      /* Sorts segment {h[i..j-1]} */

    if (n >= 2) { sort(0, n); }
    return;

    void sort(uint32_t i, uint32_t j)
      {
        while (j-i > isrt_quicksort_middle_SMALL)
          { /* Separator element is middle one: */
            uint32_t m = (i+j-1)/2;
            int32_t v = h[m];
            /* Put the first element h[i] in its place: */
            h[m] = h[i];
            /* Sift {h[i+1..j-1]} into leqs {h[i+1..r-1]} and geqs {h[s+1..j-1]}: */
            uint32_t r = i+1; uint32_t s = j-1;
            do
              { while ((r <= s) && (sgn*cmp(h[r],v) <= 0)) { r++; }
                while ((s >= r) && (sgn*cmp(v,h[s]) <= 0)) { s--; }
                if (r <= s)
                  { /* Note {r == s} is impossible here. */
                    affirm(r < s, "program error");
                    int32_t t = h[r]; h[r] = h[s]; h[s] = t; r++; s--;
                  }
              }
            while (r <= s);
            /* Now leqs = {h[i+1..s]}, geqs = {h[r..j-1]}. */
            /* Insert {v} back between leqs and geqs: */
            h[i] = h[s]; h[s] = v;
            /* Now leqs = {h[i..s-1]}, geqs = {h[r..j-1]}. */
            /* Recurse on smaller half, iterate on larger one: */
            if (s - i < j - r)
              { if (s - i >= 2) { sort(i, s); } i = r; }
            else
              { if (j - r >= 2) { sort(r, j); } j = s; }
          }
        if (j-i >= 2) 
          { /* Finish with insertion sort */
            isrt_quicksort_middle_SMALLSORT(h+i, j-i, cmp, sgn);
          }
      }
  }

void isrt_quicksort_median3(int32_t *h, uint32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  {
    auto void sort(uint32_t i, uint32_t j); 
      /* Sorts segment {h[i..j]} */

    if (n >= 2) { sort(0, n-1); }
    return;

    void sort(uint32_t i, uint32_t j)
      {
        while (j-i+1 > isrt_quicksort_median3_SMALL)
          { /* Separator element is median of {h[i],h[m],h[j-1]}: */
            if (sgn*cmp(h[i],h[j]) > 0)
              { int32_t t = h[i]; h[i] = h[j]; h[j] = t; }
            uint32_t m = (i+j)/2;
            int32_t v = h[m];
            if (sgn*cmp(h[i],v) > 0)
              { h[m] = h[i]; h[i] = v; v = h[m]; }
            else if (sgn*cmp(v,h[j]) > 0)
              { h[m] = h[j]; h[j] = v; v = h[m]; }
            /* Sift {h[i+1..j-1]} into leqs {h[i+1..r-1]} and geqs {h[s+1..j-1]}: */
            uint32_t r = i+1; uint32_t s = j-1;
            do
              { /* Now {r <= s}, {h[i..r-1] <= v}, {h[s+1..j] >= v}. */
                while ((r <= s) && (sgn*cmp(h[r],v) <= 0)) { r++; }
                while ((r <= s) && (sgn*cmp(v,h[s]) <= 0)) { s--; }
                /* Now {h[i..r-1] <= v} and {h[s+1..j] >= v} and */
                /*   if {r <= s} then {h[r] > v > h[s]}. */
                if (r <= s)
                  { int32_t t = h[r]; h[r] = h[s]; h[s] = t; r++; s--; }
                /* Now {h[i..r-1] <= v} and {h[s+1..j] >= v} */
              }
            while (r <= s);
            affirm((s >= i) && (r <= j), "bad r,s");
            /* Now leqs = {h[i..s]}, geqs = {h[r..j]}. */
            /* Recurse on smaller half, iterate on larger one: */
            if (s - i < j - r)
              { if (i < s) { sort(i, s); } i = r; }
            else
              { if (r < j) { sort(r, j); } j = s; }
          }
        if (i < j) 
          { /* Finish with insertion sort */
            isrt_quicksort_median3_SMALLSORT(h+i, j-i+1, cmp, sgn);
          }
      }
  }
  
void isrt_mergesort_pivot(int32_t *h, uint32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  {
    auto void sort(int32_t *a, int32_t *b);
    /* Sorts elements of {h} from {*a} to {*(b-1)} (inclusive).
      Assumes that {a<=b} iff there are any elements, {a>b} 
      if there are none. */
  
    auto void split (int32_t *af, int32_t *bf, int32_t *cf, int32_t **as, int32_t **bs, int32_t **piv);
    /* Given two consecutive sorted blocks {A,B} of entries in {h}, splits them
      into two consecutive sub-blocks each, respectively {A1,A2} and {B1,B2},
      and a single element {*piv} between either such that, in the order {sgn*cmp},
      (1) all elements in {A1} preced or are equivalent to the elements
      in {B2}; (2) all elements in {A2} strictly follow those in {B1};
      and (3) both pairs {A1+B2} and {A2+B1} are large compared to 
      the total {A+B}.  Moreover the parameter {*piv} is made to point to an element
      that, after {A2} and {B1} are swapped, will be in its correct place. */

    auto int32_t *locate_lower_pivot (int32_t val, int32_t *x, int32_t *y);
    /* Given two pointers {x,y} into a sorted sub-array of {h}, with
      {x<y}, returns {z} such that {*x..*(z-1)} strictly precede {val}
      in the order {sgn*cmp}, while {*z..*(y-1)} are equivalent to {val}
      or follow it in the order. */

    auto int32_t *locate_upper_pivot (int32_t val, int32_t *x, int32_t *y);
    /* Given two pointers {x,y} into a sorted sub-array of {h}, with
      {x<y}, returns {z} such that {*z..*(y-1)} strictly follow {val}
      in the order {sgn*cmp}, while {*x..*(z-1)} are equivalent to {val} or
      preced it in the order. */

    auto void flip_block (int32_t *ar, int32_t *br);
    /* Reverses the element block from {*ar} to {*br} (note: both inclusive). */

//      auto void swap_blocks_by_cycles (int32_t *as, int32_t *bs, int32_t *cs); 
//      /* Swaps the consecutive blocks starting at {as} and {bs} and ending at {cs}.
//        Uses the cycle decomposition algorithm. */

    auto void swap_blocks_by_flips (int32_t *as, int32_t *bs, int32_t *cs); 
    /* Swaps the consecutive blocks starting at {as} and {bs} and ending at {cs}.
      Uses the three-flip algorithm. */

    auto void merge_blocks(int32_t *am, int32_t *bm, int32_t *cm);
    /* Merges two consecutive sorted blocks {*am..*(bm-1)} and {*bm..*(cm-1)}. */

    if (n >= 2) { sort(h, h+n); }
    return;

    void sort(int32_t *as, int32_t *bs)
      { assert(as < bs);
        uint32_t nab = (uint32_t)(bs - as);
        if (nab <= 1) 
          { /* Do nothing. */ }
        else if (nab <= isrt_mergesort_pivot_SMALL) 
          { isrt_mergesort_pivot_SMALLSORT(as, nab, cmp, sgn); }
        else
          { int32_t *ms = as + nab/2;
            sort (as, ms);
            sort (ms, bs);
            merge_blocks(as, ms, bs);
          }
      }

    void split (int32_t *af, int32_t *bf, int32_t *cf, int32_t **as, int32_t **bs, int32_t **piv)
      { assert(bf > af);
        uint32_t na = (uint32_t)(bf - af);
        assert(cf > bf);
        uint32_t nb = (uint32_t)(cf - bf);
        if (na >= nb)
          { int32_t *tpiv = af + na/2;  /* Address of pivot before block swap. */
            /* fprintf(stderr, "  pivot h[%d] = %d\n", tpiv-a, *tpiv); */
            int32_t *bsplit = locate_lower_pivot(*tpiv, bf, cf);
            (*as) = tpiv;
            (*bs) = bsplit;
            (*piv) = tpiv + (bsplit - bf);
          }
        else
          { int32_t *tpiv = cf - nb/2 - 1;  /* Address of pivot before block swap. */
            /* fprintf(stderr, "  pivot h[%d] = %d\n", tpiv-a, *tpiv); */
            int32_t *asplit = locate_upper_pivot(*tpiv, af, bf);
            (*as) = asplit;
            (*bs) = tpiv + 1;
            (*piv) = tpiv - (bf - asplit);
          }
      }

    int32_t *locate_lower_pivot (int32_t val, int32_t *x, int32_t *y)
      { while (x != y)
          { int32_t *mid = x + (y - x)/2;
            if (cmp(val, *mid)*sgn > 0)
              { x = mid + 1; }
            else
              { y = mid; }
          }
        return x;
      }

    int32_t *locate_upper_pivot (int32_t val, int32_t *x, int32_t *y)
      { while (x != y)
          { int32_t *mid = x + (y - x)/2;
            if (cmp(*mid, val)*sgn > 0)
              { y = mid; }
            else
              { x = mid + 1; }
          }
        return x;
      }

    void flip_block (int32_t *ar, int32_t *br)
      { while (ar < br)
          { int32_t t = *ar; *ar = *br; *br = t; ar++; br--; }
      }

//      void swap_blocks_by_cycles (int32_t *as, int32_t *bs, int32_t *cs)
//        {
//          if (as == bs || bs == cs) return;
//          /* fprintf(stderr, "+ swap [%d..%d] [%d..%d]\n", as-a, bs-1-a, bs-a, cs-1-a); */
//          int32_t shift = (int32_t)(bs - as);
//          int32_t len = (int32_t)(cs - as);
//          /* fprintf(stderr, "  shift = %d\n", shift); */
//          int32_t *fold = as + (len - shift);
//          int32_t *start = as + gcd(shift, len - shift);
//          do {
//            start--;
//            int32_t *p = start; 
//            int32_t *q = p + shift;
//            int32_t val = *p;
//            /* fprintf(stderr, "  val = h[%d]\n", p-a); */
//            while (q != start)
//              { *p = *q; 
//                /* fprintf(stderr, "  h[%d] = h[%d]\n", p-a, q-a); */
//                p = q; q += shift;
//                if (p >= fold) { q -= len; }
//              }
//            *p = val;
//            /* fprintf(stderr, "  h[%d] = val\n", p-a); */
//          } while (start != as);
//          /* fprintf(stderr, "- swap [%d..%d] [%d..%d]\n", as-a, bs-1-a, bs-a, cs-1-a); */
//        }
//  

    void swap_blocks_by_flips (int32_t *as, int32_t *bs, int32_t *cs)
      {
        if (as == bs || bs == cs) return;
        /* fprintf(stderr, "+ swap [%d..%d] [%d..%d]\n", as-a, bs-1-a, bs-a, cs-1-a); */
        flip_block(as, bs-1);
        flip_block(bs, cs-1);
        flip_block(as, cs-1);
        /* fprintf(stderr, "- swap [%d..%d] [%d..%d]\n", as-a, bs-1-a, bs-a, cs-1-a); */
      }

    void merge_blocks(int32_t *am, int32_t *bm, int32_t *cm)
      { 
        if ((am < bm) && (bm < cm))
          { 
            /* fprintf(stderr, "+ merge [%d..%d] [%d..%d]\n", am-a, bm-1-a, bm-a, cm-1-a); */
            int32_t na = (int32_t)(bm - am);
            int32_t nb = (int32_t)(cm - bm);
            if ((na == 1) && (nb == 1))
              { /* Either swap or do nothing: */
                if (cmp(am[0],bm[0])*sgn > 0) 
                  { int32_t t = am[0]; am[0] = bm[0]; bm[0] = t; }
              }
            else
              { /* Split blocks {A,B} into {A1,A2,B1,B2} such that {A2,B2} must swap. */
                int32_t *as, *bs, *piv; /* Splits in the {am} and {bm} blocks. */
                split(am, bm, cm, &as, &bs, &piv);
                affirm(am <= as, "mrg bug");
                affirm(as <= bm, "mrg bug");
                affirm(bm <= bs, "mrg bug");
                affirm(bs <= cm, "mrg bug");
                affirm(as < bs, "mrg bug");
                swap_blocks_by_flips(as, bm, bs);
                affirm(as <= piv, "mrg bug");
                affirm(piv < bs, "mrg bug");
                merge_blocks(am, as, piv);
                merge_blocks(piv+1,bs,cm);
              }
            /* fprintf(stderr, "- merge [%d..%d] [%d..%d]\n", am-a, bm-1-a, bm-a, cm-1-a); */
          }
      }
  }

void isrt_mergesort_symsplit(int32_t *h, uint32_t n, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  {
    auto void sort(int32_t *a, int32_t *b);
    /* Sorts elements of {h} from {*a} to {*(b-1)} (inclusive). */
    
    if (n >= 2) { sort(h, h+n); }
    return;

    void sort(int32_t *as, int32_t *bs)
      { assert(bs > as);
        uint32_t nab = (uint32_t)(bs - as);
        if (nab <= 1) 
          { /* Do nothing. */ }
        else if (nab <= isrt_mergesort_symsplit_SMALL) 
          { isrt_mergesort_symsplit_SMALLSORT(as, nab, cmp, sgn); }
        else
          { int32_t *ms = as + nab/2;
            sort (as, ms);
            sort (ms, bs);
            imrg_merge_symsplit(as, ms, bs, cmp, sgn);
          }
      }
  }

/* Created by J. Stolfi, Unicamp, Nov/2004. */
