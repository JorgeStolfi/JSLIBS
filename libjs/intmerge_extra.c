/* See intmerge_extra.h */
/* Last edited on 2024-11-16 00:58:35 by stolfi */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <affirm.h>
#include <intmerge_extra.h>
#include <jsmath.h>

void imrg_merge_pivot(int32_t *a, int32_t *b, int32_t *c, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  {
    /* Inspired on a Java implementation by Thomas Baudel, 1998, which was
      in turn inspired on the STL C++ library. This C version fixes a
      perfomance bug, and makes different use of pointers versus integers. */
    
    auto int32_t *locate_lower_pivot (int32_t val, int32_t *x, int32_t *y);
    /* Given two pointers {x,y} into a sorted sub-array of {h}, with
      {x<y}, returns {z} such that {*x..*(z-1)} strictly precede {val}
      in the order {sgn*cmp}, while {*z..*(y-1)} are equivalent to {val}
      or follow it in the order. */

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

    auto int32_t *locate_upper_pivot (int32_t val, int32_t *x, int32_t *y);
    /* Given two pointers {x,y} into a sorted sub-array of {h}, with
      {x<y}, returns {z} such that {*z..*(y-1)} strictly follow {val}
      in the order {sgn*cmp}, while {*x..*(z-1)} are equivalent to {val} or
      preced it in the order. */

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

//      auto void flip_block (int32_t *ar, int32_t *br);
//      /* Reverses the element block from {ar} (inclusive) to {br} (exclusive). */
//  
//      void flip_block (int32_t *ar, int32_t *br)
//        { while (ar <= br)
//            { int32_t t = *ar; br--; *ar = *br; *br = t; ar++; }
//        }
//  
//      auto void swap_blocks_by_flips (int32_t *as, int32_t *bs, int32_t *cs); 
//      /* Swaps the consecutive blocks starting at {as} and {bs} and
//        ending at {cs}, with three flips. */
//  
//      void swap_blocks_by_flips (int32_t *as, int32_t *bs, int32_t *cs)
//        { flip_block(as, bs);
//          flip_block(bs, cs);
//          flip_block(as, cs);
//        }
      
    auto void swap_blocks_by_cycles (int32_t *as, int32_t *bs, int32_t *cs); 
    /* Swaps the consecutive blocks starting at {as} and {bs} and
      ending at {cs}, by decomposition into cycles. */

    void swap_blocks_by_cycles (int32_t *as, int32_t *bs, int32_t *cs)
      { if (as == bs || bs == cs) return;
        assert((uint64_t)as < (uint64_t)bs);
        /* fprintf(stderr, "+ swap [%d..%d] [%d..%d]\n", as-a, bs-1-a, bs-a, cs-1-a); */
        uint32_t shift = (uint32_t)(bs - as);
        uint32_t len = (uint32_t)(cs - as);
        /* fprintf(stderr, "  shift = %d\n", shift); */
        int32_t *fold = as + (len - shift);
        int32_t *start = as + gcd(shift, len - shift);
        do {
          start--;
          int32_t *p = start; 
          int32_t *q = p + shift;
          int32_t val = *p;
          /* fprintf(stderr, "  val = h[%d]\n", p-a); */
          while (q != start)
            { *p = *q; 
              /* fprintf(stderr, "  h[%d] = h[%d]\n", p-a, q-a); */
              p = q; q += shift;
              if (p >= fold) { q -= len; }
            }
          *p = val;
          /* fprintf(stderr, "  h[%d] = val\n", p-a); */
        } while (start != as);
        /* fprintf(stderr, "- swap [%d..%d] [%d..%d]\n", as-a, bs-1-a, bs-a, cs-1-a); */
      }

    auto void merge(int32_t *am, int32_t *bm, int32_t *cm);
    /* Merges two consecutive sorted blocks {*am..*(bm-1)} and {*bm..*(cm-1)}. */

    void merge(int32_t *am, int32_t *bm, int32_t *cm)
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
              { int32_t *as, *bs; /* Splits in the {am} and {bm} blocks. */
                int32_t *d;  /* Address of pivot after block swap. */
                if (na >= nb)
                  { int32_t ka = na/2; as = am + ka;
                    d = as;
                    /* fprintf(stderr, "  pivot h[%d] = %d\n", d-a, *d); */
                    bs = locate_lower_pivot(*d, bm, cm);
                    d = d + (bs - bm); /* After block swap */
                  }
                else
                  { int32_t kb = nb/2; bs = cm - kb;
                    d = bs - 1;
                    /* fprintf(stderr, "  pivot h[%d] = %d\n", d-a, *d); */
                    as = locate_upper_pivot(*d, am, bm);
                    d = d - (bm - as); /* After block swap */
                  }
                affirm(am <= as, "mrg bug");
                affirm(as <= bm, "mrg bug");
                affirm(bm <= bs, "mrg bug");
                affirm(bs <= cm, "mrg bug");
                affirm(as < bs, "mrg bug");
                swap_blocks_by_cycles(as, bm, bs);
                affirm(as <= d, "mrg bug");
                affirm(d < bs, "mrg bug");
                merge(am, as, d);
                merge(d+1,bs,cm);
              }
            /* fprintf(stderr, "- merge [%d..%d] [%d..%d]\n", am-a, bm-1-a, bm-a, cm-1-a); */
          }
      }

    merge(a, b, c);
  }
  
void imrg_merge_symsplit(int32_t *a, int32_t *b, int32_t *c, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  {
    auto void merge(int32_t *am, int32_t *bm, int32_t *cm);
    /* Merges two consecutive non-empty sorted blocks {A,B}, respectively from
      {*am} through {*(bm-1)} and {*bm} through {*(cm-1)}, all inclusive. */

    auto void flip_block (int32_t *ar, int32_t *br);
    /* Reverses the element block from {*ar} to {*br} (note: both inclusive). */

    auto int32_t *swap_blocks_by_flips(int32_t *as, int32_t *bs, int32_t *cs); 
    /* Swaps the consecutive blocks {A = *as..*(bs-1)} and {B = *bs..*(cs-1)}.
      Uses the three-flip algorithm.  Returns the boundary between the two 
      blocks after the swap, i.e. the address where {*as} was moved to. */

    void flip_block (int32_t *ar, int32_t *br)
      { while (ar < br)
          { int32_t t = *ar; *ar = *br; *br = t; ar++; br--; }
      }

    int32_t *swap_blocks_by_flips (int32_t *as, int32_t *bs, int32_t *cs)
      { /* fprintf(stderr, "  swap [%d..%d] [%d..%d]\n", as-a, bs-1-a, bs-a, cs-1-a); */
        if ((as < bs) && (bs < cs))
          { flip_block(as, bs-1);
            flip_block(bs, cs-1);
            flip_block(as, cs-1);
          }
        /* fprintf(stderr, "- swap [%d..%d] [%d..%d]\n", as-a, bs-1-a, bs-a, cs-1-a); */
        return as + (cs - bs);
      }

   auto void merge(int32_t *am, int32_t *bm, int32_t *cm)
      { 
        /* fprintf(stderr, "+ merge [%d..%d] [%d..%d]\n", am-a, bm-1-a, bm-a, cm-1-a); */
        int32_t na = (int32_t)(bm - am);
        int32_t nb = (int32_t)(cm - bm);
        if ((na == 1) && (nb == 1))
          { /* Only one element in each list, just compare and swap as needed: */
            if (sgn*cmp(*am,*bm) > 0) 
              { int32_t t = *am; *am = *bm; *bm = t;
                /* fprintf(stderr, "  exch %d %d\n", am-a, bm-a); */
              }
          }
        else
          { /* 
              Let {A,B} be the two blocks to be merged. We first split
              each of them into two consecutive sub-blocks, respectively
              {A1,A2} and {B1,B2}, such that, in the order {sgn*cmp}, (1)
              all elements in {A1} precede or are equivalent to the
              elements in {B2}; (2) all elements in {A2} strictly follow
              those in {B1}; and (3) both pairs {A1+B2} and {A2+B1} have
              nearly the same size.

              The split points are two addresses {as,bs} such that {am <=
              as <= bm} and {bm <= bs <= cm}. Block {A2} starts at {*as},
              and is empty if {as == bm}. Similarly {B2} starts at {*bs},
              and is empty if {bs == cm}.
            */

            int32_t *as, *bs; int32_t n;
            if (na < nb)
              { n = na; as = am; bs = cm - (nb - na + 1)/2; }
            else
              { n = nb; as = am + (na - nb + 1)/2; bs = cm; }
            affirm((am <= as) && (as+n <= bm), "bug");
            affirm((bm <= bs-n) && (bs <= cm), "bug");

            /* 
              Loop invariant: the candidate pairs of split points are
              {(as+k,bs-k)} for {k} in {0..n}. Elements {*as} through
              {*(as+n-1)} are inside {A}, and {*(bs-n)} through
              {*(bs-1)} are inside {B}. We know already that {$(as-1)}
              must precede or is equivalent to {$bs}, and {$(as+n)}
              must follow {$(bs-n-1)}, where {$p} means {*p} if elem
              {*p} is inside the corresponding block, and {-oo} or
              {+oo} otherwise. If {n} is zero, we have our split
              point. We locate it by binary search on {k}.
            */
            while (n > 0)
              { 
                affirm((am <= as) && (as+n <= bm), "bug");
                affirm((bm <= bs-n) && (bs <= cm), "bug");
                affirm((as == am) || (bs == cm) || (sgn*cmp(*(as-1),*bs) <= 0), "bug");
                affirm((as+n == bm) || (bs-n == bm) || (sgn*cmp(*(as+n),*(bs-n-1)) > 0), "bug");

                int32_t half = n/2;
                int32_t *ap = as + half, *bp = bs - half;
                affirm((am <= ap) && (ap < bm), "bug");
                affirm((bm < bp) && (bp <= cm), "bug");
                if (sgn*cmp(*ap,*(bp-1)) <= 0)
                  { as = ap + 1; bs = bp - 1; n = n - half - 1; }
                else
                  { n = half; }
              }

            /* Sanity check: */
            affirm((as >= am) && (as <= bm), "bug");
            affirm((bs >= bm) && (bs <= cm), "bug");
            affirm((as == am) || (bs == cm) || (sgn*cmp(*(as-1),*bs) <= 0), "bug");
            affirm((as == bm) || (bs == bm) || (sgn*cmp(*as,*(bs-1)) > 0), "bug");

            /* fprintf(stderr, "  split = %d %d\n", as-a, bs-a); */
            
            /* Now we just swap {A2} with {B1}, and recurse: */
            bm = swap_blocks_by_flips(as, bm, bs);
            if ((am < as) && (as < bm)) { merge(am, as, bm); }
            if ((bm < bs) && (bs < cm)) { merge(bm, bs, cm); }
          }
        /* fprintf(stderr, "- merge [%d..%d] [%d..%d]\n", am-a, bm-1-a, bm-a, cm-1-a); */
      }

    if ((a < b) && (b < c)) merge(a, b, c);
  }

/* Created by J. Stolfi, Unicamp, Nov/2004. */
