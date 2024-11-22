/* See intmerge.h */
/* Last edited on 2024-11-17 15:54:23 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <affirm.h>

#include <jsmath.h>
#include <intmerge.h>

void imrg_merge(int32_t *a, int32_t *b, int32_t *c, int32_t cmp(int32_t x, int32_t y), int32_t sgn)
  {

    auto void merge(int32_t *am, int32_t *bm, int32_t *cm);
    /* Merges two consecutive sorted blocks {*am..*(bm-1)} and {*bm..*(cm-1)}. */

    if ((a < b) && (b < c)) { merge(a, b, c); }
    return;

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

    auto void swap_blocks_by_cycles (int32_t *as, int32_t *bs, int32_t *cs); 
      /* Swaps the consecutive blocks starting at {as} and {bs} and
        ending at {cs}, by decomposition into cycles. */

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
      
//      auto void flip_block (int32_t *ar, int32_t *br);
//      /* Reverses the element block from {ar} (inclusive) to {br} (exclusive). */
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

    void swap_blocks_by_cycles (int32_t *as, int32_t *bs, int32_t *cs)
      { if ((as == bs) || (bs == cs)) return;
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

    void merge(int32_t *am, int32_t *bm, int32_t *cm)
      { 
        if ((am < bm) && (bm < cm))
          { 
            /* fprintf(stderr, "+ merge [%d..%d] [%d..%d]\n", am-a, bm-1-a, bm-a, cm-1-a); */
            uint32_t na = (uint32_t)(bm - am);
            uint32_t nb = (uint32_t)(cm - bm);
            if ((na == 1) && (nb == 1))
              { /* Either swap or do nothing: */
                if (cmp(am[0],bm[0])*sgn > 0) 
                  { int32_t t = am[0]; am[0] = bm[0]; bm[0] = t; }
              }
            else
              { int32_t *as, *bs; /* Splits in the {am} and {bm} blocks. */
                int32_t *d;  /* Address of pivot after block swap. */
                if (na >= nb)
                  { uint32_t ka = na/2; 
                    as = am + ka;
                    d = as;
                    /* fprintf(stderr, "  pivot h[%d] = %d\n", d-a, *d); */
                    bs = locate_lower_pivot(*d, bm, cm);
                    d = d + (bs - bm); /* After block swap */
                  }
                else
                  { uint32_t kb = nb/2; 
                    bs = cm - kb;
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
  }

/* Created by J. Stolfi, Unicamp, Oct/2004.
  
  Inspired on a Java implementation by Thomas Baudel, 1998, which was
  in turn inspired on the STL C++ library. This C version fixes a
  perfomance bug and makes different use of pointers versus integers.

*/
