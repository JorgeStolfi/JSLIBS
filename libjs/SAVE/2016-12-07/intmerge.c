/* See intmerge.h */
/* Last edited on 2013-10-25 01:23:39 by stolfilocal */

#include <intmerge.h>
#include <jsmath.h>

#include <stdlib.h>
#include <affirm.h>
#include <stdio.h>

void imrg_merge(int *a, int *b, int *c, int cmp(int x, int y), int sgn)
  {
    auto int *locate_lower_pivot (int val, int *x, int *y);
    /* Given two pointers {x,y} into a sorted sub-array of {h}, with
      {x<y}, returns {z} such that {*x..*(z-1)} strictly precede {val}
      in the order {sgn*cmp}, while {*z..*(y-1)} are equivalent to {val}
      or follow it in the order. */

    int *locate_lower_pivot (int val, int *x, int *y)
      { while (x != y)
          { int *mid = x + (y - x)/2;
            if (cmp(val, *mid)*sgn > 0)
              { x = mid + 1; }
            else
              { y = mid; }
          }
        return x;
      }

    auto int *locate_upper_pivot (int val, int *x, int *y);
    /* Given two pointers {x,y} into a sorted sub-array of {h}, with
      {x<y}, returns {z} such that {*z..*(y-1)} strictly follow {val}
      in the order {sgn*cmp}, while {*x..*(z-1)} are equivalent to {val} or
      preced it in the order. */

    int *locate_upper_pivot (int val, int *x, int *y)
      { while (x != y)
          { int *mid = x + (y - x)/2;
            if (cmp(*mid, val)*sgn > 0)
              { y = mid; }
            else
              { x = mid + 1; }
          }
        return x;
      }

    auto void flip_block (int *ar, int *br);
    /* Reverses the element block from {ar} (inclusive) to {br} (exclusive). */
    void flip_block (int *ar, int *br)
      { while (ar <= br)
          { int t = *ar; br--; *ar = *br; *br = t; ar++; }
      }

    auto void swap_blocks_by_flips (int *as, int *bs, int *cs); 
    /* Swaps the consecutive blocks starting at {as} and {bs} and
      ending at {cs}, with three flips. */

    void swap_blocks_by_flips (int *as, int *bs, int *cs)
      { flip_block(as, bs);
        flip_block(bs, cs);
        flip_block(as, cs);
      }
      
    auto void swap_blocks_by_cycles (int *as, int *bs, int *cs); 
    /* Swaps the consecutive blocks starting at {as} and {bs} and
      ending at {cs}, by decomposition into cycles. */

    void swap_blocks_by_cycles (int *as, int *bs, int *cs)
      { if (as == bs || bs == cs) return;
        /* fprintf(stderr, "+ swap [%d..%d] [%d..%d]\n", as-a, bs-1-a, bs-a, cs-1-a); */
        int shift = (int)(bs - as);
        int len = (int)(cs - as);
        /* fprintf(stderr, "  shift = %d\n", shift); */
        int *fold = as + (len - shift);
        int *start = as + gcd(shift, len - shift);
        do {
          start--;
          int *p = start; 
          int *q = p + shift;
          int val = *p;
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

    auto void merge(int *am, int *bm, int *cm);
    /* Merges two consecutive sorted blocks {*am..*(bm-1)} and {*bm..*(cm-1)}. */

    void merge(int *am, int *bm, int *cm)
      { 
        if ((am < bm) && (bm < cm))
          { 
            /* fprintf(stderr, "+ merge [%d..%d] [%d..%d]\n", am-a, bm-1-a, bm-a, cm-1-a); */
            int na = (int)(bm - am);
            int nb = (int)(cm - bm);
            if ((na == 1) && (nb == 1))
              { /* Either swap or do nothing: */
                if (cmp(am[0],bm[0])*sgn > 0) 
                  { int t = am[0]; am[0] = bm[0]; bm[0] = t; }
              }
            else
              { int *as, *bs; /* Splits in the {am} and {bm} blocks. */
                int *d;  /* Address of pivot after block swap. */
                if (na >= nb)
                  { int ka = na/2; as = am + ka;
                    d = as;
                    /* fprintf(stderr, "  pivot h[%d] = %d\n", d-a, *d); */
                    bs = locate_lower_pivot(*d, bm, cm);
                    d = d + (bs - bm); /* After block swap */
                  }
                else
                  { int kb = nb/2; bs = cm - kb;
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

/* Created by J. Stolfi, Unicamp, Oct/2004.
  
  Inspired on a Java implementation by Thomas Baudel, 1998, which was
  in turn inspired on the STL C++ library. This C version fixes a
  perfomance bug and makes different use of pointers versus integers.

*/
