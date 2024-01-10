/* FOI routines */

#include "foi.h"
#include "foifloat.h"
#include "foimisc.h"
#include "interval.h"
#include "iomisc.h"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#define QUICKMUL 1

#define FOI_TRACE_ADD     0
#define FOI_TRACE_MUL     0
#define FOI_TRACE_SQR     0
#define FOI_TRACE_NEG     0
#define FOI_TRACE_INV     0
#define FOI_TRACE_SQRT    0
#define FOI_TRACE_AFFINE  0
#define FOI_TRACE_CONST   0
#define FOI_TRACE_FROM_IV 0

/*** PROTOTYPES FOR INTERNAL PROCEDURES ***/

void foi_mul_naf_range(
    FOIP x, 
    FOIP y, 
    FloatP zcp, 
    FloatP errp
  );

void foi_sqr_naf_range(
    FOIP x, 
    FloatP zcp, 
    FloatP errp
  );

void foi_sqrt_approx(
    Interval xr,
    FloatP alphap,
    FloatP betap,
    FloatP gammap
  );

void foi_inv_approx(
    Interval xr,
    FloatP alphap,
    FloatP betap,
    FloatP gammap
  );
  
Float foi_throw_coef(void);

/*** ZERO ***/

static FOIHead FOI_ZERO_BDY = {
      #if MIXED
        {Zero, Zero},
      #endif
      0,
      Zero
    };

#define FOI_ZERO ((FOIP) &FOI_ZERO_BDY)

/*** NOISEVAR ID'S ***/

static VarId foi_next_id = 0; /* next unused noisevar id */

/*** STACK ALLOCATION ***/

#define FOI_STACK_SIZE 100000L

static MemP foi_stack_bot;  /* bottom-of-stack pointer. */
static MemP foi_stack_top;  /* The top of the FOI stack. */
static MemP foi_stack_lim;  /* limiting pointer for stack memory area. */

MemP foi_top (void)
  { return(foi_stack_top); }

#define FOI_VALID_FRAME(f) (((f) <= foi_stack_top) && ((f) >= foi_stack_bot))

#define FOI_IN_FRAME(p, f) (((f) <= (p)) && ((p) < foi_stack_top))

void foi_flush (MemP frame)
  { assert(FOI_VALID_FRAME(frame), "foi_flush: bad frame pointer");
    foi_stack_top = frame;
  }

FOIP foi_return (MemP frame, FOIP result)
  { assert(FOI_VALID_FRAME(frame), "foi_return: bad frame pointer");
    if (FOI_IN_FRAME((MemP) result, frame))
      { MemSize nbytes = sizeof(FOIHead);
        nbytes += (result->nterms) * sizeof(Term);
        if (((MemP) result) > ((MemP) frame))
          { memcpy((char *) frame, (char *) result, nbytes); }
        foi_stack_top = (MemP)(((char *) frame) + nbytes);
      }
    else
      { foi_stack_top = frame; }
    return((FOIP) frame);
  }

FOIP foi_alloc_head(void)
  {
    MemP old_top = foi_stack_top;
    foi_stack_top = (MemP) (((char *) foi_stack_top) + sizeof(FOIHead));
    if (((MemP) foi_stack_top) > ((MemP) foi_stack_lim))
      { fprintf(stderr, "old_top =       %p\n", old_top);
        fprintf(stderr, "foi_stack_top = %p\n", foi_stack_top);
        fprintf(stderr, "foi_stack_lim = %p\n", foi_stack_lim);
        error("foi_alloc_head: not enough space for FOI head");
      }
    return (FOIP) old_top;
  }

void foi_pop_head(void)
  { foi_stack_top = (MemP)(((char *) foi_stack_top) - sizeof(FOIHead)); }

TermP foi_alloc_term(void)
  {
    MemP old_top = foi_stack_top;
    foi_stack_top = (MemP)(((char *) foi_stack_top) + sizeof(Term));
    if (((MemP) foi_stack_top) > ((MemP) foi_stack_lim))
      { fprintf(stderr, "old_top =       %p\n", old_top);
        fprintf(stderr, "foi_stack_top = %p\n", foi_stack_top);
        fprintf(stderr, "foi_stack_lim = %p\n", foi_stack_lim);
        error("foi_alloc_term: not enough space for FOI term");
      }
    return (TermP) old_top;
  }

void foi_pop_term(void)
  { foi_stack_top = (MemP) (((char *) foi_stack_top) - sizeof(Term)); }

/*** INITIALIZATION ***/

void foi_init (void)
  {
    long dif;
    foi_stack_bot = (MemP) malloc(FOI_STACK_SIZE);
    fprintf(stderr, "foi_stack_bot = %p\n", foi_stack_bot);
    if (foi_stack_bot == (MemP) NULL)
      error("foi_init: heap not allocated");
    foi_stack_lim = (MemP) (((char *) foi_stack_bot) + FOI_STACK_SIZE);
    fprintf(stderr, "foi_stack_lim = %p\n", foi_stack_lim);
    dif = ((char *) foi_stack_lim) - ((char *) foi_stack_bot);
    fprintf(stderr, "stack size = %ld bytes\n", dif);
    foi_stack_top = foi_stack_bot;
    fprintf(stderr, "foi_stack_top = %p (%ld)\n", 
      foi_stack_top, (unsigned long) foi_stack_top
    );

    foi_next_id = 0;
  }

/*** HEAP ALLOCATION ***/

FOIP foi_heap_alloc (TermCount n)
  {
    FOIP r;
    MemSize nbytes = sizeof(FOIHead);
    nbytes += n * sizeof(Term);
    r = (FOIP) malloc(nbytes);
    if (r == (FOIP) NULL)
      { error("out of memory"); }
    return (r);
  }

FOIP foi_heap_realloc (FOIP x, TermCount n)
  { FOIP r;
    MemSize nbytes = sizeof(FOIHead);
    nbytes += n * sizeof(Term);
    r = (FOIP) realloc((void *) x, nbytes);
    if (r == (FOIP) NULL)
      { error("out of memory"); }
    return (r);
  }

void foi_heap_free (FOIP x)
  { free((void *) x); }

/*** MISCELLANEOUS ***/

void foi_print (FILE *f, FOIP x)
  {
    iv_print (f, foi_range(x));
    if (FOI_ISFULL(x))
      { int i;
        for (i=0; i<F_FMT_WIDTH; i++) putc('*', f);
        return;
      }
    else
      { TermP xp = (TermP) (x + 1);
        TermCount xn = x->nterms;

        flt_print (f, x->center);

        while (xn > 0)
          {
            if (xp->coef >= Zero)
              fprintf(f, " + ");
            else
              fprintf(f, " - ");
            flt_print(f, FABS(xp->coef)) ;
            fprintf(f, "(%d)", xp->id);
            xp++; xn--;
          }
      }
  }

void foi_move (FOIP source, FOIP destination)
  { if ( (source != destination)
      && (source != (FOIP) NULL)
      && (! (FOI_ISFULL(source)))
    ) { MemSize nbytes = sizeof(FOIHead);
        nbytes += (source->nterms) * sizeof(Term);
        memcpy((char *) destination, (char *) source, nbytes);
      }
  }

Interval foi_range(FOIP x)
  {
    if (FOI_ISFULL(x)) return (IV_FULL);
    #if MIXED
      return (x->range);
    #else
      return (foi_implicit_range(x));
    #endif
  }

Interval foi_implicit_range (FOIP x)
  {
    if (FOI_ISFULL(x))
      return (IV_FULL);
    else
      {
        Interval range;
        Float tdev = foi_sum_abs_terms((TermP)(x+1), x->nterms);
        if (tdev >= PlusInfinity) return (IV_FULL);
        ROUND_DOWN;
        range.lo = x->center - tdev;
        ROUND_UP;
        range.hi = x->center + tdev;
        IV_NORMFULL(range);
        return (range);
      }
  }

Float foi_sum_abs_terms (TermP xp, TermCount n)
  {
    Float tdev = Zero;
    ROUND_UP;
    while (n > 0)
      {
        tdev = tdev + FABS(xp->coef);
        if (tdev >= PlusInfinity) return (PlusInfinity);
        xp++; n--;
      }
    return (tdev);
  }
  
Float foi_throw_coef(void)
  {
    int coins, j;
    Float t;

    coins = random();
    if ((coins&3) == 0)
      return (Zero);
    else if ((coins&3) == 1)
      t = One;
    else
      t = flt_random();

    coins = random();
    for (j=0; j<25; j++)
      {
	switch(coins&3) 
          {
	    case 0: t *= 2.0; break;
	    case 1: t /= 2.0; break;
	    case 2: t *= 16.0; break;
	    case 3: t /= 16.0; break;
	  }
        coins >>= 2;
        if (coins < 4) coins = random();
      }
    if (t >= PlusInfinity) return (PlusInfinity);
    if ((random()&1) == 0) t = -t;
    return (t);
  }

FOIP foi_throw(int nterms)
  {
    int coins;
    coins = random();
    if ((coins&255) == 0)
      return (FOI_FULL);
    else if ((coins&63) == 0)
      return (foi_zero());
    else
      {
        MemP frame = foi_stack_top;
        FOIP z;
        TermP zp;
        TermCount zn = 0;
        Float t;
        int i;

        z = foi_alloc_head();

        #if FOI_TRACE_THROW
            fprintf(stderr, "enter foi_throw:\n");
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  nterms = %d\n", nterms);
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif

        ROUND_NEAR;
        z->center = foi_throw_coef();
	if (FABS(z->center) >= Infinity)
	  { foi_flush(frame); return(FOI_FULL); }
    
        for (i=0; i<nterms; i++)
          {
            t = foi_throw_coef();
	    if (FABS(t) >= Infinity)
	      { foi_flush(frame); return(FOI_FULL); }
            else if (t != Zero)
              {
		zp = foi_alloc_term();
		zp->coef = t;
		zp->id = i;
		zn++;
              }
          }

        /* Store term count: */
        z->nterms = zn;
        
        /* Make sure the noise symbols used do exist: */
        if (nterms > foi_next_id) foi_next_id = nterms;

        #if MIXED
        /* Compute z->range: */
        z->range = foi_implicit_range(z);
        #endif

        #if FOI_TRACE_THROW
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit foi_throw.\n");
        #endif

        return (z);
      }
  }

/*** ARITHMETIC ***/

FOIP foi_zero(void)
  { return (FOI_ZERO); }

FOIP foi_add(FOIP x, FOIP y)
  {
    if (FOI_ISFULL(x) || FOI_ISFULL(y))
      return (FOI_FULL);
    else
      {
        MemP frame = foi_stack_top;
        FOIP z = foi_alloc_head();
        TermP xp = (TermP) (x + 1);
        TermP yp = (TermP) (y + 1);
        TermP zp;
        Float err = Zero;
        TermCount xn = x->nterms;
        TermCount yn = y->nterms;
        TermCount zn = 0;

        #if FOI_TRACE_ADD
            fprintf(stderr, "enter foi_add:\n");
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  x = %p nterms = %ld\n", x, x->nterms);
            fprintf(stderr, "  y = %p nterms = %ld\n", y, y->nterms);
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif

        /* Add center values: */
        ROUND_NEAR;
        flt_add(x->center, y->center, &(z->center), &err);
        if (FABS(z->center) >= Infinity || err >= Infinity)
          { foi_flush(frame); return(FOI_FULL); }

        /* Merge noise terms: */
        while (xn > 0 || yn > 0)
          {
            zp = foi_alloc_term();
            if ((yn == 0) || ((xn != 0) && (xp->id < yp->id)))
              { *zp = *xp; xp++; xn--; zn++; }
            else if ((xn == 0) || ((yn != 0) && (yp->id < xp->id)))
              { *zp = *yp; yp++; yn--; zn++; }
            else
              {
                zp->id = xp->id;
                ROUND_NEAR;
                flt_add(xp->coef, yp->coef, &(zp->coef), &err);
                if (FABS(zp->coef) >= Infinity || err >= Infinity)
                  { foi_flush(frame); return(FOI_FULL); }
                xp++; xn--;
                yp++; yn--;
                if (zp->coef == Zero) foi_pop_term(); else zn++;
              }
          }

        /* Add rounding error term: */
        if (err != Zero)
          {
            zp = foi_alloc_term();
            zp->coef = err;
            zp->id = foi_next_id;
            foi_next_id++;
            zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if MIXED
        /* Compute z->range: */
        z->range = iv_meet(
          iv_add(x->range, y->range),
          foi_implicit_range(z)
        );
        #endif

        #if FOI_TRACE_ADD
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit foi_add.\n");
        #endif

        return (z);
      }
  }

FOIP foi_sub(FOIP x, FOIP y)
  {
    if (FOI_ISFULL(x) || FOI_ISFULL(y))
      return (FOI_FULL);
    else
      {
        MemP frame = foi_stack_top;
        FOIP z = foi_alloc_head();
        TermP xp = (TermP) (x + 1);
        TermP yp = (TermP) (y + 1);
        TermP zp;
        Float err = Zero;
        TermCount xn = x->nterms;
        TermCount yn = y->nterms;
        TermCount zn = 0;

        #if FOI_TRACE_SUB
            fprintf(stderr, "enter foi_sub:\n");
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  x = %p nterms = %ld\n", x, x->nterms);
            fprintf(stderr, "  y = %p nterms = %ld\n", y, y->nterms);
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif

        /* Subtract center values: */
        ROUND_NEAR;
        flt_sub(x->center, y->center, &(z->center), &err);
        if (FABS(z->center) >= Infinity || err >= Infinity)
          { foi_flush(frame); return(FOI_FULL); }

        /* Merge noise terms: */
        while (xn > 0 || yn > 0)
          {
            zp = foi_alloc_term();
            if ((yn == 0) || ((xn != 0) && (xp->id < yp->id)))
              { zp->id = xp->id; zp->coef =  xp->coef; xp++; xn--; zn++; }
            else if ((xn == 0) || ((yn != 0) && (yp->id < xp->id)))
              { zp->id = yp->id; zp->coef = -yp->coef; yp++; yn--; zn++; }
            else
              {
                zp->id = xp->id;
                ROUND_NEAR;
                flt_sub(xp->coef, yp->coef, &(zp->coef), &err);
                if (FABS(zp->coef) >= Infinity || err >= Infinity)
                  { foi_flush(frame); return(FOI_FULL); }
                xp++; xn--;
                yp++; yn--;
                if (zp->coef == Zero) foi_pop_term(); else zn++;
              }
          }

        /* Add rounding error term: */
        if (err != Zero)
          {
            zp = foi_alloc_term();
            zp->coef = err;
            zp->id = foi_next_id;
            foi_next_id++;
            zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if MIXED
        /* Compute z->range: */
        z->range = iv_meet(
          iv_add(x->range, y->range),
          foi_implicit_range(z)
        );
        #endif

        #if FOI_TRACE_SUB
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit foi_sub.\n");
        #endif

        return (z);
      }
  }

FOIP foi_neg(FOIP x)
  {
    if (FOI_ISFULL(x))
      return (FOI_FULL);
    else
      {
        FOIP z = foi_alloc_head();
        TermP xp = (TermP) (x + 1);
        TermP zp;
        TermCount xn = x->nterms;
        TermCount zn = 0;

        #if FOI_TRACE_NEG
            fprintf(stderr, "enter foi_neg:\n");
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  x = %p nterms = %ld\n", x, x->nterms);
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif


        /* Negate center value: */
        z->center = - (x->center);

        /* Copy and negate noise terms: */
        while (xn > 0)
          {
            zp = foi_alloc_term();
            zp->id = xp->id;
            zp->coef = - (xp->coef);
            xp++; xn--;
            zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if MIXED
        /* Compute z->range: */
        z->range = iv_neg(x->range);
        #endif

        #if FOI_TRACE_NEG
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit foi_neg.\n");
        #endif

        return (z);
      }
  }

FOIP foi_mul (FOIP x, FOIP y)
  {
    if (FOI_ISFULL(x) || FOI_ISFULL(y))
      return (FOI_FULL);
    else
      { MemP frame = foi_stack_top;
        FOIP z = foi_alloc_head();
        TermP zp;
        TermCount zn = 0;
        Float err = Zero;

        #if FOI_TRACE_MUL
            fprintf(stderr, "enter foi_mul:\n");
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  x = %p nterms = %ld\n", x, x->nterms);
            fprintf(stderr, "  y = %p nterms = %ld\n", y, y->nterms);
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif

        /*********************************************/
        /* I. Compute the affine part of the product */
        /*********************************************/

        {
          TermP xp = (TermP) (x + 1);
          TermP yp = (TermP) (y + 1);
          TermCount xn = x->nterms;
          TermCount yn = y->nterms;
          Float xt, yt;

          /* Multiply center values: */
          ROUND_NEAR;
          flt_mul(x->center, y->center, &(z->center), &err);
          if (FABS(z->center) >= Infinity || err >= Infinity)
            { foi_flush(frame); return(FOI_FULL); }

          /* Merge noise terms: */
          while (xn > 0 || yn > 0)
            {
              zp = foi_alloc_term();
              if ((yn == 0) || ((xn != 0) && (xp->id < yp->id)))
                { zp->id = xp->id;
                  ROUND_NEAR;
                  flt_mul(xp->coef, y->center, &(zp->coef), &err);
                  xp++; xn--;
                  zn++;
                }
              else if ((xn == 0) || ((yn != 0) && (yp->id < xp->id)))
                { zp->id = yp->id;
                  ROUND_NEAR;
                  flt_mul(x->center, yp->coef, &(zp->coef), &err);
                  yp++; yn--;
                  zn++;
                }
              else
                {
                  zp->id = xp->id;
                  ROUND_NEAR;
                  flt_mul(xp->coef, y->center, &xt, &err);
                  ROUND_NEAR;
                  flt_mul(x->center, yp->coef, &yt, &err);
                  if (FABS(xt) >= Infinity || FABS(yt) >= Infinity || err >= Infinity)
                    { foi_flush(frame); return(FOI_FULL); }
                  ROUND_NEAR;
                  flt_add(xt, yt, &(zp->coef), &err);
                  if (FABS(zp->coef) >= Infinity || err >= Infinity)
                    { foi_flush(frame); return(FOI_FULL); }
                  xp++; xn--;
                  yp++; yn--;
                  if (zp->coef == Zero) foi_pop_term(); else zn++;
                }
            }
        }

        /************************************/
        /* II. Approximate non-affine part. */
        /************************************/

        foi_mul_naf_range(x, y, &(z->center), &err);
        if (err >= PlusInfinity)
          { foi_flush(frame); return (FOI_FULL); }

        /* Add rounding error term: */
        if (err != Zero)
          {
            zp = foi_alloc_term();
            zp->coef = err;
            zp->id = foi_next_id;
            foi_next_id++;
            zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if MIXED
        /* Compute z->range: */
        z->range = iv_meet(
          iv_mul(x->range, y->range),
          foi_implicit_range(z)
        );
        #endif

        #if FOI_TRACE_MUL
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit foi_mul.\n");
        #endif

        return (z);
      }
  }

void foi_mul_naf_range(FOIP x, FOIP y, FloatP zcp, FloatP errp)
  {
    #if QUICKMUL
      Float xr = foi_sum_abs_terms((TermP)(x+1), x->nterms);
      Float yr = foi_sum_abs_terms((TermP)(y+1), y->nterms);
      Float xyr;
      if ((xr >= PlusInfinity) || (yr >= PlusInfinity))
        { *errp = PlusInfinity; return; }
      ROUND_UP;
      xyr = xr * yr;
      if (xyr >= PlusInfinity)
        { *errp = PlusInfinity; return; }
      *errp += xyr;
      if (*errp >= PlusInfinity)
        { *errp = PlusInfinity; return; }
      *zcp = *zcp; /* to supress the "not used" warning */
    #else
      NOT IMPLEMENTED
    #endif
  }

FOIP foi_sqr (FOIP x)
  {
    if (FOI_ISFULL(x))
      return (FOI_FULL);
    else
      { MemP frame = foi_stack_top;
        FOIP z = foi_alloc_head();
        TermP zp;
        TermCount zn = 0;
        Float err = Zero;

        #if FOI_TRACE_SQR
            fprintf(stderr, "enter foi_sqr:\n");
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  x = %p nterms = %ld\n", x, x->nterms);
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif

        /*********************************************/
        /* I. Compute the affine part of the square  */
        /*********************************************/

        {
          TermP xp = (TermP) (x + 1);
          TermCount xn = x->nterms;

          /* Square the center value: */
          ROUND_NEAR;
          flt_mul(x->center, x->center, &(z->center), &err);
          if (FABS(z->center) >= Infinity || err >= Infinity)
            { foi_flush(frame); return(FOI_FULL); }

          /* Process noise terms: */
          while (xn > 0)
            {
              zp = foi_alloc_term();
              zp->id = xp->id;
	      ROUND_NEAR;
	      flt_mul(xp->coef, x->center, &(zp->coef), &err);
	      flt_add(zp->coef, zp->coef, &(zp->coef), &err);
	      xp++; xn--;
	      zn++;
            }
        }

        /************************************/
        /* II. Approximate non-affine part. */
        /************************************/

        foi_sqr_naf_range(x, &(z->center), &err);
        if (err >= PlusInfinity)
          { foi_flush(frame); return (FOI_FULL); }

        /* Add rounding error term: */
        if (err != Zero)
          {
            zp = foi_alloc_term();
            zp->coef = err;
            zp->id = foi_next_id;
            foi_next_id++;
            zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if MIXED
        /* Compute z->range: */
        z->range = iv_meet(
          iv_sqr(x->range),
          foi_implicit_range(z)
        );
        #endif

        #if FOI_TRACE_SQR
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit foi_sqr.\n");
        #endif

        return (z);
      }
  }

void foi_sqr_naf_range(FOIP x, FloatP zcp, FloatP errp)
  {
    #if QUICKMUL
      Float xr = foi_sum_abs_terms((TermP)(x+1), x->nterms);
      Float xxr;
      Float nzc;
      if (xr >= PlusInfinity)
        { *errp = PlusInfinity; return; }
      ROUND_UP;
      xxr = xr * xr;
      if (xxr >= PlusInfinity)
        { *errp = PlusInfinity; return; }
      ROUND_NEAR;
      nzc = *zcp + Half * xxr;
      ROUND_UP;
      *errp += FMAX(nzc - *zcp, (*zcp + xxr) - nzc);
      if (*errp >= PlusInfinity)
        { *errp = PlusInfinity; return; }
      *zcp = nzc;
    #else
      NOT IMPLEMENTED
    #endif
  }

FOIP foi_inv (FOIP x)
  {
    if (FOI_ISFULL(x))
      { return (FOI_FULL); }
    else if (x->nterms == 0)
      {
        Float sx;
        Float err = Zero;
        if (x->center == Zero)
          { error("foi_inv: argument is zero"); }
        flt_inv(x->center, &sx, &err);
        if (FABS(sx) >= Infinity || err >= Infinity)
          { return(FOI_FULL); }
        return (foi_const(sx, err));
      }
    else
      {
        Interval xr = foi_range(x);
        Float alpha, beta, gamma;
        FOIP z;
        if (xr.lo <= Zero && xr.hi >= Zero)
          { return(FOI_FULL); }
        foi_inv_approx(xr, &alpha, &beta,&gamma);
        z = foi_affine(x, alpha, beta, gamma);
        #if MIXED
        /* Fix z->range: */
        z->range = iv_meet(iv_inv(x->range), z->range);
        #endif
        return(z);
      }
  }

void foi_inv_approx(
    Interval xr,
    FloatP alphap,
    FloatP betap,
    FloatP gammap
  )
  /* Computes chebyshev approx to $1/x$ in $xr$. */
  {
    Float axa, axb, da, db, dlo, dhi;
    int negative;

    #if FOI_TRACE_INV
        fprintf(stderr, "  enter foi_inv_approx:\n");
        fprintf(stderr, "    xr.lo = %15.8e\n", xr.lo);
        fprintf(stderr, "    xr.hi = %15.8e\n", xr.hi);
    #endif

    /* Assumes $xr$ doesn't contain zero: */
    negative = (xr.lo < Zero);
    if (negative)
      { Float t; t = xr.lo; xr.lo = -xr.hi; xr.hi = -t; }
    /* Compute approximate slope in interval: */
    ROUND_NEAR;
    *alphap = -(One/xr.lo)*(One/xr.hi);
    if (-(*alphap) >= MaxFloat)
      { *alphap = MaxFloat/Two; }
    /* Compute max and min of $1/x - \alpha x$ in $xr$: */
    ROUND_UP;
    if ( -(*alphap)*xr.hi*xr.hi <= One )
      { /* Difference is monotonically decreasing in $xr$. */
        ROUND_UP;
        axb = (*alphap) * xr.hi;
        ROUND_DOWN;
        axa = (*alphap) * xr.lo;
        dlo = One/xr.hi - axb;
        ROUND_UP;
        dhi = One/xr.lo - axa;
        assert((dlo <= dhi), "foi_inv: case 1");
      }
    else if ( -(*alphap)*xr.lo*xr.lo >= One )
      { /* Difference is monotonically increasing in $xr$. */
        ROUND_DOWN;
        axb = (*alphap) * xr.hi;
        ROUND_UP;
        axa = (*alphap) * xr.lo;
        dhi = One/xr.hi - axb;
        ROUND_DOWN;
        dlo = One/xr.lo - axa;
        assert((dlo <= dhi), "foi_inv: case 2");
      }
    else
      { /* Difference may be monotonic or concave in $xr$. */
        ROUND_DOWN;
        dlo = Two * sqrt(-(*alphap));
        ROUND_DOWN;
        axa = (*alphap) * xr.lo;
        axb = (*alphap) * xr.hi;
        ROUND_UP;
        da = One/xr.lo - axa;
        db = One/xr.hi - axb;
        dhi = FMAX(da, db);
        assert ((dlo <= dhi), "foi_inv: case 3");
      }
    ROUND_NEAR;
    *betap = (dhi + dlo)/Two;
    ROUND_UP;
    *gammap = FMAX(dhi - (*betap), (*betap) - dlo);
    if (negative) *betap = -(*betap);

    #if FOI_TRACE_INV
        fprintf(stderr,
          "    alpha = %15.8e  beta  = %15.8e  gamma = %15.8e\n",
          *alphap, *betap, *gammap
        );
        fprintf(stderr, "  exit foi_inv_approx.\n");
    #endif
  }

FOIP foi_sqrt(FOIP x)
  {
    if (FOI_ISFULL(x))
      { return (FOI_FULL); }
    else if (x->nterms == 0)
      {
        Float sx;
        Float err = Zero;
        if (x->center < Zero)
          { error("foi_sqrt: argument is negative"); }
        if (x->center == Zero)
          { return (FOI_ZERO); }
        flt_sqrt(x->center, &sx, &err);
        /* Overflow can't possibly have occurred, can it? */
        return (foi_const(sx, err));
      }
    else
      {
        Interval xr = foi_range(x);
        Float alpha, beta, gamma;
        FOIP z;
        if (xr.hi < Zero)
          { error("foi_sqrt: argument is negative"); }
        if (xr.hi == Zero)
          { return (FOI_ZERO); }
        foi_sqrt_approx(xr, &alpha, &beta,&gamma);
        z = foi_affine(x, alpha, beta, gamma);
        #if MIXED
        /* Fix z->range: */
        z->range = iv_meet(iv_sqrt(x->range), z->range);
        #endif
        return(z);
      }
  }

void foi_sqrt_approx(
    Interval xr,
    FloatP alphap,
    FloatP betap,
    FloatP gammap
  )
  /* Computes chebyshev approx to $\sqrt{x}$ in $xr$. */
  {
    Float sa, sb, axa, axb, da, db, dlo, dm;

    #if FOI_TRACE_SQRT
        fprintf(stderr, "  enter foi_sqrt_approx:\n");
        fprintf(stderr, "    xr.lo = %15.8e\n", xr.lo);
        fprintf(stderr, "    xr.hi = %15.8e\n", xr.hi);
    #endif

    /* I believe these computations can't overflow: */
    if (xr.lo < Zero) xr.lo = Zero;
    ROUND_DOWN;
    sa = sqrt(xr.lo);
    sb = sqrt(xr.hi);
    *alphap = One/(sa + sb);
    ROUND_UP;
    dm = (One/(*alphap))/Four;
    ROUND_UP;
    axa = (*alphap) * xr.lo;
    axb = (*alphap) * xr.hi;
    ROUND_DOWN;
    da = sa - axa;
    db = sb - axb;
    if ((dm < da) || (dm < db))
      { fprintf(stderr, "foi_sqrt: xr.lo = %14.8e\n", xr.lo);
        fprintf(stderr, "foi_sqrt: xr.hi = %14.8e\n", xr.hi);
        fprintf(stderr, "foi_sqrt: sqrt(0.5) = %14.8e\n", sqrt(0.5));
        fprintf(stderr, "foi_sqrt: sa = %14.8e\n", sa);
        fprintf(stderr, "foi_sqrt: sb = %14.8e\n", sb);
        fprintf(stderr, "foi_sqrt: alpha = %14.8e\n", *alphap);
        fprintf(stderr, "foi_sqrt: axa = %14.8e\n", dm);
        fprintf(stderr, "foi_sqrt: axb = %14.8e\n", dm);
        fprintf(stderr, "foi_sqrt: dm = %14.8e\n", dm);
        fprintf(stderr, "foi_sqrt: da = %14.8e\n", da);
        fprintf(stderr, "foi_sqrt: db = %14.8e\n", db);
        error ("foi_sqrt: assertion failure, dm/da/db");
      }
    dlo = FMIN(da, db);
    ROUND_NEAR;
    *betap = (dm + dlo)/Two;
    ROUND_UP;
    *gammap = FMAX(dm - (*betap), (*betap) - dlo);

    #if FOI_TRACE_SQRT
        fprintf(stderr,
          "    alpha = %15.8e  beta  = %15.8e  gamma = %15.8e\n",
          *alphap, *betap, *gammap
        );
        fprintf(stderr, "  exit foi_sqrt_approx.\n");
    #endif
  }

FOIP foi_affine(FOIP x, Float alpha, Float beta, Float gamma)
  {
    if (FOI_ISFULL(x)
    || (FABS(alpha) >= PlusInfinity)
    || (FABS(beta) >= PlusInfinity)
    || (FABS(gamma) >= PlusInfinity)
    )
      { return (FOI_FULL); }
    else if (alpha == Zero)
      { return (foi_const(beta, gamma)); }
    else if (gamma < Zero)
      { error ("foi_affine: gamma is negative.");
        return(FOI_FULL);
      }
    else
      { MemP frame = foi_stack_top;
        FOIP z = foi_alloc_head();
        TermP zp;
        TermCount zn = 0;
        TermP xp = (TermP) (x + 1);
        TermCount xn = x->nterms;
        Float ax;
        Float err = gamma;

        #if FOI_TRACE_AFFINE
            fprintf(stderr, "enter foi_affine:\n");
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  x = %p nterms = %ld\n", x, x->nterms);
            fprintf(stderr,
              "  alpha = %15.8e  beta  = %15.8e  gamma = %15.8e\n",
              alpha, beta, gamma
            );
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif

        ROUND_NEAR;
        flt_mul(x->center, alpha, &ax, &err);
        ROUND_NEAR;
        flt_add(ax, beta, &(z->center), &err);
        if (FABS(z->center) >= Infinity || err >= Infinity)
          { foi_flush(frame); return(FOI_FULL); }

        while (xn > 0)
          {
            zp = foi_alloc_term();
            zp->id = xp->id;
            ROUND_NEAR;
            flt_mul(xp->coef, alpha, &(zp->coef), &err);
            if (FABS(zp->coef) >= Infinity || err >= Infinity)
              { foi_flush(frame); return(FOI_FULL); }
            xp++; xn--;
            if (zp->coef == Zero) foi_pop_term(); else zn++;
          }

        /* Add rounding error term: */
        if (err != Zero)
          {
            zp = foi_alloc_term();
            zp->coef = err;
            zp->id = foi_next_id;
            foi_next_id++;
            zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if MIXED
        /* Compute z->range: */
        z->range = iv_meet(
          iv_affine(x->range, alpha, beta, gamma),
          foi_implicit_range(z)
        );
        #endif

        #if FOI_TRACE_AFFINE
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z = %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit foi_affine.\n");
        #endif

        return(z);
      }
  }

FOIP foi_const(Float c, Float err)
  {
    if ((FABS(c) >= PlusInfinity) || (FABS(err) >= PlusInfinity))
      { return (FOI_FULL); }
    else if (err < Zero)
      { error ("foi_const: err is negative"); /* and that is not human! */
        return (FOI_FULL);
      }
    else
      { FOIP z = foi_alloc_head();

        #if FOI_TRACE_CONST
            fprintf(stderr, "enter foi_const:\n");
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif

        z->center = c;
        if (err == Zero)
          {
            #if MIXED
            z->range.lo = c;
            z->range.hi = c;
            #endif
            z->nterms = 0;
          }
        else
          {
            TermP zp;
            zp = foi_alloc_term();
            zp->coef = err;
            zp->id = foi_next_id;
            foi_next_id++;
            #if MIXED
            ROUND_DOWN;
            z->range.lo = c - err;
            ROUND_UP;
            z->range.hi = c + err;
            #endif
            z->nterms = 1;
          }

        #if FOI_TRACE_CONST
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  z = %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit foi_const.\n");
        #endif

        return (z);
      }
  }

FOIP foi_from_interval(Interval x)
  {
    if (IV_ISFULL(x))
      return (FOI_FULL);
    else
      { FOIP z = foi_alloc_head();

        #if FOI_TRACE_FROM_IV
            fprintf(stderr, "enter foi_from_interval:\n");
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  x.lo = "); iv_print(stderr, x); fprintf(stderr, "\n");
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif

        if (x.lo == x.hi)
          { z->nterms = 0;
            z->center = x.lo;
          }
        else
          { TermP zp;
            Float rlo, rhi;
            z->nterms = 1;

            /* Compute center of interval: */
            ROUND_NEAR;
            z->center = (x.lo * Half) + (x.hi * Half);
            if ((z->center < x.lo) |  (z->center > x.hi))
              error("foi_from_interval: assert failed (z->center)");

            /* Append one noise term for the deviation: */
            zp = foi_alloc_term();
            ROUND_UP;
            rlo = z->center - x.lo;
            rhi = x.hi - z->center;
            if ((rlo < Zero) | (rhi < Zero))
              error("foi_from_interval: assert failed (rlo, rhi)");
            zp->coef = FMAX(rlo, rhi);
            zp->id = foi_next_id;
            foi_next_id++;
          }

        #if MIXED
        /* Set explicit range: */
        z->range = x;
        #endif

        #if FOI_TRACE_FROM_IV
            fprintf(stderr, "  stack_top = %p\n", foi_stack_top);
            fprintf(stderr, "  z = %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit foi_from_interval.\n");
        #endif

        return (z);
      }
  }

