/* FOI routines */

#include "aa.h"
#include <js.h>
#include <flt.h>
#include <ia.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#define QUICKMUL 1

#define AA_TRACE_ADD      0
#define AA_TRACE_MUL      0
#define AA_TRACE_SQR      0
#define AA_TRACE_NEG      0
#define AA_TRACE_INV      0
#define AA_TRACE_SQRT     0
#define AA_TRACE_THROW    0
#define AA_TRACE_AFFINE   0
#define AA_TRACE_CONST    0
#define AA_TRACE_FROM_INT 0

/*** PROTOTYPES FOR INTERNAL PROCEDURES ***/

void aa_mul_naf_range(
    AAP x, 
    AAP y, 
    FloatP zcp, 
    FloatP errp
  );

void aa_sqr_naf_range(
    AAP x, 
    FloatP zcp, 
    FloatP errp
  );

void aa_sqrt_approx(
    Interval xr,
    FloatP alphap,
    FloatP betap,
    FloatP gammap
  );

void aa_alt_sqrt_approx(
    Interval xr,
    FloatP alphap,
    FloatP betap,
    FloatP gammap
  );

void aa_inv_approx(
    Interval xr,
    FloatP alphap,
    FloatP betap,
    FloatP gammap
  );
  
Float aa_throw_coef(void);

/*** ZERO ***/

static AAHead AA_ZERO_BDY = {
      #if MIXED
        {Zero, Zero},
      #endif
      0,
      Zero
    };

#define AA_ZERO ((AAP) &AA_ZERO_BDY)

/*** NOISEVAR ID'S ***/

static VarId aa_next_id = 0; /* next unused noisevar id */

/*** STACK ALLOCATION ***/

#define AA_STACK_SIZE 100000L

static MemP aa_stack_bot;  /* bottom-of-stack pointer. */
static MemP aa_stack_top;  /* The top of the AA stack. */
static MemP aa_stack_lim;  /* limiting pointer for stack memory area. */

MemP aa_top (void)
  { return(aa_stack_top); }

#define AA_VALID_FRAME(f) (((f) <= aa_stack_top) && ((f) >= aa_stack_bot))

#define AA_IN_FRAME(p, f) (((f) <= (p)) && ((p) < aa_stack_top))

void aa_flush (MemP frame)
  { assert(AA_VALID_FRAME(frame), "aa_flush: bad frame pointer");
    aa_stack_top = frame;
  }

AAP aa_return (MemP frame, AAP result)
  { assert(AA_VALID_FRAME(frame), "aa_return: bad frame pointer");
    if (AA_IN_FRAME((MemP) result, frame))
      { MemSize nbytes = sizeof(AAHead);
        nbytes += (result->nterms) * sizeof(AATerm);
        if (((MemP) result) > ((MemP) frame))
          { memcpy((char *) frame, (char *) result, nbytes); }
        aa_stack_top = (MemP)(((char *) frame) + nbytes);
        return((AAP) frame);
      }
    else
      { aa_stack_top = frame;
        return(result);
      }
  }

AAP aa_alloc_head(void)
  {
    MemP old_top = aa_stack_top;
    aa_stack_top = (MemP) (((char *) aa_stack_top) + sizeof(AAHead));
    if (((MemP) aa_stack_top) > ((MemP) aa_stack_lim))
      { fprintf(stderr, "old_top =       %p\n", old_top);
        fprintf(stderr, "aa_stack_top = %p\n", aa_stack_top);
        fprintf(stderr, "aa_stack_lim = %p\n", aa_stack_lim);
        error("aa_alloc_head: not enough space for AA head");
      }
    return (AAP) old_top;
  }

void aa_pop_head(void)
  { aa_stack_top = (MemP)(((char *) aa_stack_top) - sizeof(AAHead)); }

AATermP aa_alloc_term(void)
  {
    MemP old_top = aa_stack_top;
    aa_stack_top = (MemP)(((char *) aa_stack_top) + sizeof(AATerm));
    if (((MemP) aa_stack_top) > ((MemP) aa_stack_lim))
      { fprintf(stderr, "old_top =       %p\n", old_top);
        fprintf(stderr, "aa_stack_top = %p\n", aa_stack_top);
        fprintf(stderr, "aa_stack_lim = %p\n", aa_stack_lim);
        error("aa_alloc_term: not enough space for AA term");
      }
    return (AATermP) old_top;
  }

void aa_pop_term(void)
  { aa_stack_top = (MemP) (((char *) aa_stack_top) - sizeof(AATerm)); }

/*** INITIALIZATION ***/

void aa_init (void)
  {
    long dif;
    aa_stack_bot = (MemP) malloc(AA_STACK_SIZE);
    fprintf(stderr, "aa_stack_bot = %p\n", aa_stack_bot);
    if (aa_stack_bot == (MemP) NULL)
      error("aa_init: heap not allocated");
    aa_stack_lim = (MemP) (((char *) aa_stack_bot) + AA_STACK_SIZE);
    fprintf(stderr, "aa_stack_lim = %p\n", aa_stack_lim);
    dif = ((char *) aa_stack_lim) - ((char *) aa_stack_bot);
    fprintf(stderr, "stack size = %ld bytes\n", dif);
    aa_stack_top = aa_stack_bot;
    fprintf(stderr, "aa_stack_top = %p (%ld)\n", 
      aa_stack_top, (unsigned long) aa_stack_top
    );

    aa_next_id = 0;
  }

/*** HEAP ALLOCATION ***/

AAP aa_heap_alloc (AATermCount n)
  {
    AAP r;
    MemSize nbytes = sizeof(AAHead);
    nbytes += n * sizeof(AATerm);
    r = (AAP) malloc(nbytes);
    if (r == (AAP) NULL)
      { error("out of memory"); }
    return (r);
  }

AAP aa_heap_realloc (AAP x, AATermCount n)
  { AAP r;
    MemSize nbytes = sizeof(AAHead);
    nbytes += n * sizeof(AATerm);
    r = (AAP) realloc((void *) x, nbytes);
    if (r == (AAP) NULL)
      { error("out of memory"); }
    return (r);
  }

void aa_heap_free (AAP x)
  { free((void *) x); }

/*** MISCELLANEOUS ***/

void aa_print (FILE *f, AAP x)
  {
    ia_print (f, aa_range(x));
    if (AA_ISFULL(x))
      { int i;
        for (i=0; i<F_FMT_WIDTH; i++) putc('*', f);
        return;
      }
    else
      { AATermP xp = (AATermP) (x + 1);
        AATermCount xn = x->nterms;

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

void aa_move (AAP source, AAP destination)
  { if ( (source != destination)
      && (source != (AAP) NULL)
      && (! (AA_ISFULL(source)))
    ) { MemSize nbytes = sizeof(AAHead);
        nbytes += (source->nterms) * sizeof(AATerm);
        memcpy((char *) destination, (char *) source, nbytes);
      }
  }

Interval aa_range(AAP x)
  {
    if (AA_ISFULL(x)) return (IA_FULL);
    #if MIXED
      return (x->range);
    #else
      return (aa_implicit_range(x));
    #endif
  }

Interval aa_implicit_range (AAP x)
  {
    if (AA_ISFULL(x))
      return (IA_FULL);
    else
      {
        Interval range;
        Float tdev = aa_sum_abs_terms((AATermP)(x+1), x->nterms);
        if (tdev >= PlusInfinity) return (IA_FULL);
        ROUND_DOWN;
        range.lo = x->center - tdev;
        ROUND_UP;
        range.hi = x->center + tdev;
        IA_NORMFULL(range);
        return (range);
      }
  }

Float aa_sum_abs_terms (AATermP xp, AATermCount n)
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
  
Float aa_throw_coef(void)
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

AAP aa_throw(int nterms)
  {
    int coins;
    coins = random();
    if ((coins&255) == 0)
      return (AA_FULL);
    else if ((coins&63) == 0)
      return (aa_zero());
    else
      {
        MemP frame = aa_stack_top;
        AAP z;
        AATermP zp;
        AATermCount zn = 0;
        Float t;
        int i;

        z = aa_alloc_head();

        #if AA_TRACE_THROW
            fprintf(stderr, "enter aa_throw:\n");
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  nterms = %d\n", nterms);
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif

        ROUND_NEAR;
        z->center = aa_throw_coef();
	if (FABS(z->center) >= Infinity)
	  { aa_flush(frame); return(AA_FULL); }
    
        for (i=0; i<nterms; i++)
          {
            t = aa_throw_coef();
	    if (FABS(t) >= Infinity)
	      { aa_flush(frame); return(AA_FULL); }
            else if (t != Zero)
              {
		zp = aa_alloc_term();
		zp->coef = t;
		zp->id = i;
		zn++;
              }
          }

        /* Store term count: */
        z->nterms = zn;
        
        /* Make sure the noise symbols used do exist: */
        if (nterms > aa_next_id) aa_next_id = nterms;

        #if MIXED
        /* Compute z->range: */
        z->range = aa_implicit_range(z);
        #endif

        #if AA_TRACE_THROW
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit aa_throw.\n");
        #endif

        return (z);
      }
  }

/*** ARITHMETIC ***/

AAP aa_zero(void)
  { return (AA_ZERO); }

AAP aa_add(AAP x, AAP y)
  {
    if (AA_ISFULL(x) || AA_ISFULL(y))
      return (AA_FULL);
    else
      {
        MemP frame = aa_stack_top;
        AAP z = aa_alloc_head();
        AATermP xp = (AATermP) (x + 1);
        AATermP yp = (AATermP) (y + 1);
        AATermP zp;
        Float err = Zero;
        AATermCount xn = x->nterms;
        AATermCount yn = y->nterms;
        AATermCount zn = 0;

        #if AA_TRACE_ADD
            fprintf(stderr, "enter aa_add:\n");
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
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
          { aa_flush(frame); return(AA_FULL); }

        /* Merge noise terms: */
        while (xn > 0 || yn > 0)
          {
            zp = aa_alloc_term();
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
                  { aa_flush(frame); return(AA_FULL); }
                xp++; xn--;
                yp++; yn--;
                if (zp->coef == Zero) aa_pop_term(); else zn++;
              }
          }

        /* Add rounding error term: */
        if (err != Zero)
          {
            zp = aa_alloc_term();
            zp->coef = err;
            zp->id = aa_next_id;
            aa_next_id++;
            zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if MIXED
        /* Compute z->range: */
        z->range = ia_meet(
          ia_add(x->range, y->range),
          aa_implicit_range(z)
        );
        #endif

        #if AA_TRACE_ADD
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit aa_add.\n");
        #endif

        return (z);
      }
  }

AAP aa_sub(AAP x, AAP y)
  {
    if (AA_ISFULL(x) || AA_ISFULL(y))
      return (AA_FULL);
    else
      {
        MemP frame = aa_stack_top;
        AAP z = aa_alloc_head();
        AATermP xp = (AATermP) (x + 1);
        AATermP yp = (AATermP) (y + 1);
        AATermP zp;
        Float err = Zero;
        AATermCount xn = x->nterms;
        AATermCount yn = y->nterms;
        AATermCount zn = 0;

        #if AA_TRACE_SUB
            fprintf(stderr, "enter aa_sub:\n");
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
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
          { aa_flush(frame); return(AA_FULL); }

        /* Merge noise terms: */
        while (xn > 0 || yn > 0)
          {
            zp = aa_alloc_term();
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
                  { aa_flush(frame); return(AA_FULL); }
                xp++; xn--;
                yp++; yn--;
                if (zp->coef == Zero) aa_pop_term(); else zn++;
              }
          }

        /* Add rounding error term: */
        if (err != Zero)
          {
            zp = aa_alloc_term();
            zp->coef = err;
            zp->id = aa_next_id;
            aa_next_id++;
            zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if MIXED
        /* Compute z->range: */
        z->range = ia_meet(
          ia_add(x->range, y->range),
          aa_implicit_range(z)
        );
        #endif

        #if AA_TRACE_SUB
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit aa_sub.\n");
        #endif

        return (z);
      }
  }

AAP aa_neg(AAP x)
  {
    if (AA_ISFULL(x))
      return (AA_FULL);
    else
      {
        AAP z = aa_alloc_head();
        AATermP xp = (AATermP) (x + 1);
        AATermP zp;
        AATermCount xn = x->nterms;
        AATermCount zn = 0;

        #if AA_TRACE_NEG
            fprintf(stderr, "enter aa_neg:\n");
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  x = %p nterms = %ld\n", x, x->nterms);
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif


        /* Negate center value: */
        z->center = - (x->center);

        /* Copy and negate noise terms: */
        while (xn > 0)
          {
            zp = aa_alloc_term();
            zp->id = xp->id;
            zp->coef = - (xp->coef);
            xp++; xn--;
            zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if MIXED
        /* Compute z->range: */
        z->range = ia_neg(x->range);
        #endif

        #if AA_TRACE_NEG
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit aa_neg.\n");
        #endif

        return (z);
      }
  }

AAP aa_mul (AAP x, AAP y)
  {
    if (AA_ISFULL(x) || AA_ISFULL(y))
      return (AA_FULL);
    else
      { MemP frame = aa_stack_top;
        AAP z = aa_alloc_head();
        AATermP zp;
        AATermCount zn = 0;
        Float err = Zero;

        #if AA_TRACE_MUL
            fprintf(stderr, "enter aa_mul:\n");
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
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
          AATermP xp = (AATermP) (x + 1);
          AATermP yp = (AATermP) (y + 1);
          AATermCount xn = x->nterms;
          AATermCount yn = y->nterms;
          Float xt, yt;

          /* Multiply center values: */
          ROUND_NEAR;
          flt_mul(x->center, y->center, &(z->center), &err);
          if (FABS(z->center) >= Infinity || err >= Infinity)
            { aa_flush(frame); return(AA_FULL); }

          /* Merge noise terms: */
          while (xn > 0 || yn > 0)
            {
              zp = aa_alloc_term();
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
                    { aa_flush(frame); return(AA_FULL); }
                  ROUND_NEAR;
                  flt_add(xt, yt, &(zp->coef), &err);
                  if (FABS(zp->coef) >= Infinity || err >= Infinity)
                    { aa_flush(frame); return(AA_FULL); }
                  xp++; xn--;
                  yp++; yn--;
                  if (zp->coef == Zero) aa_pop_term(); else zn++;
                }
            }
        }

        /************************************/
        /* II. Approximate non-affine part. */
        /************************************/

        aa_mul_naf_range(x, y, &(z->center), &err);
        if (err >= PlusInfinity)
          { aa_flush(frame); return (AA_FULL); }

        /* Add rounding error term: */
        if (err != Zero)
          {
            zp = aa_alloc_term();
            zp->coef = err;
            zp->id = aa_next_id;
            aa_next_id++;
            zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if MIXED
        /* Compute z->range: */
        z->range = ia_meet(
          ia_mul(x->range, y->range),
          aa_implicit_range(z)
        );
        #endif

        #if AA_TRACE_MUL
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit aa_mul.\n");
        #endif

        return (z);
      }
  }

void aa_mul_naf_range(AAP x, AAP y, FloatP zcp, FloatP errp)
  {
    #if QUICKMUL
      Float xr = aa_sum_abs_terms((AATermP)(x+1), x->nterms);
      Float yr = aa_sum_abs_terms((AATermP)(y+1), y->nterms);
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

AAP aa_sqr (AAP x)
  {
    if (AA_ISFULL(x))
      return (AA_FULL);
    else
      { MemP frame = aa_stack_top;
        AAP z = aa_alloc_head();
        AATermP zp;
        AATermCount zn = 0;
        Float err = Zero;

        #if AA_TRACE_SQR
            fprintf(stderr, "enter aa_sqr:\n");
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  x = %p nterms = %ld\n", x, x->nterms);
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif

        /*********************************************/
        /* I. Compute the affine part of the square  */
        /*********************************************/

        {
          AATermP xp = (AATermP) (x + 1);
          AATermCount xn = x->nterms;

          /* Square the center value: */
          ROUND_NEAR;
          flt_mul(x->center, x->center, &(z->center), &err);
          if (FABS(z->center) >= Infinity || err >= Infinity)
            { aa_flush(frame); return(AA_FULL); }

          /* Process noise terms: */
          while (xn > 0)
            {
              zp = aa_alloc_term();
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

        aa_sqr_naf_range(x, &(z->center), &err);
        if (err >= PlusInfinity)
          { aa_flush(frame); return (AA_FULL); }

        /* Add rounding error term: */
        if (err != Zero)
          {
            zp = aa_alloc_term();
            zp->coef = err;
            zp->id = aa_next_id;
            aa_next_id++;
            zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if MIXED
        /* Compute z->range: */
        z->range = ia_meet(
          ia_sqr(x->range),
          aa_implicit_range(z)
        );
        #endif

        #if AA_TRACE_SQR
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit aa_sqr.\n");
        #endif

        return (z);
      }
  }

void aa_sqr_naf_range(AAP x, FloatP zcp, FloatP errp)
  {
    #if QUICKMUL
      Float xr = aa_sum_abs_terms((AATermP)(x+1), x->nterms);
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

AAP aa_inv (AAP x)
  {
    if (AA_ISFULL(x))
      { return (AA_FULL); }
    else if (x->nterms == 0)
      {
        Float sx;
        Float err = Zero;
        if (x->center == Zero)
          { error("aa_inv: argument is zero"); }
        flt_inv(x->center, &sx, &err);
        if (FABS(sx) >= Infinity || err >= Infinity)
          { return(AA_FULL); }
        return (aa_const(sx, err));
      }
    else
      {
        Interval xr = aa_range(x);
        Float alpha, beta, gamma;
        AAP z;
        if (xr.lo <= Zero && xr.hi >= Zero)
          { return(AA_FULL); }
        aa_inv_approx(xr, &alpha, &beta,&gamma);
        z = aa_affine(x, alpha, beta, gamma);
        #if MIXED
        /* Fix z->range: */
        z->range = ia_meet(ia_inv(x->range), z->range);
        #endif
        return(z);
      }
  }

void aa_inv_approx(
    Interval xr,
    Float xc,
    FloatP alphap,
    FloatP betap,
    FloatP gammap
  )
  /* Computes chebyshev approx to $1/x$ in $xr$. */
  {
    Float maxabsalpha, axa, axb, da, db, dlo, dhi;
    int negative;

    #if AA_TRACE_INV
        fprintf(stderr, "  enter aa_inv_approx:\n");
        fprintf(stderr, "    xr.lo = %15.8e\n", xr.lo);
        fprintf(stderr, "    xr.hi = %15.8e\n", xr.hi);
    #endif

    /* Assumes $xr$ doesn't contain zero: */
    negative = (xr.lo < Zero);
    if (negative)
      { Float t; t = xr.lo; xr.lo = -xr.hi; xr.hi = -t; xc = -xc; }
    /* Compute approximate slope in interval: */
    ROUND_DOWN;
    *alphap = -((One/xr.hi)/xr.lo);
    maxabsalpha = MaxFloat/Two/xc;
    if (-(*alphap) >= maxabsalpha)
      { *alphap = -(maxabsalpha); }
    /* Compute max and min of $1/x - \alpha x$ in $xr$: */
    ROUND_NEAR;
    if ( -(*alphap)*xr.hi*xr.hi <= One )
      { /* Difference is monotonically decreasing in $xr$. */
        ROUND_UP;
        axb = (*alphap) * xr.hi;
        ROUND_DOWN;
        axa = (*alphap) * xr.lo;
        dlo = One/xr.hi - axb;
        ROUND_UP;
        dhi = One/xr.lo - axa;
        assert((dlo <= dhi), "aa_inv: case 1");
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
        assert((dlo <= dhi), "aa_inv: case 2");
      }
    else
      { /* Difference may be monotonic or concave in $xr$. */
        ROUND_DOWN;
        dlo = Two * sqrt(-(*alphap));
        axa = (*alphap) * xr.lo;
        axb = (*alphap) * xr.hi;
        ROUND_UP;
        da = One/xr.lo - axa;
        db = One/xr.hi - axb;
        dhi = FMAX(da, db);
        assert ((dlo <= dhi), "aa_inv: case 3");
      }
    ROUND_NEAR;
    *betap = (dhi + dlo)/Two;
    ROUND_UP;
    *gammap = FMAX(dhi - (*betap), (*betap) - dlo);
    if (negative) *betap = -(*betap);

    #if AA_TRACE_INV
        fprintf(stderr,
          "    alpha = %15.8e  beta  = %15.8e  gamma = %15.8e\n",
          *alphap, *betap, *gammap
        );
        fprintf(stderr, "  exit aa_inv_approx.\n");
    #endif
  }

AAP aa_sqrt(AAP x)
  {
    if (AA_ISFULL(x))
      { return (AA_FULL); }
    else if (x->nterms == 0)
      {
        Float sx;
        Float err = Zero;
        if (x->center < Zero)
          { error("aa_sqrt: argument is negative"); }
        if (x->center == Zero)
          { return (AA_ZERO); }
        flt_sqrt(x->center, &sx, &err);
        /* Overflow can't possibly have occurred, can it? */
        return (aa_const(sx, err));
      }
    else
      {
        Interval xr = aa_range(x);
        Float alpha, beta, gamma;
        AAP z;
        if (xr.hi < Zero)
          { error("aa_sqrt: argument is negative"); }
        if (xr.hi == Zero)
          { return (AA_ZERO); }
        aa_sqrt_approx(xr, &alpha, &beta,&gamma);
        z = aa_affine(x, alpha, beta, gamma);
        #if MIXED
        /* Fix z->range: */
        z->range = ia_meet(ia_sqrt(x->range), z->range);
        #endif
        return(z);
      }
  }

void aa_sqrt_approx(
    Interval xr,
    FloatP alphap,
    FloatP betap,
    FloatP gammap
  )
  /* Computes chebyshev approx to $\sqrt{x}$ in $xr$. */
  {
    Float sa, sb, axa, axb, da, db, dlo, dm;

    #if AA_TRACE_SQRT
        fprintf(stderr, "  enter aa_sqrt_approx:\n");
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
      { fprintf(stderr, "aa_sqrt: xr.lo = %14.8e\n", xr.lo);
        fprintf(stderr, "aa_sqrt: xr.hi = %14.8e\n", xr.hi);
        fprintf(stderr, "aa_sqrt: sqrt(0.5) = %14.8e\n", sqrt(0.5));
        fprintf(stderr, "aa_sqrt: sa = %14.8e\n", sa);
        fprintf(stderr, "aa_sqrt: sb = %14.8e\n", sb);
        fprintf(stderr, "aa_sqrt: alpha = %14.8e\n", *alphap);
        fprintf(stderr, "aa_sqrt: axa = %14.8e\n", dm);
        fprintf(stderr, "aa_sqrt: axb = %14.8e\n", dm);
        fprintf(stderr, "aa_sqrt: dm = %14.8e\n", dm);
        fprintf(stderr, "aa_sqrt: da = %14.8e\n", da);
        fprintf(stderr, "aa_sqrt: db = %14.8e\n", db);
        error ("aa_sqrt: assertion failure, dm/da/db");
      }
    dlo = FMIN(da, db);
    ROUND_NEAR;
    *betap = (dm + dlo)/Two;
    ROUND_UP;
    *gammap = FMAX(dm - (*betap), (*betap) - dlo);

    #if AA_TRACE_SQRT
        fprintf(stderr,
          "    alpha = %15.8e  beta  = %15.8e  gamma = %15.8e\n",
          *alphap, *betap, *gammap
        );
        fprintf(stderr, "  exit aa_sqrt_approx.\n");
    #endif
  }

AAP aa_aux_scale(AAP x, Float alpha, FloatP err)
  /*
    Pushes an affine form equal to $\alpha x$ on top of the AA stack,
    and returns its address.  
    
    The form does not include the error term; instead, any rounding errors 
    are added to the $*err$ parameter.  
    
    NOTE: The range of the result is NOT defined by this procedure.

    The form is created on the stack even if it is identically zero.
    On the other hand, if overflow happens in any coefficient of the
    result, or in the $*err$ variable, the routine resets the stack,
    stores PlusInfinity in $*err$, and returns $AA_FULL$. */
  {
    if (AA_ISFULL(x)
    || (FABS(alpha) >= PlusInfinity)
    )
      { return (AA_FULL); }
    else
      { MemP frame = aa_stack_top;
        AAP z = aa_alloc_head();
        AATermP zp;
        AATermCount zn = 0;
        AATermP xp = (AATermP) (x + 1);
        AATermCount xn = x->nterms;
        Float ax;
        Float err = Zero;

        #if AA_TRACE_SCALE
            fprintf(stderr, "enter aa_aux_scale:\n");
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  x = %p nterms = %ld\n", x, x->nterms);
            fprintf(stderr, "  alpha = %15.8e\n", alpha);
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif

        ROUND_NEAR;
        flt_mul(x->center, alpha, &(z->center), &err);
        if (FABS(z->center) >= Infinity || err >= Infinity)
          { aa_flush(frame); return(AA_FULL); }

        while (xn > 0)
          {
            zp = aa_alloc_term();
            zp->id = xp->id;
            ROUND_NEAR;
            flt_mul(xp->coef, alpha, &(zp->coef), &err);
            if (FABS(zp->coef) >= Infinity || err >= Infinity)
              { aa_flush(frame); return(AA_FULL); }
            xp++; xn--;
            if (zp->coef == Zero) aa_pop_term(); else zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if AA_TRACE_SCALE
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z = %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit aa_aux_scale.\n");
        #endif

        return(z);
      }
  }

AAP aa_affine(AAP x, Float alpha, Float beta, Float gamma)
  {
    if (AA_ISFULL(x)
    || (FABS(alpha) >= PlusInfinity)
    || (FABS(beta) >= PlusInfinity)
    || (FABS(gamma) >= PlusInfinity)
    )
      { return (AA_FULL); }
    else if (alpha == Zero)
      { return (aa_const(beta, gamma)); }
    else if (gamma < Zero)
      { error ("aa_affine: gamma is negative.");
        return(AA_FULL);
      }
    else
      { MemP frame = aa_stack_top;
        AAP z = aa_alloc_head();
        AATermP zp;
        AATermCount zn = 0;
        AATermP xp = (AATermP) (x + 1);
        AATermCount xn = x->nterms;
        Float ax;
        Float err = gamma;

        #if AA_TRACE_AFFINE
            fprintf(stderr, "enter aa_affine:\n");
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
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
          { aa_flush(frame); return(AA_FULL); }

        while (xn > 0)
          {
            zp = aa_alloc_term();
            zp->id = xp->id;
            ROUND_NEAR;
            flt_mul(xp->coef, alpha, &(zp->coef), &err);
            if (FABS(zp->coef) >= Infinity || err >= Infinity)
              { aa_flush(frame); return(AA_FULL); }
            xp++; xn--;
            if (zp->coef == Zero) aa_pop_term(); else zn++;
          }

        /* Add rounding error term: */
        if (err != Zero)
          {
            zp = aa_alloc_term();
            zp->coef = err;
            zp->id = aa_next_id;
            aa_next_id++;
            zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if MIXED
        /* Compute z->range: */
        z->range = ia_meet(
          ia_affine(x->range, alpha, beta, gamma),
          aa_implicit_range(z)
        );
        #endif

        #if AA_TRACE_AFFINE
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z = %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit aa_affine.\n");
        #endif

        return(z);
      }
  }

AAP aa_div_affine(AAP x, Float zeta, Float beta, Float gamma)
  {
    if (AA_ISFULL(x)
    || (FABS(zeta) <= Zero)
    || (FABS(beta) >= PlusInfinity)
    || (FABS(gamma) >= PlusInfinity)
    )
      { return (AA_FULL); }
    else if (gamma < Zero)
      { error ("aa_affine: gamma is negative.");
        return(AA_FULL);
      }
    else
      { MemP frame = aa_stack_top;
        AAP z = aa_alloc_head();
        AATermP zp;
        AATermCount zn = 0;
        AATermP xp = (AATermP) (x + 1);
        AATermCount xn = x->nterms;
        Float ax;
        Float err = gamma;

        #if AA_TRACE_AFFINE
            fprintf(stderr, "enter aa_div_affine:\n");
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  x = %p nterms = %ld\n", x, x->nterms);
            fprintf(stderr,
              "  zeta = %15.8e  beta  = %15.8e  gamma = %15.8e\n",
              zeta, beta, gamma
            );
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif

        ROUND_NEAR;
        flt_div(x->center, zeta, &ax, &err);
        ROUND_NEAR;
        flt_add(ax, beta, &(z->center), &err);
        if (FABS(z->center) >= Infinity || err >= Infinity)
          { aa_flush(frame); return(AA_FULL); }

        while (xn > 0)
          {
            zp = aa_alloc_term();
            zp->id = xp->id;
            ROUND_NEAR;
            flt_div(xp->coef, zeta, &(zp->coef), &err);
            if (FABS(zp->coef) >= Infinity || err >= Infinity)
              { aa_flush(frame); return(AA_FULL); }
            xp++; xn--;
            if (zp->coef == Zero) aa_pop_term(); else zn++;
          }

        /* Add rounding error term: */
        if (err != Zero)
          {
            zp = aa_alloc_term();
            zp->coef = err;
            zp->id = aa_next_id;
            aa_next_id++;
            zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if MIXED
        /* Compute z->range: */
        z->range = ia_meet(
          ia_div_affine(x->range, alpha, beta, gamma),
          aa_implicit_range(z)
        );
        #endif

        #if AA_TRACE_AFFINE
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z = %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit aa_div_affine.\n");
        #endif

        return(z);
      }
  }

AAP aa_const(Float c, Float err)
  {
    if ((FABS(c) >= PlusInfinity) || (FABS(err) >= PlusInfinity))
      { return (AA_FULL); }
    else if (err < Zero)
      { error ("aa_const: err is negative"); /* and that is not human! */
        return (AA_FULL);
      }
    else
      { AAP z = aa_alloc_head();

        #if AA_TRACE_CONST
            fprintf(stderr, "enter aa_const:\n");
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
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
            AATermP zp;
            zp = aa_alloc_term();
            zp->coef = err;
            zp->id = aa_next_id;
            aa_next_id++;
            #if MIXED
            ROUND_DOWN;
            z->range.lo = c - err;
            ROUND_UP;
            z->range.hi = c + err;
            #endif
            z->nterms = 1;
          }

        #if AA_TRACE_CONST
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  z = %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit aa_const.\n");
        #endif

        return (z);
      }
  }

AAP aa_from_interval(Interval x)
  {
    if (IA_ISFULL(x))
      return (AA_FULL);
    else
      { AAP z = aa_alloc_head();

        #if AA_TRACE_FROM_INT
            fprintf(stderr, "enter aa_from_interval:\n");
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  x.lo = "); ia_print(stderr, x); fprintf(stderr, "\n");
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif

        if (x.lo == x.hi)
          { z->nterms = 0;
            z->center = x.lo;
          }
        else
          { AATermP zp;
            Float rlo, rhi;
            z->nterms = 1;

            /* Compute center of interval: */
            ROUND_NEAR;
            z->center = (x.lo * Half) + (x.hi * Half);
            if ((z->center < x.lo) |  (z->center > x.hi))
              error("aa_from_interval: assert failed (z->center)");

            /* Append one noise term for the deviation: */
            zp = aa_alloc_term();
            ROUND_UP;
            rlo = z->center - x.lo;
            rhi = x.hi - z->center;
            if ((rlo < Zero) | (rhi < Zero))
              error("aa_from_interval: assert failed (rlo, rhi)");
            zp->coef = FMAX(rlo, rhi);
            zp->id = aa_next_id;
            aa_next_id++;
          }

        #if MIXED
        /* Set explicit range: */
        z->range = x;
        #endif

        #if AA_TRACE_FROM_INT
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  z = %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit aa_from_interval.\n");
        #endif

        return (z);
      }
  }

AAP aa_fix_eps(AAP x, AATermCount neps, AATerm eps[])
  {
    if (AA_ISFULL(x))
      return (AA_FULL);
    else
      {
        MemP frame = aa_stack_top;
        AAP z = aa_alloc_head();
        AATermP xp = (AATermP) (x + 1);
        AATermP ep = (AATermP) eps;
        AATermP zp;
        Float err = Zero;
        Float tmp;
        AATermCount xn = x->nterms;
        AATermCount en = neps;
        AATermCount zn = 0;

        #if AA_TRACE_ADD
            fprintf(stderr, "enter aa_fix_eps:\n");
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  x =   %p nterms = %ld\n", x, x->nterms);
            fprintf(stderr, "  eps = %p neps =   %ld\n", eps, neps);
            fprintf(stderr, "  z = %p\n", z);
            fprintf(stderr, "  ...\n");
        #endif

        z->center = x->center;

        /* Merge noise terms, and fix matches: */
        while (xn > 0)
          {
            if ((en == 0) || (xp->id < ep->id))
              { zp = aa_alloc_term();
                *zp = *xp;
                xp++; xn--; zn++; 
              }
            else if ((en != 0) && (ep->id < xp->id))
              { ep++; en--; }
            else
              {
                ROUND_NEAR;
                flt_mul(xp->coef, ep->coef, &tmp, &err);
                flt_add(z->center, tmp, &(z->center), &err);
                if (FABS(z->center) >= Infinity || err >= Infinity)
                  { aa_flush(frame); return(AA_FULL); }
                xp++; xn--;
                ep++; en--;
              }
          }

        /* Add rounding error term: */
        if (err != Zero)
          {
            zp = aa_alloc_term();
            zp->coef = err;
            zp->id = aa_next_id;
            aa_next_id++;
            zn++;
          }

        /* Store term count: */
        z->nterms = zn;

        #if MIXED
        /* Compute z->range: */
        z->range = aa_implicit_range(z);
        #endif

        #if AA_TRACE_ADD
            fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
            fprintf(stderr, "  frame = %p\n", frame);
            fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
            fprintf(stderr, "exit aa_fix_eps.\n");
        #endif

        return (z);
      }
  }

void aa_collapse_pair(AAP x, AAP y, AAP *xr, AAP *yr)
  {
    error("aa_collapse_pair: not implemented yet!");
  }

void aa_alt_sqrt_approx(
    Interval xr,
    FloatP zetap,
    FloatP betap,
    FloatP gammap
  )
  /* Computes chebyshev approx to $\sqrt{x}$ in $xr$. */
  {
    Float ra, rb, da, db, axa, axb, dlo, dhi;

    #if AA_TRACE_SQRT
	fprintf(stderr, "  enter aa_alt_sqrt_approx:\n");
	fprintf(stderr, "    xr.lo = %15.8e\n", xr.lo);
	fprintf(stderr, "    xr.hi = %15.8e\n", xr.hi);
    #endif

    /* I believe these computations can't overflow: */
    if (xr.lo < Zero) xr.lo = Zero;
    ROUND_DOWN;
    ra = sqrt(xr.lo);
    ROUND_UP;
    rb = sqrt(xr.hi);
    ROUND_NEAR;
    (*zetap) = ra + rb;
    ROUND_UP;
    axa = ra*ra/(*zetap);
    axb = rb*rb/(*zetap);
    ROUND_DOWN;
    da = ra - axa;
    db = rb - axb;
    dlo = FMIN(da, db);
    ROUND_UP;
    dhi = (*zetap)/Four;
    ROUND_NEAR;
    (*betap) = (dhi + dlo)/Two;
    ROUND_UP;
    (*gammap) = FMAX(dhi - (*betap), (*betap) - dlo);
    #if AA_TRACE_SQRT
	fprintf(stderr,
	  "    zeta = %15.8e  beta  = %15.8e  gamma = %15.8e\n",
	  *zetap, *betap, *gammap
	);
	fprintf(stderr, "  exit aa_alt_sqrt_approx.\n");
    #endif
  }

AAP aa_alt_sqrt(AAP x)
  {
    if (AA_ISFULL(x))
      { return (AA_FULL); }
    else if (x->nterms == 0)
      {
        Float sx;
        Float err = Zero;
        if (x->center < Zero)
          { error("aa_alt_sqrt: argument is negative"); }
        if (x->center == Zero)
          { return (AA_ZERO); }
        flt_sqrt(x->center, &sx, &err);
        /* Overflow can't possibly have occurred, can it? */
        return (aa_const(sx, err));
      }
    else
      {
        Interval xr = aa_range(x);
        Float zeta, beta, gamma;
        AAP z;
        if (xr.hi < Zero)
          { error("aa_alt_sqrt: argument is negative"); }
        if (xr.hi == Zero)
          { return (AA_ZERO); }
        aa_alt_sqrt_approx(xr, &zeta, &beta, &gamma);
        z = aa_div_affine(x, zeta, beta, gamma);
        #if MIXED
        /* Fix z->range: */
        z->range = ia_meet(ia_sqrt(x->range), z->range);
        #endif
        return(z);
      }
  }  
