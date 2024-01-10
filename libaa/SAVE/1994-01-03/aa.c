/* FOI routines */

#include "aa.h"
#include <js.h>
#include <flt.h>
#include <ia.h>
#include <ioprotos.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>

#define QUICKMUL 0

#define AA_TRACE_ADD      0
#define AA_TRACE_MUL      0
#define AA_TRACE_SQR      0
#define AA_TRACE_NEG      0
#define AA_TRACE_INV      0
#define AA_TRACE_SQRT     0
#define AA_TRACE_ABS      0
#define AA_TRACE_FIX      0
#define AA_TRACE_THROW    0
#define AA_TRACE_SCALE    0
#define AA_TRACE_MIX      0
#define AA_TRACE_AFFINE   0
#define AA_TRACE_JOIN     0
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
    Float xc,
    FloatP alphap,
    FloatP betap,
    FloatP gammap
  );
  
void aa_abs_approx(
    Interval xr,
    FloatP alphap,
    FloatP gammap,
    FloatP deltap
  );
  
void aa_copy_terms(
    AATermP xp, AATermCount xn, 
    AATermCount *znp
  );
  
void aa_neg_terms(
    AATermP xp, AATermCount xn, 
    AATermCount *znp 
  );
  
void aa_scale_terms(
    AATermP xp, AATermCount xn, 
    Float alpha,  
    Float zeta,
    AATermCount *znp, 
    Float *errp
  );
  
void aa_mix_terms(
    AATermP xp, AATermCount xn, Float alpha, 
    AATermP yp, AATermCount yn, Float beta, 
    Float zeta,
    AATermCount *znp, 
    Float *errp
  );
  
void aa_fuzzy_rescale_terms(
    AATermP zp,
    AATermCount *znp, 
    Float lambda, Float sigma, 
    FloatP err
  );
  
Float aa_throw_coef(void);

/*** ANYTHING ***/

#define AA_FULL ((AAP) -1)

/*** ZERO ***/

static AAHead AA_ZERO_BODY = {
      #if MIXED
        {Zero, Zero},
      #endif
      0,
      Zero
    };

#define AA_ZERO ((AAP) &AA_ZERO_BODY)

/*** PREDICATES ***/

#define AA_ISFULL(x) ((x) == AA_FULL)
#define AA_ISZERO(x) (((x)->center == Zero) && ((x)->nterms == 0))

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

AATermP aa_push_term(Float coef, VarId id)
  {
    AATermP zp = aa_alloc_term();
    zp->coef = coef;
    zp->id = id;
    return(zp);
  }
  
void aa_append_error_term (AATermCount *znp, Float err)
  { if (err != Zero)
      { AATermP zp = aa_alloc_term();
        zp->coef = err;
        zp->id = aa_next_id;
        aa_next_id++;
        (*znp)++;
      }
  }

/*** INITIALIZATION ***/

void aa_init (void)
  {
    long dif;
    
    flt_init();
    ia_init();
    
    aa_stack_bot = (MemP) malloc(AA_STACK_SIZE);
    /* fprintf(stderr, "aa_stack_bot = %p\n", aa_stack_bot); */
    if (aa_stack_bot == (MemP) NULL)
      error("aa_init: heap not allocated");
    aa_stack_lim = (MemP) (((char *) aa_stack_bot) + AA_STACK_SIZE);
    /* fprintf(stderr, "aa_stack_lim = %p\n", aa_stack_lim); */
    dif = ((char *) aa_stack_lim) - ((char *) aa_stack_bot);
    /* fprintf(stderr, "stack size = %ld bytes\n", dif); */
    aa_stack_top = aa_stack_bot;
    /* fprintf(stderr, "aa_stack_top = %p (%ld)\n", */
    /*   aa_stack_top, (unsigned long) aa_stack_top */
    /* );                                           */

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

/*** CONSTANTS AND PREDICATES ***/

AAP aa_zero(void)
  { return (AA_ZERO); }

int aa_is_zero (AAP x)
  { return (AA_ISZERO(x)); }

AAP aa_full(void)
  { return (AA_FULL); }

int aa_is_full (AAP x)
  { return (AA_ISFULL(x)); }

/*** MISCELLANEOUS ***/

void aa_print (FILE *f, AAP x)
  {
    ia_print (f, aa_range(x));
    if (AA_ISFULL(x))
      { int i;
        putc(' ', f);
        for (i=1; i<FLT_FMT_WIDTH; i++) putc('*', f);
        return;
      }
    else
      { AATermP xp = (AATermP) (x + 1);
        AATermCount xn = x->nterms;

        flt_print (f, x->center);

        while (xn > 0)
          {
            if (xp->coef >= Zero)
              fprintf(f, " +");
            else
              fprintf(f, " -");
            flt_print(f, FABS(xp->coef)) ;
            fprintf(f, "(%lu)", xp->id);
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
    if (AA_ISFULL(x)) return (ia_full());
    #if MIXED
      return (x->range);
    #else
      return (aa_implicit_range(x));
    #endif
  }

Interval aa_implicit_range (AAP x)
  {
    if (AA_ISFULL(x))
      return (ia_full());
    else
      {
        Interval range;
        Float tdev = aa_sum_abs_terms((AATermP)(x+1), x->nterms);
        if (tdev >= PlusInfinity) return (ia_full());
        ROUND_DOWN;
        range.lo = x->center - tdev;
        ROUND_UP;
        range.hi = x->center + tdev;
        ia_norm_full(&range);
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
  
Float aa_max_abs_term (AATermP xp, AATermCount n)
  {
    Float tdev = Zero;
    while (n > 0)
      {
        Float ti = FABS(xp->coef);
        if (ti > tdev) tdev = ti;
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

/*** INTERNAL ARITHMETIC ROUTINES ***/

/* 
  The routines in this section accept and return sequences of
  AATerm's. 
  
  The input AATerm sequences are given by their addresses $xp, yp$
  and term counts $xn, yn$.  The resulting sequence is pushed on top
  of the AA stack, and its length is stored into the parameter $*znp$.

  The resuting sequence does not include the extra rounding-error term; 
  instead, any rounding errors are added to the $*errp$ parameter.

  If overflow happens in any coefficient of the
  result, or in the $*errp$ variable, the procedures reset the stack,
  store PlusInfinity in $*errp$, and leave $*znp$ unchanged. */
  
void aa_copy_terms(
    AATermP xp, AATermCount xn, 
    AATermCount *znp 
  )
  /* Copies the given terms: $zp[i] = xp[i]$. */
  {
    AATermP zp;
    (*znp) = xn;
    while (xn>0)
      { zp = aa_alloc_term();
        zp->id = xp->id;
        (zp->coef) = (xp->coef);
        xp++; xn--;
      }
  }

void aa_neg_terms(
    AATermP xp, AATermCount xn, 
    AATermCount *znp 
  )
  /* Negates the given terms: $zp[i] = - xp[i]$. */
  {
    AATermP zp;
    (*znp) = xn;
    while (xn>0)
      { zp = aa_alloc_term();
        zp->id = xp->id;
        (zp->coef) = -(xp->coef);
        xp++; xn--;
      }
  }

void aa_scale_terms(
    AATermP xp, AATermCount xn, 
    Float alpha,  
    Float zeta,
    AATermCount *znp, 
    Float *errp
  )
  /* Returns $zp[i] = \alpha xp[i] / \zeta$. */
  {
    #if AA_TRACE_SCALE
      fprintf(stderr, "enter aa_scale_terms:\n");
      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
      fprintf(stderr, "  zn = %ld\n", (*znp));
      fprintf(stderr, "  xp = %p  xn = %ld  alpha = %15.8e\n ", xp, xn, alpha);
      fprintf(stderr, "  zeta  = %15.8e\n ", zeta);
      fprintf(stderr, "  ...\n");
    #endif

    (*znp) = 0;
    if (
      (FABS(alpha) >= PlusInfinity) ||
      (FABS(zeta) <= Zero) ||
      ((*errp) >= PlusInfinity)
    )
      { (*errp) = PlusInfinity; }
    else if ((FABS(alpha) == Zero) || (FABS(zeta) >= PlusInfinity))
      { /* Nothing to do */ }
    else if (alpha == zeta)
      { aa_copy_terms(xp, xn, znp); }
    else if (alpha == -zeta) 
      { aa_neg_terms(xp, xn, znp); }
    else
      { MemP frame = aa_stack_top;
        AATermCount zn = 0;
        AATermP zp;

        while (xn > 0)
          { zp = aa_alloc_term();
            zp->id = xp->id;
            ROUND_NEAR;
            flt_scale(xp->coef, alpha, zeta, &(zp->coef), errp);
            if (FABS(zp->coef) >= PlusInfinity || (*errp) >= PlusInfinity)
              { (*errp) = PlusInfinity; aa_flush(frame); return; }
            xp++; xn--;
            if (zp->coef == Zero) aa_pop_term(); else zn++;
          }

        (*znp) = zn;
      }
    #if AA_TRACE_SCALE
      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
      fprintf(stderr, "  zn = %ld\n", (*znp));
      fprintf(stderr, "exit aa_scale_terms.\n");
    #endif
  }

void aa_mix_terms(
    AATermP xp, AATermCount xn, Float alpha, 
    AATermP yp, AATermCount yn, Float beta, 
    Float zeta,
    AATermCount *znp, 
    Float *errp
  )
  /* Returns $zp[i] = (\alpha xp[i] + \beta yp[i])/\zeta$. */
  {
    #if AA_TRACE_MIX
      fprintf(stderr, "enter aa_mix_terms:\n");
      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
      fprintf(stderr, "  zn = %ld\n", (*znp));
      fprintf(stderr, "  xp = %p  xn = %ld  alpha = %15.8e\n ", xp, xn, alpha);
      fprintf(stderr, "  yp = %p  yn = %ld  beta  = %15.8e\n ", yp, yn, beta);
      fprintf(stderr, "  zeta  = %15.8e\n ", zeta);
      fprintf(stderr, "  ...\n");
    #endif

    (*znp) = 0;
    if (
      (FABS(alpha) >= PlusInfinity) ||
      (FABS(beta) >= PlusInfinity) ||
      (FABS(zeta) == Zero) ||
      ((*errp) >= PlusInfinity)
    )
      { (*errp) = PlusInfinity; }
    else if (FABS(zeta) >= PlusInfinity)
      { /* Nothing to do */ }
    else if ((xn == 0) || (FABS(alpha) == Zero))
      { aa_scale_terms (yp, yn, beta, zeta, znp, errp); }
    else if ((yn == 0) || (FABS(beta) == Zero))
      { aa_scale_terms (xp, xn, alpha, zeta, znp, errp); }
    else
      { MemP frame = aa_stack_top;
        AATermCount zn = 0;
        AATermP zp;
        Float xt, yt;

        while (xn > 0 || yn > 0)
          { zp = aa_alloc_term();
            /* Pick next term:*/
            if (yn == 0) 
              { zp->id = xp-> id; }
            else if (xn == 0)
              { zp->id = yp->id; }
            else if (xp->id < yp->id)
              { zp->id = xp->id; }
            else
              { zp->id = yp->id; }
              
            /* Compute $x$'s contribution: */
            if ((xn != 0) && (zp->id == xp->id))
              { xt = xp->coef; xp++; xn--; }
            else
              { xt = Zero; }
            
            /* Compute $y$'s contribution: */
            if ((yn != 0) && (zp->id == yp->id))
              { yt = yp->coef; yp++; yn--; }
            else
              { yt = Zero; }
              
            /* Mix the two terms: */
            ROUND_NEAR; flt_mix(xt, alpha, yt, beta, zeta, &(zp->coef), errp);
            if (FABS(zp->coef) >= Infinity || (*errp) >= Infinity)
              { (*errp) = PlusInfinity; aa_flush(frame); return; }
            if ((zp->coef) == Zero) aa_pop_term(); else zn++;
          }
        (*znp) = zn;
      }

    #if AA_TRACE_MIX
      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
      fprintf(stderr, "  zn = %ld\n", (*znp));
      fprintf(stderr, "exit aa_mix_terms.\n");
    #endif
  }

/*** ARITHMETIC ***/

AAP aa_affine(AAP x, Float alpha, Float zeta, Float gamma, Float delta)
  {
    #if AA_TRACE_AFFINE
      fprintf(stderr, "enter aa_affine\n");
      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
      fprintf(stderr, "  x = %p  nterms = %ld  alpha = %15.8e\n", x, x->nterms, alpha);
      fprintf(stderr, "  zeta  = %15.8e\n ", zeta);
      fprintf(stderr, "  gamma = %15.8e  delta = %15.8e\n", gamma, delta);
      fprintf(stderr, "  ...\n");
    #endif

    delta = FABS(delta);
    if (
      (AA_ISFULL(x)) ||
      (FABS(alpha) >= PlusInfinity) ||
      (FABS(zeta) == Zero) || 
      (FABS(gamma) >= PlusInfinity) ||
      (FABS(delta) >= PlusInfinity)
    )
      { return (AA_FULL); }
    else if ((FABS(alpha) == Zero) || (FABS(zeta) == PlusInfinity))
      { return (aa_const(gamma, delta)); }
    else
      { MemP frame = aa_stack_top;
        Float err = delta;
        AAP z = aa_alloc_head();

        /* Combine center values: */
        {
          Float zt;
          
          if (alpha == zeta)
            { zt = (x->center); }
          else if (alpha == -zeta)
            { zt = -(x->center); }
          else
            { ROUND_NEAR; flt_scale((x->center), alpha, zeta, &zt, &err); }
          
          if (gamma != Zero)
            { ROUND_NEAR; flt_add(zt, gamma, &zt, &err); }
          
          if ((FABS(zt) >= PlusInfinity) || (err >= PlusInfinity))
            { aa_flush(frame); return(AA_FULL); }
          
          (z->center) = zt;
        }
          
        /* Merge noise terms: */
        aa_scale_terms (
          (AATermP) (x + 1), (x->nterms), alpha,
          zeta,
          &(z->nterms),
          &err
        );
        if (err >= PlusInfinity) 
          { aa_flush(frame); return(AA_FULL); }
        
        /* Add error term: */
        aa_append_error_term(&(z->nterms), err);

        /* check for zero: */
        if (AA_ISZERO(z))
          { aa_flush(frame); return(AA_ZERO); }

        #if MIXED
        /* Compute z->range: */
        z->range = ia_meet(
          ia_affine(x->range, alpha, zeta, gamma, delta),
          aa_implicit_range(z)
        );
        #endif

        #if AA_TRACE_AFFINE
	  fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
	  fprintf(stderr, "  frame = %p\n", frame);
	  fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
	  fprintf(stderr, "exit aa_affine.\n");
        #endif

        return(z);
      }
  }

AAP aa_affine_2(
    AAP x, Float alpha, 
    AAP y, Float beta, 
    Float zeta, 
    Float gamma, 
    Float delta
  )
  {
    #if AA_TRACE_AFFINE
      fprintf(stderr, "enter aa_affine_2\n");
      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
      fprintf(stderr, "  x = %p  nterms = %ld  alpha = %15.8e\n", x, x->nterms, alpha);
      fprintf(stderr, "  y = %p  nterms = %ld  beta  = %15.8e\n", y, y->nterms, beta);
      fprintf(stderr, "  zeta  = %15.8e\n ", zeta);
      fprintf(stderr, "  gamma = %15.8e  delta = %15.8e\n", gamma, delta);
      fprintf(stderr, "  ...\n");
    #endif

    delta = FABS(delta);
    if (
      AA_ISFULL(x) || 
      AA_ISFULL(y) ||
      (FABS(alpha) >= PlusInfinity) ||
      (FABS(beta)  >= PlusInfinity) ||
      (FABS(zeta) == Zero) || 
      (FABS(gamma) >= PlusInfinity) ||
      (FABS(delta) >= PlusInfinity)
    )
      { return(AA_FULL); }
    else if (FABS(zeta) >= PlusInfinity)
      { return(aa_const(gamma, delta)); }
    else if (alpha == Zero)
      { return(aa_affine(y, beta, zeta, gamma, delta)); }
    else if (beta == Zero)
      { return(aa_affine(x, alpha, zeta, gamma, delta)); }
    else
      { MemP frame = aa_stack_top;
        Float err = delta;
        AAP z = aa_alloc_head();

        /* Combine center values: */
        { Float zt;
          
          ROUND_NEAR; 
          flt_mix((x->center), alpha, (y->center), beta, zeta, &zt, &err);

          if (gamma != Zero)
            { ROUND_NEAR; flt_add(zt, gamma, &zt, &err); }

          if ((FABS(zt) >= PlusInfinity) || (err >= PlusInfinity))
            { aa_flush(frame); return(AA_FULL); }
          
          (z->center) = zt;
        }
          
        /* Merge noise terms: */
        aa_mix_terms (
          (AATermP) (x + 1), (x->nterms), alpha,
          (AATermP) (y + 1), (y->nterms), beta,
          zeta,
          &(z->nterms),
          &err
        );
        if (err >= PlusInfinity) 
          { aa_flush(frame); return(AA_FULL); }
        
        /* Add error term: */
        aa_append_error_term(&(z->nterms), err);

        /* Check for zero: */
        if (AA_ISZERO(z))
          { aa_flush(frame); return(AA_ZERO); }

        #if MIXED
        /* Compute z->range: */
        z->range = ia_meet(
          ia_affine_2(x->range, alpha, y->range, beta, zeta, gamma, delta),
          aa_implicit_range(z)
        );
        #endif

        #if AA_TRACE_AFFINE
	  fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
	  fprintf(stderr, "  frame = %p\n", frame);
	  fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
	  fprintf(stderr, "exit aa_affine_2.\n");
        #endif

        return (z);
      }
  }

AAP aa_add(AAP x, AAP y)
  {
    return(aa_affine_2(x, One, y, One, One,  Zero, Zero));
  }

AAP aa_sub(AAP x, AAP y)
  {
    return(aa_affine_2(x, One, y, -One, One,  Zero, Zero));
  }

AAP aa_neg(AAP x)
  {
    return(aa_affine(x, -One, One,  Zero, Zero));
  };

AAP aa_scale(AAP x, Float alpha, Float zeta)
  {
    return(aa_affine(x, alpha, zeta, Zero, Zero));
  }

AAP aa_shift(AAP x, Float gamma)
  {
    return(aa_affine(x, One, One, gamma, Zero));
  }

AAP aa_mul (AAP x, AAP y)
  {
    #if AA_TRACE_MUL
      fprintf(stderr, "enter aa_mul:\n");
      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
      fprintf(stderr, "  x = %p nterms = %ld\n", x, x->nterms);
      fprintf(stderr, "  y = %p nterms = %ld\n", y, y->nterms);
      fprintf(stderr, "  ...\n");
    #endif

    if (AA_ISFULL(x) || AA_ISFULL(y))
      return (AA_FULL);
    else
      { MemP frame = aa_stack_top;
        AAP z = aa_alloc_head();
        Float err = Zero;

        /* Multiply center values: */
        ROUND_NEAR; flt_mul(x->center, y->center, &(z->center), &err);
        if (FABS(z->center) >= Infinity || err >= Infinity)
          { aa_flush(frame); return(AA_FULL); }


        /* Compute affine terms: */
        aa_mix_terms(
          (AATermP) (x + 1), (x->nterms), (y->center),
          (AATermP) (y + 1), (y->nterms), (x->center),
          One,
          &(z->nterms), 
          &err
        );

        if (err >= PlusInfinity)
          { aa_flush(frame); return(AA_FULL); }

        /* Account for non-affine part: */
        aa_mul_naf_range(x, y, &(z->center), &err);
        if (err >= PlusInfinity)
          { aa_flush(frame); return (AA_FULL); }

        /* Add rounding error term: */
        aa_append_error_term(&(z->nterms), err);

        /* Check for zero: */
        if (AA_ISZERO(z))
          { aa_flush(frame); return(AA_ZERO); }

        #if MIXED
        /* Compute z->range: */
        z->range = ia_meet(
          ia_mul(x->range, y->range),
          aa_implicit_range(z)
        );
        #endif

        #if AA_TRACE_MUL
	  fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
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
      else if ((xr == Zero) || (yr == Zero))
        { return; }
      ROUND_UP;
      xyr = xr * yr;
      if (xyr >= PlusInfinity)
        { *errp = PlusInfinity; return; }
      *errp += xyr;
      if (*errp >= PlusInfinity)
        { *errp = PlusInfinity; return; }
      *zcp = *zcp; /* to supress the "not used" warning */
    #else
      Float xm = aa_max_abs_term((AATermP)(x+1), x->nterms);
      Float ym = aa_max_abs_term((AATermP)(y+1), y->nterms);
      Float ur = Zero;
      Float vr = Zero;

      if ((xm == Zero) || (ym == Zero)) 
        { return; } 

      /* Computes ur = range(u), vr = range(v), where     
            u = 0.5*(x/xm + y/ym),
            v = 0.5*(x/xm - y/ym)
         ignoring center terms of x and y. */
      { AATermP xp = (AATermP)(x+1); AATermCount xn = x->nterms;
        AATermP yp = (AATermP)(y+1); AATermCount yn = y->nterms;
        ROUND_UP;
        while ((xn > 0) || (yn > 0))
          {
            Float xi, yi;
            VarId id;
            if ((xn > 0) && ((yn == 0) || (xp->id <= yp->id)))
              { id = xp->id; xi = Half * xp->coef / xm; }
            else
              { xi = Zero; }
            if ((yn > 0) && ((xn == 0) || (yp->id <= xp->id)))
              { id = yp->id; yi = Half * yp->coef / ym; }
            else
              { yi = Zero; }
            /* The following computations should not overflow: */
            ur += FABS(xi + yi);
            vr += FABS(xi + (-yi));
            if ((xn > 0) && (xp->id == id)) { xp++; xn--; }
            if ((yn > 0) && (yp->id == id)) { yp++; yn--; }
          }
      }
      
      /* Range is now xm*ym*range(u^2 - v^2);   */
      /* approximate it by xm*ym*[-vr^2, +ur^2] */
      { Float ctr, rad;
        ROUND_NEAR;
        ctr = Half * (ur + vr)*(ur - vr) * xm * ym;
        ROUND_UP;
        rad = FMAX(ur * ur * xm * ym - ctr, ctr + vr * vr * xm * ym);
        
        if ((ctr >= PlusInfinity) || (rad >= PlusInfinity))
          { *errp = PlusInfinity; return; }
        *errp += rad;
        if (*errp >= PlusInfinity)
          { *errp = PlusInfinity; return; }
        flt_add(*zcp, ctr, zcp, errp);
        if (*errp >= PlusInfinity)
          { *errp = PlusInfinity; return; }
      }
    #endif
  }

AAP aa_sqr (AAP x)
  {
    #if AA_TRACE_SQR
      fprintf(stderr, "enter aa_sqr:\n");
      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
      fprintf(stderr, "  x = %p nterms = %ld\n", x, x->nterms);
      fprintf(stderr, "  ...\n");
    #endif

    if (AA_ISFULL(x))
      return (AA_FULL);
    else
      { MemP frame = aa_stack_top;
        AAP z = aa_alloc_head();
        Float err = Zero;

        /* Scale noise terms by center value: */
        aa_scale_terms(
          (AATermP) (x + 1), (x->nterms), 
          (x->center), One,
          &(z->nterms),
          &err
        );
        if (err >= PlusInfinity)
          { aa_flush(frame); return (AA_FULL); }
        
        /* Double noise terms and error (should be be exact, barring overflow): */
        { AATermP zp = (AATermP) (z + 1);
          AATermCount zn = (z->nterms);
          ROUND_NEAR;
          while (zn > 0)
            { (zp->coef) *= Two;
              if (FABS(zp->coef) >= PlusInfinity)
                { aa_flush(frame); return (AA_FULL); }
              zp++; zn--;
            }
        }
        err *= Two;
        if (err >= PlusInfinity)
          { aa_flush(frame); return (AA_FULL); }

        /* Square the center value: */
        ROUND_NEAR; flt_mul(x->center, x->center, &(z->center), &err);
        if (FABS(z->center) >= Infinity || err >= Infinity)
          { aa_flush(frame); return(AA_FULL); }

        /* Account for non-affine part: */
        aa_sqr_naf_range(x, &(z->center), &err);
        if (err >= PlusInfinity)
          { aa_flush(frame); return (AA_FULL); }

        /* Add rounding error term: */
        aa_append_error_term(&(z->nterms), err);

        /* Check for zero: */
        if (AA_ISZERO(z))
          { aa_flush(frame); return(AA_ZERO); }

        #if MIXED
        /* Compute z->range: */
        z->range = ia_meet(
          ia_sqr(x->range),
          aa_implicit_range(z)
        );
        #endif

        #if AA_TRACE_SQR
	  fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
	  fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
	  fprintf(stderr, "exit aa_sqr.\n");
        #endif

        return (z);
      }
  }

void aa_sqr_naf_range(AAP x, FloatP zcp, FloatP errp)
  {
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
        Float alpha, gamma, delta;
        AAP z;
        if (xr.lo <= Zero && xr.hi >= Zero)
          { return(AA_FULL); }
        aa_inv_approx(xr, x->center, &alpha, &gamma, &delta);
        z = aa_affine(x, alpha, One, gamma, delta);
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
    FloatP gammap,
    FloatP deltap
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
    *gammap = (dhi + dlo)/Two;
    ROUND_UP;
    *deltap = FMAX(dhi - (*gammap), (*gammap) - dlo);
    if (negative) *gammap = -(*gammap);

    #if AA_TRACE_INV
      fprintf(stderr,
	"    alpha = %15.8e  gamma = %15.8e  delta = %15.8e\n",
	*alphap, *gammap, *deltap
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
        Float alpha, gamma, delta;
        AAP z;
        if (xr.hi < Zero)
          { error("aa_sqrt: argument is negative"); }
        if (xr.hi == Zero)
          { return (AA_ZERO); }
        aa_sqrt_approx(xr, &alpha, &gamma,&delta);
        z = aa_affine(x, alpha, One, gamma, delta);
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
    FloatP gammap,
    FloatP deltap
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
    *gammap = (dm + dlo)/Two;
    ROUND_UP;
    *deltap = FMAX(dm - (*gammap), (*gammap) - dlo);

    #if AA_TRACE_SQRT
      fprintf(stderr,
	"    alpha = %15.8e  gamma  = %15.8e  delta = %15.8e\n",
	*alphap, *gammap, *deltap
      );
      fprintf(stderr, "  exit aa_sqrt_approx.\n");
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
        Float zeta, gamma, delta;
        AAP z;
        if (xr.hi < Zero)
          { error("aa_alt_sqrt: argument is negative"); }
        if (xr.hi == Zero)
          { return (AA_ZERO); }
        aa_alt_sqrt_approx(xr, &zeta, &gamma, &delta);
        z = aa_affine(x, One, zeta, gamma, delta);
        #if MIXED
	  /* Fix z->range: */
	  z->range = ia_meet(ia_sqrt(x->range), z->range);
        #endif
        return(z);
      }
  }  

void aa_alt_sqrt_approx(
    Interval xr,
    FloatP zetap,
    FloatP gammap,
    FloatP deltap
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
    (*gammap) = (dhi + dlo)/Two;
    ROUND_UP;
    (*deltap) = FMAX(dhi - (*gammap), (*gammap) - dlo);
    #if AA_TRACE_SQRT
      fprintf(stderr,
	"    zeta = %15.8e  gamma  = %15.8e  delta = %15.8e\n",
	*zetap, *gammap, *deltap
      );
      fprintf(stderr, "  exit aa_alt_sqrt_approx.\n");
    #endif
  }

AAP aa_div (AAP x, AAP y)
  {
    #if AA_TRACE_DIV
      fprintf(stderr, "enter aa_div:\n");
      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
      fprintf(stderr, "  x = %p nterms = %ld\n", x, x->nterms);
      fprintf(stderr, "  y = %p nterms = %ld\n", y, y->nterms);
      fprintf(stderr, "  ...\n");
    #endif

    if (AA_ISFULL(x) || AA_ISFULL(y))
      { return (AA_FULL); }
    else if (AA_ISZERO(x))
      { return (AA_ZERO); }
    else if ((x->nterms == 0) && (FABS(x->center) >= One))
      { return (aa_affine(aa_inv(y), x->center, One,  Zero, Zero)); }
    else if (y->nterms == 0)
      { return (aa_affine(x, One, y->center,  Zero, Zero)); }
    else
      { Interval yr = aa_range(y);
        if ((yr.lo <= Zero) && (yr.hi >= Zero)) 
          { return(AA_FULL); }
        else
          { MemP frame = aa_stack_top;
            AAP z = aa_alloc_head();
            Float delta = Zero;
            Float err = Zero;
            Float alpha, lambda, sigma, y0;
            
            y0 = y->center;

            /* Compute approximate quotient $\alpha$ of central values: */
            ROUND_NEAR; flt_div(x->center, y->center, &alpha, &delta);
            if (FABS(alpha) >= PlusInfinity || delta >= Infinity)
              { aa_flush(frame); return(AA_FULL); }
            (z->center) = alpha;
        
            /* Compute terms of $(x-\alpha y)/y_0$: */
            aa_mix_terms(
              (AATermP) (x + 1), (x->nterms), One,
              (AATermP) (y + 1), (y->nterms), -alpha,
              y0,
              &(z->nterms),
              &delta
            );
            
            /* Compute mean value $\lambda$ and max deviation $\sigma$ of $y0/y$: */
            { Float yqlo, yqhi;
              if (y0 > Zero)
                { ROUND_DOWN; yqlo = y0 / yr.hi;
                  ROUND_UP; yqhi = y0 / yr.lo;
                }
              else
                { ROUND_DOWN; yqlo = y0 / yr.lo;
                  ROUND_UP; yqhi = y0 / yr.hi;
                }
              if (yqhi >= PlusInfinity)
                { aa_flush(frame); return(AA_FULL); }
              ROUND_NEAR; lambda = yqlo * Half + yqhi * Half;
              ROUND_UP; sigma = FMAX(yqhi - lambda, lambda - yqlo);
            }
            
            /* Multipy $z_i$ ($i=1..n$) by $\lambda \pm \sigma$: */
            aa_fuzzy_rescale_terms(
              (AATermP) (z + 1), &(z->nterms), 
              lambda, sigma, &err
            );

            /* Multiply $\pm\delta$ by $\lambda\pm\sigma$, add to $err$: */
            ROUND_UP; err = err + delta * (lambda + sigma);

            if (err >= PlusInfinity)
              { aa_flush(frame); return(AA_FULL); }

            /* Add rounding error term: */
            aa_append_error_term(&(z->nterms), err);

            /* Check for zero: */
            if (AA_ISZERO(z))
              { aa_flush(frame); return(AA_ZERO); }

            #if MIXED
	      /* Compute z->range: */
	      z->range = ia_meet(
		ia_div(x->range, y->range),
		aa_implicit_range(z)
	      );
            #endif

            #if AA_TRACE_DIV
	      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
	      fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
	      fprintf(stderr, "exit aa_div.\n");
            #endif

            return (z);
          }
      }
  }

AAP aa_abs(AAP x)
  {
    #if AA_TRACE_ABS
      fprintf(stderr, "enter aa_abs:\n");
      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
      fprintf(stderr, "  x = %p nterms = %ld\n", x, x->nterms);
      fprintf(stderr, "  ...\n");
    #endif

    if (AA_ISFULL(x))
      { return (AA_FULL); }
    else if (AA_ISZERO(x))
      { return (AA_ZERO); }
    else
      { Interval xr = aa_range(x);
        if (xr.lo >= Zero)
          { return(x); }
        else if (xr.hi <= Zero)
          { return(aa_neg(x)); }
        else
          { AAP z = aa_alloc_head();
            Float alpha, gamma, delta;
            
            aa_abs_approx(xr, &alpha, &gamma, &delta);
            z = aa_affine(x, alpha, One, gamma, delta);
            
            #if MIXED
	      /* Compute z->range: */
	      z->range = ia_meet(
		ia_abs(x->range),
		aa_implicit_range(z)
	      );
            #endif

            #if AA_TRACE_ABS
	      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
	      fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
	      fprintf(stderr, "exit aa_abs.\n");
            #endif

            return (z);
          }
      }
  }
    
void aa_abs_approx(
    Interval xr,
    FloatP alphap,
    FloatP gammap,
    FloatP deltap
  )
  /* 
    Computes Chebyshev approximation $\alpha x + \gamma$ to $\abs{x}$,
    in the zero-straddling interval $xr$. */
  { 
    Float hd, alpha, gamma, delta;

    #if AA_TRACE_ABS
      fprintf(stderr, "  enter aa_abs_approx:\n");
      fprintf(stderr, "    xr.lo = %15.8e\n", xr.lo);
      fprintf(stderr, "    xr.hi = %15.8e\n", xr.hi);
    #endif

    ROUND_NEAR; hd = xr.hi/Two - xr.lo/Two; /* Beware of overlow */
    if (hd == Zero)
      { alpha = Zero; }
    else
      { alpha = (xr.hi/Two + xr.lo/Two)/hd; }
    gamma = (One - alpha)*xr.hi; /* Shouldn't overflow */

    /* Estimate approximation error $\delta$: */
    ROUND_UP;
    delta = FMAX(
      gamma, 
      FMAX(
	(-xr.lo) + alpha*(-xr.lo) - gamma,
	xr.hi + alpha*(-xr.hi) - gamma
      )
    );
    
    *alphap = alpha;
    *gammap = gamma;
    *deltap = delta;

    #if AA_TRACE_ABS
      fprintf(stderr,
	"    alpha = %15.8e  gamma  = %15.8e  delta = %15.8e\n",
	*alphap, *gammap, *deltap
      );
      fprintf(stderr, "  exit aa_abs_approx.\n");
    #endif
  }
  
AAP aa_max(AAP x, AAP y)
  { 
    if ((AA_ISFULL(x)) || (AA_ISFULL(y)))
      { return (AA_FULL); }
    else
      { Interval xr = aa_range(x);
        Interval yr = aa_range(y);
        if (xr.hi <= yr.lo)
          { return(y); }
        else if (yr.hi <= xr.lo)
          { return(x); }
        else
          { return(
              aa_add(
                aa_affine_2(x, Half, y, Half, One, Zero, Zero),
                aa_abs(aa_affine_2(x, Half, y, -Half, One, Zero, Zero))
              )
            );
          }
      }
  }

AAP aa_min(AAP x, AAP y)
  { 
    if ((AA_ISFULL(x)) || (AA_ISFULL(y)))
      { return (AA_FULL); }
    else
      { Interval xr = aa_range(x);
        Interval yr = aa_range(y);
        if (xr.hi <= yr.lo)
          { return(x); }
        else if (yr.hi <= xr.lo)
          { return(y); }
        else
          { return(
              aa_sub(
                aa_affine_2(x, Half, y, Half, One, Zero, Zero),
                aa_abs(aa_affine_2(x, Half, y, -Half, One, Zero, Zero))
              )
            );
          }
      }
  }

void aa_fuzzy_rescale_terms(
    AATermP zp,
    AATermCount *znp, 
    Float lambda, Float sigma, 
    FloatP errp
  )
  /*
    Multiplies the $*znp$ terms starting at $*zp$ by $\lambda \pm \delta$.
    Adds to $*errp$ all rounding errors, plus the uncertainty due to $\pm \delta$.

    Assumes the given term string is at the top of the stack.  May
    compress the string, pull back the top-of-stack pointer, and
    update $*znp$, if any coefficients become zero due to underflow.

    In case of overflow, sets $*errp$ to PlusInfinity and returns. */
  {
    AATermP xp = zp;
    AATermCount xn = (*znp);
    AATermCount zn = 0;
    Float xc, zc;
    
    while (xn > 0)
      { xc = (xp->coef);
        ROUND_NEAR; flt_mul(xc, lambda, &zc, errp);
        ROUND_UP; (*errp) = (*errp) + sigma * FABS(xc);
        if (zc != Zero) 
          { (zp->coef) = zc; (zp->id) = (xp->id); zp++; zn++; }
        xp++; xn--;
      };
    (*znp) = zn;
    
    /* Pull back top-of-stack pointer. */
    assert (
      (((MemP)zp) == aa_stack_top), 
      "aa_fuzzy_rescale_terms: not on top of stack!"
    );
    aa_stack_top = ((MemP) xp);
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

AAP aa_int_const(int i)
  {
    AAP z = aa_alloc_head();
    Float ilo, ihi;

    #if AA_TRACE_CONST
      fprintf(stderr, "enter aa_int_const:\n");
      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
      fprintf(stderr, "  z = %p\n", z);
      fprintf(stderr, "  ...\n");
    #endif

    ROUND_NEAR;
    z->center = (Float) i;
    ROUND_DOWN;
    ilo = (Float) i;
    ROUND_UP; 
    ihi = (Float) i;
    if (ilo == ihi)
      { z->nterms = 0; }
    else
      {
	AATermP zp;
	zp = aa_alloc_term();
        ROUND_UP;
	zp->coef = FMAX(ihi - z->center, z->center - ilo);
	zp->id = aa_next_id;
	aa_next_id++;
	z->nterms = 1;
      }

    #if MIXED
      z->range.lo = ilo;
      z->range.hi = ihi;
    #endif

    #if AA_TRACE_CONST
      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
      fprintf(stderr, "  z = %p nterms = %ld\n", z, z->nterms);
      fprintf(stderr, "exit aa_int_const.\n");
    #endif

    return (z);
  }

AAP aa_from_interval(Interval x)
  {
    if (ia_is_full(&x))
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

AAP aa_join (AAP x, AAP y)
  { 
    #if AA_TRACE_JOIN
      fprintf(stderr, "enter aa_join_2\n");
      fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
      fprintf(stderr, "  x = %p  nterms = %ld\n", x, x->nterms);
      fprintf(stderr, "  y = %p  nterms = %ld\n", y, y->nterms);
      fprintf(stderr, "  ...\n");
    #endif

    if (AA_ISFULL(x) || AA_ISFULL(y))
      { return(AA_FULL); }
    else
      { MemP frame = aa_stack_top;
        Float err = Zero;
        AAP z = aa_alloc_head();

        /* Combine center values: */
        { Float zt, zd;
          
          Float xc = x->center;
          Float yc = y->center;
          
          ROUND_NEAR; 
          flt_mix(xc, Half, yc, Half, One, &zt, &err);
          if (xc < yc) { Float tc = xc; xc = yc; yc = tc; }
          ROUND_UP;
          zd = Half*xc + (-Half)*yc; /* Should be non-negative! */
          err += zd;

          if ((FABS(zt) >= PlusInfinity) || (err >= PlusInfinity))
            { aa_flush(frame); return(AA_FULL); }
          
          (z->center) = zt;
        }
          
        /* Merge noise terms: */
        aa_mix_terms (
          (AATermP) (x + 1), (x->nterms), Half,
          (AATermP) (y + 1), (y->nterms), Half,
          One,
          &(z->nterms),
          &err
        );

        /* Compute terms of $(y-x)/2$, add them to $err$: */
        { AATermP xp = (AATermP) (x + 1);
          AATermCount xn = x->nterms;
	  AATermP yp = (AATermP) (y + 1);
          AATermCount yn = y->nterms;
	  Float xt, yt, dt;

	  while (xn > 0 || yn > 0)
	    { /* Pick next term:*/
	      if (yn == 0) 
		{ xt = xp->coef; yt = Zero; xp++; xn--; }
	      else if (xn == 0)
		{ xt = Zero; yt = yp->coef; yp++; yn--; }
	      else if (xp->id < yp->id)
		{ xt = xp->coef; yt = Zero; xp++; xn--; }
	      else if (xp->id > yp->id)
		{ xt = Zero; yt = yp->coef; yp++; yn--; }
	      else
		{ xt = xp->coef; yt = yp->coef; 
                  xp++; xn--; yp++; yn--;
                }

	      /* Mix the two terms: */
              if (xt < yt) { Float tt = xt; xt = yt; yt = tt; }
              ROUND_UP; 
              dt = Half*xt + (-Half)*yt; /* Should be non-negative! */
              err += dt;
	    }
	}
        
        if (err >= PlusInfinity) 
          { aa_flush(frame); return(AA_FULL); }
        
        /* Add error term: */
        aa_append_error_term(&(z->nterms), err);

        /* Check for zero: */
        if (AA_ISZERO(z))
          { aa_flush(frame); return(AA_ZERO); }

        #if MIXED
	  /* Compute z->range: */
	  z->range = ia_meet(
	    ia_join(x->range, y->range),
	    aa_implicit_range(z)
	  );
        #endif

        #if AA_TRACE_JOIN
	  fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
	  fprintf(stderr, "  frame = %p\n", frame);
	  fprintf(stderr, "  z =     %p nterms = %ld\n", z, z->nterms);
	  fprintf(stderr, "exit aa_join.\n");
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
        Float err = Zero;

        #if AA_TRACE_FIX
	  fprintf(stderr, "enter aa_fix_eps:\n");
	  fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
	  fprintf(stderr, "  x =   %p nterms = %ld\n", x, x->nterms);
	  fprintf(stderr, "  eps = %p neps =   %ld\n", eps, neps);
	  fprintf(stderr, "  ...\n");
        #endif

        z->center = x->center;

        /* Merge noise terms, and fix matches: */
        {
          AATermP zp;
          Float tmp;
          AATermP xp = (AATermP) (x + 1);
          AATermP ep = (AATermP) eps;
          AATermCount xn = x->nterms;
          AATermCount en = neps;
          AATermCount zn = 0;

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
                  ROUND_NEAR;
                  flt_add(z->center, tmp, &(z->center), &err);
                  if (FABS(z->center) >= Infinity || err >= Infinity)
                    { aa_flush(frame); return(AA_FULL); }
                  xp++; xn--;
                  ep++; en--;
                }
            }
          z->nterms = zn;
        }

        /* Add rounding error term: */
        aa_append_error_term(&(z->nterms), err);

        #if MIXED
	  /* Compute z->range: */
	  z->range = aa_implicit_range(z);
        #endif

        #if AA_TRACE_FIX
	  fprintf(stderr, "  stack_top = %p\n", aa_stack_top);
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

