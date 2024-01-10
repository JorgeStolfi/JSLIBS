/* Basic Affine Arithmetic definitions and operations */

#ifndef AA_H
#define AA_H

#include <js.h>
#include <flt.h>
#include <ia.h>
#include <stdio.h>
#include <math.h>

/*** TYPES ***/

typedef
  unsigned long VarId;

typedef
  struct {
      VarId  id;
      Float  coef;
    } AATerm;

typedef
  AATerm *AATermP;

typedef
  unsigned long AATermCount;

#define MIXED 0

typedef
  struct {
      #if MIXED
        Interval range;
      #endif
      AATermCount nterms;
      Float center;
    } AAHead;

typedef
  AAHead *AAP;

/*** CONSTANTS ***/

#define AA_FULL ((AAP) -1)
#define AA_ISFULL(x) (x == AA_FULL)

/*** INITIALIZATION ***/

void aa_init (void);
  /* Initializes the AA module (stack, $aa_next_id$, etc.) */

/*** ARITHMETIC ***/

/* These routines allocate space for the result from the stack */
/* They may return a poiner to one of their arguments, or a constant */

AAP aa_zero(void);

AAP aa_add (AAP x, AAP y);

AAP aa_sub(AAP x, AAP y);

AAP aa_neg (AAP x);

AAP aa_mul (AAP x, AAP y);

AAP aa_sqr (AAP x);   /* same as aa_mul(x,x), only better */

AAP aa_inv (AAP x);   /* 1/x */

AAP aa_sqrt(AAP x);
AAP aa_alt_sqrt(AAP x);

AAP aa_affine(AAP x, Float alpha, Float beta, Float gamma); 
  /* 
    Returns $\alpha x + \beta \pm \gamma$.
    
    The $\pm \gamma$ term and all rounding errors are lumped in
    a new noise symbol. */

AAP aa_div_affine(AAP x, Float zeta, Float beta, Float gamma);
  /* 
    Like aa_affine, but computes $x / \zeta + \beta$. */

AAP aa_fix_eps(AAP x, AATermCount neps, AATerm eps[]);
  /*
    Returns an affine form obtained from $x$ by fixing the value
    of all noise symbols with $id = eps[i].id$ ($i=0,..neps-1$) to
    $eps[i].coef$.

    Any noise symbols that appear in $x$ but not in $eps$ are retained
    in the result.  A new noise symbol is added, if necessary, to
    account for rounding errors.

    The fields $eps[i].id$ must be sorted in increasing order.
    The values $eps[i].coef$ need not be in [-1 __ +1]. */

void aa_collapse_pair(AAP x, AAP y, AAP *xr, AAP *yr);
  /*
    Returns in $*xr$ and $*yr$ two affine forms that describe the same 
    quantities described by $x$ and $y$, but share at most two noise
    symbols (possibly new ones). */

/*** CONVERSION ***/

Interval aa_range (AAP x);
  /* 
    Returns the range of $x$, considering $x->range$ (if any). */

AAP aa_const(Float c, Float err);
  /* 
    Makes a new affine form with center value $c$ and
    an independent error term $err$. */

AAP aa_from_interval (Interval x);
  /*
    Makes a new affine form from the given interval,
    using a brand-new noise variable. */

/*** MISCELLANEOUS TOOLS ***/

void aa_print (FILE *f, AAP x);
  /*
    Prints out $x$ on file $f$. */

void aa_move (AAP source, AAP destination);
  /*
    Copies the affine form $*source$ to $*destination$.
    It is up to the client to ensure that there is enough space there. */

Float aa_sum_abs_terms (AATermP xp, AATermCount n);
  /* 
    Returns the sum of the absolute values of the deviations 
    $xp[0..n-1].coef$, rounded up. */

Interval aa_implicit_range (AAP x);
  /* 
    Returns the range of $x$, considering only the affine form. */

AAP aa_throw (int nterms);
  /* 
    Returns a random affine form with up to $nterms$ terms, suitable for testing. 
    The $VarId$s of the terms will be 0 through $nterms-1$, but
    some of the terms may be zero (i.e. missing).
    The client must have called $srandom(<seed>)$. */

/*** STACK ALLOCATION ***/

MemP aa_top (void);
  /*
    Returns the current top-of-stack pointer */

AAP aa_alloc_head(void);
  /*
    Reserves space at the top of the AA stack
    for inserting another AAHead; returns pointer to it. */

void aa_pop_head(void);
  /*
    Assumes an AAHead is at the top of the stack; removes it. */

AATermP aa_alloc_term(void);
  /*
    Reserves space at the top of the stack
    for inserting another AATerm; returns pointer to it. */

AATermP aa_push_term(Float coef, VarId id);
  /*
    Pushes a new AATerm on top of the stack, with given fields. */

void aa_pop_term(void);
  /*
    Assumes a AATerm is at the top of the AA stack; removes it. */

void aa_flush (MemP frame);
  /*
    Ends the logical scope started at the given $frame$
    (and all inner blocks), freeing their storage. */

AAP aa_return (MemP frame, AAP result);
  /*
    If the $result$ is not on the stack, does $aa_flush(frame)$.
    Else saves $*result$, does a $aa_flush(frame)$, then pushes
    back the saved result onto the AA stack.
    In either case, returns the result's final address. */

/*** HEAP ALLOCATION ***/

AAP aa_heap_alloc (AATermCount n);
  /*
    Allocates new space for a new affine form with $n$ terms
    (from the $faralloc$ pool).  Bombs out if it runs out of memory. */

AAP aa_heap_realloc (AAP x, AATermCount n);
  /*
    Re-allocates space for the affine form $x$ (which must have been
    allocated with $aa_heap_alloc$), for $n$ terms, and copies $x$'s
    current contents there, and frees the old storage. Returns 
    a pointer to the new area (which may be the same as $x$).
    Bombs out if it runs out of memory. */

void aa_heap_free (AAP x);
  /* 
    Releases space used by the affine form $x$ (which must have been
    allocated by $aa_heap_alloc$. */

#endif
