/* basic FOI definitions */

#ifndef FOI_H
#define FOI_H

#include "foifloat.h"
#include "foimisc.h"
#include "interval.h"
#include <stdio.h>
#include <math.h>

/*** TYPES ***/

typedef
  unsigned long VarId;

typedef
  struct {
      VarId  id;
      Float  coef;
    } Term;

typedef
  Term *TermP;

typedef
  unsigned long TermCount;

#define MIXED 0

typedef
  struct {
      #if MIXED
        Interval range;
      #endif
      TermCount nterms;
      Float center;
    } FOIHead;

typedef
  FOIHead *FOIP;

/*** CONSTANTS ***/

#define FOI_FULL ((FOIP) -1)
#define FOI_ISFULL(x) (x == FOI_FULL)

/*** INITIALIZATION ***/

void foi_init (void);
  /* Initializes the FOI module (stack, $foi_next_id$, etc.) */

/*** ARITHMETIC ***/

/* These routines allocate space for the result from the stack */
/* They may return a poiner to one of their arguments, or a constant */

FOIP foi_zero(void);

FOIP foi_add (FOIP x, FOIP y);

FOIP foi_sub(FOIP x, FOIP y);

FOIP foi_neg (FOIP x);

FOIP foi_mul (FOIP x, FOIP y);

FOIP foi_sqr (FOIP x);   /* same as foi_mul(x,x), only better */

FOIP foi_inv (FOIP x);   /* 1/x */

FOIP foi_sqrt(FOIP x);

FOIP foi_affine(FOIP x, Float alpha, Float beta, Float gamma);
  /* Returns $\alpha x + \beta$.
     If $\alpha$ and $\beta$ are only approximate,
     $\gamma$ should be an upper bound on the absolute error
     between $\alpha x + \beta$ and the intended function,
     over the range of $x$.  */

/*** CONVERSION ***/

Interval foi_range (FOIP x);
  /* Returns the range of $x$, considering $x->range$ (if any). */

FOIP foi_const(Float c, Float err);
  /* Makes a new FOI with center value $c$ and
     an independent error term $err$. */

FOIP foi_from_interval (Interval x);
  /* Makes a new FOI from the given interval,
     using a brand-new noise variable. */

/*** MISCELLANEOUS TOOLS ***/

void foi_print (FILE *f, FOIP x);
  /* Prints out $x$ on file $f$. */

void foi_move (FOIP source, FOIP destination);
  /* Copies the FOI $*source$ to $*destination$. */
  /* The client should make sure that there is enough space there. */

Float foi_sum_abs_terms (TermP xp, TermCount n);
  /* Returns the sum of the absolute values of the deviations 
     $xp[0..n-1].coef$, rounded up. */

Interval foi_implicit_range (FOIP x);
  /* Returns the range of $x$, considering only the affine form. */

FOIP foi_throw (int nterms);
  /* Returns a random FOI with up to $nterms$ terms, suitable for testing. 
     The $VarId$s of the terms will be 0 through $nterms-1$, but
     some of the terms may be zero (i.e. missing).
     The client must have called $srandom(<seed>)$. */

/*** STACK ALLOCATION ***/

MemP foi_top (void);
  /* Returns the current top-of-stack pointer */

FOIP foi_alloc_head(void);
  /* Reserves space at the top of the FOI stack */
  /* for inserting another FOIHead; returns pointer to it. */

void foi_pop_head(void);
  /* Assumes a FOIHead is at the top of the stack; removes it. */

TermP foi_alloc_term(void);
  /* Reserves space at the top of the stack */
  /* for inserting another Term; returns pointer to it. */

TermP foi_push_term(Float coef, VarId id);
  /* Pushes a new Term on top of the stack, with given fields. */

void foi_pop_term(void);
  /* Assumes a Term is at the top of the FOI stack; removes it. */

void foi_flush (MemP frame);
  /* Ends the logical scope started at the given $frame$  */
  /* (and all inner blocks), freeing their storage. */

FOIP foi_return (MemP frame, FOIP result);
  /* If the $result$ is not on the stack, does $foi_flush(frame)$. */
  /* Else saves $*result$, does a $foi_flush(frame)$, then pushes */
  /* back the saved result onto the FOI stack. */
  /* In either case, returns the result's new address. */

/*** HEAP ALLOCATION ***/

FOIP foi_heap_alloc (TermCount n);
  /* Allocates new space for a FOI with $n$ terms */
  /* (from the $faralloc$ pool). */
  /* Bombs if runs out of memory. */

FOIP foi_heap_realloc (FOIP x, TermCount n);
  /* Re-allocates space for the FOI $x$ (which must have been */
  /* allocated with $foi_heap_alloc$), for $n$ terms, and copies $x$'s */
  /* current contents there, and frees the old storage. Returns */
  /* a pointer to the new area (which may be the same as $x$). */
  /* Bombs if runs out of memory. */

void foi_heap_free (FOIP x);
  /* Releases space used by FOI $x$ (which must have been */
  /* allocated by $foi_heap_alloc$. */

#endif
