/* Basic Affine Arithmetic definitions and operations */
/* Last edited on 2023-02-18 02:30:29 by stolfi */

#ifndef aa_H
#define aa_H

#include <flt.h>
#include <ia.h>
#include <stdio.h>
#include <math.h>

/* AFFINE FORMS */

typedef unsigned long VarId; 
  /* The ID number of a /noise symbol/, a symbolic variable that
    respresents some error source. The value of a noise symbol is
    usually unknown but, by definition, lies between -1 and +1. */

typedef struct AAHead *AAP;  
  /* An {AAP} is a pointer to an affine form. The form itself consists
    of an {AAHead} and zero or more {AATerm}s, all stored in
    consecutive memory positions. The terms must be sorted by
    increasing noise symbol {id}. */

#define MIXED 1
  /* Change this to 1 to get the mixed IA/AA computation model.
    In this model, the head includes an ordinary interval. 
    All AA functions also perform the corresponding IA operation,
    and the two results (IA and AA) are used synergistically
    in subsequent operations. */

typedef unsigned long AATermCount;
  /* A count of terms in some affine form. */ 

typedef struct AAHead  /* Head of an affine form. */
  {
#if MIXED
    Interval range;      /* Ordinary interval for the quantity. */
#endif
    AATermCount nterms;  /* Number of 1st degree terms. */
    Float center;        /* Center value of affine form. */
  } AAHead;
  /* The head of an affine form contains the number of terms {nterms},
    and the central value {x_0} of the form --- the term of degree 0,
    that does not depend on any noise symbols. */

typedef struct AATerm 
  { VarId  id;   /* ID number of a noise symbol. */
    Float  coef; /* Coefficient of noise symbol. */
  } AATerm;
  /* A term of degree 1 in some affine form.  It represents the given {coef}
    times the noise variable with number {id}. */

typedef AATerm *AATermP;
  /* A pointer to some 1st-degree term of an affine form. */

/* INITIALIZATION */

void aa_init (void);
  /* Initializes the AA module (the AA stack, {aa_next_id}, etc.) */

/* CONSTANTS AND PREDICATES */

AAP aa_full(void); 
  /* Returns a pointer to the /anything/-form, a special affine form
    that represents the interval {(-oo _ +oo)} of all (finite) real
    numbers. */

int aa_is_full (AAP x); 
  /* A predicate that returns 1 (true) iff {x} represents the 
    interval {(-oo _ +oo)}. */

AAP aa_zero(void); 
  /* Retuns an affine form for the constant 0. */

AAP aa_one(void);
  /* Retuns an affine form for the constant 1. */

int aa_is_zero (AAP x); 
  /* True iff the affine form {x} represents the exact value 0. */

/* CONVERSION FROM NUMBERS AND ORDINARY INTERVALS TO AFFINE FORMS

  Each call to a procedure in this section generally returns an affine
  form that depends on a brand-new noise symbol. However, if the
  result is exactly representable as a single floating-point value,
  the returned form will have zero noise terms. */

AAP aa_const(Float c, Float err);
  /* Makes a new affine form with center value {c} and
    a noise term with coefficient {err}, depending 
    on a brand new noise symbol. */

AAP aa_int_const (int i);
  /* Makes an affine form for the integer {i}. If the integer
    is too big, the form may not be exact. */

AAP aa_from_interval (Interval x);
  /* Makes a new affine form from the given interval,
    using a brand-new noise variable. */

/* CONVERSION OF AFFINE FORMS TO ORDINARY INTERVALS */

Interval aa_range (AAP x);
  /* Returns an interval that contains all real numbers described by
    the affine form {x}.  If {MIXED} is 1, this procedure also takes into
    account the IA interval {x->range} as well. */

/* ARITHMETIC OPERATIONS ON AFFINE FORMS

  Each procedure in this section takes as arguments the addresses of
  one or more affine forms, computes a new affine form, and
  returns the adddress of this form. 
  
  Generally, a procedure with {n} arguments in this section computes
  an affine form that describes the result of applying some arithmetic
  operation or basic function to all {n}-tuples of real numbers
  descibed by the given forms. Thus, for example, {aa_add(x,y)}
  computes an affine form for {av+yv}, for all pairs {(xv,yv)}
  described by the forms {x,y} together. The resulting form tries to
  take into account any linear dependency between the two arguments
  (implied by the sharing of noise symbols) and tries to preserve that
  information in the result.  
  
  If the operation cannot be performed exactly, any rounding or
  truncation errors are combined into a new term, that depends on a
  brand-new noise symbol.
  
  The argument forms may reside anywhere in memory, and are not
  modified in any way. The procedure typically allocates space for the
  result from the private AA stack (see AA STACK ALLOCATION below). In
  some cases, however, the procedure may return a pointer to one of the
  given forms, or to some statically allocated constant form.
  Therefore, clients should never try to modify the contents of an
  affine form that has been returned by any of these procedures. */

AAP aa_add (AAP x, AAP y);
  /* Computes an affine form for {-x}. */
  
AAP aa_neg (AAP x);
  /* Computes an affine form for {-x}. */

AAP aa_sub (AAP x, AAP y);
  /* Computes an affine form for {x-y}. */

AAP aa_scale (AAP x, Float alpha, Float zeta); 
  /* Computes an affine form for {alpha*x + zeta}. */

AAP aa_shift (AAP x, Float gamma);
  /* Computes an affine form for {x + gamma}. */

AAP aa_affine(AAP x, Float alpha, Float zeta, Float gamma, Float delta);
  /* Computes an affine form for general 1st degree operation on {x},
    namely {alpha*x/zeta + gamma ± delta}. The {± delta} term and all
    rounding errors are lumped in a single terms with a new noise
    symbol. */

AAP aa_affine_2
  ( AAP x, Float alpha, 
    AAP y, Float beta, 
    Float zeta, 
    Float gamma, 
    Float delta
  );
  /* Computes an affine form for general 1st degree operation 
    on {x} and {y}, namely {(alpha x + beta y)/zeta + gamma ± delta}.
    The {± delta} term and all rounding errors are lumped in a
    single terms with a new noise symbol. */

AAP aa_mul (AAP x, AAP y);
  /* Computes an affine form for {x*y}. */

AAP aa_inv (AAP x); 
  /* Computes an affine form for {1/x}. */

AAP aa_div (AAP x, AAP y);
  /* Computes an affine form for {x/y}. */

AAP aa_sqr (AAP x); 
  /* Computes an affine form for {x^2}. Usually
    provides a tighter result than {aa_mul(x,x)}. */

AAP aa_sqrt(AAP x);
  /* Computes an affine form for {sqrt(x)}. */

AAP aa_alt_sqrt(AAP x);
  /* An alternative implementation of {aa_sqrt}. */

AAP aa_exp(AAP x);
  /* Computes an affine form for {exp(x)}. */

AAP aa_abs(AAP x);
  /* Computes an affine form for {abs(x)}. If {x} has a definite sign,
    the result is either {x} itself or {aa_neg(x)}. Otherwise the
    procedure uses a new form. */

AAP aa_max(AAP x, AAP y);
AAP aa_min(AAP x, AAP y);
  /* These procedures compute affine forms for {max(x,y)}
    and {min(x,y)}, respectively.  If {x} and {y}
    have disjoint ranges, returns either {x} or {y}.
    Otherwise returns a new form. */

AAP aa_join (AAP x, AAP y);
  /* Returns an affine form that describes an arbitrary
    convex combination of {x} and {y}; that is, a quantity that
    is known to lie between some value of {x} and some value of {y}. */

/* MANIPULATING NOISE TERMS */

AAP aa_fix_eps(AAP x, AATermCount neps, AATerm eps[]);
  /* Returns an affine form obtained from {x} by fixing the value of
    the noise symbol with {id = eps[i].id} to {eps[i].coef}, for 
    every {i} in {0..neps-1}.

    Any noise symbols that appear in {x} but not in {eps} are retained
    in the result.  A new noise symbol is added, if necessary, to
    account for rounding errors.
    
    If {x} is the full affine form, return the full affine form.

    The fields {eps[i].id} must be sorted in increasing order.
    The values {eps[i].coef} need not be in [-1 __ +1]. */

void aa_get_eps_coef(AAP x, AATermCount neps, AATerm eps[]);
  /* Extracts from the affine form {xr} the coefficients of each noise
    variable with {id = eps[i].id} and saves it into {eps[i].coef},
    for every {i} in {0..neps-1}. In particular, if {x} does not
    depend on the noise variable {eps[i].id}, stores zero into
    {eps[i].coef}.
    
    The procedure fails if {x} is the full affine form.

    The fields {eps[i].id} must be sorted in increasing order. */

void aa_collapse_pair(AAP x, AAP y, AAP *xr, AAP *yr);
  /* Returns in {*xr} and {*yr} two affine forms that describe the same 
    quantities described by {x} and {y}, but depend on at most two noise
    symbols (possibly new ones), which are shared between them. */

/* MISCELLANEOUS TOOLS */

void aa_print (FILE *f, AAP x);
  /* Prints out {x} on file {f}. */

Float aa_sum_abs_terms (AATermP xp, AATermCount n);
  /* Returns the sum of the absolute values of the deviations 
    {xp[0..n-1].coef}, rounded up. */

Float aa_max_abs_term (AATermP xp, AATermCount n);
  /* Returns the maximum absolute magnitude of the partial deviations
    {xp[0..n-1].coef}. */
  
Interval aa_implicit_range (AAP x);
  /* An interval that contains all real numbers described by the
    affine form {x}, considering only the AA information proper.
    
    Thus, if {MIXED} is 0, {aa_range(x)} is equivalent to this
    function; if {MIXED} is 1, {aa_range(x)} returns the intersection
    of {aa_implicit_range(x)} and {x->range}. */

AAP aa_throw (int nterms);
  /* Returns a random affine form with up to {nterms} terms, suitable for testing. 
    The {VarId}s of the terms will be 0 through {nterms-1}, but
    some of the terms may be zero (i.e. missing).
    The client must have called {srandom(<seed>)}. */

/* AA STACK ALLOCATION

  In order to avoid the relatively high cost of {malloc()},
  an AA procedure that needs to create new forms (such as {aa_add} and {aa_const})
  will allocate space for their result from a private pool, the 
  /AA stack/. 
  
  As the name implies, the AA stack is (mostly) managed in LIFO fashion. 
  The pool is a contiguous memory area, divided into
  two segments, `used' and `free'.  The starting address of the free 
  segment is the `top' of the stack.  Storage is allocated by incrementing
  the top-of-stack pointer, and freed by decrementing it.
  
  Clients that need to build an affine form term-by-term may use
  {aa_alloc_head}, {aa_alloc_term}, {aa_push_term}, and
  {aa_append_error_term}.  */

#define aa_STACK_SIZE 100000L

AAP aa_alloc_head(void);
  /* Reserves space at the top of the AA stack
    for inserting another AAHead; returns pointer to it. */

AATermP aa_alloc_term(void);
  /* Reserves space at the top of the stack
    for inserting another AATerm; returns pointer to it. */

AATermP aa_push_term(Float coef, VarId id);
  /* Pushes a new AATerm on top of the stack, with given fields. */

void aa_append_error_term (AATermCount *znp, Float err);
  /* If {err} is not zero, pushes an extra AATerm on the stack,
    with {err} as the coefficient and with a newly allocated noise
    symbol {id}, and increments {*znp}.  Otherwise does nothing. */

/* AA STACK DEALLOCATION
  
  There is no automatic deallocation for the AA stack. In particular,
  an AA operation never deletes its argument forms; it just pushes its
  result onto the stack. 
  
  Complicated or iterative AA computations often create many
  affine forms that become useless after a short while.
  In order reclaim the stack space used by such forms, clients
  must save the top-of-stack pointer with {aa_top} (see below) before
  the computation, and explicitly reset it with {aa_flush} at the end. 
  
  Thus, the typical paradigm for AA computations is
  
    { MemP tp = aa_top(); // Save the current position of the AA stack.
      ...                 // Code that calls {aa_const}, {aa_add}, etc.
      aa_flush(tp);       // Restore the AA stack to the original level.
    }

  One often needs to flush all temporary affine forms that
  were created during a complicated computation, except 
  one specific form that is to be retained on the stack
  for further use.  The procedure {aa_return} is intended 
  for such situations.  In particular,  typical paradigm 
  for a client procedure that returns one affine form is

    AAP myformula(AAP x, AAP y, ...)
      { MemP tp = aa_top();   // Save the current AA stack position.
        AAP r;                
        ...                   // Code that calls {aa_const}, {aa_add}, etc.
        r = ...;              // Compute the final result;
        ...                   // More AA code.
        r = aa_return(tp,r);  // Flush the AA stack keeping the form {r}.
        return r;
      }

  Clients that build their forms term-by-term may find a use for
  {aa_pop_term} and {aa_pop_head}. 
  
  */

typedef void * MemP;  /* Pointer to memory area */
typedef size_t MemSize;  /* Size in bytes of memory area */

MemP aa_top (void);
  /* Returns the current top-of-stack pointer */

void aa_flush (MemP frame);
  /* Ends the logical scope started at the given {frame}
    (and all inner blocks), freeing any AA stack storage
    allocated within their scope. */

AAP aa_return (MemP frame, AAP result);
  /* If the {result} is not on the stack, does {aa_flush(frame)}.
    Else saves {*result}, does a {aa_flush(frame)}, then pushes
    back the saved result onto the AA stack.
    In either case, returns the result's final address. */

void aa_return_n (MemP frame, int n, AAP result[]);
  /* For each {i} in {0..n-1}, checks whether the AA form that starts
    at {result[i]} lies in the AA stack above the {frame} address.
    Then copies all the AA forms that satisfy that condition to
    consecutive positions starting at the {frame} address, and updates
    the corresponding addresses {result[i]}. AA forms that do not
    satisfy the condition are left alone. Then sets the top-of-stack
    pointer to the end of those copies. */

void aa_pop_head(void);
  /* Assumes an AAHead is at the top of the stack; removes it. */

void aa_pop_term(void);
  /* Assumes a AATerm is at the top of the AA stack; removes it. */

/* OFF-STACK STORAGE OF AFFINE FORMS

  Affine forms whose useful lives are not compatible with the LIFO
  policy of the AA stack should be copied elsewhere.
  The function {aa_move} may be used for that purpose. */

void aa_move (AAP source, AAP destination);
  /* Copies the affine form {*source} to {*destination}.
    It is up to the client to ensure that there is enough space there. */

/* HEAP ALLOCATION 

  Clients who need to store affine forms in the {malloc()}
  heap may find the following procedures useful: */

AAP aa_heap_alloc (AATermCount n);
  /* Allocates new space for a new affine form with {n} terms
    from the {malloc()} pool.  Bombs out if it runs out of memory. */

AAP aa_heap_realloc (AAP x, AATermCount n);
  /* Re-allocates space for the affine form {x} (which must have been
    allocated with {aa_heap_alloc}), for {n} terms, and copies {x}'s
    current contents there, and frees the old storage. Returns 
    a pointer to the new area (which may be the same as {x}).
    Bombs out if it runs out of memory. */

void aa_heap_free (AAP x);
  /* Releases space used by the affine form {x} (which must have been
    allocated by {aa_heap_alloc}. */

#endif
