#ifndef affirm_H
#define affirm_H

/* Variants of {assert} with explicit message argument. */
/* Last edited on 2024-10-23 07:07:21 by stolfi */

/* ERRORS AND ASSERTIONS */

#define fatalerror(msg) \
  programerror((msg), __FILE__, __LINE__, __FUNCTION__)
  /* Prints string {*msg} to {stderr} and stops. */

#define affirm(test, msg) \
  do \
    { if (!(test)) programerror((msg), __FILE__, __LINE__, __FUNCTION__); } \
  while (0)
  /* Like the standard {assert}, but with an error message. If {test} is
    false, aborts after printing {*msg} and the source location.
    
    It is declared as a macro, rather than as a procedure, in order to
    avoid evaluating {msg} when {test} is true. The macro is wrapped
    in a {do {...} while(0)} in order to make it behave like a
    procedure call vis-a-vis semicolons, {else}s, etc.. This is an
    obscure but standard C programming trick. */

#define demand affirm
  /* Same as {affirm}, but intended to check the validity of
    client/user data (rather than paranoia checks). */

#define fail_test(die,msg) \
  do { \
    if (die) { fatalerror((msg)); } \
    return FALSE; \
  } while(0)
  /* If the parameter {die} is TRUE, aborts the program with message
    {msg}; else returns from the current procedure with a FALSE result
    but without any message. 
    
    This macro is useful in functions that check complicated
    conditions and invariants of data structures.
    
    The macro is wrapped in a {do {...} while(0)} in order to make it
    behave like a procedure call vis-a-vis semicolons, {else}s,
    etc..  This is an obfuscating but standard C programming trick. */

void programerror (const char *msg, const char *file, unsigned int line, const char* proc)
  __attribute__ ((noreturn));
  /* Prints {file ":" line ": (" *proc ")" *msg} to {stderr} and stops.
    Meant for use by {affirm}. */

/* NULL POINTER CHECKING */

#define notnull(p, msg) \
  checknotnull((p), (msg), __FILE__, __LINE__, __FUNCTION__)
  /* If {p == NULL}, aborts with {*msg} with the program 
    location; otherwise returns {p} itself. */
  
void *checknotnull(void *p, const char *msg, const char *file, unsigned int line, const char *proc);
  /* If {p == NULL}, prints {file ":" line ": (" *proc ")" *msg} to {stderr} and
    stops; otherwise returns {p} itself.  Meant for {notnull} below. */

/* SANE HEAP ALLOCATION */
   
#define talloc(n, T) \
  ((T*)((n) == 0 ? NULL : notnull(calloc((n), sizeof(T)), "no mem")))
  /* Allocates an array of {n} elements of type {T}, casting the result as a {T*} pointer.
    If {n} is zero, returns {NULL}.  Aborts with error if {n} is positive but the allocation
    returns {NULL} (presumably because of not enough memory).
    Equivalent to {calloc(n,sizeof(T))} except for the type casting and {NULL} checking. */
    
#define retalloc(p, n, T) \
  ((T*)((n) == 0 ? reallocarray((p), 0, sizeof(T)) : notnull(reallocarray((p), (n), sizeof(T)), "no mem")))
  /* Re-allocates the heap area pointed by {p} as an array of {n} elements of type {T}, 
    freeing the previous area (if {p} was not {NULL}), and casting the result as a {T*} pointer.
    If {n} is zero, performs {free(p)} and returns {NULL}.  Aborts with error if {n} is positive
    but the reallocation returns {NULL} (presumably because of not enough memory).
    Equivalent to {reallocarray(p,n,sizeof(T))} except for the type casting and {NULL} checking. */

#endif
