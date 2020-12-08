/* See sexp.h. */
/* Last edited on 2005-06-27 23:41:23 by stolfi */

#include <sexp.h>

#include <bool.h>
#include <affirm.h>

#include <stdlib.h>

bool_t sexp_is_empty_list(sexp_t *s)
  { return (s == NULL); }

bool_t sexp_is_proper_list(sexp_t *s)
  { return (s != NULL) && (! sexp_is_symbol(s)) && (! sexp_is_number(s)); }
  
bool_t sexp_is_atom(sexp_t *s)
  { return sexp_is_symbol(s) || sexp_is_number(s); }

bool_t sexp_is_symbol(sexp_t *s)
  { return 
     (s != NULL) && 
     ( (s->tail == SEXP_SYM_L) || 
       (s->tail == SEXP_SYM_S)
     );
  }

bool_t sexp_is_number(sexp_t *s)
  { return 
      (s != NULL) && 
      ( (s->tail == SEXP_NUM_S) || 
        (s->tail == SEXP_NUM_L) || 
        (s->tail == SEXP_NUM_E) || 
        (s->tail == SEXP_NUM_Q) || 
        (s->tail == SEXP_NUM_U) || 
        (s->tail == SEXP_NUM_D)
      );
  }

sexp_t *sexp_head(sexp_t *s)
  { demand(s != NULL, "argument is empty");
    demand(! sexp_is_atom(s), "argument is atom");
    return (sexp_t *)s->head;
  }

sexp_t *sexp_tail(sexp_t *s)
  { demand(s != NULL, "argument is empty");
    demand(! sexp_is_atom(s), "argument is atom");
    return (sexp_t *)s->tail;
  }
