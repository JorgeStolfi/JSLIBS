/* LISP-like Symbolic Expressions. */
/* Last edited on 2009-02-09 17:42:35 by stolfi */

#ifndef sexp_H
#define sexp_H

#include <bool.h>

#include <stdlib.h>
#include <stdio.h>

/*
  PURPOSE
  
  Symbolic expressions, or /sexps/ for short, are multipurpose
  mathematical structures. Sexps can be used to represent many other
  mathematical structures -- sets, algebraic or logical formulas,
  trees, graphs, etc.. Sexps are especially convenient when describing
  algorithms on such objects, because they provide efficient computer
  representations, internal and external, and an algorithmically
  complete repertoire of basic operations.
  
  Symbolic expressions were introduced by McCarthy in the 1960s,
  as part of his definition of the LISP language.
  
  DEFINITION
  
  A /symbolic expression/ or /sexp/ is either a /list/ or an /atom/. A
  list is either the /empty list/, or is a /proper list/ with a /head/
  (any sexp) and a /tail/ (another list). An atom is either a /symbol/
  (an arbitrary string of Unicode characters, other than U+0000 or
  U+FFFF) or a /number/ (a rational number).
  
  The head and tail of a sexp {E} are undefined if {E} is the empty
  list or an atom. By repeatedly taking the tail of a list {L} one
  must eventually get the empty list; the number of tail steps needed
  is the /length/ of {L}, here denoted {#L}. For any {k} in {0..#L-1},
  /element {k}/ of {L}, here denoted {L[k]}, is {L}'s head if {k=0},
  and element {k-1} of {L}'s tail otherwise. */
  
/*  
  INTERNAL REPRESENTATION */

#define SEXP_SYM_S ((void *)((char *)NULL+1))
#define SEXP_SYM_L ((void *)((char *)NULL+2))
/* Special pseudo-addresses to tag non-numeric atoms. */

#define SEXP_NUM_S ((void *)((char *)NULL+3))
#define SEXP_NUM_L ((void *)((char *)NULL+4))
#define SEXP_NUM_E ((void *)((char *)NULL+5))
#define SEXP_NUM_Q ((void *)((char *)NULL+6))
#define SEXP_NUM_U ((void *)((char *)NULL+7))
#define SEXP_NUM_D ((void *)((char *)NULL+8))
/* Special pseudo-addresses to tag numeric atoms. */

typedef struct sexp_t { void *head; void *tail; } sexp_t;
  /*
    A sexp {E} is represented internally by a pointer value of type
    {sexpr_t *}. The empty list is represented by the NULL address.
    Any other sexp {E} is represented by the address of a {sexp_t}
    record {R} in memory:

      If {E} is a list, then {E.head} and {E.tail} are its 
      head and tail, and both are {sexp_t *} values;

      If {E} is a symbol, then {E.tail} is one of the special
      addresses {SEXP_SYM_x}. If {E} is a number, then {E.tail} is one
      of the special pseudo-addresses {SEXP_NUM_x}. In either case,
      the field {E.head} (which is not a {sexp_t *} value) contains
      the representation of the symbol or number. */ 

#define DEBUG TRUE
/* Define this as TRUE to print various diagnostics. */

/* 
  FUNDAMENTAL OPERATIONS */

bool_t sexp_is_empty_list(sexp_t *s);
  /* TRUE iff {s} is the empty list. */

bool_t sexp_is_proper_list(sexp_t *s);
  /* TRUE iff {s} is a proper (non-empty) list. */
  
bool_t sexp_is_atom(sexp_t *s);
  /* TRUE iff {s} is an atom --- i.e., a symbol or number.
    False if {s} is a list, empty or not. */

bool_t sexp_is_symbol(sexp_t *s);
  /* TRUE iff {s} is a symbol. */

bool_t sexp_is_number(sexp_t *s);
  /* TRUE iff {s} is a number. */

sexp_t *sexp_head(sexp_t *s);
  /* If {s} is a proper (non-empty) list, returns its head
    (some sexp). Otherwise it is a fatal error. */

sexp_t *sexp_tail(sexp_t *s);
  /* If {s} is a proper (non-empty) list, returns its tail
    (a list). Otherwise it is a fatal error. */

/*
  EXTERNAL REPRESENTATION 
  
  The /external representation/, or /print form/, of a sexp {E} is a
  string {S} of Unicode characters that encodes E up to equality:
  
    If {E} is a list, {S} is the concatenation of the print forms of
    its elements, separated by spaces and enclosed in parentheses "("
    and ")". In particular, the print form of the empty list is "()".
    The spaces may be ASCII SP (octal 040) or any ASCII formatting
    characters (octal 011 to 015: HT, NL, VT, NP, CR). Extra spaces
    may be inserted before and after each element.
    
    If {E} is an integer number, then {S} can be its decimal
    representation. If {E} is a number which can be expressed as
    {f*10^p}, where {f} is a decimal fraction and {p} is an integer,
    then {S} can be "{F}@{P}" where {F} and {P} are any signed decimal
    representations of {f} and {p}. If {E} can be expressed as {x/y},
    for some integers {x} and {y}, then {S} can be "{X}/{Y}" where {X}
    and {Y} are any signed decimal representations of {x} and {y}. In
    all cases, each decimal representation must have at least one
    digit, at most one period, and at most one sign. The '+' sign is
    optional and meaninglesss, as well as leading zeros in the integer
    part, and trailing zeros in the fractional part (if any).
  
    If {E} is a symbol, then {S} is either the /naked form/ of {E}, 
    the /latin-quoted form/ of {E}, or the /Unicode-quoted form/
    of {E}:
    
      The naked print form consists of {E} itself. It is valid only if
      {E} consists of one or more printable ASCII characters,
      excluding space, comma, single and double quote, backslash,
      backquote, parentheses, braces, or brackets, and cannot be
      parsed as the print form of a number.
      
      The latin-quoted print form is valid for any symbol. It consists
      of {E} enclosed in ASCII double quotes, where each character {c}
      or {E} may be represented in one of these ways:
      
        (0) the character {c} itself. 
        
        (1) a backslash followed by {c} itself;
        
        (2) a backslash followed by a lower-case letter,
        't' (for ASCII HT = octal 011), 'n' (NL = 012),
        'v' (VT = 013), 'l' (NP = 014), 'r' (CR = 015),
        'e' (ESC = 033);
        
        (3) a backslash and exactly three digits on '0'..'7',
        representing the Unicode position of {c} in octal;
        
        (4) a backslash, the ascii letter 'u' or 'U', and exactly 
        four characters in '0'..'9', 'a'..'f', or 'A'..'F',
        denoting the Unicode position of {c} in hexadecimal.
        
      Choice (0), the character {c} itself, may be used only for ASCII
      space (SP, 040) or printable ISO Latin-1 characters other than
      double quote, backslash, non-breaking space (NBSP, 240), and
      soft hyphen (255). Choice (1) may be used only for double quotes
      and backslash. Choice (3) may be used for any character with
      Unicode position between 0 and 255 decimal (000 and 377 octal).
      Choice (4) may be used for any character (other than U+0000 or
      U+FFFF, which are not valid in Unicode).
      
      The Unicode-quoted print form too is valid for any symbol. It
      consists of {E}, with the following changes: (1) a backslash is
      inserted before every ISO Latin-1 backslash or ISO Latin-1
      closing guillemot; and (2) everything is enclosed in ISO Latin-1
      guillemots "«" and "»".
      
    Note that the (unique) symbol with zero characters has only the
    quoted print forms, consisting of two ASCII double quotes or an
    open-close guillemot pair. It should not be confused with the
    empty list, whose print form is "()".
    
  Note that the same expression may have several print forms. */
  
void sexp_print(FILE *wr, sexp_t *e);
  /* Writes the print form of {e} to file {wr}, using a minimum of
    spaces (ASCII SP) for separators, and no line breaks. */

void sexp_pprint(FILE *wr, sexp_t *e, int ind, int maxw, int minw, int step);
  /* Writes the print form of {e} to file {wr}. If the expression does
    ot fit in a single line, uses spaces and newlines to obtain
    suitable alignment and indantation between conecutive lines.
    
    Assumes that {wr} is at the beginning or just after a newline. The
    parameter {ind} is the indentation to use for the whole formula,
    i.e. the number of spaces to print before it. Then, if {e} is an
    atom, writes its print form. A list is printed in one line if it
    fits, else it is printed one element per line, with each element
    indented {step} extra spaces.
    
    The maximum line width is assumed to be the maximum of {maxw} and
    {ind+minw}. This limit may be exceeded by long atoms. */

#endif
