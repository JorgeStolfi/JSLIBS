#ifndef baltuple_H
#define baltuple_H

/* Enumeration of balanced integer tuples */
/* Last edited on 2008-05-15 15:24:56 by stolfi */

/* These procedures provide tools to enumerate all tuples
  of integers with specified ranges for each component
  and for the sum of all components. */

#include <bool.h>

/* INTEGER ENUMERATION
  
  These procedures enumerate integers in a specified range {vmax..vmin},
  in the  /0-symmetric order/,  defined by the sequence
  {SEQ = (0,-1,+1,-2,+2,..)}. */

int balt_ix_first(int vmin, int vmax);
  /* Returns the first value in 0-symmetric order
    that lies in the range {vmax..vmin}. */

int balt_ix_next(int v, int vmin, int vmax);
  /* Returns the next integer after {v} in 0-symmetric order that
    lies in the range {vmax..vmin}.  If there is no such integer,
    returns some integer outside that range. */
    
/* TUPLE ENUMERATION

  These procedure enumerate {n}-tuples of integers in lexicographic
  order, with specified ranges {vmax..vmin} for each element and
  {smax..smin} for the sum, where each component is compared according
  to 0-symmetric order. */

void balt_first_tuple(int n, int d[], int vmin, int vmax, int smin, int smax);
  /* Sets the tuple {d[0..n-1]} of integers to the first tuple in
    0-symmetric lexicographic order, among those sepcified by
    {vmin..vmax} and {smin..smax}.  Bombs out if there is no such {n}-tuple.  */

bool_t balt_next_tuple(int n, int d[], int vmin, int vmax, int smin, int smax);
  /* Given a tuple {d[0..n-1]} of integers, among the set specified by
    {vmin..vmax} and {smin..smax}, replaces it by the next tuple from
    that set in 0-symmetric lexicographic order. Returns TRUE if there
    is no such tuple, FALSE otherwise. */

#endif
