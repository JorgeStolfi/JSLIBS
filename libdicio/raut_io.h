#ifndef raut_io_H
#define raut_io_H

/* Last edited on 2009-10-29 20:57:57 by stolfi */

#include <stdio.h>
#include <raut.h>

/* AUTOMATON IO */
 
void raut_write(FILE *wr, raut_t *A);
  /* Writes a description of {A} to {wr}. The description is human-readable
    but designed for reading by {raut_read}. */
    
raut_t *raut_read(FILE *rd);
  /* Creates an autoamton from its description in {rd}. Assumes
    the same format used by {raut_write}. */

/* MISCELLANEOUS */

void raut_state_debug(FILE *wr, char* pref, raut_t *A, raut_state_t v, char *suff);
  /* Writes to {vr} a readable representation of the accept bit {v.ac}
     and the dag data fields of state {v}. */

#endif
