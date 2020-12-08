/* See {spmat.h} */

#define spmat_linalg_C_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2009-08-31 21:48:46 by stolfi */

#include <spmat.h>
#include <stdlib.h>
#include <stdint.h>
#include <affirm.h>

void *spmat_alloc(spmat_count_t ents, size_t esz)
  { demand(ents <= spmat_MAX_ENTS, "too many entries");
    void *e = (ents == 0 ? NULL : malloc(ents*esz));
    affirm((ents == 0) || (e != NULL), "out of mem");
    return e;
  }

void spmat_expand(void **eP, spmat_count_t *entsP, spmat_pos_t index, size_t esz)
  { if (index >= (*entsP))
      { demand(index <= spmat_MAX_POS, "index too large");
        if ((*entsP) == 0) { affirm((*eP) == NULL, "bad elem pointer"); }
        spmat_count_t ents = index + 1; /* New min entry count. */
        /* Expand the array to at least twice its size (but beware of overflow): */
        if (ents <= spmat_MAX_ENTS - (*entsP))
          { ents = ents + (*entsP); }
        else
          { ents = spmat_MAX_ENTS; }
        (*eP) = realloc((*eP), ents*esz);
        affirm((*eP) != NULL, "out of mem");
        (*entsP) = ents;
      }
  }

void spmat_trim(void **eP, spmat_count_t *entsP, spmat_count_t ents, size_t esz)
  { if (ents != (*entsP))
      { if (ents == 0)
          { free((*eP)); (*eP) = NULL; }
        else
          { demand(ents <= spmat_MAX_ENTS, "too many entries");
            (*eP) = realloc((*eP), ents*esz);
            affirm((*eP) != NULL, "out of mem");
          }
        (*entsP) = ents;
      }
  }

int spmat_compare_indices
  ( spmat_index_t arow, 
    spmat_index_t acol, 
    spmat_index_t brow, 
    spmat_index_t bcol, 
    int orow, 
    int ocol
  )
  { unsigned int zrow = (orow < 0 ? -orow : orow);
    unsigned int zcol = (ocol < 0 ? -ocol : ocol);
    demand(zrow != zcol, "ambiguous sorting criterion");
    if (zrow > zcol)
      { /* Row index is more important: */
        if (arow < brow)
          { return -orow; }
        else if (arow > brow)
          { return +orow; }
        if (zcol > 0)
          { /* Break ties by column index: */
            if (acol < bcol)
              { return -ocol; }
            else if (acol > bcol)
              { return +ocol; }
            else
              { return 0; }
          }
        else
          { return 0; }
      }
    else
      { /* Column index is more important: */
        if (acol < bcol)
          { return -ocol; }
        else if (acol > bcol)
          { return +ocol; }
        if (zrow > 0)
          { /* Break ties by row index: */
            if (arow < brow)
              { return -orow; }
            else if (arow > brow)
              { return +orow; }
            else
              { return 0; }
          }
        else
          { return 0; }
      }
  }

