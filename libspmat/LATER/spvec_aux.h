/* Auxiliary definitions for sparse vectors */
/* Last edited on 2008-03-29 14:27:04 by stolfi */

#ifndef spvec_aux_H
#define spvec_aux_H

#include <stdlib.h>
#include <bool.h>
#include <ref.h>
#include <stdint.h>

/* UNTYPED SPARSE VECTORS 
  
  The types and procedures in this interface may be used directly by clients,
  but their main purpose is to help implement the functions created by the
  macro {spvec_typeimpl} of {spvec.h}. */

/* DESCRIPTORS */

typedef struct spvec_t /* Untyped sparse vector descriptor. */
  { uint32_t nen;   /* Number of explicit entries. */
    void *en;       /* Pointer to the first explicit entry. */
  } spvec_t;
  /* An {spvec_t} {M} describes a uni-dimensional array of {M.nen} elements
    of the same size and type, stored in consecutive memory locations,
    starting at address {M.en}.  The type and size of the elements 
    must be provided by the client. */

/* PROCEDURES */

spvec_t spvec_new(int nen, size_t esz);
  /* Allocates a new sparse vector with {nen} elements of size {esz}.
    Bombs out if there is no space for the request. 
    If {nen == 0}, the result has {en == NULL}. */

void spvec_expand(spvec_t *M, uint32_t pos, size_t esz);
  /* Makes sure that element {M->en[pos]} exists, reallocating and
    copying the array {*M->en} if {pos >= M->nen}. If that happens,
    the new {M->nen} will be about twice as big as the old one. */

void spvec_trim(spvec_t *M, uint32_t nen, size_t esz);
  /* Makes sure that {M.nen == nen}, reallocating and copying 
    the  array {*M->en} if {M.nen != nen}. If {nen == 0}, the 
    result will have {en == NULL}. */

#endif
