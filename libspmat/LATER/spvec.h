/* spvec.h -- sparse vectors (one-dimensional arrays) of things */
/* Last edited on 2024-12-05 10:40:26 by stolfi */

#ifndef spvec_H
#define spvec_H

#include <stdint.h>
#include <stdlib.h>

#include <bool.h>
#include <ref.h>

#include <spvec_gen.h>

/* SPARSE VECTORS 
  
  A /sparse vector/ is a vector that contains a large proportion of
  /trivial elements/ --- elements that are equal to a client-defined
  /trivial value/. (The trivial value is usually zero in linear
  algebra, but {+oo} in graph processing, ' ' in board games, etc..)
  
  This module defines a /compressed representation/ for sparse
  vectors, namely a list of the element values together with their
  indices in the original vector. Any element not on this list is
  assumed to be trivial. The representation is handled through a
  /descriptor/, a record containing the count of entries the list, and
  a pointer to an array of (index,value) triplets.
  
  This interface also defines efficient and convenient procedures for
  extending and trimming the list of entries, as needed. Each
  expansion doubles the size of the entry list, so the total cost of
  multiple expansions is small and proportional to the size of the
  final list. */

/* DECLARING A NEW SPARSE VECTOR TYPE */

#define spvec_typedef(VEC_TYPE,PREFIX,ELEM_TYPE) \
  spvec_DEFINE_ENTRY_TYPE(VEC_TYPE,PREFIX,ELEM_TYPE); \
  vec_typedef(VEC_TYPE,PREFIX,PREFIX##_entry_t); \
  spvec_DECLARE_UNPACK(VEC_TYPE,PREFIX,ELEM_TYPE); \
  spvec_DECLARE_PACK(VEC_TYPE,PREFIX,ELEM_TYPE)
/* 
  This macro defines a new sparse vector type called {VEC_TYPE}, whose
  elements are of type {ELEM_TYPE}.
    
  This macro will also declare the data type {{PREFIX}_entry_t} (an
  (index,value) pair), and the procedures
  
    { {PREFIX}_new
      {PREFIX}_expand
      {PREFIX}_trim
      {PREFIX}_unpack
      {PREFIX}_pack
      {PREFIX}_make_desc }
      
  which are described below.
  
  For example, to define a new sparse vector type, called
  {string_spvec_t}, with elements of type {char *}, one should write
  
    { spvec_typedef(string_spvec_t, string_spvec, char *); }
    
  This macro will also declare the type {string_spvec_entry_t},
  and the procedures {string_spvec_new}, {string_spvec_expand}, etc.
  
  The names {VECTYPE} and {PREFIX} can be chosen quite arbitrarily.
  However, if {ELEM_TYPE} is some previously defined type named
  {{XXX}_t} or {XXX}, it is recommended to use {XXX}_spvec_t}
  as the {VEC_TYPE}, and {{XXX}_spvec} as the {PREFIX}. */

/* SPARSE VECTOR ENTRIES */
      
#define spvec_DEFINE_ENTRY_TYPE(VEC_TYPE,PREFIX,ELEM_TYPE) \
  typedef struct PREFIX##_entry_t { uint32_t idx; ELEM_TYPE val; } PREFIX##_entry_t
/* 
  A variable {e} of type {my_spvec_entry_t} represents an explicit
  element of the sparse vector. It contains the index {e.idx} of the
  element, and its value {e.val}. */

/* SPARSE VECTOR DESCRIPTOR */
      
#define spvec_DEFINE_VEC_TYPE(VEC_TYPE,PREFIX,ELEM_TYPE) \
  typedef struct VEC_TYPE \
    { uint32_t nen; \
      PREFIX##_entry_t *en; \
    } VEC_TYPE
/* 
  A variable {M} of type {VEC_TYPE} is a descriptor for a sparse
  vector, which has {M.nen} explicitly stored elements. The field
  {M.en} points to a heap area with space for {M.nen} explicit entries
  of type {{PREFIX}_entry_t}, namely {M.en[0..M.nen-1]}. */

/* ALLOCATION */
  
#define spvec_DECLARE_NEW(VEC_TYPE,PREFIX,ELEM_TYPE) \
  VEC_TYPE PREFIX##_new(uint32_t nen)
/*
  The call {{PREFIX}_new(nen)} allocates a new sparse vector of
  type {VEC_TYPE}, with space for {nen} explicit elements, and returns
  its descriptor. */

/* STORAGE EXPANSION */
  
#define spvec_DECLARE_EXPAND(VEC_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_expand(VEC_TYPE *M, int32_t pos)
/*
  The call {{PREFIX}_expand(&M,pos)} makes sure that the entry
  {M.en[pos]} exists. If necessary, it allocates a larger entry list
  area from the heap, copies the old entries into it, reclaims the old
  area {M.en}, sets {M.en} to the new area and {M.nen} to its size
  (which will be strictly greater than {pos}). The entry count {M.nen}
  is approximately doubled at each reallocation, to ensure total
  {O(N)} time for adding {N} entries. */

/* STORAGE TRIMMING */
  
#define spvec_DECLARE_TRIM(VEC_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_trim(VEC_TYPE *M, uint32_t nen)
/*
  The call {{PREFIX}_trim(&M,size)} reallocates the entry list {M.en}
  if necessary so that it has exactly {size} explicit entries.
  It preserves the entries {M.en[0..size-1]} and sets {M.nen = size}. */

/* UNPACKING */
  
#define spvec_DECLARE_UNPACK(VEC_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_unpack(VEC_TYPE *M, ELEM_TYPE tval, ELEM_TYPE v[], uint32_t nv)
/*
  The call {{PREFIX}_unpack(&M, tval, v, nv)} expands the sparse
  vector {M}, supplying {tval} for any omitted entry, and stores the
  result into the ordinary vector {v[0..nv-1]}.
  
  More precisely, for each explicit entry {(idx,val)) in
  {M.en[0..M.nen-1]}, the procedure sets {v[idx] = val}. All other
  elements of {v[0..nv-1]} are set to {tval}.
  
  Note that {M.nv} need not be equal to {nv}. The procedure fails with
  error if any {idx} in {M} is not in the range {0..nv-1}. */

/* PACKING */
  
#define spvec_DECLARE_PACK(VEC_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_pack(ELEM_TYPE v[], uint32_t nv, ELEM_TYPE tval, VEC_TYPE *M, uint32_t *ke)
/*
  The call {{PREFIX}_pack(v, nv, tval, &M, &ke)} gathers all
  non-trivial elements of the vector {v[0..nv-1]} and appends them to
  the list of entries of the sparse vector {M}, starting at position
  {M.en[ke]}. More precisely, for each element {v[idx]} that is not
  equal to {tval}, the procedure stores the entry {(idx,val)) to
  {M.en[ke]} and increments {ke}. 
  
  Before setting each entry {M.en[ke]}, the procedure makes sure that
  it exists, by calling {{PREFIX}_expand(&M, ke)}. */

/* DESCRIPTOR ASSEMBLY */
    
#define spvec_DECLARE_MAKE_DESC(VEC_TYPE,PREFIX,ELEM_TYPE) \
  VEC_TYPE PREFIX##_make_desc(PREFIX##_entry_t *en, uint32_t nen)
/*
  The call {{PREFIX}_make_desc(nen, en)} assembles a sparse vector
  descriptor from the given element count {nen} and
  the address {en} of the first explicit entry. The client must ensure
  that the variables {en[0..nen-1]} actually exist. The procedures
  {{PREFIX}_expand} and {{PREFIX}_trim} can be applied to the
  resulting descriptor only if the address {en} was obtained through
  {malloc}.  */

/* PROCEDURE IMPLEMENTATIONS */

#define spvec_typeimpl(VEC_TYPE,PREFIX,ELEM_TYPE) \
  spvec_IMPLEMENT_NEW(VEC_TYPE,PREFIX,ELEM_TYPE) \
  spvec_IMPLEMENT_EXPAND(VEC_TYPE,PREFIX,ELEM_TYPE) \
  spvec_IMPLEMENT_TRIM(VEC_TYPE,PREFIX,ELEM_TYPE) \
  spvec_IMPLEMENT_UNPACK(VEC_TYPE,PREFIX,ELEM_TYPE) \
  spvec_IMPLEMENT_PACK(VEC_TYPE,PREFIX,ELEM_TYPE) \
  spvec_IMPLEMENT_MAKE_DESC(VEC_TYPE,PREFIX,ELEM_TYPE) \
  extern void PREFIX##_bOgUs /* To eat the semicolon. */
/*
  This macro expands into the implementation of the sparse vector
  operations {{PROFIX}_new}, {PROFIX}_expand}, etc. It should be
  called after the corresponding call to {spvec_typedef}, which
  declares the prototypes for those procedures. */
  
#define spvec_IMPLEMENT_NEW(VEC_TYPE,PREFIX,ELEM_TYPE) \ 
  VEC_TYPE PREFIX##_new(uint32_t nen) \
    { spvec_t M = spvec_new(nen, sizeof(PREFIX##_entry_t)); \
      return (VEC_TYPE){M.nen, (ELEM_TYPE *)M.en}; \
    }

#define spvec_IMPLEMENT_EXPAND(VEC_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_expand(VEC_TYPE *M, int32_t pos) \
    { if (pos >= M->nen) \
        { spvec_expand((spvec_t *)M, pos, sizeof(PREFIX##_entry_t)); } \
    }
    
#define spvec_IMPLEMENT_TRIM(VEC_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_trim(VEC_TYPE *M, uint32_t nen) \
    { spvec_trim((spvec_t *)M, nen, sizeof(PREFIX##_entry_t)); }
    
#define spvec_IMPLEMENT_UNPACK(VEC_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_unpack(VEC_TYPE *M, ELEM_TYPE tval, ELEM_TYPE v[], uint32_t nv) \
    { uint32_t i; \
      for (i = 0; i < nv; i++) { v[i] = tval; } \
      for (i = 0; i < M->nen; i++) \
        { PREFIX##_entry_t *e = &(M->en[i]); \
          demand((e->idx >= 0) && (e->idx < nv) "invalid index"); \
          v[e->idx] = e->val; \
        } \
    }

#define spvec_IMPLEMENT_PACK(VEC_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_pack(ELEM_TYPE v[], uint32_t nv, ELEM_TYPE tval, VEC_TYPE *M, uint32_t *ke) \
    { uint32_t i; \
      for (i = 0; i < nv; i++) \
        { if (v[i] != tval) \
            { PREFIX##_expand(M, *ke); \
              PREFIX##_entry_t *e = &(M->en[*ke]); \
              e->idx = i; e->val = v[i]; \
            } \
        } \
    }

#define spvec_IMPLEMENT_MAKE_DESC(VEC_TYPE,PREFIX,ELEM_TYPE)  
  VEC_TYPE PREFIX##_desc(PREFIX##_entry_t *en, uint32_t nen) \
    { return (VEC_TYPE){ nen, en }; }
    
/* SOME USEFUL TYPED SPARSE VECTORS */

spvec_typedef(int32_t_spvec_t,    int32_t_spvec,    int32_t);
  /* Sparse Vectors of {int32_t}s. */

spvec_typedef(uint_spvec_t,   uint_spvec,   uint32_t);
  /* Sparse Vectors of {uint32_t}s. */

spvec_typedef(float_spvec_t,  float_spvec,  float);
  /* Sparse Vectors of {float}s */

spvec_typedef(double_spvec_t, double_spvec, double);
  /* Sparse Vectors of {double}s. */

spvec_typedef(bool_spvec_t,   bool_spvec,   bool_t);
  /* Sparse Vectors of {bool_t}s. */

spvec_typedef(char_spvec_t,   char_spvec,   char);
  /* Sparse Vectors of {char}s. */

spvec_typedef(ubyte_spvec_t,   ubyte_spvec, uint8_t);
  /* Sparse Vectors of {unsigned char}s. */

spvec_typedef(sbyte_spvec_t,   sbyte_spvec, int8_t);
  /* Sparse Vectors of {signed char}s. */

spvec_typedef(string_spvec_t, string_spvec, char*);
  /* Sparse Vectors of strings ({char*}s). */

spvec_typedef(ref_spvec_t, ref_spvec, ref_t);
  /* Sparse Vectors of arbitrary addresses. */

#endif
