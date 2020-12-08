/* Self-bounded vectors (one-dimensional arrays) of things */
/* Last edited on 2019-01-11 19:54:22 by jstolfi */

#ifndef vec_H
#define vec_H

#include <stdlib.h>
#include <stdint.h>

#include <bool.h>
#include <ref.h>

/* SELF-BOUNDED VECTORS

  This interface defines self-bounded vectors (unidimensional arrays)
  with elements of any specified type. Unlike standard C vectors, the
  number of elements is part of the vector's representation. Thus,
  when passing a self-bounded vector to a procedure, there is no need
  to pass its size as as a separate parameter.
  
  This interface also defines efficient and convenient procedures for
  extending and trimming vectors as needed. Expansion doubles the size
  of the vector, so the total cost of the expansion is small and
  proportional to the final size of the vector.
  
  In many contexts, these vectors can replace linked lists, with
  advantages. For example, consider this bit of code that generates an
  undetermined number of {double} values, one value at a time, and
  saves them into a self-bounded vector (a {double_vec_t}):
  
    double_vec_t x = double_vec_new(100);
    { int nx = 0;
      while(! finished(...))
        { double z = generate_next_value(...); 
          double_vec_expand(&x, nx);
          x.e[nx] = z; nx++;
        }
      double_vec_trim(&x, nx);
    }
    
  The saved values can then be printed with 
  
    for (i = 0; i < x.ne; i++) { printf("%f\n", x.e[i]); }

  */

/* DECLARING A NEW SPARSE VECTOR TYPE */

#define vec_typedef(VEC_TYPE,PREFIX,ELEM_TYPE) \
  vec_DEFINE_VEC_TYPE(VEC_TYPE,PREFIX,ELEM_TYPE); \
  vec_DECLARE_NEW(VEC_TYPE,PREFIX,ELEM_TYPE); \
  vec_DECLARE_EXPAND(VEC_TYPE,PREFIX,ELEM_TYPE); \
  vec_DECLARE_TRIM(VEC_TYPE,PREFIX,ELEM_TYPE); \
  vec_DECLARE_MAKE_DESC(VEC_TYPE,PREFIX,ELEM_TYPE)
/* 
  This macro defines a new vector type called {VEC_TYPE}, whose
  elements are of type {ELEM_TYPE}.
    
  This macro will also declare prototypes for the functions
  
    { {PREFIX}_new
      {PREFIX}_expand
      {PREFIX}_trim
      {PREFIX}_make_desc }
      
  which are described below.
  
  For example, the macro call 
  
    { vec_typedef(string_vec_t, string_vec, char *); } 
    
  will declare the type {string_vec_t}, a self-bounded vector with
  elements of type {char *}; and will also declare the procedures
  {string_vec_new}, {string_vec_expand}, etc.
  
  The names {VEC_TYPE} and {PREFIX} can be chosen quite arbitrarily.
  However, if {ELEM_TYPE} is some previously defined type named
  {{XXX}_t} or {XXX}, it is recommended to use {XXX}_vec_t}
  as the {VEC_TYPE}, and {{XXX}_vec} as the {PREFIX}. */

#define vec_typeimpl(VEC_TYPE,PREFIX,ELEM_TYPE) \
  vec_IMPLEMENT_NEW(VEC_TYPE,PREFIX,ELEM_TYPE) \
  vec_IMPLEMENT_EXPAND(VEC_TYPE,PREFIX,ELEM_TYPE) \
  vec_IMPLEMENT_TRIM(VEC_TYPE,PREFIX,ELEM_TYPE) \
  vec_IMPLEMENT_MAKE_DESC(VEC_TYPE,PREFIX,ELEM_TYPE) \
  extern void PREFIX##_bOgUs /* To eat the semicolon. */
/*
  This macro expands into the implementation of the sparse vector
  operations {{PREFIX}_new}, {PREFIX}_expand}, etc. It should be
  called after the corresponding call to {vec_typedef}, which
  declares the prototypes for those procedures. */

/* VECTOR DESCRIPTOR */
      
typedef uint32_t vec_size_t; /* Type for the number of elements in a vector. */
#define vec_MAX_SIZE (2147483647LU)
  /* Maximum allocated elements in a vector (2^32-1). */

#define vec_DEFINE_VEC_TYPE(VEC_TYPE,PREFIX,ELEM_TYPE) \
  typedef struct VEC_TYPE \
    { vec_size_t ne; \
      ELEM_TYPE *e; \
    } VEC_TYPE
/*
  Defines the type {VEC_TYPE}. A variable {v} of that type is a
  descriptor for a vector with {v.ne} elements of type {ELEM_TYPE},
  namely {v.e[0..v.ne-1]}. */

typedef int32_t vec_index_t;
#define vec_MAX_INDEX (vec_MAX_SIZE - 1)
  /* We need at least this much in order to index every element, and no more
    than this since we cannot create arrays bigger than that.
    Note that {vec_MAX_INDEX + 1} does not overflow a {vec_index_t}. */


/* ALLOCATION */
  
#define vec_DECLARE_NEW(VEC_TYPE,PREFIX,ELEM_TYPE) \
  VEC_TYPE PREFIX##_new(vec_size_t ne)
/*
  This macro declares the function {{PREFIX}_new}. The call
  {{PREFIX}_new(ne)} allocates a new vector of type {VEC_TYPE}, with
  space for {ne} elements, and returns its descriptor. */

/* STORAGE EXPANSION */
  
#define vec_DECLARE_EXPAND(VEC_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_expand(VEC_TYPE *vp, vec_index_t index)
/*
  This macro declares the function {{PREFIX}_expand}. The call
  {{PREFIX}_expand(&v,pos)} makes sure that the element {v.e[pos]}
  exists. If necessary, it allocates a larger area from the heap,
  copies the old elements of {v} into it, reclaims the old area {v.e},
  and finally sets {v.e} to the new area and {v.ne} to its size (which
  will be strictly greater than {pos}). The element count {v.ne} is
  approximately doubled at each reallocation, to ensure total {O(N)}
  time for adding {N} elements. */

/* STORAGE TRIMMING */
  
#define vec_DECLARE_TRIM(VEC_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_trim(VEC_TYPE *vp, vec_size_t ne)
/*
  This macro declares the function {{PREFIX}_trim}. The call
  {{PREFIX}_trim(&v,size)} reallocates the element area {v.e} if
  necessary so that it has exactly {size} elements. It preserves the
  elements {v.e[0..size-1]} and sets {v.ne = size}. */

/* DESCRIPTOR ASSEMBLY */
    
#define vec_DECLARE_MAKE_DESC(VEC_TYPE,PREFIX,ELEM_TYPE) \
  VEC_TYPE PREFIX##_make_desc(ELEM_TYPE *e, vec_size_t ne)
/*
  This macro declares the function {{PREFIX}_make_desc}. The call
  {{PREFIX}_make_desc(e, ne)} assembles a vector descriptor from the
  given element count {ne} and the address {e} of the first explicit
  entry. The client must ensure that the variables {e[0..ne-1]}
  actually exist. The procedures {{PREFIX}_expand} and {{PREFIX}_trim}
  can be applied to the resulting descriptor only if the address {e}
  was obtained through {malloc}. */

/* PROCEDURE IMPLEMENTATIONS */
  
#define vec_IMPLEMENT_NEW(VEC_TYPE,PREFIX,ELEM_TYPE) \
  VEC_TYPE PREFIX##_new(vec_size_t ne) \
    { void *e = vec_alloc(ne, sizeof(ELEM_TYPE)); \
      return (VEC_TYPE){ne, (ELEM_TYPE *)e}; \
    }

#define vec_IMPLEMENT_EXPAND(VEC_TYPE,PREFIX,ELEM_TYPE) \
void PREFIX##_expand(VEC_TYPE *vp, vec_index_t index) \
  { if (index >= vp->ne) \
      { vec_expand(&(vp->ne), (void**)&(vp->e), index, sizeof(ELEM_TYPE)); } \
  }

#define vec_IMPLEMENT_TRIM(VEC_TYPE,PREFIX,ELEM_TYPE)\
  void PREFIX##_trim(VEC_TYPE *vp, vec_size_t ne) \
    { vec_trim(&(vp->ne), (void**)&(vp->e), ne, sizeof(ELEM_TYPE)); }

#define vec_IMPLEMENT_MAKE_DESC(VEC_TYPE,PREFIX,ELEM_TYPE)\
  VEC_TYPE PREFIX##_make_desc(ELEM_TYPE *e, vec_size_t ne) \
    { return (VEC_TYPE){ne, e}; }
  
/* SOME USEFUL TYPED VECTORS */

vec_typedef(int_vec_t,    int_vec,     int);
  /* Vectors of {int}s. */

vec_typedef(uint_vec_t,   uint_vec,    unsigned int);
  /* Vectors of {unsigned int}s. */

vec_typedef(char_vec_t,   char_vec,    char);
  /* Vectors of {char}s. */

vec_typedef(bool_vec_t,   bool_vec,    bool_t);
  /* Vectors of {bool_t}s. */

vec_typedef(float_vec_t,  float_vec,   float);
  /* Vectors of {float}s */

vec_typedef(double_vec_t, double_vec,  double);
  /* Vectors of {double}s. */

vec_typedef(ref_vec_t,    ref_vec,     ref_t);
  /* Vectors of arbitrary addresses. */

vec_typedef(string_vec_t, string_vec,  char*);
  /* Vectors of strings ({char*}s). */

vec_typedef(int8_vec_t,   int8_vec,    int8_t);
  /* Vectors of {int8_t}s (signed bytes, {signed char}s). */

vec_typedef(int16_vec_t,  int16_vec,   int16_t);
  /* Vectors of {int16_t}s (signed half-words, {signed short int}s). */

vec_typedef(int32_vec_t,  int32_vec,   int32_t);
  /* Vectors of {int32_t}s. */

vec_typedef(int64_vec_t,  int64_vec,   int64_t);
  /* Vectors of {int64_t}s. */

vec_typedef(uint8_vec_t,  uint8_vec,   uint8_t);
  /* Vectors of {unit8_t}s (unsigned bytes, {unsigned char}s). */

vec_typedef(uint16_vec_t, uint16_vec,  uint16_t);
  /* Vectors of {unit16_t}s (unsigned half-words, {unsigned short int}s. */

vec_typedef(uint32_vec_t, uint32_vec,  uint32_t);
  /* Vectors of {uint32_t}s. */

vec_typedef(uint64_vec_t, uint64_vec,  uint64_t);
  /* Vectors of {uint64_t}s. */
/* GENERIC PROCEDURES FOR SELF-BOUNDED VECTORS 
  
  The functions in this section may be called directly by clients,
  but their main purpose is to implement the 
  functions {my_vec_new}, {my_vec_expand}, and {my_vec_trim}
  defined by {vec_typedef} and {vec_typeimpl}. */

void *vec_alloc(vec_size_t ne, size_t esz);
  /* Allocates a new storage area with space for {ne} elements of
    size {esz}. Bombs out if there is no space for the request. If
    {ne == 0}, the result is {NULL}. */

void vec_expand(vec_size_t *nep, void **ep, vec_index_t index, size_t esz);
  /* Makes sure that element {(*ep)[index]} exists, reallocating and
    copying the array {**ep} if {index >= (*nep)}. If that happens,
    sets {(*nep) to the new element count; which will be strictly
    greater than {index}, and about twice as big as the old
    {(*nep)}. */

void vec_trim(vec_size_t *nep, void **ep, vec_size_t ne, size_t esz);
  /* Makes sure that {(*nep) == ne}, reallocating and copying 
    the  array {**ep} if {(*nep) != ne}. If {ne == 0}, sets 
    {*ep} to {NULL}. */

#endif
