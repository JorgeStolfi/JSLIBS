#ifndef array_def_H
#define array_def_H
/* The ugly entrails of {array.h}. */

#define array_def_H_COPYRIGHT "Copyright © 2008 by J. Stolfi, UNICAMP"
/* Created on 2008-07-19 by J.Stolfi, UNICAMP */
/* Last edited on 2023-03-18 11:03:37 by stolfi */

/* These inclusions are necessary if this file is included or compiled on its own: */
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <array.h>

/* DEFINITION MACROS */

#define array_DEFINE_BASIC_TYPES(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  typedef ix_dim_t PREFIX##_dim_t; \
  typedef ix_size_t PREFIX##_size_t; \
  typedef ix_index_t PREFIX##_index_t; \
  typedef ix_step_t PREFIX##_step_t; \
  typedef ix_count_t PREFIX##_count_t; \
  typedef ix_pos_t PREFIX##_pos_t; \
  typedef ELEM_TYPE PREFIX##_elem_t

#define array_DEFINE_ARRAY_TYPE(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  typedef struct ARRAY_TYPE \
    { ix_descr_t ds; \
      ELEM_TYPE *e; \
    } ARRAY_TYPE

#define array_DECLARE_new(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  ARRAY_TYPE PREFIX##_new ( ix_dim_t na, ix_size_t sz[] )

#define array_DECLARE_copy(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  ARRAY_TYPE PREFIX##_copy( ARRAY_TYPE *A )

#define array_DECLARE_new_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  ARRAY_TYPE *PREFIX##_new_descr ( void )
 
#define array_DECLARE_copy_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  ARRAY_TYPE *PREFIX##_copy_descr ( ARRAY_TYPE *A )

#define array_DECLARE_free_elems(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_free_elems ( ARRAY_TYPE *A )

#define array_DECLARE_free_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_free_descr ( ARRAY_TYPE *A )

#define array_DECLARE_get_elem(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  ELEM_TYPE PREFIX##_get_elem ( ARRAY_TYPE *A, ix_index_t ix[] )

#define array_DECLARE_get_elem_position(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  PREFIX##_pos_t PREFIX##_get_elem_position ( ARRAY_TYPE *A, ix_index_t ix[] )

#define array_DECLARE_set_elem(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_set_elem ( ARRAY_TYPE *A, ix_index_t ix[], ELEM_TYPE v )

#define array_DECLARE_assign(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_assign(ARRAY_TYPE *A, ARRAY_TYPE *B )

#define array_DECLARE_get_size(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_get_size ( ARRAY_TYPE *A, ix_size_t sz[] )
  
#define array_DECLARE_check_size(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  void PREFIX##_check_size ( ARRAY_TYPE *A, ix_dim_t na, ix_size_t sz[] )

#define array_TYPEDEF(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DEFINE_BASIC_TYPES(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DEFINE_ARRAY_TYPE(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_new(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_copy(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_new_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_copy_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_free_elems(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_free_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_get_elem(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_get_elem_position(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_set_elem(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_assign(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_get_size(ARRAY_TYPE,PREFIX,ELEM_TYPE); \
  array_DECLARE_check_size(ARRAY_TYPE,PREFIX,ELEM_TYPE)

#define array_IMPLEMENT_new(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_new(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
    { ARRAY_TYPE A; \
      ix_descr_t *DA = &(A.ds); \
      (*DA) = ix_descr_from_sizes(na, sz);            \
      ix_count_t ne = ix_descr_num_positions(DA); \
      if (ne == 0) \
        { A.e = (ELEM_TYPE *)NULL; } \
      else \
        { A.e = (ELEM_TYPE *)notnull(malloc(ne*sizeof(ELEM_TYPE)), "no mem"); } \
      /* Paranoia: */ \
      (void)ix_descr_is_valid(DA, /*die:*/ TRUE); \
      return A; \
    }

#define array_IMPLEMENT_copy(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_copy(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
    { ARRAY_TYPE C = PREFIX##_new(A->ds.na, A->ds.sz); \
      /* !!! Should check for zero steps !!! */ \
      PREFIX##_assign(&C, A);                  \
      return C; \
    }

#define array_IMPLEMENT_free_elems(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_free_elems(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
    { if (A == NULL) return; \
      if (A->e != NULL) { free(A->e); A->e = NULL; } \
    }

#define array_IMPLEMENT_new_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_new_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
    { ARRAY_TYPE *A = notnull(malloc(sizeof(ARRAY_TYPE)), "no mem"); \
      A->e = NULL; \
      return A; \
    }

#define array_IMPLEMENT_copy_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_copy_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
    { ARRAY_TYPE *C = notnull(malloc(sizeof(ARRAY_TYPE)), "no mem"); \
      (*C) = (*A); \
      return C; \
    } \

#define array_IMPLEMENT_free_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_free_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
    { if (A != NULL) { free(A); } }

#define array_IMPLEMENT_get_elem(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_get_elem(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
    { ix_pos_t p = ix_descr_position(&(A->ds), ix); \
      return A->e[p]; \
    }

#define array_IMPLEMENT_get_elem_position(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_get_elem_position(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
    { ix_pos_t p = ix_descr_position(&(A->ds), ix); \
      return p; \
    }

#define array_IMPLEMENT_set_elem(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_set_elem(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
    { ix_pos_t p = ix_descr_position(&(A->ds), ix); \
      A->e[p] = v; \
    }

#define array_IMPLEMENT_get_size(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_get_size(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
    { ix_descr_t *DA = &(A->ds); \
      int32_t ia; \
      for (ia = 0; ia < DA->na; ia++) { sz[ia] = DA->sz[ia]; } \
    }

#define array_IMPLEMENT_check_size(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_check_size(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
    { ix_descr_t *DA = &(A->ds); \
      int32_t ia; \
      demand(DA->na == na, "wrong number of indices"); \
      for (ia = 0; ia < DA->na; ia++) \
        { demand(DA->sz[ia] == sz[ia], "wrong size"); } \
    }

#define array_IMPLEMENT_assign(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_DECLARE_assign(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
    { auto bool_t assign_elem ( const ix_index_t ix[], ix_pos_t pA, ix_pos_t pB, ix_pos_t pC ); \
      bool_t assign_elem ( const ix_index_t ix[], ix_pos_t pA, ix_pos_t pB, ix_pos_t pC ) \
        { A->e[pA] = B->e[pB]; return FALSE; } \
      (void)ix_descr_enum(&assign_elem, ix_order_L, FALSE, &(A->ds), &(B->ds), NULL); \
     }

#define array_TYPEIMPL(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_IMPLEMENT_new(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_IMPLEMENT_copy(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_IMPLEMENT_new_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_IMPLEMENT_copy_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_IMPLEMENT_free_elems(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_IMPLEMENT_free_descr(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_IMPLEMENT_get_elem(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_IMPLEMENT_get_elem_position(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_IMPLEMENT_set_elem(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_IMPLEMENT_assign(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_IMPLEMENT_get_size(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  array_IMPLEMENT_check_size(ARRAY_TYPE,PREFIX,ELEM_TYPE) \
  extern void PREFIX##_bOgUs /* To eat the semicolon. */

#endif
