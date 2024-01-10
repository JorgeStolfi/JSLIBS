/* See float_image.h */
/* Last edited on 2009-08-31 11:20:25 by stolfi */ 

#include <float_array.h>

#include <indexing.h>
#include <indexing_descr.h>

#include <filefmt.h>
#include <nget.h>
#include <fget.h>
#include <affirm.h>
#include <bool.h>
#include <jsmath.h>

#include <limits.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define NA float_array_NAXES

/* IMPLEMENTATIONS */

float_array_t float_array_new ( ix_dim_t na, ix_size_t sz[], ix_order_t ixor )
  { float_array_t A;
    A.ds = ix_descr_from_sizes(na, sz, ixor);
    ix_pos_count_t NE = ix_descr_num_positions(&(A.ds));
    if (NE == 0)
      { A.el = (float *)NULL; }
    else
      { A.el = (float *)notnull(malloc(NE*sizeof(float)), "no mem"); }
    /* Paranoia: */
    (void)ix_descr_is_valid(&(A.ds), /*die:*/ TRUE);
    return A;
  }
  
void float_array_free_elems ( float_array_t *A )
  { if (A == NULL) return;
    if (A->el != NULL) { free(A->el); A->el = NULL; }
  }
  
float_array_t *float_array_new_descr ( void )
  { float_array_t *A = notnull(malloc(sizeof(float_array_t)), "no mem");
    A->el = NULL;
    return A;
  }

float_array_t *float_array_copy_descr ( float_array_t *A )
  { float_array_t *C = notnull(malloc(sizeof(float_array_t)), "no mem");
    (*C) = (*A); 
    return C; 
  }

void float_array_free_descr ( float_array_t *A )
  { if (A != NULL) { free(A); } }

float float_array_get_element ( float_array_t *A, ix_index_t ix[] )
  { ix_pos_t p = ix_descr_position(&(A->ds), ix);
    return A->el[p];
  }

float *float_array_get_element_address ( float_array_t *A, ix_index_t ix[] )
  { ix_pos_t p = ix_descr_position(&(A->ds), ix);
    return &(A->el[p]);
  }

void float_array_set_element ( float_array_t *A, ix_index_t ix[], float v )
  { ix_pos_t p = ix_descr_position(&(A->ds), ix);
    A->el[p] = v;
  }

#define float_array_file_version "2009-08-22"

void float_array_write(FILE *wr, float_array_t *A, ix_order_t ixor)
  { 
    ix_descr_t *D = &(A->ds);
    
    /* Write the header line: */
    filefmt_write_header(wr, "float_array_t", float_array_file_version);

    /* Get the effective dimension {d}: */
    ix_dim_t na = D->na;
    
    /* Find the maximum size {maxsz} along any axis: */
    ix_size_t maxsz = ix_descr_max_size(D);
    int dsz = digits(maxsz);                      /* Number of digits of max size. */
    int dix = (maxsz == 0 ? 0 : digits(maxsz-1)); /* Number of digits of max index. */

    /* Write the effective dimension {na} and the respective sizes: */
    int i; 
    fprintf(wr, "axes = %d\n", na);
    fprintf(wr, "size =");
    for (i = 0; i < na; i++) { fprintf(wr, " %*llu", dsz, D->sz[i]); }
    fprintf(wr, "\n");
    fprintf(wr, "order = %d\n", ixor);
    
    /* Write the elements, in the specified order: */
    ix_index_t ix[na];
    bool_t first = TRUE;
    if (ix_descr_indices_first(D, ix))
      { ix_pos_t p = ix_descr_position(D, ix);
        do { 
          if (! first)
            { /* Print one blank line for eack trailing 0 in {ix[0..na-1]}: */
              i = na;
              while ((i > 0) && (ix[i-1] == 0)) { fprintf(wr, "\n"); i--; }
              first = FALSE;
            }
          /* Print the index tuple: */
          for (i = 0; i < na; i++) { fprintf(wr, " %*lld", dix, ix[i]); }
          /* Print the value: */
          fprintf(wr, "  %13.7e\n", A->el[p]);
        } while (! ix_descr_next(D,ix,ixor,&p));
      }

    /* Write the footer line: */
    filefmt_write_footer(wr, "float_array_t");
    fflush(wr);
  }

float_array_t float_array_read(FILE *rd)
  {
    /* Read and check the header line: */
    filefmt_read_header(rd, "float_array_t", float_array_file_version);
    
    /* Read the effective dimension {na}: */
    uint64_t d64 = nget_uint64(rd, "axes", 10); fget_eol(rd);
    demand(d64 <= float_array_MAX_AXES, "array has too many axes");
    ix_dim_t na = d64;
    
    /* Read the size vector: */
    nget_name_eq(rd, "size");
    ix_size_t sz[na];
    int i;
    for (i = 0; i < na; i++)
      { uint64_t sz64 = (i < na ? fget_uint64(rd, 10) : 1); 
        demand(sz64 <= ix_MAX_SIZE, "size too large"); 
        sz[i] = sz64;
      }
    fget_eol(rd);
    
    /* Read the index scan order: */
    ix_order_t ixor = nget_uint64(rd, "order", 10); fget_eol(rd);
    demand(ixor <= 1, "invalid element order"); 
    
    /* Allocate the array: */
    float_array_t A = float_array_new(na, sz, ixor);
    ix_descr_t *D = &(A.ds);

    /* Read the elements, in the specified order: */
    ix_index_t ix[na];
    uint64_t ix64[na];
    if (ix_descr_indices_first(D, ix))
      { ix_pos_t p = ix_descr_position(D, ix);
        do { 
          /* Skip blank lines, if any: */
          fget_skip_formatting_chars(rd);
          /* Read and check the index tuple: */
          for (i = 0; i < na; i++) 
            { ix64[i] = fget_uint64(rd, 10);
              demand(ix64[i] == ix[i], "indices invalid or out of order"); 
            }
          /* Read the value: */
          A.el[p] = fget_double(rd);
        } while (! ix_descr_next(D,ix,ixor,&p));
      }

    /* Read the footer line: */
    filefmt_read_footer(rd, "float_array_t");
    return A;
  }

