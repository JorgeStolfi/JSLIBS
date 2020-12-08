/* See {raut_io.h}. */
/* Last edited on 2019-04-09 13:02:20 by jstolfi */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <affirm.h>
#include <filefmt.h>
#include <fget.h>
#include <nget.h>
#include <bool.h>

#include <rdag.h>
#include <rdag_io.h>
#include <raut.h>
#include <raut_def.h>

#include <raut_io.h>


#define raut_FILE_TYPE "raut_t"
#define raut_FILE_VERSION "2009-10-28"
    
void raut_write(FILE *wr, raut_t *A)
  { 
    /* Write the file header: */
    filefmt_write_header(wr, raut_FILE_TYPE, raut_FILE_VERSION);

    /* Write the root state: */
    raut_state_t root = raut_root_get(A);
    fprintf(wr, "root = %u %u\n", root.ac, root.nd);

    /* Write the doc string: */
    char *doc = raut_doc_get(A);
    int ind = 0; /* Comment indentation. */
    filefmt_write_comment(wr, doc, ind, '|');

    /* Write the underlying dag: */
    rdag_t *D = raut_dag_get(A);
    rdag_write(wr, D);

    /* Write the file footer: */
    filefmt_write_footer(wr, raut_FILE_TYPE);
    fflush(wr);
  }

raut_t *raut_read(FILE *rd)
  { 
    /* Read and check the file header: */
    filefmt_read_header(rd, raut_FILE_TYPE, raut_FILE_VERSION);
    
    /* Get the root state: */
    uint32_t root_ac = nget_uint32(rd, "root", 10);
    rdag_node_t root_nd = fget_uint32(rd, 10);
    demand(root_ac <= 1, "invalid root class");
    
    /* Read the doc string: */
    char *doc = filefmt_read_comment(rd, '|');
    
    /* Get the dag: */
    rdag_t *D = rdag_read(rd);
    demand(root_nd <= rdag_node_max(D), "invalid root node");
    
    /* Put it together: */
    raut_t *A = (raut_t *)notnull(malloc(sizeof(raut_t)), "no mem");
    A->root = (raut_state_t){ .ac = root_ac, .nd = root_nd };
    A->doc = doc;
    A->D = D;
    
    /* Read and check the file footer: */
    filefmt_read_footer(rd, raut_FILE_TYPE);

    return A;
  }

void raut_state_debug(FILE *wr, char* pref, raut_t *A, raut_state_t v, char *suff)
  {
    fprintf(wr, "%s%u", pref, v.ac);
    if (v.nd == rdag_node_NULL) 
      { fprintf(wr, ":[-- -- -- --]"); }
    else
      { rdag_node_data_t dt; 
        rdag_node_data_get(raut_dag_get(A), v.nd, &dt); 
        fprintf(wr, ":[%u %u %u %u]", dt.f_link, dt.i_mark, dt.o_mark, dt.p_link);
      }
    fprintf(wr, "%s", suff);
  }

