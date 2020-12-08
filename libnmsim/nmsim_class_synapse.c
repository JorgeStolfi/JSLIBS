/* See {nmsim_class_synapse.h} */
/* Last edited on 2019-05-28 14:38:36 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <affirm.h>
#include <nget.h>
#include <fget.h>
#include <vec.h>
#include <filefmt.h>
#include <jsmath.h>
#include <jsrandom.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>

#include <nmsim_firing_func.h>
#include <nmsim_class_synapse.h>
 
nmsim_class_synapse_t *nmsim_class_synapse_new(double W_avg, double W_dev)
  {
    nmsim_class_synapse_t *sclass = notnull(malloc(sizeof(nmsim_class_synapse_t)), "no mem");
    (*sclass) = (nmsim_class_synapse_t) { .W_avg = W_avg, .W_dev = W_dev };
    return sclass;
  }

void nmsim_class_synapse_free(nmsim_class_synapse_t *sclass)
  { 
    if (sclass != NULL) { free(sclass); }
  }

void nmsim_class_synapse_write(FILE *wr, nmsim_class_synapse_t *sclass)
  {
    char *ind1 = "      ";   /* Indent of header lines. */
    char *ind2 = "        "; /* Indent of parameter lines. */
    
    fputs(ind1, wr);
    filefmt_write_header(wr, nmsim_class_synapse_FILE_TYPE, nmsim_class_synapse_VERSION);
    
    nmsim_write_double_param(wr, ind2, "W_avg", sclass->W_avg, nmsim_write_KW_PREC, TRUE,  TRUE, FALSE);
    nmsim_write_double_param(wr, ind2, "W_dev", sclass->W_dev, nmsim_write_KW_PREC, FALSE, TRUE, FALSE);
    
    fputs(ind1, wr);
    filefmt_write_footer(wr, nmsim_class_synapse_FILE_TYPE);
    fflush(wr);
  }

nmsim_class_synapse_t *nmsim_class_synapse_read(FILE *rd)
  {
    /* Read header line: */
    filefmt_read_header(rd, nmsim_class_synapse_FILE_TYPE, nmsim_class_synapse_VERSION);
    
    /* Read the parameters: */
    double W_avg  = nmsim_read_double_param(rd, "W_avg", -INF, +INF);
    double W_dev  = nmsim_read_double_param(rd, "W_dev", 0.0, +INF);
    /* Build record: */
    nmsim_class_synapse_t *sclass = nmsim_class_synapse_new(W_avg, W_dev);

    /* Read footer line: */
    filefmt_read_footer(rd, nmsim_class_synapse_FILE_TYPE);
    
    return sclass;
  }

void nmsim_class_synapse_show(FILE *wr, char *pref, nmsim_class_synapse_t *sclass, char *suff)
  { 
    if (pref != NULL) { fputs(pref, wr); }
    fprintf(wr, "{ W_avg = %+.4f W_dev = %.4f", sclass->W_avg, sclass->W_dev);
    fprintf(wr, " }");
    if (suff != NULL) { fputs(suff, wr); }
  }

nmsim_class_synapse_t* nmsim_class_synapse_throw(void)
  { 
    double W_avg = nmsim_throw_double(-10.0, +10.0);
    double W_dev = nmsim_throw_double(0.0, 0.25*fabs(W_avg) + 2.5);

    nmsim_class_synapse_t *sclass = nmsim_class_synapse_new(W_avg, W_dev);
      
    return sclass;
  }

void nmsim_class_synapse_compare(nmsim_class_synapse_t *sclass_read, nmsim_class_synapse_t *sclass_orig)
  { 
    nmsim_compare_double_param("W_avg", sclass_read->W_avg, sclass_orig->W_avg, nmsim_write_KW_PREC, TRUE, FALSE);
    nmsim_compare_double_param("W_dev", sclass_read->W_dev, sclass_orig->W_dev, nmsim_write_KW_PREC, TRUE, FALSE);
  }    

vec_typeimpl(nmsim_class_synapse_ref_vec_t,nmsim_class_synapse_ref_vec,nmsim_class_synapse_ref_t);
