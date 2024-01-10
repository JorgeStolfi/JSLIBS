/* See {nmsim_class_net.h} */
/* Last edited on 2019-06-18 11:21:59 by jstolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <affirm.h>
#include <bool.h>
#include <fget.h>
#include <nget.h>
#include <filefmt.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_class_neuron.h>
#include <nmsim_class_synapse.h>

#include <nmsim_class_net.h>

nmsim_class_net_t *nmsim_class_net_new(void)
  {
    nmsim_class_net_t *cnet = notnull(malloc(sizeof(nmsim_class_net_t)), "no mem");
    (*cnet) = (nmsim_class_net_t)
      { .nnc = 0, .nsc = 0,
        .nclass = nmsim_class_neuron_ref_vec_new(20),
        .sclass = nmsim_class_synapse_ref_vec_new(20)
      };
    return cnet;
  }

void nmsim_class_net_free(nmsim_class_net_t *cnet)
  {
    if (cnet != NULL)
      { if (cnet->nclass.e != NULL) 
          { for (nmsim_class_neuron_ix_t inc = 0; inc < cnet->nnc; inc++)
              { if (cnet->nclass.e[inc] != NULL) { free(cnet->nclass.e[inc]); } }
            free(cnet->nclass.e);
          }
        if (cnet->sclass.e != NULL) 
          { for (nmsim_class_synapse_ix_t isc = 0; isc < cnet->nsc; isc++)
              { if (cnet->sclass.e[isc] != NULL) { free(cnet->sclass.e[isc]); } }
            free(cnet->sclass.e);
          }
        free(cnet);
      }
  }
  
nmsim_class_neuron_ix_t nmsim_class_net_add_neuron_class
  ( nmsim_class_net_t *cnet, 
    nmsim_class_neuron_t *nclass
  )
  { nmsim_class_neuron_ix_t inc = cnet->nnc;
    demand(inc <= nmsim_class_neuron_ix_MAX, "too many neuron classes");
    nmsim_class_neuron_ref_vec_expand(&(cnet->nclass), inc);
    cnet->nclass.e[inc] = nclass;
    (cnet->nnc)++;
    return inc;
  }
  
nmsim_class_synapse_ix_t nmsim_class_net_add_synapse_class
  ( nmsim_class_net_t *cnet, 
    nmsim_class_synapse_t *sclass
  )
  { nmsim_class_synapse_ix_t isc = cnet->nsc;
    demand(isc <= nmsim_class_synapse_ix_MAX, "too many synapse classes");
    nmsim_class_synapse_ref_vec_expand(&(cnet->sclass), isc);
    cnet->sclass.e[isc] = sclass;
    (cnet->nsc)++;
    return isc;
  }

void nmsim_class_net_write(FILE *wr, nmsim_class_net_t *cnet, double timeStep)
  { 
    char *ind1 = "    ";    /* Indentation for whole description. */
    char *ind2 = "      ";  /* Indentation for items. */
    fputs(ind1, wr);
    filefmt_write_header(wr, nmsim_class_net_FILE_TYPE, nmsim_class_net_VERSION);
    fprintf(wr, "%sneuron_classes = %d\n", ind2, cnet->nnc);
    for (nmsim_class_neuron_ix_t inc = 0; inc < cnet->nnc; inc++)
      { fprintf(wr, "%sclass = %d\n", ind2, inc);
        nmsim_class_neuron_write(wr, cnet->nclass.e[inc], timeStep);
      }
    fprintf(wr, "%ssynapse_classes = %d\n", ind2, cnet->nsc);
    for (nmsim_class_synapse_ix_t isc = 0; isc < cnet->nsc; isc++)
      { fprintf(wr, "%sclass = %d\n", ind2, isc);
        nmsim_class_synapse_write(wr, cnet->sclass.e[isc]);
      }
    fputs(ind1, wr);
    filefmt_write_footer(wr, nmsim_class_net_FILE_TYPE);

    fflush(wr);
  }

nmsim_class_net_t *nmsim_class_net_read(FILE *rd, double timeStep)
  {
    /* Read header line: */
    filefmt_read_header(rd, nmsim_class_net_FILE_TYPE, nmsim_class_net_VERSION);
    
    /* Create the empty network description: */
    nmsim_class_net_t *cnet = nmsim_class_net_new();
    
    /* Read the neuron classes: */
    nmsim_class_neuron_count_t nnc = (nmsim_class_neuron_count_t)
      nmsim_read_int64_param(rd, "neuron_classes", 0, nmsim_class_neuron_count_MAX);
    for (nmsim_class_neuron_ix_t inc = 0; inc < nnc; inc++)
      { (void)nmsim_read_int64_param(rd, "class", inc, inc);
        nmsim_class_neuron_t *nclass = nmsim_class_neuron_read(rd, timeStep);
        nmsim_class_neuron_ix_t jnc = nmsim_class_net_add_neuron_class(cnet, nclass);
        assert(jnc == inc);
      }
    
    /* Read the synapse classes: */
    nmsim_class_synapse_count_t nsc = (nmsim_class_synapse_count_t)
      nmsim_read_int64_param(rd, "synapse_classes", 0, nmsim_class_synapse_count_MAX);
    for (nmsim_class_synapse_ix_t isc = 0; isc < nsc; isc++)
      { (void)nmsim_read_int64_param(rd, "class", isc, isc);
        nmsim_class_synapse_t *sclass = nmsim_class_synapse_read(rd);
        nmsim_class_synapse_ix_t jsc = nmsim_class_net_add_synapse_class(cnet, sclass);
        assert(jsc == isc);
      }

    /* Read footer line: */
    filefmt_read_footer(rd, nmsim_class_net_FILE_TYPE);

    return cnet;
  }

nmsim_class_net_t *nmsim_class_net_throw
  ( nmsim_class_neuron_count_t nnc, 
    nmsim_class_synapse_count_t nsc
  )
  {
    demand((nnc >= 0) && (nnc <= nmsim_class_neuron_count_MAX), "invalid neuron class count");
    demand((nsc >= 0) && (nsc <= nmsim_class_synapse_count_MAX), "invalid synapse class count");
    
    nmsim_class_net_t *cnet = nmsim_class_net_new();
    
    /* Add neuron classes: */
    for (nmsim_class_neuron_ix_t inc = 0; inc < nnc; inc++)
      { nmsim_class_neuron_t *nclass = nmsim_class_neuron_throw();
        nmsim_class_neuron_ix_t jnc = nmsim_class_net_add_neuron_class(cnet, nclass);
        assert(jnc == inc);
      }
    
    /* Add synapse classes: */
    for (nmsim_class_synapse_ix_t isc = 0; isc < nsc; isc++)
      { nmsim_class_synapse_t *sclass = nmsim_class_synapse_throw();
        nmsim_class_synapse_ix_t jsc = nmsim_class_net_add_synapse_class(cnet, sclass);
        assert(jsc == isc);
      }
    return cnet;
  }

void nmsim_class_net_compare(nmsim_class_net_t *cnet_read, nmsim_class_net_t *cnet_orig)
  { 
    /* Compare neuron classes: */
    nmsim_compare_int64_param("neuron_classes", cnet_read->nnc, cnet_orig->nnc);
    assert(cnet_read->nnc <= cnet_read->nclass.ne);
    for (nmsim_class_neuron_ix_t inc = 0; inc < cnet_orig->nnc; inc++)
      { nmsim_class_neuron_compare(cnet_read->nclass.e[inc], cnet_orig->nclass.e[inc]); }

    /* Compare synapse classes: */
    nmsim_compare_int64_param("synapse_classes", cnet_read->nsc, cnet_orig->nsc);
    assert(cnet_read->nsc <= cnet_read->sclass.ne);
    for (nmsim_class_synapse_ix_t isc = 0; isc < cnet_orig->nsc; isc++)
      {  nmsim_class_synapse_compare(cnet_read->sclass.e[isc], cnet_orig->sclass.e[isc]); }
  }    

