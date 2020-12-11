/* See {nmsim_class_net.h} */
/* Last edited on 2020-12-11 17:56:38 by jstolfi */

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

nmsim_class_net_t *nmsim_class_net_new
  ( nmsim_class_neuron_count_t nnc, 
    nmsim_class_synapse_count_t nsc
  )
  {
    nmsim_class_net_t *cnet = notnull(malloc(sizeof(nmsim_class_net_t)), "no mem");
    demand((nnc >= 1) && (nnc <= nmsim_class_neuron_count_MAX), "invalid number of neuron classes");
    demand((nsc >= 0) && (nsc <= nmsim_class_synapse_count_MAX), "invalid number of synapse classes");
    (*cnet) = (nmsim_class_net_t)
      { .nnc = nnc, 
        .nsc = nsc,
        .nclass = notnull(malloc(nnc*sizeof(nmsim_class_neuron_t *)), "no mem"),
        .sclass = (nsc == 0 ? NULL : notnull(malloc(nsc*sizeof(nmsim_class_synapse_t *)), "no mem"))
      };
    for (nmsim_class_neuron_ix_t inc = 0; inc < nnc; inc++) { cnet->nclass[inc] = NULL; }
    for (nmsim_class_synapse_ix_t isc = 0; isc < nsc; isc++) { cnet->sclass[isc] = NULL; }
    return cnet;
  }

void nmsim_class_net_free(nmsim_class_net_t *cnet)
  {
    if (cnet != NULL)
      { if (cnet->nclass != NULL) 
          { for (nmsim_class_neuron_ix_t inc = 0; inc < cnet->nnc; inc++)
              { if (cnet->nclass[inc] != NULL) { free(cnet->nclass[inc]); } }
            free(cnet->nclass);
          }
        if (cnet->sclass != NULL) 
          { for (nmsim_class_synapse_ix_t isc = 0; isc < cnet->nsc; isc++)
              { if (cnet->sclass[isc] != NULL) { free(cnet->sclass[isc]); } }
            free(cnet->sclass);
          }
        free(cnet);
      }
  }

void nmsim_class_net_write(FILE *wr, nmsim_class_net_t *cnet, double timeStep)
  { 
    char *ind1 = "    ";    /* Indentation for whole description. */
    char *ind2 = "      ";  /* Indentation for items. */
    fputs(ind1, wr);
    filefmt_write_header(wr, nmsim_class_net_FILE_TYPE, nmsim_class_net_VERSION);

    fprintf(wr, "%sneuron_classes = %d\n", ind2, cnet->nnc);
    fprintf(wr, "%ssynapse_classes = %d\n", ind2, cnet->nsc);

    for (nmsim_class_neuron_ix_t inc = 0; inc < cnet->nnc; inc++)
      { fprintf(wr, "%sneuron_class = %d\n", ind2, inc);
        nmsim_class_neuron_write(wr, cnet->nclass[inc], timeStep);
      }
    for (nmsim_class_synapse_ix_t isc = 0; isc < cnet->nsc; isc++)
      { fprintf(wr, "%ssynapse_class = %d\n", ind2, isc);
        nmsim_class_synapse_write(wr, cnet->sclass[isc]);
      }
    fputs(ind1, wr);
    filefmt_write_footer(wr, nmsim_class_net_FILE_TYPE);

    fflush(wr);
  }

nmsim_class_net_t *nmsim_class_net_read(FILE *rd, double timeStep)
  {
    /* Read header line: */
    filefmt_read_header(rd, nmsim_class_net_FILE_TYPE, nmsim_class_net_VERSION);
    
    /* Read the number of neuron classes: */
    nmsim_class_neuron_count_t nnc = (nmsim_class_neuron_count_t)
      nmsim_read_int64_param(rd, "neuron_classes", 1, nmsim_class_neuron_count_MAX);

   /* Read the number of synapse classes: */
    nmsim_class_synapse_count_t nsc = (nmsim_class_synapse_count_t)
      nmsim_read_int64_param(rd, "synapse_classes", 0, nmsim_class_synapse_count_MAX);
 
    /* Create the network description: */
    nmsim_class_net_t *cnet = nmsim_class_net_new(nnc, nsc);
    
    /* Read the neuron classes: */
    for (nmsim_class_neuron_ix_t inc = 0; inc < nnc; inc++)
      { (void)nmsim_read_int64_param(rd, "neuron_class", inc, inc);
        nmsim_class_neuron_t *nclass = nmsim_class_neuron_read(rd, timeStep);
        cnet->nclass[inc] = nclass;
      }
    
    /* Read the synapse classes: */
    for (nmsim_class_synapse_ix_t isc = 0; isc < nsc; isc++)
      { (void)nmsim_read_int64_param(rd, "synapse_class", isc, isc);
        nmsim_class_synapse_t *sclass = nmsim_class_synapse_read(rd);
        cnet->sclass[isc] = sclass;
      }

    /* Read footer line: */
    filefmt_read_footer(rd, nmsim_class_net_FILE_TYPE);

    return cnet;
  }

void nmsim_class_net_compare(nmsim_class_net_t *cnet_read, nmsim_class_net_t *cnet_orig)
  { 
    /* Compare neuron classes: */
    nmsim_compare_int64_param("neuron_classes", cnet_read->nnc, cnet_orig->nnc);
    for (nmsim_class_neuron_ix_t inc = 0; inc < cnet_orig->nnc; inc++)
      { nmsim_class_neuron_compare(cnet_read->nclass[inc], cnet_orig->nclass[inc]); }

    /* Compare synapse classes: */
    nmsim_compare_int64_param("synapse_classes", cnet_read->nsc, cnet_orig->nsc);
    for (nmsim_class_synapse_ix_t isc = 0; isc < cnet_orig->nsc; isc++)
      {  nmsim_class_synapse_compare(cnet_read->sclass[isc], cnet_orig->sclass[isc]); }
  }    

