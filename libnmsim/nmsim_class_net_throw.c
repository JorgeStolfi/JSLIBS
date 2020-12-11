/* See {nmsim_class_net_throw.h} */
/* Last edited on 2020-12-11 13:38:49 by jstolfi */

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

#include <nmsim_class_net_throw.h>

nmsim_class_net_t *nmsim_class_net_throw
  ( nmsim_class_neuron_count_t nnc, 
    nmsim_class_synapse_count_t nsc
  )
  {
    demand((nnc >= 1) && (nnc <= nmsim_class_neuron_count_MAX), "invalid neuron class count");
    demand((nsc >= 0) && (nsc <= nmsim_class_synapse_count_MAX), "invalid synapse class count");
    
    nmsim_class_net_t *cnet = nmsim_class_net_new(nnc, nsc);
    
    /* Add neuron classes: */
    for (nmsim_class_neuron_ix_t inc = 0; inc < nnc; inc++)
      { nmsim_class_neuron_t *nclass = nmsim_class_neuron_throw();
        cnet->nclass[inc] = nclass;
      }
    
    /* Add synapse classes: */
    for (nmsim_class_synapse_ix_t isc = 0; isc < nsc; isc++)
      { nmsim_class_synapse_t *sclass = nmsim_class_synapse_throw();
        cnet->sclass[isc] = sclass;
      }
    return cnet;
  }
