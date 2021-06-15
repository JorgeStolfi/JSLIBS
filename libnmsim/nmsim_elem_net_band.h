#ifndef nmsim_elem_net_band_H
#define nmsim_elem_net_band_H
 
/* Generates a simple multilayer all-to-all network of GL neurons, open or closed. */
/* Last edited on 2021-01-09 16:41:21 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>

#include <nmsim_basic.h>

#include <nmsim_class_net.h>
#include <nmsim_group_net.h>
#include <nmsim_elem_net.h>
  
nmsim_elem_net_t *nmsim_elem_net_band_make
  ( nmsim_class_net_t *cnet,
    int32_t numBands,
    nmsim_group_neuron_count_t numLayers, 
    nmsim_elem_neuron_count_t bandWidth,
    bool_t closed,
    double WMin,
    double WMax,
    double WDev
  );
  /* Creates an a simple element-level networks, suitable for testing the effect 
    of multiple parameters like synapse strength or input current. 
    All neurons and synapses are assumed to have the same class, respectively {cnet->nclass[0]} and
    {cnet->sclass[0]}.  See {nmsim_elem_net_band_make_INFO} for explanation of the parameters. */
  
#define nmsim_elem_net_band_make_INFO \
  "The basic network has {n = numLayers} layers with {m = bandWidth} neurons" \
  " in each layer.  Each neuron in layer {k+1} receives inputs from all" \
  " the {m} neurons of in layer {k}.  If {closed} is true, there each" \
  " neuron in the first layer receives inputs from all {m} neurons of" \
  " the last layer.\n" \
  "\n" \
  "  The total weight of all input synapses of each neuron is the" \
  " same, a parameter {W}.\n" \
  "\n" \
  "  The network actually consists of {numBands} independent copies" \
  " (/bands/) of this basic network, provided for the sake of simultaneous" \
  " testing with different values of the parameters or external" \
  " inputs.\n" \
  "\n" \
  "  The weights of the synapses into each neuron are drawn from a log-normal" \
  "  distribution with mean {W_avg} and deviation {W_dev}; where {W_avg} is" \
  "  the same for every neuron in each band" \
  " and varies from band to band linearly from {WMin/bandWidth} to {WMax/bandWidth}, and" \
  "  {W_dev} is {WDev*W_avg}. If {numBands} is 1, {W_avg} will be" \
  " {(WMin+WMiax}/2/bandWidth}."

#endif

