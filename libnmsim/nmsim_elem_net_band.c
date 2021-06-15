/* See {nmsim_elem_net_band.h} */
/* Last edited on 2021-01-09 20:01:30 by jstolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsrandom.h>

#include <nmsim_class_neuron.h>
#include <nmsim_class_synapse.h>
#include <nmsim_group_neuron.h>
#include <nmsim_group_synapse.h>
#include <nmsim_elem_neuron.h>
#include <nmsim_elem_synapse.h>
#include <nmsim_elem_net.h>

#include <nmsim_elem_net_band.h>

void nmsim_elem_net_band_add_one
  ( nmsim_class_net_t *cnet,
    nmsim_group_net_t *gnet,
    nmsim_elem_net_t *enet,
    int32_t iBand,
    int32_t numBands,
    nmsim_group_neuron_count_t numLayers, 
    nmsim_elem_neuron_count_t bandWidth,
    bool_t closed,
    double W_avg,
    double W_dev,
    nmsim_group_neuron_ix_t *ingNextP,
    nmsim_group_synapse_ix_t *isgNextP,
    nmsim_elem_neuron_ix_t *ineNextP,
    nmsim_elem_synapse_ix_t *iseNextP
  );
  /* Adds the groups and elements of band {iBand} to the networks {gnet} and {enet}.
    All synapses will have weight {W}.
    On entry, {*ingNextP,*isgNextP} should be the indices of the 
    next free slot in the neuron group and synapse group tables
    of {gnet}, and {*ineNextP,*iseNextP} shoul be the next free entries
    in the neuron elem and synapse elem tables of {enet}.  These
    variables are updated on exit. */

nmsim_group_neuron_t nmsim_elem_net_band_make_group_neuron
  ( nmsim_class_neuron_ix_t inc,          /* Class of neurons in population. */
    nmsim_elem_neuron_count_t nne,        /* Number of neurons in group. */
    nmsim_elem_neuron_ix_t ine_start,     /* Index of first neuron in group, or 0. */
    nmsim_group_synapse_count_t nsg_out   /* Number of synaptic bundles out of this group. */
  );
  /* Creates a neuron group with the given parameters and default values for the
    others. */

nmsim_group_synapse_t nmsim_elem_net_band_make_group_synapse
  ( nmsim_class_synapse_ix_t isc,      /* Class of synapses in synapse group. */
    nmsim_group_neuron_ix_t ing_pre,   /* Index of pre-synaptic neuron group. */
    nmsim_group_neuron_ix_t ing_pos,   /* Index of post-synaptic neuron group. */
    nmsim_elem_synapse_count_t nse  
  );
  /* Creates a synapse group with the given parameters and default values for the
    others. */
    
nmsim_elem_neuron_t nmsim_elem_net_band_make_elem_neuron
  ( nmsim_group_neuron_ix_t ing,            /* Index of neuron group. */
    nmsim_elem_synapse_ix_t nse_out,        /* Number of output synapses (may be zero). */
    nmsim_elem_synapse_ix_t ise_out_start   /* Index of first output synapse, or 0. */
  );
  /* Creates a neuron elem with the given parameters and default values for the
    others. */

nmsim_elem_synapse_t nmsim_elem_net_band_make_elem_synapse
  ( nmsim_group_synapse_ix_t isg,      /* Index of synaptic group. */
    nmsim_elem_neuron_ix_t ine_pre,    /* Index of pre-synaptic neuron. */
    nmsim_elem_neuron_ix_t ine_pos,    /* Index of post-synaptic neuron. */
    float W                            /* Resting weight (mV). A {float} to save mem. */
  );
  /* Creates a synapse elem with the given parameters and default values for the
    others. */

/* IMPLEMENTATIONS */

nmsim_elem_net_t *nmsim_elem_net_band_make
  ( nmsim_class_net_t *cnet,
    int32_t numBands,
    nmsim_group_neuron_count_t numLayers, 
    nmsim_elem_neuron_count_t bandWidth,
    bool_t closed,
    double WMin,
    double WMax,
    double WDev
  )
  { 
    /* Create the group_level network: */
    nmsim_group_neuron_count_t nng_per_band = numLayers;      /* Neuron groups per band. */
    nmsim_group_neuron_count_t nng = numBands*nng_per_band;  /* Neuron groups total. */
    
    nmsim_group_synapse_count_t nsg_per_band = (numLayers + closed - 1); /* Synapse bundles per band. */
    nmsim_group_synapse_count_t nsg = numBands*nsg_per_band;  /* Synapse bundles total. */

    nmsim_group_net_t *gnet = nmsim_group_net_new(cnet, nng, nsg);

    nmsim_elem_neuron_count_t nne_per_ng = bandWidth;           /* Neurons per neuron group. */
    nmsim_elem_neuron_count_t nse_per_sg = bandWidth*bandWidth; /* Synapses per synapse bundle. */
    
    /* Create the elem level network: */
    nmsim_elem_neuron_count_t  nne = nng*nne_per_ng;    /* Number of neuron elems. */
    nmsim_elem_synapse_count_t nse = nsg*nse_per_sg;    /* Number of synapse elems. */

    nmsim_elem_net_t *enet = nmsim_elem_net_new(gnet, nne, nse);

    /* Add the bands: */
    nmsim_group_neuron_ix_t ingNext = 0;
    nmsim_group_synapse_ix_t isgNext = 0;
    nmsim_elem_neuron_ix_t ineNext = 0;
    nmsim_elem_synapse_ix_t iseNext = 0;
    for (int32_t iBand = 0; iBand < numBands; iBand++)
      { double W; /* Total input synapse weight for each neuron. */
        if (numBands == 1)
          { W = (WMin + WMax)/2; }
        else
          { double r = ((double)iBand)/((double)numBands - 1); /* Ranges in {[0 _ 1]}. */
            W = (1-r)*WMin + r*WMax;
          }
        double W_avg = W/bandWidth;
        double W_dev = WDev*W_avg; 
        nmsim_elem_net_band_add_one
          ( cnet, gnet, enet,
            iBand, numBands,
            numLayers, bandWidth, closed,
            W_avg, W_dev,
            &ingNext, &isgNext, &ineNext, &iseNext
          );
      }
    
    assert(ingNext == nng);
    assert(isgNext == nsg);
    assert(ineNext == nne);
    assert(iseNext == nse);

    return enet;

  }

void nmsim_elem_net_band_add_one
  ( nmsim_class_net_t *cnet,
    nmsim_group_net_t *gnet,
    nmsim_elem_net_t *enet,
    int32_t iBand,
    int32_t numBands,
    nmsim_group_neuron_count_t numLayers, 
    nmsim_elem_neuron_count_t bandWidth,
    bool_t closed,
    double W_avg,
    double W_dev,
    nmsim_group_neuron_ix_t  *ingNextP,
    nmsim_group_synapse_ix_t *isgNextP,
    nmsim_elem_neuron_ix_t   *ineNextP,
    nmsim_elem_synapse_ix_t  *iseNextP
  )
  { 
    nmsim_group_neuron_count_t nng_per_band = numLayers;      /* Neuron groups per band. */

    nmsim_elem_neuron_count_t nne_per_ng = bandWidth;           /* Neurons per neuron group. */
    nmsim_elem_neuron_count_t nse_per_sg = bandWidth*bandWidth; /* Synapses per synapse bundle. */

    nmsim_group_neuron_ix_t   ingNext = (*ingNextP);
    nmsim_group_synapse_ix_t  isgNext = (*isgNextP);
    nmsim_elem_neuron_ix_t    ineNext = (*ineNextP);
    nmsim_elem_synapse_ix_t   iseNext = (*iseNextP);

    /* Record inital indices in band: */
    nmsim_group_neuron_ix_t   ing_start_band = ingNext;
    nmsim_group_neuron_ix_t   ine_start_band = ineNext;
    nmsim_group_synapse_ix_t  isg_start_band = isgNext;
    
    fprintf(stderr, "band %d  W synapse = %.4f ± %.4f", iBand, W_avg, W_dev);
    fprintf(stderr, "  W input total = %.4f ± %.4f\n", W_avg*bandWidth, W_dev*sqrt(bandWidth));
    
    /* Loop over the neuron groups of {iBand}: */
    for (nmsim_group_neuron_ix_t ing_in_band = 0; ing_in_band < nng_per_band; ing_in_band++)
      { /* Define the neuron group {ing_in_band} in band {iBand}: */
        nmsim_group_neuron_ix_t ing = ing_in_band + ing_start_band;
        nmsim_class_neuron_ix_t inc = 0;
        nmsim_elem_neuron_ix_t ine_start = ineNext;
        nmsim_group_synapse_count_t nsg_out = 1;
        gnet->ngrp[ingNext] = nmsim_elem_net_band_make_group_neuron(inc, nne_per_ng, ine_start, nsg_out);
        ingNext++;
        
        bool_t has_out = ((ing_in_band < nng_per_band-1) || closed); /* Neuron group has outputs? */
        if (has_out)
          { /* nmsim_group_synapse_ix_t isg_in_band = ing_in_band; */
            /* Define the synapse group {isg_in_band} in band {iBand}: */
            nmsim_class_synapse_ix_t isc = 0;
            nmsim_group_neuron_ix_t ing_pre_in_band = ing_in_band;
            nmsim_group_neuron_ix_t ing_pos_in_band = (ing_pre_in_band + 1) % nng_per_band;
            nmsim_group_neuron_ix_t ing_pre = ing_pre_in_band + ing_start_band;
            nmsim_group_neuron_ix_t ing_pos = ing_pos_in_band + ing_start_band;
            gnet->sgrp[isgNext] = nmsim_elem_net_band_make_group_synapse(isc, ing_pre, ing_pos, nse_per_sg);
            isgNext++;
          }
          
        /* Loop over the neurons of group {ing_in_band} in band {iBand}: */
        nmsim_elem_synapse_count_t nse_out = (has_out ? nne_per_ng : 0);
        for (nmsim_elem_neuron_ix_t ine_in_ng = 0; ine_in_ng < nne_per_ng; ine_in_ng++)
          { nmsim_elem_neuron_ix_t ine = ineNext;
            nmsim_elem_synapse_ix_t ise_out_start = iseNext;
            enet->neu[ineNext] = nmsim_elem_net_band_make_elem_neuron(ing, nse_out, ise_out_start);
            ineNext++;
            if (has_out)
              { /* Define the synapses out of this neuron: */
                nmsim_group_synapse_ix_t isg_in_band = ing_in_band;
                nmsim_group_synapse_ix_t isg = isg_in_band + isg_start_band;
                nmsim_elem_neuron_ix_t ine_pre = ine;
                /* Compute the index {ine_pos} of the first post-synaptic neuron of {ine_pre}: */
                nmsim_group_neuron_ix_t ing_pos_in_band = (isg_in_band + 1) % nng_per_band;
                nmsim_elem_neuron_ix_t ine_pos_in_band = ing_pos_in_band * nne_per_ng;
                nmsim_elem_neuron_ix_t ine_pos = ine_start_band + ine_pos_in_band;
                for (nmsim_elem_synapse_ix_t ise_neu = 0; ise_neu < nne_per_ng; ise_neu++)
                  { double W = dloggaussrand(W_avg, W_dev);
                    enet->syn[iseNext] = 
                      nmsim_elem_net_band_make_elem_synapse(isg, ine_pre, ine_pos, (float)W);
                    iseNext++;
                    ine_pos++;
                  }
              }
          }
      }

    (*ingNextP) = ingNext;
    (*isgNextP) = isgNext;
    (*ineNextP) = ineNext;
    (*iseNextP) = iseNext;
  }
 
nmsim_group_neuron_t nmsim_elem_net_band_make_group_neuron
  ( nmsim_class_neuron_ix_t inc,          /* Class of neurons in population. */
    nmsim_elem_neuron_count_t nne,        /* Number of neurons in group. */
    nmsim_elem_neuron_ix_t ine_start,     /* Index of first neuron in group, or 0. */
    nmsim_group_synapse_count_t nsg_out   /* Number of synaptic bundles out of this group. */
  )
  {
    nmsim_group_neuron_t ngrp = 
      (nmsim_group_neuron_t){ .inc = inc, .nne = nne, .ine_start = ine_start, .nsg_out = nsg_out };
    return ngrp;
  }

nmsim_group_synapse_t nmsim_elem_net_band_make_group_synapse
  ( nmsim_class_synapse_ix_t isc,      /* Class of synapses in synapse group. */
    nmsim_group_neuron_ix_t ing_pre,   /* Index of pre-synaptic neuron group. */
    nmsim_group_neuron_ix_t ing_pos,   /* Index of post-synaptic neuron group. */
    nmsim_elem_synapse_count_t nse  
  )
  { nmsim_group_synapse_t sgrp = 
      (nmsim_group_synapse_t){ .isc = isc, .ing_pre = ing_pre, .ing_pos = ing_pos, .nse = nse };
    return sgrp;
  }
nmsim_elem_neuron_t nmsim_elem_net_band_make_elem_neuron
  ( nmsim_group_neuron_ix_t ing,            /* Index of neuron group. */
    nmsim_elem_synapse_ix_t nse_out,        /* Number of output synapses (may be zero). */
    nmsim_elem_synapse_ix_t ise_out_start   /* Index of first output synapse, or 0. */
  )
  { nmsim_elem_neuron_t neu = 
      (nmsim_elem_neuron_t){ .ing = ing, .nse_out = nse_out, .ise_out_start = ise_out_start };
    return neu;
  }

nmsim_elem_synapse_t nmsim_elem_net_band_make_elem_synapse
  ( nmsim_group_synapse_ix_t isg,      /* Index of synaptic group. */
    nmsim_elem_neuron_ix_t ine_pre,    /* Index of pre-synaptic neuron. */
    nmsim_elem_neuron_ix_t ine_pos,    /* Index of post-synaptic neuron. */
    float W                            /* Resting weight (mV). A {float} to save mem. */
  )
  {
    nmsim_elem_synapse_t syn = 
      (nmsim_elem_synapse_t){ .isg = isg, .ine_pre = ine_pre, .ine_pos = ine_pos, .W = W };
    return syn;
  }

