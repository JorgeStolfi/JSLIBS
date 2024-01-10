/* See {nmsim_elem_net_sim.h} */
/* Last edited on 2019-02-13 18:30:15 by jstolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsrandom.h>

#include <nmsim_basic.h>

#include <nmsim_class_net.h>
#include <nmsim_group_net.h>
#include <nmsim_elem_net.h>

#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>

#include <nmsim_firing_func.h>

#include <nmsim_elem_net_sim.h>

/* AUXILIARY PROCEDURES */

void nmsim_elem_net_determine_firings
  ( nmsim_elem_net_t *enet, 
    double V[], /* Membrane potential of each neuron (IN/OUT). */
    bool_t X[]  /* Firing indicators in time step (OUT). */
  );
  /* Determines which neurons fire during a time step.
    
    On input, the vector {V[0..nne-1]} must describe the membrane
    potential of each neuron at some discrete time {t}; where {nne} is
    {enet.nne}.
    
    The procedure decides, based on {V[i]}, whether each neuron will
    fire between times {t} and {t+1}; and sets the firing indicator
    {X[i]} to {TRUE} or {FALSE}, accordingly.  The vector {V} is not
    modified. */

void nmsim_elem_net_compute_tot_inputs
  ( nmsim_elem_net_t *enet, 
    nmsim_step_count_t age[], /* Firing age of each neuron (IN). */
    double I[],               /* External input of each neuron  (IN). */
    bool_t X[],               /* Firing indicators in time step (IN). */
    double H[],               /* Output modulator of each neuron at time {t} (OUT). */
    double J[]                /* Total input of each neuron (OUT). */
  );
  /* Computes the total input {J[i]} of each neuron {i} in {0..nne-1} 
    in the step from some time {t} to time {t+1}; where {nne}
    is {enet.nne}.
    
    On input, {I[i]} must be the the external input received by each
    neuron {i} in {0..nne-1} between time {t} and {t+1}. 
    Specifically, the voltage increment in that time step that is due to signals
    from outside the network, such as injected current or firings in axons from sensory
    neurons.
    
    On input, {X[i]} must be true iff the neuron {i} fired within that interval.
    
    On output, {H[i]} is the {J[i]} is {I[i]} plus the increment in the membrane
    potential between those two times that is due to the firings of
    afferent neurons in that time step. It does not include the
    effects of potential recharge and potential reset after firing. */

void nmsim_elem_net_update_states
  ( nmsim_elem_net_t *enet, 
    bool_t X[],               /* Firing indicators in time step (IN). */
    double J[],               /* Total input of each neuron (IN). */
    double V[],               /* Membrane potential of each neuron (IN/OUT). */
    nmsim_step_count_t age[], /* Firing age of each neuron (IN/OUT). */
    double M[]                /* Recharge modulator of each neuron (OUT). */
  );
  /* Updates the potential {V[i]} and firing age {age[i]} of each neuron {i}
    from time {t} to time {t+1}, according to the firing indicators
    {X[i]} and the total input {J[i]}. 
    
    Namely, if {X[i]} is true, the potential {V[i]} is set to the
    reset potential and {age[i]} is set to zero; otherwise, {V[i]} is
    decayed towards the resting potential, and {age[i]} is
    incremented. In either case, the total input {J[i]} is then added
    to {V[i]}. */

/* IMPLEMENTATIONS */

void nmsim_elem_net_sim_step
  ( nmsim_elem_net_t *enet,    /* Network description. */
    nmsim_time_t t,            /* Time at start of step. */
    double I[],                /* External neuron inputs (IN). */
    double V[],                /* Neuron potentials (IN/OUT). */
    nmsim_step_count_t age[],  /* Firing ages of neurons (IN/OUT). */
    bool_t X[],                /* Firing indicator of each neuron (OUT). */
    double H[],                /* Output modulator of each neuron (OUT). */
    double M[],                /* Recharge modulator of each neuron (OUT). */
    double J[],                /* Total input of each neuron (OUT). */
    nmsim_elem_net_trace_t *etrace /* Monitored neuron traces. */
  )
  {
    nmsim_elem_neuron_count_t nne = enet->nne; /* Number of neurons in network. */
    
    /* Determine which neurons fire: */
    nmsim_elem_net_determine_firings(enet, V, X);
    if (etrace != NULL) { nmsim_elem_net_trace_set_V_X(etrace, t, nne, V, X); }
    
    /* Compute {H} modulators and the inputs {I,J} of each neuron: */
    nmsim_elem_net_compute_tot_inputs(enet, age, I, X, H, J);
    if (etrace != NULL) { nmsim_elem_net_trace_set_age_H_I_J(etrace, t, nne, age, H, I, J); }

    /* Compute {M} modulators and update potentials and ages: */
    nmsim_elem_net_update_states(enet, X, J, V, age, M);
    if (etrace != NULL) { nmsim_elem_net_trace_set_M(etrace, t, nne, M); }
  }

/* INTERNAL IMPLEMENTATIONS */

void nmsim_elem_net_determine_firings
  ( nmsim_elem_net_t *enet, 
    double V[],      /* Membrane potential of each neuron (IN/OUT). */
    bool_t X[]       /* Firing indicators in time step (OUT). */
  )
  {
    nmsim_group_net_t *gnet = enet->gnet;
    nmsim_class_net_t *cnet = gnet->cnet;
    nmsim_elem_neuron_count_t nne = enet->nne;
    
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++)
      { 
        nmsim_elem_neuron_ix_t ing = enet->neu.e[ine].ing; /* Neuron's group ix. */
        nmsim_class_neuron_ix_t inc = gnet->ngrp.e[ing].inc; /* Neuron's class ix. */
        nmsim_class_neuron_t *nclass = cnet->nclass.e[inc]; /* Neuron's class. */
        nmsim_firing_func_t *Phi = &(nclass->Phi); /* Neuron's firing funtion. */
        double pr; /* Probability of neuron firing in time step. */
        nmsim_firing_func_eval(Phi, V[ine], &pr, NULL);
        /* Decide firing: */
        X[ine] = (drandom() < pr);
      }
  }

void nmsim_elem_net_compute_tot_inputs
  ( nmsim_elem_net_t *enet, 
    nmsim_step_count_t age[], /* Firing age of each neuron at time {t} (IN). */
    double I[],               /* External input of each neuron from {t} to {t+1} (IN). */
    bool_t X[],               /* Firing indicators from {t} to {t+1} (IN). */
    double H[],               /* Output modulator of each neuron at time {t} (OUT). */
    double J[]                /* Total input of each neuron from {t} to {t+1} (OUT). */
  )
  {
    nmsim_group_net_t *gnet = enet->gnet;
    nmsim_class_net_t *cnet = gnet->cnet;
    nmsim_elem_synapse_count_t nse = enet->nse;
    nmsim_elem_neuron_count_t nne = enet->nne;
    
    /* Initialize the total inputs vector, and compute the modulators {H[]}: */
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) 
      { /* The total input starts with the external input: */
        J[ine] = I[ine];
        /* Compute the output strength modulator: */
        nmsim_elem_neuron_ix_t ing = enet->neu.e[ine].ing; /* Neuron's group ix. */
        nmsim_class_neuron_ix_t inc = gnet->ngrp.e[ing].inc; /* Neuron's class ix. */
        nmsim_class_neuron_t *nclass = cnet->nclass.e[inc]; /* Neuron's class. */
        H[ine]= nmsim_class_neuron_compute_H(nclass, age[ine]);
      }

    for (nmsim_elem_synapse_ix_t kse = 0; kse < nse; kse++)
      { 
        /* Get the synapse {kse]}: */
        nmsim_elem_synapse_t *syn = &(enet->syn.e[kse]);
        /* Get indices of its neurons: */ 
        nmsim_elem_neuron_ix_t ine_pre = syn->ine_pre; /* Pre-synaptic neuron. */
        nmsim_elem_neuron_ix_t ine_pos = syn->ine_pos; /* Post-synaptic neuron. */
        if (X[ine_pre])
          { /* The pre-synaptic neuron fired; accumulate its pulse. */
            double dV = H[ine_pre] * syn->W;
            J[ine_pos] += dV;
          }
      } 
  }

void nmsim_elem_net_update_states
  ( nmsim_elem_net_t *enet, 
    bool_t X[],               /* Firing indicators in time step (IN). */
    double J[],               /* Total input of each neuron (IN). */
    double V[],               /* Membrane potential of each neuron (IN/OUT). */
    nmsim_step_count_t age[], /* Firing age of each neuron (IN/OUT). */
    double M[]                /* Recharge modulator of each neuron (OUT). */
  )
  {
    nmsim_group_net_t *gnet = enet->gnet;
    nmsim_class_net_t *cnet = gnet->cnet;
    nmsim_elem_neuron_count_t nne = enet->nne;
    
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++)
      { 
        nmsim_elem_neuron_ix_t ing = enet->neu.e[ine].ing; /* Neuron's group ix. */
        nmsim_class_neuron_ix_t inc = gnet->ngrp.e[ing].inc; /* Neuron's class ix. */
        nmsim_class_neuron_t *nclass = cnet->nclass.e[inc]; /* Neuron's class. */
        M[ine] = nmsim_class_neuron_compute_M(nclass, age[ine]);
        if (X[ine])
          { /* Neuron fired: */
            V[ine] = nclass->V_R;
            age[ine] = 0;
          }
        else
          { /* Neuron did not fire: */
            V[ine] = nmsim_class_neuron_recharge(nclass, V[ine], M[ine]);
            age[ine]++;
          }
      }
  }
