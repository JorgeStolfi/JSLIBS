/* See {nmsim_elem_net_sim.h} */
/* Last edited on 2020-12-10 01:46:33 by jstolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <affirm.h>
#include <bool.h>
#include <jsrandom.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>

#include <nmsim_class_net.h>
#include <nmsim_group_net.h>
#include <nmsim_elem_net.h>

#include <nmsim_class_neuron.h>
#include <nmsim_group_neuron.h>
#include <nmsim_elem_neuron.h>

#include <nmsim_firing_func.h>

#include <nmsim_elem_net_sim.h>

/* AUXILIARY PROCEDURES */

void nmsim_elem_net_sim_determine_firings
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

void nmsim_elem_net_sim_compute_tot_inputs
  ( nmsim_elem_net_t *enet, 
    bool_t X[],               /* Firing indicators in time step (IN). */
    double H[],               /* Output modulator of each neuron at time {t} (IN). */
    double I[],               /* External input of each neuron (IN, mV). */
    double J[]                /* Total input of each neuron (OUT, mV). */
  );
  /* Computes the total input {J[i]} (mV) of each neuron {i} in {0..nne-1} 
    in the step from some time {t} to time {t+1}; where {nne}
    is {enet.nne}.
    
    On input, {X[i]} must be true iff the neuron {i} fired within that interval.
    
    On input, {H[i]} must be the output synapse strength modulator
    of neuron {i}.
    
    On input, {I[i]} must be the the external input received by each
    neuron {i} in {0..nne-1} between time {t} and {t+1}. Specifically,
    it must be the increment in the membrane potential of neuron {i}
    between those two times that is due to signals received from outside
    the network during that time step, such as itegrated injected
    current or the firings in axons from sensory neurons.
    
    On output, {J[i]} will be {I[i]} plus the increment in the
    membrane potential of neuron {i} between those two times that is due to the
    firings {X[j]} of afferent neurons during that time step, weighted by
    their synaptic weights {W[j-->i]} and respective modulators
    {H[j]}. The value of {J[i]} will not include the effects of
    potential recharge and potential reset after firing. */

void nmsim_elem_net_sim_update_states
  ( nmsim_elem_net_t *enet, 
    double V[],               /* Membrane potential of each neuron (IN/OUT,mV). */
    nmsim_step_count_t age[], /* Firing age of each neuron (IN/OUT). */
    double M[],               /* Recharge modulator of each neuron (IN/OUT). */
    double H[],               /* Synapse strength modulator of each neuron (IN/OUT). */
    bool_t X[],               /* Firing indicators in time step (IN). */
    double J[]                /* Total input of each neuron (IN,mV). */
  );
  /* Updates the potential {V[i]} and firing age {age[i]} of each neuron {i}
    from time {t} to time {t+1}, according to the firing indicators
    {X[i]} and the total input {J[i]} for that interval. 
    
    Namely, if {X[i]} is true, the potential {V[i]} is set to the
    reset potential, ignoring {M[i]}, and {age[i]} is set to zero; otherwise, {V[i]} is
    decayed towards the resting potential with recharge factor multiplied by {M[i]}, 
    and {age[i]} is incremented. In either case, the total input {J[i]} is then added
    to {V[i]}.  The modulators {M[i]} and {H[i]} are then
    updated to reflect the new firing age. */

/* IMPLEMENTATIONS */

void nmsim_elem_net_sim_step
  ( nmsim_elem_net_t *enet,    /* Network description. */
    nmsim_time_t t,            /* Time at start of step. */
    double V[],                /* Neuron potentials (IN/OUT,mV). */
    nmsim_step_count_t age[],  /* Firing ages of neurons (IN/OUT). */
    double M[],                /* Recharge modulator of each neuron (IN/OUT). */
    double H[],                /* Output modulator of each neuron (IN/OUT). */
    bool_t X[],                /* Firing indicator of each neuron (OUT). */
    double I[],                /* External neuron inputs (IN,mV). */
    double J[],                /* Total input of each neuron (OUT,mV). */
    nmsim_elem_net_trace_t *etrace /* Monitored neuron traces. */
  )
  {
    nmsim_elem_neuron_count_t nne = enet->nne; /* Number of neurons in network. */
    
    /* Save the trace data for time {t}: */
    if (etrace != NULL) { nmsim_elem_net_trace_set_V_age_M_H(etrace, t, nne, V, age, M, H); }
    
    /* Compute the firing indicators {X} and total inputs {J} for step {t} to {t+1}: */
    nmsim_elem_net_sim_determine_firings(enet, V, X);
    nmsim_elem_net_sim_compute_tot_inputs(enet, X, H, I, J);
    
    /* Save the tace data for the step {t} to {t+1}: */
    if (etrace != NULL) { nmsim_elem_net_trace_set_X_I_J(etrace, t, nne, X, I, J); }

    /* Update potentials, ages, and modulators for time {t+1}: */
    nmsim_elem_net_sim_update_states(enet, V, age, M, H, X, J);
  }

void nmsim_elem_net_sim_compute_modulators
  ( nmsim_elem_net_t *enet, 
    nmsim_step_count_t age[], /* Firing age of each neuron (IN). */
    double M[],               /* Recharge modulator of each neuron (OUT). */
    double H[]               /* Output modulator of each neuron at time {t} (OUT). */
  )
  {
    nmsim_group_net_t *gnet = enet->gnet;
    nmsim_class_net_t *cnet = gnet->cnet;
    nmsim_elem_neuron_count_t nne = enet->nne;
    
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) 
      { nmsim_elem_neuron_ix_t ing = enet->neu.e[ine].ing; /* Neuron's group ix. */
        nmsim_class_neuron_ix_t inc = gnet->ngrp.e[ing].inc; /* Neuron's class ix. */
        nmsim_class_neuron_t *nclass = cnet->nclass.e[inc]; /* Neuron's class. */
        M[ine] = nmsim_class_neuron_compute_M(nclass, age[ine]);
        H[ine] = nmsim_class_neuron_compute_H(nclass, age[ine]);
      }
  }

/* INTERNAL IMPLEMENTATIONS */

void nmsim_elem_net_sim_determine_firings
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
        X[ine] = (pr <= 0.0 ? FALSE : (pr >= 1.0 ? TRUE : (drandom() < pr)));
      }
  }

void nmsim_elem_net_sim_compute_tot_inputs
  ( nmsim_elem_net_t *enet, 
    bool_t X[],               /* Firing indicators from {t} to {t+1} (IN). */
    double H[],               /* Output modulator of each neuron at time {t} (IN). */
    double I[],               /* External input of each neuron from {t} to {t+1} (IN). */
    double J[]                /* Total input of each neuron from {t} to {t+1} (OUT). */
  )
  {
    nmsim_elem_synapse_count_t nse = enet->nse;
    nmsim_elem_neuron_count_t nne = enet->nne;
    
    /* Initialize the total inputs vector: */
    for (nmsim_elem_neuron_ix_t ine = 0; ine < nne; ine++) 
      { /* The total input starts with the external input: */
        J[ine] = I[ine];
      }

    for (nmsim_elem_synapse_ix_t kse = 0; kse < nse; kse++)
      { 
        /* Get the synapse {kse]}: */
        nmsim_elem_synapse_t *syn = &(enet->syn.e[kse]);
        /* Get indices of its neurons: */ 
        nmsim_elem_neuron_ix_t ine_pre = syn->ine_pre; /* Pre-synaptic neuron. */
        if (X[ine_pre])
          { /* The pre-synaptic neuron fired; accumulate its pulse. */
            nmsim_elem_neuron_ix_t ine_pos = syn->ine_pos; /* Post-synaptic neuron. */
            double dV = H[ine_pre] * syn->W;
            J[ine_pos] += dV;
          }
      } 
  }

void nmsim_elem_net_sim_update_states
  ( nmsim_elem_net_t *enet, 
    double V[],               /* Membrane potential of each neuron (IN/OUT). */
    nmsim_step_count_t age[], /* Firing age of each neuron (IN/OUT). */
    double M[],               /* Recharge modulator of each neuron (IN/OUT). */
    double H[],               /* Synapse strength modulator of each neuron (IN/OUT). */
    bool_t X[],               /* Firing indicators in time step (IN). */
    double J[]                /* Total input of each neuron (IN). */
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
        if (X[ine])
          { /* Neuron fired: */
            V[ine] = nclass->V_R;
            age[ine] = 0;
            M[ine] = nclass->M_R;
            H[ine] = nclass->H_R;
          }
        else
          { /* Neuron did not fire: */
            demand(! isnan(V[ine]), "invalid potential");
            demand(! isnan(M[ine]), "invalid recharge modulator");
            demand(! isnan(H[ine]), "invalid output modulator");
            demand(age[ine] >= 0, "invalid age");
            V[ine] = nmsim_class_neuron_recharge(nclass, V[ine], M[ine]);
            age[ine]++;
            M[ine] = 1 - (1 - M[ine])*nclass->M_mu;
            H[ine] = 1 - (1 - H[ine])*nclass->H_mu;
          }
        /* Add input: */
        V[ine] += J[ine];
      }
  }
