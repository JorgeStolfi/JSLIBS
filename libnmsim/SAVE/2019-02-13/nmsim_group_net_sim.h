#ifndef nmsim_group_net_sim_H
#define nmsim_group_net_sim_H
  
/* Group-level simulation ofa networks of Galves-Löcherbach neurons. */
/* Last edited on 2019-02-13 17:42:38 by jstolfi */

#define _GNU_SOURCE

#include <nmsim_basic.h>
#include <nmsim_group_net.h>
#include <nmsim_group_neuron_state.h>
#include <nmsim_group_net_trace.h>

void nmsim_group_net_sim_step
  ( nmsim_group_net_t *gnet,             /* Network description. */
    nmsim_time_t t,                      /* Time at start of step. */
    double Iavg[],                       /* Avg of external neuron input per population (IN). */
    double Idev[],                       /* Dev of external neuron input per population (IN). */
    nmsim_group_neuron_state_t state[],  /* Statistical state of each population (IN/OUT). */
    nmsim_group_net_trace_t *gtrace      /* Monitored neuron group traces. */
  );
  /* Simulates the evolution of a network {gnet} of {nng=gnet.nng}
    groups (populations) of GL neurons, during the time step from
    discrete time {t} to time {t+1}.
    
    On input, the element {st=state[i]}, for {i} in {0..nng-1}, must
    describe the distribution of the states of neurons in each group {i} at time {t}. 
    Specifically, {st.cs[s].V_avg} and {st.cs[s].V_dev} should contain the mean and deviation of
    the voltage distribution of the neurons in the cohort with age {s},
    and {st.cs[s].eta} must be the and fraction of the neuron group that is in that cohort.

    In addition to the synaptic inputs, each neuron in a group {i} is also
    assumed to have received, in that time step, an "external" input
    current that, by itself, would raise its potential by a certain
    amount.  That amount is a Gaussian random variable with mean {Iavg[i]}
    and deviation {Idev[i]}.  If {Iavg} is {NULL}, the external input is
    assumed to be zero for all neurons.  If {Idev} is {NULL}, the deviation is assumed
    to be zero for all neuron groups.
    
    The procedure computes the fields {st.cs[s].rho}, {st.cs[s].M}, and {st.cs[s].H} for each cohort.
    It also saves {Iavg[i]}, {Javg[i]} in {st.Iavg,st.Idev}, and computes {st.rho} (the overall 
    firing rate for the group).
    
    With those values, it updates the cohort states {st.cs[s].V_avg}, {st.cs[s].V_dev},
    and {st.cs[s].eta} for time time {t+1}.
    
    The procedure sets the field {state[i].Xavg[i]} to the estimated fraction of the neurons in 
    group {i} that are supposed to fire in the interval from {t} to {t+1}.
    
    The procedure also sets {state[i].Javg} and {state[i].Jdev} to the estimated mean and deviation 
    of the total input (synaptic inputs plus external input) for the neurons in group {i},
    between {t} and {t+1}. 
    
    The input contents of the fields {.Xavg}, {.Javg}, and {.Jdev} are ignored.
    
    If {gtrace} is not {NULL}, the procedure also stores into it 
    the states of the selected neuron groups. Specifically, for each monitored
    neuron group whose trace {trng} is in {gtrace.trng[0..gtrace.nng-1]} such 
    that {t} is in {trng.tini .. trng.tfin}, saves in the entry of {trng.st} corresponding to
    time {t} a record with a suitable summary of that state.  Note that 
    the entries of {gtrace} corresponding to time {t+1} are NOT set. */

#endif
