#ifndef nmsim_description_H
#define nmsim_description_H
 
/* Basic description of the Galves-Löcherbach neuron model. */
/* Last edited on 2020-12-06 19:13:22 by jstolfi */

/* 
  THE GL NEURON MODEL
  
  The Galves-Löcherbach (GL) neuron is a stochastic leaky
  integrate-and-fire model of a neuron's behavior. It is expected to be
  useful in the statistical simulation of large artificial neuronal
  networks with partially random topologies and synaptic weights,
  comprised of a relatively small number of homogeneous neuron
  populations.
  
  The state of a GL neuron is assumed to change at discrete /sampling
  times/, with uniform spacing {\Delta}. Each sampling time is
  identified by an integer {t} and represents the time {t\Delta},
  counted from the start of the simulation.  
  
  Between discrete times {t} and {t+1}, the neuron may /fire/ or not. If
  the neuron fires, its state returns to a fixed /reset state/ that is a
  time-invariant parameter of the model.  At the same time, pulse-like
  signals are sent to some set of /output neurons/.
  
  The time step {\Delta} must be large enough to allow us to assume that
  every spike happens entirely between to sampling times. On the other
  hand, {\Delta} must small enough to prevent a neuron from spiking
  twice in the same interval.
  
  The /firing indcator/ {X[t]} is the number of times (0 or 1) that
  neuron {i} fired between times {t} and {t+1}.
  
  The state of the neuron is its /membrane potential/ {V[t]} (a real
  number) and a non-negative /firing age/ {a[t]}. 
  
  In this library, the of membrane potential {V[t]} is taken to be its
  /resting potential/ {V_B}, when there has been no inputs and no firings for
  a long time.
  
  The firing age {a[t]} is the number of whole time steps that elapsed
  between its last firing and time {t}. That is, {a[t]} is {k} iff the
  neuron fired betwen times {t-k-1} and {t-k}, and did not fire between
  {t-k} and {t}. Thus, {a[t+1]} is zero if {X[t]} is 1, otherwise it
  is {a[t]+1}.  In some contexts, a negative value is used to indicate
  "unknown age".
  
  All other parameters of the neuron are determined by the potential and
  age, and therefore do not depend on what happened in the network
  before the time of the last firing. This is a defining property of the
  Galves-Löcherbach model.
  
  In the GL model, the firing of neuron {i} between times {t} and {t+1}
  is assumed to be dependent only on the potential {V[t]} through a
  fixed \emph{firing function} {\Phi}.
  
  Note that the state of a neuron depends only on its age {k} and on the
  inputs it received between sampling times {t-k} and {t}.
  
  If the neuron fails to fire between {t} and {t+1}, a part of its
  potential {V[t]} is ``leaked'', and only a fraction {c[t]*V[t]}
  of it is retained in the new potential {V[t+1]}.
  
  Anyway, from {t} to {t+1} the potential is incremented by a
  combination of /external inputs/ external to
  the network, and signals received from some set of /input neurons/
  through chemical synapses. 
  
  NEURON NETWORK
  
  A /Galves-Loecherbach neuron network/ is a finite set of {N} /neurons/,
  each conventionally identified by a /neuron index/ from 0 to {N-1}, 
  that are connected by /synapses/ and receive /external inputs/
  generated outside the network.
  
  Specifically, we denote by {I[i][t]} the membrane potential increment that is due to the
  external input signals received by neuron {i} between discrete times {t} and
  {t+1}.
  
  A synapse between a neuron {j} and a neuron {i} causes an
  increment in the potential of neuron {i} from time {t} to time {t+1} if, and 
  only if, neuron {j} fires in that time interval.   In that case, the increment,
  denoted by {w[j-->i][t]}, is called the /weight/ or /strength/ of the synapse for that
  time step.  By definition, it is the product of an invariant /resting weight/ or /resting strength/
  {w_B[j-->i]}, positive or negative, and an /output synapse modulator/ {H[j][t]} 
  that depends on the firing age of the respective neuron at time {t}.
  
  A synapse in the GL model thus models a chemical synapse in a biological neuronal
  network.  More precisely, it represents the path consisting of the pre-synaptic axon, a 
  chemical synaptic junction, and a path in the dendritic tree 
  leading to the post-synaptic neuron body.  
  
  The factor {H[t]} can be seen as modeling the temporary inhibition of
  chemical synapses of a neuron due to depletion of neurotransmitter
  vescicles near the synaptic junction.
  
  !!! Can we get Hebbian learning from this formulation? !!!
  
  !!! Removing {G[t]} because it seems unnecessary. 
  Multiplying {V[t]} does not make sense because it is
  relative to an arbitrary reference potential. Adding a bias to
  {V[t]} is equivalent to changing the reset potential and leak
  coefficient. !!!

  ???
    The firing indicator {X[t]} tells whether neuron {i} fired between
    times {t} and {t+1}. It is assumed to be a random variable, whose
    distribution depends only on {V[t]} and {G[t]}:

      { \Pr(X[t]=1 | V[t]) = Phi(G[t]*V[t]) }

    If the neuron {i} fires between times {t} and {t+1} (that is, if
    {X[t]=1}), the state {(V[t+1],c[t+1],G[t+1],H[t+1])} is reset to
    {(V_R,c_B,G_R,H_R)}.

    If the neuron does not fire between times {t} and {t+1} (that is, if
    {X[t]=0}), the state at {t+1} is

      { V[t+1] = V_B + c[t]*(V[t] - V_B) + I[t] + \sum_{j} (X[j][t]*H[j][t]*w[j-->i]) }

      { c[t+1] = c_B + c_mu*(c[t] - c_B) }

      { G[t+1] = 1 - G_mu*(1 - G[t]) }

      { H[t+1] = 1 - H_mu*(1 - H[t]) }  
*/ 
 
#endif
