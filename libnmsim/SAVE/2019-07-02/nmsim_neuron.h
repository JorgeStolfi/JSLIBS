#ifndef nmsim_neuron_H
#define nmsim_neuron_H
 
/* Types and functions for the Galves-Löcherbach neuron model. */
/* Last edited on 2016-12-08 15:54:41 by jstolfi */

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
  
  The /firing indcator/ {X[i,t]} is the number of times (0 or 1) that
  neuron {i} fired between times {t} and {t+1}.
  
  The state of the neuron is its /membrane potential/ {V[i,t]} (a real
  number) and a /firing age/ {a[i,t]}. 
  
  In this library, the of membrane potential {V[i,t]} is taken to be its
  /resting potential/ {V_B}, when there has been no inputs and no firings for
  a long time.
  
  The firing age {a[i,t]} is the number of whole time steps that elapsed
  between its last firing and time {t}. That is, {a[i,t]} is {k} iff the
  neuron fired betwen times {t-k-1} and {t-k}, and did not fire between
  {t-k} and {t}. Thus, {a[i,t+1]} is zero if {X[i,t]} is 1, otherwise it
  is {a[i,t]+1}.
  
  All other parameters of the neuron are determined by the potential and
  age, and therefore do not depend on what happened in the network
  before the time of the last firing. This is a defining property of the
  Galves-Löcherbach model.
  
  If the neuron fails to fire between {t} and {t+1}, a part of its
  potential {V[i,t]} is ``leaked'', and only a fraction {c[i,t]*V[i,t]}
  of it is retained in the new potential {V[i,t+1]}.
  
  Anyway, from {t} to {t+1} the potential is incremented by a
  combination of signals received from some set of /input neurons/
  through chemical synapses and gap junctions, and stimuli external to
  the network. 
  
  The external inputs are represented by a potential increment {I[i,t]}
  due the external stimuli received by neuron {i} between times {t} and
  {t+1}.
  
  A (chemical) synapse between a neuron {j} and a neuron {i} is
  described by a fixed /resting weight/ or /resting strength/
  {w_B[j-->i]}, positive or negative. The actual strength {w[j-->i,t]}
  is {H[j,t]*w_B[j-->i]*G[i,t]} where {H[j,t]} and {G[i,t]} are
  modulating factors that depend on the age of the respective neuron.
  
  The factor {H[j,t]} modulates the strength of the output synapses of
  neuron {j}. It can be seen as modeling the temporary inhibition of
  synapses due to depletion of neurotransmitter vescicles near the
  synaptic junction. The factor {G[i,t]} modulates the strength of
  all input synapses of neuron {i}.
  
  !!! Can we get Hebbian learning from this formulation? !!!
  
  !!! Taking again {G[i,t]} out of the {\Phi} function because it is
  cleaner. Multiplying {V[i,t]} does not make sense because it is
  relative to an arbitrary reference potential. Adding a bias to
  {V[i,t]} is equivalent to changing the reset potential and leak
  coefficient. Maybe multiply {G[i,t]} by the probability? !!!
  
  Similarly, a gap junction between a neuron {j} and a neuron {i} is
  described by a fixed non-negative /resting transmittance/
  {r_B[j-->i]}. The actual transmittance is {F[j,t]*r_B[j-->i]*E[i,t]}
  where {F[j,t]} and {E[i,t]} are modulating factors that depend on the
  firing ages of the two neurons.
  
  !!! Gap juntions. !!!
  
  !!! Modulation of leakage coefficient. !!!
  
  In the GL model, the firing of neuron {i} between times {t} and {t+1}
  is assumed to be dependent only on the potential {V[i,t]} through a
  fixed \emph{firing function} {\Phi}.
  
  Note that the state of a neuron depends only on its age {k} and on the
  inputs it received between sampling times {t-k} and {t}.
  
*/ 

#define _GNU_SOURCE
#include <stdio.h>

#include <nmsim_firing_func.h>

typedef struct nmsim_neuron_parms_t
  { /* Static neuron parameters: */
    double V_R;   /* Potential after firing. */
    double c_B;   /* Resting leakage factor. */
    struct nmsim_firing_func_t *Phi;  /* Firing function. */
    /* Leakage plasticity parameters: */
    double M_R;   /* Reset value of leakage modulator. */
    double M_mu;  /* Recovery factor of leakage modulator. */
    /* Chemical synapse platicity parameters: */
    double G_R;   /* Reset value of input synapse modulator. */
    double G_mu;  /* Recovery factor of output synapse modulator. */
    double H_R;   /* Reset value of input synapse modulator. */
    double H_mu;  /* Recovery factor of output synapse modulator. */
    /* Parameters of output gain dynamics: */
    /* Parameters of firing function:
  } nmsim_neuron_parms_t;
  /* Parameters of a GL neuron.
  
    The firing indicator {X[i,t]} tells whether neuron {i} fired between
    times {t} and {t+1}. It is assumed to be a random variable, whose
    distribution depends only on {V[i,t]} and {G[i,t]}:

      { \Pr(X[i,t]=1 | V[i,t]) = Phi(G[i,t]*V[i,t]) }

    If the neuron {i} fires between times {t} and {t+1} (that is, if
    {X[i,t]=1}), the state {(V[i,t+1],G[i,t+1],H[i,t+1])} is reset to
    {(V_R,G_R,H_R)}.

    If the neuron does not fire between times {t} and {t+1} (that is, if
    {X[i,t]=0}), the state at {t+1} is

      { V[i,t+1] = mu_V*V[i,t] + I[i,t] + \sum_{j} (X[j,t]*H[j,t]*w[j-->i]) }

      { G[i,t+1] = 1 - mu_G*(1 - G[i,t]) }

      { H[i,t+1] = 1 - mu_H*(1 - H[i,t]) }
    
  */

nmsim_neuron_parms_t* nmsim_neuron_parms_new
  ( double V_R,   /* Potential after firing. */
    double M_mu,  /* Potential decay factor. */
    double G_R,   /* Firing gain factor after firing. */
    double mu_G,  /* Firing gain recovery factor. */
    double H_R,   /* Output gain factor after firing. */
    double mu_H,  /* Output gain recovery factor. */
    struct nmsim_firing_func_t *Phi
  );
  /* Allocates and initializes a neuron parameter record. */
  
nmsim_neuron_parms_t* nmsim_neuron_parms_read(FILE *rd, double timeStep);
  /* Reads the parameters of a GL neuron model
    from file {rd}, in the format described 
    by {nmsim_neuron_parms_read_INFO} below. */
    
#define nmsim_neuron_parms_read_INFO \
  "The neuron description begins with a line \n" \
  "    begin " nmsim_neuron_parms_FILE_TYPE " (format of " nmsim_neuron_parms_VERSION ")\n" \
  "  followed by six lines in the format \"{NAME} = {VALUE}\",\n" \
  " where {NAME} is \"V_R\", \"tau_V\", \"G_R\", \"tau_G\",\n" \
  " \"H_R\", or \"tau_H\", in that" \
  " order; and {VALUE} is a fractional number.  Then\n" \
  " comes a line \"Phi = {PHI_NAME} {PHI_DEG} {PHI_V_M} {PHI_SIGMA}\", and" \
  " a closing line \"end nmsim_neuron_parms\"." \
  "\n" \
  "  There must be at least one space befor and after the \"=\" in the" \
  " parameter lines.  The paramter \"V_R\" is the /reset potential/," \
  " that the neuron assumes imemdiately after firing.  The parameters" \
  " \"G_R\" and \"H_R\" are the /reset gains/ for firing and" \
  " output, respectively.  The paramters \"tau_V\", \"tau_G\", and" \
  " \"tau_H\" are the charactristic times of potential decay and" \
  " gain recovery.  They are converted to the decay" \
  " factors {mu_V}, {mu_G}, and {mu_H} using the specified" \
  " time step of the simulation.\n" \
  "\n" \
  "  The parameters {PHI_NAME}, {PHI_DEG}, {PHI_V_M}, and {PHI_SIGMA} define" \
  " the firing function.\n" \
  "  " 
    
#define nmsim_neuron_parms_FILE_TYPE "nmsim_neuron_parms"
    
#define nmsim_neuron_parms_VERSION "2016-08-11"

#endif
