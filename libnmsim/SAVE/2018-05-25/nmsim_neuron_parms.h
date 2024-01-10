#ifndef nmsim_neuron_H
#define nmsim_neuron_H
 
/* Types and functions for the Galves-Löcherbach neuron model. */
/* Last edited on 2018-05-25 17:54:43 by jstolfi */

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
  number) and a /firing age/ {a[t]}. 
  
  In this library, the of membrane potential {V[t]} is taken to be its
  /resting potential/ {V_B}, when there has been no inputs and no firings for
  a long time.
  
  The firing age {a[t]} is the number of whole time steps that elapsed
  between its last firing and time {t}. That is, {a[t]} is {k} iff the
  neuron fired betwen times {t-k-1} and {t-k}, and did not fire between
  {t-k} and {t}. Thus, {a[t+1]} is zero if {X[t]} is 1, otherwise it
  is {a[t]+1}.
  
  All other parameters of the neuron are determined by the potential and
  age, and therefore do not depend on what happened in the network
  before the time of the last firing. This is a defining property of the
  Galves-Löcherbach model.
  
  If the neuron fails to fire between {t} and {t+1}, a part of its
  potential {V[t]} is ``leaked'', and only a fraction {c[t]*V[t]}
  of it is retained in the new potential {V[t+1]}.
  
  Anyway, from {t} to {t+1} the potential is incremented by a
  combination of signals received from some set of /input neurons/
  through chemical synapses and gap junctions, and stimuli external to
  the network. 
  
  The external inputs are represented by a potential increment {I[t]}
  due the external stimuli received by neuron {i} between times {t} and
  {t+1}.
  
  A (chemical) synapse between a neuron {j} and a neuron {i} is
  described by a fixed /resting weight/ or /resting strength/
  {w_B[j-->i]}, positive or negative. The actual strength {w[j-->i,t]}
  is {H[j][t]*w_B[j-->i]*G[i][t]} where {H[j][t]} and {G[i][t]} are
  factors that depend on the age of the respective neuron.
  
  The factor {H[t]} modulates the strength of the output synapses of
  neuron {j}. It can be seen as modeling the temporary inhibition of
  synapses due to depletion of neurotransmitter vescicles near the
  synaptic junction. The factor {G[t]} modulates the strength of
  all input synapses of neuron {i}.
  
  !!! Can we get Hebbian learning from this formulation? !!!
  
  !!! Taking again {G[t]} out of the {\Phi} function because it is
  cleaner. Multiplying {V[t]} does not make sense because it is
  relative to an arbitrary reference potential. Adding a bias to
  {V[t]} is equivalent to changing the reset potential and leak
  coefficient. Maybe multiply {G[t]} by the probability? !!!
  
  Similarly, a gap junction between a neuron {j} and a neuron {i} is
  described by a fixed non-negative /resting transmittance/
  {r_B[j-->i]}. The actual transmittance is {F[j][t]*r_B[j-->i]*E[i][t]}
  where {F[j][t]} and {E[i][t]} are factors that depend on the
  firing ages of the two neurons.
  
  !!! Gap juntions. !!!
  
  !!! Modulation of leakage coefficient. !!!
  
  In the GL model, the firing of neuron {i} between times {t} and {t+1}
  is assumed to be dependent only on the potential {V[t]} through a
  fixed \emph{firing function} {\Phi}.
  
  Note that the state of a neuron depends only on its age {k} and on the
  inputs it received between sampling times {t-k} and {t}.
  
*/ 

#define _GNU_SOURCE
#include <stdio.h>

#include <nmsim_firing_func.h>

typedef struct nmsim_neuron_parms_t
  {/* Potential dynamics: */
    double V_B;   /* Resting potential. */
    double V_R;   /* Potential after firing. */
    double c_B;   /* Resting leakage. */
    /* Leakage dynamics parameters: */
    double L_R;   /* Reset value of leakage modulator. */
    double L_mu;  /* Recovery factor of leakage modulator. */
    /* Chemical synapse dynamics parameters: */
    double G_R;   /* Reset value of input synapse modulator. */
    double G_mu;  /* Recovery factor of output synapse modulator. */
    double H_R;   /* Reset value of input synapse modulator. */
    double H_mu;  /* Recovery factor of output synapse modulator. */
    /* Parameters of firing function:  */
    struct nmsim_firing_func_t *Phi;  /* Firing function. */
  } nmsim_neuron_parms_t;
  /* Parameters of a GL neuron.
  
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

nmsim_neuron_parms_t *nmsim_neuron_parms_new
  ( double V_B,   
    double V_R,   
    double c_B,   
    double L_R,   
    double L_mu,  
    double G_R,   
    double G_mu,  
    double H_R,   
    double H_mu,  
    struct nmsim_firing_func_t *Phi
  );
  /* Allocates and initializes a neuron parameter record. */
  
void nmsim_neuron_parms_write(FILE *wr, nmsim_neuron_parms_t *parms);
  /* Writes the parameter record {p} to file {wr}, as described by
  {nmsim_neuron_parms_read_INFO} below. */

nmsim_neuron_parms_t* nmsim_neuron_parms_read(FILE *rd, double timeStep);
  /* Reads the parameters of a GL neuron model
    from file {rd}, in the format described 
    by {nmsim_neuron_parms_read_INFO} below. */
    
#define nmsim_neuron_parms_read_INFO \
  "The neuron description begins with a line \n" \
  "    begin " nmsim_neuron_parms_FILE_TYPE " (format of " nmsim_neuron_parms_VERSION ")\n" \
  "  followed by six lines in the format \"{NAME} = {VALUE}\",\n" \
  " where {NAME} is \"V_R\", \"V_B\", \"c_B\", \"L_R\", \"L_tau\", \"G_R\", \"G_tau\",\n" \
  " \"H_R\", or \"H_tau\", in that" \
  " order; and {VALUE} is a fractional number.  Then\n" \
  " comes a line \"Phi = {PHI_NAME} {PHI_DEG} {PHI_V_M} {PHI_SIGMA}\", and" \
  " a closing line \"end nmsim_neuron_parms\"." \
  "\n" \
  "  There must be at least one space befor and after the \"=\" in the" \
  " parameter lines.  The paramter \"V_R\" is the /reset potential/," \
  " that the neuron assumes imemdiately after firing.  The parameters" \
  " \"G_R\" and \"H_R\" are the /reset gains/ for firing and" \
  " output, respectively.  The parameters \"L_tau\", \"G_tau\", and" \
  " \"H_tau\" are the charactristic times of recovery for modulating factors.  They" \
  " are converted to the decay" \
  " factors {L_mu}, {G_mu}, and {H_mu} using the specified" \
  " time step of the simulation.\n" \
  "\n" \
  "  The parameters {PHI_NAME}, {PHI_DEG}, {PHI_V_M}, and {PHI_SIGMA} define" \
  " the firing function.\n" \
  "  " 
    
#define nmsim_neuron_parms_FILE_TYPE "nmsim_neuron_parms"
    
#define nmsim_neuron_parms_VERSION "2016-08-11"

#endif
