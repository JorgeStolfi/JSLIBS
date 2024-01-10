#ifndef nmsim_neuron_elem_H
#define nmsim_neuron_elem_H
 
/* Types and functions for the Galves-Löcherbach neuron model. */
/* Last edited on 2016-11-16 15:44:00 by stolfilocal */

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
  
  A synapse between a neuron {j} and a neuron {i} is described
  by a fixed /resting weight/ or /resting strength/ {w[j-->i]}, 
  positive or negative. 

  In the version implemented here, the state of the GL neuron number {i}
  consists of a /membrane potential/ {V[i,t]}, an /firing gain factor/ 
  {G[i,t]}, and an /output gain factor/ {H[i,t]}. 
  
  The factor {H[i,t]} modulates the strength of the output synapses of
  neuron {i}. It can be seen as modeling the temporary inhibition of
  synapses due to depletion of neurotransmitter vescicles near the
  synaptic junction. The factor {G[i,t]} modulates the neuron's
  propensity to fire; it can be seen as modeling the temporary
  inhibition of spiking due to activation or deactivation of ion
  channels in the neuron's body and in the initial segment of the axon.
  
  Between discrete times {t} and {t+1}, the neuron may /fire/ or not. If
  the neuron fires, its state returns to a fixed /reset state/ that is a
  time-invariant parameter of the model.  At the same time, pulse-like
  signals are sent to some set of /output neurons/.
  
  If the neuron fails to fire between {t} and {t+1}, its potential
  {V[i,t]} is incremented by a combination of signals received from some
  set of /input neurons/ and stimuli external to the network. However, a
  part of the old potential is ``leaked'', and only a fraction of it is
  retained in the new potential {V[i,t+1]}. At the same time, the gain
  factors {G[t]} and {H[t]} are also incremented, so as to approach the
  value 1 exponetially (as long as the neuron fails to fire).
  
  The sum of the external stimuli received by neuron {i} between
  times {t} and {t+1} are modeled as a potential increment {I[i,t]} 
  to be addded to its potential, if the  neuron does not fire.  
  
  We use {X[i,t]} to be a boolean variable that indicates whether neuron
  {i} fired between times {t} and {t+1}.distribution of each variable {X[i,t]}
  
  In the GL model, the firing of neuron {i} between times
  {t} and {t+1} is assumed to be dependent only
  on the potential {V[i,t]}, modulated by the firing gain factor
  {G[i,t]}, through a fixed \emph{firing function} {\Phi}.
  
  The /firing age/ of a neuron at time {t} is the number of whole time
  steps that elapsed between its last firing and time {t}. That is, the
  age of neuron {i} at time {t} is {k} iff the neuron fired betwen times
  {t-k-1} and {t-k}, and did not fire between {t-k} and {t}. In
  particular, the age at time {t} is zero iff the neuron fired between
  times {t-1} and {t}.
  
  Note that the state of a neuron depends only on its age {k} and on the
  inputs it received between sampling times {t-k} and {t}. This is a
  defining property of the Galves-Löcherbach model.
  
*/ 
  
  /* ??? Eliminate the 1-step resting state after firing? ??? */

#define _GNU_SOURCE
#include <stdio.h>

#include <nmsim_firing_func.h>

typedef struct nmsim_neuron_class_t
  { /* Parameters of potential dynamics: */
    double V_R;   /* Potential after firing. */
    double mu_V;  /* Potential decay factor. */
    /* Parameters of input gain dynamics: */
    double G_R;   /* Firing gain factor after firing. */
    double mu_G;  /* Firing gain recovery factor. */
    /* Parameters of output gain dynamics: */
    double H_R;   /* Output gain factor after firing. */
    double mu_H;  /* Output gain recovery factor. */
    /* Parameters of firing function: */
    struct nmsim_firing_func_t *Phi;
  } nmsim_neuron_class_t;
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

nmsim_neuron_class_t* nmsim_neuron_class_new
  ( double V_R,   /* Potential after firing. */
    double mu_V,  /* Potential decay factor. */
    double G_R,   /* Firing gain factor after firing. */
    double mu_G,  /* Firing gain recovery factor. */
    double H_R,   /* Output gain factor after firing. */
    double mu_H,  /* Output gain recovery factor. */
    struct nmsim_firing_func_t *Phi
  );
  /* Allocates and initializes a neuron parameter record. */
  
void nmsim_neuron_class_read(FILE *rd, double timeStep, nmsim_neuron_class_t*);
  /* Reads the parameters of a GL neuron model
    from file {rd}, in the format described 
    by {nmsim_neuron_class_read_INFO} below. */
    
#define nmsim_neuron_class_read_INFO \
  "The neuron description begins with a line \n" \
  "    begin " nmsim_neuron_class_FILE_TYPE " (format of " nmsim_neuron_class_VERSION ")\n" \
  "  followed by six lines in the format \"{NAME} = {VALUE}\",\n" \
  " where {NAME} is \"V_R\", \"tau_V\", \"G_R\", \"tau_G\",\n" \
  " \"H_R\", or \"tau_H\", in that" \
  " order; and {VALUE} is a fractional number.  Then\n" \
  " comes a line \"Phi = {PHI_NAME} {PHI_DEG} {PHI_V_M} {PHI_SIGMA}\", and" \
  " a closing line \"end nmsim_neuron_elem_parms\"." \
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
    
#define nmsim_neuron_class_FILE_TYPE "GL_neuron"
    
#define nmsim_neuron_class_VERSION "2016-08-11"

#endif
