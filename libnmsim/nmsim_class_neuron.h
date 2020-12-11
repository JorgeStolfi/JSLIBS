#ifndef nmsim_class_neuron_H
#define nmsim_class_neuron_H
 
/* Common parameters for a class of Galves-Löcherbach neurons. */
/* Last edited on 2020-12-11 17:51:26 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <vec.h>

#include <nmsim_basic.h>
#include <nmsim_write.h>
#include <nmsim_read.h>
#include <nmsim_compare.h>
#include <nmsim_firing_func.h>

typedef struct nmsim_class_neuron_t
  { /* Potential dynamics: */
    double V_B;   /* Resting potential (mV). */
    double V_R;   /* Potential after firing (mV). */
    double c_B;   /* Resting recharge factor. */
    /* Recharge modulator parameters: */
    double M_R;   /* Reset value of recharge modulator. */
    double M_mu;  /* Recovery factor of recharge modulator. */
    /* Synapse output strength modulator parameters: */
    double H_R;   /* Reset value of output synapse modulator. */
    double H_mu;  /* Recovery factor of output synapse modulator. */
    /* Parameters of firing function:  */
    nmsim_firing_func_t Phi;  /* Firing function. */
  } nmsim_class_neuron_t;
  /* Parameters of a GL neuron class. See {nmsim_class_neuron_INFO} for details. */
  
#define nmsim_class_neuron_INFO \
  "The parameter {V_R} is the /reset potential/," \
  " that the neuron assumes imemdiately after firing.  The parameters" \
  " {M_R} and {H_R} are the /reset values/ for the modulators" \
  " of recharge factor and " \
  " output synapse strengths, respectively.  The parameter {V_tau} is" \
  " the characteristic recharge time, apart" \
  " from recharge modulation.  {M_tau} and {H_tau} are" \
  " the characteritic recovery times for the modulating" \
  " factors {M} and {H}, respectively.\n" \
  "\n" \
  "   The parameter {V_M} is the potential at which the" \
  " firing function {\Phi} is {1/2}, and" \
  " {V_D} is the half-width of the transtion from `low' to" \
  " `high'.  The parameter {Phi_class} is the type of {\Phi}:\n" \
  "\n" \
  nmsim_firing_func_class_INFO "\n" \
  "\n" \
  "  The values of {V_tau}, {M_tau} and {H_tau} are" \
  " functionally equivalent to the parameters {c_B}, {M_mu}, and {H_mu} of the" \
  " model as usually presented, but are independent of the " \
  " simulation time step {\Delta}."

nmsim_class_neuron_t *nmsim_class_neuron_new
  ( double V_B,   
    double V_R,   
    double c_B,   
    double M_R,   
    double M_mu,  
    double H_R,   
    double H_mu,  
    struct nmsim_firing_func_t *Phi
  );
  /* Allocates and initializes a neuron parameter record. */
   
typedef nmsim_class_neuron_t *nmsim_class_neuron_ref_t;

void nmsim_class_neuron_free(nmsim_class_neuron_t *nclass);
  /* Releases the storage of {*nclass}. */
  
/* SIMULATION TOOLS */

double nmsim_class_neuron_compute_M(nmsim_class_neuron_t *nclass, nmsim_step_count_t age);
  /* Returns the recharge modulator {M} for a neuron
    of the given class and age, as defined by {nclass->{M_R,M_mu}}.
    The {age} must be zero or positive. */

double nmsim_class_neuron_compute_H(nmsim_class_neuron_t *nclass, nmsim_step_count_t age);
  /* Returns the output synapse strength modulator {H} for a neuron
    of the given class and age, as defined by {nclass->{H_R,H_mu}}.
    The {age} must be zero or positive. */
    
double nmsim_class_neuron_recharge(nmsim_class_neuron_t *nclass, double V, double M);
  /* Applies one step of the recharge process for a neuron of class
    {nclass},  currently with potential {V} and recharge modulator {M}
    at some time {t}. Namely, reduces the difference between {V} 
    and the rest potential {nclass->V_B} by the recharge factor {nclass->{c_B}}
    times the recharge modulator {M}.  Returns the modified potential. */  

void nmsim_class_neuron_throw_state
  ( nmsim_class_neuron_t *nclass, 
    double *VP, 
    nmsim_step_count_t *ageP
  );
  /* Chooses a suitable initial potential {*VP} and a non-negative initial age {*ageP}
    for a neuron of the given class.  The input values of {*VP} and {*ageP} are
    ignored. */

/* INPUT/OUTPUT */

void nmsim_class_neuron_write(FILE *wr, nmsim_class_neuron_t *nclass, double timeStep);
  /* Writes the parameter record {p} to file {wr}, as described by
    {nmsim_class_neuron_read_INFO} below.  The {timeStep} parameter (in milliseconds) is used to 
    convert {c_B}, {M_mu}, and {H_mu} to step-invariant {V_tau}, {M_tau}, and {H_tau}.
    The output takes multiple lines and ends with a newline. */

nmsim_class_neuron_t* nmsim_class_neuron_read(FILE *rd, double timeStep);
  /* Reads the parameters of a GL neuron class 
    from file {rd}, in the format described 
    by {nmsim_class_neuron_read_INFO} below.  The {timeStep} parameter (in milliseconds) is used to 
    convert step-invariant parameters {V_tau}, {M_tau}, and {H_tau} 
    from the file into the step-dependent parameters {c_B}, {M_mu}, and {H_mu}. */
    
#define nmsim_class_neuron_read_INFO \
  "The neuron parameters block begins with a line \n" \
  "    begin " nmsim_class_neuron_FILE_TYPE " (format of " nmsim_class_neuron_VERSION ")\n" \
  " followed by lines in the format \"{NAME} = {VALUE}\"," \
  " where {NAME} is \"V_R\", \"V_B\", \"V_tau\", \"M_R\", \"M_tau\",\n" \
  " \"H_R\", or \"H_tau\", or \"Phi_class\", \"V_M\", \"V_D\", in that order; and" \
  " a closing line \"end " nmsim_class_neuron_FILE_TYPE "\".\n" \
  "\n" \
  "  There must be at least one space before and after the \"=\" in the" \
  " parameter lines. The {VALUE} for \"Phi_class\" is a single uppercase" \
  " letter without quotes. For all other parameters, the {VALUE} is a" \
  " decimal fraction.  The values \"V_R\", \"V_B\", \"V_M\", and \"V_D\" are" \
  " in millivolts; those of \"V_tau\", \"M_tau\", or \"H_tau\" are" \
  " in milliseconds;  and those of \"M_R\" and \"H_R\" are" \
  " adimensional.\n" \
  "\n" \
  "  " nmsim_class_neuron_INFO
    
#define nmsim_class_neuron_FILE_TYPE "nmsim_class_neuron"
    
#define nmsim_class_neuron_VERSION "2019-01-08"

/* TESTING AND DEBUGGING */

void nmsim_class_neuron_show(FILE *wr, char *pref, nmsim_class_neuron_t *nclass, char *suff);
  /* Writes the parameters {nclass} as a compact line, with prefix {pref} 
    and suffix {suff}. */
 
nmsim_class_neuron_t* nmsim_class_neuron_throw(void);
  /* Generates a random neuron class record, for testing purposes. */

void nmsim_class_neuron_compare(nmsim_class_neuron_t *nclass_read, nmsim_class_neuron_t *nclass_orig);
  /* Compares the neuron class parameters {nclass_read} read from a file
    with the expected class parameters {nclass_orig}.  Aborts if any parameter
    does not match (within roundoff tolerance). */

vec_typedef(nmsim_class_neuron_ref_vec_t,nmsim_class_neuron_ref_vec,nmsim_class_neuron_ref_t);

#endif
   
