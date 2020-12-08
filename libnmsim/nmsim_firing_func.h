#ifndef nmsim_firing_func_H
#define nmsim_firing_func_H

#define _GNU_SOURCE
#include <stdint.h>
 
/* Firing functions for the Galves-Löcherbach neuron model. */
/* Last edited on 2020-12-04 21:04:39 by jstolfi */

typedef char nmsim_firing_func_class_t;
  /* A byte that describes the type of firing function.  The cuurrently
    defined classes are in {nmsim_firing_func_class_INFO}. */

#define nmsim_firing_func_class_INFO \
  "    'G' = integral of Gaussian distr with mean {V_M} and deviation {V_D}.\n" \
  "\n" \
  "    'L' = linear ramp from 0 at {V_M - V_D} to 1 at {V_M + V_D}.\n" \
  "\n" \
  "    'N' = Nilton Kamiji's function, 0.4 power from 0 at {V_M - V_D} to 1 at {V_M + V_D}."
     
typedef struct nmsim_firing_func_t
  { nmsim_firing_func_class_t class; /* Type of firing function. */
    double V_M;   /* Midpoint voltage (mV). */
    double V_D;   /* Voltage blur (mV). */
  } nmsim_firing_func_t;
  /* Describes a firing function {\Phi}.  
  
    The {.class} field is the general form, e.g. "G" for "Gaussian
    integral". The actual shape is defined by the fields {.V_M} and
    {V_D}.
    
    For all classes, {.V_M} is the internal potential at which the
    probability of firing is {1/2}. The precise meaning of {.V_D} depends on
    the class, but it defines the horizontal stretch of the function:
    {\Phi} will vary from "low" to "high" approximately in the range
    {V_M - V_D} to {V_M + V_D}. */

nmsim_firing_func_t nmsim_firing_func_make
  ( nmsim_firing_func_class_t class,
    double V_M,
    double V_D
  );
  /* Assembles a record of type {nmsim_firing_func_t}  with the given fields. */

void nmsim_firing_func_eval
  ( nmsim_firing_func_t *Phi,
    double V, 
    double *prP,
    double *dprP
  );
  /* Evaluates the firing function {Phi}  at the argument {V}.
    
    Stores the firing probability {Phi(V)} into {*prP}, and the derivative 
    of {Phi} at that point into {*dprP}.
    
    The procedure should return {*pr = 1/2}, for any {class}, when {V = Phi.V_M}.
    The derivative at {V = Phi.V_M} will depend on {Phi.class}.
    
    If any of {prP} and/or {dprP} is {NULL}, the corresponding 
    parameter will not be computed. */

double nmsim_firing_func_eval_inv(nmsim_firing_func_t *Phi, double pr);
  /* Inverts the firing function {Phi}. Namely, computes the potential {V} such that
    {Phi(V) = pr}.
    
    The procedure should return {V = Phi.V_M} when {*pr = 1/2}, for any {class}.
    
    Fails if {pr < 0} or {pr > 1}.  Otherwise returns 
    {-INF} or {+INF} if {Phi(V)} is greater than or less than {pr},
    respectively, for all finite {V}. */

double nmsim_firing_func_eval_gauss
  ( nmsim_firing_func_t *Phi,
    double V_avg, 
    double V_dev
  );
  /* Evaluates the firing function {Phi} at a random argument {V}
    drawn from a Gaussian distribution with mean {V_avg}
    and deviation {V_dev}. */

void nmsim_firing_func_compare(nmsim_firing_func_t *Phi_read, nmsim_firing_func_t *Phi_orig);
  /* Compares function {Phi_read} read from a file with the expected 
    function {Phi_orig}.  Fails with error message if they are 
    different classes and the parameters are too different. */

#endif
