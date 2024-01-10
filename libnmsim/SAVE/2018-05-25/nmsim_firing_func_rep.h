#ifndef nmsim_firing_func_rep_H
#define nmsim_firing_func_rep_H

#define _GNU_SOURCE
#include <stdint.h>
 
/* Internal representation of an {nmsim_firing_func_t}. */
/* Last edited on 2017-08-02 12:33:56 by jstolfi */

typedef void nmsim_firing_func_eval_proc_t
  ( double V, 
    double V_M, 
    double D_M,
    double deg,
    double *pr,
    double *dpr
  );
  /* Type of a procedure that evaluates the firing probability
    for a given membrane potential, for a certain class of firing functions. 
  
    The call {eval(V,dev,V_M,D_M,&pr,&dpr)} should return the firing
    probability {*pr=Phi(V)} and its derivative {*dpr = \partial_V
    Phi(V)} for a given membrane potential {V} (which should have been
    mutiplied by the firing modulator).
    
    The procedure should return {*pr = 1/2} and {*dpr = D_M}, for any
    {deg}, when {V = V_M}.
    
    If {dpr} is {NULL}, the derivative will not be computed.
    
    The {deg} parameter depends on the class, and may be ignored. */

typedef double nmsim_firing_func_eval_inv_proc_t
  ( double pr, 
    double deg,
    double V_M, 
    double D_M
  );
  /* Type of a procedure that computes the membrane potential needed to 
    get a certain firing probability, for a certain class of firing 
    functions.
  
    The call {eval_inv(pr,deg,V_M,D_M)} should return the membrane
    potential {V} that yields firing probability {*pr=Phi(V)} and its
    derivative {*dpr = \partial_V Phi(V)} for a given (which should
    have been mutiplied by the firing modulator).
    
    The procedure should return {*pr = 1/2} and {*dpr = D_M}, for any
    {deg}, when {V = V_M}.
    
    If {dpr} is {NULL}, the derivative will not be computed.
    
    The {deg} parameter depends on the class, and may be ignored. */

typedef struct nmsim_firing_func_t
  { nmsim_firing_func_eval_proc_t *eval;         /* Evaluates {Phi(V)}. */
    nmsim_firing_func_eval_inv_proc_t *eval_inv;  /* Evaluates {Phi^{-1}(pr)}. */
    double deg;      /* A shape parameter. */
    double V_M;      /* Midpoint reference potential. */
    double D_M;      /* The slope at the midpoint. */
    char *descr;     /* Description of function. */
  } nmsim_firing_func_t;
  /* Description of a firing function {Phi}, namely the member of 
    the class defined by {eval} and {eval_inv}, selected by the 
    parameters {deg,V_M,D_M}.
    
    The {descr} field is a string that contains a human-readable
    description of the firing function, including the values of
    the parameters {deg,V_M,D_M} when relevant. */  
    
nmsim_firing_func_t *nmsim_firing_func_new_gen
  ( nmsim_firing_func_eval_proc_t *eval, 
    nmsim_firing_func_eval_inv_proc_t *eval_inv, 
    double deg, 
    double V_M, 
    double D_M,
    char *descr
  );
  /* Creates a new {nmsim_firing_func_t} record with given fields. */ 

#endif
