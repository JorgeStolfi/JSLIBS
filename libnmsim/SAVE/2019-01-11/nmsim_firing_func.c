/* See {nmsim_firing_func.h} */
/* Last edited on 2016-08-11 16:30:06 by stolfilocal */

#define _GNU_SOURCE

#include <nmsim_firing_func.h>
typedef struct nmsim_firing_func_t
  { nmsim_firing_proc_t *eval;            /* Evaluates {Phi(V)}. */
    double V_M;      /* The midpoint potential. */
    double D_M;      /* The slope at the midpoint. */
    double deg;      /* A shape parameter. */
    char *desc;      /* Description of function. */
  } nmsim_firing_func_t;
  /* Description of a firing function {Phi}, namely the member of 
    the class defined by {eval} selected by the parameters
    {V_M,D_M,deg}.
    
    The {desc} field is a string that contains a human-readable
    description of the firing function, including the values of
    the parameters {V_M,D_M,deg} when relevant. */    

typedef void nmsim_firing_proc_t
  ( double V, 
    double V_M, 
    double D_M,
    double deg,
    double *pr,
    double *dpr
  );
  /* Type of a procedure that implements a certain class
    of firing functions. 
  
    The call {eval(V,V_M,D_M,deg,&pr,&dpr)} should return the firing
    probability {*pr=Phi(V)} and its derivative {*dpr = \partial_V
    Phi(V)} for a given membrane potential {V} (which should have been
    mutiplied by the firing gain factor).
    
    The procedure should return {*pr = 1/2} and {*dpr = D_M}, for any
    {deg}, when {V = V_M}.
    
    If {dpr} is {NULL}, the derivative will not be computed.
    
    The {deg} parameter depends on the class, and may be ignored. */

void nmsim_firing_func_eval
  ( nmsim_firing_func_t *Phi,
    double V, 
    double *pr,
    double *dpr
  )
  {
    Phi->eval(V, Phi->V_M, Phi->D_M, Phi->deg, pr, dpr);
  }
 
nmsim_firing_func_t *nmsim_firing_func_gauss(double V_M, double D_M)
  {
    char *desc = NULL
    asprintf(
    nmsim_firing_func_t *Phi = nmsim_firing_func_new(desc,proc,0,V_M,D_M);
    return Phi;
  }
