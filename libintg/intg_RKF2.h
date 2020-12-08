/* intg_RKF2.h - Runge-Kutta-Fehlberg 2nd order integrator for ODEs. */
/* Last edited on 2007-01-04 00:17:20 by stolfi */

#ifndef intg_RKF2_H
#define intg_RKF2_H

#include <vec.h>
#include <intg_gen.h>

/* 
  A second-order Runge-Kutta-Fehlberg integrator for ordinary
  differential equations.
  Created 1995/06 by R.L.Liesenfeld, C.H.Q.Forster, and J.Stolfi. 
  Converted to C 2003/08 by J.Stolfi. */

typedef struct Intg_RKF2_T
  {
    Intg_T intg;     /* Stuff inherited from the parent {Intg_T} class. */
    State s;         /* Used internally by "step" */
    Velocity v;      /* Ditto */
  } Intg_RKF2_T;
  /*
    A Runge-Kutta-Fehlberg integrator for ordinary differential equations,
    packaged as a subclass of {Intg_T}.

    The {step} method performs a second-order RKF2 method.  
    The error estimate at each step is the difference between the
    first- and second-order approximators. */

#define Intg_RKF2_TypeId "INTG.RKF2."
    
Intg_RKF2_T *Intg_RKF2_new(void);
  /* Returns a new RKF2 integrator object. */
    
#endif


