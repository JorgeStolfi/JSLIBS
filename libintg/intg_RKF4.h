/* intg_RKF4.h - RKF4 integrator for ODEs. */
/* Last edited on 2007-01-04 00:17:23 by stolfi */

#ifndef intg_RKF4_H
#define intg_RKF4_H

#include <vec.h>
#include <intg_gen.h>

/* 
  A 4th-order Runge-Kutta integrator for ordinary differential equations.
  Created 1995/06 by R.L.Liesenfeld, C.H.Q.Forster, and J.Stolfi.
  Converted to C 2003/09 by J.Stolfi. */

typedef struct Intg_RKF4_T
  {
    Intg_T intg;     /* Stuff inherited from the parent {Intg_T} class. */
    State s;         /* Used internally by "step" */
    Velocity v[5];   /* Ditto */
  } Intg_RKF4_T;
  /*
    An adaptive Runge-Kutta-Fehlberg integrator for ordinary
    differential equations, packaged as a subclass of {Intg_T}.

    The {step} method performs one step of the Runge-Kutta-Fehlberg
    method, which is theoretically accurate to fourth order.  The
    error estimate is the difference between the fourth- and
    fifth-order Runge-Kutta approximators. */

#define Intg_RKF4_TypeId "INTG.RKF4."
    
Intg_RKF4_T *Intg_RKF4_new(void);
  /* Returns a new RKF4 integrator object. */

#endif
