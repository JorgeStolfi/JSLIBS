/* intg_Euler.h - Euler integrator for ODEs. */
/* Last edited on 2007-01-04 00:17:18 by stolfi */

#ifndef intg_Euler_H
#define intg_Euler_H

#include <vec.h>
#include <intg_gen.h>

/* 
  An adaptive Euler integrator for ordinary differential equations.
  Created 1995/06 by R.L.Liesenfeld, C.H.Q.Forster, and J.Stolfi. 
  Converted to C 2003/08 by J.Stolfi. */

typedef struct Intg_Euler_T
  {
    Intg_T intg;     /* Stuff inherited from the parent {Intg_T} class. */
  } Intg_Euler_T;
  /*
    An adaptive Euler integrator for ordinary differential equations,
    packaged as a subclass of {Intg_T}.

    The {step} method performs the trivial Euler method, i.e. 
    linear extrapolation, which is only accurate to first order.  
    The error estimate at each step is the difference between the
    first- and second-order approximators. */

#define Intg_Euler_TypeId "INTG.EUL."
    
Intg_Euler_T *Intg_Euler_new(void);
  /* Returns a new Euler integrator object. */
    
#endif


