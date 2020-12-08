/* intg_gen.h - generic ODE integrator. */
/* Last edited on 2009-02-11 01:28:09 by stolfi */

#ifndef intg_gen_H
#define intg_gen_H

#include <vec.h>

/*
  This interface defines the abstract data type "Integrator.T",
  an object that numerically solves ordinary differential equations
  of the form 
  
     {ds/dt = rhs(t,s)}
  
  where {s} is an unknown function from time {t} to some 
  ``state space'' {R^n}, and {rhs} is a client-provided function.
  
  Created 1995/06 by R.L.Liesenfeld, C.H.Q.Forster, and J.Stolfi.
  Converted to C 2003/08 by J.Stolfi. */
  
/* !!! The {step} method should call the RHS for the final endpoint
  internally, instead of leving that task to the user. That
  information is required for proper estimation of the integration
  error. E.g. in the Euler step {s1 = s0 + dt*v0} the error
  depends on the second derivative, which cannot be estimated
  from the known data ({s0}, {v0}). !!! */

typedef double Time;      /* A time value. */
typedef double Coord;     /* A state coordinate. */
typedef double Dist;      /* Distance between two coordinates. */
typedef double Speed;     /* Derivative of {Coord} relative to {Time}. */

typedef double_vec_t State;     /* The system's state. */
typedef double_vec_t Error;     /* Error estimate for a {State}. */
typedef double_vec_t Velocity;  /* Derivative of {State} relative to {Time}. */

typedef bool_t Intg_RHS(Time t, State *s, Velocity *v);
  /*
    Called by the integrator to evaluate the right-hand side of the
    differential equation. If the result is {TRUE}, the integration is 
    aborted. */

/* INTEGRATOR METHODS */

typedef void OBJ; /* Actually an Intg_T or a subclass thereof. */

typedef bool_t Intg_Mth_Step
  ( OBJ *self,        /* The integrator itself. */
    Intg_RHS rhs,     /* Computes the right-hand side {F(t,s)} */
    Time ta,          /* Starting time */
    State *sa,        /* (IN) Starting state */
    Velocity *va,     /* (IN) Derivatives at starting state */
    Time tb,          /* Final time */
    State *sb,        /* (OUT) Final state */
    Error *eb         /* (OUT) Error estimate for {sb} */
  );
  /* Integrates the equation {ds/dt == rhs(t,s)} from time {ta} and state
    {sa} to time {tb}, resulting in state {sb}.

    The method must be given in {va} the velocity (derivative of state
    with respect to time) {rhs(ta,sa)}, already computed. The method
    will typically call the {rhs} procedure in order to evaluate
    {F(ti,si)} for additional times {ti} in the open range {(ta__tb)},
    and the corresponding states {si}, generated internally.

    The method will return in {eb} an estimate of the absolute
    integration error affecting the computed state {sb}. This estimate
    is usually computed by comparing two different extrapolations of
    {sb}, and may be used to select the proper step size {tb - ta}.

    The vectors {sa}, {va}, {sb}, and {eb} must be memory-disjoint.

    The return value is TRUE if the integration was aborted because 
    {rhs} returned {TRUE}.  In that case, the contents of {*sb}
    and {*eb} are undefined.

    Fixed-step integration
    ----------------------

    In most integrators, the {step} method assumes implicitly that the
    function {s(t)} belongs to some restricted space {S}, for instance
    polynomials of a certain small degree. When that assumption is not
    true, the integration error affecting the final state {sb} may
    grow explosively as {tb-ta} increases.

    Therefore, integration over ``large'' intervals (where the
    function {s(t)} cannot be approximated by functions from
    {S}) is usually performed as a sequence of smaller {step}s,
    each one starting at the time and state where the previous
    one ended.

    In the simplest ``fixed stepsize'' scheme, the desired range
    {[ta__tb]} is divided into {NSteps} equal intervals. Here is a
    typical example:

      { Time tIni = (... initial time ...);
        Time tFin == (... final time ...);
        State s0 = double_vec_new(NCoords);
        State s1 = double_vec_new(NCoords);
        Velocity v0 = double_vec_new(NCoords);
        Error e = double_vec_new(NCoords);
       
        Time t0 = tIni;
        s0 = (... initial state ...);
        for (k = 1; k <= NSteps; k++)
          { v0 = RHS(t0, s0, v0);
            double r == ((double)i)/((double)NSteps);
            Time t1 = (1 - r)*tIni + r*tFin;
            intg.step(RHS, t0,s0,v0, t1,s1, e);
            t0 = t1; s0 = s1;
          }
      }

    For a more sophisticated scheme, with variable stepsize,
    see the comments of {Intg_Mth_Adjust} below. */

typedef Time Intg_Mth_Adjust
  ( OBJ *self,             /* The integrator itself. */
    Time dt, 
    Dist error, 
    Dist tol, 
    Time dtMin, 
    Time dtMax
  );
  /* Given that a time step of size {dt} generated the specified
    integration {error}, returns the time step that should make the
    integration error approximately {tol}.

    The result is clipped to the range {[dtMin _ dtMax]}.

    This method assumes that {error} is some mathematical norm of the
    estimated error vector {eb} returned by {step}. That is,
    multiplying {eb} by some positive factor should multiply {error}
    by that same factor.

    Here is a typical example of integration with ``adaptive stepsize
    control'':

      {
        Time tIni = (... initial time ...);
        Time tFin = (... final time ...);
        State s0 = double_vec_new(NCoords);
        State s1 = double_vec_new(NCoords);
        Velocity v0 = double_vec_new(NCoords);
        Error e = double_vec_new(NCoords);
      
        Time t0 = tIni;
        s0 = (... initial state ...);
        Time dt = (... initial stepsize ...);
        while (t0 < tb)
          { RHS(t0, s0, v0);
            Time t1 = MIN(t0 + dt, tFin);
            while (1)
              { intg.step(RHS, t0,s0,v0, t1,s1, e);
                Dist error = rn_norm(e);
                if ((error < tol) || (dt < dtMin)) { break;}
                dt = intg.adjustStepSize(dt, error, tol, dtMin, dtMax);
                t1 = MIN(t0 + dt, tFin);
              }
            t0 = t1;
            s0 = s1;
         }
      }

    Note that the error estimates returned by {step} may be unreliable
    if the right-hand side behaves too wildly in the specified
    interval. Threfore, this adaptive scheme should be used with
    caution.

    In particular, if {rhs(t,s)} has discontinuities along the path from
    {(t0,s0)} to {(t1, s1)}, then it is advisable to break the
    integration into two separate steps at the discontinuity. */


/* GENERIC INTEGRATOR OBJECT */

typedef struct Intg_T
  { char *type;               /* Identifies the concrete type. */
    char *descr;              /* Printeable description of integrator. */
    Intg_Mth_Step *step;      /* Single-step integration method. */
    Intg_Mth_Adjust *adjust;  /* Stepsize adjustment method. */
  } Intg_T;
  /* A generic integrator for ordinary differential equations.
    Clients should create instances of appropriate concrete
    subclasses, such as {RKF4Integrator}. */

#define Intg_TypeId "INTG."

#endif


