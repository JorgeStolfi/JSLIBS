INTERFACE EulerIntegrator;

(* An adaptive Euler integrator for ordinary differential equations. *)
(* Created 95/06 by R.L.Liesenfeld, C.H.Q.Forster, and J.Stolfi.     *)

IMPORT Integrator;

TYPE
  T <: Public;
    (*
      An adaptive Euler integrator for ordinary differential equations.
      
      The "step" method performs the trivial Euler method, i.e. 
      linear extrapolation, which is only accurate to first order.  
      The error estimate at each step is the difference between the
      first- and second-order approximators. *)
    
  Public = Integrator.T OBJECT METHODS
    END;

END EulerIntegrator.

