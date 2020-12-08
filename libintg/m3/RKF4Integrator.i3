INTERFACE RKF4Integrator;

(* A 4th-order Runge-Kutta integrator for ordinary differential equations. *)
(* Created 95/06 by R.L.Liesenfeld, C.H.Q.Forster, and J.Stolfi.           *)

IMPORT Integrator;

TYPE
  T <: Public;
    (*
      An adaptive Runge-Kutta integrator for ordinary differential 
      equations.
      
      The "step" method performs one step of the Runge-Kutta-Fehlberg
      method, which is theoretically accurate to fourth order.  The
      error estimate is the difference between the fourth- and
      fifth-order Runge-Kutta approximators. *)
    
  Public = Integrator.T OBJECT METHODS
    END;

END RKF4Integrator.



