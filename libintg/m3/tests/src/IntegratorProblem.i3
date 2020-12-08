INTERFACE IntegratorProblem;

FROM Integrator IMPORT Time, State, Velocity;

TYPE 
  T = OBJECT
    n: CARDINAL;
    tag: TEXT;
    t0: Time := -1.0d0;
    t1: Time := +1.0d0;
    dt: Time := 0.1d0;
    dtMin: Time := 1.0d-6;
    dtMax: Time := 1.0d0;
  METHODS
    solution (t: Time; VAR s: State);
    diff (t: Time; READONLY s: State; VAR v: Velocity);
  END;
  
PROCEDURE Sample(): REF ARRAY OF T;
 (*
   Returns a sample of ODE integration problems. *)

END IntegratorProblem.
