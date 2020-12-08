INTERFACE Integrator;

(*
  This interface defines the abstract data type "Integrator.T",
  an object that numerically solves ordinary differential equations
  of the form 
  |
  |    ds/dt = F(t,s)
  |
  where "s" is an unknown function from time "t" to some 
  ``state space'' "R^n", and "F" is a client-provided function.
  
  Created 95/06 by R.L.Liesenfeld, C.H.Q.Forster, and J.Stolfi
*)

EXCEPTION Abort;

TYPE
  Time = LONGREAL;            (* A time value *)
  Coord = LONGREAL;           (* A state coordinate *)
  State = ARRAY OF Coord;     (* The system's state. *)
  Dist  = LONGREAL;           (* Distance between two states *)
  Error = ARRAY OF Coord;     (* Error estimate for a "State". *)
  Speed = LONGREAL;           (* Derivative of "Coord" relative to "Time" *)
  Velocity = ARRAY OF Speed;  (* Derivative of "State" relative to "Time". *)

  T  = OBJECT METHODS
    (*
      A generic integrator for ordinary differential equations.
      
      This is an abstract object class that implements no methods.
      Therefore, clients should not create instances of "Integrator.T"
      directly.  Instead, they should create instances of appropriate
      concrete subclasses, such as "RKF4Integrator.T". *)

      step(
          diff: DiffProc;        (* Computes the right-hand side "F(t,s)" *)
          ta: Time;              (* Starting time *)
          READONLY sa: State;    (* Starting state *)
          READONLY va: Velocity; (* Derivatives at starting state *)
          tb: Time;              (* Final time *)
          VAR sb: State;         (* OUT: Final state *)
          VAR eb: Error;         (* OUT: Error estimate for "sb" *)
        ) RAISES {Abort};
        (*
          Integrates the equation "ds/dt = F(t,s)" from time "ta"
          and state "sa" to time "tb", resulting in state "sb".
          
          The method must be given in "va" the velocity (derivative of
          state with respect to time) "F(ta,sa)", already computed.
          The method will typically call the "diff" procedure in order
          to evaluate "F(ti,si)" for additional times "ti" in the open
          range "(ta__tb)", and the corresponding states "si",
          generated internally.
          
          The method will return in "eb" an estimate of the 
          absolute integration error affecting the computed state
          "sb".  This estimate is usually computed by comparing
          two different extrapolations of "sb", and may be used
          to select the proper step size "tb - ta".
          
          The vectors "sa", "va", "sb", and "er" must be 
          memory-disjoint.
          
          The "Abort" exception is not raised by "step"; is provided
          so that "diff" may interrupt the integration at mid-step.
          
          Fixed-step integration
          ----------------------
          
          In most integrators, the "step" method assumes implicitly
          that the function "s(t)" belongs to some restricted space "S", 
          for instance polynomials of a certain small degree.
          When that assumption is not true,  the integration error
          affecting the final state "sb" may grow explosively as 
          "tb-ta" increases.  
          
          Therefore, integration over ``large'' intervals (where the
          function "s(t)" cannot be approximated by functions from
          "S") is usually performed as a sequence of smaller "step"s,
          each one starting at the time and state where the previous
          one ended.
          
          In the simplest ``fixed stepsize'' scheme, the desired range
          "[ta__tb]" is divided into "NSteps" equal intervals. 
          Here is a typical example:

            | WITH
            |   tIni = (... initial time ...),
            |   tFin = (... final time ...),
            |   s = NEW(REF State, NCoords)^,
            |   v = NEW(REF State, NCoords)^,
            |   e = NEW(REF State, NCoords)^
            | DO
            |   VAR t: LONGREAL := tIni;
            |   BEGIN
            |     s := (... initial state ...);
            |     v := (... initial velocity ...);
            |     FOR k := 1 TO NSteps DO
            |       WITH
            |         r = FLOAT(i, LONGREAL)/FLOAT(NSteps, LONGREAL),
            |         tn = (1.0d0 - r)*tIni + r*tFin
            |       DO
            |         intg.step(Diff, t,s,v, tn,s, e);
            |         t := tn;
            |         Diff(s, v); 
            |       END
            |     END
            |   END
            | END

          For a more sophisticated scheme, with variable stepsize,
          see the comments in "adjustStepSize" below. *)
         
      adjustStepSize(dt: Time; error, tol: Dist; dtMin, dtMax: Time): Time;
        (*
          Given that a time step of size "dt" generated the specified 
          integration "error", returns the time step that should 
          make the integration error approximately "tol".
          
          The result is clipped to the range "[dtMin _ dtMax]".
          
          This method assumes that "error" is some mathematical norm
          of the estimated error vector "eb" returned by "step".  That
          is, multiplying "eb" by some positive factor should multiply
          "error" by that same factor.
          
          Here is a typical example of integration with 
          ``adaptive stepsize control'':
          
            | WITH
            |   tIni = (... initial time ...),
            |   tFin = (... final time ...),
            |   s0 = NEW(REF State, NCoords)^,
            |   s1 = NEW(REF State, NCoords)^,
            |   v0 = NEW(REF State, NCoords)^,
            |   e = NEW(REF State, NCoords)^
            | DO
            |   VAR t0: LONGREAL := tIni; dt, t1: LONGREAL;
            |   BEGIN
            |     s0 := (... initial state ...);
            |     v0 := (... initial velocity ...);
            |     dt := (... initial stepsize ...);
            |     WHILE t0 < tb DO
            |       t1 := MIN(t0 + dt, tFin);
            |       LOOP
            |         intg.step(Diff, t0,s0,v0, t1,s1, e);
            |         WITH error = LRN.Norm(e) DO
            |           IF error < tol OR dt < dtMin THEN EXIT END;
            |           dt := intg.adjustStepSize(dt, error, tol, dtMin, dtMax)
            |         END
            |       END;
            |       t0 := t1;
            |       s0 := s1;
            |       Diff(t0, s0, v0)
            |     END
            |   END
            | END

          Note that the error estimates returned by "step" may be
          unreliable if the right-hand side behaves too wildly 
          in the specified interval.  Threfore, this adaptive
          scheme should be used with caution.
          
          In particular, if "F(t,s)" has discontinuities along the
          path from "(t0,s0)" to "(t1, s1)", then it is advisable to
          break the integration into two separate steps at the
          discontinuity. *)
        
      name(): TEXT;
        (*
          Returns the integrator's name, possibly with
          any internal parameters. *)
    END;

  DiffProc = PROCEDURE (t: Time; READONLY s: State; VAR v: Velocity) RAISES {Abort};
    (*
      Called by the integrator to evaluate the right-hand side of the
      differential equation. *)

END Integrator.

