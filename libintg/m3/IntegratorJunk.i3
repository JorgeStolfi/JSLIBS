INTERFACE Integrator;

  T <: Public;
    (*
      This is a semi-abstract object class that implements only the
      basic "setSize" and a generic "integrate" method.
      
      Normally, a subclass of "Integrator.T" should override only the
      "step", "adjustStepSize", and "name" methods, and inherit
      "integrate" from this type. *)
  
  Public = OBJECT METHODS

      integrate(
          VAR t: Time;
          VAR s: State;
          diff: DiffProc;
          checkState: CheckStateProc;
          checkStep: CheckStepProc := NIL;
        );
        (*
          Integrates the equation "ds/dt = F(t, s)" from time "t"
          and state "s", until a limit set by the "checkState" and
          "checkStep" procedures below.
          
          Upon entry, "t" should contain the initial time, and "s" the
          corresponding initial state. Upon exit, "t" will contain the
          final time, and "s" will contain the corresponding final
          state.
          
          Integration will proceed by a sequence of chained steps.
          The length of the first step will be "dt"; that of later
          steps is set by the client during the integration, 
          
          The integrator will call "diff(t, s, v)" several times per
          step to compute the derivative "v = ds/dt = F(t,s)" for a
          state "s" and time "t".  Beware that those times "t" are not
          always increasing, and that the states "s" may be widely off
          the true path.  The first call to "diff" will be applied to
          the initial time and state given.
          
          Stepsize control and discrete event handling
          --------------------------------------------
          
          The "checkState" and "checkStep" arguments allow the 
          client to adjust the step size based on the current
          state, and also handle ``discrete events'' such as 
          preset checktimes and discontinuous differential equations.

          The "checkState" procedure
          --------------------------
          
          The integrator will call "checkState(t, s)", for every valid
          state in the path (not including the ``tentative'' states
          mentioned above).  In particular, "checkState" will be
          applied to the first and last states of the path.
          
          The "checkState" procedure is expected to 

            (1) ``use up'' the state (print it, plot it, etc.),
            (2) adjust "s" as necessary (see below),
            (3) compute the corresponding velocity "v", and
            (4) return a final time "tf" for the next step.
          
          If the returned time "tf" is greater than "t", the next step
          will stretch from "t" to "tf".  If "tf" is equal to "t" (or
          earlier), the integration will end at time "t" and state "s".
          
          The procedure "checkState" is allowed to modify the state
          "s" before returning to the integrator.  For instance, the
          procedure may adjust "s" so as to satisfy any constraints or
          conservation laws which are implied by the original
          differential equations but could be violated in a long
          integration due to the accumulation of error.
          
          Usually, "checkState" will call "diff(t,s,v)" to compute
          the velocity.  However, if the right-hand side "F(t, s)"
          has a discontinuity at "t", then "diff"
          should compute the left limit, whereas "checkState"
          should compute the right limit.
            
          For a simple integration with fixed step "dt" and
          predetermined stopping time "tMax", the "checkState"
          procedure could simply be
          
             |   PROCEDURE CheckState(...): Time = 
             |     BEGIN
             |       ... use state ...
             |       IF t >= tMax THEN RETURN tMax END;
             |       Diff(t, s, v);
             |       RETURN MAX(t + dt, tMax)
             |     END CheckState;
          
          The "checkStep" procedure
          -------------------------
          
          If the "checkStep" procedure is not NIL, the integrator will
          also call "checkStep(ta, sa, va, tb, sb, er)" after each
          integration step, to decide ``a posteriori'' whether the
          step should indeed have been taken.
          
          The arguments of that call say that the step just taken went
          from time "ta" and state "sa" to time "tb" and state "sb";
          that the derivative at "ta" was "va"; and that the estimated
          integration error vector at "tb" for that step is "er".
          
          The "checkStep" procedure should ponder this data, 
          and then return
            
            (1) any time "tc" in "[ta _ tb)", meaning that the step
                was not acceptable, and the integration should be redone
                from "ta" to "tc"; or
            
            (2) the time "tb" itself, meaning that the step 
                was acceptable, and the system should be brought
                to state "sb".
                 
          In case (1), after the integration from "ta" to "tc"
          is redone, the "checkStep" procedure will be called again
          to inspect the truncated step, and perhaps refine the estimate
          of the final time "tc".  The client should make sure that
          this loop terminates after a finite numebr of steps.
          
          In case (2), after "checkStep" returns, the integrator 
          will call "checkState", as explained above, 
          to process the final state "sb" and define the length
          of the next step.
          
          The programmer is advised  use the "adjustStepSize" method below to 
          select the proper step size.  Here is the standard way to obtain 
          adaptive stepsize control:
          
             | WITH 
             |   intg = NEW(MyIntegrator.T).setSize(N),
             |   t = NEW(REF Time)^,
             |   s = NEW(REF State, N)^
             | DO
             | 
             |   PROCEDURE Diff(...) =
             |     BEGIN 
             |       ... compute v ...
             |     END Diff;
             | 
             |   VAR dtNext: Time := dtStart;
             | 
             |   PROCEDURE CheckState(...): Time = 
             |     BEGIN
             |       ... use state ...
             |       IF t >= tMax THEN 
             |         RETURN t
             |       ELSE
             |         ... compute v ...
             |         RETURN MIN(tMax, t + dtNext)
             |       END
             |     END CheckState;
             | 
             |   PROCEDURE CheckStep(...): Time =
             |     BEGIN
             |       WITH 
             |         dt = tb - ta,
             |         error = LRN.LInfNorm(er)
             |       DO
             |         dtNext := s.adjustStepSize(dt, error, tol, dtMax, dtMin);
             |         IF error > tol AND dt > p.dtMin THEN 
             |           RETURN ta + dtNext
             |         ELSE
             |           RETURN tb
             |         END
             |       END
             |     END CheckStep;
             |     
             |   BEGIN
             |     ... initialize "t", "s" ...
             |     s.integrate(t, s, Diff, CheckState, CheckStep);
             |     ... use final state, time ...
             |   END
             | END
             
          The "CheckStep" procedure above may use other error norms,
          e.g. "LRN.L1Norm" or "LRN.Length", as appropriate. 
        *)
        
      adjustStepSize(dt: Time; error, tol: Coord; dtMin, dtMax: Time): Time;
        (*
          Given that a time step of size "dt" generated the specified 
          integration "error", returns the time step that should 
          make the integration error approximately "tol".
          
          The result is clipped to the range "[dtMin _ dtMax]".
          
          This method assumes that "error" is some Minkowski metric
          of the estimated error vector "er" given to "checkStep";
          that is, multiplying "er" by some positive factor should
          multiply "error" by that same factor.
        *)
        
      step(
          diff: DiffProc;
          ta: Time;              (* Starting time *)
          READONLY sa: State;    (* Starting state *)
          READONLY va: Velocity; (* Derivatives at starting state *)
          tb: Time;              (* Final time *)
          VAR sb: State;         (* OUT: Final state *)
          VAR er: Error;         (* OUT: Error estimate *)
        );
        (*
          A low-level method that performs a single step of the
          integrator, extrapolating to the system's path from state "sa" 
          at time "ta" to state "sb" at time "tb". 

          The procedure must be given the velocity "va" at "sa", already
          computed.  The "diff" procedrue may be called for additional
          intermediate states, depending on the particular integrator used.
          The estimated integration error for the step is returned in "er".
          
          This method should not be called while a call to "integrate"
          is in progress. *)
        
      name(): TEXT;
        (*
          Returns the integrator's name, possibly with
          any internal parameters. *)
    END;

  DiffProc = PROCEDURE (t: Time; READONLY s: State; VAR v: Velocity);
    (*
      Called by the integrator to evaluate the right-hand side of the
      differential equation. *)

  CheckStateProc = PROCEDURE (t: Time; VAR s: State; VAR v: Velocity): Time;
    (*
      Called by the integrator for every new state along the path.
      Should adjust the state "s" if necessary, compute the velocity
      "v", and return the upper time limit for the next step. *)

  CheckStepProc = PROCEDURE (
      ta: Time;
      READONLY sa: State;
      READONLY va: Velocity; 
      tb: Time;
      READONLY sb: State;
      READONLY er: Error;
    ): Time;
    (*
      Called by the integrator after computing the step from 
      state "sa" at time "ta" to state "sb" at time "tb"..
      Should return some time "tc" in "[ta _ tb)", to cancel the
      step and try again from "ta" to "tc"; or the time "tb"
      itself, to make the step definitive. *)

END Integrator.

