MODULE Integrator;

REVEAL T = Public BRANDED OBJECT
    (* Used internally by the "integrate" method *)
    v: REF Velocity;
    sb: REF State;
    er: REF Error;
  OVERRIDES
    integrate := Integrate;
  END;

PROCEDURE Integrate(
    g: T;
    VAR t: Time;
    VAR s: State;
    diff: DiffProc;
    checkState: CheckStateProc;
    checkStep: CheckStepProc;
  ) =
  VAR tb, tf: Time;
  BEGIN
    (* Allocate work areas if needed: *)
    WITH n = NUMBER(s) DO
      IF g.sb = NIL OR NUMBER(g.sb^) # n THEN
        g.v := NEW(REF Velocity, n);
        g.sb := NEW(REF State, n);
        g.er := NEW(REF Error, n);
      END
    END;
    
    WITH
      v = g.v^,
      sb = g.sb^,
      er = g.er^
    DO
      LOOP
        (* Report state, get velocity and suggested step size: *)
        tf := checkState(t, s, v);
        IF tf <= t THEN RETURN END;
        (* Keep integrating from "t" with decreasing steps until client is happy: *)
        REPEAT
          tb := tf;
          g.step(diff, t, s, v, tb, sb, er);
          IF checkStep # NIL THEN 
            tf := checkStep(t, s, v, tf, sb, er)
          END;
        UNTIL tf >= tb;
        (* Phew, we can move on: *)
        t := tb; 
        s := sb;
      END;
    END;
  END Integrate;

BEGIN
END Integrator.
