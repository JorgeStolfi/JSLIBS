MODULE EulerIntegrator;

IMPORT Math, Integrator;
FROM Math IMPORT sqrt;
FROM Integrator IMPORT DiffProc, Time, Coord, State, Velocity, Error, Abort;

REVEAL T = Public BRANDED OBJECT
    v: REF Velocity; (* Used internally by "step" *)
    s: REF State;    (* Ditto *)
  OVERRIDES
    name := Name;
    adjustStepSize := AdjustStepSize;
    step := Step;
  END;

PROCEDURE Name(<*UNUSED*> g: T): TEXT =
  BEGIN
    RETURN "EulerIntegrator.T"
  END  Name;

PROCEDURE AdjustStepSize(
    <*UNUSED*> g: T; 
    dt: Time; 
    error, tol: Coord; 
    dtMin, dtMax: Time;
  ): Time =
  CONST MaxScale = 10.0D0;
  CONST MinScale = 1.0D0/MaxScale;
  CONST Safety   = 0.95D0;
  BEGIN
    IF error > 0.0d0 THEN
      WITH scale = Safety * sqrt(tol/error) DO
        IF scale <= MinScale  THEN dt := MinScale * dt
        ELSIF scale >= MaxScale  THEN dt := MaxScale * dt
        ELSE dt := scale * dt
        END
      END
    ELSE
      dt := MaxScale * dt;
    END;
    RETURN MIN(dtMax, MAX(dtMin, dt))
  END AdjustStepSize;

PROCEDURE Step(
    g: T;
    diff: DiffProc;
    ta: Time;
    READONLY sa: State;
    READONLY va: Velocity;
    tb: Time;
    VAR sb: State; (* Final state *)
    VAR er: Error; (* Error estimate *)
  ) RAISES {Abort} =
  CONST 
    C1a = 1.0d0/2.0d0;
    C1t = C1a;
  BEGIN
    (* Allocate work areas if needed: *)
    WITH n = NUMBER(sa) DO
      IF g.s = NIL OR NUMBER(g.s^) # n THEN
        g.s := NEW(REF State, n);
        g.v := NEW(REF Velocity, n);
      END
    END;
    
    WITH
      n = NUMBER(sa),
      dt = tb - ta,
      s = g.s^,
      v1 = g.v^
    DO
      (* Degree 1 extrapolation: *)
      FOR i := 0 TO n-1 DO
        s[i] := sa[i] + dt * C1a*va[i];
      END;
      diff(ta + C1t*dt, s, v1);

      (* Error estimation (hack --- should check the thory and do it right!): *)
      FOR i := 0 TO n-1 DO
        sb[i] := sa[i] + dt*va[i];
        er[i] := 0.5D0 * dt * (va[i] - v1[i])
      END;
    END
  END Step;

BEGIN 
END EulerIntegrator.
 
