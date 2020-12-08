MODULE RKF4Integrator;

IMPORT Math, Integrator;
FROM Math IMPORT sqrt;
FROM Integrator IMPORT DiffProc, Time, Coord, State, Velocity, Error, Abort;

REVEAL T = Public BRANDED OBJECT
    v: ARRAY [1..5] OF REF Velocity;
    s: REF State;
  OVERRIDES
    name := Name;
    adjustStepSize := AdjustStepSize;
    step := Step;
  END;

PROCEDURE Name(<*UNUSED*> ri: T): TEXT =
  BEGIN
    RETURN "RKF4Integrator.T"
  END  Name;

PROCEDURE AdjustStepSize(
    <*UNUSED*> g: T; 
    dt: Time; 
    error, tol: Coord; 
    dtMin, dtMax: Time;
  ): Time =
  CONST MaxScale = 10.0D0;
  CONST MinScale = 1.0D0/MaxScale;
  CONST Safety   = 0.840896415253D0;
  BEGIN
    IF error > 0.0d0 THEN
      WITH scale = Safety * sqrt(sqrt(tol/error)) DO
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
    ta: Time;              (* Starting time *)
    READONLY sa: State;    (* Starting state *)
    READONLY va: Velocity; (* Derivatives at starting state *)
    tb: Time;              (* Final time *)
    VAR sb: State;         (* OUT: Final state *)
    VAR er: Error;         (* OUT: Error estimate *)
  ) RAISES {Abort} =
  CONST
    C1a = 1.0d0/4.0d0;
    C1t = C1a;
    C2a = 3.0d0/32.0d0;
    C21 = 9.0d0/32.0d0;
    C2t = C2a + C21;
    C3a = 1932.0d0/2197.0d0;
    C31 = -7200.0d0/2197.0d0;
    C32 = 7296.0d0 /2197.0d0;
    C3t = C3a + C31 + C32;
    C4a = 439.0d0 /216.0d0;
    C41 = -8.0d0;
    C42 = 3680.0d0 /513.0d0;
    C43 = -845.0d0 /4104.0d0;
    C4t = C4a + C41 + C42 + C43;
    C5a = -8.0d0 /27.0d0;
    C51 = 2.0d0;
    C52 = -3544.0d0 /2565.0d0;
    C53 = 1859.0d0 /4104.0d0;
    C54 = -11.0d0 /40.0d0;
    C5t = C5a + C51 + C52 + C53 + C54;
    Cea = 1.0d0 /360.0d0;
    Ce2 = -128.0d0 /4275.0d0;
    Ce3 = -2197.0d0 /75240.0d0;
    Ce4 = 1.0d0 /50.0d0;
    Ce5 = 2.0d0 /55.0d0;
    Cba = 25.0d0 /216.0d0;
    Cb2 = 1408.0d0 /2565.0d0;
    Cb3 = 2197.0d0 /4104.0d0;
    Cb4 = -1.0d0 /5.0d0;
  BEGIN
    (* Allocate work areas if needed: *)
    WITH n = NUMBER(sa) DO 
      IF g.s = NIL OR NUMBER(g.s^) # n THEN 
        g.s := NEW(REF State, n);
        FOR i := 1 TO 5 DO 
          g.v[i] := NEW(REF Velocity, n);
        END;
      END
    END;
  
    WITH
      n = NUMBER(sa),
      dt = tb - ta,
      s = g.s^,
      v1 = g.v[1]^,
      v2 = g.v[2]^,
      v3 = g.v[3]^,
      v4 = g.v[4]^,
      v5 = g.v[5]^
    DO
      (* 
        We could do without the internal work vector "g.s" 
        by using "sb" instead.  However, if the client 
        gave us the same vector for "sa" and "sb", the 
        result would be garbage.  
        
        We could also use "v5" instead of "s" in all steps
        excet the degree 5 extrapolation, where we could use "v1". 
        But that would be too obscure...
        
        It is safer this way...
      *)
      
      (* Degree 1 extrapolation: *)
      FOR i := 0 TO n-1 DO
        s[i] := sa[i] + dt * C1a*va[i];
      END;
      diff(ta + C1t * dt, s, v1);

      (* Degree 2 extrapolation: *)
      FOR i := 0 TO n-1 DO
        s[i] := sa[i] + dt * (C2a*va[i] + C21*v1[i]);
      END;
      diff(ta + C2t * dt, s, v2);

      (* Degree 3 extrapolation: *)
      FOR i := 0 TO n-1 DO
        s[i] := sa[i] + dt * (C3a*va[i] + C31*v1[i] + C32*v2[i]);
      END;
      diff(ta + C3t * dt, s, v3);
      
      (* Degree 4 extrapolation: *)
      FOR i := 0 TO n-1 DO
        s[i] := sa[i] + dt * (C4a*va[i] + C41*v1[i] + C42*v2[i] + C43*v3[i]);
      END;
      diff(ta + C4t * dt, s, v4);
      
      (* Degree 5 extrapolation: *)
      FOR i := 0 TO n-1 DO
        s[i] := sa[i] + dt*(C5a*va[i] + C51*v1[i] + C52*v2[i] + C53*v3[i] + C54*v4[i]);
      END;
      diff(ta + C5t * dt, s, v5);
      
      (* Final extrapolation and error estimate: *)
      FOR i := 0 TO n-1 DO
        sb[i] := sa[i] + dt*(Cba*va[i] + Cb2*v2[i] + Cb3*v3[i] + Cb4*v4[i]);
        er[i] := Cea*va[i] + Ce2*v2[i] + Ce3*v3[i] + Ce4*v4[i] + Ce5*v5[i]
      END;

    END
  END Step;

BEGIN 
END RKF4Integrator.
 
