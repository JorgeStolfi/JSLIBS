MODULE IntegratorProblem;

FROM Integrator IMPORT Time, State, Velocity;
FROM Math IMPORT sqrt, sin, cos, exp, log;

PROCEDURE Sample(): REF ARRAY OF T =
  BEGIN
    WITH
      p = ARRAY OF T{
        NEW(
          T OBJECT OVERRIDES
              solution := QuadrSol;
              diff := QuadrDiff;
            END,
          tag := "quadr", n := 2
        ),
        NEW(
          T OBJECT OVERRIDES
              solution := QuartSol;
              diff := QuartDiff;
            END,
          tag := "quart", n := 2
        ),
        NEW(
          T OBJECT OVERRIDES
              solution := QuintSol;
              diff := QuintDiff;
            END,
          tag := "quint", n := 2
        ),
        NEW(
          T OBJECT OVERRIDES
              solution := SinusSol;
              diff := SinusDiff;
            END,
          tag := "sinus", n := 2
        ),
        NEW(
          T OBJECT OVERRIDES
              solution := GaussSol;
              diff := GaussDiff;
            END,
          tag := "gauss", n := 2
        ),
        NEW(
          T OBJECT OVERRIDES
              solution := TwangSol;
              diff := TwangDiff;
            END,
          tag := "twang", n := 2
        ),
        NEW(
          T OBJECT OVERRIDES
              solution := HyperSol;
              diff := HyperDiff;
            END,
          tag := "hyper", n := 2
        ),
        NEW(
          T OBJECT OVERRIDES
              solution := NotchSol;
              diff := NotchDiff;
            END,
          tag := "notch", n := 1
        )
      },
      rp = NEW(REF ARRAY OF T, NUMBER(p))
    DO
      rp^ := p;
      RETURN rp
    END
  END Sample;

(* Quadratic *)

PROCEDURE QuadrSol(<*UNUSED*> p: T; t: Time; VAR s: State) =
  BEGIN
    s[0] := t*t;
    s[1] := 2.0d0 * t;
  END QuadrSol;

PROCEDURE QuadrDiff(
    <*UNUSED*> p: T; 
    <*UNUSED*> t: Time; 
    READONLY s: State; 
    VAR v: State
  ) =
  BEGIN
    v[0] := s[1];
    v[1] := 2.0d0;
  END QuadrDiff;

(* Quartic *)

PROCEDURE QuartSol(<*UNUSED*> p: T; t: Time; VAR s: State) =
  BEGIN
    s[0] := t*(t*t*t - 1.0d0);
    s[1] := 4.0d0 * t*t*t - 1.0d0;
  END QuartSol;

PROCEDURE QuartDiff(
    <*UNUSED*> p: T; 
    t: Time; 
    READONLY s: State; 
    VAR v: Velocity
  ) =
  BEGIN
    v[0] := s[1];
    v[1] := 12.0d0*t*t;
  END QuartDiff;

(* Quintic *)

PROCEDURE QuintSol(<*UNUSED*> p: T; t: Time; VAR s: State) =
  BEGIN
    s[0] := t*(t*t*t*t - 1.0d0);
    s[1] := 5.0d0 * t*t*t*t - 1.0d0;
  END QuintSol;

PROCEDURE QuintDiff(
    <*UNUSED*> p: T; 
    t: Time; 
    READONLY s: State; 
    VAR v: Velocity
  ) =
  BEGIN
    v[0] := s[1];
    v[1] := 20.0d0*t*t*t;
  END QuintDiff;

(* Sinusoid *)

CONST SinusW = 20.0d0;

PROCEDURE SinusSol(<*UNUSED*> p: T; t: Time; VAR s: State) =
  BEGIN
    WITH
      W = SinusW,
      ct = cos(W*t),
      st = sin(W*t)
    DO
      s[0] := ct;
      s[1] := - W*st;
    END
  END SinusSol;

PROCEDURE SinusDiff(
    <*UNUSED*> p: T; 
    <*UNUSED*> t: Time; 
    READONLY s: State; 
    VAR v: Velocity
  ) =
  BEGIN
    WITH
      W = SinusW
    DO
      v[0] := s[1];
      v[1] := - W*W*s[0];
    END
  END SinusDiff;

(* Decaying sinusoid *)

CONST 
  TwangW = 20.0d0;
  TwangK = 2.0d0;

PROCEDURE TwangSol(<*UNUSED*> p: T; t: Time; VAR s: State) =
  BEGIN
    WITH
      W = TwangW,
      K = TwangK,
      et = exp(-K*t),
      ct = et * cos(W*t),
      st = et * sin(W*t)
    DO
      s[0] := ct;
      s[1] := st;
    END
  END TwangSol;

PROCEDURE TwangDiff(
    <*UNUSED*> p: T; 
    <*UNUSED*> t: Time; 
    READONLY s: State; 
    VAR v: Velocity
  ) =
  BEGIN
    WITH
      W = TwangW,
      K = TwangK
    DO
      v[0] := - K * s[0] - W * s[1];
      v[1] := - K * s[1] + W * s[0];
    END
  END TwangDiff;

(* Gaussian bell *)

CONST GaussK = 3.0d0;

PROCEDURE GaussSol(<*UNUSED*> p: T; t: Time; VAR s: State) =
  BEGIN
    WITH
      K = GaussK,
      kt = K*t,
      et = exp(-kt*kt)
    DO
      s[0] := et;
      s[1] := -2.0d0*K*K*t*et;
    END
  END GaussSol;

PROCEDURE GaussDiff(
    <*UNUSED*> p: T; 
    <*UNUSED*> t: Time; 
    READONLY s: State; 
    VAR v: Velocity
  ) =
  BEGIN
    WITH
      K = GaussK,
      r = s[1]/s[0]
    DO
      v[0] := s[1];
      v[1] := s[0]*(r*r - 2.0d0*K*K)
    END
  END GaussDiff;

(* Hyperbola *)

CONST HyperE = 0.05d0;

PROCEDURE HyperSol(<*UNUSED*> p: T; t: Time; VAR s: State) =
  BEGIN
    WITH
      E = HyperE,
      et = sqrt(t*t + E*E)
    DO
      s[0] := et;
      s[1] := t/et;
    END
  END HyperSol;

PROCEDURE HyperDiff(
    <*UNUSED*> p: T; 
    <*UNUSED*> t: Time; 
    READONLY s: State; 
    VAR v: Velocity
  ) =
  BEGIN
    WITH
      E = HyperE,
      r = 1.0d0/s[0]
    DO
      v[0] := s[1];
      v[1] := E*E * r*r*r
    END
  END HyperDiff;

(* Notches: log(2 + sin(t)) *)

CONST 
  NotchK = 1.41421356237309504880D0;
  NotchW = 4.0D0;

PROCEDURE NotchSol(<*UNUSED*> p: T; t: Time; VAR s: State) =
  BEGIN
    WITH
      K = NotchK,
      W = NotchW,
      et = log(K + sin(W*t))
    DO
      s[0] := et;
    END
  END NotchSol;

PROCEDURE NotchDiff(
    <*UNUSED*> p: T; 
    t: Time; 
    READONLY s: State; 
    VAR v: Velocity
  ) =
  BEGIN
    WITH
      W = NotchW
    DO
      v[0] := W*cos(W*t)/exp(s[0])
    END
  END NotchDiff;

BEGIN
END IntegratorProblem.
