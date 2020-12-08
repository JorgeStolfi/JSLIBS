MODULE TestIntegrators EXPORTS Main;

IMPORT Integrator, IntegratorProblem, EulerIntegrator, RKF4Integrator;
IMPORT Wr, Fmt, LRN, Text, Thread, OSError, FileWr;
FROM Integrator IMPORT State, Velocity, Time, Error, Coord;
FROM Math IMPORT exp, log;
FROM Stdio IMPORT stderr, stdout;

TYPE 
  Problem = IntegratorProblem.T;
  
  Grater = OBJECT
    i: Integrator.T;
    tag: TEXT
  END;
  
  TestResults = RECORD
    nEvals: CARDINAL;    (* Number of RHS evaluations *)
    nSteps: CARDINAL;    (* Number of steps really taken *)
    nRedos: CARDINAL;    (* Number of steps discarded because of large error *)
    tf: LONGREAL;        (* Stopping time *)
    ec, ek: LONGREAL;    (* Estimated integration error was approximately "c*dt^k" *)
    mc, mk: LONGREAL;    (* Measured integration error was approximately "c*dt^k" *)
  END;

  Stats = RECORD 
    xx, xy, xu, yu, uu: LONGREAL := 0.0d0;  (* Error/stepsize sums *)
    (* 
      Used for the least-squares fit of integration error.
      Let "dt[i]" be the length of step "i", "error[i]" be the 
      corresponding error, "x[i]" be "log(dt[i])",
      and "y[i]" be "log(error[i])".  Then 
        "xx = sum(x[i]*x[i])",
        "xy = sum(x[i]*y[i])",
        "xu = sum(x[i])",
        "yu = sum(y[i])",
        "uu = sum(1)",
      summed over all steps actually taken.
    *)
  END;

PROCEDURE DoIt() =
  BEGIN
    WITH
      p = IntegratorProblem.Sample()^,
      g = ARRAY OF Grater{
        NEW(Grater, tag := "rkfo4", i := NEW(RKF4Integrator.T)),
        NEW(Grater, tag := "euler", i := NEW(EulerIntegrator.T))
      }
    DO
      FOR i := 0 TO LAST(p) DO
        PlotTrueSolution(p[i], 512);
      END;
      WriteTestResultsHeader(stdout);
      FOR j := 0 TO LAST(g) DO
        DoStepTest(g[j], p, 1.0d-3);
        FOR i := 0 TO LAST(p) DO
          IF Text.Equal(p[i].tag, "notch") THEN 
            DoStandardTest(p[i], g[j], "1", 1.0d-1);
            DoStandardTest(p[i], g[j], "2", 1.0d-2);
            DoStandardTest(p[i], g[j], "3", 1.0d-3);
            DoStandardTest(p[i], g[j], "6", 1.0d-6);
          END
        END;
        WriteEOL(stdout);
      END
    END
  END DoIt;
 
PROCEDURE DoStandardTest(p: Problem; g: Grater; tag: TEXT; tol: Coord) =
  <* FATAL Wr.Failure, Thread.Alerted, OSError.E *>
  BEGIN
    WITH 
      name = p.tag & "-" & g.tag & "-" & tag,
      wr = FileWr.Open(name & ".plot")
    DO
      VAR 
        trFree  := AdaptiveTest(p, g, wr,  tol, force := FALSE, rStep := 1.0d0);
        trBound := AdaptiveTest(p, g, NIL, tol, force := TRUE,  rStep := 1.0d0);
      BEGIN
        trFree.mc := trBound.mc;
        trFree.mk := trBound.mk;
        WriteTestResultsLine(stdout, name, trFree);
      END;
      Wr.Close(wr)
    END;    
  END DoStandardTest;
  
PROCEDURE DoStepTest(
    g: Grater; 
    READONLY p: ARRAY OF Problem;
    tol: LONGREAL;
  ) =
  CONST NScales = 51;
        MinScale = 0.75d0;
        MaxScale = 1.50d0;
  <* FATAL Wr.Failure, Thread.Alerted, OSError.E *>
  BEGIN
    WITH 
      wr = FileWr.Open(g.tag & "-step" & ".plot"),
      tr = NEW(REF ARRAY OF TestResults, NUMBER(p))^,
      trRef = NEW(REF ARRAY OF TestResults, NUMBER(p))^
    DO
      WriteStepTestPlotHeader(wr, g.i.name(), p);
      FOR k := 1 TO NScales DO
        WITH 
          e = FLOAT(k-1, LONGREAL)/FLOAT(NScales-1, LONGREAL),
          scale = MinScale * exp(log(MaxScale/MinScale)*e)
        DO
          FOR i := 0 TO LAST(p) DO 
            tr[i] := AdaptiveTest(p[i], g, NIL, tol, force := TRUE, rStep := scale);
            IF k = 1 THEN trRef[i] := tr[i] END;
          END;
          WriteStepTestPlotLine(wr, scale, tr, trRef);
        END;
      END;        
      Wr.Close(wr)
    END;    
  END DoStepTest;
  
PROCEDURE PlotTrueSolution(p: Problem; nSteps: CARDINAL) =
  <* FATAL Wr.Failure, Thread.Alerted, OSError.E *>
  VAR t: Time;
  BEGIN
    WITH 
      wr = FileWr.Open(p.tag & "-true" & ".plot"),
      dt = (p.t1-p.t0)/FLOAT(nSteps, LONGREAL)
    DO
      WriteAdaptiveTestPlotHeader(wr, "true solution", p.dt, p.dt, p.dt, 0.0d0);
      (* Plot true solution: *)
      WITH 
        s = NEW(REF State, p.n)^,
        v = NEW(REF State, p.n)^
      DO
        FOR i := 0 TO nSteps DO
          t := p.t0 + dt * FLOAT(i, LONGREAL);
          p.solution(t, s);
          p.diff(t, s, v);
          WriteAdaptiveTestPlotLine(wr, t, s, v);
        END;
      END;
      Wr.Close(wr)
    END
  END PlotTrueSolution;

PROCEDURE AdaptiveTest(
    p: Problem;                (* Problem *)
    g: Grater;                 (* Integrator *)
    wr: Wr.T;                  (* Plot file, or NIL *)
    tol: Coord;                (* Error tolerance *)
    force: BOOLEAN := FALSE;   (* Forces path to follow the true solution *)
    rStep: LONGREAL := 1.0d0;  (* Extra stepsize factor after "g.i.adjustStepSize" *)
  ): TestResults =

  <* FATAL Wr.Failure, Thread.Alerted, Integrator.Abort *>

  VAR nEvals: CARDINAL := 0;

  PROCEDURE Diff(t: Time; READONLY s: State; VAR v: Velocity) =
    BEGIN
      INC(nEvals);
      p.diff(t, s, v)
    END Diff;

  VAR eStats: Stats;          (* Estimated error × dt statistics. *)
  VAR mStats: Stats;          (* Measured error × dt statistics. *)

  PROCEDURE GatherErrorStats(VAR s: Stats; dt, error: LONGREAL) =
    BEGIN
      IF error > 1.0d-8 AND dt > 1.0d-5 THEN
        WITH
          x = 0.5D0 * log(dt*dt + 1.0D-200),
          y = 0.5D0 * log(error*error + 1.0D-200)
        DO
          StatsGather(s, x, y)
        END;
      END
    END GatherErrorStats;

  VAR trace: BOOLEAN := FALSE;
  VAR nRedos: CARDINAL := 0;  (* Steps redone because of large error. *)
  VAR nSteps: CARDINAL := 0;  (* Steps actually taken. *)
  VAR t0: Time := p.t0; tPrev, dt, dtNext, t1: Time;
  BEGIN
    IF wr # NIL THEN
      WriteAdaptiveTestPlotHeader(wr, g.i.name(), p.dt, p.dtMin, p.dtMax, tol)
    END;
    IF trace THEN Wr.PutText(stderr, "\n") END;
    (* Plot numerical solution: *)
    WITH
      s0 = NEW(REF State, p.n)^,
      s1 = NEW(REF State, p.n)^,
      v0 = NEW(REF State, p.n)^,
      e = NEW(REF Error, p.n)^
    DO
      p.solution(t0, s0);
      Diff(t0, s0, v0);
      dt := p.dt;
      WHILE t0 < p.t1 DO
        IF wr # NIL THEN WriteAdaptiveTestPlotLine(wr, t0, s0, v0) END;
        IF force AND t0 > p.t0 THEN
          v0 := s0;
          p.solution(t0, s0);
          WITH dt = t0 - tPrev, error = LRN.Dist(v0, s0) DO
            GatherErrorStats(mStats, dt, error)
          END;
          p.diff(t0, s0, v0);
        END;
        tPrev := t0;
        t1 := MIN(t0 + dt, p.t1);
        LOOP
          g.i.step(Diff, t0,s0,v0, t1,s1, e);
          WITH error = LRN.Norm(e) DO
            IF error < tol OR dt < p.dtMin THEN 
              GatherErrorStats(eStats, dt, error);
              EXIT
            END;
            (* Ops, last step was too long, reduce it: *)
            IF trace AND nSteps < 200 THEN Trace('*') END;
            INC(nRedos);
            dtNext := g.i.adjustStepSize(rStep * dt, error, tol, p.dtMin, p.dtMax);
            <* ASSERT dtNext < rStep * dt *>
            dt := dtNext
          END
        END;
        IF trace AND nSteps < 200 THEN Trace('-') END;
        INC(nSteps);
        t0 := t1;
        s0 := s1;
        Diff(t0, s0, v0);
      END
    END;
    IF trace THEN Wr.PutText(stderr, "\n") END;
    IF wr # NIL THEN Wr.Flush(wr) END;
    VAR ek, ec, mk, mc: LONGREAL := 0.0d0;
    BEGIN
      StatsFit(eStats, a := ek, b := ec); ec := exp(ec);
      IF force THEN StatsFit(mStats, a := mk, b := mc); mc := exp(mc) END;
      RETURN TestResults{
          nSteps := nSteps, nRedos := nRedos, nEvals := nEvals,
          ec := ec, ek := ek,
          mc := mc, mk := mk,
          tf := t1
        }
    END;
  END AdaptiveTest;

PROCEDURE WriteTestResultsLine(wr: Wr.T; name: TEXT; READONLY tr: TestResults) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stdout, Fmt.Pad(name, 15, align := Fmt.Align.Left));
    Wr.PutText(wr, " ");
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(tr.nSteps), 7));
    Wr.PutText(wr, " ");
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(tr.nRedos), 7));
    Wr.PutText(wr, " ");
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(tr.nEvals), 7));
    Wr.PutText(wr, " ");
    Wr.PutText(wr, " ");
    Wr.PutText(wr, Fmt.Pad(Fmt.LongReal(tr.ec, Fmt.Style.Fix, 4), 8));
    Wr.PutText(wr, " * dt ^");
    Wr.PutText(wr, Fmt.Pad(Fmt.LongReal(tr.ek, Fmt.Style.Fix, 2), 5));
    Wr.PutText(wr, " ");
    IF tr.mc = 0.0d0 THEN
      Wr.PutText(wr, Fmt.Pad("", 20))
    ELSE
      Wr.PutText(wr, Fmt.Pad(Fmt.LongReal(tr.mc, Fmt.Style.Fix, 4), 8));
      Wr.PutText(wr, " * dt ^");
      Wr.PutText(wr, Fmt.Pad(Fmt.LongReal(tr.mk, Fmt.Style.Fix, 2), 5));
    END;
    Wr.PutText(wr, "\n");
    Wr.Flush(wr);
  END WriteTestResultsLine;

PROCEDURE WriteTestResultsHeader(wr: Wr.T) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(stdout, Fmt.Pad("TEST", 15, align := Fmt.Align.Left));
    Wr.PutText(wr, " ");
    Wr.PutText(wr, Fmt.Pad("nSteps", 7));
    Wr.PutText(wr, " ");
    Wr.PutText(wr, Fmt.Pad("nRedos", 7));
    Wr.PutText(wr, " ");
    Wr.PutText(wr, Fmt.Pad("nEvals", 7));
    Wr.PutText(wr, " ");
    Wr.PutText(wr, " ");
    Wr.PutText(wr, Fmt.Pad("estimated error", 20));
    Wr.PutText(wr, " ");
    Wr.PutText(wr, Fmt.Pad("actual error", 20));
    Wr.PutText(wr, "\n") 
  END WriteTestResultsHeader;

PROCEDURE WriteAdaptiveTestPlotLine(
    wr: Wr.T; 
    t: Time; 
    READONLY s: State; 
    READONLY v: Velocity;
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(wr, " ");
    Wr.PutText(wr, Fmt.Pad(Fmt.LongReal(t, Fmt.Style.Fix, 4), 12));
    Wr.PutText(wr, " ");
    FOR i := 0 TO LAST(s) DO 
      Wr.PutText(wr, " ");
      Wr.PutText(wr, Fmt.Pad(Fmt.LongReal(s[i], Fmt.Style.Sci, 3), 12));
    END;
    FOR i := 0 TO LAST(s) DO 
      Wr.PutText(wr, " ");
      Wr.PutText(wr, Fmt.Pad(Fmt.LongReal(v[i], Fmt.Style.Sci, 3), 12));
    END;
    Wr.PutText(wr, "\n");
  END WriteAdaptiveTestPlotLine;

PROCEDURE WriteAdaptiveTestPlotHeader(
    wr: Wr.T; 
    iName: TEXT; 
    dt, dtMin, dtMax, tol: LONGREAL;
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(wr, "#");
    Wr.PutText(wr, " integrator = ");
    Wr.PutText(wr, iName);
    Wr.PutText(wr, "\n");
    Wr.PutText(wr, "#");
    Wr.PutText(wr, " dt = ");
    Wr.PutText(wr, Fmt.LongReal(dt,    Fmt.Style.Sci, 4));
    Wr.PutText(wr, " dtMin = ");                        
    Wr.PutText(wr, Fmt.LongReal(dtMin, Fmt.Style.Sci, 4));
    Wr.PutText(wr, " dtMax = ");                        
    Wr.PutText(wr, Fmt.LongReal(dtMax, Fmt.Style.Sci, 4));
    Wr.PutText(wr, " tol = ");                          
    Wr.PutText(wr, Fmt.LongReal(tol,   Fmt.Style.Sci, 4));
    Wr.PutText(wr, "\n");
  END WriteAdaptiveTestPlotHeader;

PROCEDURE WriteStepTestPlotLine(
    wr: Wr.T; 
    scale: LONGREAL;
    READONLY tr: ARRAY OF TestResults;
    READONLY trRef: ARRAY OF TestResults;
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(wr, " ");
    Wr.PutText(wr, Fmt.Pad(Fmt.LongReal(scale, Fmt.Style.Fix, 8), 12));
    Wr.PutText(wr, " ");
    FOR i := 0 TO LAST(tr) DO 
      Wr.PutText(wr, " ");
      WITH
        n = tr[i].nEvals,
        nRef = trRef[i].nEvals,
        r = FLOAT(n, LONGREAL)/FLOAT(nRef, LONGREAL)
      DO
        Wr.PutText(wr, Fmt.Pad(Fmt.LongReal(r, style := Fmt.Style.Fix, prec := 4), 8))
      END
    END;
    Wr.PutText(wr, "\n");
  END WriteStepTestPlotLine;

PROCEDURE WriteStepTestPlotHeader(
    wr: Wr.T; 
    iName: TEXT;
    READONLY p: ARRAY OF Problem;
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(wr, "#");
    Wr.PutText(wr, " integrator = ");
    Wr.PutText(wr, iName);
    Wr.PutText(wr, "\n");
    Wr.PutText(wr, "#");
    Wr.PutText(wr, Fmt.Pad("scale", 12));
    Wr.PutText(wr, " ");
    FOR i := 0 TO LAST(p) DO
      Wr.PutText(wr, " ");
      Wr.PutText(wr, Fmt.Pad(p[i].tag, 8));
    END;
    Wr.PutText(wr, "\n");
  END WriteStepTestPlotHeader;

PROCEDURE StatsGather(VAR s: Stats; x, y: LONGREAL) =
  (* Accumulates statistics for least-squares linear fit "y = a*x + b" *)
  BEGIN
    s.xx := s.xx + x*x; s.xy := s.xy + x*y;
    s.xu := s.xu + x;   s.yu := s.yu + y;
    s.uu := s.uu + 1.0d0
  END StatsGather;

PROCEDURE StatsFit(READONLY s: Stats; VAR a, b: LONGREAL) =
  (* performs least-squares linear fit "y = a*x + b". *)
  BEGIN
    WITH
      det = MAX(s.xx*s.uu - s.xu*s.xu, 1.0d-200)
    DO
      a := (s.xy*s.uu - s.xu*s.yu)/det;
      b := (s.xx*s.yu - s.xu*s.xy)/det
    END
  END StatsFit;

PROCEDURE Trace(c: CHAR) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutChar(stderr, c)
  END Trace;

PROCEDURE WriteEOL(wr: Wr.T) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutChar(wr, '\n')
  END WriteEOL;

BEGIN
  DoIt()
END TestIntegrators.
