MODULE Util;

(* See the copyright and disclaimer note at the end of this file. *)

IMPORT Date, Time, Text, Fmt, Utime, ISOLatin1;
IMPORT Rd, FileRd, TextRd, Wr, FileWr, TextWr, OSError, Process, Thread;
FROM Stdio IMPORT stdin, stdout;

PROCEDURE GetDate(): Date.T =
  BEGIN
    RETURN Date.FromTime(Time.Now())
  END GetDate;
  
PROCEDURE GetUTCDate(): Date.T =
  BEGIN
    RETURN Date.FromTime(Time.Now(), z := Date.UTC)
  END GetUTCDate;

PROCEDURE ElapsedTime(READONLY d1, d2: Date.T): Time.T =
  <* FATAL Date.Error *>
  BEGIN
    RETURN Date.ToTime(d2) - Date.ToTime(d1)
  END ElapsedTime;
  
PROCEDURE FmtDate(READONLY d: Date.T): TEXT =
  CONST
    WeekDayName = ARRAY Date.WeekDay OF TEXT {
      "Sun","Mon","Tue","Wed","Thu", "Fri","Sat"
    };
  BEGIN
    WITH
      month = ORD(d.month)+1
    DO
      RETURN 
        WeekDayName[d.weekDay] & " " & 
        Fmt.Pad(Fmt.Int(d.year),   4) & "-" &
        Fmt.Pad(Fmt.Int(month),    2, '0') & "-" &
        Fmt.Pad(Fmt.Int(d.day),    2, '0') & " " &
        Fmt.Pad(Fmt.Int(d.hour),   2, '0') & ":" &
        Fmt.Pad(Fmt.Int(d.minute), 2, '0') & ":" &
        Fmt.Pad(Fmt.Int(d.second), 2, '0') & " " &
        d.zone & " (" &
        FmtZoneOffset(d.offset) & ")"
    END
  END FmtDate;
  
PROCEDURE FmtZoneOffset(secs: INTEGER): TEXT =
(*
  Formats a timezone offset (given in seconds) as "shh:mm",
  where "s" is the sign. *)
  VAR sign: TEXT;
  BEGIN
    IF secs < 0 THEN sign := "-" ELSE sign := "+" END;
    RETURN 
      sign &
      Fmt.Pad(Fmt.Int(secs DIV 3600), 2, '0') & ":" &
      Fmt.Pad(Fmt.Int((secs DIV 60) MOD 60), 2, '0')
  END FmtZoneOffset;
  
PROCEDURE FmtDateNum(READONLY d: Date.T): TEXT =
  BEGIN
    WITH
      month = ORD(d.month)+1
    DO
      RETURN 
        Fmt.Pad(Fmt.Int(d.year),        4, '0') & "-" &
        Fmt.Pad(Fmt.Int(month),         2, '0') & "-" &
        Fmt.Pad(Fmt.Int(d.day),         2, '0') & "-" &
        Fmt.Pad(Fmt.Int(d.hour),        2, '0') &
        Fmt.Pad(Fmt.Int(d.minute),      2, '0') &
        Fmt.Pad(Fmt.Int(d.second),      2, '0')
    END
  END FmtDateNum;
  
PROCEDURE FmtTime(t: Time.T): TEXT =
  VAR d: TEXT;
  BEGIN
    WITH
      secs = ROUND(t),
      mins = secs DIV 60,
      hurs = mins DIV 60,
      days = hurs DIV 24
    DO
      IF days = 0 THEN
        d := ""
      ELSE
        d := Fmt.Pad(Fmt.Int(days), 3) & "+"
      END;
      RETURN 
        d &
        Fmt.Pad(Fmt.Int(hurs MOD 24), 2, '0') & ":" &
        Fmt.Pad(Fmt.Int(mins MOD 60), 2, '0') & ":" &
        Fmt.Pad(Fmt.Int(secs MOD 60), 2, '0');
    END
  END FmtTime;

VAR
  StartingTimes := NEW(Utime.struct_tms_star);
  CurrentTimes  := NEW(Utime.struct_tms_star);

PROCEDURE GetExecTimes(): ExecTimes =
  VAR res: ExecTimes;
  BEGIN
    WITH retCode = Utime.times(CurrentTimes) DO
      <* ASSERT retCode = 0 *>
    END;
    WITH
      cu = CurrentTimes^.tms_utime,
      su = StartingTimes^.tms_utime,
      cs = CurrentTimes^.tms_stime,
      ss = StartingTimes^.tms_stime,
      u = cu - su,
      s = cs - ss
    DO
      res.user := FLOAT(u, LONGREAL);
      res.system := FLOAT(s, LONGREAL);
      res.total := FLOAT(u+s, LONGREAL);
      su := cu;
      ss := cs;
      RETURN res
    END
  END GetExecTimes;

PROCEDURE FmtExecTimes(READONLY t: ExecTimes): ExecTimesText =
  BEGIN
    RETURN ExecTimesText{
      user   := FmtTime(t.user),
      system := FmtTime(t.system),
      total  := FmtTime(t.total)
    }
  END FmtExecTimes;

PROCEDURE GetExecTimesText(): ExecTimesText =
  BEGIN
    RETURN FmtExecTimes(GetExecTimes())
  END GetExecTimesText;
  
PROCEDURE OpenRd(name: TEXT; err: Wr.T): Rd.T =
  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    IF Text.Empty(name) THEN
      RETURN NIL
    ELSIF Text.Equal(name, "-") THEN
      RETURN stdin
    ELSE
      TRY
        WITH rd = FileRd.Open(name) DO
          RETURN rd
        END
      EXCEPT
      | OSError.E =>
          Wr.PutText(err, "*** Util.OpenRd: error opening file = \"");
          Wr.PutText(err, name);
          Wr.PutText(err, "\"\n");
          Process.Exit(1);
          <* ASSERT FALSE *>
      END
    END
  END OpenRd;

PROCEDURE OpenWr(name: TEXT; err: Wr.T): Wr.T =
  <* FATAL Thread.Alerted, Wr.Failure *>
  BEGIN
    IF Text.Empty(name) THEN
      RETURN NIL
    ELSIF Text.Equal(name, "-") THEN
      RETURN stdout
    ELSE
      TRY
        WITH wr = FileWr.Open(name) DO
          RETURN wr
        END
      EXCEPT
      | OSError.E =>
          Wr.PutText(err, "*** Util.OpenWr: error opening file = \"");
          Wr.PutText(err, name);
          Wr.PutText(err, "\"\n");
          Process.Exit(1);
          <* ASSERT FALSE *>
      END
    END
  END OpenWr;

PROCEDURE PrintDoc(wr: Wr.T; doc: TEXT; prefix: TEXT := "|") =

  VAR rd: Rd.T := TextRd.New(doc);
  
  PROCEDURE CopyLine() RAISES {Rd.EndOfFile} =
  (*
    Copy one line from "rd" to "wr", prefixed by "prefix". 
    Supplies a final '\n' if next line exists but does not end with newline.
    Raises Rd.EndOfFile if there are no more lines in "rd". *)
    
    <* FATAL Rd.Failure, Wr.Failure, Thread.Alerted *>
    VAR c: CHAR;
    BEGIN
      c := Rd.GetChar(rd); (* If EOF here, propagate to caller *)
      Wr.PutText(wr, prefix);
      Wr.PutChar(wr, c);
      WHILE c # '\n' DO
        TRY c := Rd.GetChar(rd) EXCEPT Rd.EndOfFile => c := '\n' END;
        Wr.PutChar(wr, c)
      END
    END CopyLine;

  BEGIN
    TRY LOOP CopyLine() END EXCEPT Rd.EndOfFile => (* Ok *) END;
  END PrintDoc;

PROCEDURE ToLowerCase(t: TEXT): TEXT =
  <* FATAL Rd.EndOfFile, Rd.Failure, Wr.Failure, Thread.Alerted *>
  BEGIN
    WITH wr = NEW(TextWr.T).init(), rd = NEW(TextRd.T).init(t) DO
      WHILE NOT Rd.EOF(rd) DO Wr.PutChar(wr, ISOLatin1.Lower[Rd.GetChar(rd)]) END;
      RETURN TextWr.ToText(wr)
    END
  END ToLowerCase;

BEGIN
  WITH retCode = Utime.times(StartingTimes) DO
    <* ASSERT retCode = 0 *>
  END
END Util.

(****************************************************************************)
(* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           *)
(*                    Campinas, SP, Brazil                                  *)
(*                                                                          *)
(* Authors:                                                                 *)
(*                                                                          *)
(*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.br>         *)
(*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.br>       *)
(*   Jorge Stolfi        - CS Dept, UNICAMP <stolfi@dcc.unicamp.br>         *)
(*                                                                          *)
(* This file can be freely distributed, modified, and used for any          *)
(*   non-commercial purpose, provided that this copyright and authorship    *)
(*   notice be included in any copy or derived version of this file.        *)
(*                                                                          *)
(* DISCLAIMER: This software is offered ``as is'', without any guarantee    *)
(*   as to fitness for any particular purpose.  Neither the copyright       *)
(*   holder nor the authors or their employers can be held responsible for  *)
(*   any damages that may result from its use.                              *)
(****************************************************************************)
