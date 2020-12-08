(****************************************************************************)
(* (C) Copyright 1994 Univrsidade Estadual de Campinas (UNICAMP)           *)
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

MODULE Example EXPORTS Main;          (* Version 1.0 *)

(*
  This program writes into standard output file all the words
  contained in a automaton, in lexicographic order.

*)

(* HISTORY:

  03.Dec.94: Version 1.0 (simplified version of MaintainAutomaton
  Version 2.1) [tk]

  06.Aug.95: Converted release 3.5.3 of Modula-3.  Rewrote command
  line parsing to use ParamUtil.  Changed formatting. style, and other
  minor things.

*)

IMPORT ParseParams, Process, Wr, Rd, Text, TextWr, Fmt, Thread;
IMPORT Reduced, ParamUtil, Util, Encoding, PlainEncoding;
FROM Stdio IMPORT stdout,stderr;
FROM Text IMPORT Empty;
FROM Basics IMPORT NAT, Abort, String;

CONST
  Help =
    "Example \\\n" &
    "  -load <file>        : file with input automaton\\\n" &
    "  [ -log <file> ]     : log file, for messages and such [stderr]\\\n" &
    "  [ -help ]           : prints this message\n";

CONST
  Version = "3.5.3-1";

VAR

  arg: RECORD
      log: Wr.T;     (* Log file *)
      load: Rd.T;    (* File with input automaton *)
      cmd: TEXT;     (* The formatted command line, for documentation *)
    END;

PROCEDURE Main() =
  <* FATAL Wr.Failure, Abort, Thread.Alerted *>
  VAR lineCount: NAT := 0;
      aut: Reduced.T;
  BEGIN
    GetCommandLineArgs();

    StartedMessage("Example " & Version & " " & Util.FmtDate(Util.GetDate()));

    NL(); Wr.PutText(arg.log, arg.cmd); NL();

    StartedMessage("loading the automaton");
    aut := Reduced.Load(arg.load);
    FinishedMessage("loading the automaton");

    StartedMessage("spelling out all automaton words");
    WITH
      root = aut.Root(),
      encoding = PlainEncoding.New()
    DO

      PROCEDURE EnumAction(READONLY w: String) RAISES {} =
        <* FATAL Wr.Failure, Encoding.BadString *>
        BEGIN
          Wr.PutText(stdout, encoding.StringToText(w));
          Wr.PutChar(stdout, '\n');
          INC(lineCount);
        END EnumAction;

      BEGIN
        aut.EnumSuffs(root, EnumAction)
      END

    END;

    Wr.Flush(stdout);
    Wr.PutText(arg.log, "spelled " & Fmt.Int(lineCount) & " words\n");
    FinishedMessage("spelling out all automaton words");

    FinishedMessage("Example " & Version);
    Wr.Close(arg.log)
  END Main;

PROCEDURE GetCommandLineArgs() =
  (*
    Parses the command line parameters and stores their values in "arg".
  *)
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    TRY
      WITH
        pp = NEW(ParseParams.T).init(stderr),
        cmdWr = NEW(TextWr.T).init()
      DO
        IF pp.keywordPresent("-help") THEN
          Wr.PutText(stderr, "Usage: \n" & Help);
          Process.Exit(0);
        END;
        
        Wr.PutText(cmdWr, "Example");
        
        arg.log := ParamUtil.GetLogWr(pp, cmdWr);
        
        arg.load := ParamUtil.GetRd(pp, cmdWr, "-load");
        IF arg.load = NIL THEN
          pp.error("\"-load\" argument is required")
        END;

        pp.finish();

        Wr.PutText(cmdWr, "\n");
        arg.cmd := TextWr.ToText(cmdWr);
      END
    EXCEPT
      ParseParams.Error  =>
         Wr.PutText(stderr, "Usage:\n" & Help);
         Process.Exit(1)
    END;
  END GetCommandLineArgs;

PROCEDURE StartedMessage(m: TEXT)  =
  (*
    Prints "started m" to "arg.log", with separating line. *)
  BEGIN
    NL();
    Message("started " & m);
    EVAL Util.GetExecTimes(); (* Reset the clock for FinishedMessage *)
  END StartedMessage;

PROCEDURE FinishedMessage(m: TEXT)  =
  (*
    Prints "finished m" to "arg.log", with CPU time and separating line. *)
  BEGIN
    WITH tim = Util.GetExecTimesText().total DO
      Message("finished " & m & " (" & tim & ")")
    END;
    NL();
  END FinishedMessage;

CONST MessageWidth = 76;

PROCEDURE Message(m: TEXT)  =
  (*
    Prints a separating line to "arg.log". *)
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(arg.log,"---");
    IF Empty(m) THEN
      Wr.PutText(arg.log, "--")
    ELSE
      Wr.PutText(arg.log, " ");
      Wr.PutText(arg.log, m);
      Wr.PutText(arg.log, " ")
    END;
    FOR i:=1 TO MessageWidth - 5 - Text.Length(m) DO 
      Wr.PutChar(arg.log, '-')
    END;
    NL();
    Wr.Flush(arg.log)
  END Message;

PROCEDURE NL()  =
  (*
    Prints a blank line to "arg.log". *)
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(arg.log,"\n");
  END NL;

BEGIN
  Main()
END Example.



