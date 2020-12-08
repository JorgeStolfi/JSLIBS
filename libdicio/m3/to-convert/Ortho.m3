MODULE Ortho EXPORTS Main;

(* IC-UNICAMP "Ortho" spell checker and advisor. *)
(* See the copyright and disclaimer notice at the end of this file. *)

(* 
  This program implements spell checkers using the automaton-based
  techniques developed at IC/Unicamp by the "dicio" project.

  The program has three operating modes: 

    "-batch"  

      Reads a list of words, looks them up in the vocabulary.
      and prints the ones not found.  Optionally, suggests
      correct alternatives (for all words, or only the wrong ones).

    "-email" 

      Similar to "-batch", but assumes that the input has been
      submitted through the "ortho" e-mail spelling service.

    "-interactive" 

      Assumes the input and output are connected to the "OpenWin"
      interface written by R. Anido.

  The vocabulary is given by a basic wordlist "V0" (in the form of a
  compressed automaton) modified by auxiliary word lists "V1", "V2",...
  (plain text format, one word per line).

  Each auxiliary word list may be either ``additive'' (words that
  should be treated as correct) or ``subtractive'' (words that
  should be treated as incorrect).  These amends are applied in
  sequence to the basic vocabulary "V0".  The set of correct words
  is, therefore,

    V = ((V0 op1 V1) op2 V2) op3 V3) ...

  where "op1", "op2", ... are either set union or set difference, as
  specified by the user. 

  The "-advise" option will look for alternatives in the set "V"
  above.  "

  In "-interactive" mode, there are two additional auxiliary
  wordlists, "V+" and "V-" initially empty, respectively added and
  subtracted to the vocabulary (after all user-specified wordlists).
  In this mode, each line from "stdin" should begin with a
  one-character operation code, followed by one word, with the
  following convention:
  
    ' '  save this word for later checking
    '#'  check this word, plus all saved ones, and print result immediately
    '?'  list suggestions for this word.
    '+'  add this word to the "V+" list (and remove it from "V-").
    '-'  add this word to the "V-" list (and remove it from "V+").
    
  Error messages that are due to incorrect usage (command line
  errors, missing or malformed vocabularies, bad characters 
  on input, etc.) are always printed to "stderr", and cause 
  the program to abort with status = 1. 
  
  If "-logDir" is specified, a summary of the run is printed
  to a log file in that directory.  The file name will be the 
  date and time of the run.  The file will contain a
  printout of the command line options, the user and group id
  of the caller, the machine name, and a list of the words
  submitted to "stdin" that were not found in the vocabulary.

*)

IMPORT Wr, Rd, Fmt, Text, FileStream, Scan, 
       ParseParams, Thread,
       ReadOnlyPermDAG, Util, Adviser, AVLTextTree,
       M3toC, Upwd, Uugid;

FROM Stdio IMPORT stdin,stdout,stderr;
FROM Basics IMPORT Abort;
FROM Util IMPORT BadVersion;

CONST
  Version = "Ortho.m3 version 3.0 (97-01-19)";

CONST
  HelpText =
    "options:\n" &
    "  Ortho \\\n" &
    "  [ -help | -version | -batch | -email SENDER | -interactive ] \\\n" &
    "  -use FILE [ plus FILE | minus FILE ]... \\\n" &
    "  [ -advise [ -all ] [ -typos ] ] \\\n" &
    "  [ -logDir DIRECTORY ]\n";

CONST
  MaxAux = 4;

TYPE
  Options = RECORD
      logDir: TEXT;                      (* Log file directory, or "" *)
      mode: Mode;                        (* Processing mode *)
      sender: TEXT;                      (* Email address of sender, if "-mail" *)
      voc: TEXT;                         (* Filename of automaton *)
      aux: ARRAY [0..MaxAux-1] OF TEXT;  (* Filenames of auxiliary wordlists *)
      op: ARRAY [0..MaxAux-1] OF SetOp;  (* "op[k]" tells "aux[k]" is "+" or "-" *)
      advise: BOOLEAN;                   (* TRUE to suggest alternatives *)
      adviseAll: BOOLEAN;                (* TRUE to advise even on correct wds *)
      adviseTypos: BOOLEAN;              (* TRUE to try correcting typos too *)
      cmd: TEXT;                         (* The command line options, formatted *)    
    END;

  SetOp = { Plus, Minus, None };

TYPE
  Vocabulary = RECORD
      voc: ReadOnlyPermDAG.T;
      aux: ARRAY [0..MaxAux-1] OF AVLTextTree.T;
      op: ARRAY [0..MaxAux-1] OF SetOp;
    END;

<* FATAL Thread.Alerted, Wr.Failure *>

PROCEDURE DoIt() =
  BEGIN
    WITH
      start = Util.FmtDate(Util.GetDate()),
      o = GetCommandLineArguments(),
      log = OpenLog(o.logDir),
      v = GetVocabulary(o.ignore, o.voc, o.auxFile, o.op)
    DO
      TRY
        Initialize(log, start, o.cmd);
        CASE o.mode OF
        | Mode.Batch, Mode.Mail => 
            CheckBatch(log, v, o....);
        | Mode.Interactive =>
            CheckInteractive(log, v, o, ...)
        END;
        Finalize(log)
      FINALLY
        Wr.Close(stdout);
        Wr.Flush(stderr);
        IF log # NIL THEN Wr.Close(log) END;
      END
    END
  END DoIt;

PROCEDURE GetCommandLineArguments() =
  VAR o: Options;
      auxKey: TEXT;

  PROCEDURE NoStdin(pp: ParseParams.T; rd: Rd.T) RAISES {ParseParams.Error} =
    BEGIN
      IF rd = stdin THEN pp.error("cannot use \"-\" as wordlist") END;
    END NoStdin;

  BEGIN
    o.user := Util.GetUserId();
    o.start := Util.GetDate();
    WITH
      pp = ParseParams.New(ParseParams.T).init(stderr),
      cmdWr = NEW(TextWr.T).init()
    DO
      TRY
        IF pp.keywordPresent("-help") THEN
          pp.finish();
          Wr.PutText(stderr, HelpText);
          Wr.Flush(stderr);
          Process.Exit(0)
        ELSIF pp.keywordPresent("-version") THEN
          pp.finish();
          Wr.PutText(stderr, Version & "\n");
          Wr.Flush(stderr);
          Process.Exit(0)
        ELSIF pp.keywordPresent("-batch") THEN
          o.mode := Mode.Batch;
          o.sender := "";
          ParamUtil.PrintBoolArg(cmdWr, "-batch")
        ELSIF pp.keywordPresent("-mail") THEN
          o.mode := Mode.Mail;
          o.sender := pp.getNext();
          IF Text.Empty(o.sender) 
          OR Text.GetChar(o.sender, 0) = '-' THEN
            pp.error("bad sender address \"" & o.sender & "\"")
          END;
          ParamUtil.PrintTextArg(cmdWr, "-mail", o.sender)
        ELSIF pp.keywordPresent("-interactive") THEN
          o.mode := Mode.Interactive;
          o.sender := "";
          ParamUtil.PrintBoolArg(cmdWr, "-interactive")
        ELSE
          o.mode := Mode.Batch;
          o.sender := "";
          ParamUtil.PrintBoolArg(cmdWr, "-batch")
        END;

        o.ignore := ParamUtil.GetFileName(pp, cmdWr, "-ignore");
        NoStdin(pp, o.ignore);
        o.voc := ParamUtil.GetFileName(pp, cmdWr, "-use");
        NoStdin(pp, o.voc);
        IF o.voc = NIL THEN pp.error("\"-use\" is required") END;
        FOR k := 0 TO MaxAux-1 DO
          IF pp.textNext("plus") THEN
            o.op[k] := SetOp.Plus;  auxKey := "plus"
          ELSIF pp.textNext("minus") THEN
            o.op[k] := SetOp.Minus; auxKey := "minus"
          ELSE
            o.op[k] := SetOp.None
          END;
          IF o.op[k] = SetOp.None THEN
            o.auxFile[k] := ""
          ELSE
            WITH name = pp.getNext() DO
              EVAL ParamUtil.CheckFileName(name);
              o.aux[k] := name; NoStdin(pp, name);
              ParamUtil.PrintTextArg(cmdWr, "  " & auxKey, name)
            END;
          END
        END;

        o.advise := ParamUtil.GetBool(pp, cmdWr, "-advise");
        IF o.advise THEN
          o.adviseAll := ParamUtil.GetBool(pp, cmdWr, "-all");
          o.adviseTypos := ParamUtil.GetBool(pp, cmdWr, "-typos");
        ELSE
          o.adviseAll := FALSE;
          o.adviseTypos := FALSE
        END;
        pp.finish();
        o.cmd := TextWr.ToText(cmdWr);
      EXCEPT
        ParseParams.Error  => 
          Wr.PutText(stderr, HelpText);
          Wr.Flush(stderr);
          Process.Exit(1)
      END
    END
  END GetCommandLineArguments;

PROCEDURE GetVocabulary(
    ignore: TEXT; 
    voc: TEXT; 
    READONLY aux: ARRAY OF TEXT;
    op: ARRAY OF SetOp;
  ): Vocabulary =
  VAR v: Vocabulary;
  BEGIN
    v.ignore := ReadWordList(ignore);
    v.voc := ReadAutomaton(voc);
    FOR k := 0 TO MaxAux-1 DO
      IF NOT Text.Empty(aux[k]) THEN
        v.aux[k] := ReadWordList(aux[k]);
        <* ASSERT op[k] # SetOp.None *>
      ELSE
        v.aux[k] := NIL;
        <* ASSERT op[k] = SetOp.None *>
      END;
      v.op[k] := op[k]
    END;
  END GetVocabulary;
  
PROCEDURE ReadAutomaton(name: TEXT): ReadOnlyPermDAG.T =
  BEGIN
    WITH rd = Util.OpenRd(name, stderr) DO
      RETURN ReadOnlyPermDAG.LoadCompr(rd)
    END
  END ReadAutomaton;

PROCEDURE ReadWordList(name: TEXT): AVLTextTree.T =
  VAR buf: REF String;
  VAR nLines: CARDINAL := 0;
  BEGIN
    WITH 
      rd = Util.OpenRead(name, stderr),
      t = NEW(...)
    DO
      
      TRY
        LOOP
          WITH s = ReadString(rd, buf, name, nLines+1) DO
            INC(nLines);
            IF NUMBER(s^) > 0 THEN t.Insert(s) END
          END
        END
      EXCEPT
      | Done  => Rd.Close(rd);
      RETURN t
    END
  END ReadWordList;

PROCEDURE ReadString(
    rd: Rd.T; 
    VAR buf: REF String; 
    file: TEXT; 
    line: CARDINAL;
  ): REF String 
  RAISES {Done} =
  VAR len: CARDINAL;
  BEGIN
    TRY
      encoding.ReadString(rd, buf, len)
    EXCEPT
    | Encoding.BadChars => 
        Wr.PutText(stderr, file & ", line " & Fmt.Int(line) & ": invalid character\n");
        Wr.Flush(stderr);
        Process.Exit(1)          
    END;
    WITH s = NEW(REF String, len) DO
      s^ := SUBARRAY(s, 0, len);
      RETURN s
    END
  END ReadString;

PROCEDURE CheckBatch() RAISES{Terminate} =
  VAR buf: REF String;
  VAR len: CARDINAL;
  VAR nLines: CARDINAL := 0;
 
  PROCEDURE CheckStringAction(READONLY s: String) RAISES {Abort}  =

    PROCEDURE PrintAltAction(READONLY a: String) RAISES {Abort} =
      BEGIN
        Wr.PutChar(stdout, ' ')
        Wr.PutChar(stdout, ' ')
        encoding.PrintString(stdout, a);
        Wr.PutChar(stdout, '\n')
      END PrintAltAction;

    BEGIN
      TRY
        encoding.PrintString(stdout, s);
        Wr.PutChar(stdout, '\n')
        IF log # NIL THEN 
          Wr.PutChar(log, '!')
          Wr.PutChar(log, ' ')
          encoding.PrintString(log, s)
          Wr.PutChar(log, '\n')
        END;
        IF advise THEN
          WITH t = Adviser.Alternatives(s, v.voc, v.aux, v.op, adviseTypos) DO
            t.Enum(PrintAltAction);
          END;
        END;
      END
    END CheckStringAction;

  BEGIN
    TRY
      LOOP
        WITH s = ReadString(stdin, buf, "stdin", nLines+1) DO
          INC(nLines);
          IF NOT IsBogus(v, s^) THEN
            IF o.adviseAll OR NOT IsCorrect(v, s^) THEN
              e.Insert(s)
            END
          END
        END
      END
    EXCEPT
    | Done => 
        e.Enum(CheckStringAction);
        Wr.Flush(stdout);
    | Encoding.BadChars =>
        Wr.PutText(stderr, 
          "stdin, line " & Fmt.Int(nLines + 1) & ": invalid character\n"
        );
        Wr.Flush(stderr);
        Process.Exit(1)          
    END
  END CheckBatch;

PROCEDURE CheckInteractive() RAISES {Terminate} = 
  VAR buf: REF String;
  VAR len: CARDINAL;
  VAR nLines: CARDINAL := 0;

  PROCEDURE NextWord(): REF String =
    VAR s: REF String;
    BEGIN
      TRY
        s := ReadString(stdin, buf, ""stdin", nLines+1)
      EXCEPT 
      | Done => 
          Wr.PutText(stderr, 
            "stdin, line " & Fmt.Int(nLines + 1) & ": unexpected end-of-file\n"
          );
        Wr.Flush(stderr);
        Process.Exit(1)          
      END;
      INC(nLines);
      RETURN s
    END NextWord;

  PROCEDURE Save() =
    BEGIN
      WITH s = NextWord() DO
        IF NOT IsCorrect(v, s^) THEN
          e.Insert(s)
        END
      END
    END Save;

  PROCEDURE Check() =
    
    PROCEDURE PrintStringAction(READONLY s: String) RAISES {Abort}  =
      BEGIN
        encoding.PrintString(stdout, a);
        Wr.PutChar(stdout,'\n')
        IF log # NIL THEN 
          Wr.PutChar(log, '!')
          Wr.PutChar(log, ' ')
          encoding.PrintString(log, s)
          Wr.PutChar(log, '\n')
        END;
      END PrintStringAction;
    
    BEGIN
      Save();
      e.Enum(PrintStringAction);
      Wr.PutText(stdout, "#\n")
      Wr.Flush(stdout);
      e.DeleteAll();
    END Check;
    
  PROCEDURE Advise() =

    PROCEDURE PrintAltAction(READONLY a: String) RAISES {Abort} =
      <* FATAL Wr.Failure *>    
      BEGIN
        encoding.PrintString(stdout, a);
        Wr.PutChar(stdout, '\n')
      END PrintAltAction;

    BEGIN
      <* ASSERT e.IsEmpty() *>
      WITH s = NextWord() DO
        WITH t = Adviser.Alternatives(s, v.voc, v.aux, v.op, adviseTypos) DO
          t.Enum(PrintAltAction);
        END;
      END;
      Wr.PutText(stdout, "#\n")
      Wr.Flush(stdout);
    END Advise;

  PROCEDURE Add() =
    BEGIN
      WITH s = NextWord() DO
        v.add.Insert(s)
        v.sub.Delete(s)
      END;
    END Add;

  PROCEDURE Sub() =
    BEGIN
      WITH s = NextWord() DO
        v.add.Delete(s);
        v.sub.Insert(s)
      END;
    END Sub;

  BEGIN
    TRY
      LOOP
        c := Rd.getChar(stdin); (* OK if end-of-file here *)
        INC(nLines);
        TRY 
          IF c = '#' THEN
            Check()
          ELSIF c = '*' THEN
            Save()
          ELSIF c = '=' THEN
            Advise()
          ELSIF c = '+' THEN
            Add()
          ELSIF c = '-' THEN
            Sub()
          ELSE
            Wr.PutText(stderr, 
              "stdin, line " & Fmt.INT(nLines) & ": invalid op code\n"
            );
            Wr.Flush(stderr);
            Process.EXIT(1);
          END;
        EXCEPT 
        END
      END
    EXCEPT
    | Rd.EndOfFile => (* OK *)
    END
  END CheckInteractive;

PROCEDURE Initialize(log: Wr.T; start: TEXT; cmd: TEXT) =
  BEGIN
    IF log # NIL THEN
      WITH
        user  = Util.GetUserId(),
        group = Util.GetGroupId(),
        host  = Util.GetHostName()
      DO
        Wr.PutText(log, o.start); 
        Wr.PutChar(log, '\n')
        
        Wr.PutText(log, Version);  
        Wr.PutChar(log, '\n')
        
        Wr.PutText(log, user); 
        Wr.PutChar(log, '.')
        Wr.PutText(log, group); 
        Wr.PutChar(log, '@')
        Wr.PutText(log, host); 
        Wr.PutChar(log, '\n')
        
        Wr.PutText(log, o.cmd);
        Wr.PutChar(log, '\n')
      END;
    END;

PROCEDURE Finalize(log: Wr.T) =
  BEGIN
    IF log # NIL THEN
      WITH
        cpuTime = Util.GetExecTimesText(),
        stop = Util.FmtDate(Util.GetDate())
      DO
        Wr.PutText(log, o.stop);
        Wr.PutChar(log, '\n')
      END
    END
  END Finalize;

PROCEDURE PrintTree(wr:Wr.T; t: AVLTextTree.T)  =
  <* FATAL Abort *>

  PROCEDURE Action(READONLY s: String) =
    <* FATAL Wr.Failure *>
    BEGIN
      encoding.PrintString(wr, s);
      Wr.PutText(wr, "\n");
    END Action;

  BEGIN
    t.Enum(Action);
    Wr.Flush(wr)
  END PrintTree;

BEGIN 
  DoIt()
END Ortho.

(*
  HISTORY:

  92-08-25 (tk) Modified for compatibility with version 2.04 of 
                Modula-3; added new options: -transpose, -remove, 
                -subst, -insert

  92-09-30 (tk) Small modifications when recompiled under 2.07

  92-10-23 (tk) Inclusion OF the "-remote" option

  93-05-06 (tk) Intallation path for DCC changed to /proj/dicio

  97-01-18 (js) General rewrite for Modula-3 release 3.5.3 and 
                libm3reduced-2.  Removed all security options,
                redesigned options and user interface protocol.

*)

(****************************************************************************)
(* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           *)
(*                    Campinas, SP, Brazil                                  *)
(*                                                                          *)
(* Authors:                                                                 *)
(*                                                                          *)
(*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.ansp.br>    *)
(*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.ansp.br>  *)
(*   Jorge Stolfi        - DEC Systems Research Center <stolfi@src.dec.com> *)
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

