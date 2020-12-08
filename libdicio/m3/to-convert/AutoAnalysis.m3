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

MODULE AutoAnalysis EXPORTS Main;

(* 
  The purpose of this program is to analyze an automaton to
  identify possible errors and/or reveal the structure of the lexicon.
  
  The following command line argument is required:
  
    -load FILENAME
    
      Read the automaton from FILENAME, in .dmp format, as created by
      MaintainAutomaton.
      
  Each analyis on that automaton is specified in the command line by
  an operation code, and additional arguments specific to each
  operation, as follows:

    -prefs FILENAME
    -suffs FILENAME
    -prefSuffs FILENAME

      Write into the designated file a complete listing of the prefix
      and/or sufix sets for the states whose numbers are specified in
      the standard input. These three commands are mutually exclusive.
      Additional arguments are
      
        -maxPrefSize NAT
        -maxSuffSize NAT
        
          Maximum number of bytes of prefixes and suffixes to print 
          for each state. Default is infinity.
          
    -unprod FILENAME

      Identify and write to FILENAME the "unproductive" states: those
      which are used by too few words.  Modified by the parameters:
      
        -maxUnprod NAT
        
          A state is considered "unproductive" if it is used by 
          at most NAT words.  Default is 1.
          
        -unprodSugg FILENAME
        
          write to FILENAME a list of the words that use the "unproductive"
          states.
          
        -maxUnprodSugg NAT
        
          Truncate the list of unproductive words after NAT entries.
          Default is infinity.
          
   -prod FILENAME
    
     Lists on FILENAME the "productive states": proper states 
     for which the number of words is much greater than
     the numebr of prefixes and suffixes.
     
     Modified by:

       -minProductivity NAT 

         Minimum value of the state's "productivity", defined
         as (NP-1)(NS-1).  Default is 4.

        -maxProdStates NAT
        
          Truncate the list of productive states after NAT entries.
          Default is infinity.
          
        -maxPrefSize NAT 
        -maxSuffSize NAT
        
          When printing a productive state, truncate the prefix (resp. suffix)
          printout after NAT bytes.  Default is infinity.
          
    -similar FILENAME

      Identifies and reports in FILENAME all pairs of "similar" states
      (which have similar sets of suffixes).  The following auxiliary
      parameters control the definition of "similar states" and the
      size of the output:
      
        -maxSimDiff NAT
        
          Maximum acceptable difference (number of words) between the
          suffix sets of two similar states.  Default is 2.
          
        -maxRelSimDiff NAT
        
          Maximum acceptable difference (% relative to total suffixes)
          between the suffix sets of two similar states. Default is 100%.
          
        -minSimInter NAT
        
          Minimum acceptable number of common suffixes of two similar states.
          Default is 1.
          
        -maxSimPref NAT
        
          Only states with at most NAT prefixes are considered.
          Default is infinity.
          
        -maxSimGen NAT
        
          Considers only pairs of states that could be made equal 
          by adding at most NAT words (not suffixes).
          Default is infinity.
          
        -maxPrefSize NAT 
        -maxSuffSize NAT
        
          When printing a pair of similar states, truncate the 
          prefix (resp. suffix) list after NAT bytes. 
          Default is infinity.
          
        -simSugg FILENAME
        
          For each pair of similar states found, write to FILENAME the
          words that would need to be included in order to coalesce 
          those states.
          
        -maxSimSugg NAT
        
          When "-simSugg" is given, truncate the list of sugegstions for
          each pair after NAT words.  Default is no limit.
          
   -classes FILENAME
    
     Lists on FILENAME the "lexical classes": proper states which can be
     reached by at least a certain minimum number of prefixes from the
     root, with at least one prefix not going through another class.
     Modified by:
  
       -minClassPrefs NAT  

         Minimum number of prefixes for a state to be considered 
         a lexical class.  Default is 2.

       -radicals FILENAME

         Write to FILENAME a table of "lexical radicals" (prefixes of
         lexical classes which do not go through other lexical
         classes).  All words accepted by the automaton which do not
         go through a lexical class are listed as if they belonged to
         the "NullState" lexical class (which is never a lexical
         class, according to the definition above).

    -help

      Prints the help  message

  Other general arguments, applicable to all operations above, are:
  
    -log FILENAME
    
      Write all disgnostic messages and progress reports to FILENAME.
      (Default: stderr.)
  
  REMARK: "infinity" actually means "LAST(CARDINAL)".  For
  compatibility with previous versions, all keywords above may be
  given entierly in lower case.  Also, "-dump" and "-dumpIn" and
  "-dumpin" are equivalent to "-load"; and "-mess" is equivalent to 
  "-log".
*)


(* HISTORY:

  18.april.92: [TK] Modified for compatibility with version 2.04 of Modula-3 and the
               new interfaces

  7.oct.92:    [TK] New definition of lexical classes implemented. Improved listings
               separating true radicals from isolated ones. Prefix sets listings
               in alphabetical order (using "AddStringMaybeCrunch").

  21.oct.92:   [TK] Inclusion of "-singlePath" option.   

  11.nov.92:   [JS] Fixed bug in FindSimilarStates, that caused it to miss many 
               (most?) pairs of similar states.

  04.dec.92:   [JS] Added the "-minSimInter" parameter to provide 
               an alternative control on the minimum number
               of common suffixes required for two states to be 
               declared equivalent.

  04.dec.92:   [JS] Changed the default "-maxSimSugg" from 100 to "no limit".

  07.dec.92:   [TK] Fixed a typo in an output message.

  09.dec.92:   [TK] Included printing of sizes of suggested sets of words.

  10.dec.92:   [TK] Included "-maxSimGen"

  18.aug.95:   [JS] Converted to Modula-3 release 3.5.3.  Revised comments
                    and command line parsing.  Changed capitalizations 
                    to better match the standard Modula-3 style.
                    Changed the rounding semantics of "-maxRelSimDiff".

  27.sep.95:   [JS] Removed the "-singlePath" option, since it didn't seem
                    to have much success in the Portuguese tests
                    (see /proj/dicio/linguas/portugues/tti-analise/history.old).
                    Renamed the "-minPrefs" parameter to "-minClassPrefs".

  05.jul.97:   [JS] Added the "-prod" analysis.

*)

IMPORT 
  Rd, Wr, Text, Fmt, TextWr, Lex, ParseParams, 
  Thread, Process;

IMPORT 
  Basics, Reduced, ReducedPair, StringPrinter, Util, ParamUtil,
  IntQuickSort, Encoding, PlainEncoding;

FROM Basics IMPORT BOOL, NAT, Skip, Abort, String, RdNat, WrNat;
FROM Stdio IMPORT stderr,stdin;
FROM Reduced IMPORT State, UnitState, NullState;
FROM Text IMPORT Empty;

<* FATAL Wr.Failure, Rd.Failure, Thread.Alerted *>

CONST
  Version = "3.5.3-1";
  Help = 
    "AutoAnalysis \\\n" &
    "  -load FILE.dmp\\\n" &
    "  [ -prefs FILE.pre | -suffs FILE.suf | -prefSuffs FILE.psf ] \\\n" &
    "  [ -unprod FILE.unp \\\n" &
    "    [ -maxUnprod NAT ] \\\n" &
    "    [ -unprodSugg FILE.usg ] \\\n" &
    "    [ -maxUnprodSugg NAT ] \\\n" &
    "  ] \\\n" &
    "  [ -prod FILE.prd \\\n" &
    "    [ -minProductivity FILE.prd ] \\\n" &
    "    [ -maxProdStates NAT ] \\\n" &
    "  ] \\\n" &
    "  [ -similar FILE.sim \\\n" &
    "    [ -maxSimDiff NAT ] \\\n" &
    "    [ -maxRelSimDiff PERCENT ] \\\n" &
    "    [ -minSimInter NAT ] \\\n" &
    "    [ -maxSimPref NAT ] \\\n" &
    "    [ -maxSimGen NAT ] \\\n" &
    "    [ -simSugg FILE.ssg ] \\\n" &
    "    [ -maxSimSugg NAT ] \\\n" &
    "  ] \\\n" &
    "  [ -classes FILE.cls \\\n" &
    "    [ -minClassPrefs NAT ] \\\n" &
    "    [ -radicals FILE.rad ] \\\n" &
    "  ] \\\n" &
    "  [ -log FILE.log ] \\\n" &
    "  [ -maxPrefSize NAT ] \\\n" &
    "  [ -maxSuffSize NAT ] \n";

VAR
  encoding: Encoding.T := PlainEncoding.New();

TYPE Options = RECORD (* Command line options: *)

    (* General arguments: *)

    log: Wr.T;           (* The message file. *)
    cmd: TEXT;           (* Formatted command line arguments. *)

    load: TEXT;          (* The input automaton. *)

    (* Main operations ("" if not requested): *)

    prefs: TEXT;         (* Filename for prefix listing of specified states. *)
    suffs: TEXT;         (* Filename for prefix listing of specified states. *)
    prefSuffs: TEXT;     (* Filename for prefix/suffix listing of specified states. *)
    unprod: TEXT;        (* Filename for list of unproductive states. *)
    prod: TEXT;          (* Filename for list of productive states. *)
    similar: TEXT;       (* Filename for list of similar states. *)
    classes: TEXT;       (* Filename for list of lexical classes. *)

    (* Secondary parameters (filenames are "" if not specified): *)

    (* For "-prefs", "-suffs", "-prefSuffs", "-similar": *)
    maxPrefSize: NAT;    (* Maximum number of prefix bytes to print per state. *)
    maxSuffSize: NAT;    (* Maximum number of suffix bytes to print per state. *)

    (* For "-unprod": *)
    maxUnprod: NAT;      (* A state is unproductive if it is used by this many words. *)
    unprodSugg: TEXT;    (* Filename for listing of "strange" words. *)
    maxUnprodSugg: NAT;  (* Maximum total number of "unproductive" words to list. *)

    (* For "-prod": *)
    minProductivity: NAT;(* A state is "productive" if "(NP-1)*(NS-1)" is this or more.  *)
    maxProdStates: NAT;  (* Maximum number of productive states to list. *)

    (* For "-similar": *)
    maxSimDiff: NAT;     (* Maximum absolute difference between suffix sets. *)
    maxRelSimDiff: NAT;  (* Maximum relative difference between suffix sets. *)
    maxSimGen: NAT;      (* Maximum prefix-weighted difference between suffix sets. *)
    minSimInter: NAT;    (* Minimum number of common suffixes. *)
    maxSimPref: NAT;     (* Maximum number of prefixes. *)
    simSugg: TEXT;       (* Filename for listing of "possibly missing" words. *)
    maxSimSugg: NAT;     (* Maximum number of suggestions per state pair. *)

    (* For "-classes": *)
    minClassPrefs: NAT;  (* Minimum number of prefixes for lexical class. *)
    radicals: TEXT;      (* Filename for table of radicals. *)
  END;
    
VAR 
  log: Wr.T;

CONST
  DefaultMaxPrefSize = LAST(NAT);
  DefaultMaxSuffSize = LAST(NAT);

  DefaultMaxUnprod = 1;
  DefaultMaxUnprodSugg = LAST(NAT);

  DefaultMinProductivity = 4;
  DefaultMaxProdStates = LAST(NAT);

  DefaultMaxSimDiff = 2;
  DefaultMaxRelSimDiff = 100;  (* in % *)
  DefaultMinSimInter = 1;
  DefaultMaxSimPref = LAST(NAT);
  DefaultmaxSimSugg = LAST(NAT);
  DefaultMaxSimGen:  NAT = LAST(NAT);

  DefaultMinClassPrefs = 2;

PROCEDURE Main() =
  BEGIN

    WITH arg = GetCommandLineArguments() DO
      log := arg.log;

      StartedMessage("AutoAnalysis " & Version & " " & Util.FmtDate(Util.GetDate())); 

      DoShowArgs(arg.cmd);

      WITH
        aut = DoLoad(arg.load)
      DO
        IF NOT Empty(arg.prefSuffs) THEN
          ListPrefSuffs(
            arg.prefSuffs, aut, "-prefSuffs",
            arg.maxPrefSize, arg.maxSuffSize
          )
        END;
        IF NOT Empty(arg.suffs) THEN
          ListPrefSuffs(arg.suffs, aut, "-suffs", 0, arg.maxSuffSize)
        END;
        IF NOT Empty(arg.prefs) THEN
          ListPrefSuffs(arg.prefs, aut, "-prefs", arg.maxPrefSize, 0)
        END;

        IF NOT Empty(arg.unprod) THEN
          FindUnproductiveStates(
            arg.unprod, arg.unprodSugg, aut, 
            maxUnprod := arg.maxUnprod, 
            maxUnprodSugg := arg.maxUnprodSugg
          )
        END;

        IF NOT Empty(arg.prod) THEN
          FindProductiveStates(
            arg.prod, aut, 
            minProductivity := arg.minProductivity, 
            maxProdStates := arg.maxProdStates,
            maxPrefSize := arg.maxPrefSize, 
            maxSuffSize := arg.maxSuffSize
          )
        END;

        IF NOT Empty(arg.similar) THEN 
          FindSimilarStates(
            arg.similar, arg.simSugg, aut, 
            maxPrefSize := arg.maxPrefSize, 
            maxSuffSize := arg.maxSuffSize, 
            maxSimDiff := arg.maxSimDiff,
            maxRelSimDiff := arg.maxRelSimDiff,
            maxSimGen := arg.maxSimGen,
            minSimInter := arg.minSimInter,
            maxSimPref := arg.maxSimPref,
            maxSimSugg := arg.maxSimSugg
          )
        END;

        (* Lexical classes: *)
        IF NOT Empty(arg.classes) THEN 
          WITH
            nRadicals = NEW(REF ARRAY OF NAT, aut.Root()+1)^
          DO
            ComputeClasses(aut, nRadicals, minClassPrefs := arg.minClassPrefs);
            OutputClasses(
              arg.classes, aut, nRadicals, 
              onlyOneRadical := TRUE, tag := "classes"
            );
            IF NOT Empty(arg.radicals) THEN
              OutputClasses(
                arg.radicals, aut, nRadicals,
                onlyOneRadical := FALSE, tag := "radicals"
              )
            END;
          END
        END;

      END;
      FinishedMessage("AutoAnalysis " & Version);
      Wr.Close(log)
    END;
  END Main;
   
PROCEDURE GetCommandLineArguments(): Options =
  VAR arg: Options;
  BEGIN
    WITH 
      pp = NEW(ParseParams.T).init(stderr),
      cmdWr = NEW(TextWr.T).init()
    DO
      TRY
        IF pp.keywordPresent("-help") THEN 
          Wr.PutText(stderr, "Usage: ");
          Wr.PutText(stderr, Help);
          Process.Exit(0);
          <* ASSERT FALSE *>
        END;
        
        Wr.PutText(cmdWr, "AutoAnalysis");
        
        arg.log := ParamUtil.GetLogWr(pp, cmdWr);

        arg.load := ParamUtil.GetFileName(pp, cmdWr, "-load");
        IF Empty(arg.load) THEN 
          arg.load := ParamUtil.GetFileName(pp, cmdWr, "-dump");
        END;
        IF Empty(arg.load) THEN 
          arg.load := ParamUtil.GetFileName(pp, cmdWr, "-dumpIn");
        END;
        
        arg.prefs     := ParamUtil.GetFileName(pp, cmdWr, "-prefs");
        arg.suffs     := ParamUtil.GetFileName(pp, cmdWr, "-suffs");
        arg.prefSuffs := ParamUtil.GetFileName(pp, cmdWr, "-prefSuffs");
        arg.unprod    := ParamUtil.GetFileName(pp, cmdWr, "-unprod");
        arg.prod      := ParamUtil.GetFileName(pp, cmdWr, "-prod");
        arg.similar   := ParamUtil.GetFileName(pp, cmdWr, "-similar");
        arg.classes   := ParamUtil.GetFileName(pp, cmdWr, "-classes");
  
        IF  Empty(arg.prefs)
        AND Empty(arg.suffs)
        AND Empty(arg.prefSuffs)
        AND Empty(arg.unprod)
        AND Empty(arg.prod)
        AND Empty(arg.similar)
        AND Empty(arg.classes)
        THEN
          pp.error("no operation specified")
        END;
        
        IF ORD(NOT Empty(arg.prefs))
        +  ORD(NOT Empty(arg.suffs))
        +  ORD(NOT Empty(arg.prefSuffs)) > 1 THEN
          pp.error("commands \"-prefs\", \"-suffs\", and \"-prefSufss\"" &
            " are mutually exclusive")
        END;
        
        IF NOT Empty(arg.prefs) 
        OR NOT Empty(arg.suffs) 
        OR NOT Empty(arg.prefSuffs) 
        OR NOT Empty(arg.similar) 
        OR NOT Empty(arg.prod) 
        THEN
          arg.maxPrefSize  := ParamUtil.GetInt(pp, cmdWr, "-maxPrefSize",
            DefaultMaxPrefSize,   1, LAST(NAT)
          );
          arg.maxSuffSize  := ParamUtil.GetInt(pp, cmdWr, "-maxSuffSize",
            DefaultMaxSuffSize,   1, LAST(NAT)
          );
        END;
        
        IF NOT Empty(arg.unprod) THEN
          arg.unprodSugg := ParamUtil.GetFileName(pp, cmdWr, "-unprodSugg");
          arg.maxUnprod := ParamUtil.GetInt(pp, cmdWr, "-maxUnprod",
            DefaultMaxUnprod,     1, LAST(NAT)
          );
          arg.maxUnprodSugg := ParamUtil.GetInt(pp, cmdWr, "-maxUnprodSugg",
            DefaultMaxUnprodSugg, 1, LAST(NAT)
          );
        END;

        IF NOT Empty(arg.prod) THEN
          arg.minProductivity := ParamUtil.GetInt(pp, cmdWr, "-minProductivity",
            DefaultMinProductivity, 1, LAST(NAT)
          );
          arg.maxProdStates := ParamUtil.GetInt(pp, cmdWr, "-maxProdStates",
            DefaultMaxProdStates, 1, LAST(NAT)
          );
        END;

        IF NOT Empty(arg.similar) THEN
          arg.simSugg := ParamUtil.GetFileName(pp, cmdWr, "-simSugg");
          arg.maxSimDiff    := ParamUtil.GetInt(pp, cmdWr, "-maxSimDiff", 
            DefaultMaxSimDiff,    1, LAST(NAT)
          );
          arg.maxRelSimDiff := ParamUtil.GetInt(pp, cmdWr, "-maxRelSimDiff",
            DefaultMaxRelSimDiff, 1, 100
          );
          arg.minSimInter   := ParamUtil.GetInt(pp, cmdWr, "-minSimInter",
            DefaultMinSimInter,   0, LAST(NAT)
          );
          arg.maxSimPref    := ParamUtil.GetInt(pp, cmdWr, "-maxSimPref",
            DefaultMaxSimPref,    1, LAST(NAT)
          );
          arg.maxSimSugg    := ParamUtil.GetInt(pp, cmdWr, "-maxSimSugg",
            DefaultmaxSimSugg,    1, LAST(NAT)
          );
          arg.maxSimGen     := ParamUtil.GetInt(pp, cmdWr, "-maxSimGen",
            DefaultMaxSimGen,     1, LAST(NAT)
          );
        END;
        
        IF NOT Empty(arg.classes) THEN
          arg.minClassPrefs := ParamUtil.GetInt(pp, cmdWr, "-minClassPrefs",
            DefaultMinClassPrefs, 1, LAST(NAT)
          );
          arg.radicals := ParamUtil.GetFileName(pp, cmdWr, "-radicals");
        END;
        
        pp.finish();
        arg.cmd := TextWr.ToText(cmdWr);
        RETURN arg
      EXCEPT
      | ParseParams.Error =>
          Wr.PutText(stderr, "Usage: ");
          Wr.PutText(stderr, Help);
          Process.Exit(1);
          <* ASSERT FALSE *>
      END
    END
  END GetCommandLineArguments;

PROCEDURE DoShowArgs(cmd: TEXT) =
  BEGIN
    NL();
    Wr.PutText(log, cmd);
    NL();
    NL();
  END DoShowArgs;

PROCEDURE DoLoad(fileName: TEXT): Reduced.T =
  VAR aut: Reduced.T;
  BEGIN
    StartedMessage("loading automaton");
    aut := Reduced.Load(Util.OpenRd(fileName, "-load"));
    Util.PrintDoc(log, aut.doc, "  |");
    FinishedMessage("loading automaton");
    RETURN aut
  END DoLoad;
  
PROCEDURE ListPrefSuffs(
    fileName: TEXT;
    aut: Reduced.T;
    tag: TEXT;
    maxPrefSize: NAT;
    maxSuffSize: NAT;
  ) =

  PROCEDURE BadStateMessage(wr: Wr.T; msg: TEXT) =
    BEGIN
      Wr.PutText(wr, "*** ");
      Wr.PutText(wr, msg);
      Wr.PutText(wr, " \n");
      Wr.Flush(wr);
    END BadStateMessage;

  VAR
    st: State;

    (* Automaton used to sort the prefixes, and its root: *)
    autP: Reduced.T := Reduced.New(10000);
    rootP: State := autP.Root();

  BEGIN
    StartedMessage("writing " & tag & " file");
    WITH
      wr = Util.OpenWr(fileName, tag),
      prefpr = StringPrinter.New(wr, encoding, maxChars := maxPrefSize),
      suffpr = StringPrinter.New(wr, encoding, maxChars := maxSuffSize),
      root = aut.Root()
    DO
      LOOP
        st := RdNat(stdin);
        WrNat(log, st); NL();
        IF st = NullState OR st > root THEN 
          BadStateMessage(log, "invalid state");
        ELSE
          Wr.PutChar(wr, '\n');
          Wr.PutText(wr, "State ");
          WrNat(wr, st);
          Wr.PutText(wr, ";   ");
          WrNat(wr, aut.NPrefs(st));
          Wr.PutText(wr, " prefixes,  ");
          WrNat(wr, aut.NSuffs(st));
          Wr.PutText(wr, " suffixes\n\n");
          IF maxPrefSize > 0 THEN
            (* Enumerate prefixes that lead to "st".  To get them sorted, *)
            (* build an automaton "autP" and then enumerate it. *)
            PROCEDURE PrefixActionProc (* : Reduced PrefixAction*) (
                READONLY w: String;
              ) =
              BEGIN
                WITH 
                  nw = NUMBER(w),
                  nl = nw - 1,
                  rs = NEW(REF String, nw)
                DO
                  FOR i := 0 TO nl DO rs^[i] := w[nl-i] END;
                  rootP := autP.AddStringMaybeCrunch(rootP, rs^);
                END
              END PrefixActionProc;

            <* FATAL Basics.Abort *>
            BEGIN
              prefpr.Reset();
              autP.SetRoot(NullState);
              autP.Crunch();
              rootP := NullState;
              Wr.PutText(wr, "Prefixes:\n\n{ ");
              (* The test below will not be necessary when Reduced.m3 is fixed. *)
              IF aut.NPrefs(st) > 0 THEN 
                aut.EnumPrefs(st, PrefixActionProc);
                autP.PrintSuffs(rootP, prefpr)
              END;
              Wr.PutText(wr, " }\n");
              Wr.Flush(wr)
            END
          END;
          IF maxSuffSize > 0 THEN
            Wr.PutText(wr,"\nSuffixes:\n\n{ ");
            suffpr.Reset();
            aut.PrintSuffs(st, suffpr);
            Wr.PutText(wr, " }\n");
            Wr.Flush(wr)
          END
        END;
        Lex.Skip(stdin);
        IF Rd.EOF(stdin) THEN EXIT END 
      END;
      Wr.Close(wr)
    END;
    FinishedMessage("writing " & tag & " file");
  END ListPrefSuffs;

PROCEDURE FindUnproductiveStates(
    stateFileName: TEXT;
    suggFileName: TEXT;
    aut: Reduced.T;
    maxUnprod: NAT;
    maxUnprodSugg: NAT;
  ) =
  VAR
    stateCount: NAT := 0;
    suggCount: NAT := 0;
    
  <* FATAL Basics.Abort *>
  BEGIN
    StartedMessage("looking for unproductive states");
    WITH
      stWr = Util.OpenWr(stateFileName, "-unprod"),
      suggWr = Util.OpenWr(suggFileName, "-unprodSugg"),
      uppr = StringPrinter.New(stWr, encoding, maxChars := LAST(NAT)),
      root = aut.Root()
    DO
      FOR st := UnitState TO root DO
        WITH
          ns = aut.NSuffs(st),
          np = aut.NPrefs(st),
          nw = ns*np
        DO
          IF nw > 0 AND nw <= maxUnprod THEN
            IF stateCount >= maxUnprodSugg THEN
              Wr.PutText(stWr,"\n(more...)\n");
              EXIT 
            END;
            INC(stateCount);
            Wr.PutText(stWr, Fmt.Pad(Fmt.Int(st), 7));
            Wr.PutText(stWr, " np = ");
            Wr.PutText(stWr, Fmt.Int(np));
            Wr.PutText(stWr, " ns = ");
            Wr.PutText(stWr, Fmt.Int(ns));
            Wr.PutText(stWr, "  { ");
            uppr.Reset();
            aut.PrintPrefs(st, uppr);
            Wr.PutText(stWr," }:{ ");
            aut.PrintSuffs(st, uppr);
            Wr.PutText(stWr," }\n");
            Wr.Flush(stWr);
            IF suggWr # NIL THEN

              PROCEDURE PrefixAction(READONLY pw: String ) =

                PROCEDURE SuffixAction(READONLY sw: String) =
                  BEGIN
                    INC(suggCount);
                    FOR j := LAST(pw) TO 0 BY -1 DO 
                      encoding.PrintLetter(suggWr,pw[j])
                    END;
                    FOR j := 0 TO LAST(sw) DO 
                      encoding.PrintLetter(suggWr,sw[j])
                    END;
                    Wr.PutChar(suggWr,'\n')
                  END SuffixAction;

                BEGIN
                  aut.EnumSuffs(st, SuffixAction)
                END PrefixAction;

              BEGIN
                aut.EnumPrefs(st, PrefixAction);
                Wr.Flush(suggWr)
              END
            END
          END
        END (* WITH *)
      END (* FOR *);
      Wr.Close(stWr);
      Wr.PutText(log, 
        Fmt.Pad(Fmt.Int(stateCount), 9) & " unproductive states\n"
      );
      IF suggWr # NIL THEN 
        Wr.Close(suggWr);
        Wr.PutText(log,
          Fmt.Pad(Fmt.Int(suggCount), 9) & " strange words (with repetitions) listed\n"
        )
      END;
    END (* WITH *);
    FinishedMessage("looking for unproductive states");
  END FindUnproductiveStates;
  
PROCEDURE FindProductiveStates(
    stateFileName: TEXT;
    aut: Reduced.T;
    minProductivity: NAT; (* Min value of (NS-1)*(NP-1) *)
    maxProdStates: NAT;   (* Max states to print *)
    maxPrefSize: NAT;     (* Maximum number of prefix bytes to print per state *)
    maxSuffSize: NAT;     (* Maximum number of suffix bytes to print per state *)
  ) =
  VAR
    stateCount: NAT := 0;
  BEGIN
    StartedMessage("looking for productive states");
    WITH
      stWr = Util.OpenWr(stateFileName, "-prod"),
      prefPr = StringPrinter.New(stWr, encoding, maxChars := maxPrefSize),
      suffPr = StringPrinter.New(stWr, encoding, maxChars := maxSuffSize),
      root = aut.Root()
    DO
      Wr.PutText(stWr, "  state  nprefs  nsuffs  nwords  prodty  prefs/suffs\n");
      Wr.PutText(stWr, "------- ------- ------- ------- -------  -----------------\n");
      FOR st := UnitState TO root DO
        WITH
          ns = aut.NSuffs(st),
          np = aut.NPrefs(st),
          nw = ns*np,
          py = (ns - 1)*(np - 1)
        DO
          IF py > minProductivity THEN
            IF stateCount >= maxProdStates THEN
              Wr.PutText(stWr, "\n(more...)\n");
              EXIT 
            END;
            INC(stateCount);
            Wr.PutText(stWr, Fmt.Pad(Fmt.Int(st), 7));
            Wr.PutText(stWr, Fmt.Pad(Fmt.Int(np), 8));
            Wr.PutText(stWr, Fmt.Pad(Fmt.Int(ns), 8));
            Wr.PutText(stWr, Fmt.Pad(Fmt.Int(nw), 8));
            Wr.PutText(stWr, Fmt.Pad(Fmt.Int(py), 8));
            Wr.PutText(stWr, "  { ");
            prefPr.Reset();
            aut.PrintPrefs(st, prefPr);
            Wr.PutText(stWr," }:{ ");
            suffPr.Reset();
            aut.PrintSuffs(st, prefPr);
            Wr.PutText(stWr," }\n");
            Wr.Flush(stWr);
          END
        END (* WITH *)
      END (* FOR *);
      Wr.Close(stWr);
      Wr.PutText(log, 
        Fmt.Pad(Fmt.Int(stateCount), 9) & " productive states\n"
      );
    END (* WITH *);
    FinishedMessage("looking for productive states");
  END FindProductiveStates;
  
PROCEDURE FindSimilarStates(
    stateFileName: TEXT;
    suggFileName: TEXT;
    aut: Reduced.T;
    maxPrefSize: NAT;    (* Maximum number of prefix bytes to print per state *)
    maxSuffSize: NAT;    (* Maximum number of suffix bytes to print per state *)
    maxSimDiff: NAT;     (* Maximum absolute difference between suffix sets. *)
    maxRelSimDiff: NAT;  (* Maximum relative difference between suffix sets. *)
    maxSimGen: NAT;      (* Maximum prefix-weighted difference between suffix sets. *)
    minSimInter: NAT;    (* Minimum number of common suffixes *)
    maxSimPref: NAT;     (* Maximum number of prefixes *)
    maxSimSugg: NAT;     (* Maximum number of suggestions per state pair *)
  ) =
  VAR
    dst: NAT;     (* | Suffs(s) - Suffs(t) | *)
    dts: NAT;     (* | Suffs(t) - Suffs(s) | *)
    
    simf: BOOL;   (* TRUE iff states s, t are similar *)

    simCount: NAT := 0;     (* Counts pairs of similar states *)
    simSuggCount: NAT := 0; (* Counts suggested words *)

    nGood: NAT := 0; (* Number of "good" states to compare *)
    
  PROCEDURE MaxUncovered(nSuffs: NAT): NAT =
    (* 
      Given NSuffs(s), computes how many of those suffixes
      may be left uncovered by a similar state. *)
    BEGIN
      WITH
        maxRDiff = CEILING(FLOAT(nSuffs)*100.0/FLOAT(maxRelSimDiff))
      DO
        RETURN MIN(MIN(maxSimDiff, maxRDiff), nSuffs - minSimInter)
      END
    END MaxUncovered;


  BEGIN (* FindSimilarStates *)
    StartedMessage("Computing similar states");
    WITH
      simWr  = Util.OpenWr(stateFileName, "-similar"),
      suggWr = Util.OpenWr(suggFileName,  "-simSugg"),
      sppr = StringPrinter.New(simWr, encoding, maxChars := maxPrefSize),
      sspr = StringPrinter.New(simWr, encoding, maxChars := maxSuffSize),
      root = aut.Root(),
      
      nSuffGoodRef = NEW(REF ARRAY OF INTEGER, root+1),
      stateGoodRef = NEW(REF ARRAY OF INTEGER, root+1),

      autPair = ReducedPair.New(aut, aut)
    DO

      PROCEDURE OutputStatePair(s, t: State; stint: NAT) =
        BEGIN
          WITH
            stp = ReducedPair.State{s, t},
            tsp = ReducedPair.State{t, s},
            nspref = aut.NPrefs(s),
            nssuff = aut.NSuffs(s),
            ntpref = aut.NPrefs(t),
            ntsuff = aut.NSuffs(t),
            dst = nssuff - stint,
            dts = ntsuff - stint,
            stsug = nspref*dts,
            tssug = ntpref*dst,
            totsug = stsug+tssug
          DO
            Wr.PutText(simWr, "\n\nStates:  " & 
              " s = " & Fmt.Int(s) &
              " (" & Fmt.Int(nspref) & "," & Fmt.Int(nssuff) & ")  " &
              " t = " & Fmt.Int(t) &  
              " (" & Fmt.Int(ntpref) & "," & Fmt.Int(ntsuff) & ")  " &
              "  " & Fmt.Int(totsug) & " suggestions" &
              "\n"
            );
            Wr.PutText(simWr, "\nPref(s)\n{ ");
            sppr.Reset();
            aut.PrintPrefs(s, sppr);
            Wr.PutText(simWr, " }\n");
            Wr.PutText(simWr, "\nPref(t)\n{ ");
            sppr.Reset();
            aut.PrintPrefs(t, sppr);
            Wr.PutText(simWr, " }\n");
            Wr.PutText(simWr, "\nSuff(s) \\ Suff(t)  (");
            Wr.PutText(simWr, Fmt.Int(dst));
            Wr.PutText(simWr, ", ");
            Wr.PutText(simWr, Fmt.Int(stsug));
            Wr.PutText(simWr, ")\n{ ");
            sspr.Reset();
            autPair.PrintSuffs(stp, ReducedPair.BoolOp.Diff, sspr);
            Wr.PutText(simWr, " }\n");
            Wr.PutText(simWr, "\nSuff(t) \\ Suff(s)  (");
            Wr.PutText(simWr, Fmt.Int(dts));
            Wr.PutText(simWr, ", ");
            Wr.PutText(simWr, Fmt.Int(tssug));
            Wr.PutText(simWr, ")\n{ ");
            sspr.Reset();
            autPair.PrintSuffs(tsp, ReducedPair.BoolOp.Diff, sspr);
            Wr.PutText(simWr, " }\n");
            Wr.PutText(simWr, "\nSuff(t) & Suff(s)  (");
            Wr.PutText(simWr, Fmt.Int(stint));
            Wr.PutText(simWr, ")\n{ ");
            sspr.Reset();
            autPair.PrintSuffs(tsp, ReducedPair.BoolOp.Inter, sspr);
            Wr.PutText(simWr," }\n");
            Wr.Flush(simWr);
          END
        END OutputStatePair;

      PROCEDURE OutputSuggestions(s, t: State) =
        VAR pairSuggCount: NAT := 0;
        BEGIN
          WITH
            st = ARRAY BOOL OF State{s, t}
          DO
            TRY
              FOR which := FALSE TO TRUE DO

                PROCEDURE PrefixAction(READONLY pw: String) RAISES {Abort} =

                  PROCEDURE SuffixAction(READONLY sw: String) RAISES {Abort} =
                    BEGIN
                      IF pairSuggCount >= maxSimSugg THEN
                        Wr.PutText(suggWr, "\n(more...)\n");
                        RAISE Abort
                      END;
                      INC(pairSuggCount);
                      FOR j := LAST(pw) TO 0 BY -1 DO 
                        encoding.PrintLetter(suggWr, pw[j])
                      END;
                      FOR j := 0 TO LAST(sw) DO 
                        encoding.PrintLetter(suggWr, sw[j])
                      END;
                      Wr.PutChar(suggWr,'\n');
                      INC(simSuggCount)
                    END SuffixAction;

                  BEGIN
                    autPair.EnumSuffs(
                      ReducedPair.State{st[NOT which], st[which]}, 
                      ReducedPair.BoolOp.Diff, 
                      SuffixAction
                    )
                  END PrefixAction;

                BEGIN
                  Wr.PutText(suggWr, "Pref(");
                  Wr.PutText(suggWr, Fmt.Int(st[which]));
                  Wr.PutText(suggWr, ") x Suff(");
                  Wr.PutText(suggWr, Fmt.Int(st[NOT which]));
                  Wr.PutText(suggWr, ")\\Suff(");
                  Wr.PutText(suggWr, Fmt.Int(st[which]));
                  Wr.PutText(suggWr, ")\n");
                  pairSuggCount := 0;
                  aut.EnumPrefs(st[which], PrefixAction);
                  Wr.PutChar(suggWr, '\n');
                  Wr.Flush(suggWr)
                END
              END
            EXCEPT
              Abort  =>  Wr.Flush(suggWr) (* OK *)
            END
          END;
        END OutputSuggestions;

      BEGIN
        (* Collect "good" states that are worth comparing: *)
        WITH
          nSuffGood = nSuffGoodRef^,
          stateGood = stateGoodRef^
        DO
          FOR s := UnitState TO root - 1 DO 
            WITH nspref = aut.NPrefs(s), nssuff = aut.NSuffs(s) DO
              IF nssuff > 0 
              AND nspref > 0 
              AND nspref <= maxSimPref
              AND nssuff >= minSimInter
              THEN
                nSuffGood[nGood] := nssuff;
                stateGood[nGood] := s;
                INC(nGood)
              END
            END
          END
        END;
        
        WITH
          nSuffGood = SUBARRAY(nSuffGoodRef^, 0, nGood),
          stateGood = SUBARRAY(stateGoodRef^, 0, nGood)
        DO
          (* Sort "good" states by incresing number of suffixes: *)
          VAR t := stateGoodRef;
          BEGIN
            IntQuickSort.QuickSort(nSuffGood, t)
          END;
          
          (* Test all pairs of "good" states: *)
          FOR i := 0 TO nGood - 1 DO
            WITH
              s = stateGood[i],
              nssuff = nSuffGood[i],
              nspref = aut.NPrefs(s),
              stMaxDiff = MaxUncovered(nssuff)
            DO
              <* ASSERT nssuff > 0 *>
              <* ASSERT nspref > 0 *>
              <* ASSERT nspref <= maxSimPref *>
              <* ASSERT nssuff >= minSimInter *>
              FOR j := i + 1 TO nGood - 1 DO
                WITH
                  t = stateGood[j],
                  ntsuff = nSuffGood[j],
                  ntpref = aut.NPrefs(t),
                  tsMaxDiff = MaxUncovered(ntsuff)
                DO
                  <* ASSERT ntsuff >= nssuff *>
                  IF ntsuff - nssuff > tsMaxDiff THEN EXIT END;
                  IF tsMaxDiff = 0 AND stMaxDiff = 0 THEN EXIT END;
                  WITH
                    stp = ReducedPair.State{s, t},
                    tsp = ReducedPair.State{t, s}
                  DO
                    simf := TRUE; (* in principle... *)
                    TRY
                      dst := autPair.NSuffs(
                        stp, 
                        ReducedPair.BoolOp.Diff, 
                        limit := stMaxDiff + 1
                      );
                      simf := simf AND dst <= stMaxDiff;
                      IF simf THEN
                        dts := autPair.NSuffs(
                          tsp,
                          ReducedPair.BoolOp.Diff,
                          limit := tsMaxDiff + 1
                        );
                        <* ASSERT (nssuff - dst) = (ntsuff - dts) *>
                        <* ASSERT dst > 0 OR dts > 0 *>
                        simf := simf AND dts <= tsMaxDiff
                      END
                    EXCEPT
                      Abort => simf := FALSE;
                    END;
                    IF simf THEN
                      WITH
                        stsug = nspref*dts,
                        tssug = ntpref*dst,
                        totsug = stsug+tssug
                      DO
                        IF totsug <= maxSimGen THEN
                          INC(simCount);
                          OutputStatePair(s, t, nssuff - dst);
                          IF suggWr # NIL THEN
                            OutputSuggestions(s, t)
                          END
                        END
                      END
                    END;
                  END
                END
              END
            END
          END
        END
      END;
      Wr.Close(simWr);
      Wr.PutText(log, Fmt.Pad(Fmt.Int(simCount), 9));
      Wr.PutText(log, " pairs of similar states found\n");
      IF suggWr # NIL THEN 
        Wr.Close(suggWr);
        Wr.PutText(log, Fmt.Pad(Fmt.Int(simSuggCount), 9));
        Wr.PutText(log, " words (with repetitions) suggested\n");
      END
    END;
    FinishedMessage("computing similar states");
  END FindSimilarStates;

PROCEDURE ComputeClasses(
    aut: Reduced.T;
    VAR nRadicals: ARRAY OF NAT;
    minClassPrefs: NAT;
  ) =
  (*
    For each proper state "s": if "s" is a lexical class, sets
    "nRadicals[s]" to the number of prefixes of "s" that do not go
    though another lexical class; otherwise sets "nRadicals" to 0. *)
  
  VAR
    lexClassCount: NAT := 0;
    
  PROCEDURE ProcessRadicalPrefix(
      <*UNUSED*> len: NAT;
      s: State;
      <*UNUSED*> f: BOOL;
    ) RAISES {Skip} =
    (*
      Called for every path from the root that does not go through any lexical
      class.  Checks whether the final state "s" of the path is a lexical class;
      if so, bumps its radical count, and prevents the enumeration
      from proceeding through "s". *)
    BEGIN
      IF aut.NPrefs(s) >= minClassPrefs THEN
        (* "s" is a lexical class. *)
        WITH nrs = nRadicals[s] DO
          IF nrs = 0 THEN INC(lexClassCount) END;
          INC(nrs)
        END;
        (* Do not go past this state: *)
        RAISE Skip
      END
    END ProcessRadicalPrefix;

  <* FATAL Basics.Abort *>
  BEGIN
    StartedMessage("computing lexical classes");
    WITH
      root = aut.Root()
    DO
      FOR s := NullState TO root DO nRadicals[s] := 0 END;
      aut.EnumPaths(root, ProcessRadicalPrefix, NIL, NIL, NIL);
      Wr.PutText(log, Fmt.Pad(Fmt.Int(lexClassCount), 7));
      Wr.PutText(log, " lexical classes found\n");
    END;
    FinishedMessage("computing lexical classes")
  END ComputeClasses;
  
PROCEDURE PrintAligned(wr: Wr.T; t: TEXT; width: NAT) =
  BEGIN
    WITH len = Text.Length(t) DO
      FOR k := 1 TO width-len DO Wr.PutChar(wr, ' ') END;
    END;
    Wr.PutText(wr, t)
  END PrintAligned;

PROCEDURE PrintAlignedNat(wr: Wr.T; x: NAT; width: NAT) =
  VAR t: NAT := x; len: NAT := 1;
  BEGIN
    WHILE t > 9 DO INC(len); t := t DIV 10 END;
    FOR k := 1 TO width-len DO Wr.PutChar(wr, ' ') END;
    Basics.WrNat(wr, x)
  END PrintAlignedNat;

PROCEDURE PrintAlignedReal(wr: Wr.T; x: REAL; width, decs: NAT) =
  VAR mul: NAT := 1;
  BEGIN
    (* Scale up to integer: *)
    FOR i := 1 TO decs DO mul := mul*10 END;
    WITH
      ax = ABS(x)*FLOAT(mul)
    DO
      IF ax - 0.5 >= FLOAT(LAST(INTEGER)) THEN 
        (* Give up - use Modula-3 routines *)
        Wr.PutText(wr, Fmt.Pad(Fmt.Real(x, prec := decs, style := Fmt.Style.Fix), width))
      ELSE
        VAR tx: NAT := ROUND(ax); 
            len: NAT := 0;
            digs: ARRAY [0..30] OF CHAR;
            minus: BOOL := tx # 0 AND x < 0.0;
        BEGIN
          (* Break into digits: *)
          WHILE tx > 9 DO 
            digs[len] := VAL(ORD('0') + tx MOD 10, CHAR);
            INC(len);
            tx := tx DIV 10
          END;
          digs[len] := VAL(ORD('0') + tx, CHAR);
          INC(len);
          (* Print them out: *)
          IF minus THEN digs[len] := '-'; INC(len) END;
          FOR k := 1 TO width-len-1 DO Wr.PutChar(wr, ' ') END;
          WHILE len > decs DO DEC(len); Wr.PutChar(wr, digs[len]) END;
          Wr.PutChar(wr, '.');
          WHILE len > 0 DO DEC(len); Wr.PutChar(wr, digs[len]) END;
        END
      END
    END
  END PrintAlignedReal;

PROCEDURE OutputClasses(
    fileName: TEXT;
    aut: Reduced.T;
    READONLY nRadicals: ARRAY OF NAT;
    onlyOneRadical: BOOL;  (* TRUE prints only one radical per class *)
    tag: TEXT;             (* "classes" or "radicals" *)
  ) =
  
  BEGIN
    StartedMessage("writing " & tag & " file");
    WITH
      wr = Util.OpenWr(fileName, "-" & tag),
      pr = StringPrinter.New(wr, encoding, maxChars := 100, rightMargin := 120),
      classSeen = NEW(REF ARRAY OF BOOL, NUMBER(nRadicals))^
    DO

      PROCEDURE PrintHeader1() =
        BEGIN
          PrintAligned(wr,  "State",  6);
          PrintAligned(wr, "NPrefs",  8);
          PrintAligned(wr, "NSuffs",  8);
          PrintAligned(wr, "NWords",  8);
          PrintAligned(wr, "NRadcs",  8);
          PrintAligned(wr, "NClWds",  8);
          PrintAligned(wr,    "NW/",  8);
          PrintAligned(wr,   "NCW/",  8);
          Wr.PutText(wr, " ");
          Wr.PutText(wr, Fmt.Pad("FirstPrefix", 15, align := Fmt.Align.Left));
          Wr.PutText(wr, " ");
          Wr.PutText(wr, "Suffixes");
          Wr.PutText(wr, "\n");
          Wr.Flush(wr)
        END PrintHeader1;
        
      PROCEDURE PrintHeader2() =
        BEGIN
          PrintAligned(wr,        "",  6);
          PrintAligned(wr,        "",  8);
          PrintAligned(wr,        "",  8);
          PrintAligned(wr,        "",  8);
          PrintAligned(wr,        "",  8);
          PrintAligned(wr,        "",  8);
          PrintAligned(wr, "NS+NP-1",  8);
          PrintAligned(wr,      "NW",  8);
          Wr.PutText(wr, " ");
          Wr.PutText(wr, Fmt.Pad("", 15, align := Fmt.Align.Left));
          Wr.PutText(wr, " ");
          Wr.PutText(wr, "");
          Wr.PutText(wr, "\n");
          Wr.Flush(wr)
        END PrintHeader2;

      PROCEDURE PrintHeader3() =
        BEGIN
          PrintAligned(wr,  "-----",  6);
          PrintAligned(wr,  "-------",  8);
          PrintAligned(wr,  "-------",  8);
          PrintAligned(wr,  "-------",  8);
          PrintAligned(wr,  "-------",  8);
          PrintAligned(wr,  "-------",  8);
          PrintAligned(wr,  "-------",  8);
          PrintAligned(wr,  "-------",  8);
          Wr.PutText(wr, " ");
          Wr.PutText(wr, "---------------");
          Wr.PutText(wr, " ");
          Wr.PutText(wr, "---------------");
          Wr.PutText(wr, "\n");
          Wr.Flush(wr)
        END PrintHeader3;

      PROCEDURE PrintRadical(
          READONLY w: String;
          s: State;
        ) =
        (*
          State "s" is a proper class, and "w" is one of its prefixes. *)
        BEGIN
          WITH
            NP = aut.NPrefs(s),
            NS = aut.NSuffs(s),
            NW = NP * NS,
            NR = nRadicals[s],
            NCW = NR * NS,
            eff1 = FLOAT(NW)/FLOAT(NS+NP-1),
            eff2 = FLOAT(NCW)/FLOAT(NW)
          DO
            PrintAlignedNat(wr, s,    6);
            PrintAlignedNat(wr, NP,   8);
            PrintAlignedNat(wr, NS,   8);
            PrintAlignedNat(wr, NW,   8);
            PrintAlignedNat(wr, NR,   8);
            PrintAlignedNat(wr, NCW,  8);
            PrintAlignedReal(wr, eff1, 8, 4);
            PrintAlignedReal(wr, eff2, 8, 4);
            Wr.PutText(wr, " ");
            FOR i := 0 TO LAST(w) DO
              encoding.PrintLetter(wr, w[i])
            END;
            FOR i := LAST(w) + 1 TO 15 DO
              Wr.PutChar(wr,' ')
            END;
            Wr.PutText(wr, " { ");
            pr.Reset();
            aut.PrintSuffs(s, pr);
            Wr.PutText(wr," }\n");
            Wr.Flush(wr)
          END
        END PrintRadical;
        
      PROCEDURE PrintIsolatedWord(
          READONLY w: String;
          s: State;
        ) =
        (*
          State "s" is final but is not a proper class, and "w" is 
          a string that reaches "s" without going through
          a proper class. *)
        BEGIN
          WITH
            NP = aut.NPrefs(s),
            NS = 1,
            NW = NP * NS,
            NR = 0,
            NCW = NR * NS,
            eff1 = FLOAT(NW)/FLOAT(NS+NP-1),
            eff2 = FLOAT(NCW)/FLOAT(NW)
          DO
            PrintAlignedNat(wr, NullState,    8);
            PrintAlignedNat(wr, NP,   8);
            PrintAlignedNat(wr, NS,   8);
            PrintAlignedNat(wr, NW,   8);
            PrintAlignedNat(wr, NR,   8);
            PrintAlignedNat(wr, NCW,  8);
            PrintAlignedReal(wr, eff1, 8, 4);
            PrintAlignedReal(wr, eff2, 8, 4);
            Wr.PutText(wr, " ");
            FOR i := 0 TO LAST(w) DO
              encoding.PrintLetter(wr, w[i])
            END;
            FOR i := LAST(w) + 1 TO 15 DO
              Wr.PutChar(wr,' ')
            END;
            Wr.PutText(wr, " { () }");
            Wr.Flush(wr)
          END
        END PrintIsolatedWord;
        
      PROCEDURE VisitPrefix(
          READONLY w: String;
          s: State;
          <*UNUSED*> f: BOOL;
        ) RAISES{Skip} =
        BEGIN
          IF nRadicals[s] > 0 THEN
            IF NOT classSeen[s] OR NOT onlyOneRadical THEN
              PrintRadical(w, s)
            END;
            classSeen[s] := TRUE;
            RAISE Skip
          ELSIF aut.Final(s) AND NOT onlyOneRadical THEN
            PrintIsolatedWord(w, s)
          END
        END VisitPrefix;

      <* FATAL Basics.Abort *>
      BEGIN
        PrintHeader1();
        PrintHeader2();
        PrintHeader3();
        FOR s := NullState TO aut.Root() DO classSeen[s] := FALSE END;
        aut.EnumStrings(aut.Root(), VisitPrefix, NIL);
      END;
      Wr.Close(wr);
    END;
    FinishedMessage("writing  " & tag & " file");
  END OutputClasses;

PROCEDURE StartedMessage(m: TEXT)  =
  (*
    Prints "started m" to "log", with separating line. *)
  BEGIN
    NL();
    Message("started " & m);
    EVAL Util.GetExecTimes(); (* Reset the clock for FinishedMessage *)
  END StartedMessage;

PROCEDURE FinishedMessage(m: TEXT)  =
  (*
    Prints "finished m" to "log", with CPU time and separating line. *)
  BEGIN
    WITH tim = Util.GetExecTimesText().total DO
      Message("finished " & m & " (" & tim & ")")
    END;
    NL();
  END FinishedMessage;

CONST MessageWidth = 76;

PROCEDURE Message(m: TEXT)  =
  (*
    Prints a separating line to "log". *)
  <* FATAL Wr.Failure *>
  BEGIN
    Wr.PutText(log,"===");
    IF Empty(m) THEN
      Wr.PutText(log, "==")
    ELSE
      Wr.PutText(log, " ");
      Wr.PutText(log, m);
      Wr.PutText(log, " ")
    END;
    FOR i:=1 TO MessageWidth - 5 - Text.Length(m) DO 
      Wr.PutChar(log, '=')
    END;
    NL();
  END Message;

PROCEDURE NL()  =
  (* Prints a blank line to "log". *)
  <* FATAL Wr.Failure *>
  BEGIN
    Wr.PutText(log,"\n");
    Wr.Flush(log)
  END NL;

BEGIN
  Main()
END AutoAnalysis.
