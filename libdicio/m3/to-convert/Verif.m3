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

UNSAFE MODULE Verif EXPORTS Main;       (* Version 2.8 -- 23/Oct/92 *)

(* HISTORY:

25.august.92 (tk): Modified for compatibility with version 2.04 of 
                   Modula-3; added new options: -transpose, -remove, 
                   -subst, -insert
30.sept.92 (tk):   Small modifications when recompiled under 2.07

23.oct.92 (tk):    Inclusion OF the "-remote" option
06.May.93 (tk):    Intallation path for DCC changed to /proj/dicio
*)

(* This  program is used   to  implement   speller  checkers --  alpha
version. Includes some  special security features  in order to protect
the set of  words used  for   the automaton.  Runs only  at prescribed
installations and   within fixed directories.   The module  is  unsafe
because   of  some imported   system   modules.  In its  OpenWin  mode
communicates with an  interface written by  R. Anido. Under batch mode
works as a filter.

The program reads  the automaton  from a   file called  "cdados"  (see
/AutFileName/),  a set of  additional  words (i.e.  which  should have
been  included   in the automaton) from  a  file  called "extras" (see
/ExtrasFileName/) and a  set  of words to  be  considered  wrong (i.e.
which  should not have been included  in the automaton) from  the file
"tirar" (see /IgnoreFileName/).

The user has several options specified below, including  the inclusion
of up to four specific vocabularies of his own.

*)

(*

COMMAND LINE PARAMETERS (with defaults whenever applicable):
============================================================

-help                : prints the help the message
-1 FileName          : user's vocabulary # 1
-2 FileName          : user's vocabulary # 2
-3 FileName          : user's vocabulary # 3
-4 FileName          : user's vocabulary # 4
-keyword FileName    : keywords to be ignored
-transpose           : try alternatives by transposing two consecutive 
                       letters (one pair at a time) -- applied to words
                       of minimal length specified in the module /Adviser/
-remove              : try alternatives by removing one letter at a time 
                       -- applied to words of minimal length specified in 
                       the module /Adviser/
-subst               : try alternatives by substituting one letter at a time 
                       -- applied to words of minimal length specified in 
                       the module /Adviser/
-insert              : try alternatives by inserting one letter at a time 
                       between two consecutive letters -- applied to words 
                       of minimal length specified in the module /Adviser/
-s SiteCode          : specifies site code (at present DCC=10 and IME=20 -- 
                       see constants)
-m ModeCode          : specifies running mode: 0 -- batch, 1 -- OpenWin
                       interface 
-allwrong            : treats all input words as having been mispelled
-alt                 : when in batch mode, produces alternatives for misspelled
                       words
-t FileName          : file with secret code to allow testing (see TEST)
-d FileName          : when in test mode allows specification of the path where
                       to read the automaton FROM
                       
-remote String       : it is assumed that the program is run in batch mode through mail;
                       /String/ should be the sender's e-mail information (address AND name);
                       log information will be recorded in the file /rlog/ instead of /log/.
                       
*)

  IMPORT Wr, Rd, Fmt, Text, FileStream, Scan, 
         ParseParams, Thread,
         ReadOnlyPermDAG, Util, Adviser, AVLTextTree,
         M3toC, Upwd, Uugid;
  
  FROM Stdio IMPORT stdin,stdout,stderr;
  FROM Basics IMPORT Abort;
  FROM Util IMPORT BadVersion;
  
  <* FATAL Thread.Alerted *>

  EXCEPTION 
    Terminate; 
    Done;
  
  CONST
  
    Version = "versão 2.8 (23/out/92)";
  
    MaxPathLen = 40;
    ObLen = 5;
    
    DCC = 10;
    IME = 20;
    TEST = "42447481";  (* Guess what is that?! *)
    
    BATCH = 0;
    OPENWIN = 1;
    
    CodeOne = "0";
    CodeMany = "1";
    CodeInsert = "2";
    
    LOG = "log";
    XLOG = "xlog";
    RLOG = "rlog";
    
  VAR
    site: [DCC .. IME] := DCC;
    mode: [BATCH .. OPENWIN] := BATCH;
    remote: TEXT;
    
    AltFlag: BOOLEAN;
    RemoteFlag,
    TestFlag,
    AllWrongFlag,
    TransposeFlag,
    RemoveFlag,
    SubstFlag,
    InsertFlag,
    KeywordFlag: BOOLEAN := FALSE;
    
    KeywordFileName: TEXT := "";
    
    
    LocalTree := AVLTextTree.New();
    IgnoreTree := AVLTextTree.New();
    KeywordTree := AVLTextTree.New();
    ErrorTree := AVLTextTree.New();
    AltTree: AVLTextTree.T;
    
    LocalFlag1,
    LocalFlag2,
    LocalFlag3,
    LocalFlag4: BOOLEAN;
    
    LocalFileName1,
    LocalFileName2,
    LocalFileName3,
    LocalFileName4: TEXT := "";
    
    AutFileName: TEXT := "cdados";
    LogFileName: TEXT;
    ExtrasFileName: TEXT := "extras";
    IgnoreFileName: TEXT := "tirar";

    path: ARRAY [0..MaxPathLen-1] OF CHAR;
    object: ARRAY [0..ObLen-1] OF CHAR; 
    
    pathlen: CARDINAL;
    
    PUBPathName,
    notPUBPathName: TEXT;
    
    logwr: Wr.T; 
    obrd: Rd.T;
    
    InitTime,
    FinalTime,
    ExecTime: TEXT;
    
    WordCount,
    ErrorCount,
    NotFoundCount,
    FileCount: CARDINAL := 0;
    
    user: ADDRESS;
    username: TEXT :="";
    blanks: TEXT := "            ";
    
    Aut  : ReadOnlyPermDAG.T;
    Root : ReadOnlyPermDAG.State;

(* ---------------------------------------------------------------------- *)

  PROCEDURE Main() =
  <* FATAL Wr.Failure *>
  BEGIN
    TRY
    
      Initialize();
      
      GetOptions();
      
      SetUp();
      
      Check();

      Finalize()

    EXCEPT
    | Terminate => Wr.Flush(stderr); Wr.Flush(stdout)  
    END
  END Main;

  PROCEDURE Initialize() RAISES{Terminate} =
  <* FATAL Wr.Failure *>
  BEGIN
    TRY
      InitTime := Util.DateTime();
      EVAL Util.ExecTimesTexts().TotalTime;
      user := (Upwd.getpwuid(Uugid.getuid())).pw_name;
      IF user=NIL THEN
        Wr.PutText(stderr,"\nUsuário desconhecido\n");
        RAISE Terminate 
      END;

      username := M3toC.StoT(user);
      username := username & 
                  Text.Sub(blanks,0,Text.Length(blanks)-Text.Length(username));
    EXCEPT
    | Abort => Wr.PutText(stderr,"\nProblemas com utilitários\n");
               RAISE Terminate
    END
  END Initialize;
  
  PROCEDURE GetOptions() RAISES{Terminate} =
  <* FATAL Wr.Failure *>
  BEGIN
      TRY
        ParseParams.BeginParsing(stderr);
        IF ParseParams.KeywordPresent("-help") THEN RAISE Scan.BadFormat END;
        LocalFlag1 := ParseParams.KeywordPresent("-1");
        IF LocalFlag1 THEN LocalFileName1 := ParseParams.GetNext() END;
        LocalFlag2 := ParseParams.KeywordPresent("-2");
        IF LocalFlag2 THEN LocalFileName2 := ParseParams.GetNext() END;
        LocalFlag3 := ParseParams.KeywordPresent("-3");
        IF LocalFlag3 THEN LocalFileName3 := ParseParams.GetNext() END;
        LocalFlag4 := ParseParams.KeywordPresent("-4");
        IF LocalFlag4 THEN LocalFileName4 := ParseParams.GetNext() END;
        KeywordFlag := ParseParams.KeywordPresent("-keyword");
        IF KeywordFlag THEN KeywordFileName := ParseParams.GetNext() END;
        RemoteFlag := ParseParams.KeywordPresent("-remote");
        IF RemoteFlag THEN remote := ParseParams.GetNext() END;
        TransposeFlag := ParseParams.KeywordPresent("-transpose");
        RemoveFlag := ParseParams.KeywordPresent("-remove");
        SubstFlag := ParseParams.KeywordPresent("-subst");
        InsertFlag := ParseParams.KeywordPresent("-insert");
        IF NOT ParseParams.KeywordPresent("-s") THEN RAISE Scan.BadFormat END;
        site := ParseParams.GetNextInt(DCC,IME);
        IF NOT ParseParams.KeywordPresent("-m") THEN RAISE Scan.BadFormat END;
        mode := ParseParams.GetNextInt(BATCH,OPENWIN);
        AllWrongFlag := ParseParams.KeywordPresent("-allwrong");
        IF mode=OPENWIN AND AllWrongFlag THEN RAISE Scan.BadFormat END;
        IF ParseParams.KeywordPresent("-t") THEN 
          TRY
            WITH
              tf=FileStream.OpenRead(ParseParams.GetNext())
            DO
              IF Text.Equal(TEST,Rd.GetLine(tf)) THEN
                TestFlag := TRUE;
                IF ParseParams.KeywordPresent("-d") THEN
                  AutFileName := ParseParams.GetNext()
                END
              ELSE
                RAISE Terminate
              END
            END
          EXCEPT
            Terminate,
            Rd.EndOfFile,
            Rd.Failure => Wr.PutText(stderr,"\nUso indevido do programa\n");
                          RAISE Terminate
          END
        END;
        IF mode=BATCH  THEN
          AltFlag := ParseParams.KeywordPresent("-alt")
        END;
        IF RemoteFlag AND (site#DCC OR mode#BATCH) THEN
          Wr.PutText(stderr,"\nUso indevido do programa\n");
           RAISE Terminate
        END;
        ParseParams.EndParsing();
      EXCEPT
        Rd.Failure,
        Scan.BadFormat  => 
            Wr.PutText(stderr,"\nUso:\n");
            Wr.PutText(stderr, "     ProgramName   [-1 LocDic1] [-2 LocDic2] [-3 LocDic3] [-4 LocDic4]\n");
            Wr.PutText(stderr, "                   [-keyword] [-transpose] [-remove] [-subst] [-insert] [-alt]\n");
            Wr.PutText(stderr, "                   -s site -m mode [remote String] [-allwrong] [-t TestFile] [-d automatondir]\n");
            RAISE Terminate
      END;
  END GetOptions;

  PROCEDURE SetUp() RAISES{Terminate} =
  <* FATAL Wr.Failure *>
  BEGIN
      IF TestFlag THEN
        Wr.PutText(stderr,"\n" & "Automaton: " & AutFileName &
                          "\n" & "Local1: " & LocalFileName1 &
                          "\n" & "Local2: " & LocalFileName2 &
                          "\n" & "Local3: " & LocalFileName3 &
                          "\n" & "Local4: " & LocalFileName4 &
                          "\n" & "Keywords: " & KeywordFileName & "\n")
      END;

      IF mode#OPENWIN AND NOT RemoteFlag THEN
        Wr.PutText(stderr,
                      "\nVerificador Ortográfico -- " & Version &
                      "\nUso exclusivo em pesquisa e ensino" &
                      "\nAutores: T.Kowaltowski, C.L.Lucchesi, J.Stolfi");
        Wr.PutText(stderr,"\nLeia LEIAME\n\n")
      END;


      (* This part makes it somwhat more difficult to fiddle with 
      directory names directly in the object file!  *)

      CASE site OF
      
        DCC => pathlen := 12; 
               path[0] :=  '/';
               path[1] :=  'p';
               path[2] :=  'r';
               path[3] :=  'o';
               path[4] :=  'j';
               path[5] :=  '/';
               path[6] :=  'd';
               path[7] :=  'i';
               path[8] :=  'c';
               path[9] :=  'i';
               path[10] := 'o';
               path[11] := '/' 

      | IME => pathlen := 21;
               path[0] :=  '/';
               path[1] :=  'h';
               path[2] :=  'o';
               path[3] :=  'm';
               path[4] :=  'e';
               path[5] :=  '/';
               path[6] :=  's';
               path[7] :=  'p';
               path[8] :=  'e';
               path[9] :=  'c';
               path[10] := 'm';
               path[11] := 'a';
               path[12] := 'c';
               path[13] := '/';
               path[14] := 't';
               path[15] := 'o';
               path[16] := 'm';
               path[17] := 'a';
               path[18] := 's';
               path[19] := 'z';
               path[20] := '/'
               
      ELSE      Wr.PutText(stderr,"\nUso indevido do programa\n");
                RAISE Terminate

      END;

      object[0] := 'V';
      object[1] := 'e';
      object[2] := 'r';
      object[3] := 'i';
      object[4] := 'f';

      PUBPathName := Text.FromChars(SUBARRAY(path,0,pathlen)) & "PUB/";
      notPUBPathName := PUBPathName & "notPUB/";
      
      IF TestFlag THEN
        Wr.PutText(stderr,"\n" & "PUB: " & PUBPathName &
                          "\n" & "notPUB: " & notPUBPathName & "\n")
      END;
      
      TRY
      
        IF TestFlag THEN
          Wr.PutText(stderr,"\nCarregando o autômato");
          Aut  := ReadOnlyPermDAG.LoadCompr(FileStream.OpenRead(AutFileName))
        ELSE
          Aut  := ReadOnlyPermDAG.LoadCompr(FileStream.OpenRead(notPUBPathName & AutFileName))
        END;
      
        Root := Aut.Root();

        (* just to make sure nobody fiddled too much! *)
        
        IF NOT TestFlag THEN
          obrd := FileStream.OpenRead(notPUBPathName & Text.FromChars(object));
          Rd.Close(obrd)
        END;

        LoadWordFile(notPUBPathName & ExtrasFileName,LocalTree,"\nArquivo original não encontrado;"
                               & "\nProvavelmente tentativa de uso indevido\n");

        LoadWordFile(notPUBPathName & IgnoreFileName,IgnoreTree,"\nArquivo original não encontrado;"
                               & "\nProvavelmente tentativa de uso indevido\n");

        IF KeywordFlag THEN
          LoadWordFile(KeywordFileName,KeywordTree,"\nProblemas com arquivo de palavras reservadas\n")
        END;

        LoadLocalDic(LocalFileName1,LocalFlag1);
        LoadLocalDic(LocalFileName2,LocalFlag2);
        LoadLocalDic(LocalFileName3,LocalFlag3);
        LoadLocalDic(LocalFileName4,LocalFlag4);

      EXCEPT
      
        Rd.Failure   => Wr.PutText(stderr,"\nArquivos originais não encontrados;"
                                & "\nProvavelmente tentativa de uso indevido\n");
                        RAISE Terminate
      | BadVersion,
        Rd.EndOfFile => Wr.PutText(stderr,"\nProblemas com o autômato.\n");
                        RAISE Terminate

      END;
    
  END SetUp;

  PROCEDURE Check() RAISES{Terminate} =
  BEGIN
    CASE mode OF
      BATCH    => CheckBatch()
    | OPENWIN  => CheckOpen()
    END
  END Check;

  PROCEDURE CheckBatch() RAISES{Terminate} =
  <* FATAL Wr.Failure *>
  VAR
    txt: TEXT;
  BEGIN
    IF TestFlag THEN
      Wr.PutText(stderr,"\nEntering BATCH mode\n")
    END;
    TRY
      LOOP
        txt := NextWord();
        IF Text.Length(txt)>0 THEN
          INC(WordCount);
          IF AllWrongFlag OR 
             NOT ((Aut.Accepts(Root,txt) AND (NOT IgnoreTree.Search(txt))) 
                  OR  LocalTree.Search(txt) OR KeywordTree.Search(txt)) THEN
            INC(NotFoundCount);
            ErrorTree.Insert(txt)
          END
        END
      END
    EXCEPT
      Done  =>
            TRY
              ErrorTree.Enum(BatchAction);
              Wr.Flush(stdout);
              RETURN
            EXCEPT
              Wr.Failure,
              Abort => Wr.PutText(stderr,"\nProblemas na enumeração\n");
                            RAISE Terminate
              END
    END
  END CheckBatch;

  PROCEDURE CheckOpen() RAISES {Terminate} = 
  <* FATAL Wr.Failure *>
  TYPE
    Status = {ONE,MANY,INSERT,UNDEF};
  VAR
    txt: TEXT;
    stat: Status := Status.UNDEF;
  BEGIN
    TRY
      IF TestFlag THEN
        Wr.PutText(stderr,"\nEntering OPENWIN mode\n")
      END;
      LOOP
        TRY 
          LOOP 
            txt := NextWord();
            IF Text.Equal(txt,CodeOne) THEN
              stat := Status.ONE; EXIT
            ELSIF Text.Equal(txt,CodeMany) THEN
              stat := Status.MANY; EXIT
            ELSIF Text.Equal(txt,CodeInsert) THEN
              stat := Status.INSERT; EXIT
            END
          END;
          CASE stat OF
            Status.UNDEF => <* ASSERT stat#Status.UNDEF *>  (* Impossible *)
          | Status.ONE  =>
                IF TestFlag THEN
                  Wr.PutText(stderr,"\nEntering ONE word\n")
                END;
                txt := NextWord();
                IF Text.Length(txt)>0 THEN
                  INC(WordCount);
                  IF AllWrongFlag OR 
                     NOT ((Aut.Accepts(Root,txt) AND (NOT IgnoreTree.Search(txt))) 
                           OR  LocalTree.Search(txt)) THEN
                    INC(NotFoundCount);
                    Wr.PutChar(stdout,'0')
                  ELSE
                    Wr.PutChar(stdout,'1')
                  END;
                  Adviser.Alternatives(txt,Aut,LocalTree,AltTree,TransposeFlag,RemoveFlag,SubstFlag,InsertFlag,FALSE);
                  CleanTree(AltTree,IgnoreTree);
                  EnumAlternatives(AltTree,txt);
                  Wr.PutChar(stdout,'#')
                ELSE  (* interactive version always wants an answer! *)
                  Wr.PutText(stdout,"0#")
                END;
                stat := Status.UNDEF;
(*                IF TestFlag THEN  *)
                  Wr.Flush(stdout)
(*                END   *)
          | Status.MANY  =>
                IF TestFlag THEN
                  Wr.PutText(stderr,"\nEntering MANY words\n")
                END;
                INC(FileCount);
                LOOP
                  txt := NextWord();
                  IF Text.Equal(txt,"#") THEN EXIT END;
                  IF Text.Length(txt)>0 THEN
                    INC(WordCount);
                    IF AllWrongFlag OR 
                       NOT ((Aut.Accepts(Root,txt) AND (NOT IgnoreTree.Search(txt))) 
                             OR  LocalTree.Search(txt) OR KeywordTree.Search(txt)) THEN
                      INC(NotFoundCount);
                      ErrorTree.Insert(txt)
                    END
                  END
                END;
                ErrorTree.Enum(ManyAction);
                ErrorTree := AVLTextTree.New();
                Wr.PutChar(stdout,'#');
                Wr.Flush(stdout);
                stat := Status.UNDEF
          | Status.INSERT  =>
                IF TestFlag THEN
                  Wr.PutText(stderr,"\nEntering INSERT\n")
                END;
                txt := NextWord();
                IF Text.Length(txt)>0 THEN
                  LocalTree.Insert(txt)
                END;
                stat := Status.UNDEF
          END (* CASE stat *)
        EXCEPT
            Done  => IF stat#Status.UNDEF THEN
                       Wr.PutText(stderr,"\nProblemas com status\n");
                       RAISE Terminate 
                     ELSE 
                       RETURN 
                     END
        END (* TRY *)
      END (* LOOP *)
    EXCEPT
         Wr.Failure => Wr.PutText(stderr,"\nProblemas na saída\n");
                       RAISE Terminate
      |  Abort => Wr.PutText(stderr,"\nProblemas nas alternativas\n");
                  RAISE Terminate
    END
  END CheckOpen;
 
  PROCEDURE LoadWordFile(lf: TEXT; VAR (*IO*) t: AVLTextTree.T; mess: TEXT) 
                        RAISES {Terminate} =
  <* FATAL Wr.Failure *>                        
  BEGIN
    IF TestFlag THEN
      Wr.PutText(stderr,"\nCarregando: " & lf & "\n")
    END;
    TRY
      WITH lfrd=FileStream.OpenRead(lf) DO
        TRY
          LOOP
            WITH
              line = Util.Trim(Rd.GetLine(lfrd))
            DO
              IF Text.Length(line)>0 THEN 
                t.Insert(line);
              END
            END
          END
        EXCEPT
          Rd.EndOfFile  => Rd.Close(lfrd);
                           IF TestFlag THEN
                             PrintTree(stderr,t)
                           END;
                           RETURN
        END
      END
    EXCEPT
      Rd.Failure  =>
                Wr.PutText(stderr,"\n" & mess & "\n");
                RAISE Terminate
    END;
  END LoadWordFile;
  
  PROCEDURE PrintTree(trwr:Wr.T; t: AVLTextTree.T)  =
  <* FATAL Abort *>
    PROCEDURE Action(READONLY w: TEXT) =
    <* FATAL Wr.Failure *>
    BEGIN
      Wr.PutText(trwr,w & "\n");
    END Action;
  BEGIN
    t.Enum(Action)
  END PrintTree;
  
  PROCEDURE LoadLocalDic(fn: TEXT; fl: BOOLEAN) RAISES {Terminate} = 
  BEGIN
    IF fl THEN
      LoadWordFile(fn,LocalTree,"\nProblemas com arquivos locais\n")
    END
  END LoadLocalDic;
  
  PROCEDURE NextWord(): TEXT RAISES {Done, Terminate} =
  <* FATAL Wr.Failure *>
  BEGIN
    TRY
      RETURN Util.Trim(Rd.GetLine(stdin))
    EXCEPT
      Rd.EndOfFile   =>  RAISE Done
    | Rd.Failure     => 
              Wr.PutText(stderr,"\nNão consegue ler um arquivo\n\n");
              RAISE Terminate
    END
  END NextWord;

  PROCEDURE BatchAction(READONLY w: TEXT) RAISES {Abort}  =
  <* FATAL Wr.Failure *>    
    PROCEDURE AltAction(READONLY u: TEXT) RAISES {Abort} =
    BEGIN
      TRY
        Wr.PutText(stdout,u);
        Wr.PutChar(stdout,' ')
      EXCEPT
        Wr.Failure => Wr.PutText(stderr,"\nProblemas na escrita de alternativas\n");
                      RAISE Abort 
      END
    END AltAction;
  BEGIN (* BatchAction *)
    TRY
      INC(ErrorCount);
      Wr.PutText(stdout,w);
      IF AltFlag THEN
        Adviser.Alternatives(w,Aut,LocalTree,AltTree,TransposeFlag,RemoveFlag,SubstFlag,
                             InsertFlag,FALSE);
        CleanTree(AltTree,IgnoreTree);
        IF AltTree.Size()>0 THEN
          Wr.PutText(stdout," ( ");
          AltTree.Enum(AltAction);
          Wr.PutChar(stdout,')')
        END;
      END;
      Wr.PutChar(stdout,'\n')
    EXCEPT
      Wr.Failure,
      Abort  => Wr.PutText(stderr,"\nProblemas na busca ou escrita de alternativas\n");
                RAISE Abort
    END
  END BatchAction;

  PROCEDURE EnumAlternatives(tr: AVLTextTree.T; txt: TEXT)  =
  <* FATAL Abort *>
    PROCEDURE Action(READONLY w: TEXT ) =
    <* FATAL Wr.Failure *>
    BEGIN
      IF NOT Text.Equal(txt,w) THEN
        Wr.PutText(stdout,w & "\n")
      END
    END Action;
  BEGIN
    tr.Enum(Action)
  END EnumAlternatives;
  
  PROCEDURE CleanTree(VAR tr: AVLTextTree.T; deltr: AVLTextTree.T) =
  <* FATAL Abort *>
    PROCEDURE Action(READONLY w: TEXT) =
    BEGIN
      tr.Delete(w)
    END Action;
  BEGIN
    deltr.Enum(Action)  
  END CleanTree;
  
  PROCEDURE ManyAction(READONLY w: TEXT)  RAISES{Abort} =
  <* FATAL Wr.Failure *>
  BEGIN
    TRY
      INC(ErrorCount);
      Wr.PutText(stdout,w & "\n")
    EXCEPT
      Wr.Failure,
      Abort  => Wr.PutText(stderr,"\nProblemas na busca ou escrita de alternativas\n");
                RAISE Abort
    END
  END ManyAction;
  
  PROCEDURE Finalize() RAISES{Terminate} =
  <* FATAL Wr.Failure *>
  BEGIN
    TRY
      ExecTime := Util.ExecTimesTexts().TotalTime;

      CASE mode OF
        BATCH   => Wr.PutText(stderr,"\n" & Fmt.Int(WordCount) & " palavras lidas; " &
                          Fmt.Int(ErrorCount) & " palavras não encontradas (" &
                          Fmt.Int(NotFoundCount) & " ocorrências)\n" & ExecTime & "s\n");
                   IF RemoteFlag THEN 
                     LogFileName := RLOG;
                     username := Fmt.Pad(text:=remote,length:=50,align:=Fmt.Align.Left)
                   ELSE 
                     LogFileName := LOG 
                   END;
                   
      | OPENWIN => LogFileName := XLOG
      END;

      FinalTime := Util.DateTime();

      IF NOT TestFlag THEN
        logwr := FileStream.OpenAppend(notPUBPathName & LogFileName);

        Wr.PutText(logwr, username & InitTime & "  " & FinalTime &
                          " " & ExecTime & 
                          Fmt.Pad(Fmt.Int(WordCount),9) & 
                          Fmt.Pad(Fmt.Int(NotFoundCount),9) & 
                          Fmt.Pad(Fmt.Int(ErrorCount),8) & Fmt.Pad(Fmt.Int(FileCount),3) & "\n");
        Wr.Close(logwr)
      END;

      Wr.Flush(stderr);

    EXCEPT

      Abort      => Wr.PutText(stderr,"\nProblemas com utilitários\n");
                    RAISE Terminate
    | Wr.Failure => Wr.PutText(stderr,"\nProblemas na saída final\n");
                    RAISE Terminate
    END
    
  END Finalize;
  
BEGIN 

  Main()

END Verif.

