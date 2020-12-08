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

MODULE IdealClasses EXPORTS Main;      (* Version 2.3 *)

(* HISTORY: (T.K.)

18.Aug.92: Redone, starting from the older version FMatchClasses,
           for compatibility with version 2.04 of Modula-3 and the 
           new interfaces 
03.Sep.92: Minor adjustments and additional options  
14.Sep.92: Three separate criteria for matching 
09.Apr.93: Modified input format for classes 
23.Apr.93: Added option -representative 
08.Jan.94: Started working again -- fixing bugs 
12.Jan.94: Apparently works!

*)              

(* Program to identify  states close to  given classes.  These classes
are pre-established and not necessarily  those produced by the program
'AutoAnalysis'.   Classes description   is  given  as  an  input  file
(ClassesFile) specifying the set of endings for each class; its format
is described  below. Each  state  is matched against  all classes  and
identified  with a class  subset  (within the subset, all  classes are
disjoint), according to criteria described below.


------------------------------- FILE FORMAT --------------------------

Rmk.:  Any line starting  with a  '#'   before a class description  is
considered a comment.

----------------------------------------------------------------------

Class <class name>: <mcd1> <mcd2> (<repr>)
{ <ending>, <ending>, ..., <ending>}

Class <class name>: <mcd1> <mcd2> (<repr>)
{ <ending>, <ending>, ..., <ending>}

...

Class <class name>: <mcd1> <mcd2> (<repr>)
{ <ending>, <ending>, ..., <ending>}

-----------------------------END FILE FORMAT --------------------------

where:  <class  name> is any string not   containing ':';  <mcd1> is a
cardinal  denoting the maximum number  of endings  that can be missing
for a given state to be considered matched to a given class; <mcd2> is
a cardinal  denoting the maximum number  of endings corresponding to a
state and not in the class;  <repr> is  a representative ending of the
class and must occur among the endings -- its presence is optional and
required only  if the  -representative  option is   specified. Newline
characters are  ignored.   The maximum  number  of classes is  64 (see
constant "MaxNClasses").

Given  the  classes g0,g1,...,gn (i.e.  sets of  endings), the program
determines, for each state s, a  subset G(s) of  these classes through
the following criteria:

        Ideal(s) = SG(s) union Rest(s)
        
        Rest(s) = UNION a & Ideal(s*a)
                    a

  where  SG(s) =  union   g
                 g in G(s)
      
         and s*a is the state reached from s by the letter a (if it exists).
         
  A set (class) g is included in G(s) only if:

    (1)    | g\Suff(s) | <= MaxDiff1(g) |    and
        
    (2)    | (Suff(s)\Rest(s))\g | <= MaxDiff2(g) and
    
    (3)    | Suffs(s)\Rest(s)\  UNION   g |   <=  MaxDiff
                               g in G(s) 
    (4)    | Suffs(s)\Rest(s)\  UNION   g |   <=  MaxRelDiff*|Suffs(s)|
                               g in G(s) 

Each criterion can be selected independently through  command  line  options:
"-criterion1",   "-criterion2",   "-criterion3   MaxDiff"   and  "-criterion4
MaxRelDiff".

The program produces an output file (MatchesFile) which is a listing of
words preceded by '*', '+' or '-':

        *P word      : the word belongs to both the input automaton and to
                       Ideal(Root) matching an ending from class P
        +PQ... word  : the word belongs to the input automaton but not to 
                       Ideal(Root); the last matched state within the 
                       automaton matches classes PQ... (in this case a
                       hyphen separates the matching prefix from the ending);
        -P word      : the word belongs to  Ideal(Root) with an ending
                       from class  P but not to the input automaton

where  P,Q,   ...   are  the  <class  prefixes>  (this  are characters
associated by the program  with each  class  starting with "A"  -- see
class listings in  the output messages file.   The  output of  each of
these  sets  of words is  controlled  through appropriate command line
options. 

Example  of  a  classes file  for   the   Portuguese regular (at least
phonetically) verbs. Notice that  the first class denotes the standard
first regular conjugation, whereas other prefixes denote its variants.

------------------------------ EXAMPLE ---------------------------------

#
# First conjugation with some variants
# 

Class First conjugation -ar: 10 5

{ o, as, a,  amos, ais, am, ava, avas,  ava, ávamos, áveis,  avam, ei,
aste, ou, amos,  astes, aram, ara,  aras,  ara, áramos,  áreis,  aram,
arei,  arás, ará, aremos,  areis, arão, aria,  arias,  aria,  aríamos,
aríeis, ariam, e,  es, e, emos, eis,  em,  asse, asses, asse, ássemos,
ásseis, assem, ar,  ares, ar, armos, ardes, arem,  e, a, e,  emos, ai,
em, e, es, e, emos, eis, em, ando, ado, ados, ada, adas }

Class First conjugation variant -gar: 5 1

{ go, gas,  ga, gamos, gais, gam,  gava, gavas, gava, gávamos, gáveis,
gavam, guei,  gaste,   gou, gamos,  gastes,  garam, gara, garas, gara,
gáramos, gáreis, garam, garei, garás,  gará,  garemos, gareis,  garão,
garia, garias,   garia, garíamos,  garíeis,  gariam, gue,  gues,  gue,
guemos, gueis, guem, gasse, gasses,  gasse, gássemos, gásseis, gassem,
gar, gares, gar, garmos, gardes,   garem, gue, ga, gue, guemos,   gai,
guem, gue, gues,  gue, guemos, gueis, guem,  gando, gado, gados, gada,
gadas }

Class First conjugation variant -car: 5 1

{ co, cas, ca, camos,  cais, cam, cava,  cavas, cava, cávamos, cáveis,
cavam, quei,  caste, cou,  camos,  castes,  caram,  cara, caras, cara,
cáramos, cáreis, caram,  carei, carás, cará, caremos,  careis,  carão,
caria, carias,  caria,  caríamos,  caríeis,  cariam, que,  ques,  que,
quemos, queis, quem, casse, casses,  casse, cássemos, cásseis, cassem,
car,  cares,  car, carmos, cardes, carem,  que, ca,  que, quemos, cai,
quem, que, ques, que, quemos, queis,  quem, cando,  cado, cados, cada,
cadas }




----------------------------END EXAMPLE ---------------------------------

*)

(*

COMMAND LINE PARAMETERS (with defaults whenever applicable):
============================================================

-help                : prints the help  message
-mess  FileName      : messages file name (default: stdout)
-dump FileName       : file containing the automaton to be analyzed
-classes FileName    : file for classes description data (default: stdin)
-matches FileName    : file for results of matching (default: stdout)
-criterion1          : apply criterion (1) for matching (default: off)
-criterion2          : apply criterion (2) for matching (default: off)
-criterion3 NAT      : apply criterion (3) for matching (default: no limit)
-criterion4 NAT      : apply criterion (4) for matching (default: 100 (%))
-inter               : include in program output words marked by * 
                       (intersection)
-plus                : include in program output words marked by + (those
                       not in the ideal)
-minus               : include in program output words marked by - (those 
                       not in the automaton)
-representative      : under the -inter and -minus options inlude only 
                       words terminating with the representative ending of 
                       the class (default: FALSE);
-report NAT          : report interval for some of the more time consuming
                       procedures (default: no limit)
-showmatches         : shows classes matched by each state (default: FALSE)
-radicals FileName   : file for radicals, ie, sets of prefixes of all the 
                       states that have been matched against at least one class 
                       (including the class information) (default: no output)
     
(Actually, "no limit" means LAST(NAT).)

*)


  IMPORT Rd, Wr, Text, Fmt, FileRd, FileWr, OSError, Scan, ParseParams,
         Thread, Basics, Reduced, ReducedPair, Encoding, Lex, Char, Convert,
         PlainEncoding, StringPrinter, Util;
         
  FROM Basics IMPORT BOOL, NAT, Letter, String, Abort, Done, CopyString, WrNat; 
  FROM Stdio IMPORT stderr,stdin, stdout; 
  FROM Reduced IMPORT State, NullState, UnitState;
  
  EXCEPTION Terminate; MyAbort(TEXT);
  
  <* FATAL Thread.Alerted *>
  
  CONST
  
    Version = "Version 2.3";
    MaxLen = 100;
    MaxNClasses = 64;
    NoMatch: StateClasses = NIL;
    ClassAutDagSize = 50000;
  
    ClassPref = ARRAY [0..MaxNClasses-1] OF CHAR{
        'A','B','C','D','E','F','G','H','I','J','K','L','M',
        'N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
        'a','b','c','d','e','f','g','h','i','j','k','l','m',
        'n','o','p','q','r','s','t','u','v','w','x','y','z',
        '0','1','2','3','4','5','6','7','8','9','$','&'};
                
  TYPE
    MatchRecord = RECORD
        class: [0..MaxNClasses-1];
        next: REF MatchRecord
      END;
    StateClasses = REF MatchRecord;
    
    StringRecord = RECORD
        string:  REF String;
        next: REF StringRecord
      END;
      
    ClassStateRecord = RECORD
        state: State;
        next: REF ClassStateRecord
      END;
  
  VAR

    meswr: Wr.T := stderr;
    classrd: Rd.T := stdin;
    matchwr: Wr.T := stdout;
    dummywr: Wr.T;

    mespr: StringPrinter.T;
    
    encoding: Encoding.T := PlainEncoding.New();

    MessFlag,
    DumpFlag,
    ClassesFlag,
    RadFlag,
    Crit1Flag,
    Crit2Flag,
    InterFlag,
    PlusFlag,
    MinusFlag,
    ReprFlag,
    ShowFlag,
    MatchesFlag: BOOL;

    MessFileName,
    DumpFileName,
    ClassesFileName,
    RadFileName,
    MatchesFileName: TEXT;

    ClassAut: Reduced.T := Reduced.New(ClassAutDagSize);
    Aut: Reduced.T;
    ACPair: ReducedPair.T;
    
    Root: State;
  
    NClasses: NAT;
    
    Classes := ARRAY [0..MaxNClasses-1] OF State{NullState, ..};
    SubClasses := ARRAY [0..MaxNClasses-1] OF 
                    ARRAY [0..MaxNClasses-1] OF 
                         BOOL{ARRAY [0..MaxNClasses-1] OF BOOL{FALSE, ..}, ..};
    nSClasses: ARRAY [0..MaxNClasses-1] OF NAT;
    ClassStrings: REF ARRAY OF REF StringRecord;
    MaxDiff1: ARRAY [0..MaxNClasses-1] OF NAT;
    MaxDiff2: ARRAY [0..MaxNClasses-1] OF NAT;
    ClassReprString: ARRAY [0..MaxNClasses-1] OF REF String;
    ClassComment: ARRAY [0..MaxNClasses-1] OF TEXT;
    BestClasses: REF ARRAY OF StateClasses;
    ClassStates := ARRAY [0..MaxNClasses-1] OF REF ClassStateRecord{NIL, ..};

    MaxDiff: NAT := LAST(NAT);
    MaxRelDiff: NAT := 100;  (* % *)
    NMStates: NAT := 0;
    NMatches: NAT :=0;
    MatchCounts := ARRAY [0..MaxNClasses-1] OF NAT{0, ..};
    ClassMatchCounts := ARRAY [0..MaxNClasses-1] OF NAT{0, ..};
    
    CurrentClassPref: CHAR;
    CurrentClassReprString: REF String;
    
    ReportInterval: NAT := LAST(NAT);
    BestMatchReportCount,
    GenIdealCount,
    GenIdealComplementCount: NAT := 0;
    
(* ---------------------------------------------------------------------- *)

  PROCEDURE Main() = 
  <* FATAL Wr.Failure, Abort *>
  BEGIN
    
    TRY
      EVAL Util.ExecTimesTexts();
      
      GetOptions();
      
      mespr := StringPrinter.New(wr:=meswr,
                                 encoding:=encoding,
                                 leftMargin:=4);
      dummywr :=  FileWr.Open("DummyFile");

      GetAutomaton();
      
      GetClasses();
      
      ComputeSubClasses();
      
      BestClasses := NEW(REF ARRAY OF StateClasses, Root+1);
      
      ACPair := ReducedPair.New(Aut,ClassAut);
      
      FindBestMatches();
      
      GenIdeal();
      
      GenIdealComplement();
      
      GenRadicals();
      
      Wr.Flush(matchwr);
      Wr.PutText(meswr,"\n\nEnd of processing\n\n" & Util.DateTime() & "\n");
      Wr.Close(meswr);
      Wr.Flush(dummywr);

    EXCEPT
       Wr.Failure  =>  Wr.PutText(meswr,"\nProblems in writing a file.\n")
     | Rd.Failure,
       Rd.EndOfFile=>  Wr.PutText(meswr,"\nProblems in reading a file.\n")
     | Abort       =>  Wr.PutText(meswr,"\nProblems with utility functions.\n")
     | Terminate   =>  Wr.PutText(meswr,"\n\nEnd of processing\n\n" & Util.DateTime() & "\n");
                       Wr.Close(meswr);
    END;
    
  END Main;


  PROCEDURE GetOptions() RAISES {Terminate,Wr.Failure,Rd.Failure,Abort} =
  BEGIN
    TRY
      WITH pp = NEW(ParseParams.T).init(stderr) DO
      IF pp.keywordPresent("-help") THEN RAISE Scan.BadFormat END;
      Util.GetFlagAndTextParam(meswr,"-mess",MessFlag,MessFileName);
      IF MessFlag THEN
        meswr := FileWr.Open(MessFileName)
      END;
      Wr.PutText(meswr,"\n" & Util.DateTime() & "\n\n");
      Wr.PutText(meswr,
          "\nProcessing \'IdealClasses\' (" & Version & 
          ") with the following options:");

      Wr.PutText(meswr,
          "\n==================================================================\n");
      Util.GetFlagAndTextParam(meswr,"-dump",DumpFlag,DumpFileName);
      Util.GetFlagAndTextParam(meswr,"-classes",ClassesFlag,ClassesFileName);
      Util.GetFlagAndTextParam(meswr,"-matches",MatchesFlag,MatchesFileName);
      Util.GetFlagAndTextParam(meswr,"-radicals",RadFlag,RadFileName);
      Util.GetFlagParam(meswr,"-criterion1",Crit1Flag);
      Util.GetFlagParam(meswr,"-criterion2",Crit2Flag);
      Util.GetNatParam(meswr,"-criterion3",MaxDiff);
      Util.GetNatParam(meswr,"-criterion4",MaxRelDiff);
      Util.GetFlagParam(meswr,"-inter",InterFlag);
      Util.GetFlagParam(meswr,"-plus",PlusFlag);
      Util.GetFlagParam(meswr,"-minus",MinusFlag);
      Util.GetFlagParam(meswr,"-representative",ReprFlag);
      Util.GetFlagParam(meswr,"-showmatches",ShowFlag);
      Util.GetNatParam(meswr,"-report",ReportInterval);
      pp.finish();
      Wr.PutChar(meswr,'\n');
      IF NOT DumpFlag   THEN
        Wr.PutText(stderr,"Missing automaton (-dump).\n");
        RAISE Scan.BadFormat
      END;
      IF  ClassesFlag  THEN
        classrd := FileRd.Open(ClassesFileName)
      END;
      IF MatchesFlag   THEN
        matchwr := FileWr.Open(MatchesFileName)
      END;
      IF (NOT Crit1Flag) AND (NOT Crit2Flag) 
          AND (MaxDiff=LAST(NAT)) AND (MaxRelDiff=100) THEN
        Wr.PutText(meswr,"\n\nWARNING: No matching criterion specified!\n\n")
      END;
      ProcessingTime();
      Wr.Flush(meswr)
    EXCEPT
      Scan.BadFormat  =>
           Wr.PutText(stderr, "\nUsage:\n");     
           Wr.PutText(stderr, "     ProgramName   -dump DumpFile [-mess MessageFileName (stderr)]\n");
           Wr.PutText(stderr, "                   [-classes ClassesFileName (stdin)]\n");
           Wr.PutText(stderr, "                   [-matches MatchesFileName (stdout)] \n");
           Wr.PutText(stderr, "                   [-radicals RadicalsFileName (none)] \n");
           Wr.PutText(stderr, "                   [-inter (FALSE)] [-plus (FALSE)] \n");
           Wr.PutText(stderr, "                   [-minus (FALSE)] [-representative (FALSE)]\n");
           Wr.PutText(stderr, "                   [-criterion1 (FALSE)] [-criterion2 (FALSE)]\n");
           Wr.PutText(stderr, "                   [-criterion3 (no limit)] [-criterion4 (100%)]\n");
           Wr.PutText(stderr, "                   [-showmatches (FALSE)] [-report (no limit)]\n");
           RAISE Terminate
    | Wr.Failure      =>
           Wr.PutText(stderr,"\nProblems with the message file\n");
           RAISE Terminate
    END
  END GetOptions;

  
  PROCEDURE GetAutomaton() RAISES {Rd.Failure,Wr.Failure,Abort} =
  BEGIN 
     PhaseStartMessage("Reading the automaton");
     Aut := Reduced.Load(FileRd.Open(DumpFileName));
     Wr.PutText(meswr,"\nAutomaton loaded\nComment:\n" & Aut.doc & "\n");
     Root := Aut.Root();
     ProcessingTime()
  END GetAutomaton;


  PROCEDURE GetClasses() RAISES {Terminate,Wr.Failure} =
  <* FATAL Abort *>

     CONST
       empty: TEXT = "";
       commentchars = SET OF CHAR{'#'};
       classword: TEXT = "class";
       Classword: TEXT = "Class";
       classNameChars = Char.All - SET OF CHAR{':'};
       skipChars = Lex.Blanks + SET OF CHAR{':'};
       openClassChars = SET OF CHAR{'{'};
       closeClassChars = SET OF CHAR{'}'};
       separators = Lex.Blanks + SET OF CHAR{','};
       openReprChars = SET OF CHAR{'('};
       closeReprChars = SET OF CHAR{')'};
       endingChars = Char.All - (separators + closeClassChars+closeReprChars); 
                           (* Allow any accents! *)
       
       
     VAR 
       token: TEXT;
       classNumber: CARDINAL := 0;
       classPrefix: CHAR;
       classPrefixTxt: TEXT;
       reprLen: NAT;
       reprString := NEW(REF String, 100);
       
     PROCEDURE NextString (* Reduced.NextStringProc *)
                         (VAR (* IO *) s: REF String; 
                          VAR (* OUT *) len: NAT;
                          VAR (* OUT *) add: BOOL) RAISES {Done,Abort} =
     <* FATAL Wr.Failure *>
       VAR
         line: TEXT;
         
     BEGIN
       add := TRUE;
       TRY
         Lex.Skip(classrd,separators);
         line := Lex.Scan(classrd,endingChars);
         IF Text.Equal(line,empty) THEN
            RAISE Done
         END;
         encoding.TextToString(classPrefixTxt & line,s,len)
       EXCEPT
       | Encoding.BadChar(ch) =>
                 Wr.PutText(meswr,"\nIllegal classes input char:" & 
                 Text.FromChar(ch) & "\n\n");
                 RAISE Abort
       | Encoding.BadChars =>
                 Wr.PutText(meswr,"\nIllegal classes input chars\n\n");
                 RAISE Abort
       | Rd.Failure,
         Rd.EndOfFile => 
                 Wr.PutText(meswr,"\nCannot read classes input file\n\n");
                 RAISE Abort
       END
     END NextString;
       
   BEGIN (* GetClasses *)
     PhaseStartMessage("Getting classes");
     TRY
       LOOP
         LOOP
           token := Lex.Scan(classrd,commentchars);
           IF  Text.Equal(token,empty) THEN EXIT END;
             EVAL Rd.GetLine(classrd)
           END;
         Lex.Skip(classrd,Lex.Blanks);
         token := Lex.Scan(classrd,Char.Letters);
         IF Text.Equal(token,empty) THEN
           RAISE Done
         ELSIF NOT (Text.Equal(token,classword) OR Text.Equal(token,Classword)) THEN
           RAISE MyAbort("Wrong class identification")
         END;
         ClassComment[classNumber] := Lex.Scan(classrd,classNameChars);
         classPrefix := ClassPref[classNumber];
         classPrefixTxt := Text.FromChar(classPrefix);
         Lex.Skip(classrd,skipChars);
         MaxDiff1[classNumber] := Lex.Int(classrd);
         Lex.Skip(classrd,separators);
         MaxDiff2[classNumber] := Lex.Int(classrd);
         Lex.Skip(classrd,Lex.Blanks);
         ClassReprString[classNumber] := NIL;
         token := Lex.Scan(classrd,openReprChars);
         IF ReprFlag AND Text.Equal(token,empty) THEN
           RAISE MyAbort("Missing representative ending") 
         ELSIF NOT Text.Equal(token,empty)  THEN
           encoding.TextToString(Lex.Scan(classrd,endingChars),
                                 reprString,reprLen);
           ClassReprString[classNumber] := NEW(REF String,reprLen);
           WITH crs=ClassReprString[classNumber]^ DO
             FOR i:=0 TO reprLen-1 DO crs[i] := reprString^[i] END
           END;
           IF Text.Equal(Lex.Scan(classrd,closeReprChars),empty) THEN 
             RAISE MyAbort("Representative ending not properly closed")
           END
         END;
         Lex.Skip(classrd,Lex.Blanks);
         IF Text.Equal(Lex.Scan(classrd,openClassChars),empty) THEN 
           RAISE MyAbort("Class not properly open")
         END;
         ClassAut.Build(next:=NextString,
                        wr:=dummywr,
                        reportInterval:=LAST(NAT),
                        flagRedundant:=FALSE);
         IF Text.Equal(Lex.Scan(classrd,closeClassChars),empty) THEN 
           RAISE MyAbort("Class not properly closed")
         END;
         INC(classNumber)
       END
     EXCEPT
     | Lex.Error       => Wr.PutText(meswr,"Bad input file format: Lex.Error\n");
                          RAISE Terminate
     | Convert.Failed  => Wr.PutText(meswr,"Bad input file format: Convert.failed\n");
                          RAISE Terminate
     | MyAbort(mess)   => Wr.PutText(meswr,mess & "\n");
                          RAISE Terminate
     | Rd.Failure      => Wr.PutText(meswr,"Bad input file format: Rd.Failure\n");
                          RAISE Terminate
     | Rd.EndOfFile    => Wr.PutText(meswr,"Bad input file format: Rd.EndOfFile\n");
                          RAISE Terminate
     | Abort => Wr.PutText(meswr,"Problems in building the class automaton\n");
                RAISE Terminate
     | Encoding.BadChars 
                       => Wr.PutText(meswr,"Bad input file format: Encoding.BadChars\n");
                          RAISE Terminate  
     | Done =>  NClasses := classNumber
     END;

     ClassStrings := NEW(REF ARRAY OF REF StringRecord,NClasses);
 
     WITH 
       cr = ClassAut.Root() 
     DO
       FOR i := 0 TO NClasses-1 DO
         WITH 
           cli = Classes[i],
           clsi = ClassStrings[i]
         DO
         
           PROCEDURE EnumClassSuffsAction(READONLY w: String) RAISES {Abort} =
           BEGIN
             clsi := NEW(REF StringRecord, string:=CopyString(w), next:=clsi)
           END EnumClassSuffsAction;
           
           BEGIN
             TRY
               cli := ClassAut.Succ(cr,ORD(ClassPref[i]));
               clsi := NIL;
               nSClasses[i] := ClassAut.NSuffs(cli);
               ClassAut.EnumSuffs(cli,EnumClassSuffsAction);
               IF cli=NullState OR clsi=NIL THEN 
                 RAISE MyAbort("Problems with the classes automaton.") 
               END;
               IF ReprFlag AND 
                       NOT ClassAut.Accepts(cli,ClassReprString[i]^) THEN
                 RAISE MyAbort("A class representative not in the class")
               END
             EXCEPT
               MyAbort(mess) => 
                 Wr.PutText(meswr, mess & "\n");
                 RAISE Terminate 
             END
           END;
         END;
       END
     END;
 
     FOR c := 0 TO NClasses-1 DO
       PrintClass(meswr,mespr,c)
     END;
     ProcessingTime()
   END GetClasses;


  PROCEDURE ComputeSubClasses() RAISES{Wr.Failure} =
  <* FATAL Abort *>
    PROCEDURE RepSubClass(x,y: [0..MaxNClasses-1]; diff: NAT) 
                        RAISES {Wr.Failure} =
    BEGIN
      Wr.PutText(meswr,"Class ");
      WrNat(meswr,x);
      Wr.PutText(meswr," subclass of class ");
      WrNat(meswr,y);
      Wr.PutText(meswr," (");
      WrNat(meswr,diff);
      Wr.PutText(meswr,")\n")
    END RepSubClass;
  BEGIN
    PhaseStartMessage("Computing SubClasses");
    WITH
      ccp=ReducedPair.New(ClassAut,ClassAut)
    DO
      FOR c := 0 TO NClasses-1 DO
        <* ASSERT NOT SubClasses[c,c] *>
        FOR d := c+1 TO NClasses-1  DO
          TRY
            EVAL ccp.NSuffs(ReducedPair.State{Classes[c],Classes[d]},
                            ReducedPair.BoolOp.Diff,limit:=1);
            SubClasses[c,d] := TRUE;
            RepSubClass(c,d,nSClasses[d]-nSClasses[c])
          EXCEPT
            Abort => <* ASSERT NOT SubClasses[c,d] *>
          END;
          TRY
            EVAL ccp.NSuffs(ReducedPair.State{Classes[d],Classes[c]},
                            ReducedPair.BoolOp.Diff,limit:=1);
            SubClasses[d,c] := TRUE;
            RepSubClass(d,c,nSClasses[c]-nSClasses[d]);
          EXCEPT
            Abort => <* ASSERT NOT SubClasses[d,c] *>
          END
        END        
      END
    END;
    ProcessingTime()
  END ComputeSubClasses;

  
  PROCEDURE FindBestMatches() RAISES{Wr.Failure} =
  <* FATAL Abort *>
  BEGIN
    <* ASSERT BestClasses#NIL *>
    PhaseStartMessage("Computing Best Matches");
    FOR s := UnitState TO Root DO
      IF IsState(s) THEN
        BestClasses[s] := BestMatch(s);
      END;
    END;
    Wr.PutText(meswr,"\nTotal number of matches: " & 
                          Fmt.Int(NMatches));
    Wr.PutText(meswr,"\nTotal of matched states: " &
                          Fmt.Int(NMStates));
    Wr.PutText(meswr,"\n\nMatches per state distribution:\n\n");
    FOR i := 0 TO NClasses-1 DO
      IF MatchCounts[i]>0 THEN
        Wr.PutText(meswr,Fmt.Pad(Fmt.Int(i),6) &
                              " matches: " &
                              Fmt.Pad(Fmt.Int(MatchCounts[i]),6) &
                              " states\n")
      END 
    END;
    Wr.PutText(meswr,"\n\nNumber of matches for each class:\n\n");
    FOR i := 0 TO NClasses-1 DO
      Wr.PutText(meswr,"      ");
      Wr.PutChar(meswr,ClassPref[i]);
      Wr.PutChar(meswr,' ');
      WrNat(meswr,ClassMatchCounts[i]);
      Wr.PutChar(meswr,'\n')
    END;
    IF ShowFlag THEN
      Wr.PutText(meswr,"\n\nMatches for each state:\n\n");
      FOR s:=UnitState TO Root DO 
        VAR bcls := BestClasses[s];
        BEGIN
          IF bcls#NIL THEN
            Wr.PutText(meswr,Fmt.Pad(Fmt.Int(s),7) & ": ");
            REPEAT
              Wr.PutChar(meswr,ClassPref[bcls^.class]);
              Wr.PutChar(meswr,' ');
              bcls := bcls^.next
            UNTIL bcls=NIL;
            Wr.PutChar(meswr,'\n')
          END
        END
      END
    END;
    ProcessingTime()
  END FindBestMatches;
  
  PROCEDURE BestMatch(s: State): StateClasses RAISES{Wr.Failure} =
  <* FATAL Encoding.BadString *>
    VAR
      dcs1: NAT;
      dcs2: NAT;
      dcs3: NAT := 0;
      sClass: StateClasses := NoMatch;
      aux: StateClasses;
      cs: REF StringRecord;
      len: NAT;
      md3: NAT := MIN(MaxDiff,
                      CEILING((FLOAT(MaxRelDiff)/100.0)*FLOAT(Aut.NSuffs(s))));
  BEGIN
    ProgressReport(BestMatchReportCount);
    FOR c := 0 TO NClasses-1 DO
      WITH 
        md1 = MaxDiff1[c],
        md2 = MaxDiff2[c],
        cr  = Classes[c]
      DO
        dcs1 := nSClasses[c];
        dcs2 := 0;
        cs := ClassStrings[c];
        <* ASSERT cs#NIL *>
        IF Crit1Flag THEN
          REPEAT
            IF Aut.Accepts(s,cs.string^) THEN DEC(dcs1) END;
            cs := cs.next
          UNTIL cs=NIL OR dcs1<=md1
        ELSE
          dcs1 := 0
        END;
        IF dcs1<=md1 AND Crit2Flag THEN
          PROCEDURE EnumAction(READONLY w: String) RAISES {Abort} =
          BEGIN
            IF (NOT ClassAut.Accepts(cr,w)) AND (NOT InRest(w,s))  THEN
              INC(dcs2);
              IF  dcs2>md2 THEN 
                RAISE Abort  
              END
            END
          END EnumAction;
          BEGIN
            TRY
              Aut.EnumSuffs(s,EnumAction)
            EXCEPT
              Abort => (* OK *)
            END;
          END
        END;
        IF (Crit1Flag <= (dcs1<=md1)) AND (Crit2Flag <= (dcs2<=md2)) THEN
          sClass := NEW(REF MatchRecord,class:=c,next:=sClass);
          INC(NMatches)
        END       
      END (* WITH *);
    END (* FOR c *);
    len := CleanSubClasses(sClass);
    IF MaxDiff<LAST(NAT) OR MaxRelDiff<100 THEN
      PROCEDURE EnumAction(READONLY w: String) RAISES {Abort} =
        VAR notfound: BOOL := TRUE;
      BEGIN
        aux := sClass;
        WHILE aux#NIL AND notfound DO
          WITH
            cr = Classes[aux^.class]
          DO
            IF ClassAut.Accepts(cr,w) OR InRest(w,s) THEN
              notfound := FALSE
            END
          END;
          aux := aux^.next
        END;
        IF notfound THEN
          INC(dcs3); 
          IF dcs3>md3 THEN RAISE Abort END
        END
      END EnumAction;

      BEGIN
        TRY
          Aut.EnumSuffs(s,EnumAction)
        EXCEPT
          Abort => sClass := NoMatch; len := 0
        END        
      END
    END;
    aux := sClass;
    WHILE aux#NIL DO
      INC(ClassMatchCounts[aux^.class]);
      aux := aux^.next
    END;
    INC(MatchCounts[len]);
    IF len#0 THEN INC(NMStates) END;
    RETURN sClass
  END BestMatch;
  

  PROCEDURE CleanSubClasses(scl: StateClasses): NAT =
    VAR
      py, x, y: StateClasses;
      len: NAT := 0;
  BEGIN
    x := scl;
    WHILE x#NIL DO
      WITH
        c=x.class
      DO
        INC(len);
        py := x;
        y := x.next;
        WHILE y#NIL DO
          WITH
            d=y.class
          DO
            IF SubClasses[d,c] THEN
              py.next := y.next
            END
          END;
          py := py.next;
          IF py=NIL THEN EXIT END;
          y := py.next
        END;
        x := x.next
      END
    END;
    RETURN len
  END CleanSubClasses;
 
  PROCEDURE GenIdeal() RAISES{Wr.Failure} =
  <* FATAL Encoding.BadString, Abort *>
    VAR InterCount, MinusCount: NAT := 0;
    PROCEDURE EnumIdealAction(READONLY w: String) RAISES {Abort} =
    <* FATAL Wr.Failure *>
      VAR 
        mark: CHAR;
        flag := FALSE;
    BEGIN
      IF Aut.Accepts(Root,w)
         (* Quick and very dirty! *)
            AND (NOT ReprFlag OR StringSuffix(CurrentClassReprString^,w)) THEN
        IF InterFlag THEN
          mark := '*';
          flag := TRUE;
          INC(InterCount)
        ELSIF MinusFlag THEN
          mark := '-';
          flag := TRUE;
          INC(MinusCount)
        END
      END;
      IF flag THEN
        Wr.PutChar(matchwr,mark);
        Wr.PutChar(matchwr,CurrentClassPref);
        Wr.PutChar(matchwr,' ');
        encoding.PrintString(matchwr,w);
        Wr.PutChar(matchwr,'\n');
        ProgressReport(GenIdealCount,10)
      END
    END EnumIdealAction;
  BEGIN
    PhaseStartMessage("Generating ideals ('*' and '-')");
    IF InterFlag OR MinusFlag THEN
      EnumIdeal(Root,EnumIdealAction)
    END;
    Wr.PutText(meswr,Fmt.Int(InterCount) & " * words written\n");
    Wr.PutText(meswr,Fmt.Int(MinusCount) & " - words written\n");
    ProcessingTime()
  END GenIdeal;

  PROCEDURE GenIdealComplement() RAISES{Wr.Failure} =
  <* FATAL Abort, Encoding.BadString *>
    VAR PlusCount: NAT := 0;

    PROCEDURE EnumSuffsAction(READONLY w: String) RAISES {Abort} =
    <* FATAL Wr.Failure *>
      VAR 
        lastMatch: StateClasses := NoMatch;
        len: NAT := 0;
        idtxt: TEXT := "";  
    BEGIN
      IF NOT InIdeal(w,Root,lastMatch,len) AND lastMatch#NoMatch THEN
        REPEAT
          WITH
            c=lastMatch^.class
          DO
            idtxt := idtxt & Text.FromChar(ClassPref[c]);
            lastMatch := lastMatch^.next
          END
        UNTIL lastMatch=NoMatch;
        INC(PlusCount);
        Wr.PutChar(matchwr,'+');
        Wr.PutText(matchwr,idtxt & " ");
        encoding.PrintString(matchwr,SUBARRAY(w,0,len));
        WITH
          n=NUMBER(w)
        DO
          IF len<n THEN
            Wr.PutChar(matchwr,'-');
            encoding.PrintString(matchwr,SUBARRAY(w,len,n-len))
          END
        END;
        Wr.PutChar(matchwr,'\n');
        ProgressReport(GenIdealComplementCount,10);
      END
    END EnumSuffsAction;

  BEGIN
    PhaseStartMessage("Generating ideals complement ('+')");
    IF PlusFlag THEN
      Aut.EnumSuffs(Root,EnumSuffsAction)
    END;
    Wr.PutText(meswr,"\n" & Fmt.Int(PlusCount) & " + words written\n");
    ProcessingTime()
  END GenIdealComplement;

  PROCEDURE GenRadicals() RAISES {Wr.Failure} =
  <* FATAL Abort *>
    VAR 
      radwr: Wr.T;
      radpr,
      auxradpr: StringPrinter.T;
      radCount: NAT := 0;
  BEGIN
    PhaseStartMessage("Generating radicals");
    IF RadFlag THEN
      radwr := FileWr.Open(RadFileName);
      radpr := StringPrinter.New(wr:=radwr,
                                 encoding:=encoding,
                                 empty:="|empty|",
                                 leftMargin:=4);
      auxradpr := StringPrinter.New(wr:=radwr,
                                 encoding:=encoding,
                                 empty:="|empty|",
                                 rightMargin:=LAST(NAT),
                                 sep := "\n");
      Wr.PutText(radwr,"Automaton file: " & DumpFileName & "\n");
      Wr.PutText(radwr,"Automaton comment:\n");
      Wr.PutText(radwr,Aut.doc & "\n");
      
      FOR s:=Root TO UnitState BY -1 DO 
        VAR 
          bcls := BestClasses[s];
        BEGIN
          WHILE bcls#NIL DO
            WITH csbc=ClassStates[bcls^.class] DO
               csbc := NEW(REF ClassStateRecord,state:=s,next:=csbc)
            END;
            bcls := bcls^.next
          END
        END
      END;

      FOR c := 0 TO NClasses-1 DO
        VAR 
          strlist: REF StringRecord := NIL;
          csl: REF ClassStateRecord := ClassStates[c];
          AuxAut: Reduced.T;

        PROCEDURE EnumPrefAction (* Reduced.PrefixAction*)
                                (READONLY w: String) RAISES {Abort} =
        BEGIN
          WITH nw=NUMBER(w),
               nl=nw-1,
               rs=NEW(REF String,nw)
          DO
            FOR i := 0 TO nl DO
              rs^[i] := w[nl-i];
            END;
            strlist := NEW(REF StringRecord,string:=rs,next:=strlist)
          END
        END EnumPrefAction;
        
        PROCEDURE NextString (* Reduced.NextStringProc *)
                            (VAR (*IO*) s: REF String; 
                             VAR (*OUT*) len: NAT;
                             VAR (*OUT*) add: BOOL) RAISES {Done, Abort} =
        BEGIN
          IF strlist=NIL THEN RAISE Done END;
          add := TRUE;
          s := strlist^.string;
          strlist := strlist^.next;
          len := NUMBER(s^)
        END NextString;

        BEGIN
          IF csl#NIL THEN
            radpr.Reset();
            PrintClass(radwr,radpr,c);
            <* ASSERT strlist=NIL *>
            REPEAT
              Aut.EnumPrefs(csl^.state,EnumPrefAction);
              csl := csl^.next
            UNTIL csl=NIL;
            AuxAut := Reduced.New(ClassAutDagSize);
            AuxAut.Build(next:=NextString,
                         wr:=dummywr,
                         reportInterval:=LAST(NAT),
                         flagRedundant:=FALSE);
            WITH ns=AuxAut.NSuffs(AuxAut.Root()) DO
              Wr.PutText(radwr,"\n<<<Radicals: " & Fmt.Int(ns) & "\n");
              auxradpr.Reset();
              AuxAut.PrintSuffs(AuxAut.Root(),auxradpr);
              Wr.PutText(radwr,"\n>>>\n");
              INC(radCount,ns)
            END
          END
        END
      END; (* FOR c *)
      Wr.Flush(radwr)      
    END; (* IF RadFlag *)
    Wr.PutText(meswr,"\n" & Fmt.Int(radCount) & " radicals written\n");
    ProcessingTime()
  END GenRadicals;

  PROCEDURE EnumIdeal(s: State;
                    Action: Reduced.SuffixAction) =
    VAR  
      str: ARRAY[0..MaxLen-1] OF Letter;
      len: NAT := 0;

    PROCEDURE ClassAction(READONLY w: String) RAISES {Abort} =
      VAR 
        CurrLen: NAT:= len;
    BEGIN
      FOR i := 0 TO LAST(w) DO
        str[len] := w[i];
        INC(len)
      END;
      Action(SUBARRAY(str,0,len));
      len := CurrLen
    END ClassAction;

    PROCEDURE EnumRest(s: State) =
      VAR 
        t: State;
    BEGIN
      t := s;
      WHILE Aut.HasArcs(t) DO
        WITH arc=Aut.Last(t) DO
          str[len] := arc.letter;
          INC(len);
          DoEnumIdeal(arc.dest);
          DEC(len);
          t := Aut.Rest(t)
        END
      END;
      RETURN
    END EnumRest;

    PROCEDURE DoEnumIdeal(s: State) =
    <* FATAL Abort *>
      VAR
        bcls: StateClasses;
    BEGIN
      bcls := BestClasses[s];
      WHILE bcls#NIL DO
        WITH
          c=bcls.class
        DO
          CurrentClassPref := ClassPref[c];
          CurrentClassReprString := ClassReprString[c];
          ClassAut.EnumSuffs(Classes[c],ClassAction);
        END;
        bcls := bcls.next
      END;
      EnumRest(s)
    END DoEnumIdeal;

  BEGIN (* EnumIdeal *)
    DoEnumIdeal(s)
  END EnumIdeal;


  PROCEDURE InIdeal(READONLY w: String;
                   s :State;
                   VAR (*INOUT*) lastMatch: StateClasses;
                   VAR (*INOUT*) len: NAT): BOOL =
    VAR 
      localLen: NAT := 0;
  
    PROCEDURE DoInIdeal(READONLY w: String; s :State): BOOL =
      VAR
        bcls: StateClasses := BestClasses[s];
    BEGIN
      IF bcls#NIL THEN lastMatch := bcls; len := localLen END;
      WHILE bcls#NIL DO
        IF ClassAut.Accepts(Classes[bcls.class],w) THEN
          RETURN TRUE
        END;
        bcls := bcls.next
      END;
      INC(localLen);
      RETURN DoInRest(w,s)
    END DoInIdeal;

    PROCEDURE DoInRest(READONLY w: String; s: State): BOOL =
      VAR 
        t: State;
    BEGIN
      WITH n=NUMBER(w) DO
        IF n=0 THEN RETURN FALSE END;
        t := s;
        LOOP
          IF NOT Aut.HasArcs(t) THEN RETURN FALSE END;
          WITH  arc=Aut.Last(t) DO
            IF arc.letter=w[0] THEN 
              RETURN DoInIdeal(SUBARRAY(w,1,n-1),arc.dest)
            END;
            t := Aut.Rest(t)
          END
        END
      END
    END DoInRest;

  BEGIN
    RETURN DoInIdeal(w,s)
  END InIdeal;


  PROCEDURE InRest(READONLY w: String; s: State): BOOL =
    VAR 
      t: State;
      
    PROCEDURE DoInIdeal(READONLY w: String; s :State): BOOL =
      VAR
        bcls: StateClasses := BestClasses[s];
    BEGIN
      WHILE bcls#NIL DO
        IF ClassAut.Accepts(Classes[bcls.class],w) THEN
          RETURN TRUE
        END;
        bcls := bcls.next
      END;
      RETURN InRest(w,s)
    END DoInIdeal;
      
  BEGIN
    WITH n=NUMBER(w) DO
      IF n=0 THEN RETURN FALSE END;
      t := s;
      LOOP
        IF NOT Aut.HasArcs(t) THEN RETURN FALSE END;
        WITH  arc=Aut.Last(t) DO
          IF arc.letter=w[0] THEN 
            RETURN DoInIdeal(SUBARRAY(w,1,n-1),arc.dest)
          END;
          t := Aut.Rest(t)
        END
      END
    END
  END InRest;
  
  PROCEDURE IsState(s: State): BOOL =
  BEGIN
    RETURN Aut.NPrefs(s)>0
  END IsState;

  PROCEDURE PrintClass(wr: Wr.T; 
                       pr: StringPrinter.T;
                       c: [0..MaxNClasses-1]) RAISES {Wr.Failure} =
  <* FATAL Encoding.BadString *>
  BEGIN
    WITH cr = Classes[c] DO
      Wr.PutText(wr,"\n\n#### Class No. ");
      WrNat(wr,c);
      Wr.PutText(wr, "  (");
      Wr.PutChar(wr,ClassPref[c]);                     
      Wr.PutChar(wr,':'); 
      Wr.PutText(wr,ClassComment[c]);
      Wr.PutText(wr,")\n     ");
      WrNat(wr,nSClasses[c]);
      Wr.PutText(wr," suffixes  limits = ");
      WrNat(wr,MaxDiff1[c]);
      Wr.PutText(wr,", ");
      WrNat(wr,MaxDiff2[c]);
      IF ReprFlag THEN
        Wr.PutText(wr,"   representative:  ");
        encoding.PrintString(wr,ClassReprString[c]^)
      END;
      Wr.PutText(wr,"\n\n  { ");
      pr.Reset();
      ClassAut.PrintSuffs(cr,pr);
      Wr.PutText(wr," }\n\n");
      Wr.Flush(wr)
    END;
  END PrintClass;

  PROCEDURE ProcessingTime() RAISES{Wr.Failure,Abort} = 
  BEGIN
    Wr.PutText(meswr,"\nProcessing time for this phase: " & 
                     Util.ExecTimesTexts().TotalTime & "\n");
    Wr.Flush(meswr)                     
  END ProcessingTime;

  PROCEDURE PhaseStartMessage(mess: TEXT) RAISES {Wr.Failure} =
  BEGIN 
    Wr.PutText(meswr,"\n\n" & mess & "\n");
    FOR i:=1 TO Text.Length(mess) DO Wr.PutChar(meswr,'=') END;
    Wr.PutText(meswr,"\n\n");
    Wr.Flush(meswr)
  END PhaseStartMessage;

  PROCEDURE ProgressReport(VAR count: NAT; factor: NAT:=1) =
  <* FATAL Wr.Failure *>
  BEGIN
    INC(count);
    IF (count MOD ReportInterval*factor)=0 THEN
      WrNat(meswr,count);
      Wr.PutChar(meswr,'\n');
      Wr.Flush(meswr)
    END
  END ProgressReport;

  PROCEDURE StringSuffix(READONLY v,w: String): BOOL =
  BEGIN
    WITH
      lv=LAST(v),
      lw=LAST(w)
    DO
      IF lw>=lv THEN
        FOR i:=0 TO lv DO
          IF v[lv-i]#w[lw-i] THEN RETURN FALSE END
        END;
        RETURN TRUE
      ELSE
        RETURN FALSE 
      END
    END
  END StringSuffix;

BEGIN
  Main()
END IdealClasses.
