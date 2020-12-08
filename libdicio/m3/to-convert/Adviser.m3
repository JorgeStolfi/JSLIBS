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

MODULE Adviser;

 IMPORT Char, Text, ReadOnlyPermDAG, AVLTextTree;
 
 FROM Char IMPORT NUL;
 FROM ReadOnlyPermDAG IMPORT State;
 FROM Basics IMPORT Abort;
 
EXCEPTION NoTransition;

 CONST
    (* Minimum length to try some searches *)
    minTransp  = 4;
    minRemove = 4;
    minSubst  = 4;
    minInsert = 5;
  
    (* Sounds: using chars in [@-Z] -- assumes they do not appear in words! *)
    
    FirstSound = '@';
    LastSound = 'Z';
    
    SoundA                 = '@'; (* all accented forms of 'a' 
                                     will appear as 'a' itself ! *)
    
    SoundCCcX              = 'A';
    SoundCCed              = 'B';
    SoundCedCced           = 'C';
    SoundCedS              = 'D';
    SoundCedScedSsCcedPced = 'E';
    SoundChSXZ             = 'F';
    SoundChX               = 'G';
    SoundCScSsXcXCc        = 'H';
    SoundEI                = 'I'; (* includes accented forms of 'e' and 'i' *)
    SoundF                 = 'J';
    SoundGJ                = 'K';
    SoundLOU               = 'L';
    SoundLU                = 'M';
    SoundMN                = 'N';
    SoundOU                = 'O'; (* includes accented forms of 'o' and 'u' *)
    SoundS                 = 'P';
    SoundSXZ               = 'Q';
    SoundScXc              = 'R';
    SoundSCSc              = 'S'; 
    SoundX                 = 'T';
              
    (* optional sounds *)
    
    XSoundC                = 'U';
    XSoundEI               = 'V';
    XSoundH                = 'W';    
    XSoundP                = 'X';
    XSoundPC               = 'Y';
    XSoundU                = 'Z';
                      
    Sounds = SET OF CHAR{'@' .. 'Z'};                              
                              
    Letters = SET OF CHAR{'a' .. 'z','ç'};
    InsertionLetters = ARRAY [0 .. 26] OF CHAR{'a','b','c','d','e','f','g','h','i',
               'j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','ç'};
                                  
    MaxAutoState = 27;
    MaxAutoAction = 56;
    MaxSubstLetters = 3;
    MaxSubst = 8;
    
    NoSub     = Subst{NUL,..};
    
  TYPE
    
    ChArr = ARRAY OF CHAR;
    
    AutoState = [0 .. MaxAutoState];
    AutoAction = [0 .. MaxAutoAction];

    AutoActionRec = RECORD
        NextState: AutoState;
        Action: AutoAction
      END;
    
    Subst = ARRAY [0 .. MaxSubstLetters-1] OF CHAR;
    ArrSubst = ARRAY [0..MaxSubst-1] OF  Subst;
     
  VAR
  
    SoundAutomaton:  ARRAY AutoState,CHAR OF AutoActionRec;
    StripAccents := ARRAY CHAR OF CHAR{Char.NUL,..};
    SoundSubst: ARRAY [FirstSound .. LastSound] OF  ArrSubst;
    
      (* Allocate only once ! *)
      
    InWord := NEW(REF ChArr,100);
    InWord2 := NEW(REF ChArr,100);
    OutWord := NEW(REF ChArr,100);
    
    
  PROCEDURE InitAlternatives() =
    CONST
      SubNUL    = Subst{NUL,..};
      SubC      = Subst{'c',NUL,..};
      SubCc     = Subst{'c','c',NUL,..};
      SubCced   = Subst{'c','ç',NUL,..};
      SubCed    = Subst{'ç',NUL,..};
      SubCh     = Subst{'c','h',NUL,..};
      SubF      = Subst{'f',NUL,..};
      SubG      = Subst{'g',NUL,..};
      SubH      = Subst{'h',NUL,..};
      SubJ      = Subst{'j',NUL,..};
      SubL      = Subst{'l',NUL,..};
      SubM      = Subst{'m',NUL,..};
      SubN      = Subst{'n',NUL,..};
      SubO      = Subst{'o',NUL,..};
      SubP      = Subst{'p',NUL,..};
      SubU      = Subst{'u',NUL,..};
      SubPced   = Subst{'p','ç',NUL,..};
      SubS      = Subst{'s',NUL,..};
      SubSc     = Subst{'s','c',NUL,..};
      SubSced   = Subst{'s','ç',NUL,..};
      SubSs     = Subst{'s','s',NUL,..};
      SubX      = Subst{'x',NUL,..};
      SubXc     = Subst{'x','c',NUL,..};
      SubZ      = Subst{'z',NUL,..};
      
  BEGIN
    
    SoundSubst[SoundA]            := ArrSubst{Subst{'a',NUL,..},Subst{'à',NUL,..},Subst{'á',NUL,..},
                                              Subst{'â',NUL,..},Subst{'ã',NUL,..},NoSub,..};
    SoundSubst[SoundCCcX]         := ArrSubst{SubC,SubCc,SubX,NoSub,..};
    SoundSubst[SoundCCed]         := ArrSubst{SubC,SubCed,NoSub,..};
    SoundSubst[SoundCedCced]      := ArrSubst{SubCed,SubCced,NoSub,..};
    SoundSubst[SoundCedS]         := ArrSubst{SubCed,SubS,NoSub,..};
    SoundSubst[SoundCedScedSsCcedPced] := ArrSubst{SubCed,SubSced,SubSs,SubCced,SubPced,NoSub,..};
    SoundSubst[SoundChSXZ]        := ArrSubst{SubCh,SubX,SubZ,NoSub,..};
    SoundSubst[SoundChX]          := ArrSubst{SubCh,SubX,NoSub,..};
    SoundSubst[SoundCScSsXcXCc]   := ArrSubst{SubC,SubSc,SubSs,SubXc,SubX,SubCc,NoSub,..};
    SoundSubst[SoundEI]           := ArrSubst{Subst{'e',NUL,..},Subst{'é',NUL,..},Subst{'ê',NUL,..},
                                               Subst{'i',NUL,..},Subst{'í',NUL,..},NoSub,..};
    SoundSubst[SoundF]            := ArrSubst{SubF,NoSub,..};      
    SoundSubst[SoundGJ]           := ArrSubst{SubG,SubJ,NoSub,..};
    SoundSubst[SoundLOU]          := ArrSubst{SubL,SubO,SubU,NoSub,..};
    SoundSubst[SoundLU]           := ArrSubst{SubL,SubU,NoSub,..};
    SoundSubst[SoundMN]           := ArrSubst{SubM,SubN,NoSub,..};
    SoundSubst[SoundOU]           := ArrSubst{Subst{'o',NUL,..},Subst{'ó',NUL,..},Subst{'ô',NUL,..},
                                               Subst{'õ',NUL,..},Subst{'u',NUL,..},Subst{'ú',NUL,..},
                                               Subst{'ü',NUL,..},NoSub,..};
    SoundSubst[SoundS]            := ArrSubst{SubS,NoSub,..};
    SoundSubst[SoundSXZ]          := ArrSubst{SubS,SubX,SubZ,NoSub,..};
    SoundSubst[SoundScXc]         := ArrSubst{SubSc,SubXc,NoSub,..};
    SoundSubst[SoundSCSc]         := ArrSubst{SubS,SubC,SubSc,NoSub,..};
    SoundSubst[SoundX]            := ArrSubst{SubX,NoSub,..};
                           
    SoundSubst[XSoundC]           := ArrSubst{SubNUL,SubC,NoSub,..};
    SoundSubst[XSoundEI]          := ArrSubst{SubNUL,Subst{'e',NUL,..},Subst{'é',NUL,..},
                                               Subst{'ê',NUL,..},Subst{'i',NUL,..},
                                               Subst{'í',NUL,..},NoSub,..};
    SoundSubst[XSoundH]           := ArrSubst{SubNUL,SubH,NoSub,..};  
    SoundSubst[XSoundP]           := ArrSubst{SubNUL,SubP,NoSub,..};
    SoundSubst[XSoundPC]          := ArrSubst{SubNUL,SubP,SubC,NoSub,..};
    SoundSubst[XSoundU]           := ArrSubst{SubNUL,SubU,NoSub,..};

    
    
  END InitAlternatives;

  PROCEDURE InitSoundex() =
  BEGIN
  
    FOR c := 'a' TO  'z'  DO StripAccents[c] := c END;
    StripAccents['à'] := 'a';
    StripAccents['á'] := 'a';
    StripAccents['â'] := 'a';
    StripAccents['ã'] := 'a';
    StripAccents['é'] := 'e';
    StripAccents['ê'] := 'e';
    StripAccents['í'] := 'i';
    StripAccents['ó'] := 'o';
    StripAccents['ô'] := 'o';
    StripAccents['õ'] := 'o';
    StripAccents['ú'] := 'u';
    StripAccents['ü'] := 'u';
    StripAccents['ç'] := 'ç';
  
    FOR s:=0 TO MaxAutoState DO
      FOR c:=FIRST(CHAR) TO LAST(CHAR) DO 
          SoundAutomaton[s,c] := AutoActionRec{NextState:=0,Action:=0} END
    END;
    
    (* '#' is the word end marker *)
    
    (* State 0 *)
    (* inicio,fim,h^,- *)

    SoundAutomaton[ 0,'#'] := AutoActionRec{NextState:= 0,Action:= 2};
    SoundAutomaton[ 0,'a'] := AutoActionRec{NextState:= 1,Action:=16};
    SoundAutomaton[ 0,'b'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 0,'c'] := AutoActionRec{NextState:= 3,Action:= 1};
    SoundAutomaton[ 0,'d'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 0,'e'] := AutoActionRec{NextState:= 1,Action:=16};
    SoundAutomaton[ 0,'f'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 0,'g'] := AutoActionRec{NextState:= 4,Action:= 1};
    SoundAutomaton[ 0,'h'] := AutoActionRec{NextState:= 0,Action:=51};
    SoundAutomaton[ 0,'i'] := AutoActionRec{NextState:= 1,Action:=16};
    SoundAutomaton[ 0,'j'] := AutoActionRec{NextState:= 5,Action:= 1};
    SoundAutomaton[ 0,'k'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 0,'l'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 0,'m'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 0,'n'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 0,'o'] := AutoActionRec{NextState:= 1,Action:=16};
    SoundAutomaton[ 0,'p'] := AutoActionRec{NextState:=21,Action:= 1};
    SoundAutomaton[ 0,'q'] := AutoActionRec{NextState:=14,Action:= 2};
    SoundAutomaton[ 0,'r'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 0,'s'] := AutoActionRec{NextState:= 6,Action:= 1};
    SoundAutomaton[ 0,'t'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 0,'u'] := AutoActionRec{NextState:= 1,Action:=16};
    SoundAutomaton[ 0,'v'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 0,'w'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 0,'x'] := AutoActionRec{NextState:= 2,Action:=17};
    SoundAutomaton[ 0,'y'] := AutoActionRec{NextState:= 1,Action:= 2};
    SoundAutomaton[ 0,'z'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 0,'ç'] := AutoActionRec{NextState:= 6,Action:= 1};


    (* State 1 *)
    (* ...V^ *)

    SoundAutomaton[ 1,'#'] := AutoActionRec{NextState:= 0,Action:= 2};
    SoundAutomaton[ 1,'a'] := AutoActionRec{NextState:= 1,Action:= 2};
    SoundAutomaton[ 1,'b'] := AutoActionRec{NextState:= 7,Action:= 2};
    SoundAutomaton[ 1,'c'] := AutoActionRec{NextState:= 8,Action:= 1};
    SoundAutomaton[ 1,'d'] := AutoActionRec{NextState:= 7,Action:= 2};
    SoundAutomaton[ 1,'e'] := AutoActionRec{NextState:= 1,Action:= 7};
    SoundAutomaton[ 1,'f'] := AutoActionRec{NextState:= 7,Action:= 2};
    SoundAutomaton[ 1,'g'] := AutoActionRec{NextState:= 4,Action:= 1};
    SoundAutomaton[ 1,'h'] := AutoActionRec{NextState:= 2,Action:=51};
    SoundAutomaton[ 1,'i'] := AutoActionRec{NextState:= 1,Action:= 7};
    SoundAutomaton[ 1,'j'] := AutoActionRec{NextState:= 5,Action:= 1};
    SoundAutomaton[ 1,'k'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 1,'l'] := AutoActionRec{NextState:= 9,Action:= 1};
    SoundAutomaton[ 1,'m'] := AutoActionRec{NextState:=10,Action:= 1};
    SoundAutomaton[ 1,'n'] := AutoActionRec{NextState:=10,Action:= 1};
    SoundAutomaton[ 1,'o'] := AutoActionRec{NextState:= 1,Action:= 7};
    SoundAutomaton[ 1,'p'] := AutoActionRec{NextState:=21,Action:= 1};
    SoundAutomaton[ 1,'q'] := AutoActionRec{NextState:=14,Action:= 2};
    SoundAutomaton[ 1,'r'] := AutoActionRec{NextState:=11,Action:= 2};
    SoundAutomaton[ 1,'s'] := AutoActionRec{NextState:=12,Action:= 1};
    SoundAutomaton[ 1,'t'] := AutoActionRec{NextState:= 7,Action:=18};
    SoundAutomaton[ 1,'u'] := AutoActionRec{NextState:=18,Action:= 1};
    SoundAutomaton[ 1,'v'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 1,'w'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 1,'x'] := AutoActionRec{NextState:=13,Action:= 1};
    SoundAutomaton[ 1,'y'] := AutoActionRec{NextState:= 1,Action:= 2};
    SoundAutomaton[ 1,'z'] := AutoActionRec{NextState:=12,Action:= 1};
    SoundAutomaton[ 1,'ç'] := AutoActionRec{NextState:=15,Action:= 1};


    (* State 2 *)
    (* ...C^ *)

    SoundAutomaton[ 2,'#'] := AutoActionRec{NextState:= 0,Action:= 2};
    SoundAutomaton[ 2,'a'] := AutoActionRec{NextState:= 1,Action:= 2};
    SoundAutomaton[ 2,'b'] := AutoActionRec{NextState:= 2,Action:=11};
    SoundAutomaton[ 2,'c'] := AutoActionRec{NextState:= 3,Action:= 1};
    SoundAutomaton[ 2,'d'] := AutoActionRec{NextState:= 2,Action:=11};
    SoundAutomaton[ 2,'e'] := AutoActionRec{NextState:= 1,Action:= 7};
    SoundAutomaton[ 2,'f'] := AutoActionRec{NextState:= 2,Action:=11};
    SoundAutomaton[ 2,'g'] := AutoActionRec{NextState:= 4,Action:= 1};
    SoundAutomaton[ 2,'h'] := AutoActionRec{NextState:= 2,Action:= 1};
    SoundAutomaton[ 2,'i'] := AutoActionRec{NextState:= 1,Action:= 7};
    SoundAutomaton[ 2,'j'] := AutoActionRec{NextState:= 5,Action:= 1};
    SoundAutomaton[ 2,'k'] := AutoActionRec{NextState:= 2,Action:=11};
    SoundAutomaton[ 2,'l'] := AutoActionRec{NextState:= 2,Action:=11};
    SoundAutomaton[ 2,'m'] := AutoActionRec{NextState:= 2,Action:=11};
    SoundAutomaton[ 2,'n'] := AutoActionRec{NextState:= 2,Action:=11};
    SoundAutomaton[ 2,'o'] := AutoActionRec{NextState:= 1,Action:= 7};
    SoundAutomaton[ 2,'p'] := AutoActionRec{NextState:=27,Action:=11};
    SoundAutomaton[ 2,'q'] := AutoActionRec{NextState:=14,Action:= 2};
    SoundAutomaton[ 2,'r'] := AutoActionRec{NextState:= 2,Action:=11};
    SoundAutomaton[ 2,'s'] := AutoActionRec{NextState:=17,Action:= 1};
    SoundAutomaton[ 2,'t'] := AutoActionRec{NextState:= 2,Action:=11};
    SoundAutomaton[ 2,'u'] := AutoActionRec{NextState:= 1,Action:= 7};
    SoundAutomaton[ 2,'v'] := AutoActionRec{NextState:= 2,Action:=11};
    SoundAutomaton[ 2,'w'] := AutoActionRec{NextState:= 2,Action:=11};
    SoundAutomaton[ 2,'x'] := AutoActionRec{NextState:= 2,Action:=11};
    SoundAutomaton[ 2,'y'] := AutoActionRec{NextState:= 1,Action:= 2};
    SoundAutomaton[ 2,'z'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[ 2,'ç'] := AutoActionRec{NextState:=17,Action:= 1};


    (* State 3 *)
    (* ...C^c, ^c *)

    SoundAutomaton[ 3,'#'] := AutoActionRec{NextState:= 0,Action:= 4};
    SoundAutomaton[ 3,'a'] := AutoActionRec{NextState:= 1,Action:=53};
    SoundAutomaton[ 3,'b'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[ 3,'c'] := AutoActionRec{NextState:= 3,Action:= 1};
    SoundAutomaton[ 3,'d'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[ 3,'e'] := AutoActionRec{NextState:= 1,Action:=19};
    SoundAutomaton[ 3,'f'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[ 3,'g'] := AutoActionRec{NextState:= 4,Action:= 3};
    SoundAutomaton[ 3,'h'] := AutoActionRec{NextState:= 2,Action:=17};
    SoundAutomaton[ 3,'i'] := AutoActionRec{NextState:= 1,Action:=19};
    SoundAutomaton[ 3,'j'] := AutoActionRec{NextState:= 5,Action:= 3};
    SoundAutomaton[ 3,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[ 3,'l'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[ 3,'m'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[ 3,'n'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[ 3,'o'] := AutoActionRec{NextState:= 1,Action:=53};
    SoundAutomaton[ 3,'p'] := AutoActionRec{NextState:=27,Action:= 3};
    SoundAutomaton[ 3,'q'] := AutoActionRec{NextState:=14,Action:= 2};
    SoundAutomaton[ 3,'r'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[ 3,'s'] := AutoActionRec{NextState:=17,Action:= 3};
    SoundAutomaton[ 3,'t'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[ 3,'u'] := AutoActionRec{NextState:= 1,Action:=53};
    SoundAutomaton[ 3,'v'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[ 3,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[ 3,'x'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[ 3,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[ 3,'z'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[ 3,'ç'] := AutoActionRec{NextState:=17,Action:= 3};


    (* State 4 *)
    (* ...^g *)

    SoundAutomaton[4,'#'] := AutoActionRec{NextState:= 0,Action:= 4};
    SoundAutomaton[4,'a'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[4,'b'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[4,'c'] := AutoActionRec{NextState:= 3,Action:= 3};
    SoundAutomaton[4,'d'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[4,'e'] := AutoActionRec{NextState:= 1,Action:=20};
    SoundAutomaton[4,'f'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[4,'g'] := AutoActionRec{NextState:= 4,Action:= 1};
    SoundAutomaton[4,'h'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[4,'i'] := AutoActionRec{NextState:= 1,Action:=20};
    SoundAutomaton[4,'j'] := AutoActionRec{NextState:= 5,Action:= 3};
    SoundAutomaton[4,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[4,'l'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[4,'m'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[4,'n'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[4,'o'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[4,'p'] := AutoActionRec{NextState:=27,Action:= 3};
    SoundAutomaton[4,'q'] := AutoActionRec{NextState:=14,Action:= 4};
    SoundAutomaton[4,'r'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[4,'s'] := AutoActionRec{NextState:=17,Action:= 3};
    SoundAutomaton[4,'t'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[4,'u'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[4,'v'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[4,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[4,'x'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[4,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[4,'z'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[4,'ç'] := AutoActionRec{NextState:=17,Action:= 3};


    (* State 5 *)
    (* ...^j *)

    SoundAutomaton[5,'#'] := AutoActionRec{NextState:= 0,Action:= 4};
    SoundAutomaton[5,'a'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[5,'b'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[5,'c'] := AutoActionRec{NextState:= 3,Action:= 3};
    SoundAutomaton[5,'d'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[5,'e'] := AutoActionRec{NextState:= 1,Action:=20};
    SoundAutomaton[5,'f'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[5,'g'] := AutoActionRec{NextState:= 4,Action:= 3};
    SoundAutomaton[5,'h'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[5,'i'] := AutoActionRec{NextState:= 1,Action:=20};
    SoundAutomaton[5,'j'] := AutoActionRec{NextState:= 5,Action:= 1};
    SoundAutomaton[5,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[5,'l'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[5,'m'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[5,'n'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[5,'o'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[5,'p'] := AutoActionRec{NextState:=27,Action:= 3};
    SoundAutomaton[5,'q'] := AutoActionRec{NextState:=14,Action:= 4};
    SoundAutomaton[5,'r'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[5,'s'] := AutoActionRec{NextState:=17,Action:= 3};
    SoundAutomaton[5,'t'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[5,'u'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[5,'v'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[5,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[5,'x'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[5,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[5,'z'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[5,'ç'] := AutoActionRec{NextState:=17,Action:= 3};


    (* State 6 *)
    (* ^s,^ç *)

    SoundAutomaton[6,'#'] := AutoActionRec{NextState:= 0,Action:= 4};
    SoundAutomaton[6,'a'] := AutoActionRec{NextState:= 1,Action:=37};
    SoundAutomaton[6,'b'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[6,'c'] := AutoActionRec{NextState:= 3,Action:= 1};
    SoundAutomaton[6,'d'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[6,'e'] := AutoActionRec{NextState:= 1,Action:=19};
    SoundAutomaton[6,'f'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[6,'g'] := AutoActionRec{NextState:= 4,Action:= 3};
    SoundAutomaton[6,'h'] := AutoActionRec{NextState:= 2,Action:=17};
    SoundAutomaton[6,'i'] := AutoActionRec{NextState:= 1,Action:=19};
    SoundAutomaton[6,'j'] := AutoActionRec{NextState:= 5,Action:= 3};
    SoundAutomaton[6,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[6,'l'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[6,'m'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[6,'n'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[6,'o'] := AutoActionRec{NextState:= 1,Action:=37};
    SoundAutomaton[6,'p'] := AutoActionRec{NextState:=27,Action:= 3};
    SoundAutomaton[6,'q'] := AutoActionRec{NextState:=14,Action:= 4};
    SoundAutomaton[6,'r'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[6,'s'] := AutoActionRec{NextState:= 6,Action:= 1};
    SoundAutomaton[6,'t'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[6,'u'] := AutoActionRec{NextState:= 1,Action:=37};
    SoundAutomaton[6,'v'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[6,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[6,'x'] := AutoActionRec{NextState:= 2,Action:=17};
    SoundAutomaton[6,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[6,'z'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[6,'ç'] := AutoActionRec{NextState:= 6,Action:= 1};


    (* State 7 *)
    (* ...V[bdft]^ *)

    SoundAutomaton[7,'#'] := AutoActionRec{NextState:= 0,Action:= 2};
    SoundAutomaton[7,'a'] := AutoActionRec{NextState:= 1,Action:= 2};
    SoundAutomaton[7,'b'] := AutoActionRec{NextState:= 7,Action:=12};
    SoundAutomaton[7,'c'] := AutoActionRec{NextState:= 3,Action:=21};
    SoundAutomaton[7,'d'] := AutoActionRec{NextState:= 7,Action:=12};
    SoundAutomaton[7,'e'] := AutoActionRec{NextState:=19,Action:= 1};
    SoundAutomaton[7,'f'] := AutoActionRec{NextState:= 7,Action:=12};
    SoundAutomaton[7,'g'] := AutoActionRec{NextState:= 4,Action:=21};
    SoundAutomaton[7,'h'] := AutoActionRec{NextState:= 2,Action:= 1};
    SoundAutomaton[7,'i'] := AutoActionRec{NextState:=19,Action:= 1};
    SoundAutomaton[7,'j'] := AutoActionRec{NextState:= 5,Action:=21};
    SoundAutomaton[7,'k'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[7,'l'] := AutoActionRec{NextState:= 2,Action:=22};
    SoundAutomaton[7,'m'] := AutoActionRec{NextState:= 2,Action:=22};
    SoundAutomaton[7,'n'] := AutoActionRec{NextState:= 2,Action:=22};
    SoundAutomaton[7,'o'] := AutoActionRec{NextState:= 1,Action:= 7};
    SoundAutomaton[7,'p'] := AutoActionRec{NextState:=27,Action:= 1};
    SoundAutomaton[7,'q'] := AutoActionRec{NextState:=14,Action:=22};
    SoundAutomaton[7,'r'] := AutoActionRec{NextState:= 2,Action:=22};
    SoundAutomaton[7,'s'] := AutoActionRec{NextState:=17,Action:= 1};
    SoundAutomaton[7,'t'] := AutoActionRec{NextState:= 7,Action:=12};
    SoundAutomaton[7,'u'] := AutoActionRec{NextState:= 1,Action:= 7};
    SoundAutomaton[7,'v'] := AutoActionRec{NextState:= 2,Action:=22};
    SoundAutomaton[7,'w'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[7,'x'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[7,'y'] := AutoActionRec{NextState:= 1,Action:= 2};
    SoundAutomaton[7,'z'] := AutoActionRec{NextState:= 2,Action:=22};
    SoundAutomaton[7,'ç'] := AutoActionRec{NextState:=17,Action:= 1};


    (* State 8 *)
    (* ...V^c *)

    SoundAutomaton[8,'#'] := AutoActionRec{NextState:= 0,Action:= 4};
    SoundAutomaton[8,'a'] := AutoActionRec{NextState:= 1,Action:=53};
    SoundAutomaton[8,'b'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[8,'c'] := AutoActionRec{NextState:=16,Action:= 1};
    SoundAutomaton[8,'d'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[8,'e'] := AutoActionRec{NextState:= 1,Action:=23};
    SoundAutomaton[8,'f'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[8,'g'] := AutoActionRec{NextState:= 4,Action:= 3};
    SoundAutomaton[8,'h'] := AutoActionRec{NextState:= 2,Action:=17};
    SoundAutomaton[8,'i'] := AutoActionRec{NextState:= 1,Action:=23};
    SoundAutomaton[8,'j'] := AutoActionRec{NextState:= 5,Action:= 3};
    SoundAutomaton[8,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[8,'l'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[8,'m'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[8,'n'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[8,'o'] := AutoActionRec{NextState:= 1,Action:=53};
    SoundAutomaton[8,'p'] := AutoActionRec{NextState:=27,Action:= 3};
    SoundAutomaton[8,'q'] := AutoActionRec{NextState:=14,Action:= 2};
    SoundAutomaton[8,'r'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[8,'s'] := AutoActionRec{NextState:=20,Action:= 1};
    SoundAutomaton[8,'t'] := AutoActionRec{NextState:= 2,Action:=24};
    SoundAutomaton[8,'u'] := AutoActionRec{NextState:= 1,Action:=53};
    SoundAutomaton[8,'v'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[8,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[8,'x'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[8,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[8,'z'] := AutoActionRec{NextState:=20,Action:= 1};
    SoundAutomaton[8,'ç'] := AutoActionRec{NextState:=26,Action:= 1};


    (* State 9 *)
    (* ...V^l *)

    SoundAutomaton[9,'#'] := AutoActionRec{NextState:= 0,Action:=28};
    SoundAutomaton[9,'a'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[9,'b'] := AutoActionRec{NextState:= 2,Action:=28};
    SoundAutomaton[9,'c'] := AutoActionRec{NextState:= 3,Action:=39};
    SoundAutomaton[9,'d'] := AutoActionRec{NextState:= 2,Action:=28};
    SoundAutomaton[9,'e'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[9,'f'] := AutoActionRec{NextState:= 2,Action:=28};
    SoundAutomaton[9,'g'] := AutoActionRec{NextState:= 4,Action:=39};
    SoundAutomaton[9,'h'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[9,'i'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[9,'j'] := AutoActionRec{NextState:= 5,Action:=39};
    SoundAutomaton[9,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[9,'l'] := AutoActionRec{NextState:= 9,Action:= 1};
    SoundAutomaton[9,'m'] := AutoActionRec{NextState:= 2,Action:=28};
    SoundAutomaton[9,'n'] := AutoActionRec{NextState:= 2,Action:=28};
    SoundAutomaton[9,'o'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[9,'p'] := AutoActionRec{NextState:=27,Action:=39};
    SoundAutomaton[9,'q'] := AutoActionRec{NextState:=14,Action:=28};
    SoundAutomaton[9,'r'] := AutoActionRec{NextState:= 2,Action:=28};
    SoundAutomaton[9,'s'] := AutoActionRec{NextState:=17,Action:= 3};
    SoundAutomaton[9,'t'] := AutoActionRec{NextState:= 2,Action:=28};
    SoundAutomaton[9,'u'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[9,'v'] := AutoActionRec{NextState:= 2,Action:=28};
    SoundAutomaton[9,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[9,'x'] := AutoActionRec{NextState:= 2,Action:=25};
    SoundAutomaton[9,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[9,'z'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[9,'ç'] := AutoActionRec{NextState:=17,Action:= 3};


    (* State 10 *)
    (* ...V^[mn] *)

    SoundAutomaton[10,'#'] := AutoActionRec{NextState:= 0,Action:=27};
    SoundAutomaton[10,'a'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[10,'b'] := AutoActionRec{NextState:= 2,Action:=27};
    SoundAutomaton[10,'c'] := AutoActionRec{NextState:= 3,Action:=26};
    SoundAutomaton[10,'d'] := AutoActionRec{NextState:= 2,Action:=27};
    SoundAutomaton[10,'e'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[10,'f'] := AutoActionRec{NextState:= 2,Action:=27};
    SoundAutomaton[10,'g'] := AutoActionRec{NextState:= 4,Action:=26};
    SoundAutomaton[10,'h'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[10,'i'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[10,'j'] := AutoActionRec{NextState:= 5,Action:=26};
    SoundAutomaton[10,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[10,'l'] := AutoActionRec{NextState:= 2,Action:=27};
    SoundAutomaton[10,'m'] := AutoActionRec{NextState:=10,Action:=10};
    SoundAutomaton[10,'n'] := AutoActionRec{NextState:=10,Action:=10};
    SoundAutomaton[10,'o'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[10,'p'] := AutoActionRec{NextState:=27,Action:=26};
    SoundAutomaton[10,'q'] := AutoActionRec{NextState:=14,Action:=27};
    SoundAutomaton[10,'r'] := AutoActionRec{NextState:= 2,Action:=27};
    SoundAutomaton[10,'s'] := AutoActionRec{NextState:=17,Action:=26};
    SoundAutomaton[10,'t'] := AutoActionRec{NextState:= 2,Action:=27};
    SoundAutomaton[10,'u'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[10,'v'] := AutoActionRec{NextState:= 2,Action:=27};
    SoundAutomaton[10,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[10,'x'] := AutoActionRec{NextState:=22,Action:=26};
    SoundAutomaton[10,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[10,'z'] := AutoActionRec{NextState:= 2,Action:=27};
    SoundAutomaton[10,'ç'] := AutoActionRec{NextState:=17,Action:=26};


    (* State 11 *)
    (* ...Vr^ *)

    SoundAutomaton[11,'#'] := AutoActionRec{NextState:= 0,Action:= 2};
    SoundAutomaton[11,'a'] := AutoActionRec{NextState:= 1,Action:= 2};
    SoundAutomaton[11,'b'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[11,'c'] := AutoActionRec{NextState:= 3,Action:= 1};
    SoundAutomaton[11,'d'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[11,'e'] := AutoActionRec{NextState:= 1,Action:= 7};
    SoundAutomaton[11,'f'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[11,'g'] := AutoActionRec{NextState:= 4,Action:= 1};
    SoundAutomaton[11,'h'] := AutoActionRec{NextState:= 2,Action:= 1};
    SoundAutomaton[11,'i'] := AutoActionRec{NextState:= 1,Action:= 7};
    SoundAutomaton[11,'j'] := AutoActionRec{NextState:= 5,Action:= 1};
    SoundAutomaton[11,'k'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[11,'l'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[11,'m'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[11,'n'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[11,'o'] := AutoActionRec{NextState:= 1,Action:= 7};
    SoundAutomaton[11,'p'] := AutoActionRec{NextState:=27,Action:= 1};
    SoundAutomaton[11,'q'] := AutoActionRec{NextState:=14,Action:= 2};
    SoundAutomaton[11,'r'] := AutoActionRec{NextState:=23,Action:= 1};
    SoundAutomaton[11,'s'] := AutoActionRec{NextState:=17,Action:= 1};
    SoundAutomaton[11,'t'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[11,'u'] := AutoActionRec{NextState:= 1,Action:= 7};
    SoundAutomaton[11,'v'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[11,'w'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[11,'x'] := AutoActionRec{NextState:= 2,Action:=17};
    SoundAutomaton[11,'y'] := AutoActionRec{NextState:= 1,Action:= 2};
    SoundAutomaton[11,'z'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[11,'ç'] := AutoActionRec{NextState:=17,Action:= 1};

    (* State 12 *)
    (* ...V^[sz] *)

    SoundAutomaton[12,'#'] := AutoActionRec{NextState:= 0,Action:=29};
    SoundAutomaton[12,'a'] := AutoActionRec{NextState:= 1,Action:=29};
    SoundAutomaton[12,'b'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[12,'c'] := AutoActionRec{NextState:=24,Action:= 1};
    SoundAutomaton[12,'d'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[12,'e'] := AutoActionRec{NextState:= 1,Action:=30};
    SoundAutomaton[12,'f'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[12,'g'] := AutoActionRec{NextState:= 4,Action:=31};
    SoundAutomaton[12,'h'] := AutoActionRec{NextState:= 2,Action:=17};
    SoundAutomaton[12,'i'] := AutoActionRec{NextState:= 1,Action:=30};
    SoundAutomaton[12,'j'] := AutoActionRec{NextState:= 5,Action:=31};
    SoundAutomaton[12,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[12,'l'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[12,'m'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[12,'n'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[12,'o'] := AutoActionRec{NextState:= 1,Action:=30};
    SoundAutomaton[12,'p'] := AutoActionRec{NextState:=27,Action:=31};
    SoundAutomaton[12,'q'] := AutoActionRec{NextState:=14,Action:=29};
    SoundAutomaton[12,'r'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[12,'s'] := AutoActionRec{NextState:=15,Action:= 1};
    SoundAutomaton[12,'t'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[12,'u'] := AutoActionRec{NextState:= 1,Action:=30};
    SoundAutomaton[12,'v'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[12,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[12,'x'] := AutoActionRec{NextState:=13,Action:= 1};
    SoundAutomaton[12,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[12,'z'] := AutoActionRec{NextState:=12,Action:= 1};
    SoundAutomaton[12,'ç'] := AutoActionRec{NextState:=15,Action:= 1};

    (* State 13 *)
    (* ...V^x *)

    SoundAutomaton[13,'#'] := AutoActionRec{NextState:= 0,Action:=29};
    SoundAutomaton[13,'a'] := AutoActionRec{NextState:= 1,Action:=32};
    SoundAutomaton[13,'b'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[13,'c'] := AutoActionRec{NextState:=24,Action:= 1};
    SoundAutomaton[13,'d'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[13,'e'] := AutoActionRec{NextState:= 1,Action:=33};
    SoundAutomaton[13,'f'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[13,'g'] := AutoActionRec{NextState:= 4,Action:=31};
    SoundAutomaton[13,'h'] := AutoActionRec{NextState:= 2,Action:= 3};
    SoundAutomaton[13,'i'] := AutoActionRec{NextState:= 1,Action:=33};
    SoundAutomaton[13,'j'] := AutoActionRec{NextState:= 5,Action:=31};
    SoundAutomaton[13,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[13,'l'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[13,'m'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[13,'n'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[13,'o'] := AutoActionRec{NextState:= 1,Action:=33};
    SoundAutomaton[13,'p'] := AutoActionRec{NextState:=27,Action:=31};
    SoundAutomaton[13,'q'] := AutoActionRec{NextState:=14,Action:=29};
    SoundAutomaton[13,'r'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[13,'s'] := AutoActionRec{NextState:=24,Action:= 1};
    SoundAutomaton[13,'t'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[13,'u'] := AutoActionRec{NextState:= 1,Action:=33};
    SoundAutomaton[13,'v'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[13,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[13,'x'] := AutoActionRec{NextState:=13,Action:= 1};
    SoundAutomaton[13,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[13,'z'] := AutoActionRec{NextState:=12,Action:= 1};
    SoundAutomaton[13,'ç'] := AutoActionRec{NextState:=15,Action:= 1};


    (* State 14 *)
    (* ...q^ *)

    SoundAutomaton[14,'#'] := AutoActionRec{NextState:= 0,Action:= 2};
    SoundAutomaton[14,'a'] := AutoActionRec{NextState:= 1,Action:=52};
    SoundAutomaton[14,'b'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[14,'c'] := AutoActionRec{NextState:= 3,Action:= 1};
    SoundAutomaton[14,'d'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[14,'e'] := AutoActionRec{NextState:= 1,Action:=52};
    SoundAutomaton[14,'f'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[14,'g'] := AutoActionRec{NextState:= 4,Action:= 1};
    SoundAutomaton[14,'h'] := AutoActionRec{NextState:= 2,Action:= 1};
    SoundAutomaton[14,'i'] := AutoActionRec{NextState:= 1,Action:=52};
    SoundAutomaton[14,'j'] := AutoActionRec{NextState:= 5,Action:= 1};
    SoundAutomaton[14,'k'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[14,'l'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[14,'m'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[14,'n'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[14,'o'] := AutoActionRec{NextState:= 1,Action:=52};
    SoundAutomaton[14,'p'] := AutoActionRec{NextState:=27,Action:= 2};
    SoundAutomaton[14,'q'] := AutoActionRec{NextState:=14,Action:= 1};
    SoundAutomaton[14,'r'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[14,'s'] := AutoActionRec{NextState:=17,Action:= 1};
    SoundAutomaton[14,'t'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[14,'u'] := AutoActionRec{NextState:= 1,Action:=52};
    SoundAutomaton[14,'v'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[14,'w'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[14,'x'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[14,'y'] := AutoActionRec{NextState:= 1,Action:= 2};
    SoundAutomaton[14,'z'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[14,'ç'] := AutoActionRec{NextState:=17,Action:= 1};


    (* State 15 *)
    (* ...V^[ç,ss,sç,xç] *)

    SoundAutomaton[15,'#'] := AutoActionRec{NextState:= 0,Action:=29};
    SoundAutomaton[15,'a'] := AutoActionRec{NextState:= 1,Action:=34};
    SoundAutomaton[15,'b'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[15,'c'] := AutoActionRec{NextState:=24,Action:= 1};
    SoundAutomaton[15,'d'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[15,'e'] := AutoActionRec{NextState:= 1,Action:=23};
    SoundAutomaton[15,'f'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[15,'g'] := AutoActionRec{NextState:= 4,Action:=31};
    SoundAutomaton[15,'h'] := AutoActionRec{NextState:= 2,Action:=17};
    SoundAutomaton[15,'i'] := AutoActionRec{NextState:= 1,Action:=23};
    SoundAutomaton[15,'j'] := AutoActionRec{NextState:= 5,Action:=31};
    SoundAutomaton[15,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[15,'l'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[15,'m'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[15,'n'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[15,'o'] := AutoActionRec{NextState:= 1,Action:=35};
    SoundAutomaton[15,'p'] := AutoActionRec{NextState:=27,Action:=31};
    SoundAutomaton[15,'q'] := AutoActionRec{NextState:=14,Action:=29};
    SoundAutomaton[15,'r'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[15,'s'] := AutoActionRec{NextState:=15,Action:= 1};
    SoundAutomaton[15,'t'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[15,'u'] := AutoActionRec{NextState:= 1,Action:=35};
    SoundAutomaton[15,'v'] := AutoActionRec{NextState:= 2,Action:=29};
    SoundAutomaton[15,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[15,'x'] := AutoActionRec{NextState:=15,Action:= 1};
    SoundAutomaton[15,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[15,'z'] := AutoActionRec{NextState:=15,Action:= 1};
    SoundAutomaton[15,'ç'] := AutoActionRec{NextState:=15,Action:= 1};


    (* State 16 *)
    (* ...V^cc *)

    SoundAutomaton[16,'#'] := AutoActionRec{NextState:= 0,Action:= 4};
    SoundAutomaton[16,'a'] := AutoActionRec{NextState:= 1,Action:=53};
    SoundAutomaton[16,'b'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[16,'c'] := AutoActionRec{NextState:=16,Action:= 1};
    SoundAutomaton[16,'d'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[16,'e'] := AutoActionRec{NextState:= 1,Action:=54};
    SoundAutomaton[16,'f'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[16,'g'] := AutoActionRec{NextState:= 4,Action:= 3};
    SoundAutomaton[16,'h'] := AutoActionRec{NextState:= 2,Action:=17};
    SoundAutomaton[16,'i'] := AutoActionRec{NextState:= 1,Action:=54};
    SoundAutomaton[16,'j'] := AutoActionRec{NextState:= 5,Action:= 3};
    SoundAutomaton[16,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[16,'l'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[16,'m'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[16,'n'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[16,'o'] := AutoActionRec{NextState:= 1,Action:=53};
    SoundAutomaton[16,'p'] := AutoActionRec{NextState:=27,Action:= 3};
    SoundAutomaton[16,'q'] := AutoActionRec{NextState:=14,Action:= 2};
    SoundAutomaton[16,'r'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[16,'s'] := AutoActionRec{NextState:=20,Action:= 1};
    SoundAutomaton[16,'t'] := AutoActionRec{NextState:= 2,Action:=24};
    SoundAutomaton[16,'u'] := AutoActionRec{NextState:= 1,Action:=53};
    SoundAutomaton[16,'v'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[16,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[16,'x'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[16,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[16,'z'] := AutoActionRec{NextState:=20,Action:= 1};
    SoundAutomaton[16,'ç'] := AutoActionRec{NextState:=26,Action:= 1};


    (* State 17 *)
    (* ...C^[sç] *)

    SoundAutomaton[17,'#'] := AutoActionRec{NextState:= 0,Action:=37};
    SoundAutomaton[17,'a'] := AutoActionRec{NextState:= 1,Action:=36};
    SoundAutomaton[17,'b'] := AutoActionRec{NextState:= 2,Action:=37};
    SoundAutomaton[17,'c'] := AutoActionRec{NextState:=25,Action:= 1};
    SoundAutomaton[17,'d'] := AutoActionRec{NextState:= 2,Action:=37};
    SoundAutomaton[17,'e'] := AutoActionRec{NextState:= 1,Action:=19};
    SoundAutomaton[17,'f'] := AutoActionRec{NextState:= 2,Action:=37};
    SoundAutomaton[17,'g'] := AutoActionRec{NextState:= 4,Action:= 3};
    SoundAutomaton[17,'h'] := AutoActionRec{NextState:= 2,Action:=17};
    SoundAutomaton[17,'i'] := AutoActionRec{NextState:= 1,Action:=19};
    SoundAutomaton[17,'j'] := AutoActionRec{NextState:= 5,Action:= 3};
    SoundAutomaton[17,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[17,'l'] := AutoActionRec{NextState:= 2,Action:=37};
    SoundAutomaton[17,'m'] := AutoActionRec{NextState:= 2,Action:=37};
    SoundAutomaton[17,'n'] := AutoActionRec{NextState:= 2,Action:=37};
    SoundAutomaton[17,'o'] := AutoActionRec{NextState:= 1,Action:=36};
    SoundAutomaton[17,'p'] := AutoActionRec{NextState:=27,Action:=31};
    SoundAutomaton[17,'q'] := AutoActionRec{NextState:=14,Action:=37};
    SoundAutomaton[17,'r'] := AutoActionRec{NextState:= 2,Action:=37};
    SoundAutomaton[17,'s'] := AutoActionRec{NextState:=17,Action:= 1};
    SoundAutomaton[17,'t'] := AutoActionRec{NextState:= 2,Action:=37};
    SoundAutomaton[17,'u'] := AutoActionRec{NextState:= 1,Action:=36};
    SoundAutomaton[17,'v'] := AutoActionRec{NextState:= 2,Action:=37};
    SoundAutomaton[17,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[17,'x'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[17,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[17,'z'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[17,'ç'] := AutoActionRec{NextState:=17,Action:= 1};


    (* State 18 *)
    (* ...V^u *)

    SoundAutomaton[18,'#'] := AutoActionRec{NextState:= 0,Action:=28};
    SoundAutomaton[18,'a'] := AutoActionRec{NextState:= 1,Action:=38};
    SoundAutomaton[18,'b'] := AutoActionRec{NextState:= 7,Action:=55};
    SoundAutomaton[18,'c'] := AutoActionRec{NextState:= 8,Action:=56};
    SoundAutomaton[18,'d'] := AutoActionRec{NextState:= 7,Action:=55};
    SoundAutomaton[18,'e'] := AutoActionRec{NextState:= 1,Action:=38};
    SoundAutomaton[18,'f'] := AutoActionRec{NextState:= 7,Action:=55};
    SoundAutomaton[18,'g'] := AutoActionRec{NextState:= 4,Action:=56};
    SoundAutomaton[18,'h'] := AutoActionRec{NextState:= 2,Action:=56};
    SoundAutomaton[18,'i'] := AutoActionRec{NextState:= 1,Action:=38};
    SoundAutomaton[18,'j'] := AutoActionRec{NextState:= 5,Action:=56};
    SoundAutomaton[18,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[18,'l'] := AutoActionRec{NextState:= 2,Action:=38};
    SoundAutomaton[18,'m'] := AutoActionRec{NextState:=10,Action:=56};
    SoundAutomaton[18,'n'] := AutoActionRec{NextState:=10,Action:=56};
    SoundAutomaton[18,'o'] := AutoActionRec{NextState:= 1,Action:=38};
    SoundAutomaton[18,'p'] := AutoActionRec{NextState:=21,Action:=56};
    SoundAutomaton[18,'q'] := AutoActionRec{NextState:=14,Action:=55};
    SoundAutomaton[18,'r'] := AutoActionRec{NextState:=11,Action:=55};
    SoundAutomaton[18,'s'] := AutoActionRec{NextState:=12,Action:=56};
    SoundAutomaton[18,'t'] := AutoActionRec{NextState:= 7,Action:=55};
    SoundAutomaton[18,'u'] := AutoActionRec{NextState:=18,Action:= 1};
    SoundAutomaton[18,'v'] := AutoActionRec{NextState:= 2,Action:=55};
    SoundAutomaton[18,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[18,'x'] := AutoActionRec{NextState:=13,Action:=56};
    SoundAutomaton[18,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[18,'z'] := AutoActionRec{NextState:=12,Action:=56};
    SoundAutomaton[18,'ç'] := AutoActionRec{NextState:=15,Action:=56};


    (* State 19 *)
    (* ...V[bdfpt]^[ei] *)

    SoundAutomaton[19,'#'] := AutoActionRec{NextState:= 0,Action:=40};
    SoundAutomaton[19,'a'] := AutoActionRec{NextState:= 1,Action:=40};
    SoundAutomaton[19,'b'] := AutoActionRec{NextState:= 2,Action:=22};
    SoundAutomaton[19,'c'] := AutoActionRec{NextState:= 8,Action:=21};
    SoundAutomaton[19,'d'] := AutoActionRec{NextState:= 2,Action:=22};
    SoundAutomaton[19,'e'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[19,'f'] := AutoActionRec{NextState:= 2,Action:=22};
    SoundAutomaton[19,'g'] := AutoActionRec{NextState:= 4,Action:=21};
    SoundAutomaton[19,'h'] := AutoActionRec{NextState:= 1,Action:=49};
    SoundAutomaton[19,'i'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[19,'j'] := AutoActionRec{NextState:= 5,Action:=21};
    SoundAutomaton[19,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[19,'l'] := AutoActionRec{NextState:= 9,Action:=21};
    SoundAutomaton[19,'m'] := AutoActionRec{NextState:=10,Action:=49};
    SoundAutomaton[19,'n'] := AutoActionRec{NextState:=10,Action:=49};
    SoundAutomaton[19,'o'] := AutoActionRec{NextState:= 1,Action:=40};
    SoundAutomaton[19,'p'] := AutoActionRec{NextState:=21,Action:=21};
    SoundAutomaton[19,'q'] := AutoActionRec{NextState:=14,Action:=22};
    SoundAutomaton[19,'r'] := AutoActionRec{NextState:= 2,Action:=22};
    SoundAutomaton[19,'s'] := AutoActionRec{NextState:=12,Action:=21};
    SoundAutomaton[19,'t'] := AutoActionRec{NextState:= 2,Action:=22};
    SoundAutomaton[19,'u'] := AutoActionRec{NextState:=18,Action:=49};
    SoundAutomaton[19,'v'] := AutoActionRec{NextState:= 2,Action:=22};
    SoundAutomaton[19,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[19,'x'] := AutoActionRec{NextState:=13,Action:=49};
    SoundAutomaton[19,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[19,'z'] := AutoActionRec{NextState:=12,Action:=49};
    SoundAutomaton[19,'ç'] := AutoActionRec{NextState:=15,Action:=21};


    (* State 20 *)
    (* ...V^[cs,cz] *)

    SoundAutomaton[20,'#'] := AutoActionRec{NextState:= 0,Action:=42};
    SoundAutomaton[20,'a'] := AutoActionRec{NextState:= 1,Action:=42};
    SoundAutomaton[20,'b'] := AutoActionRec{NextState:= 2,Action:=42};
    SoundAutomaton[20,'c'] := AutoActionRec{NextState:=24,Action:= 1};
    SoundAutomaton[20,'d'] := AutoActionRec{NextState:= 2,Action:=42};
    SoundAutomaton[20,'e'] := AutoActionRec{NextState:= 1,Action:=54};
    SoundAutomaton[20,'f'] := AutoActionRec{NextState:= 2,Action:=42};
    SoundAutomaton[20,'g'] := AutoActionRec{NextState:= 4,Action:= 5};
    SoundAutomaton[20,'h'] := AutoActionRec{NextState:=20,Action:= 1};
    SoundAutomaton[20,'i'] := AutoActionRec{NextState:= 1,Action:=54};
    SoundAutomaton[20,'j'] := AutoActionRec{NextState:= 5,Action:= 5};
    SoundAutomaton[20,'k'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[20,'l'] := AutoActionRec{NextState:= 2,Action:=42};
    SoundAutomaton[20,'m'] := AutoActionRec{NextState:= 2,Action:=42};
    SoundAutomaton[20,'n'] := AutoActionRec{NextState:= 2,Action:=42};
    SoundAutomaton[20,'o'] := AutoActionRec{NextState:= 1,Action:=42};
    SoundAutomaton[20,'p'] := AutoActionRec{NextState:=27,Action:=31};
    SoundAutomaton[20,'q'] := AutoActionRec{NextState:=14,Action:=42};
    SoundAutomaton[20,'r'] := AutoActionRec{NextState:= 2,Action:=42};
    SoundAutomaton[20,'s'] := AutoActionRec{NextState:=20,Action:= 1};
    SoundAutomaton[20,'t'] := AutoActionRec{NextState:= 2,Action:=42};
    SoundAutomaton[20,'u'] := AutoActionRec{NextState:= 1,Action:=42};
    SoundAutomaton[20,'v'] := AutoActionRec{NextState:= 2,Action:=42};
    SoundAutomaton[20,'w'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[20,'x'] := AutoActionRec{NextState:=13,Action:= 1};
    SoundAutomaton[20,'y'] := AutoActionRec{NextState:= 1,Action:= 6};
    SoundAutomaton[20,'z'] := AutoActionRec{NextState:=20,Action:= 1};
    SoundAutomaton[20,'ç'] := AutoActionRec{NextState:=26,Action:= 1};


    (* State 21 *)
    (* ^p,...V^p *)

    SoundAutomaton[21,'#'] := AutoActionRec{NextState:= 0,Action:= 4};
    SoundAutomaton[21,'a'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[21,'b'] := AutoActionRec{NextState:= 7,Action:=45};
    SoundAutomaton[21,'c'] := AutoActionRec{NextState:= 3,Action:=46};
    SoundAutomaton[21,'d'] := AutoActionRec{NextState:= 7,Action:=45};
    SoundAutomaton[21,'e'] := AutoActionRec{NextState:=19,Action:= 3};
    SoundAutomaton[21,'f'] := AutoActionRec{NextState:= 7,Action:=45};
    SoundAutomaton[21,'g'] := AutoActionRec{NextState:= 4,Action:=46};
    SoundAutomaton[21,'h'] := AutoActionRec{NextState:= 7,Action:=48};
    SoundAutomaton[21,'i'] := AutoActionRec{NextState:=19,Action:= 3};
    SoundAutomaton[21,'j'] := AutoActionRec{NextState:= 5,Action:=46};
    SoundAutomaton[21,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[21,'l'] := AutoActionRec{NextState:= 2,Action:=45};
    SoundAutomaton[21,'m'] := AutoActionRec{NextState:= 2,Action:=45};
    SoundAutomaton[21,'n'] := AutoActionRec{NextState:= 2,Action:=45};
    SoundAutomaton[21,'o'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[21,'p'] := AutoActionRec{NextState:=21,Action:= 1};
    SoundAutomaton[21,'q'] := AutoActionRec{NextState:=14,Action:=45};
    SoundAutomaton[21,'r'] := AutoActionRec{NextState:= 2,Action:=45};
    SoundAutomaton[21,'s'] := AutoActionRec{NextState:=17,Action:= 3};
    SoundAutomaton[21,'t'] := AutoActionRec{NextState:= 2,Action:=44};
    SoundAutomaton[21,'u'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[21,'v'] := AutoActionRec{NextState:= 2,Action:=45};
    SoundAutomaton[21,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[21,'x'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[21,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[21,'z'] := AutoActionRec{NextState:= 2,Action:=45};
    SoundAutomaton[21,'ç'] := AutoActionRec{NextState:=17,Action:=50};


    (* State 22 *)
    (* ...V[mn]^x *)

    SoundAutomaton[22,'#'] := AutoActionRec{NextState:= 0,Action:= 4};
    SoundAutomaton[22,'a'] := AutoActionRec{NextState:= 1,Action:=43};
    SoundAutomaton[22,'b'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[22,'c'] := AutoActionRec{NextState:= 3,Action:= 3};
    SoundAutomaton[22,'d'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[22,'e'] := AutoActionRec{NextState:= 1,Action:=43};
    SoundAutomaton[22,'f'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[22,'g'] := AutoActionRec{NextState:= 4,Action:= 3};
    SoundAutomaton[22,'h'] := AutoActionRec{NextState:= 2,Action:= 3};
    SoundAutomaton[22,'i'] := AutoActionRec{NextState:= 1,Action:=43};
    SoundAutomaton[22,'j'] := AutoActionRec{NextState:= 5,Action:= 3};
    SoundAutomaton[22,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[22,'l'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[22,'m'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[22,'n'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[22,'o'] := AutoActionRec{NextState:= 1,Action:=43};
    SoundAutomaton[22,'p'] := AutoActionRec{NextState:=27,Action:= 3};
    SoundAutomaton[22,'q'] := AutoActionRec{NextState:=14,Action:= 4};
    SoundAutomaton[22,'r'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[22,'s'] := AutoActionRec{NextState:=17,Action:= 1};
    SoundAutomaton[22,'t'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[22,'u'] := AutoActionRec{NextState:= 1,Action:=43};
    SoundAutomaton[22,'v'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[22,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[22,'x'] := AutoActionRec{NextState:=22,Action:= 1};
    SoundAutomaton[22,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[22,'z'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[22,'ç'] := AutoActionRec{NextState:=17,Action:= 1};


    (* State 23 *)
    (* ...Vr^r *)

    SoundAutomaton[23,'#'] := AutoActionRec{NextState:= 0,Action:= 2};
    SoundAutomaton[23,'a'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[23,'b'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[23,'c'] := AutoActionRec{NextState:= 3,Action:= 1};
    SoundAutomaton[23,'d'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[23,'e'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[23,'f'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[23,'g'] := AutoActionRec{NextState:= 4,Action:= 1};
    SoundAutomaton[23,'h'] := AutoActionRec{NextState:=23,Action:= 1};
    SoundAutomaton[23,'i'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[23,'j'] := AutoActionRec{NextState:= 5,Action:= 1};
    SoundAutomaton[23,'k'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[23,'l'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[23,'m'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[23,'n'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[23,'o'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[23,'p'] := AutoActionRec{NextState:=27,Action:= 1};
    SoundAutomaton[23,'q'] := AutoActionRec{NextState:=14,Action:= 2};
    SoundAutomaton[23,'r'] := AutoActionRec{NextState:=23,Action:= 1};
    SoundAutomaton[23,'s'] := AutoActionRec{NextState:=17,Action:= 1};
    SoundAutomaton[23,'t'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[23,'u'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[23,'v'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[23,'w'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[23,'x'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[23,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[23,'z'] := AutoActionRec{NextState:= 2,Action:= 2};
    SoundAutomaton[23,'ç'] := AutoActionRec{NextState:=17,Action:= 1};


    (* State 24 *)
    (* ...V^[sc,xc,zc] *)

    SoundAutomaton[24,'#'] := AutoActionRec{NextState:= 0,Action:= 6};
    SoundAutomaton[24,'a'] := AutoActionRec{NextState:= 1,Action:=47};
    SoundAutomaton[24,'b'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[24,'c'] := AutoActionRec{NextState:=24,Action:= 1};
    SoundAutomaton[24,'d'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[24,'e'] := AutoActionRec{NextState:= 1,Action:=23};
    SoundAutomaton[24,'f'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[24,'g'] := AutoActionRec{NextState:= 4,Action:= 5};
    SoundAutomaton[24,'h'] := AutoActionRec{NextState:= 2,Action:=17};
    SoundAutomaton[24,'i'] := AutoActionRec{NextState:= 1,Action:=23};
    SoundAutomaton[24,'j'] := AutoActionRec{NextState:= 5,Action:= 5};
    SoundAutomaton[24,'k'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[24,'l'] := AutoActionRec{NextState:= 2,Action:=47};
    SoundAutomaton[24,'m'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[24,'n'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[24,'o'] := AutoActionRec{NextState:= 1,Action:=47};
    SoundAutomaton[24,'p'] := AutoActionRec{NextState:=27,Action:= 5};
    SoundAutomaton[24,'q'] := AutoActionRec{NextState:=14,Action:=29};
    SoundAutomaton[24,'r'] := AutoActionRec{NextState:= 2,Action:=47};
    SoundAutomaton[24,'s'] := AutoActionRec{NextState:=20,Action:= 1};
    SoundAutomaton[24,'t'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[24,'u'] := AutoActionRec{NextState:= 1,Action:=47};
    SoundAutomaton[24,'v'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[24,'w'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[24,'x'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[24,'y'] := AutoActionRec{NextState:= 1,Action:= 6};
    SoundAutomaton[24,'z'] := AutoActionRec{NextState:=20,Action:= 1};
    SoundAutomaton[24,'ç'] := AutoActionRec{NextState:=15,Action:= 1};


    (* State 25 *)
    (* ...C^sc *)

    SoundAutomaton[25,'#'] := AutoActionRec{NextState:= 0,Action:= 6};
    SoundAutomaton[25,'a'] := AutoActionRec{NextState:= 1,Action:= 6};
    SoundAutomaton[25,'b'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[25,'c'] := AutoActionRec{NextState:=25,Action:= 1};
    SoundAutomaton[25,'d'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[25,'e'] := AutoActionRec{NextState:= 1,Action:=19};
    SoundAutomaton[25,'f'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[25,'g'] := AutoActionRec{NextState:= 4,Action:= 5};
    SoundAutomaton[25,'h'] := AutoActionRec{NextState:= 2,Action:=17};
    SoundAutomaton[25,'i'] := AutoActionRec{NextState:= 1,Action:=19};
    SoundAutomaton[25,'j'] := AutoActionRec{NextState:= 5,Action:= 5};
    SoundAutomaton[25,'k'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[25,'l'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[25,'m'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[25,'n'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[25,'o'] := AutoActionRec{NextState:= 1,Action:= 9};
    SoundAutomaton[25,'p'] := AutoActionRec{NextState:=27,Action:= 5};
    SoundAutomaton[25,'q'] := AutoActionRec{NextState:=14,Action:=29};
    SoundAutomaton[25,'r'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[25,'s'] := AutoActionRec{NextState:=17,Action:= 1};
    SoundAutomaton[25,'t'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[25,'u'] := AutoActionRec{NextState:= 1,Action:= 9};
    SoundAutomaton[25,'v'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[25,'w'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[25,'x'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[25,'y'] := AutoActionRec{NextState:= 1,Action:= 6};
    SoundAutomaton[25,'z'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[25,'ç'] := AutoActionRec{NextState:=17,Action:= 1};


    (* State 26 *)
    (* ...V^cç *)

    SoundAutomaton[26,'#'] := AutoActionRec{NextState:= 0,Action:= 6};
    SoundAutomaton[26,'a'] := AutoActionRec{NextState:= 1,Action:=41};
    SoundAutomaton[26,'b'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[26,'c'] := AutoActionRec{NextState:= 8,Action:= 1};
    SoundAutomaton[26,'d'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[26,'e'] := AutoActionRec{NextState:= 1,Action:=23};
    SoundAutomaton[26,'f'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[26,'g'] := AutoActionRec{NextState:= 4,Action:= 5};
    SoundAutomaton[26,'h'] := AutoActionRec{NextState:= 2,Action:=17};
    SoundAutomaton[26,'i'] := AutoActionRec{NextState:= 1,Action:=23};
    SoundAutomaton[26,'j'] := AutoActionRec{NextState:= 5,Action:= 5};
    SoundAutomaton[26,'k'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[26,'l'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[26,'m'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[26,'n'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[26,'o'] := AutoActionRec{NextState:= 1,Action:=41};
    SoundAutomaton[26,'p'] := AutoActionRec{NextState:=27,Action:= 5};
    SoundAutomaton[26,'q'] := AutoActionRec{NextState:=14,Action:=29};
    SoundAutomaton[26,'r'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[26,'s'] := AutoActionRec{NextState:=20,Action:= 1};
    SoundAutomaton[26,'t'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[26,'u'] := AutoActionRec{NextState:= 1,Action:=41};
    SoundAutomaton[26,'v'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[26,'w'] := AutoActionRec{NextState:= 2,Action:= 6};
    SoundAutomaton[26,'x'] := AutoActionRec{NextState:=13,Action:= 1};
    SoundAutomaton[26,'y'] := AutoActionRec{NextState:= 1,Action:= 6};
    SoundAutomaton[26,'z'] := AutoActionRec{NextState:=20,Action:= 1};
    SoundAutomaton[26,'ç'] := AutoActionRec{NextState:=26,Action:= 1};


    (* State 27 *)
    (* ...C^p *)

    SoundAutomaton[27,'#'] := AutoActionRec{NextState:= 0,Action:= 4};
    SoundAutomaton[27,'a'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[27,'b'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[27,'c'] := AutoActionRec{NextState:= 3,Action:= 3};
    SoundAutomaton[27,'d'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[27,'e'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[27,'f'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[27,'g'] := AutoActionRec{NextState:= 4,Action:= 3};
    SoundAutomaton[27,'h'] := AutoActionRec{NextState:= 2,Action:=48};
    SoundAutomaton[27,'i'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[27,'j'] := AutoActionRec{NextState:= 5,Action:= 3};
    SoundAutomaton[27,'k'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[27,'l'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[27,'m'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[27,'n'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[27,'o'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[27,'p'] := AutoActionRec{NextState:=27,Action:= 1};
    SoundAutomaton[27,'q'] := AutoActionRec{NextState:=14,Action:= 4};
    SoundAutomaton[27,'r'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[27,'s'] := AutoActionRec{NextState:=17,Action:= 3};
    SoundAutomaton[27,'t'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[27,'u'] := AutoActionRec{NextState:= 1,Action:= 8};
    SoundAutomaton[27,'v'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[27,'w'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[27,'x'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[27,'y'] := AutoActionRec{NextState:= 1,Action:= 4};
    SoundAutomaton[27,'z'] := AutoActionRec{NextState:= 2,Action:= 4};
    SoundAutomaton[27,'ç'] := AutoActionRec{NextState:=17,Action:= 3};
    
  END InitSoundex;


  PROCEDURE Soundex(READONLY InWord: REF ChArr; InLen: CARDINAL; 
                    OutWord: REF ChArr)  =
                    
    CONST
      SetEI = SET OF CHAR{'e','i'};
      SetOU = SET OF CHAR{'o','u'};
    
    VAR
      CurrState: AutoState := 0;
      InPos: CARDINAL := 0;
      CleanWord := NEW(REF ChArr,InLen+1);
      OutLen: CARDINAL;
      
    PROCEDURE Last() =
    BEGIN
      OutWord[OutLen] := CleanWord[InPos];
      INC(OutLen);
    END Last;
  
    PROCEDURE BLast() =
    BEGIN
      OutWord[OutLen] := CleanWord[InPos-1];
      INC(OutLen);
    END BLast;
  
    PROCEDURE BBLast() =
    BEGIN
      OutWord[OutLen] := CleanWord[InPos-2];
      INC(OutLen);
    END BBLast;
  
    PROCEDURE CondLast() =
    BEGIN
      WITH
        ch=CleanWord[InPos],
        och=OutWord[OutLen]
      DO
        IF ch IN SetEI THEN
          och := SoundEI
        ELSIF ch IN SetOU THEN
          och := SoundOU
        ELSE
          och := ch
        END
      END;
      INC(OutLen)
    END CondLast;
      
    PROCEDURE InsertSound(s: CHAR) =
    BEGIN
      OutWord[OutLen] := s;
      INC(OutLen)
    END InsertSound;
      
    PROCEDURE InsertSoundLast(s: CHAR) =
    BEGIN
      InsertSound(s);
      Last()
    END InsertSoundLast;
    
    PROCEDURE CondInsertSoundLast(s: CHAR) =
    BEGIN
      InsertSound(s);
      CondLast()
    END CondInsertSoundLast;
        
      
  BEGIN (* Soundex *)
  
    OutLen := 0;
    FOR i:=0 TO InLen-1 DO
      CleanWord^[i] := StripAccents[InWord^[i]];
    END;
    CleanWord^[InLen] := '#';
    
    WHILE InPos <= InLen DO
      WITH
        actrec=SoundAutomaton[CurrState,CleanWord^[InPos]],
        nst=actrec.NextState,
        act=actrec.Action
      DO
        CASE act OF 
           0  => <* ASSERT FALSE *>  (* impossible *)
        |  1  => (* Do nothing *)
        |  2  => Last()
        |  3  => BLast()
        |  4  => BLast(); Last()
        |  5  => BBLast(); BLast()
        |  6  => BBLast(); BLast(); Last()
        |  7  => CondLast()
        |  8  => BLast(); CondLast()
        |  9  => BBLast(); BLast(); CondLast()
        | 10  => IF CleanWord[InPos-1]#CleanWord[InPos] THEN BLast() END;
        | 11  => IF CleanWord[InPos-1]#CleanWord[InPos] THEN Last() END;
        | 12  => IF CleanWord[InPos-1]#CleanWord[InPos] THEN 
                   InsertSoundLast(XSoundEI)
                 END
        | 13  => <* ASSERT FALSE *> (* for future use *)
        | 14  => <* ASSERT FALSE *> (* for future use *)
        | 15  => <* ASSERT FALSE *> (* for future use *)
        | 16  => IF OutLen=0 THEN 
                   CondInsertSoundLast(XSoundH) 
                 ELSE
                   CondLast()
                 END
        | 17  => InsertSound(SoundChX)
        | 18  => InsertSoundLast(XSoundPC)
        | 19  => InsertSound(SoundSCSc); InsertSound(SoundEI)
        | 20  => InsertSound(SoundGJ); InsertSound(SoundEI)
        | 21  => InsertSound(XSoundEI)
        | 22  => InsertSoundLast(XSoundEI)
        | 23  => InsertSound(SoundCScSsXcXCc); InsertSound(SoundEI)
        | 24  => InsertSoundLast(XSoundC)
        | 25  => BLast(); InsertSound(SoundChX)
        | 26  => InsertSound(SoundMN)
	| 27  => InsertSoundLast(SoundMN) 
        | 28  => InsertSoundLast(SoundLU)
        | 29  => InsertSoundLast(SoundSXZ)
        | 30  => CondInsertSoundLast(SoundSXZ)
        | 31  => InsertSound(SoundSXZ)
        | 32  => InsertSoundLast(SoundChSXZ)
        | 33  => CondInsertSoundLast(SoundChSXZ)
        | 34  => InsertSoundLast(SoundCedScedSsCcedPced)
        | 35  => CondInsertSoundLast(SoundCedScedSsCcedPced)
        | 36  => CondInsertSoundLast(SoundCedS)
        | 37  => CondInsertSoundLast(SoundS)
        | 38  => CondInsertSoundLast(SoundOU)
        | 39  => InsertSound(SoundLU)
        | 40  => CondInsertSoundLast(SoundEI)
        | 41  => CondInsertSoundLast(SoundCedCced)
        | 42  => CondInsertSoundLast(SoundX)
        | 43  => CondInsertSoundLast(SoundChX)
        | 44  => InsertSound(XSoundP); InsertSoundLast(XSoundEI)
        | 45  => BLast(); InsertSoundLast(XSoundEI)            
        | 46  => BLast(); InsertSound(XSoundEI)        
        | 47  => CondInsertSoundLast(SoundScXc)
        | 48  => InsertSound(SoundF)
        | 49  => InsertSound(SoundEI)
        | 50  => InsertSound(XSoundP)
        | 51  => InsertSound(XSoundH)
        | 52  => CondInsertSoundLast(XSoundU)
        | 53  => CondInsertSoundLast(SoundCCed)
        | 54  => CondInsertSoundLast(SoundCCcX)
        | 55  => InsertSoundLast(SoundLOU)
        | 56  => InsertSound(SoundLOU)
        END (* case *);
        
      CurrState := nst
      
      END (* with *);
      
      INC(InPos)
      
    END (* while *);
    
    <* ASSERT OutWord[OutLen-1]='#' *>
    
    (* quick fix for the accented forms of 'a' not present in the original 
       automaton! *)
       
    FOR i := 0 TO OutLen-2 DO
      IF OutWord[i]='a' THEN OutWord[i] := SoundA END
    END
    
  END Soundex;
  
  PROCEDURE AlternativeSpellings(Aut: ReadOnlyPermDAG.T; Tree: AVLTextTree.T;
                                 READONLY InWord: REF ChArr; 
                                 VAR (*IO*) AltTree: AVLTextTree.T) RAISES{Abort}=
                                 
    VAR
      root := Aut.Root();
      TreeRoot := Tree.RootState();
      OutWord := NEW(REF ChArr,100);
      OutLen: CARDINAL := 0;
      
    PROCEDURE OutputWord() =
    BEGIN 
      AltTree.Insert(Text.FromChars(SUBARRAY(OutWord^,0,OutLen)));
    END OutputWord;
    
    PROCEDURE TryAutSubst(InPos: CARDINAL; LastState: State)  RAISES{Abort} =
      VAR
        final: BOOLEAN;
        CurrState: State;
        LocLastState: State := LastState;
        LocOutLen: CARDINAL := OutLen;
        chars: REF ARRAY OF CHAR;
        
      PROCEDURE StepLetter(cch: CHAR)  RAISES {NoTransition,Abort} =
        VAR
          chars: REF ARRAY OF CHAR;
      BEGIN
        Aut.ExpandState(LocLastState,final,chars);
        FOR i := 0 TO LAST(chars^) DO
          IF cch=chars^[i] THEN
            OutWord[OutLen] := cch;
            INC(OutLen);
            LocLastState := Aut.MakeTransition(LocLastState,cch);
            RETURN 
          END
        END;
        RAISE NoTransition
      END StepLetter;

    BEGIN (* TryAutSubst *)
      Aut.ExpandState(LastState,final,chars);
      WITH
        ch=InWord[InPos]
      DO 
        IF ch='#' THEN 
          IF final THEN
            OutputWord()
          END;
          RETURN
        END;
        IF ch IN Letters  THEN 
          FOR i := 0 TO LAST(chars^) DO
            IF ch=chars^[i] THEN
              OutWord[OutLen] := ch;
              INC(OutLen);
              CurrState := Aut.MakeTransition(LastState,ch);
              TryAutSubst(InPos+1,CurrState);
              OutLen := LocOutLen;
              EXIT
            END;
          END
        ELSIF ch IN Sounds THEN 
          WITH
            arrsub=SoundSubst[ch]
          DO
            FOR s := 0 TO MaxSubst-1 DO
              WITH
                sub=arrsub[s]
              DO
                IF s>0 AND sub=NoSub THEN RETURN  END;
                TRY
                  FOR j := 0 TO MaxSubstLetters-1 DO
                    WITH
                      cch=sub[j]
                    DO
                      IF cch=NUL THEN EXIT END; 
                      StepLetter(cch)
                    END
                  END;
                  TryAutSubst(InPos+1,LocLastState)
                EXCEPT 
                  NoTransition => (* OK *)
                END
              END;
              OutLen := LocOutLen;
              LocLastState := LastState
            END
          END
        ELSE <* ASSERT FALSE *> (* Impossible *)
        END;
      END
    END TryAutSubst;

    PROCEDURE TryTreeSubst(InPos: CARDINAL; 
                           LastState: AVLTextTree.AVLTextTreeState)  =
    
      VAR
        CurrState: AVLTextTree.AVLTextTreeState;
        LocLastState: AVLTextTree.AVLTextTreeState := LastState;
        LocOutLen: CARDINAL := OutLen;
        
      PROCEDURE StepLetter(cch: CHAR)  RAISES {NoTransition} =
        VAR
          aux: AVLTextTree.AVLTextTreeState;
      BEGIN
        aux := Tree.MakeTransition(LocLastState,cch);
        IF NOT AVLTextTree.IsNullState(aux) THEN
          OutWord[OutLen] := cch;
          INC(OutLen);
          LocLastState := aux;
          RETURN
        ELSE 
          RAISE NoTransition
        END
      END StepLetter;

    BEGIN (* TryTreeSubst *)
      WITH
        ch=InWord[InPos]
      DO 
        IF ch='#' THEN 
          IF Tree.IsFinal(LastState) THEN
            OutputWord()
          END;
          RETURN
        END;
        IF ch IN Letters THEN
          CurrState := Tree.MakeTransition(LastState,ch);
          IF NOT AVLTextTree.IsNullState(CurrState) THEN
            OutWord[OutLen] := ch;
            INC(OutLen);
            TryTreeSubst(InPos+1,CurrState);
            OutLen := LocOutLen;
          END
        ELSIF ch IN Sounds THEN 
          WITH
            arrsub=SoundSubst[ch]
          DO
            FOR s := 0 TO MaxSubst-1 DO
              WITH
                sub=arrsub[s]
              DO
                IF s>0 AND sub=NoSub THEN RETURN  END;
                TRY
                  FOR j := 0 TO MaxSubstLetters-1 DO
                    WITH
                      cch=sub[j]
                    DO
                      IF cch=NUL THEN EXIT END; 
                      StepLetter(cch)
                    END
                  END;
                  TryTreeSubst(InPos+1,LocLastState)
                EXCEPT 
                  NoTransition => (* OK *)
                END
              END;
              OutLen := LocOutLen;
              LocLastState := LastState
            END
          END
        ELSE <* ASSERT FALSE *> (* Impossible *)
        END;
      END
    END TryTreeSubst;

  BEGIN
    TryAutSubst(0,root);
    TryTreeSubst(0,TreeRoot)
  END AlternativeSpellings;                                 

  PROCEDURE Alternatives(READONLY w: TEXT;
                             Aut: ReadOnlyPermDAG.T;
                             LocalTree: AVLTextTree.T;
                             VAR (*IO*) AltTree: AVLTextTree.T;
                             Transpose,
                             Remove,
                             Subst,
                             Insert,
                             ForceAlternatives: BOOLEAN) RAISES{Abort} =
  BEGIN
    AltTree := AVLTextTree.New();
    WITH n = Text.Length(w) DO
      Text.SetChars(InWord^,w);
      Soundex(InWord,n,OutWord);
      AlternativeSpellings(Aut,LocalTree,OutWord,AltTree);
      IF  (n>=minTransp AND Transpose)
          OR ForceAlternatives THEN  (* Try transposing two consecutive letters *)
        FOR i := 1 TO n-3 DO
          SUBARRAY(InWord2^,0,n) := SUBARRAY(InWord^,0,n);
          InWord2[i] := InWord[i+1];
          InWord2[i+1] := InWord[i];
          Soundex(InWord2,n,OutWord);
          AlternativeSpellings(Aut,LocalTree,OutWord,AltTree);
        END
      END;
      IF  (n>=minRemove AND Remove) OR
           ForceAlternatives THEN  (* Try removing one letter at a time *)
        FOR i := 1 TO n-2 DO
          SUBARRAY(InWord2^,0,i) := SUBARRAY(InWord^,0,i);
          SUBARRAY(InWord2^,i,n-i-1) := SUBARRAY(InWord^,i+1,n-i-1);
          Soundex(InWord2,n-1,OutWord);
          AlternativeSpellings(Aut,LocalTree,OutWord,AltTree);
        END
      END;
      IF (n>=minSubst AND Subst) OR 
          ForceAlternatives THEN (* Try substitutions *)
        FOR i := 1 TO n-2 DO
          SUBARRAY(InWord2^,0,n) := SUBARRAY(InWord^,0,n);
          FOR k := FIRST(InsertionLetters) TO LAST(InsertionLetters) DO
            InWord2[i] := InsertionLetters[k];
            Soundex(InWord2,n,OutWord);
            AlternativeSpellings(Aut,LocalTree,OutWord,AltTree);
          END
        END
      END;
      IF (n>=minInsert AND Insert) OR 
          ForceAlternatives THEN (* Try insertions *)
        FOR i := 1 TO n-1 DO
          SUBARRAY(InWord2^,0,i) := SUBARRAY(InWord^,0,i);
          SUBARRAY(InWord2^,i+1,n-i) := SUBARRAY(InWord^,i,n-i);
          FOR k := FIRST(InsertionLetters) TO LAST(InsertionLetters) DO
            InWord2[i] := InsertionLetters[k];
            Soundex(InWord2,n+1,OutWord);
            AlternativeSpellings(Aut,LocalTree,OutWord,AltTree);
          END
        END
      END      
    END
  END Alternatives;                             
  


BEGIN (* Initilize module *)
    InitSoundex();
    InitAlternatives()
END Adviser.
