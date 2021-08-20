MODULE AVLTextTree;

IMPORT Text;

FROM Basics IMPORT Abort;

  REVEAL
  
    T = AVLPublic BRANDED OBJECT
      Tree: AVLRef;
      Nodes: CARDINAL
    OVERRIDES
      Insert := Insert;
      Search := Search;
      Delete := Delete;
      Size := Size;
      Enum := Enum;
      MakeTransition := AVLMakeTransition;
      IsFinal := AVLIsFinal;
      RootState := RootState;
    END;
    
    AVLTextTreeState = BRANDED REF RECORD
        Tree: AVLRef;
        Position: CARDINAL
      END;
      
  VAR 
    AVLTextTreeNullState: AVLTextTreeState;
      
  TYPE 
  
    BalFactor = [-1..+1];
  
    AVLRef = REF AVLRec;
    
    AVLRec = RECORD
        bal: BalFactor;
        left,right: AVLRef;
        info:  TEXT
      END;
      
    Result = {Smaller,Equal,Greater};
 
  PROCEDURE Insert(tr: T;  
                   READONLY w: TEXT) =
    VAR DummyBool: BOOLEAN;
    PROCEDURE DoInsert(VAR (*IO*) t: AVLRef;
                       VAR (*OUT*) h: BOOLEAN) =
      VAR t1,t2: AVLRef;
    BEGIN
      IF t=NIL  THEN
        t := NEW(AVLRef);
        t^:= AVLRec{bal:=0,left:=NIL,right:=NIL,info:=w};
        h := TRUE;
        INC(tr.Nodes);
        RETURN 
      ELSE
        CASE CompText(w,t^.info)  OF
          Result.Smaller =>  
                  DoInsert(t^.left,h); 
                  IF h THEN (* left branch has grown *)
                    CASE t^.bal OF
                       1 => t^.bal := 0;  h := FALSE
                    |  0 => t^.bal := -1
                    | -1 =>  (* rebalance *)
                         t1 := t^.left;
                         IF t1^.bal=-1 THEN (* single LL *)
                           t^.left := t1^.right; t1^.right := t;
                           t^.bal := 0; t := t1
                         ELSE (* double LR *)
                           t2 := t1^.right;
                           t1^.right := t2^.left; t2^.left := t1;
                           t^.left := t2^.right; t2^.right := t;
                           IF t2^.bal=-1 THEN 
                             t^.bal := 1 
                           ELSE
                             t^.bal := 0
                           END;
                           IF t2^.bal=+1 THEN
                             t1^.bal := -1
                           ELSE
                             t1^.bal := 0;
                           END;
                           t := t2
                         END;
                         t^.bal := 0; h := FALSE;
                    END
                  END;
                  RETURN
        |  Result.Equal =>
                  h := FALSE;
                  RETURN
        |  Result.Greater =>  
                  DoInsert(t^.right,h);
                  IF h THEN (* right branch has grown *)
                    CASE t^.bal OF
                      -1 => t^.bal := 0;  h := FALSE
                    |  0 => t^.bal := +1
                    |  1 =>  (* rebalance *)
                           t1 := t^.right;
                           IF t1^.bal=+1 THEN (* single RR *)
                             t^.right := t1^.left; t1^.left := t;
                             t^.bal := 0; t := t1
                           ELSE (* double RL *)
                             t2 := t1^.left;
                             t1^.left := t2^.right; t2^.right := t1;
                             t^.right := t2^.left; t2^.left := t;
                             IF t2^.bal=+1 THEN 
                               t^.bal := -1 
                             ELSE
                               t^.bal := 0
                             END;
                             IF t2^.bal=-1 THEN
                               t1^.bal := +1
                             ELSE
                               t1^.bal := 0;
                             END;
                             t := t2
                           END;
                           t^.bal := 0; h := FALSE;
                    END
                  END;
                  RETURN
        END
      END
    END DoInsert;
  BEGIN (* Insert*)
    DoInsert(tr.Tree,DummyBool)
  END Insert;

(*
  PROCEDURE Delete(<* UNUSED *> tr: T; <* UNUSED *> READONLY w: TEXT): BOOLEAN =
  BEGIN
    <* ASSERT FALSE *> (* Delete not implemented yet! *)
  END Delete;
*)

  PROCEDURE Delete(tr: T; READONLY w: TEXT) =
    VAR
      DummyBool: BOOLEAN;
    PROCEDURE DoDelete(VAR (*IO*) t: AVLRef;
                       VAR (*OUT*) h: BOOLEAN) =
      VAR
        q: AVLRef;

      PROCEDURE BalanceLeft(VAR t: AVLRef; VAR h: BOOLEAN) =
        VAR 
          t1,t2: AVLRef;
          b1,b2: BalFactor;
      BEGIN
        CASE t^.bal OF
          -1 => t^.bal := 0;
        |  0 => t^.bal := +1; h := FALSE;
        | +1 => (* rebalance *)
                t1 := t^.right; b1 := t1^.bal;
                IF b1>=0 THEN (* single RR *)
                  t^.right := t1^.left; t1^.left := t;
                  IF b1=0 THEN
                    t^.bal := +1; t1^.bal := -1; h := FALSE 
                  ELSE 
                    t^.bal := 0; t1^.bal := 0
                  END;
                  t := t1
                ELSE (* double RL *)
                  t2 := t1^.left; b2 := t2^.bal;
                  t1^.left := t2^.right; t2^.right := t1;
                  t^.right := t2^.left; t2^.left := t;
                  IF b2=+1 THEN t^.bal := -1 ELSE t^.bal := 0 END;
                  IF b2=-1 THEN t1^.bal := +1 ELSE t1^.bal := 0 END;
                  t := t2; t2^.bal := 0
                END;
        END 
      END BalanceLeft;

      PROCEDURE BalanceRight(VAR t: AVLRef; VAR h: BOOLEAN) =
        VAR 
          t1,t2: AVLRef;
          b1,b2: BalFactor;
      BEGIN
        CASE t^.bal OF
           1 => t^.bal := 0;
        |  0 => t^.bal := -1; h := FALSE;
        | -1 => (* rebalance *)
                t1 := t^.left; b1 := t1^.bal;
                IF b1<=0 THEN (* single LL *)
                  t^.left := t1^.right; t1^.right := t;
                  IF b1=0 THEN
                    t^.bal := -1; t1^.bal := +1; h := FALSE 
                  ELSE 
                    t^.bal := 0; t1^.bal := 0
                  END;
                  t := t1
                ELSE (* double LR *)
                  t2 := t1^.right; b2 := t2^.bal;
                  t1^.right := t2^.left; t2^.left := t1;
                  t^.left := t2^.right; t2^.right := t;
                  IF b2=-1 THEN t^.bal := +1 ELSE t^.bal := 0 END;
                  IF b2=+1 THEN t1^.bal := -1 ELSE t1^.bal := 0 END;
                  t := t2; t2^.bal := 0
                END;
        END 
      END BalanceRight;
      
      PROCEDURE Del(VAR r: AVLRef; VAR h: BOOLEAN) =
      BEGIN
        IF r^.right#NIL THEN
          Del(r^.right,h);
          IF h THEN BalanceRight(r,h) END;
        ELSE
          q^.info := r^.info; r := r^.left; h := TRUE 
        END
      END Del;

    BEGIN (* DoDelete *)
      IF t=NIL THEN h := FALSE; RETURN END;
      CASE CompText(w,t^.info) OF
        Result.Smaller =>
                DoDelete(t^.left,h);
                IF h THEN BalanceLeft(t,h) END
      | Result.Greater =>
                DoDelete(t^.right,h);
                IF h THEN
                  BalanceRight(t,h)
                END
      | Result.Equal =>  (* delete *)
                q := t;
                IF q^.right=NIL THEN
                  t :=q^.left; h := TRUE 
                ELSIF q^.left=NIL THEN
                  t := q^.right; h := TRUE 
                ELSE
                  Del(q^.left,h);
                  IF h THEN BalanceLeft(t,h) END 
                END
      END 
    END DoDelete;
    
  BEGIN (* Delete *)
    DoDelete(tr.Tree,DummyBool)
  END Delete;

  PROCEDURE Size(tr: T): CARDINAL =
  BEGIN
    RETURN tr.Nodes
  END Size;
  
  
  PROCEDURE Search(tr: T; READONLY w: TEXT): BOOLEAN =
    PROCEDURE DoSearch(t: AVLRef): BOOLEAN =
    BEGIN
      IF t=NIL  THEN
        RETURN FALSE
      ELSE
        CASE CompText(w,t^.info)  OF
          Result.Smaller   =>  RETURN DoSearch(t^.left)
        | Result.Equal     =>  RETURN TRUE
        | Result.Greater   =>  RETURN DoSearch(t^.right)
        END
      END
    END DoSearch;
  BEGIN
    RETURN DoSearch(tr.Tree)
  END Search;

  PROCEDURE New(): T =
  BEGIN
    RETURN NEW(T,Tree:=NIL,Nodes:=0)
  END New;


  PROCEDURE Enum(tr: T; 
                            Action: TextActionProc;
                            reverse: BOOLEAN := FALSE) RAISES {Abort} =
    PROCEDURE DoEnum(t: AVLRef) RAISES {Abort} =
    BEGIN
      TRY
        IF t#NIL THEN
          CASE reverse OF
            FALSE  => DoEnum(t^.left);
                      Action(t^.info);
                      DoEnum(t^.right)
          | TRUE   => DoEnum(t^.right);
                      Action(t^.info);
                      DoEnum(t^.left)
          END
        END
      EXCEPT
        Abort => RAISE Abort
      END
    END DoEnum;
  BEGIN
    TRY
      DoEnum(tr.Tree)
    EXCEPT
      Abort => RAISE Abort
    END
  END Enum;
   
  PROCEDURE AVLMakeTransition(<* UNUSED *> tr: T; 
                              s: AVLTextTreeState; 
                              ch: CHAR): AVLTextTreeState =
    VAR 
      t := s.Tree;
      p := s.Position;
      tinf := t^.info;

    PROCEDURE DoMake(lrt:AVLRef): AVLTextTreeState =
    BEGIN
      IF lrt=NIL THEN RETURN AVLTextTreeNullState END;
      WITH
        lrtinf=lrt^.info,
        lrtl=Text.Length(lrtinf)
      DO
        FOR i := 0 TO MIN(p-1,lrtl-1) DO
          WITH
            tch=Text.GetChar(tinf,i),
            lch=Text.GetChar(lrtinf,i)
          DO
            IF tch<lch THEN
              RETURN DoMake(lrt^.left)
            ELSIF tch>lch THEN
              RETURN DoMake(lrt^.right)
            END
          END
        END;
        IF p<lrtl THEN
          WITH
            nch=Text.GetChar(lrtinf,p)
          DO
            IF ch=nch THEN
              RETURN NEW(AVLTextTreeState,Tree:=lrt,Position:=p+1)
            ELSIF ch<nch THEN
              RETURN DoMake(lrt^.left)
            ELSE
              RETURN DoMake(lrt^.right)
            END
          END
        ELSE
          RETURN DoMake(lrt^.left)
        END
      END
    END DoMake;
    
  BEGIN (* AVLMakeTransition*)
    RETURN DoMake(t)
  END AVLMakeTransition;

  
  PROCEDURE AVLIsFinal(<* UNUSED *> tr: T; s:AVLTextTreeState): BOOLEAN =
    VAR
      t := s.Tree;
      p := s.Position;
      tinf := t^.info;
      
    PROCEDURE DoFinal(lrt: AVLRef): BOOLEAN =
    BEGIN
      IF lrt=NIL THEN RETURN FALSE END;
      WITH
        lrtinf=lrt^.info,
        lrtl=Text.Length(lrtinf)
      DO
        FOR i := 0 TO MIN(p-1,lrtl-1) DO
          WITH
            tch=Text.GetChar(tinf,i),
            lch=Text.GetChar(lrtinf,i)
          DO
            IF tch<lch THEN
              RETURN DoFinal(lrt^.left)
            ELSIF tch>lch THEN
              RETURN DoFinal(lrt^.right)
            END
          END
        END;
        IF p=lrtl THEN 
          RETURN TRUE 
        ELSIF p<lrtl THEN
          RETURN DoFinal(lrt^.left)
        ELSE
          RETURN DoFinal(lrt^.right)
        END
      END
    END DoFinal;
  BEGIN (* AVLIsFinal *)
    RETURN DoFinal(t)
  END AVLIsFinal;
      
  PROCEDURE RootState(tr: T): AVLTextTreeState =
  BEGIN
    RETURN NEW(AVLTextTreeState,Tree:=tr.Tree,Position:=0)
  END RootState;
   
  PROCEDURE IsNullState(s: AVLTextTreeState): BOOLEAN =
  BEGIN
    RETURN s.Tree=NIL
  END IsNullState;
   
  PROCEDURE CompText(READONLY x,y: TEXT): Result =
    VAR ind : INTEGER := 0;
  BEGIN
    WITH
      lx=Text.Length(x),
      ly=Text.Length(y),
      lxy=MIN(lx,ly)-1
    DO
      LOOP
        IF ind>lxy THEN 
          IF lx<ly THEN
            RETURN Result.Smaller
          ELSIF lx=ly  THEN
            RETURN Result.Equal
          ELSE
            RETURN Result.Greater
          END
        ELSE
          WITH
            chx=Text.GetChar(x,ind),
            chy=Text.GetChar(y,ind)
          DO
            IF chx<chy THEN
              RETURN Result.Smaller
            ELSIF  chx>chy THEN
              RETURN Result.Greater
            ELSE
              INC(ind)
            END
          END
        END
      END
    END
  END CompText;
  
BEGIN
(* Initialize AVL null state *)

AVLTextTreeNullState := NEW(AVLTextTreeState,Tree:=NIL,Position:=0)

END AVLTextTree.
