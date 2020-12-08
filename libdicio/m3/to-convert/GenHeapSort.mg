GENERIC MODULE GenHeapSort(X, Y);

(* See the copyright and disclaimer note at the end of this file. *)

PROCEDURE Sort(
    VAR (*IO*) s: ARRAY OF X.T; 
    VAR (*IO*) ss: REF ARRAY OF Y.T;
  ) = 
  VAR  q, r, mx: CARDINAL;
       cc,smx: X.T;
       ccss: Y.T;
  BEGIN

    WITH  
      n = NUMBER(s),
      tss = ss#NIL
    DO
      (* 1. Build heap with largest element at s[0] *)
      FOR  p := 0  TO n-1  DO
        cc := s[p]; 
        IF tss THEN ccss := ss[p] END;
        r := p;
        LOOP
          IF  r=0  THEN  EXIT  END;
          q := (r-1) DIV 2;
          IF X.Compare(s[q],cc) = -1  THEN
            s[r] := s[q]; 
            IF tss THEN ss[r] := ss[q] END;
          ELSE
            EXIT
          END;
          r := q
        END;
        s[r] := cc; 
        IF tss THEN ss[r] := ccss END
      END;
      (* 2. Remove elements from heap and insert at end *)
      FOR  p := n-1  TO 1  BY -1  DO
        (* save s[p] *)
        cc := s[p]; 
        IF tss THEN ccss := ss[p] END;
        (* Move largest heap element to pos[p] *)
        s[p] := s[0]; 
        IF tss THEN ss[p] := ss[0] END;
        (* Insert cc in remaining heap, from root down  *)
        q := 0;
        LOOP
          s[q] := cc; 
          IF tss THEN ss[q] := ccss END;
          r := 2*q+1;
          (* Find largest among cc, s[LEFT(q)], s[RIGHT(q)]  *)
          mx := q; smx := cc;
          IF (r<p) AND (X.Compare(s[r],smx) = +1)  THEN
            mx := r; smx := s[r]
          END;
          INC(r);
          IF  (r<p) AND (X.Compare(s[r],smx) = +1)  THEN
            mx := r; smx := s[r]
          END;
          (*  See who won  *)
          IF  mx=q  THEN
            (*  Stop here  *)
            EXIT
          ELSE
            (*  Promote child and advance  *)
            s[q] := s[mx];  
            IF tss THEN ss[q] := ss[mx] END; 
            q := mx
          END
        END
      END
    END
  END Sort;

BEGIN
END GenHeapSort.

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
 
