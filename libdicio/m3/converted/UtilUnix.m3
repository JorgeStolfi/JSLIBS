(* Last edited on 1999-06-05 20:50:38 by stolfi *) 

UNSAFE MODULE UtilUnix EXPORTS Util;

(* See the copyright and disclaimer note at the end of this file. *)

IMPORT M3toC, Unix, Uugid, Upwd, Ugrp, Fmt;

PROCEDURE GetUserName(effective: BOOLEAN := TRUE): TEXT =
  VAR uid: INTEGER;
  BEGIN
    IF effective THEN 
      uid := Uugid.geteuid()
    ELSE
      uid := Uugid.getuid()
    END;
    WITH user = (Upwd.getpwuid(uid)).pw_name DO
      IF user = NIL THEN
        RETURN Fmt.Int(uid)
      ELSE
        RETURN M3toC.CopyStoT(user)
      END
    END
  END GetUserName;

PROCEDURE GetGroupName(effective: BOOLEAN := TRUE): TEXT =
  VAR gid: INTEGER;
  BEGIN
    IF effective THEN 
      gid := Uugid.getegid()
    ELSE
      gid := Uugid.getgid()
    END;
    WITH group = (Ugrp.getgrgid(gid)).gr_name DO
      IF group = NIL THEN
        RETURN Fmt.Int(gid)
      ELSE
        RETURN M3toC.CopyStoT(group)
      END
    END
  END GetGroupName;

PROCEDURE GetHostName(): TEXT = 
  CONST MaxHostNameChars = 256;
  VAR buf: ARRAY [0..MaxHostNameChars+1] OF CHAR;
  BEGIN
    WITH res = Unix.gethostname(ADR(buf[0]), NUMBER(buf)-1) DO
      IF res # 0 THEN RETURN "" END;
      buf[LAST(buf)] := '\000';
      RETURN M3toC.CopyStoT(buf)
    END;
  END GetHostName;

BEGIN
END UtilUnix.

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
