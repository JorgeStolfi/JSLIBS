/* Procedures of {raut.h} that do not require access to the internal rep. */
/* Last edited on 2009-10-29 21:21:58 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>

#include <affirm.h>
#include <bool.h>

#include <rdag.h>
#include <raut.h>

/* Only for {raut_state_debug}: */
#include <raut_io.h> 
#include <rdag_io.h> 

bool_t raut_state_accepts(raut_t *A, raut_state_t u, uint32_t nstr, rdag_symbol_t str[])
  {
    demand(u.ac <= 1, "invalid accept bit");
    demand(u.nd <= raut_node_max(A), "invalid state node");
    int i;
    for (i = 0; (i < nstr) && (u.nd != rdag_node_NULL); i++)
      { u = raut_arc_follow(A, u, str[i]); }
    if (i < nstr)
      { return FALSE; }
    else
      { assert(i == nstr);
        return u.ac;
      }
  }

raut_state_t raut_string_add(raut_t *A, raut_state_t u, uint32_t nstr, rdag_symbol_t str[])
  {
    demand(u.ac <= 1, "invalid accept bit");
    demand(u.nd <= raut_node_max(A), "invalid state node");
    
    bool_t debug = FALSE;
    
    auto raut_state_t do_add(raut_state_t v, uint32_t k);
      /* Returns a state {v'} such that {Suff(v')} is {Suff(v)} plus
         the string {str[k..nstr-1]}. */
         
    raut_state_t do_add(raut_state_t v, uint32_t k)
      {
        if (debug) 
          { fprintf(stderr, "  %*sadding ", 2*k, "");
            rdag_string_iso_latin_1_write(stderr, nstr-k, str);
            fprintf(stderr, " to state (%u,%u)\n", v.ac, v.nd);
          }
        raut_state_t v1; /* Resulting state. */
        if (k >= nstr)
          { /* Make the state {v} accept the empty string too: */
            v1 = (raut_state_t){ .ac = 1, .nd = v.nd };
          }
        else
          { /* Get the state {w0} that accepts {str[k] \rquot Suff(v)}: */
            raut_state_t w0 = raut_arc_follow(A, v, str[k]);
            /* Get the state that accepts {Suff(w0) \union str[k+1..nstr-1]}: */
            raut_state_t w1 = do_add(w0, k+1);
            if ((w1.ac == w0.ac) && (w1.nd == w0.nd))
              { /* Alreay in there, no change: */ return v; }
            /* Now replace {str[k] \concat Suff(w0)} by {str[k] \concat Suff(w1)}: */
            if (debug) 
              { fprintf(stderr, "  %*sset arc (%u,%u) --[%u]--> (%u,%u)\n", 2*k, "", v.ac, v.nd, str[k], w1.ac, w1.nd); }
            v1 = raut_arc_set(A, v, str[k], w1);
          }
        if (debug) 
          { fprintf(stderr, "  %*sresult is state (%u,%u)", 2*k, "", v1.ac, v1.nd);
            raut_state_debug(stderr, " = ", A, v1, "\n");
          }
        return v1;
      }
      
    return do_add(u, 0);
  }

raut_state_t raut_string_remove(raut_t *A, raut_state_t u, uint32_t nstr, rdag_symbol_t str[])
  {
    demand(u.ac <= 1, "invalid accept bit");
    demand(u.nd <= raut_node_max(A), "invalid state node");
    
    auto raut_state_t do_remove(raut_state_t v, uint32_t k);
      /* Returns a state {v'} such that {Suff(v')} is {Suff(v)} minus
         the string {str[k..nstr-1]}. */
         
    raut_state_t do_remove(raut_state_t v, uint32_t k)
      {
        if (k >= nstr)
          { /* Make the state {v} reject the empty string: */
            return (raut_state_t){ .ac = 0, .nd = v.nd };
          }
        else
          { /* Get the state {w0} that accepts {str[k] \rquot Suff(v)}: */
            raut_state_t w0 = raut_arc_follow(A, v, str[k]);
            /* Get the state that accepts {Suff(w0) \setminus str[k+1..nstr-1]}: */
            raut_state_t w1 = do_remove(w0, k+1);
            if ((w1.ac == w0.ac) && (w1.nd == w0.nd))
              { /* Wasn't there, no change: */ return v; }
            /* Replace {str[k] \concat Suff(w0)} by {str[k] \concat Suff(w1)} in {Suff(v)}: */
            raut_state_t v1 = raut_arc_set(A, v, str[k], w1);
            return v1;
          }
      }
      
    return do_remove(u, 0);
  }

//CONST
//    NullLetter = FIRST(rdag_basics.h.rdag_symbol_t);
//    (*
//      A transition from "s" through NullLetter to raut_state_NULL is used
//      to indicate that "s" is final. *)
//  
//  PROCEDURE New(size: POS := 1): T =
//    BEGIN
//      WITH
//        dag = DAG.New(size)
//      DO
//        RETURN FromDAG(dag, doc := "", root := raut_state_NULL)
//      END
//    END
//  
//  PROCEDURE Copy(
//      aut: T;
//      size: POS := 1;
//      VAR (*OUT*) map: REF ARRAY OF raut_state_t;
//    ): T RAISES {raut_error_FULL} =
//    BEGIN
//      WITH
//        newAut = New(size),
//        newRoot = raut_copy_states(
//          from := aut,
//          to := newAut,
//          s := aut.raut_root(),
//          map := (*OUT*) map
//        )
//      DO
//        newAut.raut_set_root(newRoot);
//        RETURN newAut
//      END
//    END
//  
//  PROCEDURE raut_copy_states(
//      from: T;
//      to: T;
//      s: raut_state_t;
//      VAR (*IO*) map: REF ARRAY OF raut_state_t;
//    ): raut_state_t RAISES {raut_error_FULL} =
//    BEGIN
//      RETURN DAG.Copy(from.dag, to.dag, s, map)
//    END
//  
//  (************************)
//  (* PRIMITIVE METHODS    *)
//  (************************)
//  
//  PROCEDURE raut_has_arcs(<*UNUSED*> aut: T; s: raut_state_t): BOOL =
//    BEGIN
//      RETURN s # raut_state_NULL AND s # raut_state_UNIT
//    END
//  
//  PROCEDURE rdag_last(aut: T; s: raut_state_t): Arc =
//    BEGIN
//      <* ASSERT s # raut_state_NULL AND s # raut_state_UNIT *>
//      WITH a = aut.dag.rdag_last(s) DO
//        RETURN Arc{symbol := a.rd, dest := a.dest}
//      END
//    END
//  
//  PROCEDURE rdag_rest(aut: T; s: raut_state_t): raut_state_t =
//    BEGIN
//      <* ASSERT s # raut_state_NULL AND s # raut_state_UNIT *>
//      RETURN aut.dag.rdag_rest(s)
//    END
//  
//  PROCEDURE Append(aut:T; s: raut_state_t; symbol: rdag_symbol_t; dest: raut_state_t): raut_state_t RAISES {raut_error_FULL} =
//    BEGIN
//      IF s # raut_state_NULL THEN
//        <* ASSERT symbol > aut.dag.rdag_last(s).rd *>
//      END;
//      IF dest = raut_state_NULL THEN
//        RETURN s
//      ELSE
//        WITH
//          oldMax = aut.raut_max_state(),
//          t = aut.dag.Append(DAG.Arc{rd := symbol, wr := 0, dest := dest}, s)
//        DO
//          IF t > oldMax THEN
//            (* raut_state_t is truly new, must compute aut.nSuffs[t], aut.nSuffLetters[t] *)
//            WITH ns = aut.nSuffs^, nl = aut.nSuffLetters^ DO
//              ns[t] := ns[dest] + ns[s];
//              nl[t] := nl[dest] + ns[dest] + nl[s];
//            END;
//          END;
//          RETURN t
//        END
//      END
//    END
//  
//  PROCEDURE raut_root(aut: T): raut_state_t =
//    BEGIN
//      RETURN aut.root
//    END
//  
//  PROCEDURE raut_set_root(aut: T; s: raut_state_t) =
//    BEGIN
//      aut.root := s
//    END
//  
//  PROCEDURE raut_max_state(aut: T): POS =
//    BEGIN
//      RETURN aut.dag.raut_max_state()
//    END
//  
//  PROCEDURE raut_max_alloc_state(aut: T): POS =
//    BEGIN
//      RETURN aut.dag.raut_max_alloc_state()
//    END
//  
//  PROCEDURE raut_expand(aut: T; newSize: uint32_t := 0) RAISES {raut_error_FULL} =
//    BEGIN
//      aut.dag.raut_expand(newSize);
//      ComputeNSuffs(aut);
//    END
//  
//  PROCEDURE raut_crunch(aut: T; keep: REF ARRAY OF raut_state_t := NIL) =
//    VAR nKeep: uint32_t;
//    BEGIN
//      DiscardPrefixData(aut);
//      IF keep = NIL THEN nKeep := 0 ELSE nKeep := NUMBER(keep^) END;
//      WITH
//        roots = NEW(REF ARRAY OF raut_state_t, nKeep+1)^
//      DO
//        roots[0] := aut.root;
//        IF keep # NIL THEN SUBARRAY(roots, 1, nKeep) := keep^ END;
//        aut.dag.raut_crunch(roots);
//        aut.root := roots[0];
//        IF keep # NIL THEN keep^ := SUBARRAY(roots, 1, nKeep) END;
//      END;
//      (* Recreate unit state, which may have been crunched out: *)
//      MakeUnitState(aut.dag);
//      (* Recompute suffix counts: *)
//      ComputeNSuffs(aut);
//    END
//  
//  (******************************)
//  (* DERIVED METHODS            *)
//  (******************************)
//  
//  PROCEDURE raut_is_final(aut: T; s: raut_state_t): BOOL =
//    VAR t: raut_state_t := s;
//    BEGIN
//      WHILE t # raut_state_NULL AND t # raut_state_UNIT DO t := aut.rdag_rest(t) END;
//      RETURN t = raut_state_UNIT
//    END
//  
//  PROCEDURE raut_step(aut: T; s: raut_state_t; symbol: rdag_symbol_t): raut_state_t =
//    VAR t: raut_state_t := s;
//    BEGIN
//      IF s = raut_state_NULL THEN RETURN raut_state_NULL END;
//      WHILE t # raut_state_NULL AND t # raut_state_UNIT DO
//        WITH a = aut.rdag_last(t) DO
//          IF a.symbol = symbol THEN
//            RETURN a.dest
//          ELSIF a.symbol < symbol THEN
//            RETURN raut_state_NULL
//          ELSE
//            t := aut.rdag_rest(t)
//          END
//        END
//      END;
//      RETURN raut_state_NULL
//    END
//  
//  PROCEDURE rdag_out_deg(aut: T; s: raut_state_t): uint32_t =
//    VAR t: raut_state_t := s; n: uint32_t := 0;
//    BEGIN
//      WHILE t # raut_state_NULL AND t # raut_state_UNIT DO INC(n); t := aut.rdag_rest(t) END;
//      RETURN n
//    END
//  
//  PROCEDURE raut_in_deg(aut: T; s: raut_state_t): uint32_t =
//    BEGIN
//      IF aut.root = raut_state_NULL OR s >= aut.root THEN
//        RETURN 0
//      ELSE
//        ComputePrefixData(aut);
//        <* ASSERT aut.nPrefs # NIL *>
//        RETURN aut.rdag.rdag_out_deg(aut.rev[s])
//      END
//    END
//  
//  PROCEDURE rdag_first(aut: T; s: raut_state_t): Arc =
//    VAR p, t: raut_state_t;
//    BEGIN
//      <* ASSERT s # raut_state_NULL AND s # raut_state_UNIT *>
//      t := s;
//      REPEAT p := t; t := aut.rdag_rest(t) UNTIL t = raut_state_NULL OR t = raut_state_UNIT;
//      RETURN aut.rdag_last(p)
//    END
//  
//  PROCEDURE raut_set_arc(aut: T; s: raut_state_t; symbol: rdag_symbol_t; dest: raut_state_t): raut_state_t RAISES {raut_error_FULL} =
//  
//    PROCEDURE DoSet(t: raut_state_t): raut_state_t RAISES {raut_error_FULL} =
//      BEGIN
//        IF t = raut_state_NULL OR t = raut_state_UNIT THEN
//          RETURN aut.Append(t, symbol, dest)
//        ELSE
//          WITH a = aut.rdag_last(t) DO
//            IF a.symbol < symbol THEN
//              RETURN aut.Append(t, symbol, dest)
//            ELSIF a.symbol = symbol THEN
//              IF a.dest = dest THEN
//                RETURN t
//              ELSE
//                RETURN aut.Append(aut.rdag_rest(t), symbol, dest)
//              END;
//            ELSE (* a.symbol > symbol *)
//              WITH
//                rest = aut.rdag_rest(t),
//                newRest = DoSet(rest)
//              DO
//                IF rest = newRest THEN
//                  RETURN t
//                ELSE
//                  RETURN aut.Append(newRest, a.symbol, a.dest)
//                END
//              END
//            END
//          END
//        END
//      END
//  
//    BEGIN
//      RETURN DoSet(s)
//    END
//  
//  PROCEDURE raut_set_final(aut: T; s: raut_state_t; final: BOOL): raut_state_t RAISES {raut_error_FULL} =
//  
//    PROCEDURE DoSet(t: raut_state_t): raut_state_t RAISES {raut_error_FULL} =
//      BEGIN
//        IF t = raut_state_NULL OR t = raut_state_UNIT THEN
//          IF final THEN 
//            RETURN raut_state_UNIT 
//          ELSE
//            RETURN raut_state_NULL
//          END
//        ELSE
//          WITH 
//            a = aut.rdag_last(t),
//            rest = aut.rdag_rest(t),
//            newRest = DoSet(rest)
//          DO
//            IF newRest = rest THEN 
//              RETURN t
//            ELSE
//              RETURN aut.Append(newRest, a.symbol, a.dest)
//            END
//          END
//        END
//      END
//  
//    BEGIN
//      RETURN DoSet(s)
//    END
//  
//  PROCEDURE raut_walk(aut: T; s: raut_state_t; READONLY w: String): raut_state_t =
//    VAR t: raut_state_t := s;
//    BEGIN
//      FOR i := 0 TO LAST(w) DO
//        <* ASSERT w[i] # NullLetter *>
//        IF t = raut_state_NULL THEN
//          RETURN raut_state_NULL
//        ELSE
//          t := aut.raut_step(t, w[i])
//        END
//      END;
//      RETURN t
//    END
//  
//  PROCEDURE raut_accepts(aut: T; s: raut_state_t; READONLY w: String): BOOL =
//    BEGIN
//      RETURN aut.raut_is_final(aut.raut_walk(s, w))
//    END
//  
//  PROCEDURE raut_rank(aut: T; s: raut_state_t; READONLY w: String; reverse: BOOL := FALSE): uint32_t =
//    BEGIN
//      IF reverse THEN
//        RETURN RankFromBottom(aut, s, w)
//      ELSE
//        RETURN RankFromTop(aut, s, w)
//      END
//    END
//  
//  PROCEDURE RankFromTop(aut: T; s: raut_state_t; READONLY w: String): uint32_t =
//  (*
//    Computes raut_rank(aut, s, w, reverse := FALSE).
//    *)
//    VAR rank: uint32_t := 0;
//        i: uint32_t := 0;
//        t: raut_state_t := s;
//    BEGIN
//      WHILE t # raut_state_NULL AND i <= LAST(w) DO
//        (* Tally the empty string, if it is a suffix of "t": *)
//        IF aut.raut_is_final(t) THEN INC(rank) END;
//        WITH x = w[i] DO
//          (* Tally all suffixes of "t" starting with symbols less than "x": *)
//          <* ASSERT x # NullLetter *>
//          IF x > FIRST(rdag_symbol_t) THEN
//            rank := rank + AddNSuffsOfChildren(aut, t, lo := FIRST(rdag_symbol_t), hi := x-1)
//          END;
//          (* Now add the rank of "w[i+1..]" in raut_suff(raut_step(t, x)): *)
//          t := aut.raut_step(t, x);
//          INC(i);
//        END
//      END;
//      RETURN rank
//    END
//  
//  PROCEDURE RankFromBottom(aut: T; s: raut_state_t; READONLY w: String): uint32_t =
//  (*
//    Computes raut_rank(aut, s, w, reverse := TRUE).
//    *)
//    VAR rank: uint32_t := 0;
//        i: uint32_t := 0;
//        t: raut_state_t := s;
//    BEGIN
//      WHILE t # raut_state_NULL AND i <= LAST(w) DO
//        WITH x = w[i] DO
//          (* Tally all suffixes of "t" starting with symbols greater than "x": *)
//          <* ASSERT x # NullLetter *>
//          IF x < LAST(rdag_symbol_t) THEN
//            rank := rank + AddNSuffsOfChildren(aut, t, lo := x+1, hi := LAST(rdag_symbol_t))
//          END;
//          (* Now add the reverse rank of "w[i+1..]" in raut_suff(raut_step(t, x)): *)
//          t := aut.raut_step(t, x);
//          INC(i);
//        END
//      END;
//      IF t # raut_state_NULL THEN
//        (* Tally all suffixes of "t" that are strictly grater than (): *)
//        rank := rank + aut.raut_num_suffs(t);
//        IF aut.raut_is_final(t) THEN DEC(rank) END
//      END;
//      RETURN rank
//    END
//  
//  PROCEDURE AddNSuffsOfChildren(aut: T; s: raut_state_t; lo, hi: rdag_symbol_t): uint32_t =
//  (*
//    Returns the sum of raut_num_suffs(raut_step(s, x)) for all "x" in the range [lo..hi]. *)
//    VAR sum: uint32_t := 0;
//        t: raut_state_t := s;
//    BEGIN
//      WHILE aut.raut_has_arcs(t) DO
//        WITH a = aut.rdag_last(t) DO
//          IF a.symbol < lo THEN
//            RETURN sum
//          ELSIF a.symbol <= hi THEN
//            INC(sum, aut.raut_num_suffs(a.dest))
//          END;
//        END;
//        t := aut.rdag_rest(t)
//      END;
//      RETURN sum
//    END
//  
//  PROCEDURE raut_add_string(aut: T; s: raut_state_t; READONLY w: String): raut_state_t RAISES {raut_error_FULL} =
//    BEGIN
//      RETURN AddSubString(aut, s, w, TRUE)
//    END
//  
//  PROCEDURE raut_sub-string(aut: T; s: raut_state_t; READONLY w: String): raut_state_t RAISES {raut_error_FULL} =
//    BEGIN
//      RETURN AddSubString(aut, s, w, FALSE)
//    END raut_sub-string;
//  
//  PROCEDURE AddSubString(
//      aut: T;
//      s: raut_state_t;
//      READONLY w: String;
//      add: BOOL;
//    ): raut_state_t RAISES {raut_error_FULL} =
//  (*
//    Does raut_add_string(aut, s, w) or raut_sub-string(aut, s, w),
//    depending on "add". *)
//  
//    PROCEDURE DoAddSub(i: uint32_t; t: raut_state_t): raut_state_t RAISES {raut_error_FULL} =
//    (*
//      Returns AddSubString(aut, t, SUBARRAY(w, i), add). *)
//      BEGIN
//        IF i = NUMBER(w) THEN
//          RETURN aut.raut_set_final(t, add)
//        ELSE
//          <* ASSERT w[i] # NullLetter *>
//          RETURN aut.raut_set_arc(t, w[i], DoAddSub(i+1, aut.raut_step(t, w[i])))
//        END
//      END
//  
//    BEGIN
//      RETURN DoAddSub(0, s)
//    END
//  
//  PROCEDURE raut_enum_out_arcs(aut: T; s: raut_state_t; action: rdag_arc_action_t) RAISES {rdag_disp_STOP} =
//    VAR i: uint32_t := 0;
//  
//    PROCEDURE DoEnum(r: raut_state_t) RAISES {rdag_disp_SKIP, rdag_disp_STOP} =
//    (*
//      Does raut_enum_out_arcs on a given prefix "r" of the arcs out of "s".
//      Raises rdag_disp_SKIP, rdag_disp_STOP iff "action" does so. *)
//      BEGIN
//        IF r = raut_state_UNIT OR r = raut_state_NULL THEN
//          (* Ok *)
//        ELSE
//          DoEnum(aut.rdag_rest(r));
//          WITH a = aut.rdag_last(r) DO
//            action(i, a);
//            INC(i)
//          END
//        END
//      END
//  
//    BEGIN
//      TRY DoEnum(s) EXCEPT rdag_disp_SKIP => (* Ok *) END;
//    END
//  
//  PROCEDURE raut_add_string_maybe_crunch(
//      aut: T; 
//      VAR s: raut_state_t;
//      READONLY w: String; 
//      keep: REF States := NIL;
//    ): raut_state_t =
//    BEGIN
//      RETURN AddSubStringMaybeCrunch(aut, s, w, TRUE, keep)
//    END
//  
//  PROCEDURE raut_sub_string_maybe_crunch(
//      aut: T; 
//      VAR s: raut_state_t; 
//      READONLY w: String; 
//      keep: REF States := NIL;
//    ): raut_state_t =
//    BEGIN
//      RETURN AddSubStringMaybeCrunch(aut, s, w, FALSE, keep)
//    END
//  
//  PROCEDURE AddSubStringMaybeCrunch(
//      aut: T; 
//      VAR s: raut_state_t; 
//      READONLY w: String; 
//      add: BOOL;
//      keep: REF States := NIL;
//    ): raut_state_t =
//  
//    PROCEDURE CrunchAndExpandIt() =
//    (*
//      Does aut.raut_crunch(s & keep), and perhaps aut.raut_expand(). *)
//      VAR nKeep: CARDINAL := 0;
//      BEGIN
//        IF keep # NIL THEN nKeep := NUMBER(keep^) END;
//        WITH 
//          tKeep = NEW(REF States, nKeep + 1) 
//        DO
//          tKeep[0] := s;
//          IF keep # NIL THEN SUBARRAY(tKeep^, 1, nKeep) := keep^ END;
//          aut.raut_crunch(tKeep);
//          IF keep # NIL THEN keep^ := SUBARRAY(tKeep^, 1, nKeep) END;
//          s := tKeep[0];
//        END;
//        IF aut.raut_max_state() * 4 >= aut.raut_max_alloc_state() * 3 THEN 
//          ExpandIt()
//        END;
//      END
//  
//    PROCEDURE ExpandIt() =
//    (*
//      Expands the automaton: *)
//      BEGIN
//        WITH
//          oldSize = aut.raut_max_alloc_state() + 1,
//          newSize = ComputeNewSize(oldSize)
//        DO
//          <* FATAL rdag_basics.h.raut_error_FULL *>
//          BEGIN
//            aut.raut_expand(newSize);
//          END;
//        END;
//      END
//  
//    PROCEDURE ComputeNewSize(oldSize: uint32_t): uint32_t =
//    (*
//      Chooses a suitable new size for automatic expansion.
//      *)
//      BEGIN
//        RETURN
//          MIN(
//            MIN(
//                  1 + oldSize + oldSize,
//              10001 + oldSize + oldSize DIV 2
//            ),
//            MIN(
//              20001 + oldSize + oldSize DIV 3,
//              50001 + oldSize + oldSize DIV 4
//            )
//          )
//      END
//  
//    BEGIN
//      LOOP
//        TRY
//          RETURN AddSubString(aut, s, w, add);
//        EXCEPT 
//          raut_error_FULL => CrunchAndExpandIt();
//        END;
//      END;
//    END
//  
//  PROCEDURE raut_enum_in_arcs(aut: T; s: raut_state_t; action: rdag_arc_action_t) RAISES {rdag_disp_STOP} =
//    VAR i: uint32_t := 0;
//        rr: DAG.raut_state_t;
//    BEGIN
//      IF s = raut_state_NULL AND s >= aut.root THEN
//        RETURN
//      ELSE
//        ComputePrefixData(aut);
//        rr := aut.rev[s];
//        i := 0;
//        TRY
//          WITH rdag = aut.rdag, dir = aut.dir^ DO
//            WHILE rr # DAG.raut_state_NULL DO
//              WITH a = rdag.rdag_last(rr) DO
//                action(i, Arc{symbol := a.rd, dest := dir[a.dest]});
//                INC(i);
//                rr := rdag.rdag_rest(rr)
//              END;
//            END;
//          END
//        EXCEPT
//        | rdag_disp_SKIP => (* Ok *)
//        END;
//      END;
//    END
//  
//  PROCEDURE rdag_enum_paths(
//      aut: T;
//      s: raut_state_t;
//      enter: raut_state_action_t := NIL;
//      push: raut_path_action_t := NIL;
//      pop: raut_path_action_t := NIL;
//      exit: raut_state_action_t := NIL;
//    ) RAISES {rdag_disp_STOP} =
//  
//    VAR len: uint32_t := 0;
//  
//    PROCEDURE DoEnumPaths(t: raut_state_t) RAISES {rdag_disp_STOP} =
//    (*
//      Does rdag_enum_paths starting at "t", a generic sucessor
//      of "s". Assumes "t" is not raut_state_NULL. *)
//  
//      VAR final: BOOL; (* raut_is_final(t); set by "EnumRest" at the end of the recursion. *)
//          i: uint32_t := 0; (* Arc index (not including the NullLetter, if any) *)
//  
//      PROCEDURE EnumRest(r: raut_state_t) RAISES {rdag_disp_SKIP, rdag_disp_STOP} =
//      (*
//        Calls "enter" on "t" and enumerates a given prefix "r" of the arcs out of "t".
//        Also sets the local variable "final" of DoEnumStates(t)) to raut_is_final(t).
//        Raises rdag_disp_SKIP iff "enter" or "pop" raised "rdag_disp_SKIP". *)
//  
//        BEGIN
//          IF r = raut_state_UNIT OR r = raut_state_NULL THEN
//            final := (r = raut_state_UNIT);
//            IF enter # NIL THEN enter(len, t, final) END;
//          ELSE
//            EnumRest(aut.rdag_rest(r));
//            WITH a = aut.rdag_last(r) DO
//              TRY
//                IF push # NIL THEN push(len, t, i, a) END;
//                INC(len);
//                DoEnumPaths(a.dest);
//                DEC(len);
//              EXCEPT
//                rdag_disp_SKIP => (* Ok *)
//              END;
//              IF pop # NIL THEN pop(len, t, i, a) END;
//              INC(i);
//            END
//          END
//        END
//  
//      BEGIN
//        <* ASSERT t # raut_state_NULL *>
//        TRY EnumRest(t) EXCEPT rdag_disp_SKIP => (* Ok *) END;
//        IF exit # NIL THEN
//          <* FATAL rdag_disp_SKIP *>
//          BEGIN
//            exit(len, t, final) (* Shouldn't raise "rdag_disp_SKIP" *)
//          END
//        END
//      END
//  
//    BEGIN
//      IF s # raut_state_NULL THEN DoEnumPaths(s) END;
//    END
//  
//  PROCEDURE raut_enum_strings(
//      aut: T;
//      s: raut_state_t;
//      enter: raut_string-action_t := NIL;
//      exit: raut_string-action_t := NIL;
//    ) RAISES {rdag_disp_STOP} =
//  
//    VAR rw: REF String := NEW(REF String, 100);
//  
//    PROCEDURE PathEnter (* : raut_state_action_t *) (
//        len: uint32_t;
//        s: raut_state_t;
//        final: BOOL;
//      ) RAISES {rdag_disp_SKIP, rdag_disp_STOP} =
//      BEGIN
//        IF enter # NIL THEN enter(SUBARRAY(rw^, 0, len), s, final) END
//      END
//  
//    PROCEDURE PathPush (* : raut_path_action_t *) (
//        len: uint32_t;
//        <*UNUSED*> org: raut_state_t;
//        <*UNUSED*> i: uint32_t;
//        arc: Arc;
//      ) RAISES {} =
//      BEGIN
//        rdag_basics.h.ExpandString(rw, len + 1);
//        rw[len] := arc.symbol;
//      END
//  
//    PROCEDURE PathPop (* : raut_path_action_t *) (
//        len: uint32_t;
//        <*UNUSED*> org: raut_state_t;
//        <*UNUSED*> i: uint32_t;
//        arc: Arc;
//      ) RAISES {} =
//      BEGIN
//        <* ASSERT rw[len] = arc.symbol *>
//      END
//  
//    PROCEDURE PathExit (* : raut_state_action_t *) (
//        len: uint32_t;
//        s: raut_state_t;
//        final: BOOL;
//      ) RAISES {rdag_disp_SKIP, rdag_disp_STOP} =
//      BEGIN
//        IF exit # NIL THEN exit(SUBARRAY(rw^, 0, len), s, final) END
//      END
//  
//    BEGIN
//      aut.rdag_enum_paths(
//        s,
//        enter := PathEnter,
//        push := PathPush,
//        pop := PathPop,
//        exit := PathExit
//      );
//    END
//  
//  PROCEDURE raut_enum_states(
//      aut: T;
//      READONLY base: ARRAY OF raut_state_t;
//      substates: BOOL := FALSE;
//      enter: raut_state_action_t := NIL;
//      exit: raut_state_action_t := NIL;
//    ) RAISES {rdag_disp_STOP} =
//  
//    VAR maxState: raut_state_t;
//    BEGIN
//      (* Computes maximum reachable state: *)
//      maxState := raut_state_NULL;
//      FOR i := 0 TO LAST(base) DO maxState := MAX(maxState, base[i]) END;
//  
//      WITH
//        len = NEW(REF ARRAY OF uint32_t, maxState + 1)^
//      DO
//  
//        (* Initialize path lengths: *)
//        FOR s := 0 TO maxState DO len[s] := LAST(uint32_t) END;
//        FOR i := 0 TO LAST(base) DO len[base[i]] := 0 END;
//  
//        (* rdag_first pass: scan states from highest to lowest,
//        propagating "len" and "enter"ing all reachable states: *)
//  
//        FOR s := maxState TO 1 BY -1 DO
//          WITH
//            ls = len[s],
//            final = aut.raut_is_final(s)
//          DO
//            IF ls < LAST(uint32_t) THEN
//              TRY
//                IF enter # NIL THEN enter(ls, s, final) END;
//                IF substates THEN
//                  IF s # raut_state_UNIT THEN
//                    WITH d = aut.rdag_last(s).dest, r = aut.rdag_rest(s) DO
//                      len[d] := MIN(len[d], ls + 1);
//                      len[r] := MIN(len[r], ls + 1)
//                    END
//                  END
//                ELSE
//                  VAR t := s;
//                  BEGIN
//                    WHILE t # raut_state_UNIT AND t # raut_state_NULL DO
//                      WITH d = aut.rdag_last(t).dest DO
//                        len[d] := MIN(len[d], ls + 1)
//                      END;
//                      t := aut.rdag_rest(t)
//                    END
//                  END
//                END;
//              EXCEPT
//              | rdag_disp_SKIP => (* Ignore any arcs out of "s" *)
//              END
//            END
//          END
//        END;
//  
//        (* Second pass: scan states from low to high, "exit"ing all
//        reachable ones: *)
//  
//        FOR s := 1 TO maxState DO
//          WITH
//            ls = len[s],
//            final = aut.raut_is_final(s)
//          DO
//            IF ls < LAST(uint32_t) THEN
//              TRY
//                IF exit # NIL THEN exit(ls, s, final) END;
//              EXCEPT
//              | rdag_disp_SKIP => <* ASSERT FALSE *>
//              END;
//            END
//          END
//        END;
//      END
//    END
//  
//  CONST
//    DumpHeader = "Reduced.rdag_dump (format of 91-12-21)";
//  
//  PROCEDURE rdag_dump(wr: FILE *; aut: T) =
//    <* FATAL Wr.Failure, Thread.Alerted *>
//    BEGIN
//      Wr.PutText(wr, "Begin " & DumpHeader); Wr.PutChar(wr, '\n');
//      DumpDoc(wr, aut.doc);
//      Wr.PutText(wr, "root = " & Fmt.Int(aut.root) & "\n");
//      DAG.rdag_dump(wr, aut.dag);
//      Wr.PutText(wr, "End " & DumpHeader); Wr.PutChar(wr, '\n');
//      Wr.Flush(wr);
//    END
//    
//  EXCEPTION  (* Syntax errors in dump file: *)
//    MissingFinalNewLine;
//    InvalidHeader;
//    InvalidFooter;
//  
//  PROCEDURE Load(rd: FILE *; minSize: POS := 1): T =
//    <* FATAL Rd.Failure, Rd.EndOfFile, Thread.Alerted *>
//    <* FATAL InvalidHeader, InvalidFooter *>
//    BEGIN
//      WITH hdr = Rd.GetLine(rd) DO
//        IF NOT Text.Equal(hdr, "Begin " & DumpHeader) THEN RAISE InvalidHeader END;
//      END;
//      WITH
//        doc = LoadDoc(rd),
//        root = ReadParam(rd, "root = "),
//        dag = DAG.Load(rd, minSize)
//      DO
//        WITH hdr = Rd.GetLine(rd) DO
//          IF NOT Text.Equal(hdr, "End " & DumpHeader) THEN RAISE InvalidFooter END
//        END;
//        RETURN FromDAG(dag, doc, root)
//      END;
//    END
//    
//  PROCEDURE Print(
//      wr: FILE *;
//      aut: T;
//      e: rdag_encoding.h.T;
//    ) =
//    CONST Finalrdag_code.h = ARRAY BOOL OF CHAR {' ', '*'};
//  
//    PROCEDURE PrintState (* : raut_state_action_t *) (
//        <*UNUSED*> len: uint32_t;
//        s: raut_state_t;
//        final: BOOL;
//      ) RAISES {} =
//      <* FATAL Wr.Failure, Thread.Alerted *>
//      BEGIN
//        WrNat(wr, s);
//        Wr.PutChar(wr, ' ');
//        Wr.PutChar(wr, Finalrdag_code.h[final]);
//        Wr.PutChar(wr, '\n');
//        <* FATAL rdag_disp_SKIP, rdag_disp_STOP *>
//        BEGIN
//          aut.raut_enum_out_arcs(s, PrintArc);
//        END
//      END
//  
//    PROCEDURE PrintArc (* : rdag_arc_action_t *) (
//        <*UNUSED*> i: uint32_t;
//        arc: Arc;
//      ) RAISES {} =
//      <* FATAL Wr.Failure, Thread.Alerted *>
//      BEGIN
//        Wr.PutChar(wr, ' ');
//        Wr.PutChar(wr, ' ');
//        e.PrintLetter(wr, arc.symbol);
//        Wr.PutChar(wr, ' ');
//        Wr.PutChar(wr, '-');
//        Wr.PutChar(wr, '>');
//        Wr.PutChar(wr, ' ');
//        WrNat(wr, arc.dest);
//        Wr.PutChar(wr, '\n');
//      END
//  
//    <* FATAL Wr.Failure, Thread.Alerted, rdag_disp_STOP *>
//    BEGIN
//      Wr.PutText(wr, "root = " & Fmt.Int(aut.root) & "\n");
//      Wr.PutChar(wr, '\n');
//      aut.raut_enum_states(base := ARRAY OF raut_state_t{aut.root}, enter := PrintState);
//      Wr.Flush(wr);
//    END
//  
//  PROCEDURE raut_num_states(aut: T; READONLY base: ARRAY OF raut_state_t): uint32_t =
//    VAR nStates: uint32_t := 0;
//  
//    PROCEDURE CountState (* : raut_state_action_t *) (
//        <*UNUSED*> len: uint32_t;
//        <*UNUSED*> s: raut_state_t;
//        <*UNUSED*> final: BOOL;
//      ) RAISES {} =
//      BEGIN
//        INC(nStates);
//      END
//  
//    <* FATAL rdag_disp_STOP *>
//    BEGIN
//      aut.raut_enum_states(base := base, enter := CountState);
//      RETURN nStates;
//    END
//  
//  PROCEDURE raut_num_arcs(aut: T; READONLY base: ARRAY OF raut_state_t): uint32_t =
//    VAR nArcs: uint32_t := 0;
//  
//    PROCEDURE CountArcs (* : raut_state_action_t *) (
//        <*UNUSED*> len: uint32_t;
//        s: raut_state_t;
//        <*UNUSED*> final: BOOL;
//      ) RAISES {} =
//      BEGIN
//        INC(nArcs, aut.rdag_out_deg(s));
//      END
//  
//    <* FATAL rdag_disp_STOP *>
//    BEGIN
//      aut.raut_enum_states(base := base, enter := CountArcs);
//      RETURN nArcs;
//    END
//  
//  PROCEDURE raut_count(aut: T; READONLY base: ARRAY OF raut_state_t): raut_counts_t =
//  
//    VAR nStates: uint32_t := 0;
//        nArcs: uint32_t := 0;
//        nFinals: uint32_t := 0;
//        nSubStates: uint32_t := 0;
//        nStrings: uint32_t := 0;
//        nLetters: uint32_t := 0;
//  
//    PROCEDURE CountState (* : raut_state_action_t *) (
//        <*UNUSED*> len: uint32_t;
//        s: raut_state_t;
//        final: BOOL;
//      ) RAISES {} =
//      BEGIN
//        INC(nStates);
//        INC(nArcs, aut.rdag_out_deg(s));
//        IF final THEN INC(nFinals) END;
//      END
//  
//    PROCEDURE CountSubState (* : raut_state_action_t *) (
//        <*UNUSED*> len: uint32_t;
//        <*UNUSED*> s: raut_state_t;
//        <*UNUSED*> final: BOOL;
//      ) RAISES {} =
//      BEGIN
//        INC(nSubStates);
//      END
//  
//    <* FATAL rdag_disp_STOP *>
//    BEGIN
//      (*
//        Note that "nStates" counts only the DAG states reachable by "raut_step"
//        chains, not those reachable by "rdag_last-rdag_rest" chains.
//  
//        Note also that "nArcs" is the sum of the outdegre of all
//        "raut_step"-reachable states of the rautuced_t; thus, DAG arcs that
//        are shared by more than one reachable state are counted
//        more than once.
//      *)
//      aut.raut_enum_states(base := base, enter := CountState);
//      aut.raut_enum_states(base := base, substates := TRUE, enter := CountSubState);
//      nStrings := 0;
//      nLetters := 0;
//      FOR i := 0 TO LAST(base) DO
//        nStrings := nStrings + aut.raut_num_suffs(base[i]);
//        nLetters := nLetters + aut.raut_num_suff_letters(base[i]);
//      END;
//      RETURN raut_counts_t{
//        strings := nStrings,
//        symbols := nLetters,
//        states := nStates,
//        substates := nSubStates,
//        arcs := nArcs,
//        finals := nFinals
//      }
//    END
//  
//  PROCEDURE raut_num_prefs(aut: T; s: raut_state_t): uint32_t =
//    BEGIN
//      IF s = raut_state_NULL OR s > aut.root THEN
//        RETURN 0
//      ELSIF s = aut.root THEN
//        RETURN 1
//      ELSE
//        ComputePrefixData(aut);
//        <* ASSERT aut.nPrefs # NIL *>
//        WITH np = aut.nPrefs^ DO
//          <* ASSERT LAST(np) >= aut.root *>
//          RETURN np[s]
//        END;
//      END
//    END
//  
//  PROCEDURE raut_num_suffs(aut: T; s: raut_state_t): uint32_t =
//    BEGIN
//      IF s = raut_state_NULL THEN
//        RETURN 0
//      ELSIF s = raut_state_UNIT THEN
//        RETURN 1
//      ELSE
//        <* ASSERT aut.nSuffs # NIL *>
//        WITH ns = aut.nSuffs^ DO
//          <* ASSERT LAST(ns) >= s *>
//          RETURN ns[s]
//        END;
//      END;
//    END
//  
//  PROCEDURE raut_num_pref_letters(aut: T; s: raut_state_t): uint32_t =
//    BEGIN
//      IF s = raut_state_NULL OR s > aut.root THEN
//        (* No paths from "aut.root" to "s":  *)
//        RETURN 0
//      ELSIF s = aut.root THEN
//        (* One path from "aut.root" to "s", of length 0: *)
//        RETURN 0
//      ELSE
//        ComputePrefixData(aut);
//        <* ASSERT aut.nPrefLetters # NIL *>
//        WITH nl = aut.nPrefLetters^ DO
//          <* ASSERT LAST(nl) >= aut.root *>
//          RETURN nl[s]
//        END;
//      END
//    END
//  
//  PROCEDURE raut_num_suff_letters(aut: T; s: raut_state_t): uint32_t =
//    BEGIN
//      IF s = raut_state_NULL OR s = raut_state_UNIT THEN
//        RETURN 0
//      ELSE
//        <* ASSERT aut.nSuffLetters # NIL *>
//        WITH nl = aut.nSuffLetters^ DO
//          <* ASSERT LAST(nl) >= s *>
//          RETURN nl[s]
//        END;
//      END;
//    END
//  
//  PROCEDURE raut_enum_prefs(aut: T; s: raut_state_t; action: raut_prefix_action_t) RAISES {rdag_disp_STOP} =
//    VAR rs: REF String := NEW(REF String, 100);
//  
//    PROCEDURE PrefixPush (* : DAG.raut_path_action_t *) (
//        len: uint32_t;
//        <*UNUSED*> org: raut_state_t;
//        <*UNUSED*> i: uint32_t;
//        a: DAG.Arc
//      ) RAISES {rdag_disp_SKIP, rdag_disp_STOP} =
//      BEGIN
//        IF a.rd = NullLetter THEN
//          action(SUBARRAY(rs^, 0, len));
//          RAISE rdag_disp_SKIP
//        ELSE
//          rdag_basics.h.ExpandString(rs, len+1);
//          rs[len] := a.rd
//        END
//      END
//  
//    BEGIN
//      IF s = raut_state_NULL OR s > aut.root THEN
//        RETURN
//      ELSIF s = aut.root THEN
//        action(String{});
//      ELSE
//        ComputePrefixData(aut);
//        aut.rdag.rdag_enum_paths(aut.rev[s], push := PrefixPush)
//      END;
//    END
//  
//  PROCEDURE raut_enum_suffs(aut: T; s: raut_state_t; action: raut_suffix_action_t) RAISES {rdag_disp_STOP} =
//  
//    PROCEDURE raut_suffix_action_t (* : raut_string-action_t *) (
//         READONLY w: String;
//         <*UNUSED*> dest: raut_state_t;
//         final: BOOL;
//      ) RAISES {rdag_disp_STOP} =
//      BEGIN
//        IF final THEN action(w) END
//      END
//  
//    BEGIN
//      aut.raut_enum_strings(s, enter := raut_suffix_action_t)
//    END
//  
//  PROCEDURE raut_print_prefs(
//      aut: T;
//      s: raut_state_t;
//      spr: rdag_string_printer.h.T;
//    ) =
//  
//    PROCEDURE PrintPrefix (* : raut_prefix_action_t *) (
//        READONLY w: String
//      ) RAISES {rdag_disp_STOP} =
//      BEGIN
//        spr.PutString(w, rev := TRUE)
//      END
//  
//    BEGIN
//      TRY aut.raut_enum_prefs(s, action := PrintPrefix) EXCEPT rdag_disp_STOP => (* Ok *) END;
//      spr.Reset();
//    END
//  
//  PROCEDURE raut_print_suffs(
//      aut: T;
//      s: raut_state_t;
//      spr: rdag_string_printer.h.T;
//    ) =
//  
//    PROCEDURE PrintSuffix (* : raut_suffix_action_t *) (
//        READONLY w: String
//      ) RAISES {rdag_disp_STOP} =
//      BEGIN
//        spr.PutString(w, rev := FALSE)
//      END
//  
//    BEGIN
//      TRY aut.raut_enum_suffs(s, action := PrintSuffix) EXCEPT rdag_disp_STOP => (* Ok *) END;
//      spr.Reset();
//    END
//  
//  PROCEDURE raut_first_prefix(aut: T; s: raut_state_t; e: rdag_encoding.h.T): TEXT =
//  
//    VAR wr := TextWr.New();
//  
//    PROCEDURE FP(t: DAG.raut_state_t) =
//    (*
//      Prints into "wr" the string of the first path from "t"
//      to an initial state (recognized by its first arc being labelled
//      with NullLetter). *)
//      <* FATAL Wr.Failure, Thread.Alerted *>
//      BEGIN
//        <* ASSERT t # raut_state_NULL *>
//        WITH a = aut.rdag.rdag_first(t) DO
//          IF a.rd = NullLetter THEN
//            RETURN
//          ELSE
//            FP(a.dest);
//            e.PrintLetter(wr, a.rd);
//          END;
//        END;
//      END
//  
//    BEGIN
//      <* ASSERT s # raut_state_NULL *>
//      <* ASSERT s <= aut.root *>
//      ComputePrefixData(aut);
//      FP(aut.rev[s]);
//      RETURN TextWr.ToText(wr)
//    END
//  
//  PROCEDURE raut_first_suffix(aut: T; s: raut_state_t; e: rdag_encoding.h.T): TEXT =
//    VAR t: raut_state_t := s;
//    <* FATAL Wr.Failure, Thread.Alerted *>
//    BEGIN
//      WITH wr = TextWr.New() DO
//        LOOP
//          <* ASSERT t # raut_state_NULL *>
//          WITH a = aut.dag.rdag_first(t) DO
//            IF a.rd = NullLetter THEN RETURN TextWr.ToText(wr) END;
//            e.PrintLetter(wr, a.rd);
//            t := a.dest
//          END
//        END;
//      END;
//    END
//  
//  PROCEDURE raut_full_mark(aut: T; s: raut_state_t; e: rdag_encoding.h.T; sep: TEXT := ":"): TEXT =
//    BEGIN
//      RETURN
//        aut.raut_first_prefix(s, e) & sep & aut.raut_first_suffix(s, e)
//    END
//    
//  VAR (*CONST*) Defaultrdag_encoding.h: rdag_encoding.h.T := Plainrdag_encoding.h.New();
//  
//  PROCEDURE raut_build(
//      aut: T;
//      next: raut_next_string_proc_t;        (* Client input procedure *)
//      wr: FILE * := NIL;             (* Writer for progress report, etc: *)
//      e: rdag_encoding.h.T := NIL;        (* rdag_symbol_t/CHAR encoding for printout *)
//      reportInterval: POS := 1000; (* Print a report every this many input strings *)
//      flagRedundant: BOOL := TRUE; (* TRUE to print warnings on redundant operations *)
//    ) RAISES {rdag_disp_STOP} =
//  
//    VAR
//      nStrings: ARRAY BOOL OF uint32_t := ARRAY OF uint32_t {0, 0};
//      nLetters: ARRAY BOOL OF uint32_t := ARRAY OF uint32_t {0, 0};
//  
//    VAR
//      lastCrunchMaxState: raut_state_t := raut_state_NULL;  (* aut.raut_max_state() after last GC run *)
//  
//    PROCEDURE CrunchIt() =
//    (*
//      Does aut.raut_crunch(), printing status reports. *)
//      <* FATAL Wr.Failure, Thread.Alerted *>
//      BEGIN
//        WITH size = aut.raut_max_alloc_state() + 1 DO
//          IF wr # NIL THEN
//            IF NOT TimeToReport() THEN PrintStatusReport() END;
//            Wr.PutText(wr, "    * (crunching, alloc = " &  Fmt.Int(size) & "...");
//            Wr.Flush(wr);
//          END;
//          aut.raut_crunch(pSt);
//          IF wr # NIL THEN 
//            Wr.PutText(wr, ")\n");
//            PrintStatusReport();
//          END;
//        END;
//        lastCrunchMaxState := aut.raut_max_state();
//      END
//  
//    PROCEDURE ExpandIt() =
//    (*
//      Expands the automaton, printing some noise: *)
//      <* FATAL Wr.Failure, Thread.Alerted *>
//      BEGIN
//        WITH
//          oldSize = aut.raut_max_alloc_state() + 1,
//          newSize = ComputeNewSize(oldSize)
//        DO
//          IF wr # NIL THEN
//            IF NOT TimeToReport() THEN PrintStatusReport() END;
//            Wr.PutText(wr, "    * (expanding from ");
//            Wr.PutText(wr, Fmt.Int(oldSize));
//            Wr.PutText(wr, " to ");
//            Wr.PutText(wr, Fmt.Int(newSize));
//            Wr.PutText(wr, "...");
//            Wr.Flush(wr);
//          END;
//          <* FATAL rdag_basics.h.raut_error_FULL *>
//          BEGIN
//            aut.raut_expand(newSize);
//          END;
//          IF wr # NIL THEN 
//            Wr.PutText(wr, ")\n");
//            PrintStatusReport();
//          END;
//        END;
//      END
//  
//    PROCEDURE ComputeNewSize(oldSize: uint32_t): uint32_t =
//    (*
//      Chooses a suitable new size for automatic expansion.
//      *)
//      BEGIN
//        RETURN
//          MIN(
//            MIN(
//                  1 + oldSize + oldSize,
//              10001 + oldSize + oldSize DIV 2
//            ),
//            MIN(
//              20001 + oldSize + oldSize DIV 3,
//              50001 + oldSize + oldSize DIV 4
//            )
//          )
//      END
//  
//    VAR
//      pN: uint32_t := 0;                              (* Number of pending raut_set_arc actions *)
//      pSt: REF States := NEW(REF States, 100);   (* raut_state_t arguments for pending actions *)
//      pLet: REF String := NEW(REF String, 100);  (* rdag_symbol_t arguments for pending SetArcs *)
//        (*
//          For efficiency reasons, the strings returned by "next" are not added or
//          deleted right away.  Instead, we keep a list "pa[0..pN-1]" of "pending
//          raut_set_arc actions" to be performed later.
//          
//          The "i"th pending raut_set_arc action is described by the pair /pSt[i],
//          pLet[i]". The action consists in computing a new state "new[i]/ that 
//          is like "pSt[i]" except that under the symbol
//          "pLet[i]" it goes to the state "new[i+1]". As a special case,
//          "new[pN]" is defined to be "pSt[pN]" itself, and "pLet[pN]" is
//          not used.
//          
//          Thus, to flush out all pending actions we must repeat
//          
//  |          pSt[pN-1] := aut.raut_set_arc(pSt[pN-1], pLet[pN-1], s);
//  |          DEC(pN)
//  
//          until pN = 0, and finally set the automaton's root to the
//          resuting state "pSt[0]". *)
//      
//    PROCEDURE DoSetFinal(s: raut_state_t; final: BOOL): raut_state_t =
//    (*
//      Computes aut.raut_set_final(s, final), but expands "aut" if necessary
//      (instead of raising "raut_error_FULL"). *)
//      BEGIN
//        LOOP
//          TRY
//            RETURN aut.raut_set_final(s, final);
//          EXCEPT 
//            raut_error_FULL => ExpandIt();
//          END;
//        END;
//      END
//  
//    PROCEDURE DoOnePendingAction() =
//    (*
//      Performs the last pending "raut_set_arc" action. Expands if necessary. 
//      *)
//      BEGIN
//        <* ASSERT pN > 0 *>
//        WITH 
//          st = pSt^, let = pLet^, 
//          n1 = pN-1 
//        DO
//          LOOP
//            TRY
//              st[n1] := aut.raut_set_arc(st[n1], let[n1], st[pN]);
//              pN := n1;
//              RETURN
//            EXCEPT 
//              raut_error_FULL => ExpandIt();
//            END;
//          END;
//        END
//      END
//  
//    PROCEDURE FlagRedundant(READONLY w: String; add: BOOL) =
//      <* FATAL Wr.Failure, Thread.Alerted *>
//      BEGIN
//        Wr.PutText(wr, "  ** redundant command: ");
//        PrintCommand(w, add);
//        Wr.PutChar(wr, '\n');
//      END
//  
//    PROCEDURE AddOrSubIt(READONLY w: String; add: BOOL := TRUE) =
//    (*
//      Adds or deletes "w", crunching and/or expanding if necessary:
//      *)
//      BEGIN
//        (* raut_crunch the automaton, if it looks worth it: *)
//        WITH alloc = aut.raut_max_alloc_state() DO
//          IF aut.raut_max_state() + 30*NUMBER(w) >= alloc
//          AND (aut.raut_max_state() - lastCrunchMaxState) >= alloc DIV 10
//          THEN
//            (* We are close to the allocated size,
//              and about 10% of the current dag states were
//              created since the last CG run. Better run GC again... *)
//            CrunchIt();
//          END;
//        END;
//  
//        (* Now add the string, and expand if doesn't fit: *)
//        VAR np: CARDINAL := 0;
//        BEGIN
//          (* rdag_disp_SKIP pending actions that match the symbols of w: *)
//          WITH maxp = MIN(pN, NUMBER(w)), let = pLet^ DO
//            WHILE np < maxp AND w[np] = let[np] DO INC(np) END
//          END;
//          
//          (* Flush any remaining actions: *)
//          IF pN > np THEN
//            REPEAT DoOnePendingAction() UNTIL pN <= np
//          END;
//          
//          (* Stack new actions corresponding to remaining symbols of w: *)
//          rdag_basics.h.ExpandString(pLet, NUMBER(w));
//          ExpandStates(pSt, NUMBER(w) + 1);
//          <* ASSERT pN = np *>
//          WITH st = pSt^, let = pLet^ DO
//            FOR i := np TO LAST(w) DO
//              let[i] := w[i];
//              st[i+1] := aut.raut_step(st[i], w[i]);
//            END;
//            pN := NUMBER(w);
//            IF add = aut.raut_is_final(pSt[pN]) THEN
//              (* Command was superflous *)
//              IF wr # NIL AND flagRedundant THEN
//                FlagRedundant(w, add)
//              END;
//            ELSE
//              pSt[pN] := DoSetFinal(pSt[pN], add)
//            END;
//          END;
//        END;
//      END
//  
//    VAR
//      refw: REF String := NEW(REF String, 100);  (* String buffer *)
//      len: uint32_t;                                  (* Length of string *)
//      add: BOOL;                                 (* Add/delete flag for refw^ *)
//      
//    PROCEDURE PrintStatusReport() =
//      <* FATAL Wr.Failure, Thread.Alerted *>
//      BEGIN
//        Wr.PutText(wr, "    * ");
//        WITH 
//          ns = nStrings[FALSE] + nStrings[TRUE],
//          nl = nLetters[FALSE] + nLetters[TRUE]
//        DO
//          Wr.PutText(wr, Fmt.Pad(Fmt.Int(ns), 8) &  " strings ");
//          Wr.PutText(wr, Fmt.Pad(Fmt.Int(nl), 8) &  " symbols ");
//        END;
//        Wr.PutText(wr, Fmt.Pad(Fmt.Int(aut.raut_max_state()), 8) & " dag states  ");
//        PrintCommand(SUBARRAY(refw^, 0, len), add);
//        Wr.PutText(wr, "\n");
//        Wr.Flush(wr);
//      END
//  
//    PROCEDURE PrintCommand(READONLY w: String; add: BOOL) =
//      <* FATAL Wr.Failure, Thread.Alerted *>
//      BEGIN
//        Wr.PutChar(wr, ARRAY BOOL OF CHAR{'-', '+'}[add]);
//        Wr.PutChar(wr, ' ');
//        TRY e.PrintString(wr, w) EXCEPT rdag_encoding.h.BadString => (*IGNORE*) END;
//      END
//  
//    PROCEDURE PrintFinalReport() =
//      <* FATAL Wr.Failure, Thread.Alerted *>
//      BEGIN
//        Wr.PutText(wr, "\n");
//  
//        Wr.PutText(wr, "Input statistics:\n");
//        Wr.PutText(wr, "\n");
//  
//        Wr.PutText(wr, "add:  ");
//        Wr.PutText(wr, Fmt.Pad(Fmt.Int(nStrings[TRUE]), 8) & " strings ");
//        Wr.PutText(wr, Fmt.Pad(Fmt.Int(nLetters[TRUE]), 8) & " symbols ");
//        Wr.PutText(wr, "\n");
//  
//        Wr.PutText(wr, "sub:  ");
//        Wr.PutText(wr, Fmt.Pad(Fmt.Int(nStrings[FALSE]), 8) & " strings ");
//        Wr.PutText(wr, Fmt.Pad(Fmt.Int(nLetters[FALSE]), 8) & " symbols ");
//        Wr.PutText(wr, "\n");
//  
//        Wr.PutText(wr, "------");
//        Wr.PutText(wr, "--------");
//        Wr.PutText(wr, "---------");
//        Wr.PutText(wr, "--------");
//        Wr.PutText(wr, "---------");
//        Wr.PutText(wr, "\n");
//        WITH 
//          ns = nStrings[FALSE] + nStrings[TRUE],
//          nl = nLetters[FALSE] + nLetters[TRUE]
//        DO
//          Wr.PutText(wr, "tot:  ");
//          Wr.PutText(wr, Fmt.Pad(Fmt.Int(ns), 8) & " strings ");
//          Wr.PutText(wr, Fmt.Pad(Fmt.Int(nl), 8) & " symbols ");
//        END;
//        Wr.PutText(wr, "\n");
//        Wr.PutText(wr, "\n");
//        WITH ct = aut.raut_count(ARRAY OF raut_state_t{aut.raut_root()}) DO
//          raut_print_counts(wr, ct)
//        END;
//      END
//  
//    PROCEDURE TimeToReport(): BOOL =
//      BEGIN
//        RETURN (nStrings[FALSE] + nStrings[TRUE]) MOD reportInterval = 0
//      END
//    
//    BEGIN (*raut_build*)
//      pN := 0;
//      pSt[0] := aut.raut_root();
//      IF wr # NIL THEN
//        IF e = NIL THEN e := Defaultrdag_encoding.h END;
//      END;
//      TRY
//        TRY
//          LOOP
//            next((*IO*) refw, (*OUT*) len, (*OUT*) add);
//            WITH w = SUBARRAY(refw^, 0, len) DO
//              AddOrSubIt(w, add);
//              INC(nStrings[add]);
//              INC(nLetters[add], NUMBER(w));
//              IF wr # NIL AND TimeToReport() THEN 
//                PrintStatusReport(); 
//              END;
//            END
//          END
//        EXCEPT
//        | Done => (* Ok *)
//        END;
//      FINALLY (* Normally, or in case of "rdag_disp_STOP" *)
//        WHILE pN > 0 DO DoOnePendingAction() END;
//        aut.raut_set_root(pSt[0])
//      END;
//      CrunchIt();
//      IF wr # NIL THEN PrintFinalReport() END;
//    END
//  
//  (************************)
//  (* AUXILIARY OPERATIONS *)
//  (************************)
//  
//  PROCEDURE FromDAG(dag: rdag_t; doc: TEXT; root: DAG.raut_state_t): T =
//  (*
//    Bulds a rautuced_t given the underlying DAG and the root state.
//    The DAG must satisfy the implementation conventions about the
//    use of NullLetter and the ordering of arc labels. *)
//    BEGIN
//      WITH
//        aut = NEW(T,
//          doc := doc,
//          dag := dag,
//          root := root,
//          nSuffs := NIL,
//          nSuffLetters := NIL,
//          prefRoot := raut_state_NULL,
//          rdag := NIL,
//          rev := NIL,
//          dir := NIL,
//          nPrefs := NIL,
//          nPrefLetters := NIL
//        )
//      DO
//        MakeUnitState(dag);
//        ComputeNSuffs(aut);
//        RETURN aut;
//      END
//    END
//  
//  PROCEDURE MakeUnitState(dag: rdag_t) =
//  (*
//    Creates the "unit" state in the "dag", if necesary: *)
//    <* FATAL rdag_basics.h.raut_error_FULL *>
//    BEGIN
//      WITH 
//        unit = dag.Append(
//          last := DAG.Arc{rd := 0, wr := 0, dest := DAG.raut_state_NULL},
//          rest := DAG.raut_state_NULL
//        )
//      DO
//        <* ASSERT unit = raut_state_UNIT *>
//      END
//    END
//  
//  PROCEDURE DiscardPrefixData(aut: T) =
//    BEGIN
//      aut.prefRoot := raut_state_NULL;
//      aut.rdag := NIL;
//      aut.rev := NIL;
//      aut.dir := NIL;
//      aut.nPrefs := NIL;
//      aut.nPrefLetters := NIL;
//    END
//  
//  PROCEDURE ComputePrefixData(aut: T) =
//  (*
//    If the prefix tables of "aut" ("rdag", "rev", "dir" and "nPrefs")
//    are missing or out of date (e.g., because the root has changed),
//    recomputes them for the current root.
//    *)
//    BEGIN
//      IF aut.root = raut_state_NULL OR aut.prefRoot = aut.root THEN
//        (* Current prefix data is still OK *)
//        RETURN
//      ELSE
//        DiscardPrefixData(aut);
//        ComputeReverseDAG(aut);
//        ComputeNPrefs(aut);
//        aut.prefRoot := aut.root;
//      END
//    END
//  
//  PROCEDURE ComputeReverseDAG(aut: T) =
//  (*
//    Computes "aut.rdag", "aut.rev", "aut.dir" from "aut.dag" and "aut.root".
//    *)
//    VAR maxRev: DAG.raut_state_t;
//    <* FATAL rdag_basics.h.raut_error_FULL *>
//    BEGIN
//      <* ASSERT aut.root # raut_state_NULL *>
//      (*
//        Note that the reverse DAG may be bigger than "aut.dag",
//        because one arc of "aut" that is shared by two
//        reachable states will give rise to two distinct reverse arcs.
//        *)
//      WITH
//        totArcs = aut.raut_count(ARRAY OF raut_state_t{aut.root}).arcs,
//        rdag = DAG.New(size := totArcs + 1),
//        rev = NEW(REF ARRAY OF DAG.raut_state_t, aut.root+1),
//        r = rev^
//      DO
//        FOR t := 0 TO aut.root DO r[t] := raut_state_NULL END;
//        (* Create one transition from the root state to raut_state_NULL with NullLetter: *)
//        r[aut.root] := rdag.Append(
//          last := DAG.Arc{rd := NullLetter, wr := NullLetter, dest := raut_state_NULL},
//          rest := raut_state_NULL
//        );
//  
//        (*
//          For every state "s" reachable from "root", and every
//          proper successor "t" of "s" in "aut", add to "rdag" an arc from
//          "rev[t]" to "rev[s]", and update rev[t] accordingly.
//          Note that it is important to process the states from high to low,
//          so that rev[s] is stable by the time we add the reverse arcs into it. *)
//  
//        maxRev := raut_state_NULL;
//        FOR s := aut.root TO 1 BY -1 DO
//          IF r[s] # raut_state_NULL THEN
//            (* "s" is reachable from the root state. *)
//  
//            (* Update "maxRev": *)
//            maxRev := MAX(maxRev, r[s]);
//  
//            (* Enumerate its outgoing arcs, and add them to the reverse dag: *)
//            VAR t: raut_state_t := s;
//            BEGIN
//              WHILE aut.raut_has_arcs(t) DO
//                WITH a = aut.rdag_last(t) DO
//                  r[a.dest] := rdag.Append(
//                    last := DAG.Arc{rd := a.symbol, wr := NullLetter, dest := r[s]},
//                    rest := r[a.dest]
//                  );
//                  t := aut.rdag_rest(t)
//                END;
//              END;
//            END;
//          END;
//        END;
//  
//        (* Store results in "aut": *)
//        aut.rdag := rdag;
//        aut.rev := rev;
//  
//        (* raut_build table from reversed state to direct state: *)
//        WITH
//          dir = NEW(REF ARRAY OF raut_state_t, maxRev + 1),
//          d = dir^
//        DO
//          FOR t := 0 TO maxRev DO d[t] := raut_state_NULL END;
//          FOR s := 0 TO aut.root DO
//            WITH t = r[s] DO
//              IF t # raut_state_NULL THEN
//                <* ASSERT d[t] = raut_state_NULL *>
//                d[t] := s
//              END
//            END
//          END;
//          aut.dir := dir
//        END;
//      END;
//    END
//  
//  PROCEDURE ComputeNPrefs(aut: T) =
//  (*
//    Computes the tables "nPrefs" and "nPrefLetters", for the current root state,
//    from 
//    *)
//    VAR t: raut_state_t;
//    BEGIN
//      <* ASSERT aut.root # raut_state_NULL *>
//      (* Allocates the vector if necessary. *)
//      WITH
//        minSize = aut.root + 1
//      DO
//        IF aut.nPrefs = NIL OR NUMBER(aut.nPrefs^) < minSize THEN
//          aut.nPrefs := NEW(REF ARRAY OF uint32_t, minSize)
//        END;
//        IF aut.nPrefLetters = NIL OR NUMBER(aut.nPrefLetters^) < minSize THEN
//          aut.nPrefLetters := NEW(REF ARRAY OF uint32_t, minSize)
//        END;
//      END;
//  
//      WITH
//        np = aut.nPrefs^,
//        nl = aut.nPrefLetters^,
//        maxState = aut.root
//      DO
//        FOR s := raut_state_NULL TO maxState DO 
//          np[s] := 0; nl[s] := 0 
//        END;
//        IF maxState > raut_state_NULL THEN np[maxState] := 1 END;
//        FOR s := maxState TO 1 BY -1 DO
//          WITH  nps = np[s], nls = nl[s] DO
//            IF nps # 0 THEN
//              t := s;
//              WHILE aut.raut_has_arcs(t)  DO
//                WITH
//                  d = (aut.rdag_last(t)).dest
//                DO
//                  <* ASSERT d < s *>
//                  INC(np[d], nps);
//                  INC(nl[d], nls + nps);
//                END;
//                t := aut.rdag_rest(t)
//              END
//            END
//          END;
//        END
//      END
//    END
//  
//  PROCEDURE ComputeNSuffs(aut: T) =
//  (*
//    Computes the tables "nSuffs" and "nSuffLetters", from "aut.dag". *)
//    BEGIN
//      (* Allocates the vectors if necessary. *)
//      WITH
//        minSize = aut.raut_max_alloc_state() + 1
//      DO
//        IF aut.nSuffs = NIL OR NUMBER(aut.nSuffs^) < minSize THEN
//          aut.nSuffs := NEW(REF ARRAY OF uint32_t, minSize)
//        END;
//        IF aut.nSuffLetters = NIL OR NUMBER(aut.nSuffLetters^) < minSize THEN
//          aut.nSuffLetters := NEW(REF ARRAY OF uint32_t, minSize)
//        END;
//      END;
//  
//      (* Computes suffix counts and sizes, by a postorder scan *)
//      WITH
//        ns = aut.nSuffs^,
//        nl = aut.nSuffLetters^,
//        maxState = aut.raut_max_state()
//      DO
//        <* ASSERT maxState >= raut_state_UNIT *>
//        ns[raut_state_NULL] := 0;  nl[raut_state_NULL] := 0;
//        ns[raut_state_UNIT] := 1;  nl[raut_state_UNIT] := 0;
//        FOR s := 2 TO maxState DO
//          WITH
//            d  = (aut.rdag_last(s)).dest,
//            r  = aut.rdag_rest(s)
//          DO
//            <* ASSERT (d < s) AND (r < s) *>
//            ns[s] := ns[d] + ns[r];
//            nl[s] := nl[d] + ns[d] + nl[r]
//          END;
//        END;
//        (* Clear out remaining counts, just for tidiness: *)
//        FOR i := maxState + 1 TO LAST(ns) DO ns[i] := 0; nl[i] := 0 END;
//      END;
//    END
//  
//  CONST DocPrefix: CHAR = '|';
//    
//  PROCEDURE DumpDoc(wr: FILE *; doc: TEXT) =
//  (* 
//    Writes the given "doc" text to "wr", with a DocPrefix
//    in front of every line.  Supplies a final '\n' if the text is 
//    non-empty but does not end with newline. *)
//    
//    VAR rd: FILE * := TextRd.New(doc);
//    
//    PROCEDURE CopyLine() RAISES {Rd.EndOfFile} =
//    (*
//      Copy one line from "rd" to "wr", prefixed by DocPrefix. 
//      Supplies a final '\n' if next line exists but does not end with newline.
//      Raises Rd.EndOfFile if there are no more lines in "rd". *)
//      
//      <* FATAL Rd.Failure, Wr.Failure, Thread.Alerted *>
//      VAR c: CHAR;
//      BEGIN
//        c := Rd.GetChar(rd); (* If EOF here, propagate to caller *)
//        Wr.PutChar(wr, DocPrefix);
//        Wr.PutChar(wr, c);
//        WHILE c # '\n' DO
//          TRY c := Rd.GetChar(rd) EXCEPT Rd.EndOfFile => c := '\n' END;
//          Wr.PutChar(wr, c)
//        END
//      END
//  
//    BEGIN
//      TRY LOOP CopyLine() END EXCEPT Rd.EndOfFile => (* Ok *) END;
//    END
//  
//  PROCEDURE LoadDoc(rd: FILE *): TEXT =
//  (*
//    Reads zero or more lines from "rd" that begin with DocPrefix, strips the
//    leading DocPrefix of each line, and returns the concatenation of those lines
//    as a single TEXT with embedded and terminating newline chars. *)
//  
//    VAR wr: FILE * := TextWr.New();
//  
//    PROCEDURE CopyLine() RAISES {Rd.EndOfFile} =
//    (*
//      Copy one DocPrefix line from "rd" to "wr", removing the DocPrefix
//      but leaving the final (mandatory) newline.
//      Raises Rd.EndOfFile if "rd" is exhausted or the next char is 
//      not DocPrefix. *)
//      <* FATAL Rd.Failure, Wr.Failure, Thread.Alerted *>
//      <* FATAL MissingFinalNewLine *>
//      VAR c: CHAR;
//      BEGIN
//        c := Rd.GetChar(rd); (* If EOF here, propagate to caller *)
//        IF c # DocPrefix THEN Rd.UnGetChar(rd); RAISE Rd.EndOfFile END;
//        REPEAT
//          TRY c := Rd.GetChar(rd) EXCEPT Rd.EndOfFile => RAISE MissingFinalNewLine END;
//          Wr.PutChar(wr, c)
//        UNTIL c = '\n'
//      END
//  
//    VAR twr: TextWr.T;
//        txt: TEXT;
//    BEGIN
//      TRY LOOP CopyLine() END EXCEPT Rd.EndOfFile => (* Ok *) END;
//      twr := wr;
//      txt := TextWr.ToText(twr);
//      RETURN txt
//    END  
//  
//  PROCEDURE raut_print_counts(wr: FILE *; READONLY ct: raut_counts_t) =
//    <* FATAL Wr.Failure, Thread.Alerted *>
//    BEGIN
//      Wr.PutText(wr, " strings  symbols    states   finals     arcs  sub-sts lets/arc\n");
//      Wr.PutText(wr, "-------- --------  -------- -------- -------- -------- --------\n");
//      Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.strings),   8));
//      Wr.PutChar(wr, ' ');
//      Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.symbols),   8));
//      Wr.PutChar(wr, ' ');
//      Wr.PutChar(wr, ' ');
//      Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.states),    8));
//      Wr.PutChar(wr, ' ');
//      Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.finals),    8));
//      Wr.PutChar(wr, ' ');
//      Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.arcs),      8));
//      Wr.PutChar(wr, ' ');
//      Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.substates), 8));
//      Wr.PutChar(wr, ' ');
//      Wr.PutText(wr, 
//        Fmt.Pad(Fmt.Real(FLOAT(ct.symbols)/FLOAT(ct.arcs), Fmt.Style.Fix, 3), 8)
//      );
//      Wr.PutChar(wr, '\n');
//    END
//    
//  PROCEDURE ExpandStates(VAR s: REF States; n: CARDINAL) =
//    BEGIN
//      WITH nold = NUMBER(s^) DO
//        IF nold < n THEN
//          WITH
//            r = NEW(REF States, MAX(n + 10, 2*nold))
//          DO
//            SUBARRAY(r^, 0, nold) := s^;
//            s := r
//          END
//        END
//      END
//    END
//  
//  BEGIN
//    <* ASSERT raut_state_NULL = DAG.raut_state_NULL *>
//    <* ASSERT NullLetter < FIRST(rdag_symbol_t) *>
//  END Reduced.
