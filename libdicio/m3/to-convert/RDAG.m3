MODULE DAG;

(* See the copyright and disclaimer note at the end of this file. *)

IMPORT Word, Text, Rd, Wr, Thread;
FROM Basics IMPORT NAT;
FROM Basics IMPORT Full, Skip, Abort, RdNat, ReadParam, WrNat;

CONST
  StatePBits = 24;    (* Num bits to use for packed state numbers *)
  StatePMax = Word.Shift(1, StatePBits) - 1; (* Maximum packed state number *)

TYPE
  StateP = BITS StatePBits FOR [0..StatePMax];

  Entry = RECORD  (* Data for one state/arc *)
      dest: StateP;      (* Destination of last arc out of this state *)
      rd: Symbol;        (* rd-label of last arc out of this state *)
      rest: StateP;      (* State with same out-arcs except last *)
      wr: Symbol;        (* wr-label of last arc out of this state *)
      next: StateP;      (* Next state with same hash value *)
      foo: Symbol := 0;  (* A filler field *)
    END;

  Entries = ARRAY OF Entry;  (* Entry 0 (NullState) is not used *)

REVEAL
  T = Public BRANDED OBJECT
      maxState: NAT;          (* Number of largest-numbered created state *)
      e: REF Entries;         (* Arc/state data. Entries "[1..maxState]" are in use. *)
      h: REF ARRAY OF State;  (* Hash table for "Append" *)
    OVERRIDES
      MaxState  := MaxState;
      Last      := Last;
      Rest      := Rest;
      Append    := Append;
      First     := First;
      OutDeg    := OutDeg;
      EnumPaths := EnumPaths;
      MaxAllocState := MaxAllocState;
      Expand    := Expand;
      Crunch    := Crunch;
    END;

PROCEDURE MaxState(dag: T): State =
  BEGIN
    RETURN dag.maxState
  END MaxState;

PROCEDURE MaxAllocState(dag: T): NAT =
  BEGIN
    RETURN LAST(dag.e^)
  END MaxAllocState;

PROCEDURE Last(dag: T; s: State): Arc =
  BEGIN
    <* ASSERT s > 0 AND s <= dag.maxState *>
    WITH es = dag.e[s] DO
      RETURN Arc{rd := es.rd, wr := es.wr, dest := es.dest}
    END
  END Last;

PROCEDURE Rest(dag: T; s: State): State =
  BEGIN
    <* ASSERT s > 0 AND s <= dag.maxState *>
    RETURN dag.e[s].rest
  END Rest;

PROCEDURE OutDeg(dag: T; s: State): NAT =
  VAR t: State := s; n: NAT := 0;
  BEGIN
    WHILE t # NullState DO INC(n); t := dag.Rest(t) END;
    RETURN n
  END OutDeg;

PROCEDURE First(dag: T; s: State): Arc =
  VAR p, t: State;
  BEGIN
    <* ASSERT s # NullState *>
    t := s;
    REPEAT p := t; t := dag.Rest(t);  UNTIL t = NullState;
    RETURN dag.Last(p)
  END First;

PROCEDURE Append(dag: T; last: Arc; rest: State): State RAISES {Full} =
  BEGIN
    <* ASSERT last.dest <= dag.maxState *>
    <* ASSERT rest <= dag.maxState *>
    (* See if state already exists: *)
    WITH
      e = dag.e^,
      h = dag.h^,
      hash = AppendHash(
        Arc{rd := last.rd, wr := last.wr, dest := last.dest},
        rest,
        NUMBER(h)
      )
    DO
      (* Check if the required state already exists, using the hash table: *)
      VAR s: State := h[hash];
          s1, s2: State := NullState;
      BEGIN
        LOOP
          IF s = NullState THEN EXIT END;
          WITH es = e[s] DO
            IF es.rest = rest
            AND es.dest = last.dest
            AND es.rd = last.rd
            AND es.wr = last.wr
            THEN
              EXIT
            END
          END;
          s2 := s1; s1 := s; s := e[s].next
        END;
        IF s # NullState THEN
          IF s1 # NullState THEN
            (* Swap s with its predecessor (is this worth doing?) *)
            IF s2 = NullState THEN h[hash] := s ELSE e[s2].next := s END;
            e[s1].next := e[s].next;
            e[s].next := s1
          END;
          RETURN s
        END
      END;

      (* No luck, must create a new entry: *)
      IF dag.maxState + 1 > LAST(dag.e^) THEN
        (* Really bad luck, ran out of space: *)
        RAISE Full
      END;

      WITH
        s = dag.maxState + 1
      DO
        dag.maxState := s;
        WITH es = e[s] DO
          es.rd := last.rd;
          es.wr := last.wr;
          es.dest := last.dest;
          es.rest := rest;
          es.next := h[hash];
        END;
        h[hash] := s;
        RETURN s
      END
    END
  END Append;

PROCEDURE AppendHash(last: Arc; rest: State; size: NAT): NAT =
(*
  Computes a hash value in "[0..size-1]" for a state whose last arc
  is "last" and whose "Rest" is "rest".
  *)
  CONST
    (* The following are standard random numbers: 8-*)
    rdM = 418 - 1;
    destM = 4615;
    wrM = 1003;
    restM = 480317;
  BEGIN
    WITH
      rdH = Word.Times(last.rd, rdM),
      wrH = Word.Times(last.wr, wrM),
      destH = Word.Times(last.dest, destM),
      restH = Word.Times(rest, restM)
    DO
      RETURN Word.Mod(Word.Plus(Word.Plus(rdH, wrH), Word.Plus(destH, restH)), size)
    END
  END AppendHash;

PROCEDURE Expand(dag: T; newSize: NAT) RAISES {Full} =
  BEGIN
    IF newSize < dag.maxState THEN RAISE Full END;
    WITH
      e = dag.e,
      h = dag.h,
      eNew = NEW(REF Entries, newSize),
      hNew = NEW(REF ARRAY OF State, HashTableSize(newSize))
    DO
      WITH eNew = eNew^, e = e^ DO
        FOR i := 0 TO dag.maxState DO
          eNew[i] := e[i];
          eNew[i].next := NullState; (* Hash lists must be rebuilt *)
        END;
      END;
      e := eNew;
      h := hNew;
      Rehash(dag);
    END
  END Expand;

PROCEDURE HashTableSize(numEntries: NAT): NAT =
(*
  Computes a suitable hash table size for a DAG with "numEntries"
  allocated entries.
  *)
  BEGIN
    RETURN Word.Or(1, numEntries DIV 2)
  END HashTableSize;

PROCEDURE Rehash(dag: T) =
  BEGIN
    WITH e = dag.e^, h = dag.h^ DO
      (* Rehash states into new hash table *)
      FOR hash := 0 TO LAST(h) DO h[hash] := NullState END;
      FOR i := 1 TO dag.maxState DO
        WITH ei = e[i] DO
          <* ASSERT ei.rest < i *>
          <* ASSERT ei.dest < i *>
          WITH
            hash = AppendHash(
              Arc{rd := ei.rd, wr := ei.wr, dest := ei.dest},
              ei.rest,
              NUMBER(h)
            )
          DO
            ei.next := h[hash];
            h[hash] := i
          END
        END
      END
    END
  END Rehash;

PROCEDURE EnumPaths(
    dag: T;
    s: State;
    enter: StateAction := NIL;
    push, pop: ArcAction := NIL;
    exit: StateAction := NIL;
  ) RAISES {Abort} =

  VAR len: NAT := 0;

  PROCEDURE DoEnumPaths(t: State) RAISES {Abort} =
  (*
    Does EnumPaths starting at "t", a generic sucessor
    of "s". Assumes "t" is not "NullState". *)

    VAR i: NAT := 0;

    PROCEDURE EnumRest(r: State) RAISES {Skip, Abort} =
    (*
      Calls "enter" on "t" and enumerates a given prefix
      "r" of the arcs out of "t". *)

      BEGIN
        IF r = NullState THEN
          IF enter # NIL THEN enter(len, t) END;
        ELSE
          EnumRest(dag.Rest(r));
          WITH a = dag.Last(r) DO
            TRY
              IF push # NIL THEN push(len, t, i, a) END;
              INC(len);
              DoEnumPaths(a.dest);
              DEC(len);
            EXCEPT
              Skip => (* Ok *)
            END;
            IF pop # NIL THEN pop(len, t, i, a) END;
            INC(i);
          END;
        END
      END EnumRest;

    BEGIN (* DoEnumPaths *)
      <* ASSERT t # NullState *>
      TRY EnumRest(t) EXCEPT Skip => (* Ok *) END;
      IF exit # NIL THEN
        <* FATAL Skip *>
        BEGIN
          exit(len, t) (* Shouldn't raise "Skip" *)
        END
      END;
    END DoEnumPaths;

  BEGIN
    <* ASSERT s <= dag.maxState *>
    DoEnumPaths(s)
  END EnumPaths;

PROCEDURE New(size: NAT): T =
  BEGIN
    WITH
      e = NEW(REF Entries, size+1),
      h = NEW(REF ARRAY OF State, HashTableSize(NUMBER(e^))),
      dag = NEW(T,
         e := e,
         h := h,
         maxState := 0
       )
    DO
      (* initialize hash table *)
      WITH h = h^ DO
        FOR i := 0 TO LAST(h) DO h[i] := NullState END
      END;
      (* Initialize entries, just in case: *)
      WITH e = e^ DO
        FOR i := 0 TO LAST(e) DO
          e[i] := Entry{
            rd := 0,
            dest := 0,
            wr := 0,
            rest := 0,
            next := 0
          }
        END
      END;
      (* Note: "e[0]" should never be used. *)
      RETURN dag
    END
  END New;

PROCEDURE Copy(
    from: T;
    to: T;
    s: State;
    VAR (*IO*) map: REF ARRAY OF State;
  ): State RAISES {Full} =
  BEGIN
    (* We could do a linear scan instead of a recursive enumeration, *)
    (*   but then time would be proportional to the total DAG size, *)
    (*   instead of size of reachable part. *)

    (* Allocate the "map" array, if not given: *)
    IF map = NIL OR LAST(map^) < s THEN
      WITH
        newmap = NEW(REF ARRAY OF State, s + 1),
        nm = newmap^
      DO
        IF map = NIL THEN
          FOR i := 0 TO LAST(nm) DO nm[i] := NullState END;
        ELSE
          WITH m = map^ DO
            FOR i := 0 TO LAST(m) DO nm[i] := m[i] END;
            FOR i := NUMBER(m) TO LAST(nm) DO nm[i] := NullState END;
          END
        END;
        map := newmap
      END
    END;

    (* Now copy all uncopied nodes: *)
    WITH m = map^ DO

      PROCEDURE DoCopy(t: State): State RAISES {Full} =
      (*
        Does Copy starting at "t", a generic sucessor
        of "s": copies Last(t).dest, then copies Rest(t),
        and finally creates the copy of "t", saving it in "m".
        *)

        BEGIN
          IF t = NullState OR m[t] # NullState THEN RETURN m[t] END;
          WITH
            a = from.Last(t),
            r = from.Rest(t),
            nt = to.Append(
              last := Arc{rd := a.rd, wr := a.wr, dest := DoCopy(a.dest)},
              rest := DoCopy(r)
            )
          DO
            m[t] := nt;
            RETURN nt
          END
        END DoCopy;

      BEGIN
        RETURN DoCopy(s)
      END
    END
  END Copy;

PROCEDURE Crunch(dag: T; VAR (*IO*) root: ARRAY OF State) =
  VAR maxOldState: State;
      maxNewState: State;
  BEGIN
    (*
      Crunch is done in in four passes:

      First pass: unmark all states ("mark[s] := FALSE").

      Second pass: mark all roots as reachable ("mark[s] := TRUE"),
      then scan the states from high to low, mark the descendants of
      every marked state found.  (This works because the dag is
      acyclic and the states are topologically sorted).

      Third pass: scan all entries from 1 up, and move each marked
      node "s" to the lowest unused position, updating its "dest" and
      "rest" pointers, and saving that position in "map[s]".  Also
      update the root pointers.

      Fourth pass: rebuild the hash table for the compacted nodes.

      We actually store the marks "mark[s]" and the pointers "map[s]"
      in the "next" field of entry "s", since the hash table is not
      needed until step four.
    *)

    WITH e = dag.e^ DO

      (* Unmark all nodes *)
      FOR i := 0 TO dag.maxState DO e[i].next := NullState END;

      (* Mark root nodes, and remember highest one: *)
      (* (convention: a node is marked iff it has "next[t]=t".) *)
      maxOldState := 0;
      FOR i := 0 TO LAST(root) DO
        WITH t = root[i] DO
          IF t # NullState THEN
            maxOldState := MAX(maxOldState, t);
            e[t].next := t
          END
        END
      END;

      (* Scan from "maxOldState" down to "1", marking the descendants 
         of marked nodes: *)
      FOR t := maxOldState TO 1 BY -1 DO
        WITH et = e[t] DO
          IF et.next = t THEN
            (* "t" is reachable *)
            WITH d = et.dest + 0 DO
              <* ASSERT d < t *>
              IF d # NullState THEN e[d].next := d END
            END;
            WITH r = et.rest + 0 DO
              <* ASSERT r < t *>
              IF r # NullState THEN e[r].next := r END
            END;
          END
        END;
      END;

      (* Move every reachable state to the lowest free place
         "map[t]" (in "next[t]") *)
      e[NullState].next := NullState;
      maxNewState := NullState;
      FOR t := 1 TO maxOldState DO
        WITH
          et = e[t]
        DO
          IF et.next = t THEN
            WITH
              dest = et.dest + 0,  (* "+ 0" is a hack to force "dest" to be value, *)
              rest = et.rest + 0,  (* ("dest"/"rest" are packed, can't be designators) *)
              newdest = e[dest].next + 0,
              newrest = e[rest].next + 0,
              newt = maxNewState + 1,
              enewt = e[newt]
            DO
              <* ASSERT newt <= t *>
              <* ASSERT newdest <= dest *>
              <* ASSERT newrest <= rest *>
              enewt.dest := newdest;
              enewt.rest := newrest;
              enewt.rd := et.rd;
              enewt.wr := et.wr;
              et.next := newt;
              maxNewState := newt;
            END
          END
        END
      END;

      (* Update the root pointers and dag.maxState: *)
      FOR i := 0 TO LAST(root) DO root[i] := e[root[i]].next END;
      dag.maxState := maxNewState;

      (* Rehash nodes *)
      Rehash(dag)
    END
  END Crunch;

CONST
  DumpHeader = "DAG.Dump (format of 91-12-16)";

PROCEDURE Dump(wr: Wr.T; dag: T) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(wr, "Begin " & DumpHeader); Wr.PutChar(wr, '\n');
    Wr.PutText(wr, "max state = ");
    WrNat(wr, dag.maxState);
    Wr.PutChar(wr, '\n');
    WITH e = dag.e^ DO
      FOR i := 1 TO dag.maxState DO
        WITH ei = e[i] DO
          WrNat(wr, i);
          Wr.PutChar(wr, ' ');
          WrNat(wr, ei.rd);
          Wr.PutChar(wr, ' ');
          WrNat(wr, ei.wr);
          Wr.PutChar(wr, ' ');
          WrNat(wr, ei.dest);
          Wr.PutChar(wr, ' ');
          WrNat(wr, ei.rest);
          Wr.PutChar(wr, '\n');
        END;
      END;
    END;
    Wr.PutText(wr, "End " & DumpHeader); Wr.PutChar(wr, '\n');
    Wr.Flush(wr);
  END Dump;

PROCEDURE Load(rd: Rd.T; minSize: NAT := 0): T =
  <* FATAL Rd.Failure, Rd.EndOfFile, Thread.Alerted *>
  BEGIN
    WITH hdr = Rd.GetLine(rd) DO
      <* ASSERT Text.Equal(hdr, "Begin " & DumpHeader) *>
    END;
    WITH
      maxState = ReadParam(rd, "max state = "),
      dag = New(MAX(maxState, minSize)),
      e = dag.e^
    DO
      FOR i := 1 TO maxState DO
        WITH
          ii = RdNat(rd),
          rdLet = RdNat(rd),
          wrLet = RdNat(rd),
          dest = RdNat(rd),
          rest = RdNat(rd)
        DO
          <* ASSERT ii = i *>
          <* ASSERT rest < i *>
          <* ASSERT dest < i *>
          e[i] := Entry{
            rd := rdLet,
            wr := wrLet,
            dest := dest,
            rest := rest,
            next := 0
          };
          WITH c = Rd.GetChar(rd)DO
            <* ASSERT c = '\n' *>
          END;
        END;
      END;
      WITH hdr = Rd.GetLine(rd) DO
        <* ASSERT Text.Equal(hdr, "End " & DumpHeader) *>
      END;
      dag.maxState := maxState;
      Rehash(dag);
      RETURN dag
    END;
  END Load;

BEGIN
END DAG.

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
