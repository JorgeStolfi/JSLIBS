

FROM RDAG.i3

  It follows that an "RDAG.T" is reduced, in the following sense.  For
  any arc "e" of the "RDAG.T" with origin "s", define its "extended
  label" as being the pair "(c, i, x)", where "c = Class(s)" "x" is
  the label of "e", and "i" is the position of "e" in the list of arcs
  out of "s".  Then any directed path in the DAG from a fixed vertex
  can be identified with the "extended string" spelled by that path,
  namely the sequence of the extended labels of its arcs.  Let
  "XStrings(s)" be the set of extended strings spelled by all paths
  starting from node "s".  The unique representation invariant
  mentioned above implies that "XStrings(s) = XStrings(t)" iff "s = t".
  
  NULL NODE: Every "RDAG.T" has a distinguished node, the
  "NullNode", which has no outgoing arcs.

  INTRINSIC ARC ORDER: for every node "s" that has outgoing arcs, there is
  one distinguished arc "Last(s)" out of "s", and another node
  "Rest(s)", such that the set of arcs leaving the latter are those
  leaving "s", minus the arc "Last(s)".
  
  The methods "Last" and "Rest" define the ``intrinsic ordering'' for the
  arcs out "s".  Note that "Rest" moves {\em backwards} in this ordering.
    
  SUBNODES: by definition, a "subnode" of a node "s" of the DAG is a node
  "t" whose transitions are an initial prefix of the transitions
  out of "s", in intrinsic order (or, in other words, that can be obtained
  from "s" by zero or more applications of the "Rest" method).
  In particular, the "NullNode" is a subnode of every node.

  If a "DAG.T" includes a node "s", it automatically includes all
  subnodes of "s".

  TOPOLOGICAL ORDERING OF THE NODES: The nodes in a "DAG.T" are
  identified by numbers, assigned as they are created.  The numerical
  order of the nodes is compatible with the DAG arcs and the
  subnode relation; that is, every node is greater (numerically)
  than any of its proper subnodes, and of the destination of any of
  its outgoing arcs.

  STANDARD PATH ORDERS: The intrinsic ordering of the arcs out of each node
  induces a ``standard'' ordering of the set of all paths out of a fixed node "s".
  In this ordering, all paths that leave "s" through the arc "Last(s)" occur after 
  those that leave through other arcs. 
  
  There are actually three standard path orderings, differing on where
  the paths that {\em end} at "s" appear: either before all the paths
  that {\em go through} "s" (``pre-order''), or after them
  (``end-order''), or both before and after them (``double order'').



FROM Basics.i3:

PROCEDURE WrNat(wr: Wr.T; v: NAT);
(*
  Writes "v" to "wr", as one or more decimal digits.
  Equivalent to Wr.PutText(wr, Fmt.Int(v)), but allocates no storage. *)

PROCEDURE RdNat(rd: Rd.T): INT;
(*
  Reads a NAT (one of more decimal digits) from "rd",
  until the first non-digit, ignoring leading spaces
  (but not newlines, tabs, etc.). *)

PROCEDURE ReadParam(rd: Rd.T; label: TEXT): NAT;
(*
  Reads a line with the given label follewed by a NAT and a newline,
  returns the NAT. *)

PROCEDURE ReadParam(rd: Rd.T; label: TEXT): NAT =
  <* FATAL Rd.EndOfFile, Rd.Failure, Thread.Alerted *>
  BEGIN
    FOR i := 0 TO Text.Length(label) - 1 DO
      WITH c = Rd.GetChar(rd) DO
        <* ASSERT c = Text.GetChar(label, i) *>
      END;
    END;
    WITH i = RdNat(rd), c = Rd.GetChar(rd) DO
      <* ASSERT c = '\n' *>
      RETURN i
    END;
  END ReadParam;

CONST
  MaxExp10 = 9;  (* Max "e" such that 10**e <= LAST(NAT) *)
  P10 = ARRAY [0..MaxExp10] OF NAT{
    1,
    10,
    100,
    1000,
    10000,
    100000,
    1000000,
    10000000,
    100000000,
    1000000000
  };
  Digits = SET OF CHAR {'0' .. '9'};

PROCEDURE RdNat(rd: Rd.T): INT =
  VAR v: NAT := 0;
      c: CHAR;
  <* FATAL Rd.EndOfFile, Rd.Failure, Thread.Alerted *>
  BEGIN
    REPEAT
      c := Rd.GetChar(rd);
    UNTIL c # ' ';
    <* ASSERT c IN Digits *>
    WHILE c = '0' DO c := Rd.GetChar(rd) END;
    WHILE c IN Digits DO
      v := 10 * v + (ORD(c) - ORD('0'));
      c := Rd.GetChar(rd);
    END;
    Rd.UnGetChar(rd);
    RETURN v
  END RdNat;

PROCEDURE WrNat(wr: Wr.T; v: NAT) =
  VAR e: INT := MaxExp10;
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    IF v = 0 THEN
      Wr.PutChar(wr, '0');
      RETURN
    END;
    WHILE P10[e] > v DO e := e-1 END;
    WHILE e >= 0 DO
      WITH d = P10[e], q = v DIV d DO
        Wr.PutChar(wr, VAL(ORD('0') + q, CHAR));
        v := v - q * d;
        e := e - 1
      END
    END
  END WrNat;

PROCEDURE Copy(
    dag: T;
    from: BinGraph.T;
    s: Node;
    VAR (*IO*) map: ARRAY OF Node;
  ): Node RAISES {Full} =
  BEGIN
    WITH
      t = BinGraph.T.Copy(dag, from, s, map)
    DO
      CheckInvariants(dag, 0, dag.nNodes);
      RETURN t
    END
  END Copy;

PROCEDURE Crunch(
    dag: T; 
    READONLY root: ARRAY OF Node; 
    VAR (*OUT*) map: ARRAY OF Node;
  ) =
  BEGIN
    NARROW(dag, BinGraph.T).Crunch(dag, root, map)
  END Crunch;



FROM DAG.SET:

PROCEDURE Set(dag: T; s: State; r, w: Symbol; r, l: State) RAISES {Full} =
  VAR oldClass: Symbol; oldNStates: NAT;
  BEGIN
    IF dag.e = NIL OR s > LAST(dag.e^) THEN RAISE Full END;
    oldNStates := dag.nStates;
    IF s < oldNStates THEN
      oldClass := dag.e[s].w
    END;
    BinGraph.T.Set(dag, s, r, w, r, l);
    WITH
      e = dag.e^,
      es = e[s]
    DO
      IF l = s THEN
        (* Must be a sink: *)
        <* ASSERT r = s *>
        <* ASSERT r = 0 *>
      ELSE
        (* Must be topologically sorted,and class-homogeneous: *)
        <* ASSERT r < s *>
        <* ASSERT l < s *>
        <* ASSERT e[l].w = w *>
      END;
      IF s < oldNStates-1 AND oldClass # w THEN 
        CheckInvariants(dag, s+1, oldNStates-1)
      END;
    END
  END SET;
  


FROM BinGraph.Copy
    (* Allocate the "map" array, if not given: *)
    IF map = NIL OR LAST(map^) < s THEN
      WITH
        newmap = NEW(REF ARRAY OF State, s + 1),
        nm = newmap^
      DO
        IF map = NIL THEN
          FOR i := 0 TO LAST(nm) DO nm[i] := NoState END;
        ELSE
          WITH m = map^ DO
            FOR i := 0 TO LAST(m) DO nm[i] := m[i] END;
            FOR i := NUMBER(m) TO LAST(nm) DO nm[i] := NoState END;
          END
        END;
        map := newmap
      END
    END;

    (* Now copy all uncopied nodes: *)


FROM ReducedPair.m3:

(*

CONST
  NullLetter = FIRST(Basics.Symbol);

PROCEDURE LastLetter(aut: Reduced.T; s: Reduced.State): Basics.Symbol =
(*
  The aut.Last(s).symbol, or NullLetter if not aut.HasArcs(s). *)
  BEGIN
    IF s = RedUnit OR s = RedNull THEN
      RETURN NullLetter
    ELSE
      RETURN aut.Last(s).symbol
    END
  END LastLetter;

*)

PROCEDURE SetRep(t:Ctypes.long): Time =
  CONST
    sec = 60; (* ticks *)
    min = 60*sec;
    ho = 60*min;
    da = 24*ho;
  VAR
    res: UtilTime;
  BEGIN
    res.day    := t DIV da;
    res.hour   := t MOD da  DIV ho;
    res.minute := t MOD ho  DIV min;
    res.second := t MOD min DIV sec;
    RETURN res
  END SetRep;

    MonthName = ARRAY Date.Month OF TEXT {
      "Jan","Feb","Mar","Apr","May","Jun",
      "Jul","Aug","Sep","Oct","Nov","Dec"
    };
