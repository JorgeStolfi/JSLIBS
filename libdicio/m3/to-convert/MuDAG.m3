PROCEDURE Sink(dag: T; class: Symbol): Node RAISES {Full} =
  BEGIN
    IF dag.e = NIL OR dag.nNodes >= NUMBER(dag.e^) THEN RAISE Full END;
    WITH 
      e = dag.e^, 
      s = dag.nNodes + 0
    DO
      INC(dag.nNodes);
      WITH es = e[s] DO
        es.rSym := 0;
        es.lSym := class;
        es.r := s;
        es.l := s
      END;
      RETURN s
    END
  END Sink;

PROCEDURE Append(dag: T; rest: Node; label: Symbol; dest: Node): Node RAISES {Full} =
  BEGIN
    IF dag.e = NIL OR dag.nNodes >= NUMBER(dag.e^) THEN RAISE Full END;
    WITH 
      e = dag.e^, 
      s = dag.nNodes + 0
    DO
      <* ASSERT dest < s *>
      <* ASSERT rest < s *>
      INC(dag.nNodes);
      WITH es = e[s] DO
        es.rSym := label;
        es.r := dest;
        es.l := rest;
        es.lSym := e[rest].lSym
      END;
      RETURN s
    END
  END Append;

