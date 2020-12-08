PROCEDURE Sink(dag: T; class: Symbol): Node RAISES {Full} ==
  {
    if ((dag.e == NULL) || (dag.nNodes >= NUMBER(dag.e^))){ RAISE Full ;};
    with (
      e == dag.e^, 
      s == dag.nNodes + 0
   ){
      INC(dag.nNodes);
      with (es == e[s]){
        es.rSym = 0;
        es.lSym = class;
        es.r = s;
        es.l = s
      ;};
      return s
    ;}
  ;} Sink;

PROCEDURE Append(dag: T; rest: Node; label: Symbol; dest: Node): Node RAISES {Full} ==
  {
    if ((dag.e == NULL) || (dag.nNodes >= NUMBER(dag.e^))){ RAISE Full ;};
    with (
      e == dag.e^, 
      s == dag.nNodes + 0
   ){
      assert(dest < s );
      assert(rest < s );
      INC(dag.nNodes);
      with (es == e[s]){
        es.rSym = label;
        es.r = dest;
        es.l = rest;
        es.lSym = e[rest].lSym
      ;};
      return s
    ;}
  ;} Append;

