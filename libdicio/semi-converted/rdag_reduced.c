#ifndef _H
#define _H


/* See the copyright and disclaimer note at the end of this file. */

IMPORT Rd, Wr, TextWr, TextRd, Fmt, Text, Thread;
IMPORT Basics, DAG, StringPrinter, Encoding, PlainEncoding;
FROM Basics IMPORT NAT, POS, BOOL;
FROM Basics IMPORT Done, Full, Skip, Abort;
FROM Basics IMPORT WrNat, ReadParam;

REVEAL
  T == Public BRANDED OBJECT

      dag: DAG.T;
      /*
        The DAG.T that is being used to implement this
        reduced automaton ("aut").

        Each state of "aut" is represented by a state
        of "dag", whose out-arcs all have distinct "rd"-labels
        and zero "wr"-labels, and are sorted by decreasing "rd"-label.
        A state of "aut" is final iff the corresponding state
        of "dag" has an outgoing arc with "rd==NullLetter"; these
        special arcs all lead to the NullState of "dag".

        These conventions, plus the "state uniqueness" invariant
        of a DAG.T, ensure (I hope) that the automaton "aut" is reduced
        in the usual sense. */

      State root;
      /*
        The current root state. May be NullState. */

      /*************************************************************************/
      /* PREFIX/SUFFIX DATA                                                    */
      /*************************************************************************/

      /*
        The folowing fields are auxiliary tables used by methods that deal
        with prefixes and suffixes.

        The "nSuffs" "nSufLetters" arrays are always allocated, big enough, 
        and up-to-date.

        The prefix-related fields ("rdag", "rev", "dir", "nPrefs", and 
        "nPrefLetters") are (re)computed only when needed (and if "aut.root" 
        is not NullState). If "prefRoot" is NullState, then those fields are
        garbage; if "prefRoot" is a proper state, those fields are valid
        provided that the current root is "prefRoot".

        The reverse DAG "rdag" has one arc from rev[t] to rev[s]
        whenever the automaton has a proper arc from "s" to "t".
        The virtual arcs of "aut" that go to NullState, including
        the invisible arcs labelled with NullLetter that denote final states
        of "aut", have no correspondents in "rdag". On the other hand, "rdag"
        has one arc labelled with "NullLetter" going from rev[prefRoot]
        to NullState. */

      nSuffs: REF ARRAY OF NAT;       /* Number of suffix strings per state */
      nSuffLetters: REF ARRAY OF NAT; /* Total length of suffix strings per state */

      State prefRoot;     /* Root used to compute "rdag, rev, dir, nPrefs" */

      rdag: DAG.T;                    /* The DAG with reverse transitions */
      rev: REF ARRAY OF DAG.State;    /* Maps states of self to states of "rdag" */
      dir: REF ARRAY OF State;        /* Inverse of "rev", or NullState if undef */
      nPrefs: REF ARRAY OF NAT;       /* Number of prefix strings for each state */
      nPrefLetters: REF ARRAY OF NAT; /* Total length of prefix strings per state */

    OVERRIDES

      MaxState = MaxState;
      Final = Final;
      HasArcs = HasArcs;
      Last = Last;
      Append = Append;
      Rest = Rest;
      Root = Root;
      SetRoot = SetRoot;

      OutDeg = OutDeg;
      InDeg = InDeg;
      First = First;
      R = R;
      SetArc = SetArc;
      SetFinal = SetFinal;
      Walk = Walk;
      Accepts = Accepts;
      Rank = Rank;
      AddString = AddString;
      SubString = SubString;
      AddStringMaybeCrunch = AddStringMaybeCrunch;
      SubStringMaybeCrunch = SubStringMaybeCrunch;

      EnumOutArcs = EnumOutArcs;
      EnumInArcs = EnumInArcs;
      EnumPaths = EnumPaths;
      EnumStrings = EnumStrings;
      EnumStates = EnumStates;

      NPrefs = NPrefs;
      NSuffs = NSuffs;
      NPrefLetters = NPrefLetters;
      NSuffLetters = NSuffLetters;
      EnumPrefs = EnumPrefs;
      EnumSuffs = EnumSuffs;
      FirstPrefix = FirstPrefix;
      FirstSuffix = FirstSuffix;
      FullLabel = FullLabel;
      PrintPrefs = PrintPrefs;
      PrintSuffs = PrintSuffs;
      
      Build = Build;

      NStates = NStates;
      NArcs = NArcs;
      Count = Count;

      MaxAllocState = MaxAllocState;
      Expand = Expand;
      Crunch = Crunch;
    ;};

CONST
  NullLetter == FIRST(Basics.Symbol);
  /*
    A transition from "s" through NullLetter to NullState is used
    to indicate that "s" is final. */

PROCEDURE New(size: POS = 1): T ==
  {
    with (
      dag == DAG.New(size)
   ){
      return FromDAG(dag, doc = "", root = NullState)
    ;}
  ;} New;

PROCEDURE Copy(
    T aut;
    size: POS = 1;
    VAR /*OUT*/ map: REF ARRAY OF State;
  ): T RAISES {Full} ==
  {
    with (
      newAut == New(size),
      newRoot == CopyStates(
        from = aut,
        to = newAut,
        s = aut.Root(),
        map = /*OUT*/ map
      )
   ){
      newAut.SetRoot(newRoot);
      return newAut
    ;}
  ;} Copy;

PROCEDURE CopyStates(
    T from;
    T to;
    State s;
    VAR /*IO*/ map: REF ARRAY OF State;
  ): State RAISES {Full} ==
  {
    return DAG.Copy(from.dag, to.dag, s, map)
  ;} CopyStates;

/************************/
/* PRIMITIVE METHODS    */
/************************/

PROCEDURE HasArcs(<*UNUSED); aut: T; s: State): BOOL ==
  {
    return s!=NullState)  AND  AND  (s!=UnitState
  ;} HasArcs;

PROCEDURE Last(aut: T; s: State): Arc ==
  {
    assert(s!=NullState)  AND  AND  (s!=UnitState );
    with (a == aut.dag.Last(s)){
      return Arc{symbol = a.rd, dest = a.dest}
    ;}
  ;} Last;

PROCEDURE Rest(aut: T; s: State): State ==
  {
    assert(s!=NullState)  AND  AND  (s!=UnitState );
    return aut.dag.Rest(s)
  ;} Rest;

PROCEDURE Append(aut:T; s: State; symbol: Symbol; dest: State): State RAISES {Full} ==
  {
    if ((s!=NullState)){
      assert(symbol > aut.dag.Last(s).rd );
    ;};
    if ((dest == NullState)){
      return s
    }else{
      with (
        oldMax == aut.MaxState(),
        t == aut.dag.Append(DAG.Arc{rd = symbol, wr = 0, dest = dest}, s)
     ){
        if ((t > oldMax)){
          /* State is truly new, must compute aut.nSuffs[t], aut.nSuffLetters[t] */
          with (ns == aut.nSuffs^, nl == aut.nSuffLetters^){
            ns[t] = ns[dest] + ns[s];
            nl[t] = nl[dest] + ns[dest] + nl[s];
          ;};
        ;};
        return t
      ;}
    ;}
  ;} Append;

PROCEDURE Root(aut: T): State ==
  {
    return aut.root
  ;} Root;

PROCEDURE SetRoot(aut: T; s: State) ==
  {
    aut.root = s
  ;} SetRoot;

PROCEDURE MaxState(aut: T): POS ==
  {
    return aut.dag.MaxState()
  ;} MaxState;

PROCEDURE MaxAllocState(aut: T): POS ==
  {
    return aut.dag.MaxAllocState()
  ;} MaxAllocState;

PROCEDURE Expand(aut: T; newSize: NAT = 0) RAISES {Full} ==
  {
    aut.dag.Expand(newSize);
    ComputeNSuffs(aut);
  ;} Expand;

PROCEDURE Crunch(aut: T; keep: REF ARRAY OF State = NULL) ==
  NAT *nKeep;
  {
    DiscardPrefixData(aut);
    if ((keep == NULL)){ nKeep = 0 }else{ nKeep = NUMBER(keep^) ;};
    with (
      roots == NEW(REF ARRAY OF State, nKeep+1)^
   ){
      roots[0] = aut.root;
      if ((keep!=NULL)){ SUBARRAY(roots, 1, nKeep) = keep^ ;};
      aut.dag.Crunch(roots);
      aut.root = roots[0];
      if ((keep!=NULL)){ keep^ = SUBARRAY(roots, 1, nKeep) ;};
    ;};
    /* Recreate unit state, which may have been crunched out: */
    MakeUnitState(aut.dag);
    /* Recompute suffix counts: */
    ComputeNSuffs(aut);
  ;} Crunch;

/******************************/
/* DERIVED METHODS            */
/******************************/

PROCEDURE Final(aut: T; s: State): BOOL ==
  VAR t: State = s;
  {
    while (t!=NullState)  AND  AND  (t!=UnitState){ t = aut.Rest(t) ;};
    return t == UnitState
  ;} Final;

PROCEDURE R(aut: T; s: State; symbol: Symbol): State ==
  VAR t: State = s;
  {
    if ((s == NullState)){ return NullState ;};
    while (t!=NullState)  AND  AND  (t!=UnitState){
      with (a == aut.Last(t)){
        if ((a.symbol == symbol)){
          return a.dest
        }else if ((a.symbol < symbol)){
          return NullState
        }else{
          t = aut.Rest(t)
        ;}
      ;}
    ;};
    return NullState
  ;} R;

PROCEDURE OutDeg(aut: T; s: State): NAT ==
  VAR t: State = s; n: NAT = 0;
  {
    while (t!=NullState)  AND  AND  (t!=UnitState){ INC(n); t = aut.Rest(t) ;};
    return n
  ;} OutDeg;

PROCEDURE InDeg(aut: T; s: State): NAT ==
  {
    if ((aut.root == NullState) || (s >= aut.root)){
      return 0
    }else{
      ComputePrefixData(aut);
      assert(aut.nPrefs!=NULL );
      return aut.rdag.OutDeg(aut.rev[s])
    ;}
  ;} InDeg;

PROCEDURE First(aut: T; s: State): Arc ==
  State *p, t;
  {
    assert(s!=NullState)  AND  AND  (s!=UnitState );
    t = s;
    REPEAT p = t; t = aut.Rest(t) UNTIL t == NullState) || (t == UnitState;
    return aut.Last(p)
  ;} First;

PROCEDURE SetArc(aut: T; s: State; symbol: Symbol; dest: State): State RAISES {Full} ==

  PROCEDURE DoSet(t: State): State RAISES {Full} ==
    {
      if ((t == NullState) || (t == UnitState)){
        return aut.Append(t, symbol, dest)
      }else{
        with (a == aut.Last(t)){
          if ((a.symbol < symbol)){
            return aut.Append(t, symbol, dest)
          }else if ((a.symbol == symbol)){
            if ((a.dest == dest)){
              return t
            }else{
              return aut.Append(aut.Rest(t), symbol, dest)
            ;};
          }else{ /* a.symbol > symbol */
            with (
              rest == aut.Rest(t),
              newRest == DoSet(rest)
           ){
              if ((rest == newRest)){
                return t
              }else{
                return aut.Append(newRest, a.symbol, a.dest)
              ;}
            ;}
          ;}
        ;}
      ;}
    ;} DoSet;

  {
    return DoSet(s)
  ;} SetArc;

PROCEDURE SetFinal(aut: T; s: State; final: BOOL): State RAISES {Full} ==

  PROCEDURE DoSet(t: State): State RAISES {Full} ==
    {
      if ((t == NullState) || (t == UnitState)){
        if ((final)){ 
          return UnitState 
        }else{
          return NullState
        ;}
      }else{
        with (
          a == aut.Last(t),
          rest == aut.Rest(t),
          newRest == DoSet(rest)
       ){
          if ((newRest == rest)){ 
            return t
          }else{
            return aut.Append(newRest, a.symbol, a.dest)
          ;}
        ;}
      ;}
    ;} DoSet;

  {
    return DoSet(s)
  ;} SetFinal;

PROCEDURE Walk(aut: T; s: State; READONLY w: String): State ==
  VAR t: State = s;
  {
    for (i = 0 TO LAST(w)){
      assert(w[i]!=NullLetter );
      if ((t == NullState)){
        return NullState
      }else{
        t = aut.R(t, w[i])
      ;}
    ;};
    return t
  ;} Walk;

PROCEDURE Accepts(aut: T; s: State; READONLY w: String): BOOL ==
  {
    return aut.Final(aut.Walk(s, w))
  ;} Accepts;

PROCEDURE Rank(aut: T; s: State; READONLY w: String; reverse: BOOL = FALSE): NAT ==
  {
    if ((reverse)){
      return RankFromBottom(aut, s, w)
    }else{
      return RankFromTop(aut, s, w)
    ;}
  ;} Rank;

PROCEDURE RankFromTop(aut: T; s: State; READONLY w: String): NAT ==
/*
  Computes Rank(aut, s, w, reverse = FALSE).
  */
  VAR rank: NAT = 0;
      i: NAT = 0;
      t: State = s;
  {
    while (t!=NullState)  AND  AND  (i <= LAST(w)){
      /* Tally the empty string, if it is a suffix of "t": */
      if ((aut.Final(t))){ INC(rank) ;};
      with (x == w[i]){
        /* Tally all suffixes of "t" starting with symbols less than "x": */
        assert(x!=NullLetter );
        if ((x > FIRST(Symbol))){
          rank = rank + AddNSuffsOfChildren(aut, t, lo = FIRST(Symbol), hi = x-1)
        ;};
        /* Now add the rank of "w[i+1..]" in Suff(R(t, x)): */
        t = aut.R(t, x);
        INC(i);
      ;}
    ;};
    return rank
  ;} RankFromTop;

PROCEDURE RankFromBottom(aut: T; s: State; READONLY w: String): NAT ==
/*
  Computes Rank(aut, s, w, reverse = TRUE).
  */
  VAR rank: NAT = 0;
      i: NAT = 0;
      t: State = s;
  {
    while (t!=NullState)  AND  AND  (i <= LAST(w)){
      with (x == w[i]){
        /* Tally all suffixes of "t" starting with symbols greater than "x": */
        assert(x!=NullLetter );
        if ((x < LAST(Symbol))){
          rank = rank + AddNSuffsOfChildren(aut, t, lo = x+1, hi = LAST(Symbol))
        ;};
        /* Now add the reverse rank of "w[i+1..]" in Suff(R(t, x)): */
        t = aut.R(t, x);
        INC(i);
      ;}
    ;};
    if ((t!=NullState)){
      /* Tally all suffixes of "t" that are strictly grater than (): */
      rank = rank + aut.NSuffs(t);
      if ((aut.Final(t))){ DEC(rank) ;}
    ;};
    return rank
  ;} RankFromBottom;

PROCEDURE AddNSuffsOfChildren(aut: T; s: State; lo, hi: Symbol): NAT ==
/*
  Returns the sum of NSuffs(R(s, x)) for all "x" in the range [lo..hi]. */
  VAR sum: NAT = 0;
      t: State = s;
  {
    while (aut.HasArcs(t)){
      with (a == aut.Last(t)){
        if ((a.symbol < lo)){
          return sum
        }else if ((a.symbol <= hi)){
          INC(sum, aut.NSuffs(a.dest))
        ;};
      ;};
      t = aut.Rest(t)
    ;};
    return sum
  ;} AddNSuffsOfChildren;

PROCEDURE AddString(aut: T; s: State; READONLY w: String): State RAISES {Full} ==
  {
    return AddSubString(aut, s, w, TRUE)
  ;} AddString;

PROCEDURE SubString(aut: T; s: State; READONLY w: String): State RAISES {Full} ==
  {
    return AddSubString(aut, s, w, FALSE)
  ;} SubString;

PROCEDURE AddSubString(
    T aut;
    State s;
    String READONLY w;
    BOOL add;
  ): State RAISES {Full} ==
/*
  Does AddString(aut, s, w) or SubString(aut, s, w),
  depending on "add". */

  PROCEDURE DoAddSub(i: NAT; t: State): State RAISES {Full} ==
  /*
    Returns AddSubString(aut, t, SUBARRAY(w, i), add). */
    {
      if ((i == NUMBER(w))){
        return aut.SetFinal(t, add)
      }else{
        assert(w[i]!=NullLetter );
        return aut.SetArc(t, w[i], DoAddSub(i+1, aut.R(t, w[i])))
      ;}
    ;} DoAddSub;

  {
    return DoAddSub(0, s)
  ;} AddSubString;

PROCEDURE EnumOutArcs(aut: T; s: State; action: ArcAction) RAISES {Abort} ==
  VAR i: NAT = 0;

  PROCEDURE DoEnum(r: State) RAISES {Skip, Abort} ==
  /*
    Does EnumOutArcs on a given prefix "r" of the arcs out of "s".
    Raises Skip, Abort iff "action" does so. */
    {
      if ((r == UnitState) || (r == NullState)){
        /* Ok */
      }else{
        DoEnum(aut.Rest(r));
        with (a == aut.Last(r)){
          action(i, a);
          INC(i)
        ;}
      ;}
    ;} DoEnum;

  {
    TRY DoEnum(s) EXCEPT Skip ==> /* Ok */ ;};
  ;} EnumOutArcs;

PROCEDURE AddStringMaybeCrunch(
    T aut; 
    State *s;
    String READONLY w; 
    keep: REF States = NULL;
  ): State ==
  {
    return AddSubStringMaybeCrunch(aut, s, w, TRUE, keep)
  ;} AddStringMaybeCrunch;

PROCEDURE SubStringMaybeCrunch(
    T aut; 
    State *s; 
    String READONLY w; 
    keep: REF States = NULL;
  ): State ==
  {
    return AddSubStringMaybeCrunch(aut, s, w, FALSE, keep)
  ;} SubStringMaybeCrunch;

PROCEDURE AddSubStringMaybeCrunch(
    T aut; 
    State *s; 
    String READONLY w; 
    BOOL add;
    keep: REF States = NULL;
  ): State ==

  PROCEDURE CrunchAndExpandIt() ==
  /*
    Does aut.Crunch(s & keep), and perhaps aut.Expand(). */
    VAR nKeep: unsigned = 0;
    {
      if ((keep!=NULL)){ nKeep = NUMBER(keep^) ;};
      with (
        tKeep == NEW(REF States, nKeep + 1) 
     ){
        tKeep[0] = s;
        if ((keep!=NULL)){ SUBARRAY(tKeep^, 1, nKeep) = keep^ ;};
        aut.Crunch(tKeep);
        if ((keep!=NULL)){ keep^ = SUBARRAY(tKeep^, 1, nKeep) ;};
        s = tKeep[0];
      ;};
      if ((aut.MaxState() * 4 >= aut.MaxAllocState() * 3)){ 
        ExpandIt()
      ;};
    ;} CrunchAndExpandIt;

  PROCEDURE ExpandIt() ==
  /*
    Expands the automaton: */
    {
      with (
        oldSize == aut.MaxAllocState() + 1,
        newSize == ComputeNewSize(oldSize)
     ){
        <* FATAL Basics.Full );
        {
          aut.Expand(newSize);
        ;};
      ;};
    ;} ExpandIt;

  PROCEDURE ComputeNewSize(oldSize: NAT): NAT ==
  /*
    Chooses a suitable new size for automatic expansion.
    */
    {
      return
        MIN(
          MIN(
                1 + oldSize + oldSize,
            10001 + oldSize + oldSize DIV 2
          ),
          MIN(
            20001 + oldSize + oldSize DIV 3,
            50001 + oldSize + oldSize DIV 4
          )
        )
    ;} ComputeNewSize;

  {
    while (1){
      TRY
        return AddSubString(aut, s, w, add);
      EXCEPT 
        Full ==> CrunchAndExpandIt();
      ;};
    ;};
  ;} AddSubStringMaybeCrunch;

PROCEDURE EnumInArcs(aut: T; s: State; action: ArcAction) RAISES {Abort} ==
  VAR i: NAT = 0;
      rr: DAG.State;
  {
    if ((s == NullState)  AND  AND  (s >= aut.root)){
      return
    }else{
      ComputePrefixData(aut);
      rr = aut.rev[s];
      i = 0;
      TRY
        with (rdag == aut.rdag, dir == aut.dir^){
          while (rr!=DAG.NullState){
            with (a == rdag.Last(rr)){
              action(i, Arc{symbol = a.rd, dest = dir[a.dest]});
              INC(i);
              rr = rdag.Rest(rr)
            ;};
          ;};
        ;}
      EXCEPT
      | Skip ==> /* Ok */
      ;};
    ;};
  ;} EnumInArcs;

PROCEDURE EnumPaths(
    T aut;
    State s;
    enter: StateAction = NULL;
    push: PathAction = NULL;
    pop: PathAction = NULL;
    exit: StateAction = NULL;
  ) RAISES {Abort} ==

  VAR len: NAT = 0;

  PROCEDURE DoEnumPaths(t: State) RAISES {Abort} ==
  /*
    Does EnumPaths starting at "t", a generic sucessor
    of "s". Assumes "t" is not NullState. */

    BOOL *final; /* Final(t); set by "EnumRest" at the end of the recursion. */
        i: NAT = 0; /* Arc index (not including the NullLetter, if any) */

    PROCEDURE EnumRest(r: State) RAISES {Skip, Abort} ==
    /*
      Calls "enter" on "t" and enumerates a given prefix "r" of the arcs out of "t".
      Also sets the local variable "final" of DoEnumStates(t)) to Final(t).
      Raises Skip iff "enter" or "pop" raised "Skip". */

      {
        if ((r == UnitState) || (r == NullState)){
          final = (r == UnitState);
          if ((enter!=NULL)){ enter(len, t, final) ;};
        }else{
          EnumRest(aut.Rest(r));
          with (a == aut.Last(r)){
            TRY
              if ((push!=NULL)){ push(len, t, i, a) ;};
              INC(len);
              DoEnumPaths(a.dest);
              DEC(len);
            EXCEPT
              Skip ==> /* Ok */
            ;};
            if ((pop!=NULL)){ pop(len, t, i, a) ;};
            INC(i);
          ;}
        ;}
      ;} EnumRest;

    {
      assert(t!=NullState );
      TRY EnumRest(t) EXCEPT Skip ==> /* Ok */ ;};
      if ((exit!=NULL)){
        <* FATAL Skip );
        {
          exit(len, t, final) /* Shouldn't raise "Skip" */
        ;}
      ;}
    ;} DoEnumPaths;

  {
    if ((s!=NullState)){ DoEnumPaths(s) ;};
  ;} EnumPaths;

PROCEDURE EnumStrings(
    T aut;
    State s;
    enter: StringAction = NULL;
    exit: StringAction = NULL;
  ) RAISES {Abort} ==

  VAR rw: REF String = NEW(REF String, 100);

  PROCEDURE PathEnter /* : StateAction */ (
      NAT len;
      State s;
      BOOL final;
    ) RAISES {Skip, Abort} ==
    {
      if ((enter!=NULL)){ enter(SUBARRAY(rw^, 0, len), s, final) ;}
    ;} PathEnter;

  PROCEDURE PathPush /* : PathAction */ (
      NAT len;
      <*UNUSED); org: State;
      <*UNUSED); i: NAT;
      Arc arc;
    ) RAISES {} ==
    {
      Basics.ExpandString(rw, len + 1);
      rw[len] = arc.symbol;
    ;} PathPush;

  PROCEDURE PathPop /* : PathAction */ (
      NAT len;
      <*UNUSED); org: State;
      <*UNUSED); i: NAT;
      Arc arc;
    ) RAISES {} ==
    {
      assert(rw[len] == arc.symbol );
    ;} PathPop;

  PROCEDURE PathExit /* : StateAction */ (
      NAT len;
      State s;
      BOOL final;
    ) RAISES {Skip, Abort} ==
    {
      if ((exit!=NULL)){ exit(SUBARRAY(rw^, 0, len), s, final) ;}
    ;} PathExit;

  {
    aut.EnumPaths(
      s,
      enter = PathEnter,
      push = PathPush,
      pop = PathPop,
      exit = PathExit
    );
  ;} EnumStrings;

PROCEDURE EnumStates(
    T aut;
    READONLY base: ARRAY OF State;
    substates: BOOL = FALSE;
    enter: StateAction = NULL;
    exit: StateAction = NULL;
  ) RAISES {Abort} ==

  State *maxState;
  {
    /* Computes maximum reachable state: */
    maxState = NullState;
    for (i = 0 TO LAST(base)){ maxState = MAX(maxState, base[i]) ;};

    with (
      len == NEW(REF ARRAY OF NAT, maxState + 1)^
   ){

      /* Initialize path lengths: */
      for (s = 0 TO maxState){ len[s] = LAST(NAT) ;};
      for (i = 0 TO LAST(base)){ len[base[i]] = 0 ;};

      /* First pass: scan states from highest to lowest,
      propagating "len" and "enter"ing all reachable states: */

      for (s = maxState TO 1 BY -1){
        with (
          ls == len[s],
          final == aut.Final(s)
       ){
          if ((ls < LAST(NAT))){
            TRY
              if ((enter!=NULL)){ enter(ls, s, final) ;};
              if ((substates)){
                if ((s!=UnitState)){
                  with (d == aut.Last(s).dest, r == aut.Rest(s)){
                    len[d] = MIN(len[d], ls + 1);
                    len[r] = MIN(len[r], ls + 1)
                  ;}
                ;}
              }else{
                VAR t = s;
                {
                  while (t!=UnitState)  AND  AND  (t!=NullState){
                    with (d == aut.Last(t).dest){
                      len[d] = MIN(len[d], ls + 1)
                    ;};
                    t = aut.Rest(t)
                  ;}
                ;}
              ;};
            EXCEPT
            | Skip ==> /* Ignore any arcs out of "s" */
            ;}
          ;}
        ;}
      ;};

      /* Second pass: scan states from low to high, "exit"ing all
      reachable ones: */

      for (s = 1 TO maxState){
        with (
          ls == len[s],
          final == aut.Final(s)
       ){
          if ((ls < LAST(NAT))){
            TRY
              if ((exit!=NULL)){ exit(ls, s, final) ;};
            EXCEPT
            | Skip ==> assert(FALSE );
            ;};
          ;}
        ;}
      ;};
    ;}
  ;} EnumStates;

CONST
  DumpHeader == "Reduced.Dump (format of 91-12-21)";

PROCEDURE Dump(wr: Wr.T; aut: T) ==
  <* FATAL Wr.Failure, Thread.Alerted );
  {
    Wr.PutText(wr, "Begin " & DumpHeader); Wr.PutChar(wr, '\n');
    DumpDoc(wr, aut.doc);
    Wr.PutText(wr, "root == " & Fmt.Int(aut.root) & "\n");
    DAG.Dump(wr, aut.dag);
    Wr.PutText(wr, "End " & DumpHeader); Wr.PutChar(wr, '\n');
    fflush(wr);
  ;} Dump;
  
EXCEPTION  /* Syntax errors in dump file: */
  MissingFinalNewLine;
  InvalidHeader;
  InvalidFooter;

PROCEDURE Load(rd: Rd.T; minSize: POS = 1): T ==
  <* FATAL Rd.Failure, Rd.EndOfFile, Thread.Alerted );
  <* FATAL InvalidHeader, InvalidFooter );
  {
    with (hdr == Rd.GetLine(rd)){
      if ((NOT Text.Equal(hdr, "Begin " & DumpHeader))){ RAISE InvalidHeader ;};
    ;};
    with (
      doc == LoadDoc(rd),
      root == ReadParam(rd, "root == "),
      dag == DAG.Load(rd, minSize)
   ){
      with (hdr == Rd.GetLine(rd)){
        if ((NOT Text.Equal(hdr, "End " & DumpHeader))){ RAISE InvalidFooter ;}
      ;};
      return FromDAG(dag, doc, root)
    ;};
  ;} Load;
  
PROCEDURE Print(
    wr: Wr.T;
    T aut;
    e: Encoding.T;
  ) ==
  CONST FinalCode == ARRAY BOOL OF CHAR {' ', '*'};

  PROCEDURE PrintState /* : StateAction */ (
      <*UNUSED); len: NAT;
      State s;
      BOOL final;
    ) RAISES {} ==
    <* FATAL Wr.Failure, Thread.Alerted );
    {
      WrNat(wr, s);
      Wr.PutChar(wr, ' ');
      Wr.PutChar(wr, FinalCode[final]);
      Wr.PutChar(wr, '\n');
      <* FATAL Skip, Abort );
      {
        aut.EnumOutArcs(s, PrintArc);
      ;}
    ;} PrintState;

  PROCEDURE PrintArc /* : ArcAction */ (
      <*UNUSED); i: NAT;
      Arc arc;
    ) RAISES {} ==
    <* FATAL Wr.Failure, Thread.Alerted );
    {
      Wr.PutChar(wr, ' ');
      Wr.PutChar(wr, ' ');
      e.PrintLetter(wr, arc.symbol);
      Wr.PutChar(wr, ' ');
      Wr.PutChar(wr, '-');
      Wr.PutChar(wr, '>');
      Wr.PutChar(wr, ' ');
      WrNat(wr, arc.dest);
      Wr.PutChar(wr, '\n');
    ;} PrintArc;

  <* FATAL Wr.Failure, Thread.Alerted, Abort );
  {
    Wr.PutText(wr, "root == " & Fmt.Int(aut.root) & "\n");
    Wr.PutChar(wr, '\n');
    aut.EnumStates(base = ARRAY OF State{aut.root}, enter = PrintState);
    fflush(wr);
  ;} Print;

PROCEDURE NStates(aut: T; READONLY base: ARRAY OF State): NAT ==
  VAR nStates: NAT = 0;

  PROCEDURE CountState /* : StateAction */ (
      <*UNUSED); len: NAT;
      <*UNUSED); s: State;
      <*UNUSED); final: BOOL;
    ) RAISES {} ==
    {
      INC(nStates);
    ;} CountState;

  <* FATAL Abort );
  {
    aut.EnumStates(base = base, enter = CountState);
    return nStates;
  ;} NStates;

PROCEDURE NArcs(aut: T; READONLY base: ARRAY OF State): NAT ==
  VAR nArcs: NAT = 0;

  PROCEDURE CountArcs /* : StateAction */ (
      <*UNUSED); len: NAT;
      State s;
      <*UNUSED); final: BOOL;
    ) RAISES {} ==
    {
      INC(nArcs, aut.OutDeg(s));
    ;} CountArcs;

  <* FATAL Abort );
  {
    aut.EnumStates(base = base, enter = CountArcs);
    return nArcs;
  ;} NArcs;

PROCEDURE Count(aut: T; READONLY base: ARRAY OF State): Counts ==

  VAR nStates: NAT = 0;
      nArcs: NAT = 0;
      nFinals: NAT = 0;
      nSubStates: NAT = 0;
      nStrings: NAT = 0;
      nLetters: NAT = 0;

  PROCEDURE CountState /* : StateAction */ (
      <*UNUSED); len: NAT;
      State s;
      BOOL final;
    ) RAISES {} ==
    {
      INC(nStates);
      INC(nArcs, aut.OutDeg(s));
      if ((final)){ INC(nFinals) ;};
    ;} CountState;

  PROCEDURE CountSubState /* : StateAction */ (
      <*UNUSED); len: NAT;
      <*UNUSED); s: State;
      <*UNUSED); final: BOOL;
    ) RAISES {} ==
    {
      INC(nSubStates);
    ;} CountSubState;

  <* FATAL Abort );
  {
    /*
      Note that "nStates" counts only the DAG states reachable by "R"
      chains, not those reachable by "Last-Rest" chains.

      Note also that "nArcs" is the sum of the outdegre of all
      "R"-reachable states of the Reduced.T; thus, DAG arcs that
      are shared by more than one reachable state are counted
      more than once.
    */
    aut.EnumStates(base = base, enter = CountState);
    aut.EnumStates(base = base, substates = TRUE, enter = CountSubState);
    nStrings = 0;
    nLetters = 0;
    for (i = 0 TO LAST(base)){
      nStrings = nStrings + aut.NSuffs(base[i]);
      nLetters = nLetters + aut.NSuffLetters(base[i]);
    ;};
    return Counts{
      strings = nStrings,
      symbols = nLetters,
      states = nStates,
      substates = nSubStates,
      arcs = nArcs,
      finals = nFinals
    }
  ;} Count;

PROCEDURE NPrefs(aut: T; s: State): NAT ==
  {
    if ((s == NullState) || (s > aut.root)){
      return 0
    }else if ((s == aut.root)){
      return 1
    }else{
      ComputePrefixData(aut);
      assert(aut.nPrefs!=NULL );
      with (np == aut.nPrefs^){
        assert(LAST(np) >= aut.root );
        return np[s]
      ;};
    ;}
  ;} NPrefs;

PROCEDURE NSuffs(aut: T; s: State): NAT ==
  {
    if ((s == NullState)){
      return 0
    }else if ((s == UnitState)){
      return 1
    }else{
      assert(aut.nSuffs!=NULL );
      with (ns == aut.nSuffs^){
        assert(LAST(ns) >= s );
        return ns[s]
      ;};
    ;};
  ;} NSuffs;

PROCEDURE NPrefLetters(aut: T; s: State): NAT ==
  {
    if ((s == NullState) || (s > aut.root)){
      /* No paths from "aut.root" to "s":  */
      return 0
    }else if ((s == aut.root)){
      /* One path from "aut.root" to "s", of length 0: */
      return 0
    }else{
      ComputePrefixData(aut);
      assert(aut.nPrefLetters!=NULL );
      with (nl == aut.nPrefLetters^){
        assert(LAST(nl) >= aut.root );
        return nl[s]
      ;};
    ;}
  ;} NPrefLetters;

PROCEDURE NSuffLetters(aut: T; s: State): NAT ==
  {
    if ((s == NullState) || (s == UnitState)){
      return 0
    }else{
      assert(aut.nSuffLetters!=NULL );
      with (nl == aut.nSuffLetters^){
        assert(LAST(nl) >= s );
        return nl[s]
      ;};
    ;};
  ;} NSuffLetters;

PROCEDURE EnumPrefs(aut: T; s: State; action: PrefixAction) RAISES {Abort} ==
  VAR rs: REF String = NEW(REF String, 100);

  PROCEDURE PrefixPush /* : DAG.PathAction */ (
      NAT len;
      <*UNUSED); org: State;
      <*UNUSED); i: NAT;
      a: DAG.Arc
    ) RAISES {Skip, Abort} ==
    {
      if ((a.rd == NullLetter)){
        action(SUBARRAY(rs^, 0, len));
        RAISE Skip
      }else{
        Basics.ExpandString(rs, len+1);
        rs[len] = a.rd
      ;}
    ;} PrefixPush;

  {
    if ((s == NullState) || (s > aut.root)){
      return
    }else if ((s == aut.root)){
      action(String{});
    }else{
      ComputePrefixData(aut);
      aut.rdag.EnumPaths(aut.rev[s], push = PrefixPush)
    ;};
  ;} EnumPrefs;

PROCEDURE EnumSuffs(aut: T; s: State; action: SuffixAction) RAISES {Abort} ==

  PROCEDURE SuffixAction /* : StringAction */ (
       String READONLY w;
       <*UNUSED); dest: State;
       BOOL final;
    ) RAISES {Abort} ==
    {
      if ((final)){ action(w) ;}
    ;} SuffixAction;

  {
    aut.EnumStrings(s, enter = SuffixAction)
  ;} EnumSuffs;

PROCEDURE PrintPrefs(
    T aut;
    State s;
    spr: StringPrinter.T;
  ) ==

  PROCEDURE PrintPrefix /* : PrefixAction */ (
      READONLY w: String
    ) RAISES {Abort} ==
    {
      spr.PutString(w, rev = TRUE)
    ;} PrintPrefix;

  {
    TRY aut.EnumPrefs(s, action = PrintPrefix) EXCEPT Abort ==> /* Ok */ ;};
    spr.Reset();
  ;} PrintPrefs;

PROCEDURE PrintSuffs(
    T aut;
    State s;
    spr: StringPrinter.T;
  ) ==

  PROCEDURE PrintSuffix /* : SuffixAction */ (
      READONLY w: String
    ) RAISES {Abort} ==
    {
      spr.PutString(w, rev = FALSE)
    ;} PrintSuffix;

  {
    TRY aut.EnumSuffs(s, action = PrintSuffix) EXCEPT Abort ==> /* Ok */ ;};
    spr.Reset();
  ;} PrintSuffs;

PROCEDURE FirstPrefix(aut: T; s: State; e: Encoding.T): char *==

  VAR wr = TextWr.New();

  PROCEDURE FP(t: DAG.State) ==
  /*
    Prints into "wr" the string of the first path from "t"
    to an initial state (recognized by its first arc being labelled
    with NullLetter). */
    <* FATAL Wr.Failure, Thread.Alerted );
    {
      assert(t!=NullState );
      with (a == aut.rdag.First(t)){
        if ((a.rd == NullLetter)){
          return
        }else{
          FP(a.dest);
          e.PrintLetter(wr, a.rd);
        ;};
      ;};
    ;} FP;

  {
    assert(s!=NullState );
    assert(s <= aut.root );
    ComputePrefixData(aut);
    FP(aut.rev[s]);
    return TextWr.ToText(wr)
  ;} FirstPrefix;

PROCEDURE FirstSuffix(aut: T; s: State; e: Encoding.T): char *==
  VAR t: State = s;
  <* FATAL Wr.Failure, Thread.Alerted );
  {
    with (wr == TextWr.New()){
      while (1){
        assert(t!=NullState );
        with (a == aut.dag.First(t)){
          if ((a.rd == NullLetter)){ return TextWr.ToText(wr) ;};
          e.PrintLetter(wr, a.rd);
          t = a.dest
        ;}
      ;};
    ;};
  ;} FirstSuffix;

PROCEDURE FullLabel(aut: T; s: State; e: Encoding.T; sep: char *= ":"): char *==
  {
    return
      aut.FirstPrefix(s, e) & sep & aut.FirstSuffix(s, e)
  ;} FullLabel;
  
VAR /*CONST*/ DefaultEncoding: Encoding.T = PlainEncoding.New();

PROCEDURE Build(
    T aut;
    NextStringProc next;        /* Client input procedure */
    wr: Wr.T = NULL;             /* Writer for progress report, etc: */
    e: Encoding.T = NULL;        /* Symbol/CHAR encoding for printout */
    reportInterval: POS = 1000; /* Print a report every this many input strings */
    flagRedundant: BOOL = TRUE; /* TRUE to print warnings on redundant operations */
  ) RAISES {Abort} ==

  VAR
    nStrings: ARRAY BOOL OF NAT = ARRAY OF NAT {0, 0};
    nLetters: ARRAY BOOL OF NAT = ARRAY OF NAT {0, 0};

  VAR
    lastCrunchMaxState: State = NullState;  /* aut.MaxState() after last GC run */

  PROCEDURE CrunchIt() ==
  /*
    Does aut.Crunch(), printing status reports. */
    <* FATAL Wr.Failure, Thread.Alerted );
    {
      with (size == aut.MaxAllocState() + 1){
        if ((wr!=NULL)){
          if ((NOT TimeToReport())){ PrintStatusReport() ;};
          Wr.PutText(wr, "    * (crunching, alloc == " &  Fmt.Int(size) & "...");
          fflush(wr);
        ;};
        aut.Crunch(pSt);
        if ((wr!=NULL)){ 
          Wr.PutText(wr, ")\n");
          PrintStatusReport();
        ;};
      ;};
      lastCrunchMaxState = aut.MaxState();
    ;} CrunchIt;

  PROCEDURE ExpandIt() ==
  /*
    Expands the automaton, printing some noise: */
    <* FATAL Wr.Failure, Thread.Alerted );
    {
      with (
        oldSize == aut.MaxAllocState() + 1,
        newSize == ComputeNewSize(oldSize)
     ){
        if ((wr!=NULL)){
          if ((NOT TimeToReport())){ PrintStatusReport() ;};
          Wr.PutText(wr, "    * (expanding from ");
          Wr.PutText(wr, Fmt.Int(oldSize));
          Wr.PutText(wr, " to ");
          Wr.PutText(wr, Fmt.Int(newSize));
          Wr.PutText(wr, "...");
          fflush(wr);
        ;};
        <* FATAL Basics.Full );
        {
          aut.Expand(newSize);
        ;};
        if ((wr!=NULL)){ 
          Wr.PutText(wr, ")\n");
          PrintStatusReport();
        ;};
      ;};
    ;} ExpandIt;

  PROCEDURE ComputeNewSize(oldSize: NAT): NAT ==
  /*
    Chooses a suitable new size for automatic expansion.
    */
    {
      return
        MIN(
          MIN(
                1 + oldSize + oldSize,
            10001 + oldSize + oldSize DIV 2
          ),
          MIN(
            20001 + oldSize + oldSize DIV 3,
            50001 + oldSize + oldSize DIV 4
          )
        )
    ;} ComputeNewSize;

  VAR
    pN: NAT = 0;                              /* Number of pending SetArc actions */
    pSt: REF States = NEW(REF States, 100);   /* State arguments for pending actions */
    pLet: REF String = NEW(REF String, 100);  /* Symbol arguments for pending SetArcs */
      /*
        For efficiency reasons, the strings returned by "next" are not added or
        deleted right away.  Instead, we keep a list "pa[0..pN-1]" of "pending
        SetArc actions" to be performed later.
        
        The "i"th pending SetArc action is described by the pair /pSt[i],
        pLet[i]". The action consists in computing a new state "new[i]/ that 
        is like "pSt[i]" except that under the symbol
        "pLet[i]" it goes to the state "new[i+1]". As a special case,
        "new[pN]" is defined to be "pSt[pN]" itself, and "pLet[pN]" is
        not used.
        
        Thus, to flush out all pending actions we must repeat
        
|          pSt[pN-1] = aut.SetArc(pSt[pN-1], pLet[pN-1], s);
|          DEC(pN)

        until pN == 0, and finally set the automaton's root to the
        resuting state "pSt[0]". */
    
  PROCEDURE DoSetFinal(s: State; final: BOOL): State ==
  /*
    Computes aut.SetFinal(s, final), but expands "aut" if necessary
    (instead of raising "Full"). */
    {
      while (1){
        TRY
          return aut.SetFinal(s, final);
        EXCEPT 
          Full ==> ExpandIt();
        ;};
      ;};
    ;} DoSetFinal;

  PROCEDURE DoOnePendingAction() ==
  /*
    Performs the last pending "SetArc" action. Expands if necessary. 
    */
    {
      assert(pN > 0 );
      with (
        st == pSt^, let == pLet^, 
        n1 == pN-1 
     ){
        while (1){
          TRY
            st[n1] = aut.SetArc(st[n1], let[n1], st[pN]);
            pN = n1;
            return
          EXCEPT 
            Full ==> ExpandIt();
          ;};
        ;};
      ;}
    ;} DoOnePendingAction;

  PROCEDURE FlagRedundant(READONLY w: String; add: BOOL) ==
    <* FATAL Wr.Failure, Thread.Alerted );
    {
      Wr.PutText(wr, "  ** redundant command: ");
      PrintCommand(w, add);
      Wr.PutChar(wr, '\n');
    ;} FlagRedundant;

  PROCEDURE AddOrSubIt(READONLY w: String; add: BOOL = TRUE) ==
  /*
    Adds or deletes "w", crunching and/or expanding if necessary:
    */
    {
      /* Crunch the automaton, if it looks worth it: */
      with (alloc == aut.MaxAllocState()){
        if ((aut.MaxState() + 30*NUMBER(w) >= alloc
       )  AND  AND  ((aut.MaxState() - lastCrunchMaxState) >= alloc DIV 10
       )){
          /* We are close to the allocated size,
            and about 10% of the current dag states were
            created since the last CG run. Better run GC again... */
          CrunchIt();
        ;};
      ;};

      /* Now add the string, and expand if doesn't fit: */
      VAR np: unsigned = 0;
      {
        /* Skip pending actions that match the symbols of w: */
        with (maxp == MIN(pN, NUMBER(w)), let == pLet^){
          while (np < maxp)  AND  AND  (w[np] == let[np]){ INC(np) ;}
        ;};
        
        /* Flush any remaining actions: */
        if ((pN > np)){
          REPEAT DoOnePendingAction() UNTIL pN <= np
        ;};
        
        /* Stack new actions corresponding to remaining symbols of w: */
        Basics.ExpandString(pLet, NUMBER(w));
        ExpandStates(pSt, NUMBER(w) + 1);
        assert(pN == np );
        with (st == pSt^, let == pLet^){
          for (i = np TO LAST(w)){
            let[i] = w[i];
            st[i+1] = aut.R(st[i], w[i]);
          ;};
          pN = NUMBER(w);
          if ((add == aut.Final(pSt[pN]))){
            /* Command was superflous */
            if ((wr!=NULL)  AND  AND  (flagRedundant)){
              FlagRedundant(w, add)
            ;};
          }else{
            pSt[pN] = DoSetFinal(pSt[pN], add)
          ;};
        ;};
      ;};
    ;} AddOrSubIt;

  VAR
    refw: REF String = NEW(REF String, 100);  /* String buffer */
    NAT len;                                  /* Length of string */
    BOOL add;                                 /* Add/delete flag for refw^ */
    
  PROCEDURE PrintStatusReport() ==
    <* FATAL Wr.Failure, Thread.Alerted );
    {
      Wr.PutText(wr, "    * ");
      with (
        ns == nStrings[FALSE] + nStrings[TRUE],
        nl == nLetters[FALSE] + nLetters[TRUE]
     ){
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(ns), 8) &  " strings ");
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(nl), 8) &  " symbols ");
      ;};
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(aut.MaxState()), 8) & " dag states  ");
      PrintCommand(SUBARRAY(refw^, 0, len), add);
      Wr.PutText(wr, "\n");
      fflush(wr);
    ;} PrintStatusReport;

  PROCEDURE PrintCommand(READONLY w: String; add: BOOL) ==
    <* FATAL Wr.Failure, Thread.Alerted );
    {
      Wr.PutChar(wr, ARRAY BOOL OF CHAR{'-', '+'}[add]);
      Wr.PutChar(wr, ' ');
      TRY e.PrintString(wr, w) EXCEPT Encoding.BadString ==> /*IGNORE*/ ;};
    ;} PrintCommand;

  PROCEDURE PrintFinalReport() ==
    <* FATAL Wr.Failure, Thread.Alerted );
    {
      Wr.PutText(wr, "\n");

      Wr.PutText(wr, "Input statistics:\n");
      Wr.PutText(wr, "\n");

      Wr.PutText(wr, "add:  ");
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(nStrings[TRUE]), 8) & " strings ");
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(nLetters[TRUE]), 8) & " symbols ");
      Wr.PutText(wr, "\n");

      Wr.PutText(wr, "sub:  ");
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(nStrings[FALSE]), 8) & " strings ");
      Wr.PutText(wr, Fmt.Pad(Fmt.Int(nLetters[FALSE]), 8) & " symbols ");
      Wr.PutText(wr, "\n");

      Wr.PutText(wr, "------");
      Wr.PutText(wr, "--------");
      Wr.PutText(wr, "---------");
      Wr.PutText(wr, "--------");
      Wr.PutText(wr, "---------");
      Wr.PutText(wr, "\n");
      with (
        ns == nStrings[FALSE] + nStrings[TRUE],
        nl == nLetters[FALSE] + nLetters[TRUE]
     ){
        Wr.PutText(wr, "tot:  ");
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(ns), 8) & " strings ");
        Wr.PutText(wr, Fmt.Pad(Fmt.Int(nl), 8) & " symbols ");
      ;};
      Wr.PutText(wr, "\n");
      Wr.PutText(wr, "\n");
      with (ct == aut.Count(ARRAY OF State{aut.Root()})){
        PrintCounts(wr, ct)
      ;};
    ;} PrintFinalReport;

  PROCEDURE TimeToReport(): BOOL ==
    {
      return (nStrings[FALSE] + nStrings[TRUE]) MOD reportInterval == 0
    ;} TimeToReport;
  
  { /*Build*/
    pN = 0;
    pSt[0] = aut.Root();
    if ((wr!=NULL)){
      if ((e == NULL)){ e = DefaultEncoding ;};
    ;};
    TRY
      TRY
        while (1){
          next(/*IO*/ refw, /*OUT*/ len, /*OUT*/ add);
          with (w == SUBARRAY(refw^, 0, len)){
            AddOrSubIt(w, add);
            INC(nStrings[add]);
            INC(nLetters[add], NUMBER(w));
            if ((wr!=NULL)  AND  AND  (TimeToReport())){ 
              PrintStatusReport(); 
            ;};
          ;}
        ;}
      EXCEPT
      | Done ==> /* Ok */
      ;};
    FINALLY /* Normally, or in case of "Abort" */
      while (pN > 0){ DoOnePendingAction() ;};
      aut.SetRoot(pSt[0])
    ;};
    CrunchIt();
    if ((wr!=NULL)){ PrintFinalReport() ;};
  ;} Build;

/************************/
/* AUXILIARY OPERATIONS */
/************************/

PROCEDURE FromDAG(dag: DAG.T; doc: char *; root: DAG.State): T ==
/*
  Bulds a Reduced.T given the underlying DAG and the root state.
  The DAG must satisfy the implementation conventions about the
  use of NullLetter and the ordering of arc labels. */
  {
    with (
      aut == NEW(T,
        doc = doc,
        dag = dag,
        root = root,
        nSuffs = NULL,
        nSuffLetters = NULL,
        prefRoot = NullState,
        rdag = NULL,
        rev = NULL,
        dir = NULL,
        nPrefs = NULL,
        nPrefLetters = NULL
      )
   ){
      MakeUnitState(dag);
      ComputeNSuffs(aut);
      return aut;
    ;}
  ;} FromDAG;

PROCEDURE MakeUnitState(dag: DAG.T) ==
/*
  Creates the "unit" state in the "dag", if necesary: */
  <* FATAL Basics.Full );
  {
    with (
      unit == dag.Append(
        last = DAG.Arc{rd = 0, wr = 0, dest = DAG.NullState},
        rest = DAG.NullState
      )
   ){
      assert(unit == UnitState );
    ;}
  ;} MakeUnitState;

PROCEDURE DiscardPrefixData(aut: T) ==
  {
    aut.prefRoot = NullState;
    aut.rdag = NULL;
    aut.rev = NULL;
    aut.dir = NULL;
    aut.nPrefs = NULL;
    aut.nPrefLetters = NULL;
  ;} DiscardPrefixData;

PROCEDURE ComputePrefixData(aut: T) ==
/*
  If the prefix tables of "aut" ("rdag", "rev", "dir" and "nPrefs")
  are missing or out of date (e.g., because the root has changed),
  recomputes them for the current root.
  */
  {
    if ((aut.root == NullState) || (aut.prefRoot == aut.root)){
      /* Current prefix data is still OK */
      return
    }else{
      DiscardPrefixData(aut);
      ComputeReverseDAG(aut);
      ComputeNPrefs(aut);
      aut.prefRoot = aut.root;
    ;}
  ;} ComputePrefixData;

PROCEDURE ComputeReverseDAG(aut: T) ==
/*
  Computes "aut.rdag", "aut.rev", "aut.dir" from "aut.dag" and "aut.root".
  */
  VAR maxRev: DAG.State;
  <* FATAL Basics.Full );
  {
    assert(aut.root!=NullState );
    /*
      Note that the reverse DAG may be bigger than "aut.dag",
      because one arc of "aut" that is shared by two
      reachable states will give rise to two distinct reverse arcs.
      */
    with (
      totArcs == aut.Count(ARRAY OF State{aut.root}).arcs,
      rdag == DAG.New(size = totArcs + 1),
      rev == NEW(REF ARRAY OF DAG.State, aut.root+1),
      r == rev^
   ){
      for (t = 0 TO aut.root){ r[t] = NullState ;};
      /* Create one transition from the root state to NullState with NullLetter: */
      r[aut.root] = rdag.Append(
        last = DAG.Arc{rd = NullLetter, wr = NullLetter, dest = NullState},
        rest = NullState
      );

      /*
        For every state "s" reachable from "root", and every
        proper successor "t" of "s" in "aut", add to "rdag" an arc from
        "rev[t]" to "rev[s]", and update rev[t] accordingly.
        Note that it is important to process the states from high to low,
        so that rev[s] is stable by the time we add the reverse arcs into it. */

      maxRev = NullState;
      for (s = aut.root TO 1 BY -1){
        if ((r[s]!=NullState)){
          /* "s" is reachable from the root state. */

          /* Update "maxRev": */
          maxRev = MAX(maxRev, r[s]);

          /* Enumerate its outgoing arcs, and add them to the reverse dag: */
          VAR t: State = s;
          {
            while (aut.HasArcs(t)){
              with (a == aut.Last(t)){
                r[a.dest] = rdag.Append(
                  last = DAG.Arc{rd = a.symbol, wr = NullLetter, dest = r[s]},
                  rest = r[a.dest]
                );
                t = aut.Rest(t)
              ;};
            ;};
          ;};
        ;};
      ;};

      /* Store results in "aut": */
      aut.rdag = rdag;
      aut.rev = rev;

      /* Build table from reversed state to direct state: */
      with (
        dir == NEW(REF ARRAY OF State, maxRev + 1),
        d == dir^
     ){
        for (t = 0 TO maxRev){ d[t] = NullState ;};
        for (s = 0 TO aut.root){
          with (t == r[s]){
            if ((t!=NullState)){
              assert(d[t] == NullState );
              d[t] = s
            ;}
          ;}
        ;};
        aut.dir = dir
      ;};
    ;};
  ;} ComputeReverseDAG;

PROCEDURE ComputeNPrefs(aut: T) ==
/*
  Computes the tables "nPrefs" and "nPrefLetters", for the current root state,
  from 
  */
  State *t;
  {
    assert(aut.root!=NullState );
    /* Allocates the vector if necessary. */
    with (
      minSize == aut.root + 1
   ){
      if ((aut.nPrefs == NULL) || (NUMBER(aut.nPrefs^) < minSize)){
        aut.nPrefs = NEW(REF ARRAY OF NAT, minSize)
      ;};
      if ((aut.nPrefLetters == NULL) || (NUMBER(aut.nPrefLetters^) < minSize)){
        aut.nPrefLetters = NEW(REF ARRAY OF NAT, minSize)
      ;};
    ;};

    with (
      np == aut.nPrefs^,
      nl == aut.nPrefLetters^,
      maxState == aut.root
   ){
      for (s = NullState TO maxState){ 
        np[s] = 0; nl[s] = 0 
      ;};
      if ((maxState > NullState)){ np[maxState] = 1 ;};
      for (s = maxState TO 1 BY -1){
        with (nps == np[s], nls == nl[s]){
          if ((nps!=0)){
            t = s;
            while (aut.HasArcs(t) ){
              with (
                d == (aut.Last(t)).dest
             ){
                assert(d < s );
                INC(np[d], nps);
                INC(nl[d], nls + nps);
              ;};
              t = aut.Rest(t)
            ;}
          ;}
        ;};
      ;}
    ;}
  ;} ComputeNPrefs;

PROCEDURE ComputeNSuffs(aut: T) ==
/*
  Computes the tables "nSuffs" and "nSuffLetters", from "aut.dag". */
  {
    /* Allocates the vectors if necessary. */
    with (
      minSize == aut.MaxAllocState() + 1
   ){
      if ((aut.nSuffs == NULL) || (NUMBER(aut.nSuffs^) < minSize)){
        aut.nSuffs = NEW(REF ARRAY OF NAT, minSize)
      ;};
      if ((aut.nSuffLetters == NULL) || (NUMBER(aut.nSuffLetters^) < minSize)){
        aut.nSuffLetters = NEW(REF ARRAY OF NAT, minSize)
      ;};
    ;};

    /* Computes suffix counts and sizes, by a postorder scan */
    with (
      ns == aut.nSuffs^,
      nl == aut.nSuffLetters^,
      maxState == aut.MaxState()
   ){
      assert(maxState >= UnitState );
      ns[NullState] = 0;  nl[NullState] = 0;
      ns[UnitState] = 1;  nl[UnitState] = 0;
      for (s = 2 TO maxState){
        with (
          d  == (aut.Last(s)).dest,
          r  == aut.Rest(s)
       ){
          assert((d < s))  AND  AND  ((r < s) );
          ns[s] = ns[d] + ns[r];
          nl[s] = nl[d] + ns[d] + nl[r]
        ;};
      ;};
      /* Clear out remaining counts, just for tidiness: */
      for (i = maxState + 1 TO LAST(ns)){ ns[i] = 0; nl[i] = 0 ;};
    ;};
  ;} ComputeNSuffs;

CONST DocPrefix: CHAR == '|';
  
PROCEDURE DumpDoc(wr: Wr.T; doc: char *) ==
/* 
  Writes the given "doc" text to "wr", with a DocPrefix
  in front of every line.  Supplies a final '\n' if the text is 
  non-empty but does not end with newline. */
  
  VAR rd: Rd.T = TextRd.New(doc);
  
  PROCEDURE CopyLine() RAISES {Rd.EndOfFile} ==
  /*
    Copy one line from "rd" to "wr", prefixed by DocPrefix. 
    Supplies a final '\n' if next line exists but does not end with newline.
    Raises Rd.EndOfFile if there are no more lines in "rd". */
    
    <* FATAL Rd.Failure, Wr.Failure, Thread.Alerted );
    CHAR *c;
    {
      c = Rd.GetChar(rd); /* If EOF here, propagate to caller */
      Wr.PutChar(wr, DocPrefix);
      Wr.PutChar(wr, c);
      while (c!='\n'){
        TRY c = Rd.GetChar(rd) EXCEPT Rd.EndOfFile ==> c = '\n' ;};
        Wr.PutChar(wr, c)
      ;}
    ;} CopyLine;

  {
    TRY while (1){CopyLine() ;} EXCEPT Rd.EndOfFile ==> /* Ok */ ;};
  ;} DumpDoc;

PROCEDURE LoadDoc(rd: Rd.T): char *==
/*
  Reads zero or more lines from "rd" that begin with DocPrefix, strips the
  leading DocPrefix of each line, and returns the concatenation of those lines
  as a single char *with embedded and terminating newline chars. */

  VAR wr: Wr.T = TextWr.New();

  PROCEDURE CopyLine() RAISES {Rd.EndOfFile} ==
  /*
    Copy one DocPrefix line from "rd" to "wr", removing the DocPrefix
    but leaving the final (mandatory) newline.
    Raises Rd.EndOfFile if "rd" is exhausted or the next char is 
    not DocPrefix. */
    <* FATAL Rd.Failure, Wr.Failure, Thread.Alerted );
    <* FATAL MissingFinalNewLine );
    CHAR *c;
    {
      c = Rd.GetChar(rd); /* If EOF here, propagate to caller */
      if ((c!=DocPrefix)){ Rd.UnGetChar(rd); RAISE Rd.EndOfFile ;};
      REPEAT
        TRY c = Rd.GetChar(rd) EXCEPT Rd.EndOfFile ==> RAISE MissingFinalNewLine ;};
        Wr.PutChar(wr, c)
      UNTIL c == '\n'
    ;} CopyLine;

  VAR twr: TextWr.T;
      char *txt;
  {
    TRY while (1){CopyLine() ;} EXCEPT Rd.EndOfFile ==> /* Ok */ ;};
    twr = wr;
    txt = TextWr.ToText(twr);
    return txt
  ;} LoadDoc;  

PROCEDURE PrintCounts(wr: Wr.T; READONLY ct: Counts) ==
  <* FATAL Wr.Failure, Thread.Alerted );
  {
    Wr.PutText(wr, " strings  symbols    states   finals     arcs  sub-sts lets/arc\n");
    Wr.PutText(wr, "-------- --------  -------- -------- -------- -------- --------\n");
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.strings),   8));
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.symbols),   8));
    Wr.PutChar(wr, ' ');
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.states),    8));
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.finals),    8));
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.arcs),      8));
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, Fmt.Pad(Fmt.Int(ct.substates), 8));
    Wr.PutChar(wr, ' ');
    Wr.PutText(wr, 
      Fmt.Pad(Fmt.Real(FLOAT(ct.symbols)/FLOAT(ct.arcs), Fmt.Style.Fix, 3), 8)
    );
    Wr.PutChar(wr, '\n');
  ;} PrintCounts;
  
PROCEDURE ExpandStates(VAR s: REF States; n: unsigned) ==
  {
    with (nold == NUMBER(s^)){
      if ((nold < n)){
        with (
          r == NEW(REF States, MAX(n + 10, 2*nold))
       ){
          SUBARRAY(r^, 0, nold) = s^;
          s = r
        ;}
      ;}
    ;}
  ;} ExpandStates;

{
  assert(NullState == DAG.NullState );
  assert(NullLetter < FIRST(Symbol) );
;} Reduced.

/***********************************************************************/
/* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           */
/*                    Campinas, SP, Brazil                                  */
/*                                                                          */
/* Authors:                                                                 */
/*                                                                          */
/*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.br>         */
/*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.br>       */
/*   Jorge Stolfi        - CS Dept, UNICAMP <stolfi@dcc.unicamp.br>         */
/*                                                                          */
/* This file can be freely distributed, modified, and used for any          */
/*   non-commercial purpose, provided that this copyright and authorship    */
/*   notice be included in any copy or derived version of this file.        */
/*                                                                          */
/* DISCLAIMER: This software is offered ``as is'', without any guarantee    */
/*   as to fitness for any particular purpose.  Neither the copyright       */
/*   holder nor the authors or their employers can be held responsible for  */
/*   any damages that may result from its use.                              */
/****************************************************************************/
 
