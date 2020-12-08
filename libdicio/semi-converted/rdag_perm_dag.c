#ifndef _H
#define _H


IMPORT Code, DAG, Huffman, Util;
IMPORT OldFileFmt, FPut, NPut, FGet, NGet;
IMPORT FileRd, FileWr, Fmt, Time, Date, Rd, Text, Wr, Thread, OSError;

FROM Basics IMPORT Full;

FROM ReadOnlyPermDAG IMPORT 
  ComprHeader, ComprTrailer, DefaultCommentText, ExtendedLetter,
  FinalBitChar, InexistentLetter, Letter, State;

FROM Stdio IMPORT stderr;

<* FATAL Thread.Alerted, OSError.E, Rd.Failure, Wr.Failure );

CONST
  PermDAGFormat : unsigned ==  4;
  BackupPeriod  : Time.T == 3600.0d0; /* 3600 sec == 1 hour */

TYPE
  Arc                   == DAG.Arc;

REVEAL
  T == Public BRANDED OBJECT
      char *comment            ;
      unsigned format             ;
      char *textperm           ;
      charperm           : REF ARRAY OF CHAR;
      chartolettertable  : ARRAY CHAR OF ExtendedLetter;
      root               : REF ARRAY OF State;
      dag                : DAG.T;
    OVERRIDES
      Accepts            = Accepts;
      AddSub             = AddSub;
      Conta              = Conta;
      Comment            = Comment;
      Crunch             = Crunch;
      Dump               = Dump;
      DumpCompr          = DumpCompr;
      DumpFixedBin       = DumpFixedBin;
      ExpandState        = ExpandState;
      Fold               = Fold;
      IncRoots           = IncRoots;
      IsFinal            = IsFinal;
      MakeTransition     = MakeTransition;
      Prm2Red            = Prm2Red;
      Root               = Root;
      SortPrm            = SortPrm;
      Spell              = Spell;
      Statistics         = Statistics;
      UnFold             = UnFold
    ;};

TYPE
  AddressClass           == {Nul, Seq, Jmp};

CONST
  PermDAGFileType        : char *== "PermDAG.Dump";
  PermDAGFileVersion     : char *== "92-08-24";
  ReducedFileType        : char *== "Reduced.Dump";
  ReducedFileVersion     : char *== "91-12-21";
  FormatPrefix           : char *== "Format";
  RootPrefix             : char *== "PermDag.Root LAST";
VAR
  PermPrefix             : char *= "Permutation (" &
    Text.FromChar(FinalBitChar) & " == Char.NUL)";

CONST
  CharSize               : unsigned == BITSIZE(CHAR);
VAR
  BooleanSize            : unsigned = Code.BooleanSize();
  WordSize               : unsigned = Code.WordSize();
  WordSizeSize           : unsigned = Code.WordSizeSize() - 1;/* excludes 0 */

PROCEDURE Accepts(permdag: T; s: State; word: char *): BOOLEAN ==
  {
    with (len == Text.Length(word),
         w   == word & Text.FromChar(FinalBitChar),
         dag == permdag.dag,
         t   == permdag.chartolettertable
        ){
      for (i = 0 TO len){
        with (c == t[Text.GetChar(w, i)]){
          while (s!=0)  AND  AND  (dag.Last(s).rd!=c){
            s = dag.Rest(s)
          ;}
        ;};
        if ((s == 0)){ return FALSE ;};
        s = dag.Last(s).dest
      ;};
      return s == 0
    ;}
  ;} Accepts;
  
PROCEDURE AddSub(permdag: T; add: BOOLEAN; word: char *; level : unsigned = 0)
    RAISES {Full} ==
  VAR
    unsigned ind ;
    lst : unsigned = Text.Length(word);
    wrd : REF ARRAY OF Letter = NEW(REF ARRAY OF Letter, 1 + lst);
    
  PROCEDURE VisitAdd(VAR s: State) RAISES {Full} ==
    VAR
      arc  : Arc = Arc{rd = 0, dest = 0};
      rest : State = 0;
    {
      if ((ind > lst)){ assert(s == 0 ); return ;};
      with (let == wrd[ind],
           dag == permdag.dag
          ){
        if ((s == 0)){ arc.rd = let
        }else{ arc  = dag.Last(s); rest = dag.Rest(s)
        ;};
        if ((arc.rd == let)){ INC(ind); VisitAdd(arc.dest)
        }else{ VisitAdd(rest)
        ;};
        s = dag.Append(arc, rest)
      ;};
    ;} VisitAdd;
    
  PROCEDURE VisitSub(VAR s: State) RAISES {Full} ==
    VAR
      arc  : Arc = Arc{rd = 0, dest = 0};
      rest : State = 0;
    {
      if ((s == 0)){ return ;};
      assert(ind <= lst );
      with (let == wrd[ind],
           dag == permdag.dag
          ){
        arc = dag.Last(s); rest = dag.Rest(s);
        if ((arc.rd == let)){
          INC(ind);
          VisitSub(arc.dest);
          if ((arc.dest == 0)){ s = rest; return ;};
        }else{
          VisitSub(rest)
        ;};
        s = dag.Append(arc, rest)
      ;}
    ;} VisitSub;
    
  {
    for (i = 0 TO lst - 1){
      with (c == Text.GetChar(word, i)){
        assert(c!=FinalBitChar );
        wrd[i] = CharToLetter(permdag, c)
      ;}
    ;};
    wrd[lst] = CharToLetter(permdag, FinalBitChar);
    assert(level <= LAST(permdag.root^) );
    while (1){
      TRY
        ind = 0;
        if ((add)){ VisitAdd(permdag.root[level])
        }else{ VisitSub(permdag.root[level])
        ;};
        EXIT
      EXCEPT
        Full ==> 
          permdag.dag.Expand(2 + 11 * permdag.dag.MaxAllocState() DIV 10);
          Crunch(permdag)
      ;}
    ;}
  ;} AddSub;
  
PROCEDURE BuildPermDagTables(PermDag: T) ==
  {
    with (t == PermDag.textperm){
      PermDag.charperm = NEW(REF ARRAY OF CHAR, Text.Length(t));
      Text.SetChars(PermDag.charperm^, t);
      for (i = 0 TO LAST(PermDag.charperm^)){
        PermDag.chartolettertable[PermDag.charperm[i]] = i
      ;}
    ;}
  ;} BuildPermDagTables;

PROCEDURE CharToLetter(permdag: T; c: CHAR): Letter ==
  {
    with (e == permdag.chartolettertable[c]){
      if ((e == InexistentLetter)){
        with (textperm == permdag.textperm,
             charperm == permdag.charperm
            ){
          e = Text.Length(textperm);
          textperm = textperm & Text.FromChar(c);
          charperm = NEW(REF ARRAY OF CHAR, 1 + e);
          Text.SetChars(charperm^, textperm)
        ;}
      ;};
      return e
    ;}
  ;} CharToLetter;
  
PROCEDURE ClassifyAddress(s: State; a: State): AddressClass ==
  {
    if ((a == 0    )){ return AddressClass.Nul ;};
    if ((a == s - 1)){ return AddressClass.Seq ;};
                      return AddressClass.Jmp /* 1 <= a <= s - 2 */
  ;} ClassifyAddress;
  
PROCEDURE Comment(permdag: T; comment: char *) ==
  {
    permdag.comment = comment
  ;} Comment;
  
PROCEDURE Conta(permdag: T; wr: Wr.T) ==
  TYPE
    Linha   == {i,n};
    Coluna  == {seq, pul, fim};
  CONST
    ltext   : ARRAY Linha OF char *== ARRAY Linha OF char *{"i", "n"};
    ctext   : ARRAY Coluna OF char *== ARRAY Coluna OF char *{"seq", "pul", "fim"};
  VAR
    nl      : unsigned = 0;
    nd      : unsigned = 0;
    Linha linha   ;
    Coluna coluna  ;
    matriz  : ARRAY Linha OF ARRAY Coluna OF unsigned = ARRAY Linha OF
      ARRAY Coluna OF unsigned{ARRAY Coluna OF unsigned{0, ..}, ..};
  {
    with (perm    == permdag.charperm,
         dag     == permdag.dag,
         dagsize == dag.MaxState()
        ){
      for (s = 1 TO dagsize){
        with (arc  == dag.Last(s), rest == dag.Rest(s)){
          if ((perm[arc.rd] == FinalBitChar)){
            INC(nl);                       /* @ */
            linha = Linha.i
          }else{
            INC(nl); INC(nd);              /* r d */
            if ((arc.dest == s - 1)){
              linha = Linha.i
            }else{
              linha = Linha.n
            ;}
          ;};
          if ((rest == 0)){
            INC(nl);                       /* | */
            coluna = Coluna.fim
          }else if ((rest == s - 1)){
            coluna = Coluna.seq
          }else{
            INC(nl); INC(nd);              /* > d */
            coluna = Coluna.pul
          ;}
        ;};
        INC(matriz[linha, coluna])
      ;}
    ;};
    Wr.PutText(wr, "nl == " & Fmt.Int(nl) & " nd == " & Fmt.Int(nd) & 
      " total == " & Fmt.Int(nl + 2 * nd) & "\n");
    for (l = FIRST(Linha) TO LAST(Linha)){
      for (c = FIRST(Coluna) TO LAST(Coluna)){
        Wr.PutText(wr, ltext[l] & " " & ctext[c] & " " & 
          Fmt.Int(matriz[l,c]) & "\n")
      ;}
    ;}
  ;} Conta;

PROCEDURE Copy(frompermdag, topermdag: T) ==
  {
    topermdag.format            = frompermdag.format;
    topermdag.comment           = frompermdag.comment;
    topermdag.textperm          = frompermdag.textperm;
    topermdag.charperm          = frompermdag.charperm;
    topermdag.chartolettertable = frompermdag.chartolettertable;
    topermdag.root              = frompermdag.root;
    topermdag.dag               = frompermdag.dag
  ;} Copy;
  
PROCEDURE Crunch(permdag: T) ==
  {
    permdag.dag.Crunch(permdag.root^)
  ;} Crunch;
  
PROCEDURE Dump(permdag: T; wr: Wr.T) ==
  {
    OldFileFmt.WriteHeader(wr, PermDAGFileType, PermDAGFileVersion);
    
    NPut.Int(wr, FormatPrefix, permdag.format);
    FPut.EOL(wr);
    
    OldFileFmt.WriteComment(wr, permdag.comment, '|');
    
    /* Careful - blanks are significant here! */
    Wr.PutText(wr, PermPrefix);
    Wr.PutText(wr, " == ");
    Wr.PutText(wr, permdag.textperm);
    FPut.EOL(wr);
    
    with (
      root == permdag.root,
      Last == LAST(root^)
   ){
      NPut.Int(wr, RootPrefix, Last);
      FPut.Colon(wr);
      FPut.Int(wr, root[0]); 
      for (i = 1 TO Last){ FPut.Space(wr); FPut.Int(wr, root[i]) ;}
    ;};
    FPut.EOL(wr);
    
    DAG.Dump  (wr, permdag.dag);
    
    OldFileFmt.WriteFooter(wr, PermDAGFileType, PermDAGFileVersion);
  ;} Dump;

PROCEDURE DumpCompr(permdag: T; wr: Wr.T) RAISES {Huffman.OverFlow} ==
  VAR
    code         : Code.T = Code.NewCode();
    firststate   : State = 1;
    State laststate    ;
    tail         : char *= "";
    f            : Huffman.FrequencyTable;
    h            : Huffman.T;
  {
    Wr.PutText(wr, ComprHeader); FPut.EOL(wr);
    code.InitWrite(wr);
    code.Write(WordSizeSize, permdag.format);
    with (com == permdag.comment,
          len == Text.Length(com)
         ){
      code.Write(WordSize, len);
      for (i = 0 TO len - 1){
        code.Write(CharSize, ORD(Text.GetChar(com, i)))
      ;}
    ;};
    with (dag      == permdag.dag,
         m        == dag.MaxState(),
         msize    == Code.BitSize(m),
         permlast == LAST(permdag.charperm^),
         rdsize   == Code.BitSize(permlast),
         nulrd    == permdag.chartolettertable[FinalBitChar]
       ){
      code.Write(CharSize, rdsize);
      code.Write(rdsize, permlast);
      for (c = FIRST(CHAR) TO LAST(CHAR)){
        with (t == permdag.chartolettertable[c],
             codet == t!=InexistentLetter
            ){
          code.Write(BooleanSize, ORD(codet));
          if ((codet)){
            code.Write(rdsize, t)
          ;}
        ;}
      ;};
      code.Write(WordSizeSize, msize - 1);
      with (root == permdag.root){
        code.Write(WordSizeSize, NUMBER(root^));
        for (i = 0 TO LAST(root^)){
          code.Write(msize, root[i])
        ;}
      ;};

      f = NEW(Huffman.FrequencyTable, 1 + permlast);
      for (c = 0 TO permlast){ f[c] = 0 ;};
      for (s = 1 TO m){ INC(f[dag.Last(s).rd]) ;};
      h = Huffman.Tbl2Huffman(f);
      h.Dump(code);
      
      code.Write(msize, m);
      laststate = MIN(m, 2);
      for (statesize = 0 TO msize){
        for (s = firststate TO laststate){
          with (arc      == dag.Last(s),
               rd       == 0 + arc.rd,
               dest     == arc.dest,
               class    == ClassifyAddress(s, dest),
               codedest == class == AddressClass.Jmp
              ){
            assert((rd == nulrd) == (class == AddressClass.Nul) );
            h.Compress(VAL(rd, CHAR), code); /* code.Write(rdsize, rd); */
            code.Write(BooleanSize, ORD(codedest));
            if ((codedest)){
              code.Write(statesize, dest - 1)
            ;}
          ;};
          with (rest     == dag.Rest(s),
               class    == ClassifyAddress(s, rest),
               seq      == class == AddressClass.Seq,
               coderest == class == AddressClass.Jmp
              ){
            code.Write(BooleanSize, ORD(seq));
            if ((NOT seq)){
              code.Write(BooleanSize, ORD(coderest));
              if ((coderest)){
                code.Write(statesize, rest - 1)
              ;}
            ;}
          ;}
        ;};
        if ((laststate == m)){
          code.Fin(tail);
          EXIT
        }else{
          firststate = 1 + laststate;
          laststate = MIN(m, 2 * laststate - 1)
        ;}
      ;}
    ;};
    Wr.PutText(wr, ComprTrailer); FPut.EOL(wr);
  ;} DumpCompr;
  
PROCEDURE DumpFixedBin(permdag: T; wr: Wr.T) ==
  VAR
    code         : Code.T = Code.NewCode();
  {
    Wr.PutText(wr, ComprHeader); FPut.EOL(wr);
    code.InitWrite(wr);
    code.Write(WordSizeSize, permdag.format);
    with (com == permdag.comment,
          len == Text.Length(com)
         ){
      code.Write(WordSize, len);
      for (i = 0 TO len - 1){
        code.Write(CharSize, ORD(Text.GetChar(com, i)))
      ;}
    ;};
    with (dag      == permdag.dag,
         m        == dag.MaxState(),
         msize    == Code.BitSize(m),
         permlast == LAST(permdag.charperm^),
         rdsize   == Code.BitSize(permlast)
         /* , nulrd    == permdag.chartolettertable[FinalBitChar] */
       ){
      code.Write(CharSize, rdsize);
      code.Write(rdsize, permlast);
      for (c = FIRST(CHAR) TO LAST(CHAR)){
        with (t == permdag.chartolettertable[c],
             codet == t!=InexistentLetter
            ){
          code.Write(BooleanSize, ORD(codet));
          if ((codet)){
            code.Write(rdsize, t)
          ;}
        ;}
      ;};
      code.Write(WordSizeSize, msize - 1);
      with (root == permdag.root){
        code.Write(WordSizeSize, NUMBER(root^));
        for (i = 0 TO LAST(root^)){
          code.Write(msize, root[i])
        ;}
      ;};

      code.Write(msize, m);
      for (s = 1 TO m){
        with (arc      == dag.Last(s)
            ){
          code.Write(rdsize, arc.rd);
          code.Write(msize,  arc.dest);
          code.Write(msize,  dag.Rest(s))
        ;}
      ;}
    ;};
    Wr.PutText(wr, ComprTrailer); FPut.EOL(wr);
  ;} DumpFixedBin;
  
PROCEDURE ExpandState(    permdag : T; 
                            State s       ;
                        BOOLEAN *final   ;
                        VAR chars   : REF ARRAY OF CHAR) ==
  VAR
    fullchars : ARRAY CHAR OF CHAR;
    nchars    : CHAR = FIRST(CHAR);
  {
    final = FALSE;
    with (dag  == permdag.dag,
         perm == permdag.charperm
        ){
      while (s > 0){
        with (arc == dag.Last(s),
             c   == perm[arc.rd]
            ){
          if ((c == FinalBitChar)  AND  AND  (arc.dest == 0)){
            final = TRUE
          }else{
            fullchars[nchars] = c;
            INC(nchars)
          ;}
        ;};
        s = dag.Rest(s)
      ;}
    ;};
    with (n == ORD(nchars)){
      chars  = NEW(REF ARRAY OF CHAR, n);
      chars^ = SUBARRAY(fullchars, 0, n)
    ;}
  ;} ExpandState;
  
PROCEDURE Fold(permdag     : T;
                 char *backupname  ;
                 unsigned backuplevel ;
                 restore     : BOOLEAN
                ) RAISES {Full, Huffman.OverFlow} ==
  TYPE
    ExtendedState    == [-1..LAST(State)];
    StateTable       == ARRAY OF ExtendedState;
    StateDescr       == RECORD
        outdegree    : unsigned = 0;
        descr        : REF StateTable;
      ;};
    CardTable        == ARRAY OF unsigned;
  CONST
    flag             : ExtendedState  == -1;
  VAR
    alarmtime    : Time.T;
    StateDescr bottomdescr  ;
    count        : REF StateTable;
    currenttime  : Time.T;
    unsigned dagnumber    ;
    unsigned dagsize      ;
    equiv        : REF ARRAY OF unsigned;
    extradag     : DAG.T;
    indegree     : REF ARRAY OF unsigned;
    initdescr    : REF StateTable;
    unsigned initialsize  ;
    State lasts        ;
    matches      : REF StateTable;
    unsigned permlast     ;
    unsigned permlength   ;
    rdstack      : REF CardTable;
    StateDescr statedescr   ;
    workpermdag  : T = New();
    zeroroot     : REF ARRAY OF State = NEW(REF ARRAY OF State, 0);
  
  PROCEDURE Visit(READONLY stovisit: State) ==
    VAR
      State s              ;
      unsigned imax           ;
    {
      statedescr.outdegree  = 0; statedescr.descr^  = initdescr^;
      bottomdescr.outdegree = 0; bottomdescr.descr^ = initdescr^;
      with (e       == equiv[stovisit], 
           dag     == permdag.dag,
           workdag == workpermdag.dag,
           staout  == statedescr.outdegree,
           botout  == bottomdescr.outdegree
          ){
        s = stovisit;
        REPEAT
          INC(staout);
          with (arc == dag.Last(s)){
            statedescr.descr[arc.rd] = equiv[arc.dest];
            rdstack[permlength - staout] = arc.rd
          ;};
          s = dag.Rest(s)
        UNTIL s == 0;
        SUBARRAY(rdstack^, 0, staout) = 
          SUBARRAY(rdstack^, permlength - staout, staout);
        for (t = 1 TO workdag.MaxState()){
          with (mr  == matches[workdag.Rest(t)],
               m   == matches[t],
               arc == workdag.Last(t)
              ){
            if ((mr == flag) || (statedescr.descr[arc.rd]!=arc.dest)){
              m = flag
            }else{
              m = 1 + mr;
              if ((botout < m)){
                botout = m;
                e = t;
                if ((botout == staout)){ return ;}
              ;}
            ;}
          ;}
        ;};
        if ((botout == 0)){ e = 0 ;};
        s = e;
        while (s!=0){
          with (rd == workdag.Last(s).rd + 0, d == statedescr.descr[rd]){
            bottomdescr.descr[rd] = d;
            d = flag;
            DEC(staout);
            for (i = 0 TO staout){
              if ((rdstack[i] == rd)){
                SUBARRAY(rdstack^, i, staout - i) = 
                  SUBARRAY(rdstack^, 1 + i, staout - i);
                EXIT
              ;}
            ;}
          ;};
          s = workdag.Rest(s)
        ;};
        while (staout > 1){
          count^ = initdescr^;
          for (t = 1 TO dag.MaxState()){
            with (m   == matches[t],
                 arc == dag.Last(t)
                ){
              m = matches[dag.Rest(t)];
              if ((bottomdescr.descr[arc.rd] == equiv[arc.dest])){
                INC(m)
              ;};
              if ((m == botout)  AND  AND  (indegree[t] > 1)  AND  AND  (t > stovisit)){
                s = t;
                REPEAT
                  with (arc == dag.Last(s)){
                    if ((statedescr.descr[arc.rd] == equiv[arc.dest])){
                      INC(count[arc.rd])
                    ;}
                  ;};
                  s = dag.Rest(s)
                UNTIL s == 0
              ;}
            ;}
          ;};
          imax = 0;
          for (i = 1 TO staout){
            if ((count[rdstack[i]] > count[rdstack[imax]])){ imax = i ;}
          ;};
          if ((count[rdstack[imax]] == flag)){ EXIT ;};
          with (cmax == rdstack[imax],
               d    == statedescr.descr[cmax]
              ){
            TRY
              e = workdag.Append(Arc{rd = cmax, dest = d}, e)
            EXCEPT
              Full ==> assert(FALSE );
            ;};
            bottomdescr.descr[cmax] = d;
            d = flag
          ;};
          INC(botout);
          DEC(staout);
          SUBARRAY(rdstack^, imax, staout - imax) =
            SUBARRAY(rdstack^, 1 + imax, staout - imax)
        ;};
        for (i = 0 TO staout - 1){
          with (r == rdstack[i]){
            TRY
              e = workdag.Append(Arc{rd = r, dest = statedescr.descr[r]}, e)
            EXCEPT
              Full ==> assert(FALSE );
            ;}
          ;}
        ;}
      ;}
    ;} Visit;

  {
    currenttime = Time.Now();
    if ((restore)){
      FoldRestore(permdag, initialsize, lasts, equiv, workpermdag, backupname,
              backuplevel);
      alarmtime    = currenttime + BackupPeriod
    }else{
      initialsize  = permdag.dag.MaxState();
      lasts        = 0;
      equiv        = NEW(REF ARRAY OF unsigned, 1 + initialsize);
      equiv[0]     = 0; 
      workpermdag.textperm           = permdag.textperm;
      workpermdag.charperm           = permdag.charperm;
      workpermdag.chartolettertable  = permdag.chartolettertable;
      alarmtime    = currenttime
    ;};
    dagsize      = permdag.dag.MaxState();
    dagnumber    = 1 + dagsize;
    workpermdag.dag.Expand(dagnumber);
    indegree     = NEW(REF ARRAY OF unsigned, dagnumber);
    matches      = NEW(REF StateTable, dagnumber);
    matches[0]   = 0;
    
    permlast     = LAST(permdag.charperm^);
    permlength   = 1 + permlast;
    initdescr    = NEW(REF StateTable, permlength);
    for (i = 0 TO permlast){ initdescr[i] = flag ;}; 
    bottomdescr  = StateDescr{
      outdegree = 0, descr = NEW(REF StateTable, permlength)};
    statedescr   = StateDescr{
      outdegree = 0, descr = NEW(REF StateTable, permlength)};
    count        = NEW(REF StateTable, permlength);
    rdstack      = NEW(REF CardTable,  permlength);
    REPEAT
      dagsize = permdag.dag.MaxState();
      Wr.PutText(stderr, Fmt.Int(dagsize) & "\n");
      fflush(stderr);
      indegree[0] = 0;
      for (s = 1 TO dagsize){
        indegree[s] = 0;
        INC(indegree[permdag.dag.Last(s).dest], 2);
        INC(indegree[permdag.dag.Rest(s)])
      ;};
      for (i = 0 TO LAST(permdag.root^)){
        INC(indegree[permdag.root[i]], 2)
      ;};
      for (s = 1 + lasts TO dagsize){
        equiv[s] = 0;
        if ((indegree[s] > 1)){
          Visit(s);
          currenttime = Time.Now();
          if ((currenttime >= alarmtime)){
            FoldBackUp(permdag, initialsize, s, equiv, workpermdag,
                       backupname, backuplevel);
            alarmtime = currenttime + BackupPeriod
          ;}
        ;}
      ;};
      with (root == permdag.root){
        for (i = 0 TO LAST(root^)){
          with (r == root[i]){
            r = equiv[r]
          ;}
        ;}
      ;};
      permdag.dag.Crunch(zeroroot^);
      extradag = permdag.dag; 
      permdag.dag = workpermdag.dag;
      workpermdag.dag = extradag;
      lasts = 0
    UNTIL dagsize == permdag.dag.MaxState();
    Wr.PutText(stderr, "\ninitial == " & Fmt.Int(initialsize) &
      " final == " & Fmt.Int(dagsize) & " ==> " & 
      Fmt.Int(ROUND(FLOAT(dagsize) * 100.0 / FLOAT(initialsize))) &
      "%\n");
    fflush(stderr)
  ;} Fold;
  
PROCEDURE FoldBackUp(    permdag     : T;
                           unsigned initialsize ;
                           State lasts       ;
                           equiv       : REF ARRAY OF unsigned;
                           T workpermdag ;
                           char *backupname  ;
                       VAR backuplevel : unsigned
                      ) RAISES {Huffman.OverFlow}==
  VAR
    date      : Date.T;
    code      : Code.T   = Code.NewCode();
    char *tail      ;
  {
    date = Util.GetDate();
    Wr.PutText(stderr, Util.FmtDate(date) & " backup == " &
      Fmt.Int(backuplevel) & ": " & 
      Fmt.Int(ROUND(100.0 * FLOAT(lasts) / FLOAT(permdag.dag.MaxState()))) & 
      "%");
    fflush(stderr);
    with (
      size == Code.BitSize(permdag.dag.MaxState()),
      backupwriter == FileWr.Open( 
      backupname & Fmt.Int(backuplevel)) 
   ){
      permdag.DumpCompr(backupwriter);
      code.InitWrite(backupwriter);
      code.Write(WordSize, initialsize);
      code.Write(size, lasts);
      for (s = 0 TO lasts){
        code.Write(size, equiv[s])
      ;};
      tail = "";
      code.Fin(tail);
      workpermdag.DumpCompr(backupwriter);
      fclose(backupwriter);
      Wr.PutText(stderr, ": done\n");
      fflush(stderr)
    ;};
    backuplevel = 1 - backuplevel
  ;} FoldBackUp;
  
PROCEDURE FoldRestore(    permdag     : T;
                        unsigned *initialsize ;
                        State *lasts       ;
                        VAR equiv       : REF ARRAY OF unsigned;
                            T workpermdag ;
                            char *backupname  ;
                        VAR backuplevel : unsigned
                       ) RAISES {Full} ==
  VAR
    code : Code.T   = Code.NewCode();
    char *tail ;
  <* FATAL Rd.EndOfFile );
  {
    with (backupreader == FileRd.Open(
           backupname & Fmt.Int(backuplevel))){
      Copy(InCodeCompr(backupreader), permdag);
      code.InitRead(backupreader);
      initialsize = code.Read(WordSize);
      with (size       == permdag.dag.MaxState(),
           sizelength == Code.BitSize(size)){
        lasts = code.Read(sizelength);
        equiv = NEW(REF ARRAY OF unsigned, 1 + size);
        for (s = 0 TO lasts){
          equiv[s] = code.Read(sizelength)
        ;}
      ;};
      code.Fin(tail);
      Copy(InCodeCompr(backupreader), workpermdag);
      Rd.Close(backupreader)
    ;};
    backuplevel =  1 - backuplevel
  ;} FoldRestore;
  
PROCEDURE InCodeCompr(rd: Rd.T): T 
    RAISES {Full} ==
  VAR
    permdag      : T = New();
    firststate   : State = 1;
    State laststate    ;
    State dest         ;
    State rest         ;
    
    code         : Code.T = Code.NewCode();
    Letter rdsize       ;
    unsigned msize        ;
    char *tail         ;
    h            : Huffman.T;
    unsigned rdVar        ;
  <* FATAL Rd.EndOfFile );
  {
    FGet.Match(rd, ComprHeader); FGet.EOL(rd);
    code.InitRead(rd);
    with (format == permdag.format,
         com    == permdag.comment){
      format = code.Read(WordSizeSize);
      if ((format >= 1)){
        com = "";
        for (i = 1 TO code.Read(WordSize)){
          com = com & Text.FromChar(VAL(code.Read(CharSize), CHAR))
        ;}
      }else{
        com = DefaultCommentText
      ;}
    ;};
    rdsize = code.Read(CharSize);
    with (charperm   == permdag.charperm){
      charperm = NEW(REF ARRAY OF CHAR, 1 + code.Read(rdsize));
      for (c = FIRST(CHAR) TO LAST(CHAR)){
        with (codet == VAL(code.Read(BooleanSize), BOOLEAN),
             t     == permdag.chartolettertable[c]
            ){
          if ((codet)){
            t = code.Read(rdsize);
            charperm[t] = c
          }else{
            t = InexistentLetter
          ;}
        ;}
      ;};
      permdag.textperm = Text.FromChars(charperm^)
    ;};   
    msize = 1 + code.Read(WordSizeSize);
    with (root == permdag.root){
      root = NEW(REF ARRAY OF State, code.Read(WordSizeSize));
      for (i = 0 TO LAST(root^)){
        root[i] = code.Read(msize)
      ;}
    ;};

    if ((permdag.format >= 4)){ 
      h = Huffman.Load(code, NUMBER(permdag.charperm^));
    ;};

    with (m     == code.Read(msize),
         nulrd == permdag.chartolettertable[FinalBitChar],
         dag   == permdag.dag 
        ){
      dag = DAG.New(m);
      laststate = MIN(m, 2);
      for (statesize = 0 TO msize){
        for (s = firststate TO laststate){
          if ((permdag.format >= 4 
         )){ rdVar = ORD(h.Uncompress(code))
          }else{ rdVar = code.Read(rdsize)
          ;};
          if ((VAL(code.Read(BooleanSize), BOOLEAN))){
            dest = 1 + code.Read(statesize)
          }else if ((rdVar == nulrd)){
            dest = 0
          }else{
            dest = s - 1
          ;};
          if ((VAL(code.Read(BooleanSize), BOOLEAN))){
            rest = s - 1
          }else if ((VAL(code.Read(BooleanSize), BOOLEAN))){
            rest = 1 + code.Read(statesize)
          }else{
            rest = 0
          ;};
          EVAL dag.Append(Arc{rd = rdVar, dest = dest}, rest)
        ;};
        if ((laststate == m)){
          code.Fin(tail);
          EXIT
        }else{
          firststate = 1 + laststate;
          laststate = MIN(m, 2 * laststate - 1)
        ;}
      ;}
    ;};
    if ((permdag.format >= 2)){
      if ((permdag.format >= 3)){
        tail = ""
      ;};
      
    ;};
    permdag.format = PermDAGFormat;

    /* Check trailer: */
    TRY
      with (line == Rd.GetLine(rd)){ tail = tail & line ;}
    EXCEPT 
      Rd.EndOfFile ==> /*OK*/
    ;};
    assert(Text.Equal(tail, ComprTrailer) );

    return permdag
  ;} InCodeCompr;
  
PROCEDURE IncRoots(permdag: T) ==
  {
    with (root    == permdag.root,
         n       == NUMBER(root^),
         newroot == NEW(REF ARRAY OF State, 1 + n)
        ){
      SUBARRAY(newroot^, 0, n) = root^;
      root = newroot;
      root[n] = 0
    ;}
  ;} IncRoots;
  
PROCEDURE IsFinal(permdag: T; s: State): BOOLEAN ==
  VAR
    BOOLEAN final ;
    chars : REF ARRAY OF CHAR;
  {
    ExpandState(permdag, s, final, chars);
    return final
  ;} IsFinal;
  
PROCEDURE Load(rd: Rd.T): T ==
  <* FATAL Rd.EndOfFile );
  VAR permdag: T = New();
  {
    OldFileFmt.ReadHeader(rd, PermDAGFileType, PermDAGFileVersion);
    
    permdag.format = NGet.Int(rd, FormatPrefix);
    FGet.EOL(rd);
    
    permdag.comment = OldFileFmt.ReadComment(rd, '|');
    
    /* Careful - blanks are significant here! */
    FGet.Match(rd, PermPrefix); 
    FGet.Match(rd, " == ");
    permdag.textperm = Rd.GetLine(rd);
    BuildPermDagTables(permdag);
    
    with (
      nRoots == NGet.Int(rd, RootPrefix) + 1,
      r == NEW(REF ARRAY OF State, nRoots)
   ){
      FGet.Colon(rd);
      for (i = 0 TO nRoots-1){ r[i] = FGet.Int(rd) ;};
      permdag.root = r
    ;};
    FGet.EOL(rd);
    
    permdag.dag = DAG.Load(rd, 0);
    
    OldFileFmt.ReadFooter(rd, PermDAGFileType, PermDAGFileVersion);
    permdag.format = PermDAGFormat;
    return permdag
  ;} Load;
  
PROCEDURE LoadCompr(rd: Rd.T): T RAISES {Full} ==
  {
    return InCodeCompr(rd)
  ;} LoadCompr;
  
PROCEDURE MakeTransition(permdag: T; s: State; c: CHAR): State ==
  {
    with (dag    == permdag.dag,
         letter == permdag.chartolettertable[c]
        ){
      while (s > 0){
        with (arc == dag.Last(s)){
          if ((arc.rd == letter)){
            return arc.dest
          ;}
        ;};
        s = dag.Rest(s)
      ;}
    ;};
    return s
  ;} MakeTransition;
  
PROCEDURE New(): T ==
  {
    with (permdag == NEW(T,
      format            = PermDAGFormat,
      comment           = DefaultCommentText,
      textperm          = Text.FromChar(FinalBitChar),
      charperm          = NEW(REF ARRAY OF CHAR, 1),
      chartolettertable = ARRAY CHAR OF ExtendedLetter{InexistentLetter, ..},
      root              = NEW(REF ARRAY OF State, 1),
      dag               = DAG.New(0))
   ){
      permdag.charperm[0] = FinalBitChar;
      permdag.chartolettertable [FinalBitChar] = 0;
      permdag.root[0] = 0;
      return permdag
    ;}
  ;} New;
  
EXCEPTION BadPermDAG;
  
PROCEDURE Prm2Red(permdag: T; wr: Wr.T) RAISES {Full} ==
  unsigned *letter  ;
  <* FATAL BadPermDAG );
  {
    OldFileFmt.WriteHeader(wr, ReducedFileType, ReducedFileVersion);
    
    OldFileFmt.WriteComment(wr, permdag.comment, '|');
    
    with (dag    == permdag.dag,
         perm   == permdag.charperm,
         root   == permdag.root,
         m      == permdag.dag.MaxState(),
         reddag == DAG.New(m)
        ){
      assert(LAST(root^) == 0 );
      for (s = 1 TO m){
        with (arc  == dag.Last(s),
             c    == perm[arc.rd],
             dest == arc.dest,
             rest == dag.Rest(s)
            ){
          if ((c == FinalBitChar)  AND  AND  (dest == 0)){
            if ((rest == 0)){
              letter = 0
            }else{
              RAISE BadPermDAG
            ;}
          }else{
            letter = ORD(c)
          ;};
          EVAL reddag.Append(Arc{rd = letter, dest = dest}, rest);
          if ((rest!=0)  AND  AND  (reddag.Last(rest).rd >= letter)){
            RAISE BadPermDAG
          ;};
        ;}
      ;};
      
      NPut.Int(wr, "root", root[0]); FPut.EOL(wr);
      DAG.Dump (wr, reddag);
    ;};
    OldFileFmt.WriteFooter(wr, ReducedFileType, ReducedFileVersion);
  ;} Prm2Red;
  
PROCEDURE Red2Prm(rd: Rd.T): T RAISES {Full} ==
  VAR
    permdag   : T = New();
    usedchars : SET OF CHAR = SET OF CHAR{};
    CHAR equivchar ;
  {
    OldFileFmt.ReadHeader(rd, ReducedFileType, ReducedFileVersion);
    
    permdag.comment = OldFileFmt.ReadComment(rd, '|');
    
    with (
      redRoot == NGet.Int(rd, "root"),
      r == NEW(REF ARRAY OF State, 1)
   ){
      r[0] = redRoot;
      permdag.root = r
    ;};
    FGet.EOL(rd);
    
    with (dag    == permdag.dag,
         t      == permdag.textperm,
         reddag == DAG.Load(rd),
         m      == reddag.MaxState()
        ){
      dag = DAG.New(m);
      for (s = 1 TO m){
        with (rd == reddag.Last(s).rd + 0){
          if ((rd!=0)){
            usedchars = usedchars + SET OF CHAR{VAL(rd, CHAR)}
          ;}
        ;}
      ;};
      assert(NOT (FinalBitChar IN usedchars) );
      for (c = FIRST(CHAR) TO LAST(CHAR)){
        if ((c IN usedchars)){
          t = t & Text.FromChar(c)
        ;}
      ;};
      BuildPermDagTables(permdag);
      for (s = 1 TO m){
        with (arc == reddag.Last(s),
             rd  == arc.rd + 0
            ){
          if ((rd == 0)){
            equivchar = FinalBitChar
          }else{
            equivchar = VAL(arc.rd, CHAR)
          ;};
          EVAL dag.Append(
            Arc{rd   = permdag.chartolettertable[equivchar],
                dest = arc.dest}, 
            reddag.Rest(s))
        ;}
      ;}
    ;};
    OldFileFmt.ReadFooter(rd, ReducedFileType, ReducedFileVersion);
    return permdag
  ;} Red2Prm;

PROCEDURE Root(permdag: T; level : unsigned = 0): State ==
  {
    with (root == permdag.root){
      assert(level <= LAST(root^) );
      return root[level]
    ;}
  ;} Root;
  
PROCEDURE SortPrm(permdag: T) RAISES {Full} ==
  VAR
    Present       : ARRAY CHAR OF BOOLEAN = ARRAY CHAR OF BOOLEAN {FALSE, ..};
    WorkPermDag   : T = New();
  {
    with (dag      == permdag.dag,
         m        == dag.MaxState(),
         charperm == permdag.charperm,
         workdag  == WorkPermDag.dag
   ){
      WorkPermDag.comment = permdag.comment;
      WorkPermDag.root    = NEW(REF ARRAY OF State, NUMBER(permdag.root^));
      WorkPermDag.root^   = permdag.root^;
      workdag.Expand(1 + m);
      for (s = 1 TO m){
        Present[charperm[dag.Last(s).rd]] = TRUE
      ;};
      Present[FinalBitChar] = FALSE;
      with (t == WorkPermDag.textperm){
        for (c = FIRST(CHAR) TO LAST(CHAR)){
          if (( Present[c])){ t = t & Text.FromChar(c) ;}
        ;}
      ;};
      BuildPermDagTables(WorkPermDag);
      with (table == WorkPermDag.chartolettertable){
        TRY
          for (s = 1 TO m){
            with (arc == dag.Last(s), r == arc.rd + 0, d == arc.dest){
              EVAL workdag.Append(Arc{rd = table[charperm[r]], dest = d}, 
                                  dag.Rest(s))
            ;}
          ;}
        EXCEPT
          Full ==> assert(FALSE );
        ;}
      ;}
    ;};
    Copy(WorkPermDag, permdag)
  ;} SortPrm;
  
PROCEDURE Spell(permdag: T; wr: Wr.T; level: unsigned = 0) ==
  VAR
    Flag          : State = 1 + permdag.dag.MaxState();
    Word          : ARRAY [1..100] OF CHAR;
    WordLength    : unsigned = 0;

  PROCEDURE Visit(s: State) ==
    VAR
      StateDescr : ARRAY CHAR OF State = ARRAY CHAR OF State{Flag, ..};
    {
      with (dag == permdag.dag,
           perm == permdag.charperm
          ){
        while (s!=0){
          with (arc == dag.Last(s), c == perm[arc.rd]){
            if ((c == FinalBitChar)){
              Wr.PutString(wr, SUBARRAY(Word, 0, WordLength));
              Wr.PutChar(wr, '\n')
            }else{
              StateDescr[c] = arc.dest
            ;}
          ;};
          s = dag.Rest(s)
        ;};
        INC(WordLength);
        with (w == Word[WordLength]){
          for (c = FIRST(CHAR) TO LAST(CHAR)){
            with (dest == StateDescr[c]){
              if ((dest!=Flag)){
                w = c;
                Visit(dest)
              ;}
            ;}
          ;}
        ;};
        DEC(WordLength)
      ;}
    ;} Visit;

  {
    with (root == permdag.root){
      assert(level <= LAST(root^) );
      Visit(root[level])
    ;}
  ;} Spell;
  
PROCEDURE Statistics(permdag: T; wr: Wr.T) ==
  TYPE
    TypeCounter == RECORD
        BOOLEAN Visited  ;
        BOOLEAN IsState  ;
        BOOLEAN IsFinal  ;
        unsigned Arcs     ;
        unsigned Letters  ;
        Words    : unsigned
      ;};
  VAR
    States       : unsigned = 0;
    Finals       : unsigned = 0;
    Arcs         : unsigned = 0;
    Letters      : unsigned = 0;
    Words        : unsigned = 0;
    Counter      : REF ARRAY OF TypeCounter;
    dag          : DAG.T = permdag.dag;
    perm         : REF ARRAY OF CHAR = permdag.charperm;

  PROCEDURE Visit(s: State) ==
    {
      with (c == Counter[s]){
        if ((NOT c.Visited)){
          c.Visited = TRUE;
          with (arc == dag.Last(s),
               rd  == perm[arc.rd],
               dest == arc.dest,
               rest == dag.Rest(s)
              ){
            Visit(dest);
            Visit(rest);
            Counter[dest].IsState = TRUE;
            with (d == Counter[dest],
                 r == Counter[rest]
             ){
              c.Arcs     = r.Arcs;
              c.Letters  = d.Letters + r.Letters;
              c.Words    = d.Words + r.Words;
              if ((rd == FinalBitChar)  AND  AND  (dest == 0)){
                c.IsFinal = TRUE;
                INC(c.Words)
              }else{
                INC(c.Arcs);
                c.IsFinal = Counter[rest].IsFinal;
                INC(c.Letters, d.Words)
              ;}
            ;}
          ;}
        ;}
      ;}
    ;} Visit;

  {
    with (m    == dag.MaxState(),
         root == permdag.root
     ){
      Counter = NEW(REF ARRAY OF TypeCounter, 1 + m);
      Counter[0] = TypeCounter{Visited = TRUE, IsState = FALSE,
        IsFinal = FALSE, Arcs = 0, Letters = 0, Words   = 0};
      for (i = 1 TO m){
        Counter[i].Visited = FALSE
      ;};
      for (i = 0 TO LAST(root^)){
        with (r == root[i]){
          Counter[r].IsState = TRUE;
          Visit(r)
        ;}
      ;};
      for (i = 1 TO m){
        with (c == Counter[i]){
          if ((c.IsState)){
            INC(States);
            if ((c.IsFinal)){ INC(Finals) ;};
            INC(Arcs, c.Arcs)
          ;}
        ;}
      ;};
      for (i = 0 TO LAST(root^)){
        with (c == Counter[root[i]]){
          INC(Letters, c.Letters);
          INC(Words, c.Words)
        ;}
      ;};
      Wr.PutText(wr,
        Fmt.F("%9s %9s %9s %9s", "strings", "letters", "states", "finals") &
        Fmt.F("%9s %9s %9s", "arcs", "sub-sts", "lets/arc") &"\n" &
        Fmt.F("%9s %9s %9s %9s",
          "--------", "--------", "--------", "--------") &
        Fmt.F("%9s %9s %9s", "--------", "--------", "--------") &"\n" &
        Fmt.F("%9s %9s %9s %9s",
          Fmt.Int(Words), Fmt.Int(Letters), Fmt.Int(States), Fmt.Int(Finals)) &
        Fmt.F("%9s %9s %9s", Fmt.Int(Arcs), Fmt.Int(m),
          Fmt.Real(FLOAT(Letters) / FLOAT(Arcs),
            prec = 3, style = Fmt.Style.Fix))
        & "\n");
      fflush(wr)
    ;}
  ;} Statistics;
  
PROCEDURE UnFold(permdag: T) RAISES {Full} ==
  VAR
    dagsize : State              = permdag.dag.MaxState();
    newdag  : DAG.T              = DAG.New(dagsize);
    equiv   : REF ARRAY OF State = NEW(REF ARRAY OF State, 1 + dagsize);
    
    PROCEDURE Visit(s: State): State RAISES {Full} ==
    TYPE
      CharIndices == [ORD(FIRST(CHAR))..ORD(LAST(CHAR))];
    VAR
      UsedChars : SET OF CharIndices          = SET OF CharIndices{};
      Dests     : ARRAY CharIndices OF State;
    {
      with (e    == equiv[s],
           dag  == permdag.dag
          ){
        if ((e!=0) || (s == 0)){ return e ;};
        REPEAT
          with (arc == dag.Last(s)){
            UsedChars = UsedChars + SET OF CharIndices{arc.rd};
            Dests[arc.rd]  = arc.dest
          ;};
          s = dag.Rest(s)
        UNTIL s == 0;
        for (c = FIRST(CharIndices) TO LAST(CharIndices)){
          if ((c IN UsedChars)){
            while (1){
              TRY
                e = newdag.Append(Arc{rd = c, dest = Visit(Dests[c])}, e);
                EXIT
              EXCEPT
                Full ==> newdag.Expand(2 + 11 * newdag.MaxAllocState() DIV 10)
              ;}
            ;}
          ;}
        ;};
        return e
      ;}
    ;} Visit;
  
  {
    for (s = 0 TO dagsize){
      equiv[s] = 0
    ;};
    for (i = 0 TO LAST(permdag.root^)){
      with (r == permdag.root[i]){
        r = Visit(r)
      ;}
    ;};
    permdag.dag = newdag
  ;} UnFold;

{
;} PermDAG.
