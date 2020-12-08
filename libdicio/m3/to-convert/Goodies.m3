MODULE Goodies EXPORTS Main;

IMPORT Code, Huffman, PermDAG, ReadOnlyPermDAG, Util, ParamUtil; 
IMPORT FileRd, TextWr, Fmt, ParseParams, Text, Rd, Wr, OSError;
IMPORT Time, Process, Thread;

FROM Basics IMPORT Full;

FROM Stdio IMPORT stderr, stdin, stdout;

CONST
  UsageText = 
    "usage: Goodies \\\n" &
    "  [ -new \\\n" &
    "  | { -add | -sub } -wordrd FILENAME [ -level ROOTNUM ] \\\n" &
    "  | -spell [ -level ROOTNUM ] \\\n" &
    "  | -comment -comment TEXT \\\n" &
    "  | -compress | -fixedbin | -uncompress \\\n" &
    "  | -fold -backup FILENAME \\\n" &
    "  | -rfold -backup FILENAME [ -level BACKUPLEVEL ] \\\n" &
    "  | -unfold \\\n" &
    "  | -incroots \\\n" &
    "  | -prm2red | -red2prm \\\n" &
    "  | -sortprm \\\n" &
    "  | -crunch \\\n" &
    "  | -conta \\\n" &
    "  | -stats \\\n" &
    "  | -viewbits \\\n" &
    "  | -test \\\n" &
    "  ]\n";

TYPE 
  Parms = RECORD
      task: Task;             (* The operation to be performed *)
      backupfilename: TEXT;   (* For "-fold" and "-rfold" *)
      backuplevel: CARDINAL;  (* For "-rfold" *)
      rootnum: CARDINAL;      (* For "-spell", "-add", "-sub" *)
      wordrd: Rd.T;           (* For "-add" and "-sub" *)
      comment: TEXT;          (* For "-comment" *)
      cmd: TEXT;              (* Formatted command line *)
    END;

  Task = { 
      none, add, compress, comment, conta, 
      crunch, fixedbin, fold, incroots, new, 
      prm2red, red2prm, rfold, sortprm, spell,
      stats, sub, test, uncompress, unfold,
      viewbits
    };

CONST
  TaskName = ARRAY Task OF TEXT{
      "none", "add", "compress", "comment", "conta",
      "crunch", "fixedbin", "fold", "incroots", "new",
      "prm2red", "red2prm", "rfold", "sortprm", "spell",
      "stats", "sub", "test", "uncompress", "unfold",
      "viewbits"
    };
    
  Separator = "=====================================================";

<* FATAL Thread.Alerted, Rd.Failure, Wr. Failure *>

PROCEDURE Main() =
  BEGIN
    InitialMessage(stderr);
    WITH
      o = GetParms()
    DO
      Wr.PutText(stderr, "task: ");
      Wr.PutText(stderr, TaskName[o.task]);
      Wr.PutText(stderr, "\n");
      TRY
        CASE o.task OF 
        | Task.new        => New();
        | Task.add        => AddSub(TRUE,  o.wordrd, o.rootnum);
        | Task.sub        => AddSub(FALSE, o.wordrd, o.rootnum);
        | Task.spell      => Spell(o.rootnum);
        | Task.comment    => Comment(o.comment);
        | Task.compress   => Compress();
        | Task.fixedbin   => FixedBin();
        | Task.uncompress => Uncompress();
        | Task.fold       => Fold(o.backupfilename);
        | Task.rfold      => RestoreFold(o.backupfilename, o.backuplevel);
        | Task.unfold     => UnFold();
        | Task.incroots   => IncRoots();
        | Task.prm2red    => Prm2Red();
        | Task.red2prm    => Red2Prm();
        | Task.sortprm    => SortPrm();
        | Task.crunch     => Crunch();
        | Task.conta      => Conta();
        | Task.stats      => Stats();
        | Task.viewbits   => ViewBits();
        | Task.test       => Test();
        | Task.none       => <* ASSERT FALSE *>
        END
      EXCEPT 
      | Full =>
          Error("not enough space in DAG -- aborted")
      END
    END;
    Wr.Flush(stdout);
    FinalMessage(stderr);
    Wr.Flush(stderr)
  END Main;

PROCEDURE GetParms(): Parms =
  VAR o: Parms;
  BEGIN
    TRY
      WITH 
        cw = NEW(TextWr.T).init(),
        pp = NEW(ParseParams.T).init(stderr)
      DO
      
        PROCEDURE ParseAddSubParms() RAISES{ParseParams.Error} =
          BEGIN
            o.wordrd := ParseReader("-wordrd");
            o.rootnum := ParamUtil.GetInt(pp, cw, "-level", min := 0, default := 0);
          END ParseAddSubParms;
          
        PROCEDURE ParseFileName(key: TEXT): TEXT RAISES{ParseParams.Error} =
          BEGIN
            WITH f = ParamUtil.GetFileName(pp, cw, key) DO
              IF Text.Empty(f) THEN 
                pp.error("\"" & key & "\" not specified!")
              ELSIF Text.Equal(f, "-") THEN
                pp.error("\"" & key & "\" cannot be \"-\"!")
              END;
              RETURN f
            END
          END ParseFileName;

        PROCEDURE ParseReader(key: TEXT): Rd.T RAISES{ParseParams.Error} =
          BEGIN
            WITH rd = ParamUtil.GetRd(pp, cw, key) DO
              IF rd = NIL THEN 
                pp.error("\"" & key & "\" not specified!")
              ELSIF rd = stdin THEN
                pp.error("\"" & key & "\" cannot be standard input!")
              END;
              RETURN rd
            END
          END ParseReader;

        BEGIN  
          Wr.PutText(cw, "Goodies");
          IF ParamUtil.GetBool(pp, cw, "-new") THEN
            o.task := Task.new;
          ELSIF ParamUtil.GetBool(pp, cw, "-add") THEN
            ParseAddSubParms();
            o.task := Task.add;
          ELSIF ParamUtil.GetBool(pp, cw, "-sub") THEN
            ParseAddSubParms();
            o.task := Task.sub;
          ELSIF ParamUtil.GetBool(pp, cw, "-comment") THEN
            o.comment := ParamUtil.GetText(pp, cw, "-comment", default := "");
            o.task := Task.comment
          ELSIF ParamUtil.GetBool(pp, cw, "-compress") THEN
            o.task := Task.compress
          ELSIF ParamUtil.GetBool(pp, cw, "-fixedbin") THEN
            o.task := Task.fixedbin
          ELSIF ParamUtil.GetBool(pp, cw, "-uncompress") THEN
            o.task := Task.uncompress
          ELSIF ParamUtil.GetBool(pp, cw, "-conta") THEN
            o.task := Task.conta
          ELSIF ParamUtil.GetBool(pp, cw, "-crunch") THEN
            o.task := Task.crunch
          ELSIF ParamUtil.GetBool(pp, cw, "-fold") THEN
            o.backupfilename := ParseFileName("-backup");
            o.task := Task.fold
          ELSIF ParamUtil.GetBool(pp, cw, "-incroots") THEN
            o.task := Task.incroots
          ELSIF ParamUtil.GetBool(pp, cw, "-prm2red") THEN
            o.task := Task.prm2red
          ELSIF ParamUtil.GetBool(pp, cw, "-red2prm") THEN
            o.task := Task.red2prm
          ELSIF ParamUtil.GetBool(pp, cw, "-rfold") THEN
            o.backupfilename := ParseFileName("-backup");
            o.backuplevel := 
              ParamUtil.GetInt(pp, cw, "-level", min := 0, max := 1, default := 1);
            o.task := Task.rfold
          ELSIF ParamUtil.GetBool(pp, cw, "-sortprm") THEN
            o.task := Task.sortprm
          ELSIF ParamUtil.GetBool(pp, cw, "-spell") THEN
            o.rootnum := ParamUtil.GetInt(pp, cw, "-level", default := 0);
            o.task := Task.spell
          ELSIF ParamUtil.GetBool(pp, cw, "-stats") THEN
            o.task := Task.stats
          ELSIF ParamUtil.GetBool(pp, cw, "-unfold") THEN
            o.task := Task.unfold
          ELSIF ParamUtil.GetBool(pp, cw, "-viewbits") THEN
            o.task := Task.viewbits
          ELSIF ParamUtil.GetBool(pp, cw, "-test") THEN
            o.task := Task.test
          ELSE
            o.task := Task.none
          END;
          IF o.task = Task.none THEN 
            pp.error("no operation specified!")
          END;
          pp.finish();
        END;

        o.cmd := TextWr.ToText(cw);
      END
    EXCEPT
    | ParseParams.Error => 
        Wr.PutText(stderr, UsageText);
        Process.Exit(1);
    END;
    RETURN o
  END GetParms;
  
PROCEDURE InitialMessage(wr: Wr.T) =
  BEGIN
    Wr.PutText(wr, Separator);
    Wr.PutText(wr, "\n");
    Wr.PutText(wr, Util.FmtDate(Util.GetDate()));
    Wr.PutText(wr, ": Goodies started\n");
    Wr.Flush(stderr);
  END InitialMessage;

PROCEDURE FinalMessage(wr: Wr.T) =
  BEGIN
    Wr.PutText(wr, Util.FmtDate(Util.GetDate()));
    Wr.PutText(wr, ": Goodies finished.\n");
    Wr.Flush(stderr);
  END FinalMessage;
  
PROCEDURE Error(msg: TEXT) = 
  BEGIN
    Wr.PutText(stderr, "\n");
    Wr.PutText(stderr, msg);
    Wr.PutText(stderr, "\n");
    Wr.Flush(stderr);
    Process.Exit(1);
  END Error;
  
TYPE
  LoadProc = PROCEDURE (rd: Rd.T): PermDAG.T RAISES {Full};
  
PROCEDURE New() =
  BEGIN
    WITH permdag = PermDAG.New() DO
      permdag.Dump(stdout)
    END
  END New;

PROCEDURE AddSub(
    add: BOOLEAN := TRUE; 
    wordrd: Rd.T;
    rootnum: CARDINAL;
  ) RAISES {Full} =
  BEGIN
    WITH permdag = Load(PermDAG.Load) DO
      TRY
        LOOP
          WITH word = Trim(Rd.GetLine(wordrd)) DO
            permdag.AddSub(add, word, rootnum)
          END
        END
      EXCEPT
        Rd.EndOfFile => permdag.Crunch()
      END;
      permdag.Dump(stdout)
    END
  END AddSub;

PROCEDURE Trim(t: TEXT): TEXT =
  VAR ini, lim: CARDINAL;
  BEGIN
    WITH len = Text.Length(t) DO
      ini := 0; lim := len;
      WHILE ini < lim AND Text.GetChar(t, ini) = ' '   DO INC(ini) END;
      WHILE ini < lim AND Text.GetChar(t, lim-1) = ' ' DO DEC(lim) END;
      IF ini = 0 AND lim = len THEN
        RETURN t
      ELSE
        RETURN Text.Sub(t, ini, lim-ini)
      END
    END;
  END Trim;

PROCEDURE Compress() RAISES {Full} =
  BEGIN
    TRY
      WITH permdag = Load(PermDAG.Load) DO
        permdag.DumpCompr(stdout)
      END
    EXCEPT
      Huffman.OverFlow => Error("Huffman code too long")
    END
  END Compress;

PROCEDURE FixedBin() RAISES {Full} =
  BEGIN
    WITH permdag = Load(PermDAG.Load) DO
      permdag.DumpFixedBin(stdout)
    END
  END FixedBin;

PROCEDURE Uncompress() RAISES {Full} =
  BEGIN
    WITH permdag = Load(PermDAG.LoadCompr) DO
      permdag.Dump(stdout)
    END
  END Uncompress;

PROCEDURE Comment(comment: TEXT) RAISES {Full} =
  BEGIN
    WITH permdag = Load(PermDAG.Load) DO
      permdag.Comment(comment);
      permdag.Dump(stdout)
    END
  END Comment;

PROCEDURE Conta() RAISES {Full} =
  BEGIN
    WITH permdag = Load(PermDAG.Load) DO
      permdag.Conta(stdout)
    END
  END Conta;

PROCEDURE Crunch() RAISES {Full} =
  BEGIN
    WITH permdag = Load(PermDAG.Load) DO
      permdag.Crunch();
      permdag.Dump(stdout)
    END
  END Crunch;

PROCEDURE Fold(backupfilename: TEXT) RAISES {Full} =
  BEGIN
    WITH permdag = Load(PermDAG.Load) DO
      TRY
        permdag.Fold(backupfilename, backuplevel := 1, restore := FALSE)
      EXCEPT
        Huffman.OverFlow => Error("Huffman code too long")
      END;
      permdag.Dump(stdout)
    END
  END Fold;

PROCEDURE IncRoots() RAISES {Full} =
  BEGIN
    WITH permdag = Load(PermDAG.Load) DO
      permdag.IncRoots();
      permdag.Dump(stdout)
    END
  END IncRoots;

PROCEDURE Prm2Red() RAISES {Full} =
  BEGIN
    WITH permdag = PermDAG.Load(stdin) DO
      permdag.Prm2Red(stdout)
    END
  END Prm2Red;

PROCEDURE Red2Prm() RAISES {Full} =
  BEGIN
    WITH permdag = Load(PermDAG.Red2Prm) DO
      permdag.Dump(stdout)
    END
  END Red2Prm;

PROCEDURE RestoreFold(backupfilename: TEXT; backuplevel: CARDINAL) RAISES {Full} =
  BEGIN
    WITH permdag = PermDAG.New() DO
      TRY
        permdag.Fold(backupfilename, backuplevel, TRUE)
      EXCEPT
        Huffman.OverFlow => Error("Huffman code too long")
      END;
      permdag.Dump(stdout)
    END
  END RestoreFold;

PROCEDURE SortPrm() RAISES {Full} =
  VAR
    permdag : PermDAG.T;
  BEGIN
    permdag := Load(PermDAG.Load);
    permdag.SortPrm();
    permdag.Dump(stdout)
  END SortPrm;

PROCEDURE Spell(rootnum: CARDINAL) RAISES {Full} =
  BEGIN
    WITH permdag = Load(PermDAG.Load) DO
      permdag.Spell(stdout, rootnum)
    END
  END Spell;

PROCEDURE Stats() RAISES {Full} =
  BEGIN
    WITH permdag = Load(PermDAG.Load) DO 
      permdag.Statistics(stdout)
    END
  END Stats;

PROCEDURE UnFold() RAISES {Full} =
  BEGIN
    WITH permdag = Load(PermDAG.Load) DO
      permdag.UnFold();
      permdag.Dump(stdout)
    END
  END UnFold;

PROCEDURE ViewBits() =
  BEGIN
    Code.ViewBits(stdin, stdout)
  END ViewBits;

PROCEDURE Test() RAISES {Full} =
  CONST 
    MsgWidth    : CARDINAL = 20;
    CmpFileName : TEXT = "frances.fld.cmp";
    FldFileName : TEXT = "frances.fld";

  VAR before    : Time.T := Time.Now();

  PROCEDURE PutText(msg: TEXT) =
    BEGIN
      Wr.PutText(stdout, Fmt.Pad(msg, MsgWidth, align := Fmt.Align.Left));
      WITH 
        now = Time.Now(),
        elapsed = Util.FmtTime(now - before),
        exec = Util.GetExecTimesText()
      DO
        Wr.PutText(stdout, 
          Fmt.F("clk = %s usr = %s sys = %s tot = %s\n", 
            elapsed, exec.user, exec.system, exec.total
          )
        );
        Wr.PutText(stdout, "\n");
        Wr.Flush(stdout);
      END
    END PutText;

  <* FATAL OSError.E *>
  BEGIN
    PutText("início");

    WITH
      rd = FileRd.Open(CmpFileName)
    DO
      TRY
        LOOP
          EVAL Rd.GetChar(rd)
        END
      EXCEPT
        Rd.EndOfFile => Rd.Close(rd);
      END;
      PutText("getchar only")
    END;

    WITH
      rd = FileRd.Open(FldFileName)
    DO
      EVAL PermDAG.Load(rd);
      Rd.Close(rd);
      PutText("permdagload");
    END;

    WITH
      rd = FileRd.Open(CmpFileName)
    DO
      EVAL PermDAG.LoadCompr(rd);
      Rd.Close(rd);
      PutText("load comp")
    END;

    WITH
      rd = FileRd.Open(CmpFileName)
    DO
      EVAL ReadOnlyPermDAG.LoadCompr(rd);
      Rd.Close(rd);
      PutText("readonly load comp")
    END;
  END Test;

PROCEDURE Load(p: LoadProc): PermDAG.T RAISES {Full} =
  BEGIN
    RETURN p(stdin)
  END Load;
  
BEGIN
  Main()
END Goodies.
