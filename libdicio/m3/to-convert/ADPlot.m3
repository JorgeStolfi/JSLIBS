(****************************************************************************)
(* (C) Copyright 1992 Universidade Estadual de Campinas (UNICAMP)           *)
(*                    Campinas, SP, Brazil                                  *)
(*                                                                          *)
(* Authors:                                                                 *)
(*                                                                          *)
(*   Tomasz Kowaltowski  - CS Dept, UNICAMP <tomasz@dcc.unicamp.ansp.br>    *)
(*   Claudio L. Lucchesi - CS Dept, UNICAMP <lucchesi@dcc.unicamp.ansp.br>  *)
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

MODULE ADPlot;

IMPORT Wr, Fmt, Text, PSTools, Thread;
FROM Stdio IMPORT stderr;

CONST
  XScale = 1.0D0/(XMax - XMin);
  YScale = 1.0D0/(YMax - YMin);

CONST
  HSize = (XMax - XMin) / 0.0254D0 * 72.0D0;
  HMin = (8.5D0 * 72.0D0 - HSize) * 0.50D0;
  HMax = HMin + HSize;
  
  VSize = (YMax - YMin) / 0.0254D0 * 72.0D0;
  VMin = (11.0D0 * 72.0D0 - VSize) * 0.75D0;
  VMax = VMin + VSize;
  
PROCEDURE BeginFile (f: Wr.T) =
  BEGIN
    PSTools.BeginFile (f);
  END BeginFile;

PROCEDURE EndFile (f: Wr.T; npages: CARDINAL := LAST(CARDINAL)) =
  BEGIN
    PSTools.EndFile (f, npages);
  END EndFile;
  
PROCEDURE FIP(x: INTEGER; w: CARDINAL): TEXT =
  BEGIN
    RETURN Fmt.Pad(Fmt.Int(x), w)
  END FIP;

<*UNUSED*>
PROCEDURE FL(x: LONGREAL; d: CARDINAL): TEXT =
  BEGIN
    RETURN Fmt.LongReal(x, d, Fmt.Style.Flo)
  END FL;

PROCEDURE FLP(x: LONGREAL; d, w: CARDINAL): TEXT =
  BEGIN
    RETURN Fmt.Pad(Fmt.LongReal(x, d, Fmt.Style.Flo), w)
  END FLP;

PROCEDURE SubChar(txt: TEXT; a, b: CHAR): TEXT =
  BEGIN
    WITH i = Text.FindChar(txt, a) DO
      IF i = -1 THEN
        RETURN txt
      ELSE
        RETURN 
          Text.Sub(txt, 0, i) & 
          Text.FromChar(b) & 
          Text.Sub(txt, i+1, Text.Length(txt) - i - 1)
      END
    END
  END SubChar;

PROCEDURE EL(x: LONGREAL; d: CARDINAL := 6): TEXT =
  BEGIN
    RETURN SubChar(Fmt.LongReal(x, d, Fmt.Style.Sci), 'D', 'E')
  END EL;

PROCEDURE FR(x: REAL; d: CARDINAL): TEXT =
  BEGIN
    RETURN Fmt.Real(x, d, Fmt.Style.Flo)
  END FR;

PROCEDURE FRP(x: REAL; d, w: CARDINAL): TEXT =
  BEGIN
    RETURN Fmt.Pad(Fmt.Real(x, d, Fmt.Style.Flo), w)
  END FRP;

PROCEDURE BeginPage (
    f: Wr.T;
    page: CARDINAL;
    xn, yn: CARDINAL
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    WITH
      fxn = FLOAT(xn, LONGREAL),
      fyn = FLOAT(yn, LONGREAL)
    DO
      PSTools.BeginPage (f, page);
      Wr.PutText(f, "60 dict begin\n");
      Wr.PutText(f, "gsave\n");

      Wr.PutText(f, "% Round joints and caps:\n");
      Wr.PutText(f, "1 setlinecap 1 setlinejoin\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Black thin lines:\n");
      Wr.PutText(f, "0 setlinewidth 0 setgray [ ] 0 setdash\n");
      Wr.PutText(f, "\n");

      PSTools.SetCoords(f, PSTools.Axis.H, 0.0D0, 1.0D0, HMin, HMax);

      Wr.PutText(f, "/xmin 0.0 def  % min plottable x\n");
      Wr.PutText(f, "/xmax 1.0 def  % max plottable x\n");
      Wr.PutText(f, "/xn " & Fmt.Int(xn) & " def  % grid cells along x axis\n");
      Wr.PutText(f, "/xstep " & EL(1.0D0/fxn) & " def  % x-size of grid cell\n");

      PSTools.SetCoords (f, PSTools.Axis.V, 0.0D0, 1.0D0, VMin, VMax);

      Wr.PutText(f, "/ymin 0.0 def  % min plottable y\n");
      Wr.PutText(f, "/ymax 1.0 def  % max plottable y\n");
      Wr.PutText(f, "/yn " & Fmt.Int(yn) & " def  % grid cells along y axis\n");
      Wr.PutText(f, "/ystep " & EL(1.0D0/fyn) & " def  % y-size of grid cell\n");

      <* ASSERT (HMax - HMin) = (VMax - VMin) *>

      Wr.PutText(f, "% Units of measure:\n");
      Wr.PutText(f, "/pt " & EL(1.0D0/(VMax - VMin)) & " def\n");
      Wr.PutText(f, "/in pt 72.0 mul def \n");
      Wr.PutText(f, "/mm pt 72.0 25.4 div mul def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Set clipping path to boundary of plot area:\n");
      Wr.PutText(f, "newpath\n");
      Wr.PutText(f, "  xmin ymin moveto\n");
      Wr.PutText(f, "  xmax ymin lineto\n");
      Wr.PutText(f, "  xmax ymax lineto\n");
      Wr.PutText(f, "  xmin ymax lineto\n");
      Wr.PutText(f, "  xmin ymin lineto\n");
      Wr.PutText(f, "clip\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Caption text cursor:\n");
      Wr.PutText(f, "/captionx xmin def\n");
      Wr.PutText(f, "/captiony ymin def\n");
      Wr.PutText(f, "\n");
      
      SetCaptionFontSize (f, 10.0);

      SetLabelFontSize (f, 8.0);

      Wr.PutText(f, "% Rectangle draw operator:\n");
      Wr.PutText(f, "%   /xlo/ /xhi/ /ylo/ /yhi/ recd --> \n");
      Wr.PutText(f, "/recd\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  gsave\n");
      Wr.PutText(f, "    newpath\n");
      Wr.PutText(f, "      3 index 2 index moveto\n");
      Wr.PutText(f, "      2 index 2 index lineto\n");
      Wr.PutText(f, "      2 index 1 index lineto\n");
      Wr.PutText(f, "      3 index 1 index lineto\n");
      Wr.PutText(f, "      pop pop pop pop\n");
      Wr.PutText(f, "      closepath\n");
      Wr.PutText(f, "     stroke\n");
      Wr.PutText(f, "  grestore\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Rectangle fill operator:\n");
      Wr.PutText(f, "%   /xlo/ /xhi/ /ylo/ /yhi/ /gray/ recf --> \n");
      Wr.PutText(f, "/recf\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  gsave\n");
      Wr.PutText(f, "    setgray\n");
      Wr.PutText(f, "    newpath\n");
      Wr.PutText(f, "      3 index 2 index moveto\n");
      Wr.PutText(f, "      2 index 2 index lineto\n");
      Wr.PutText(f, "      2 index 1 index lineto\n");
      Wr.PutText(f, "      3 index 1 index lineto\n");
      Wr.PutText(f, "      pop pop pop pop\n");
      Wr.PutText(f, "      closepath\n");
      Wr.PutText(f, "    fill\n");
      Wr.PutText(f, "  grestore\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Rectangle fill and stroke operator:\n");
      Wr.PutText(f, "%   /xlo/ /xhi/ /ylo/ /yhi/ /gray/ recfd --> \n");
      Wr.PutText(f, "/recfd\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  gsave\n");
      Wr.PutText(f, "    5 1 roll\n");
      Wr.PutText(f, "    newpath\n");
      Wr.PutText(f, "      3 index 2 index moveto\n");
      Wr.PutText(f, "      2 index 2 index lineto\n");
      Wr.PutText(f, "      2 index 1 index lineto\n");
      Wr.PutText(f, "      3 index 1 index lineto\n");
      Wr.PutText(f, "      pop pop pop pop\n");
      Wr.PutText(f, "      closepath\n");
      Wr.PutText(f, "      gsave setgray fill grestore\n");
      Wr.PutText(f, "    stroke\n");
      Wr.PutText(f, "  grestore\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Circle fill operator:\n");
      Wr.PutText(f, "%   /x/ /y/ /radius/ /gray/ cirf --> \n");
      Wr.PutText(f, "/cirf\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  gsave\n");
      Wr.PutText(f, "    setgray\n");
      Wr.PutText(f, "    newpath\n");
      Wr.PutText(f, "      0 360 arc\n");
      Wr.PutText(f, "      closepath\n");
      Wr.PutText(f, "    fill\n");
      Wr.PutText(f, "  grestore\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Circle draw operator:\n");
      Wr.PutText(f, "%   /x/ /y/ /radius/ cird --> \n");
      Wr.PutText(f, "/cird\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  gsave\n");
      Wr.PutText(f, "    newpath\n");
      Wr.PutText(f, "      0 360 arc\n");
      Wr.PutText(f, "      closepath\n");
      Wr.PutText(f, "    stroke\n");
      Wr.PutText(f, "  grestore\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Circle fill and draw operator:\n");
      Wr.PutText(f, "%   /x/ /y/ /radius/ /gray/ cirfd --> \n");
      Wr.PutText(f, "/cirfd\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  gsave\n");
      Wr.PutText(f, "    4 1 roll\n");
      Wr.PutText(f, "    newpath\n");
      Wr.PutText(f, "      0 360 arc\n");
      Wr.PutText(f, "      closepath\n");
      Wr.PutText(f, "      gsave setgray fill grestore\n");
      Wr.PutText(f, "    stroke\n");
      Wr.PutText(f, "  grestore\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Cell fill operator:\n");
      Wr.PutText(f, "%   /xi/ /yi/ celf --> \n");
      Wr.PutText(f, "/celf\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  3 1 roll \n");
      Wr.PutText(f, "  exch dup \n");
      Wr.PutText(f, "  xstep mul xmin add exch 1 add xstep mul xmin add\n");
      Wr.PutText(f, "  3 2 roll dup\n");
      Wr.PutText(f, "  ystep mul ymin add exch 1 add ystep mul ymin add\n");
      Wr.PutText(f, "  5 4 roll \n");
      Wr.PutText(f, "  recf\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Segment draw operator:\n");
      Wr.PutText(f, "%   /xa/ /ya/ /xb/ /yb/ segd --> \n");
      Wr.PutText(f, "/segd\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  gsave\n");
      Wr.PutText(f, "    newpath\n");
      Wr.PutText(f, "      moveto\n");
      Wr.PutText(f, "      lineto\n");
      Wr.PutText(f, "    stroke\n");
      Wr.PutText(f, "  grestore\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Curve draw operator:\n");
      Wr.PutText(f, "%   /xa/ /ya/  /xb/ /yb/  /xc/ /yc/  /xd/ /yd/ arcd --> \n");
      Wr.PutText(f, "/arcd\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  gsave\n");
      Wr.PutText(f, "    newpath\n");
      Wr.PutText(f, "      8 -2 roll moveto curveto\n");
      Wr.PutText(f, "    stroke\n");
      Wr.PutText(f, "  grestore\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Triangle fill operator:\n");
      Wr.PutText(f, "%   /xa/ /ya/ /xb/ /yb/ /xc/ /yc/ /gray/ trif --> \n");
      Wr.PutText(f, "/trif\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  gsave\n");
      Wr.PutText(f, "    setgray\n");
      Wr.PutText(f, "    newpath\n");
      Wr.PutText(f, "      moveto\n");
      Wr.PutText(f, "      lineto\n");
      Wr.PutText(f, "      lineto\n");
      Wr.PutText(f, "      closepath\n");
      Wr.PutText(f, "    fill\n");
      Wr.PutText(f, "  grestore\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Draw an X-value grid line:\n");
      Wr.PutText(f, "%   /x/ xgrd --> \n");
      Wr.PutText(f, "/xgrd\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  gsave\n");
      Wr.PutText(f, "    newpath\n");
      Wr.PutText(f, "      dup ymin moveto\n");
      Wr.PutText(f, "      ymax lineto\n");
      Wr.PutText(f, "    stroke\n");
      Wr.PutText(f, "  grestore\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Draw an Y-value grid line:\n");
      Wr.PutText(f, "%   /y/ ygrd --> \n");
      Wr.PutText(f, "/ygrd\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  gsave\n");
      Wr.PutText(f, "    newpath\n");
      Wr.PutText(f, "      dup xmin exch moveto\n");
      Wr.PutText(f, "      xmax exch lineto\n");
      Wr.PutText(f, "    stroke\n");
      Wr.PutText(f, "  grestore\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Label printing operator:\n");
      Wr.PutText(f, "%   /str/ /xa/ /ya/ /xc/ /yc/ lbsh --> \n");
      Wr.PutText(f, "/lbsh\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  labelfont setfont\n");
      Wr.PutText(f, "  newpath moveto\n");
      Wr.PutText(f, "    % --- str, xa, ya\n");
      Wr.PutText(f, "  gsave 2 index false charpath flattenpath pathbbox grestore\n");
      Wr.PutText(f, "    % --- str, xa, ya, lox, loy, hix, hiy\n");
      Wr.PutText(f, "  3 index 3 index currentpoint \n");
      Wr.PutText(f, "    % --- str, xa, ya, lox, loy, hix, hiy, lox, loy, cx, cy\n");
      Wr.PutText(f, "  exch 4 1 roll exch sub\n");
      Wr.PutText(f, "  3 1 roll sub exch\n");
      Wr.PutText(f, "    % --- str, xa, ya, lox, loy, hix, hiy, cx-lox, cy-loy\n");
      Wr.PutText(f, "  rmoveto\n");
      Wr.PutText(f, "    % --- str, xa, ya, lox, loy, hix, hiy\n");
      Wr.PutText(f, "  exch 4 1 roll exch sub \n");
      Wr.PutText(f, "  3 1 roll sub exch\n");
      Wr.PutText(f, "    % --- str, xa, ya, dx, dy\n");
      Wr.PutText(f, "  exch 4 1 roll mul -1 mul\n");
      Wr.PutText(f, "  3 1 roll mul -1 mul exch\n");
      Wr.PutText(f, "    % --- str, -dx*xa, -dy*ya\n");
      Wr.PutText(f, "  rmoveto\n");
      Wr.PutText(f, "    % --- str\n");
      Wr.PutText(f, "  show\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Operator to move to new caption line:\n");
      Wr.PutText(f, "%   nl --> \n");
      Wr.PutText(f, "/nl\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  /captiony captiony captiondy sub def\n");
      Wr.PutText(f, "  captionx captiony moveto\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.PutText(f, "% Operator to print caption string at CP without clipping:\n");
      Wr.PutText(f, "%   /s/ shw --> \n");
      Wr.PutText(f, "/shw\n");
      Wr.PutText(f, "{\n");
      Wr.PutText(f, "  captionfont setfont\n"); 
      Wr.PutText(f, "  gsave initclip show grestore\n");
      Wr.PutText(f, "} def\n");
      Wr.PutText(f, "\n");

      Wr.Flush(f);
    END
  END BeginPage;

PROCEDURE EndPage(f: Wr.T) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(f, "grestore\n");
    Wr.PutText(f, "% Now we are back to the standard coord system.\n");
    Wr.PutText(f, "\n");

    Wr.PutText(f, "end\n");
    Wr.PutText(f, "\n");
    PSTools.EndPage (f);
  END EndPage;

PROCEDURE SetCaptionFontSize (f: Wr.T; size: REAL) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(f, "% Caption font setup:\n");
    Wr.PutText(f, "/captiondy " & FR(size, 3) & " pt mul def\n");
    Wr.PutText(f, "/captionfont\n");
    Wr.PutText(f, "  /Courier findfont captiondy scalefont\n");
    Wr.PutText(f, "def\n");
    Wr.PutText(f, "\n");
    Wr.Flush(f);
  END SetCaptionFontSize;
  
PROCEDURE AddCaption (f: Wr.T; txt: TEXT) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(f, "nl ");
    PSTools.PutText(f, txt, ") shw\nnl (");
    Wr.PutText(f, " shw\n");
    Wr.Flush(f);
  END AddCaption;

PROCEDURE BeginSection (f: Wr.T; title: TEXT) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(f, "% " & title & "\n");
    Wr.PutText(stderr, "[" & title & "]\n");
    Wr.Flush(f);
  END BeginSection;

PROCEDURE EndSection (f: Wr.T) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(f, "\n");
    Wr.Flush(f);
  END EndSection;

PROCEDURE DrawFrame (f: Wr.T) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    BeginSection (f, "Draw frame around plot area");
    Wr.PutText(f, "gsave\n");
    Wr.PutText(f, "% Assumes xmax, xmin, ymax, ymin are defined.\n");
    Wr.PutText(f, "  initclip\n");
    Wr.PutText(f, "  newpath\n");
    Wr.PutText(f, "  xmin ymin moveto\n");
    Wr.PutText(f, "  xmax ymin lineto\n");
    Wr.PutText(f, "  xmax ymax lineto\n");
    Wr.PutText(f, "  xmin ymax lineto\n");
    Wr.PutText(f, "  xmin ymin lineto\n");
    Wr.PutText(f, "  closepath stroke\n");
    Wr.PutText(f, "grestore\n");
    EndSection (f);
  END DrawFrame;

PROCEDURE SetPen (
    f: Wr.T;
    gray: REAL;
    width: REAL;
    dashlength: REAL;
    dashspace: REAL
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(f, FR(gray, 3) & " setgray\n");
    Wr.PutText(f, "mm " & FR(width, 3) & " mul setlinewidth\n");
    IF dashlength = 0.0 OR dashspace = 0.0 THEN
      Wr.PutText(f, "[ ] 0 setdash\n")
    ELSE
      Wr.PutText(f,
        "[ mm " & FR(dashlength, 6) & " mul " & 
        "  mm " & FR(dashspace, 6) & " mul ] 0 setdash\n"
      );
    END;
    Wr.PutText(f, "\n");
    Wr.Flush(f);
  END SetPen;

PROCEDURE DrawSegment (f: Wr.T; xa, ya, xb, yb: LONGREAL) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    WITH
      psxa = XScale * (xa - XMin),
      psya = YScale * (ya - YMin),
      psxb = XScale * (xb - XMin),
      psyb = YScale * (yb - YMin)
    DO
      Wr.PutText(f,
        FLP(psxa, 4, 6) & " " & FLP(psya, 4, 6) & "  " & 
        FLP(psxb, 4, 6) & " " & FLP(psyb, 4, 6) & " segd\n"
      );
    END;
    Wr.Flush(f);
  END DrawSegment;  

PROCEDURE DrawCurve (f: Wr.T; xa, ya, xb, yb, xc, yc, xd, yd: LONGREAL) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    WITH
      psxa = XScale * (xa - XMin),
      psya = YScale * (ya - YMin),
      psxb = XScale * (xb - XMin),
      psyb = YScale * (yb - YMin),
      psxc = XScale * (xc - XMin),
      psyc = YScale * (yc - YMin),
      psxd = XScale * (xd - XMin),
      psyd = YScale * (yd - YMin)
    DO
      Wr.PutText(f,
        FLP(psxa, 4, 6) & " " & FLP(psya, 4, 6) & "  " & 
        FLP(psxb, 4, 6) & " " & FLP(psyb, 4, 6) & "  " & 
        FLP(psxc, 4, 6) & " " & FLP(psyc, 4, 6) & "  " & 
        FLP(psxd, 4, 6) & " " & FLP(psyd, 4, 6) & " arcd\n"
      );
    END;
    Wr.Flush(f);
  END DrawCurve;  
  
PROCEDURE AuxDrawRectangle(
    f: Wr.T; 
    xlo, xhi, ylo, yhi: LONGREAL;
    gray: REAL;
    operator: TEXT
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    WITH
      psxlo = XScale * (xlo - XMin),
      psxhi = XScale * (xhi - XMin),
      psylo = YScale * (ylo - YMin),
      psyhi = YScale * (yhi - YMin)
    DO
      Wr.PutText(f, 
        FLP(psxlo, 4, 6) & " " & FLP(psxhi, 4, 6) & "  " & 
        FLP(psylo, 4, 6) & " " & FLP(psyhi, 4, 6)
      );
      IF gray >= 0.0 THEN Wr.PutText(f, "  " & FR(gray, 4)) END;
      Wr.PutText(f, " " & operator & "\n");
      Wr.Flush(f);
    END
  END AuxDrawRectangle;

PROCEDURE DrawRectangle (f: Wr.T; xlo, xhi, ylo, yhi: LONGREAL) =
  BEGIN
    AuxDrawRectangle (f, xlo, xhi, ylo, yhi, -1.0, "recd");
  END DrawRectangle;

PROCEDURE FillRectangle (
    f: Wr.T; 
    xlo, xhi, ylo, yhi: LONGREAL;
    gray: REAL
  ) =
  BEGIN
    AuxDrawRectangle (f, xlo, xhi, ylo, yhi, gray, "recf");
  END FillRectangle;

PROCEDURE FillAndDrawRectangle (
    f: Wr.T;
    xlo, xhi, ylo, yhi: LONGREAL;
    gray: REAL
  ) =
  BEGIN
    AuxDrawRectangle (f, xlo, xhi, ylo, yhi, gray, "recfd");
  END FillAndDrawRectangle;

PROCEDURE AuxDrawCircle (
    f: Wr.T; 
    xc, yc, radius: LONGREAL;
    gray: REAL;
    operator: TEXT
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    WITH
      psxc = XScale * (xc - XMin),
      psyc = YScale * (yc - YMin),
      psradius = YScale * radius
    DO
      Wr.PutText(f, 
        FLP(psxc, 4, 6) & " " & FLP(psyc, 4, 6) & "  " & 
        FLP(psradius, 4, 6)
      );
      IF gray >= 0.0 THEN Wr.PutText(f, "  " & FR(gray, 4)) END;
      Wr.PutText(f, " " & operator & "\n");
      Wr.Flush(f);
    END
  END AuxDrawCircle;

PROCEDURE FillCircle (f: Wr.T; xc, yc, radius: LONGREAL; gray: REAL) =
  BEGIN
    AuxDrawCircle (f, xc, yc, radius, gray, "cirf");
  END FillCircle;
  
PROCEDURE DrawCircle (f: Wr.T; xc, yc, radius: LONGREAL) =
  BEGIN
    AuxDrawCircle (f, xc, yc, radius, -1.0, "cird");
  END DrawCircle;
  
PROCEDURE FillAndDrawCircle (f: Wr.T; xc, yc, radius: LONGREAL; gray: REAL) =
  BEGIN
    AuxDrawCircle (f, xc, yc, radius, gray, "cirfd");
  END FillAndDrawCircle;

PROCEDURE FillTriangle (
    f: Wr.T;
    xa, ya, xb, yb, xc, yc: LONGREAL;
    gray: REAL
  ) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    WITH
      psxa = XScale * (xa - XMin),
      psya = YScale * (ya - YMin),
      psxb = XScale * (xb - XMin),
      psyb = YScale * (yb - YMin),
      psxc = XScale * (xc - XMin),
      psyc = YScale * (yc - YMin)
    DO
      Wr.PutText(f,
        FLP(psxa, 4, 6) & " " & FLP(psya, 4, 6) & "  " & 
        FLP(psxb, 4, 6) & " " & FLP(psyb, 4, 6) & "  " & 
        FLP(psxc, 4, 6) & " " & FLP(psyc, 4, 6) & "  " &
        FRP(gray, 4, 6) & " trif\n"
      );
      Wr.Flush(f);
    END
  END FillTriangle;

PROCEDURE FillGridCell (f: Wr.T; xi, yi: CARDINAL; gray: REAL) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(f, 
      FIP(xi, 3) & " " & FIP(yi, 3) & "  " & FRP(gray, 3, 4) & " celf\n"
    );
    Wr.Flush(f);
  END FillGridCell;

PROCEDURE DrawCoordLine (f: Wr.T; axis: Axis; coord: LONGREAL) =
  VAR pscoord: LONGREAL;
      op: TEXT;
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    IF axis = Axis.X THEN
      pscoord := XScale * (coord - XMin); op := "xgrd"
    ELSE
      pscoord := YScale * (coord - YMin); op := "ygrd"
    END;
    Wr.PutText(f, FLP(pscoord, 4, 6) & " " & op & "\n");
  END DrawCoordLine;

PROCEDURE DrawGridLines (f: Wr.T) = 
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(f, "% Grid lines:\n");
    Wr.PutText(f, "gsave\n");
    Wr.PutText(f, "  initclip\n");
    Wr.PutText(f, "  0 1 xn {\n");
    Wr.PutText(f, "    xstep mul xMin add xgrd\n");
    Wr.PutText(f, "  } for\n");
    Wr.PutText(f, "  0 1 yn {\n");
    Wr.PutText(f, "    ystep mul yMin add ygrd\n");
    Wr.PutText(f, "  } for\n");
    Wr.PutText(f, "grestore\n");
    Wr.PutText(f, "\n");
    Wr.Flush(f);
  END DrawGridLines;
  
PROCEDURE SetLabelFontSize (f: Wr.T; size: REAL) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    Wr.PutText(f, "% Label font setup:\n");
    Wr.PutText(f, "/labelfont\n");
    Wr.PutText(f, "  /Courier findfont " & FR(size, 3) & " pt mul scalefont\n");
    Wr.PutText(f, "def\n");
    Wr.PutText(f, "\n");
    Wr.Flush(f);
  END SetLabelFontSize;
  
PROCEDURE PutLabel (f: Wr.T; label: TEXT; x, y: LONGREAL; xAlign, yAlign: REAL) =
  <* FATAL Wr.Failure, Thread.Alerted *>
  BEGIN
    PSTools.PutText(f, label, "\\267");
    WITH
      psx = XScale * (x - XMin),
      psy = YScale * (y - YMin)
    DO
      Wr.PutText(f, " ");
      Wr.PutText(f, 
        FRP(xAlign, 3, 5) & " " & FRP(yAlign, 3, 5) & "  " &
        FLP(psx, 4, 6) & " " & FLP(psy, 4, 6) & "  lbsh\n"
      );
    END;
    Wr.Flush(f);
  END PutLabel;

BEGIN
END ADPlot.


