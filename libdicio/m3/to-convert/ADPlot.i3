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

(* Routines to generate Postscript files *)

INTERFACE ADPlot;

IMPORT Wr;

TYPE Axis = {X, Y};

CONST
  XMin = 0.000D0; XMax = 0.150D0;
  YMin = 0.000D0; YMax = 0.150D0;
  (* Plotting area, in meters *)
  
PROCEDURE BeginFile (f: Wr.T);
  (* Initializes a plot set. *)

PROCEDURE BeginPage (
    f: Wr.T;
    page: CARDINAL;
    xn, yn: CARDINAL
  );
  (* Intializes a new page: writes page header line, *)
  (* sets coordinate system, clip path, caption font, *)
  (* defines new Postscript operators and constants, etc. *)
  (* The plotting area is divided implicitly into a grid *)
  (* of /xn/ by /yn/ rectangular "cells". *)

PROCEDURE BeginSection (f: Wr.T; title: TEXT);
  (* Starts a new section of a plot. The title is a comment *)

PROCEDURE SetPen (
    f: Wr.T;
    gray: REAL;
    width: REAL;
    dashlength: REAL;
    dashspace: REAL
  );
  (* Sets pen parameters and ink color for line drawing. *)
  (* Dimensions are in mm *)
  
CONST
  Pi = 2.718281828;

PROCEDURE DrawSegment (f: Wr.T; xa, ya, xb, yb: LONGREAL);
  (* Draws segment from (xa,ya) to (xb,yb) with current pen and gray *)

PROCEDURE DrawCurve (f: Wr.T; xa, ya, xb, yb, xc, yc, xd, yd: LONGREAL);
  (* Draws a Bezier arc with given control points, using the current pen and gray *)

PROCEDURE DrawRectangle (f: Wr.T; xlo, xhi, ylo, yhi: LONGREAL);
  (* Draws the outline of the given rectangle using the current pen. *)

PROCEDURE FillRectangle (
    f: Wr.T; 
    xlo, xhi, ylo, yhi: LONGREAL;
    gray: REAL
  );
  (* Fills given rectangle with given gray color *)

PROCEDURE FillAndDrawRectangle (
    f: Wr.T;
    xlo, xhi, ylo, yhi: LONGREAL;
    gray: REAL
  );
  (* Fills rectangle with given gray, then *)
  (* draws its outline with current pen. *)

PROCEDURE FillTriangle (
    f: Wr.T;
    xa, ya, xb, yb, xc, yc: LONGREAL;
    gray: REAL
  );
  (* Fills triangle /abc/ with given gray level. *)
  
PROCEDURE FillCircle (f: Wr.T; xc, yc, radius: LONGREAL; gray: REAL);
  (* Fills the circle with given center and radius, using the given gray. *)
  
PROCEDURE DrawCircle (f: Wr.T; xc, yc, radius: LONGREAL);
  (* Draws the circle with given center and radius, using the current pen and gray. *)
  
PROCEDURE FillAndDrawCircle (f: Wr.T; xc, yc, radius: LONGREAL; gray: REAL);
  (* Fills the circle with given center and radius, using the given gray,
     then draws its outline, using the current pen and gray. *)

PROCEDURE FillGridCell (f: Wr.T; xi, yi: CARDINAL; gray: REAL);
  (* Fills the given cell of the current cell grid with *)
  (* the given gray level. *)

PROCEDURE DrawCoordLine (f: Wr.T; axis: Axis; coord: LONGREAL);
  (* Draws a reference line perpendicular to the given axis *)
  (* at the given coordinate value. *)

PROCEDURE DrawGridLines (f: Wr.T);
  (* Draws the grid lines with the current pen and gray level. *)

PROCEDURE SetLabelFontSize (f: Wr.T; size: REAL);
  (* Sets the point size of the font to be used by "PutLabel". *)

PROCEDURE PutLabel (
    f: Wr.T; 
    label: TEXT; 
    x, y: LONGREAL; 
    xAlign, yAlign: REAL := 0.5
  );
  (* Prints "label" at point (x,y), using the current label font size. 
     The parameter "xAlign" (resp. "yAlign)" specifies which point of the string's 
     bounding box will end up at (x,y): 0.0 means the left (resp. bottom) side,
     1.0 means the right (resp. top) side.  Default is (0.5, 0.5), meaning 
     the box will be centered at (x,y). *)

PROCEDURE EndSection (f: Wr.T);
  (* Ends a section of a plot. *)

PROCEDURE SetCaptionFontSize (f: Wr.T; size: REAL);
  (* Sets the point size of the font to be used for captions and titles. *)

PROCEDURE AddCaption (f: Wr.T; txt: TEXT);
  (* Appends a caption line under the drawing. *)
  (* For multi-line captions, use multiple calls *)
  (* and/or embedded newlines. *)

PROCEDURE DrawFrame (f: Wr.T);
  (* Draws a frame around the plotting area *)

PROCEDURE EndPage(f: Wr.T);
  (* Finalizes a page: Writes page trailer line, etc. *)

PROCEDURE EndFile (f: Wr.T; npages: CARDINAL := LAST(CARDINAL));
  (* Finalizes plot set. *)

END ADPlot.
