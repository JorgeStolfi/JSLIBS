/* See epswr.h */
/* Last edited on 2022-10-20 06:54:34 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <affirm.h>
#include <bool.h>
#include <jsfile.h>
#include <jstime.h>
#include <jsstring.h>

#include <epswr.h>
#include <epswr_def.h>
#include <epswr_vis.h>
#include <epswr_font_list.h>

#include <epswr_dev.h>

/* INTERNAL PROTOTYPES */

void epswr_dev_write_file_header(FILE *wr, double hSize, double vSize);
  /* Writes to {wr} the file preamble for an encapsulated Postscript file.
    
    The preamble defines some auxiliary Postscript operators and constants.
    It sets the line capt and join style to "round", and the Postscript 
    variables {hSize,vSize} to the specified value. The figure's bounding box (for the "%%BoundingBox"
    comment) is set to {[0 _ hSize] x [0 _ vSize]}, 
    The parameters {hSize} and {vSize} must be positive exact 
    integers, in points. 
    
    !!! Change them to {uint32_t}. !!!
    
    The procedure does not set any other Postscript variables, such as
    the label and text fonts, the plotting rectange, the fill color,
    etc.. Therefore, some of the auxiliary Postscript operators may not
    work. Other procedures below must be called to define those
    variables. */

void epswr_dev_write_file_trailer(FILE *wr, int32_t nfonts, char **fonts);
  /* Writes to {wr} the file postamble for an encapsulated Postscript file,
    namely the "%%Trailer" comment followed by the "%%DocumentFonts:
    {F}" (where {F} is the given list of font names). */

void epswr_dev_text_line
  ( epswr_figure_t *eps, 
    const char *line,
    const char *strut,
    double hAbs, double vAbs, 
    double rot, 
    bool_t clipped, 
    double hAlign, double vAlign,
    bool_t fill, bool_t draw
  );
  /* Prints the first line of {line} to {eps},
    with the current font and font size.
    
    The point of text's bounding box specified by {hAlign,vAlign} 
    will be at absolute Device coordinates {(hAbs,vAbs)},
    and the text will be rotated by {rot} degrees counterclockwise
    around that point.  The vertical extent of the bounding box
    will be that of the string {strut}, which is otherwise ignored.
    
    The end of the {line} is assumed to be at the first character
    that is null or a line break ('\000', '\012', or '\015').
    Tabs are printed as two spaces. The interpretation of any
    other characters depends on the font.
    
    If {clipped} is true, the text (after rotation) will be clipped 
    to the current plotting area; otherwise it may extend over the 
    whole figure. */ 

void epswr_dev_write_font_list(FILE *psFile, int32_t nfonts, char **fonts);
  /* Writes the structured comment "%%DocumentFonts: {F}" 
    where {F} is the given list of font names. */

void epswr_dev_write_font_set_cmds(FILE *wr, const char *psname, const char *font, double size);
  /* Writes to Postscript file {wr} the command that sets the Postscript variable 
    whose name is {psname} to the specified font scaled to the specified size.  The
    {psname} must be "labelfont" or "textfont". */
    
/* MISCELLANEOUS WRITING PROCEDURES */

void epswr_dev_write_proc_defs(FILE *psFile);
  /* Writes to {psFile} the definition of the operators assumed
    by other plotting and text writing commands. */

void epswr_dev_write_color(FILE *psFile, double *clr);
  /* Writes the color value {clr[0..2]} to the Postscript file,
    with 3 decimals, preceded by spaces.  */

void epswr_dev_compute_joining_arc_radius
  ( double ptx0, double pty0, 
    double ptx1, double pty1, 
    double ptx2, double pty2, 
    double rad,
    double *arad
  );
  /* Checks whether the corner between points {pt0=(ptx0,pty0)},
    {pt1=(ptx1,pty1)}, and {pt2=(ptx2,pty2)} can be rounded with an
    arc of radius {rad}. If so, sets {*arad} to {rad}. Otherwise
    sets {*arad} to zero. */

/* INTERNAL FILE WRITING PROCEDURES */

/* The file preamble (produced by both {epswr_write_file_header})
  creates a Postscript dictionary
  "eps$dict" that contains definitions for special Postscript
  operators and some state variables, such as the plotting rectangle. */

void epswr_dev_write_file_header
  ( FILE *wr, 
    double hSize, 
    double vSize
  )
  { fprintf(wr, "%%!PS-Adobe-3.0 EPS-2.0\n");

    char *date = today();
    fprintf(wr, "%%%%CreationDate: %s\n", date);
    free(date);

    demand(hSize == floor(hSize) && hSize > 0.0, "invalid hSize"); 
    demand(vSize == floor(vSize) && vSize > 0.0, "invalid hSize"); 
    fprintf(wr, "%%%%BoundingBox: %.0f %.0f %.0f %.0f\n", 0.0, 0.0, hSize, vSize);

    fprintf(wr, "%%%%Pages: 1\n");
    fprintf(wr, "%%%%DocumentFonts: (atend)\n");
    fprintf(wr, "%%%%EndComments\n");
    
    fprintf(wr, 
      "/epswr$dict 6400 dict def \n"
      "epswr$dict begin\n"
      "%% Round joints and caps:\n"
      "1 setlinecap 1 setlinejoin\n"
      "\n"
    );
    epswr_dev_write_proc_defs(wr);
  }

void epswr_dev_write_file_trailer(FILE *wr, int32_t nFonts, char **fonts)
  { fprintf(wr, "%%%%Trailer\n" );
    epswr_dev_write_font_list(wr, nFonts, fonts);
  }

void epswr_dev_write_font_list(FILE *wr, int32_t nFonts, char **fonts)
  { int32_t i;
    fprintf(wr, "%%%%DocumentFonts:");
    for (i = 0; i < nFonts; i++)
      { fprintf(wr, " %s", fonts[i]); }
    fprintf(wr, "\n");
  }

void epswr_dev_write_window_set_cmds
  ( FILE *wr,
    double hMin, double hMax,
    double vMin, double vMax
  )
  { 
    fprintf(wr, "/hMin %f def   %% min plottable x\n", hMin);
    fprintf(wr, "/hMax %f def   %% max plottable x\n", hMax);

    fprintf(wr, "/vMin %f def   %% min plottable y\n", vMin);
    fprintf(wr, "/vMax %f def   %% max plottable y\n", vMax);

    fprintf(wr, 
      "%% Set clipping path to boundary of plot area:\n"
      "initclip\n"
      "newpath\n"
      "  hMin vMin moveto\n"
      "  hMax vMin lineto\n"
      "  hMax vMax lineto\n"
      "  hMin vMax lineto\n"
      "  hMin vMin lineto\n"
      "clip\n"
      "\n"
    );
    fflush(wr);
  }

void epswr_dev_write_label_font_set_cmds(FILE *wr, const char *font, double size)
  { epswr_dev_write_font_set_cmds(wr, "labelFont", font, size); }

void epswr_dev_write_text_font_set_cmds(FILE *wr, const char *font, double size)
  { epswr_dev_write_font_set_cmds(wr, "textFont", font, size); }

void epswr_dev_write_font_set_cmds(FILE *wr, const char *psname, const char *font, double size)
  { fprintf(wr, 
      "/%s\n /%s findfont %.3f pt mul scalefont def\n",
      psname, font, size
    );
  }  

void epswr_dev_write_color(FILE *wr, double *clr)
  { fprintf(wr, "  %5.3f %5.3f %5.3f", clr[0], clr[1], clr[2]); }

void epswr_dev_write_fill_color_set_cmds(FILE *wr, double *fc)
  { epswr_dev_write_color(wr, fc);
    fprintf(wr, " sfc\n");
  }

void epswr_dev_write_ps_string(FILE *wr, const char *text)
  { const char *p = text;
    putc('(', wr);
    while((*p != 0) && (*p != '\n') && (*p != '\r'))
      { if ((*p == '(') || (*p == ')') || (*p == '\\'))
          { putc('\\', wr); putc(*p, wr); }
        else if ((*p >= ' ') && (*p <= '~'))
          { putc(*p, wr); }
        else if (*p == '\t')
          { putc(' ', wr); putc(' ', wr); }
        else
          { fprintf(wr, "\\%03o", *p); }
        p++;
      }
    fprintf(wr, ")");
  }

void epswr_dev_write_proc_defs(FILE *wr)
  { 
    /* Global constants and variables: */
    
    fprintf(wr, 
      "%% Units of measure:\n"
      "/pt 1.0 def\n"
      "/in pt 72.0 mul def \n"
      "/mm pt 72.0 25.4 div mul def\n"
      "\n"
    );
  
    fprintf(wr, 
      "%% Segment draw operator:\n"
      "%%   {xa} {ya} {xb} {yb} segd --> \n"
      "/segd\n"
      "{ gsave\n"
      "    newpath\n"
      "    moveto\n"
      "    lineto\n"
      "    stroke\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fprintf(wr, 
      "%% Curve draw operator:\n"
      "%%   {xa} {ya}  {xb} {yb}  {xc} {yc}  {xd} {yd} arcd --> \n"
      "/arcd\n"
      "{ gsave\n"
      "    newpath\n"
      "      8 -2 roll moveto curveto\n"
      "    stroke\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fprintf(wr, 
      "%% Draw an X-value grid line:\n"
      "%%   {x} xgrd --> \n"
      "/xgrd\n"
      "{ gsave\n"
      "  initclip\n"
      "  newpath\n"
      "    dup vMin moveto\n"
      "    vMax lineto\n"
      "    stroke\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fprintf(wr, 
      "%% Draw an Y-value grid line:\n"
      "%%   {y} ygrd --> \n"
      "/ygrd\n"
      "{ gsave\n"
      "  initclip\n"
      "  newpath\n"
      "    dup hMin exch moveto\n"
      "    hMax exch lineto\n"
      "    stroke\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fprintf(wr,
      "%% Draw a frame around the current plot window:\n"
      "%%   frame --> \n"
      "/wframe\n"
      "{ gsave\n"
      "    initclip\n"
      "    newpath\n"
      "    hMin vMin moveto\n"
      "    hMax vMin lineto\n"
      "    hMax vMax lineto\n"
      "    hMin vMax lineto\n"
      "    hMin vMin lineto\n"
      "    closepath stroke\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );
    
    /* Setting the fill color - variables {fcR,fcG,fcB}: */
    
    fprintf(wr, 
      "%% Set fill color operator:\n"
      "%%   {R} {G} {B} sfc --> \n"
      "/sfc\n"
      "{ /fcB exch def\n"
      "  /fcG exch def\n"
      "  /fcR exch def\n"
      "} bind def\n"
      "\n"
    );

    /* Figures: */
    
    fprintf(wr, 
      "%% Fill and/or stroke operator:\n"
      "%%   {draw} {fill} fs --> \n"
      "/fs\n"
      "{ dup 2 index add\n"
      "  %% --- draw fill draw+fill\n"
      "  dup 0 eq\n"
      "  %% --- draw fill draw+fill neither\n"
      "    { pop pop pop newpath }\n"
      "    {\n"
      "      %% --- draw fill draw+fill\n"
      "      2 eq\n"
      "        { gsave fcR fcG fcB setrgbcolor fill grestore stroke pop pop }\n"
      "        { 1 eq { fcR fcG fcB setrgbcolor fill } if\n"
      "          1 eq { stroke } if\n"
      "        }\n"
      "      ifelse\n"
      "    }\n"
      "  ifelse\n"
      "} bind def\n"
      "\n"
    );
    
    fprintf(wr, 
      "%% Even-odd fill and/or stroke operator:\n"
      "%%   {draw} {fill} eofs --> \n"
      "/eofs\n"
      "{ dup 2 index add\n"
      "  %% --- draw fill draw+fill\n"
      "  dup 0 eq\n"
      "  %% --- draw fill draw+fill neither\n"
      "    { pop pop pop newpath }\n"
      "    {\n"
      "      %% --- draw fill draw+fill\n"
      "      2 eq\n"
      "        { gsave fcR fcG fcB setrgbcolor eofill grestore stroke pop pop }\n"
      "        { 1 eq { fcR fcG fcB setrgbcolor eofill } if\n"
      "          1 eq { stroke } if\n"
      "        }\n"
      "      ifelse\n"
      "    }\n"
      "  ifelse\n"
      "} bind def\n"
      "\n"
    );
    
    fprintf(wr, 
      "%% Generic fill and/or stroke operator:\n"
      "%%   {draw} {fill} {eo} anyfs --> \n"
      "/anyfs\n"
      "{ 1 eq\n"
      "    { eofs }\n"
      "    { fs }\n"
      "  ifelse\n"
      "} bind def\n"
      "\n"
    );
    
    fprintf(wr, 
      "%% Triangle operator:\n"
      "%%   {draw} {fill} {xa} {ya} {xb} {yb} {xc} {yc} tri --> \n"
      "/tri\n"
      "{ gsave\n"
      "    newpath\n"
      "    moveto\n"
      "    lineto\n"
      "    lineto\n"
      "    closepath\n"
      "    fs\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );
    
    /* For smooth triangle shading we use {shfill} with {ShadingType = 4}. */
    /* The expected order of fields in the {DataSource} is: */
    /* [ 0 xa ya Ra Ga Ba  0 xb yb Rb Gb Bb  0 xc yc Rc Gc Bc ] */

    fprintf(wr, 
      "%% Gouraud-shaded triangle fill operator\n"
      "%%   (always fills, never strokes):\n"
      "%%   {xa} {ya} {Ra} {Ga} {Ba} \\\n"
      "%%   {xb} {yb} {Rb} {Gb} {Bb} \\\n"
      "%%   {xc} {yc} {Rc} {Gc} {Bc} gstrif --> \n"
      "/gstrif\n"
      "{ [ 16 1 roll\n"
      "    0 6 1 roll 16 6 roll\n"
      "    0 6 1 roll 17 6 roll\n"
      "    0 6 1 roll 18 6 roll\n"
      "  ]\n"
      "  gsave\n"
      "    newpath\n"
      "    <<\n"
      "      /DataSource 3 2 roll\n"
      "      /ShadingType 4\n"
      "      /ColorSpace [ /DeviceRGB ]\n"
      "    >>\n"
      "    shfill\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );
    
    /* For smooth quadrilateral shading we use {shfill} with {ShadingType = 7}. */
    /* The expected order of fields in the {DataSource} is: */
    /* [ 0 x00 y00 */
    /*     x10 y10 */
    /*     x20 y20 */
    /*     x30 y30 */
    /*     x31 y31 */
    /*     x32 y32 */
    /*     x33 y33 */
    /*     x23 y23 */
    /*     x13 y13 */
    /*     x03 y03 */
    /*     x02 y02 */
    /*     x01 y01 */
    /*     x11 y11 */
    /*     x21 y21 */
    /*     x22 y22 */
    /*     x12 y12 */
    /*     R00 G00 B00 */
    /*     R30 G30 B30 */
    /*     R33 G33 B33 */
    /*     R03 G03 B03 ] */
    
    /* The parameters to {gsquadf} and their indices on the PS stack are: */
    /*  x00 [19] */
    /*  y00 [18] */
    /*  R00 [17] */
    /*  G00 [16] */
    /*  B00 [15] */
    /*  x03 [14] */
    /*  y03 [13] */
    /*  R03 [12] */
    /*  G03 [11] */
    /*  B03 [10] */
    /*  x30 [09] */
    /*  y30 [08] */
    /*  R30 [07] */
    /*  G30 [06] */
    /*  B30 [05] */
    /*  x33 [04] */
    /*  y33 [03] */
    /*  R33 [02] */
    /*  G33 [01] */
    /*  B33 [00] */
    
    fprintf(wr, 
      "%% Bilinear Gouraud-shaded quadrilateral fill operator\n"
      "%%   (always fills, never strokes):\n"
      "%%   {x00} {y00} {R00} {G00} {B00} \\\n"
      "%%   {x03} {y03} {R03} {G03} {B03} \\\n"
      "%%   {x30} {y30} {R30} {G30} {B30} \\\n"
      "%%   {x33} {y33} {R33} {G33} {B33} gsquadf --> \n"
      "/ixm1\n"
      "{ index\n"
      "  21 1 roll\n"
      "} bind def\n"
      "\n"
      "/ixm21\n"
      "{ 1 add index exch\n"
      "  1 add index 2 mul\n"
      "  add 3 div\n"
      "  21 1 roll\n"
      "} bind def\n"
      "\n"
      "/ixm4221\n"
      "{ 3 add index 4 1 roll\n"
      "  3 add index 2 mul 4 1 roll\n"
      "  3 add index 2 mul 4 1 roll\n"
      "  3 add index 4 mul\n"
      "  add add add 9 div\n"
      "  21 1 roll\n"
      "} bind def\n"
      "\n"
      "/gsquadf\n"
      "{ [ 0 22 2 roll\n"
      "    19 ixm1 \n"     /* x00 = x00 */
      "    18 ixm1 \n"     /* y00 = y00 */
      "    19 09 ixm21 \n" /* x10 = (2*x00 + x30)/3 */
      "    18 08 ixm21 \n" /* y10 = (2*y00 + y30)/3 */
      "    09 19 ixm21 \n" /* x20 = (2*x30 + x00)/3 */
      "    08 18 ixm21 \n" /* y20 = (2*y30 + y00)/3 */
      "    09 ixm1 \n"     /* x30 = x30 */
      "    08 ixm1 \n"     /* y30 = y30 */
      "    09 04 ixm21 \n" /* x31 = (2*x30 + x33)/3 */
      "    08 03 ixm21 \n" /* y31 = (2*y30 + y33)/3  */
      "    04 09 ixm21 \n" /* x32 = (2*x33 + x30)/3 */
      "    03 08 ixm21 \n" /* y32 = (2*y33 + y30)/3 */
      "    04 ixm1 \n"     /* x33 = x33  */
      "    03 ixm1 \n"     /* y33 = y33  */
      "    04 14 ixm21 \n" /* x23 = (2*x33 + x03)/3  */
      "    03 13 ixm21 \n" /* y23 = (2*y33 + y03)/3  */
      "    14 04 ixm21 \n" /* x13 = (2*x03 + x33)/3 */
      "    13 03 ixm21 \n" /* y13 = (2*y03 + y33)/3 */
      "    14 ixm1 \n"     /* x03 = x03 */
      "    13 ixm1 \n"     /* y03 = y03 */
      "    14 19 ixm21 \n" /* x02 = (2*x03 + x00)/3 */
      "    13 18 ixm21 \n" /* y02 = (2*y03 + y00)/3 */
      "    19 14 ixm21 \n" /* x01 = (2*x00 + x03)/3 */
      "    18 13 ixm21 \n" /* y01 = (2*y00 + y03)/3 */
      "    19 14 09 04 ixm4221 \n" /* x11 = (4*x00 + 2*x03 + 2*x30 + x33)/9 */
      "    18 13 08 03 ixm4221 \n" /* y11 = (4*y00 + 2*y03 + 2*y30 + y33)/9 */
      "    09 19 04 14 ixm4221 \n" /* x21 = (4*x30 + 2*x00 + 2*x33 + x03)/9 */
      "    08 18 03 13 ixm4221 \n" /* y21 = (4*y30 + 2*y00 + 2*y33 + y03)/9 */
      "    04 14 09 19 ixm4221 \n" /* x22 = (4*x33 + 2*x03 + 2*x30 + x00)/9  */
      "    03 13 08 18 ixm4221 \n" /* y22 = (4*y33 + 2*y03 + 2*y30 + y00)/9 */
      "    14 19 04 09 ixm4221 \n" /* x12 = (4*x03 + 2*x00 + 2*x33 + x30)/9  */
      "    13 18 03 08 ixm4221 \n" /* y12 = (4*y03 + 2*y00 + 2*y33 + y30)/9  */
      "    17 ixm1 \n" /* R00 = R00  */
      "    16 ixm1 \n" /* G00 = G00  */
      "    15 ixm1 \n" /* B00 = B00  */
      "    07 ixm1 \n" /* R30 = R30  */
      "    06 ixm1 \n" /* G30 = G30  */
      "    05 ixm1 \n" /* B30 = B30 */
      "    02 ixm1 \n" /* R33 = R33  */
      "    01 ixm1 \n" /* G33 = G33  */
      "    00 ixm1 \n" /* B33 = B33  */
      "    12 ixm1 \n" /* R03 = R03  */
      "    11 ixm1 \n" /* G03 = G03  */
      "    10 ixm1 \n" /* B03 = B03  */
      "    20 { pop } repeat\n"
      "  ]\n"
      "  gsave\n"
      "    newpath\n"
      "    <<\n"
      "      /DataSource 3 2 roll\n"
      "      /ShadingType 7\n"
      "      /ColorSpace [ /DeviceRGB ]\n"
      "    >>\n"
      "    shfill\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );
    
    fprintf(wr, 
      "%% Rectangle operator:\n"
      "%%   {draw} {fill} {xlo} {xhi} {ylo} {yhi} rec --> \n"
      "/rec\n"
      "{ gsave\n"
      "  newpath\n"
      "    3 index 2 index moveto\n"
      "    2 index 2 index lineto\n"
      "    2 index 1 index lineto\n"
      "    3 index 1 index lineto\n"
      "    pop pop pop pop\n"
      "    closepath\n"
      "    fs\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );
    
    fprintf(wr, 
      "%% Polygon operator:\n"
      "%%   {draw} {fill} {eo} {cls} {p[n-1]} {p[n-2]} ... {p[0]} {n} pol --> \n"
      "%% where each {p[i]} is 2 coords. Note reverse order of points. \n"
      "/pol\n"
      "{ gsave\n"
      "    newpath\n"
      "    3 1 roll\n"
      "    %% --- draw fill eo cls p[n-1] ... p[1] n p[0]\n"
      "    moveto\n"
      "    %% --- draw fill eo cls p[n-1] ... p[1] n\n"
      "    1 sub 1 1  3 2 roll {\n"
      "      %% --- draw fill eo cls p[n-1] ... p[i] i\n"
      "      pop lineto\n"
      "    } for\n"
      "    %% --- draw fill eo cls\n"
      "    1 eq { closepath } if\n"
      "    %% --- draw fill eo\n"
      "    anyfs\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fprintf(wr, 
      "%% Rounded polygon operator:\n"
      "%%   {draw} {fill} {eo} {cls} {p[0]} {r[0]} {p[n-1]} {r[n-1]}... {p[0]} {r[0]} {n} {pm} cirpol --> \n"
      "%% where each {p[i]} and {pm} are 2 coords, each {r[i]} one radius. \n"
      "%% Note the reverse order of points. \n"
      "/cirpol\n"
      "{ gsave\n"
      "    newpath\n"
      "    moveto\n"
      "    %% --- draw fill eo cls p[0] r[0] p[n-1] r[n-1]... p[0] r[0] n\n"
      "    1 1  3 2 roll {\n"
      "      %% --- draw fill eo cls p[0] r[0] p[n-1] r[n-1]... p[i+1] r[i+1] p[i] r[i] i\n"
      "      pop 5 index 5 index\n"
      "      %% --- draw fill eo cls p[0] r[0] p[n-1] r[n-1]... p[i+1] r[i+1] p[i] r[i] p[i+1]\n"
      "      3 2 roll\n"
      "      %% --- draw fill eo cls p[0] r[0] p[n-1] r[n-1]... p[i+1] r[i+1] p[i] p[i+1] r[i]\n"
      "      dup 0 eq\n"
      "        { pop pop pop lineto }\n"
      "        { arcto pop pop pop pop }\n"
      "      ifelse\n"
      "    } for\n"
      "    %% --- draw fill eo cls p[0] r[0]\n"
      "    pop pop pop\n"
      "    %% --- draw fill eo cls\n"
      "    1 eq { closepath } if\n"
      "    %% --- draw fill eo\n"
      "    anyfs\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fprintf(wr, 
      "%% Bezier polygon operator:\n"
      "%%   {draw} {fill} {eo} {cls} {a[n-1]} ... {a[0]} {n} bzpol --> \n"
      "%% where each {a[i]} is 4 points, each 2 coords: \n"
      "%%   {p[1]} {p[2]} {p[3]} {p[0]}.\n"
      "%% Note the reverse order of arcs and the funny order of the 4 points. \n"
      "/bzpol\n"
      "{ gsave\n"
      "    newpath\n"
      "    3 1 roll\n"
      "    %% --- draw fill eo cls a[n-1] ... a[1] p[1] p[2] p[3] n p[0]\n"
      "    moveto\n"
      "    %% --- draw fill eo cls a[n-1] ... a[1] p[1] p[2] p[3] n\n"
      "    7 1 roll\n"
      "    %% --- draw fill eo cls a[n-1] ... a[1] n p[1] p[2] p[3]\n"
      "    curveto\n"
      "    %% --- draw fill eo cls a[n-1] ... a[1] n\n"
      "    1 sub 1 1  3 2 roll {\n"
      "      %% --- draw fill eo cls a[n-1] ... a[i] i\n"
      "      pop\n"
      "      lineto\n"
      "      curveto\n"
      "    } for\n"
      "    %% --- draw fill eo cls\n"
      "    1 eq { closepath } if\n"
      "    %% --- draw fill eo\n"
      "    anyfs\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fprintf(wr, 
      "%% Circle operator:\n"
      "%%   {draw} {fill} {x} {y} {rad} cir --> \n"
      "/cir\n"
      "{ gsave\n"
      "    newpath\n"
      "    0 360 arc\n"
      "    closepath\n"
      "    fs\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );
    
    fprintf(wr, 
      "%% Lune (lens) operator:\n"
      "%%   {draw} {fill} {xc} {yc} {rad} {tilt} lun --> \n"
      "/lun\n"
      "{ gsave\n"
      "    4 2 roll\n"
      "    %% --- draw fill rad tilt xc yc\n"
      "    translate\n"
      "    %% --- draw fill rad tilt\n"
      "    rotate\n"
      "    %% --- draw fill rad \n"
      "    newpath\n"
      "    dup\n"
      "    %% --- draw fill rad rad\n"
      "    dup neg 0\n"
      "    %% --- draw fill rad rad -rad 0\n"
      "    3 2 roll 2 mul\n"
      "    %% --- draw fill rad -rad 0  2*rad\n"
      "    -60 60 arc\n"
      "    %% --- draw fill rad\n"
      "    dup 0\n"
      "    %% --- draw fill rad rad 0\n"
      "    3 2 roll 2 mul\n"
      "    %% --- draw fill rad 0 2*rad\n"
      "    120 240 arc\n"
      "    closepath\n"
      "    fs\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fprintf(wr, 
      "%% Pie slice operator:\n"
      "%%   {draw} {fill} {xc} {yc} {rad} {start} {stop} pie --> \n"
      "/pie\n"
      "{ gsave\n"
      "    newpath\n"
      "    4 index 4 index moveto arc\n"
      "    closepath\n"
      "    fs\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fprintf(wr, 
      "%% Grid cell operator:\n"
      "%%   {draw} {fill} {ih} {nh} {iv} {nv} cel --> \n"
      "/cel\n"
      "{ 1 index 1 add 1 index exch\n"      /* --- draw fill ih nh iv nv nv iv+1 */
      "  vMax vMin sub mul\n"               /* --- draw fill ih nh iv nv (iv+1)*(vMax-vMin) */
      "  exch div vMin add\n"               /* --- draw fill ih nh iv nv vhi */
      "  3 1 roll exch\n"                   /* --- draw fill ih nh vhi nv iv */
      "  vMax vMin sub mul\n"               /* --- draw fill ih nh vhi nv iv*(vMax-vMin) */
      "  exch div vMin add\n"               /* --- draw fill ih nh vhi vlo */
      "  exch 4 2 roll\n"                   /* --- draw fill vlo vhi ih nh */
      "\n"
      "  1 index 1 add 1 index exch\n"      /* --- draw fill vlo vhi ih nh nh ih+1 */
      "  hMax hMin sub mul\n"               /* --- draw fill vlo vhi ih nh nh (ih+1)*(hMax-hMin) */
      "  exch div hMin add\n"               /* --- draw fill vlo vhi ih nh hhi */
      "  3 1 roll exch\n"                   /* --- draw fill vlo vhi hhi nh ih */
      "  hMax hMin sub mul\n"               /* --- draw fill vlo vhi hhi nh ih*(hMax-hMin) */
      "  exch div hMin add\n"               /* --- draw fill vlo vhi hhi hlo */
      "  exch 4 2 roll\n"                   /* --- draw fill hlo hhi vlo vhi */
      "  rec\n"
      "} bind def\n"
      "\n"
    );

    /* Labels and text: */
    
    /* The {strfs} and {ncstrfs} operators uses the current font and font size. */
    
    fprintf(wr, 
      "%% Generic clipped label/text draw and fill operator:\n"
      "%%   {draw} {fill} {str} {str1} {xa} {ya} {rot} {xc} {yc} strfs --> \n"
      "/strfs\n"
      "{ newpath moveto\n"
      "    %% --- draw, fill, str, str1, xa, ya, rot\n"
      "  5 1 roll\n"
      "    %% --- draw, fill, rot, str, str1, xa, ya\n"
      "  gsave 3 index false charpath flattenpath pathbbox grestore\n"
      "    %% --- draw, fill, rot, str, str1, xa, ya, lox, loy, hix, hiy\n"
      "  pop exch pop\n"
      "    %% --- draw, fill, rot, str, str1, xa, ya, lox, hix\n"
      "  gsave 5 -1 roll false charpath flattenpath pathbbox grestore\n"
      "    %% --- draw, fill, rot, str, xa, ya, lox, hix lox1 loy1 hix1 hiy1\n"
      "  exch pop 3 -1 roll pop\n"
      "    %% --- draw, fill, rot, str, xa, ya, lox, hix loy1 hiy1\n"
      "  3 1 roll exch 3 -1 roll\n"
      "    %% --- draw, fill, rot, str, xa, ya, lox, loy1 hix hiy1\n"
      "  \n"
      "  3 index 3 index currentpoint\n"
      "    %% --- draw, fill, rot, str, xa, ya, lox, loy1, hix, hiy1, lox, loy1, cx, cy\n"
      "  exch\n"
      "    %% --- draw, fill, rot, str, xa, ya, lox, loy1, hix, hiy1, lox, loy1, cy, cx\n"
      "  4 1 roll\n"
      "    %% --- draw, fill, rot, str, xa, ya, lox, loy1, hix, hiy1, cx, lox, loy1, cy\n"
      "  exch\n"
      "    %% --- draw, fill, rot, str, xa, ya, lox, loy1, hix, hiy1, cx, lox, cy, loy1\n"
      "  sub\n"
      "    %% --- draw, fill, rot, str, xa, ya, lox, loy1, hix, hiy1, cx, lox, cy-loy1\n"
      "  3 1 roll\n"
      "    %% --- draw, fill, rot, str, xa, ya, lox, loy1, hix, hiy1, cy-loy1, cx, lox\n"
      "  sub\n"
      "    %% --- draw, fill, rot, str, xa, ya, lox, loy1, hix, hiy1, cy-loy1, cx-lox\n"
      "  exch\n"
      "    %% --- draw, fill, rot, str, xa, ya, lox, loy1, hix, hiy1, cx-lox, cy-loy1\n"
      "  10 -1 roll (C) pstack pop\n"
      "    %% --- draw, fill, str, xa, ya, lox, loy, hix1, hiy1, cx-lox, cy-loy1, rot\n"
      "  gsave\n"
      "  rotate\n"
      "    %% --- draw, fill, str, xa, ya, lox, loy1, hix, hiy1, cx-lox, cy-loy1\n"
      "  rmoveto\n"
      "    %% --- draw, fill, str, xa, ya, lox, loy1, hix, hiy1\n"
      "  exch\n"
      "    %% --- draw, fill, str, xa, ya, lox, loy1, hiy1, hix\n"
      "  4 1 roll\n"
      "    %% --- draw, fill, str, xa, ya, hix, lox, loy1, hiy1\n"
      "  exch\n"
      "    %% --- draw, fill, str, xa, ya, hix, lox, hiy1, loy1\n"
      "  sub \n"
      "    %% --- draw, fill, str, xa, ya, hix, lox, dy1\n"
      "  3 1 roll\n"
      "    %% --- draw, fill, str, xa, ya, dy1, hix, lox\n"
      "  sub\n"
      "    %% --- draw, fill, str, xa, ya, dy1, dx\n"
      "  exch (D) pstack pop\n"
      "    %% --- draw, fill, str, xa, ya, dx, dy1\n"
      "  exch 4 1 roll mul -1 mul\n"
      "  3 1 roll mul -1 mul exch (E) pstack pop\n"
      "    %% --- draw, fill, str, -dx*xa, -dy1*ya\n"
      "  rmoveto\n"
      "    %% --- draw, fill, str\n"
      "  true charpath\n"
      "    %% --- draw, fill\n"
      "  fs\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fprintf(wr, 
      "%% Generic unclipped label/text draw and fill operator:\n"
      "%%   {draw} {fill} {str} {str1} {xa} {ya} {rot} {xc} {yc} ncstrfs --> \n"
      "/ncstrfs\n"
      "{ gsave\n"
      "    initclip\n"
      "    strfs\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fflush(wr);
  }

epswr_figure_t *epswr_dev_new_figure
  ( FILE *wr,
    double hSize,   /* Figure width (in pt). */
    double vSize,   /* Figure height (in pt). */
    bool_t verbose  /* TRUE to print diagnostics. */
  )
  {
    epswr_figure_t *eps = (epswr_figure_t *)notnull(malloc(sizeof(epswr_figure_t)), "out of mem");
    
    eps->verbose = verbose;
    eps->wr = wr;

    /* Set the total figure size: */
    epswr_check_param("hSize", hSize, 1.0, epswr_MAX_SIZE);
    epswr_check_param("vSize", vSize, 1.0, epswr_MAX_SIZE);

    if (eps->verbose)
      { fprintf(stderr, "epswr_dev_new_figure:\n");
        fprintf(stderr, "  total figure size = %6.1f × %6.1f (pt)\n", hSize, vSize);
      }

    eps->hSize = hSize;
    eps->vSize = vSize;
    
    /* Write the Postscript preamble and definitions: */
    epswr_dev_write_file_header(eps->wr, eps->hSize, eps->vSize);

    /* Clear font list: */
    eps->nFonts = 0;
    eps->fonts = NULL;
    epswr_font_list_initialize(&(eps->nFonts), &(eps->fonts));

    /* Initialize the Device coords of the plot windows: */
    epswr_dev_set_window(eps, 0, hSize, 0, vSize, FALSE);

    /* Initialize the plotting state, just in case the client forgets to set it up: */
    epswr_dev_set_fill_color(eps, 0.500,0.500,0.500);
    epswr_dev_set_pen(eps, 0.000,0.000,0.000, 1.0, 0.0, 0.0);
    epswr_dev_set_label_font(eps, "Courier", 10.0); 
    epswr_dev_set_text_font(eps, "Courier", 14.0); 
    epswr_dev_set_text_geometry(eps, 0, hSize, 0, hSize, 0.0); 
    
    fflush(eps->wr);
    
    return eps;
  }

void epswr_dev_end_figure(epswr_figure_t *eps)
  { 
    demand(eps->wr != NULL, "figure already closed?");
    epswr_dev_write_file_trailer(eps->wr, eps->nFonts, eps->fonts);
    epswr_font_list_free(&(eps->nFonts), &(eps->fonts));
    fclose(eps->wr);
    eps->wr = NULL;
    free(eps);
  }

void epswr_dev_get_figure_size
  ( epswr_figure_t *eps,  /* Picture stream. */
    double *hSizeP,        /* OUT: Total width of figure (in pt). */
    double *vSizeP         /* OUT: Total height of figure (in pt). */
  )
  { *hSizeP = eps->hSize;
    *vSizeP = eps->vSize;
  }

void epswr_dev_set_window
  ( epswr_figure_t *eps,
    double hMin, double hMax,
    double vMin, double vMax,
    bool_t relative
  )
  { if (relative) 
      { /* Displace the given window by {eps->(hMin,vMin)}: */
        hMin += eps->hMin;  hMax += eps->hMin;
        vMin += eps->vMin;  vMax += eps->vMin;
      }
      
    demand((0 <= hMin) && (hMin < hMax) && (hMax <= eps->hSize), "invalid {hMin,hMax}");
    demand((0 <= vMin) && (vMin < vMax) && (vMax <= eps->vSize), "invalid {vMin,vMax}");
    
    /* Save the window parameters in the {epswr_figure_t} record: */
    epswr_def_set_device_window(eps, hMin, hMax, vMin, vMax, "epswr_dev_set_window");
    epswr_def_undefine_client_window(eps);
    
    /* Write the window setup commands to the Postcript file: */
    epswr_dev_write_window_set_cmds(eps->wr, hMin, hMax, vMin, vMax);
    fflush(eps->wr);
  }
  
void epswr_dev_get_window
  ( epswr_figure_t *eps,
    double *hMinP, double *hMaxP,
    double *vMinP, double *vMaxP
  )
  {
    (*hMinP) = eps->hMin;
    (*hMaxP) = eps->hMax;
    (*vMinP) = eps->vMin;
    (*vMaxP) = eps->vMax;
  }

void epswr_dev_shrink_window
  ( epswr_figure_t *eps, 
    double dhMin, double dhMax, 
    double dvMin, double dvMax
  )
  { double hMin = eps->hMin + dhMin;
    double hMax = eps->hMax - dhMax;
    double vMin = eps->vMin + dvMin;
    double vMax = eps->vMax - dvMax;
    epswr_dev_set_window(eps, hMin, hMax, vMin, vMax, FALSE);
  }

void epswr_dev_set_window_to_grid_cell
  ( epswr_figure_t *eps, 
    double hMin, double hMax, int32_t ih, int32_t nh, 
    double vMin, double vMax, int32_t iv, int32_t nv
  )
  { demand((0 <= ih) && (ih < nh), "invalid grid column spec {ih,nh}");
    demand((0 <= iv) && (iv < nv), "invalid grid column spec {iv,nv}");
    double rhMin = ((double)ih)/((double)nh), rhMax = ((double)ih+1)/((double)nh);
    double hMinCell = (1 - rhMin)*hMin + rhMin*hMax;
    double hMaxCell = (1 - rhMax)*hMin + rhMax*hMax;
    double rvMin = ((double)iv)/((double)nv), rvMax = ((double)iv+1)/((double)nv);
    double vMinCell = (1 - rvMin)*vMin + rvMin*vMax;
    double vMaxCell = (1 - rvMax)*vMin + rvMax*vMax;
    epswr_dev_set_window(eps, hMinCell, hMaxCell, vMinCell, vMaxCell, FALSE);
  }

void epswr_dev_set_pen
  ( epswr_figure_t *eps,
    double R, double G, double B,
    double pswidth,
    double psdashLength,
    double psdashSpace
  )
  { FILE *wr = eps->wr;
    demand((!isnan(R)) && (!isnan(G)) && (!isnan(B)), "invalid pen color");
    if (R < 0.0) { R = 0.0; } else if (R > 1.0) { R = 1.0; }
    if (G < 0.0) { G = 0.0; } else if (G > 1.0) { G = 1.0; }
    if (B < 0.0) { B = 0.0; } else if (B > 1.0) { B = 1.0; }
    fprintf(wr, "%5.3f %5.3f %5.3f setrgbcolor\n", R, G, B);
    fprintf(wr, "%.3f setlinewidth\n", pswidth);
    if ((psdashLength == 0.0) | (psdashSpace == 0.0))
      { fprintf(wr, " [ ] 0 setdash\n"); }
    else
      { fprintf(wr,
          " [ %.3f %.3f ] 0 setdash\n",
          psdashLength, psdashSpace
        );
      }
  }

void epswr_dev_segment
  ( epswr_figure_t *eps,
    double psxa, double psya,
    double psxb, double psyb
  )
  { if (eps->verbose) 
      { fprintf(stderr, "segment: (%.3f %.3f) --> (%.3f %.3f)\n",  psxa, psya, psxb, psyb); }
    if (epswr_segment_is_invisible(eps, psxa, psya, psxb, psyb)) 
      { return; }
    FILE *wr = eps->wr;
    fprintf(wr,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxa, psya, psxb, psyb
    );
  }
  
void epswr_dev_curve
  ( epswr_figure_t *eps,
    double psxa, double psya,
    double psxb, double psyb,
    double psxc, double psyc,
    double psxd, double psyd
  )
  { if (epswr_curve_is_invisible(eps, psxa, psya, psxb, psyb, psxc, psyc, psxd, psyd))
      { return; }
    FILE *wr = eps->wr;
    fprintf(wr, "%6.1f %6.1f  %6.1f %6.1f  %6.1f %6.1f  %6.1f %6.1f arcd\n",
      psxa, psya, psxb, psyb, psxc, psyc, psxd, psyd
    );
  }

void epswr_dev_coord_line
  ( epswr_figure_t *eps, 
    epswr_axis_t axis, 
    double pspos
  )
  { FILE *wr = eps->wr;
    fprintf(wr, "%6.1f %sgrd\n", pspos, (axis == epswr_axis_HOR ? "x" : "y"));
  }

void epswr_dev_frame (epswr_figure_t *eps)
  { FILE *wr = eps->wr;
    fprintf(wr, "wframe\n");
  }

void epswr_dev_set_fill_color(epswr_figure_t *eps, double R, double G, double B)
  { if (isnan(R) || (fabs(R) == INFINITY) || (R < 0.0)) { R = G = B = -1.0; }
    double *fc = eps->fillColor;
    if ((R != fc[0]) || (G != fc[1]) || (B != fc[2]))
      { FILE *wr = eps->wr;
        fc[0] = R; fc[1] = G; fc[2] = B;
        epswr_dev_write_fill_color_set_cmds(wr, fc); 
      }
  }

void epswr_dev_rectangle
  ( epswr_figure_t *eps,
    double psxlo, double psxhi,
    double psylo, double psyhi,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    if (eps->verbose) 
      { fprintf(stderr, "rectangle: [ %.3f _ %.3f] x [%.3f _ %.3f]\n",  psxlo, psxhi, psylo, psyhi); }
    if (epswr_rectangle_is_invisible(eps, psxlo, psxhi, psylo, psyhi))
      { return; }
    FILE *wr = eps->wr;
    fprintf(wr, "%d %d %6.1f %6.1f  %6.1f %6.1f",  draw, fill, psxlo, psxhi, psylo, psyhi);
    fprintf(wr, " rec\n");
    fflush(wr);
  }
  
void epswr_dev_triangle
  ( epswr_figure_t *eps,
    double psxa, double psya,
    double psxb, double psyb,
    double psxc, double psyc,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    if (epswr_triangle_is_invisible(eps, psxa, psya, psxb, psyb, psxc, psyc))
      { return; }
    FILE *wr = eps->wr;
    fprintf(wr,
      "%d %d %6.1f %6.1f  %6.1f %6.1f  %6.1f %6.1f ",
      draw, fill, psxa, psya, psxb, psyb, psxc, psyc
    );
    fprintf(wr, " tri\n");
    fflush(wr);
  }
  
void epswr_dev_quadrilateral
  ( epswr_figure_t *eps,
    double psx00, double psy00,
    double psx01, double psy01,
    double psx10, double psy10,
    double psx11, double psy11,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psx[4], psy[4];
    psx[0] = psx00; psy[0] = psy00;
    psx[1] = psx01; psy[1] = psy01;
    psx[2] = psx11; psy[2] = psy11;
    psx[3] = psx10; psy[3] = psy10;
    epswr_dev_polygon(eps, TRUE, psx, psy, 4, fill, draw, TRUE);
  }
  
/* !!! Polygons should do the right thing with dashes. !!! */
   
void epswr_dev_polygon
  ( epswr_figure_t *eps,
    bool_t closed,
    double psx[], double psy[],
    int32_t n,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if (! closed) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    if (! epswr_polygon_is_invisible(eps, psx, psy, n))
      { /* Plot them: */
        FILE *wr = eps->wr;
        fprintf(wr, "%d %d %d %d", draw, fill, evenOdd, closed);
        if (n > 6) { fprintf(wr, "\n"); }
        /* Write the sides in the reverse order: */
        int32_t nplin = 0; /* Number of points in current line. */
        for (int32_t i = n-1; i >= 0; i--)
          { if (nplin >= 6) { fprintf(wr, "\n"); nplin = 0; }
            fprintf(wr, "  %6.1f %6.1f", psx[i], psy[i]);
          }
        fprintf(wr, "  %d", n);
        fprintf(wr, " %s\n", "pol");
        fflush(wr);
      }
  }
  
void epswr_dev_rounded_polygon
  ( epswr_figure_t *eps,
    bool_t closed,
    double psx[], double psy[],
    int32_t n,
    double psrad,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if (! closed) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    if (! epswr_polygon_is_invisible(eps, psx, psy, n))
      { /* Compute the radii to use at each corner: */
        double *arad = (double *)malloc(n*sizeof(double));
        for (int32_t i = 0; i<n; i++)
          { int32_t j = (i + 1) % n;
            int32_t k = (i + 2) % n;
            /* Adjust the rounding radius at corner {j = i+1}: */
            if ((! closed) && ((j == 0) || (j == n-1)))
              { arad[j] = 0; }
            else
              { epswr_dev_compute_joining_arc_radius
                  ( psx[i],psy[i], psx[j],psy[j], psx[k],psy[k], psrad, &(arad[j]) );
              }
          }
        /* Output the plot command: */
        FILE *wr = eps->wr;
        fprintf(wr, "%d %d %d %d", draw, fill, evenOdd, closed);
        if (n > 3) { fprintf(wr, "\n"); }
        /* Write the corners and radii in reverse order: */
        int32_t nplin = 0; /* Number of points in current line. */
        for (int32_t ii = n; ii >= 0; ii--)
          { if (nplin >= 3) { fprintf(wr, "\n"); nplin = 0; }
            int32_t i = ii % n;
            fprintf(wr, "  %6.1f %6.1f %8.3f", psx[i], psy[i], arad[i]);
            nplin++;
          }
        fprintf(wr, "\n");
        /* Write the number of points and the starting point: */
        fprintf(wr, "  %d", n);
        double pmdx = (closed ? (psx[0]+psx[n-1])/2 : psx[0]); 
        double pmdy = (closed ? (psy[0]+psy[n-1])/2 : psy[0]); 
        fprintf(wr, "  %6.1f %6.1f", pmdx, pmdy);
        fprintf(wr, " %s\n", "cirpol");
        fflush(wr);
        free(arad);
      }
  }
    
void epswr_dev_bezier_polygon
  ( epswr_figure_t *eps,
    bool_t closed,
    double psx[], double psy[],
    int32_t n,
    bool_t fill, bool_t draw,
    bool_t evenOdd
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if (! closed) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    int32_t np = 4*n; /* Number of points. */
    /* 
      The following test assumes that if the straight polygon with
      those {np} vertices is invisible, the Bézier polygon is
      invisible too. This is true for simple visibility tests (such as
      bounding box), but not for more precise ones. To be correct in
      any case, we should checl the convex hull of those points
      instead.
    */
    if (! epswr_polygon_is_invisible(eps, psx, psy, np))
      { /* Plot it: */
        FILE *wr = eps->wr;
        fprintf(wr, "%d %d %d %d", draw, fill, evenOdd, closed);
        fprintf(wr, "\n");
        /* Write the arcs in the reverse order: */
        for (int32_t i = n-1; i >= 0; i--)
          { int32_t k0 = 4*i;              /* Start of arc number {i} */
            for (int32_t j = 0; j < 4; j++)
              { int32_t j1 = (j + 1) % 4;
                double psxi = psx[k0+j1];
                double psyi = psy[k0+j1];
                fprintf(wr, "  %6.1f %6.1f", psxi, psyi);
              }
            fprintf(wr, "\n");
          }
        fprintf(wr, "  %d", n);
        fprintf(wr, " %s\n", "bzpol");
        fflush(wr);
      }
  }
  
void epswr_dev_circle
  ( epswr_figure_t *eps,
    double psxc, double psyc, double psrad,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    if (epswr_circle_is_invisible(eps, psxc, psyc, psrad))
      { return; }
    FILE *wr = eps->wr;
    fprintf(wr, "%d %d %6.1f %6.1f  %6.2f", draw, fill, psxc, psyc, psrad); 
    fprintf(wr, " cir\n");
    fflush(wr);
  }

void epswr_dev_lune
  ( epswr_figure_t *eps,
    double psxc, double psyc, double psrad, double pstilt,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    if (epswr_lune_is_invisible(eps, psxc, psyc, psrad, pstilt))
      { return; }
    FILE *wr = eps->wr;
    fprintf(wr, "%d %d %6.1f %6.1f  %6.1f %6.2f",
      draw, fill, psxc, psyc, psrad, pstilt); 
    fprintf(wr, " lun\n");
    fflush(wr);
  }
  
void epswr_dev_slice
  ( epswr_figure_t *eps,
    double psxc, double psyc, double psrad,
    double start, double stop,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    if (epswr_slice_is_invisible(eps, psxc, psyc, psrad, start, stop))
      { return; }
    FILE *wr = eps->wr;
    fprintf(wr, "%d %d %6.1f %6.1f  %6.1f  %6.2f %6.2f",
      draw, fill, psxc, psyc, psrad, start, stop); 
    fprintf(wr, " pie\n");
    fflush(wr);
  }

void epswr_dev_dot
  ( epswr_figure_t *eps,
    double psxc, double psyc, double psrad,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    if (epswr_circle_is_invisible(eps, psxc, psyc, psrad))
      { return; }
    FILE *wr = eps->wr;
    fprintf(wr, "%d %d %6.1f %6.1f  %6.1f", draw, fill, psxc, psyc, psrad); 
    fprintf(wr, " cir\n");
    fflush(wr);
    
  }
  
void epswr_dev_tic
  ( epswr_figure_t *eps, 
    epswr_axis_t axis, 
    double psxc, double psyc, 
    double psticSize,
    double align 
  )
  {
    double psxa, psya, psxb, psyb;
    if (axis == epswr_axis_HOR)
      { psxa = psxb = psxc; 
        psya = psyc - align*psticSize; 
        psyb = psyc + (1-align)*psticSize; 
      }
    else if (axis == epswr_axis_VER)
      { psxa = psxc - align*psticSize; 
        psxb = psxc + (1-align)*psticSize;
        psya = psyb = psyc; 
      }
    else
      { affirm(FALSE, "invalid axis"); psxa = psxb = psya = psyb = 0.0; }
    epswr_dev_segment(eps, psxa, psya, psxb, psyb);
  }
  
void epswr_dev_cross
  ( epswr_figure_t *eps, 
    double psxc, double psyc, double psrad, bool_t diag,
    bool_t draw
  )
  { if (! draw) { return; }
    if (epswr_rectangle_is_invisible(eps, psxc-psrad, psxc+psrad, psyc-psrad, psyc+psrad))
      { return; }
    double psxA, psyA, psxB, psyB, psxC, psyC, psxD, psyD;
    if (diag)
      { double d = psrad*0.7071067812;
        psxA = psxc - d; psyA = psyc - d;
        psxB = psxc + d; psyB = psyc + d;
        psxC = psxc - d; psyC = psyc + d;
        psxD = psxc + d; psyD = psyc - d;
      }
    else
      { psxA = psxc - psrad; psyA = psyc;
        psxB = psxc + psrad; psyB = psyc;
        psxC = psxc;         psyC = psyc - psrad;
        psxD = psxc;         psyD = psyc + psrad;
      }
    FILE *wr = eps->wr;
    fprintf(wr,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxA, psyA, psxB, psyB
    );
    fprintf(wr,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxC, psyC, psxD, psyD
    );
  }
  
void epswr_dev_asterisk
  ( epswr_figure_t *eps, 
    double psxc, double psyc, double psrad,
    bool_t draw
  )
  { if (!draw) { return; }
    double ct = 0.92387953 * psrad;
    double st = 0.38268343 * psrad;
    if (epswr_rectangle_is_invisible(eps, psxc-ct, psxc+ct, psyc-ct, psyc+ct))
      { return; }
    double psxA = psxc - st; double psyA = psyc - ct;
    double psxB = psxc + st; double psyB = psyc + ct;
    double psxC = psxc - ct; double psyC = psyc - st;
    double psxD = psxc + ct; double psyD = psyc + st;
    FILE *wr = eps->wr;
    fprintf(wr,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxA, psyA, psxB, psyB
    );
    fprintf(wr,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxA, psyB, psxB, psyA
    );
    fprintf(wr,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxC, psyC, psxD, psyD
    );
    fprintf(wr,
      "%6.1f %6.1f  %6.1f %6.1f segd\n",
      psxC, psyD, psxD, psyC
    );
  }

void epswr_dev_square
  ( epswr_figure_t *eps,
    double psxc, double psyc, double psrad,
    bool_t fill, bool_t draw
  )
  { double d = 0.7071067812 * psrad;
    double psxlo = psxc - d;
    double psxhi = psxc + d;
    double psylo = psyc - d;
    double psyhi = psyc + d;
    epswr_dev_rectangle(eps, psxlo, psxhi, psylo, psyhi, fill, draw);
  }

void epswr_dev_diamond
  ( epswr_figure_t *eps, 
    double psxc, double psyc,
    double psxRad, double psyRad,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    double psxlo = psxc - psxRad;
    double psxhi = psxc + psxRad;
    double psylo = psyc - psyRad;
    double psyhi = psyc + psyRad;
    if (epswr_rectangle_is_invisible(eps, psxlo, psxhi, psylo, psyhi))
      { return; }
    FILE *wr = eps->wr;
    bool_t evenOdd = TRUE;
    bool_t closed = TRUE;
    fprintf(wr, "%d %d %d %d", draw, fill, evenOdd, closed);
    /* Write the sides in the reverse order: */
    fprintf(wr, "  %6.1f %6.1f", psxlo, psyc); 
    fprintf(wr, "  %6.1f %6.1f", psxc,  psylo); 
    fprintf(wr, "  %6.1f %6.1f", psxhi, psyc); 
    fprintf(wr, "  %6.1f %6.1f", psxc,  psyhi); 
    fprintf(wr, "  %d", 4);
    fprintf(wr, " %s\n", "pol");
    fflush(wr);
  }
    
void epswr_dev_arrowhead 
  ( epswr_figure_t *eps,
    double psxa, double psya, double psxb, double psyb,
    double pswidth, double pslength, 
    double fraction,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    
    /* Unit direction vector: */
    double dx = psxb - psxa;
    double dy = psyb - psya;
    double d = sqrt(dx*dx + dy*dy);
    dx /= d; dy /= d;
    
    /* Arrow dimensions in Device units: */
    double psw = pswidth/2.0;
    double psh = pslength;
    
    /* Corners of triangle: */
    double noitcarf = 1.0 - fraction;
    double psxt = noitcarf*psxa + fraction*psxb;
    double psyt = noitcarf*psya + fraction*psyb;
    double psxu = psxt - psh * dx + psw * dy;
    double psyu = psyt - psh * dy - psw * dx;
    double psxv = psxt - psh * dx - psw * dy;
    double psyv = psyt - psh * dy + psw * dx;
    
    if (epswr_triangle_is_invisible(eps, psxt, psyt, psxu, psyu, psxv, psyv))
      { return; }
    FILE *wr = eps->wr;
    fprintf(wr, "%d %d %6.1f %6.1f  %6.1f %6.1f  %6.1f %6.1f",
      draw, fill, psxt, psyt, psxu, psyu, psxv, psyv); 
    fprintf(wr, " tri\n");
    fflush(wr);
  }
    
/* GRID LINES AND GRID CELLS */

void epswr_dev_grid_lines(epswr_figure_t *eps, int32_t nh, int32_t nv)
  {
    double hMin, hMax, vMin, vMax;
    epswr_dev_get_window(eps, &hMin, &hMax, &vMin, &vMax);
    FILE *wr = eps->wr;
    for (int32_t ih = 0; ih<=nh; ih++)
      { double r = ((double)ih)/((double)nh);
        double h = (1 - r)*hMin + r*hMax;
        fprintf(wr, "%6.1f xgrd\n", h);
      }
    for (int32_t iv = 0; iv<=nv; iv++)
      { double r = ((double)iv)/((double)nv);
        double v = (1 - r)*vMin + r*vMax;
        fprintf(wr, "%6.1f ygrd\n", v);
      }
  }

void epswr_dev_grid_cell
  ( epswr_figure_t *eps, 
    int32_t ih, int32_t nh,
    int32_t iv, int32_t nv,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    FILE *wr = eps->wr;
    fprintf(wr, "%d %d   %3d %3d  %3d %3d", draw, fill, ih, nh, iv, nv);
    fprintf(wr, " cel\n");
    fflush(wr);
  }

void epswr_dev_set_label_font(epswr_figure_t *eps, const char *font, double size)
  { epswr_font_list_add(font, &(eps->nFonts), &(eps->fonts));
    eps->labelFont = font;
    eps->labelFontSize = size;
    FILE *wr = eps->wr;
    epswr_dev_write_label_font_set_cmds(wr, font, size);
  }

void epswr_dev_label
  ( epswr_figure_t *eps, 
    const char *text, 
    const char *strut, 
    double psx, double psy, 
    double rot,
    bool_t clipped,
    double hAlign, double vAlign,
    bool_t fill, bool_t draw
  )
  { if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    if (eps->verbose)
      { fprintf(stderr, "label: \"%s\" at (%.3f,%.3f) rot = %.1f\n", text, psx, psy, rot); }

    FILE *wr = eps->wr;
    fprintf(wr, "labelFont setfont\n");
    epswr_dev_text_line(eps, text, strut, psx, psy, rot, clipped, hAlign, vAlign, fill, draw);
  }

/* RUNNING TEXT */

void epswr_dev_set_text_geometry
  ( epswr_figure_t *eps, 
    double hMin, double hMax, 
    double vMin, double vMax,
    double rot
  )
  {
    eps->hCtrText = (hMin + hMax)/2;
    eps->vCtrText = (vMin + vMax)/2;
    eps->hSizeText = hMax - hMin;
    eps->rotText = rot;
    eps->vTopText = vMax - eps->vCtrText;
  }

void epswr_dev_set_text_font(epswr_figure_t *eps, const char *font, double size)
  { epswr_font_list_add(font, &(eps->nFonts), &(eps->fonts));
    eps->textFont = font;
    eps->textFontSize = size;
    FILE *wr = eps->wr;
    epswr_dev_write_text_font_set_cmds(wr, font, size);
  }

void epswr_dev_text
  ( epswr_figure_t *eps, 
    const char *text, 
    bool_t clipped,
    double hAlign, 
    bool_t fill, bool_t draw
  )
  { 
    bool_t debug = FALSE;
    if (eps->fillColor[0] < 0.0) { fill = FALSE; }
    if ((!draw) && (!fill)) { return; }
    
    FILE *wr = eps->wr;
    fprintf(wr, "textFont setfont\n");
    /* Loop on text lines */
    double dv = eps->textFontSize; /* Assumed baseline spacing. */
    double rot = eps->rotText;
    double ang = rot*M_PI/180; /* Rotation angle in radians. */
    double ca = cos(ang), sa = sin(ang);
    double vAlign = 1.00;
    char *line = (char *)text; /* Start of next line of text. */
    while(TRUE)
      { /* Find the end {pend} of the line: */
        char *pend = line;
        while (((*pend) != 0) && ((*pend) != '\n') && ((*pend) != '\r')) { pend++; }
        /* Compute the vert coord of the line bottom in the unrotated text area, rel to center: */
        if (pend > line)
          { /* !!! Should account for font depth. !!! */
            /* Compute the vert coord of text bottom relto center, unrotated. */
            double vRefRaw = eps->vTopText; 
            /* Compute the horiz coord of the line's ref point rel to center, unrotated: */
            double hRefRaw = (hAlign - 0.5)*eps->hSizeText;
            if (debug) { fprintf(stderr, "raw = ( %.3f %.3f )", hRefRaw, vRefRaw); }
            /* Rotate the ref point around the center: */
            double hRef = eps->hCtrText + ca*hRefRaw - sa*vRefRaw;
            double vRef = eps->vCtrText + sa*hRefRaw + ca*vRefRaw;
            if (debug) { fprintf(stderr, " rot = ( %.3f %.3f )\n", hRef, vRef); }
            if (debug) { epswr_dev_dot(eps, hRef, vRef, 0.5, TRUE, FALSE); }
            epswr_dev_text_line(eps, line, "Rg", hRef, vRef, rot, clipped, hAlign, vAlign, fill, draw);
          }
        /* Update the text "cursor" position: */
        eps->vTopText = eps->vTopText - dv;
        /* Prepare for next line: */
        line = pend;
        /* Skip line break (ASCII CR, LF, or CR-LF), if any: */
        if ((*line) == '\r') { line++; }
        if ((*line) == '\n') { line++; }
        /* Are we done? */
        if ((*line) == '\000') { return; }
      }
  }

void epswr_dev_text_line
  ( epswr_figure_t *eps, 
    const char *text,
    const char *strut,
    double hAbs, double vAbs, 
    double rot, 
    bool_t clipped, 
    double hAlign, double vAlign,
    bool_t fill, bool_t draw
  )
  { FILE *wr = eps->wr;
    fprintf(wr, "%d %d ", draw, fill);
    epswr_dev_write_ps_string(wr, text);
    epswr_dev_write_ps_string(wr, strut);
    fprintf(wr, "  %5.3f %5.3f  %7.3f",  hAlign, vAlign, rot);
    fprintf(wr, "  %6.1f %6.1f", hAbs, vAbs);
    fprintf(wr, "  %s\n", (clipped ? "strfs" : "ncstrfs"));
    fflush(wr);
  }

/* MISCELLANEOUS */ 

void epswr_dev_set_verbose(epswr_figure_t *eps, const bool_t verbose)
  { eps->verbose = verbose; }

void epswr_dev_comment(epswr_figure_t *eps, const char *title)
  {
    FILE *wr = eps->wr;
    affirm(wr != NULL, "no wr");
    fprintf(eps->wr, "\n%% [%s]\n", title);
    if (eps->verbose) { fprintf(stderr, "[%s]\n", title); }
    fflush(wr);
  }

void epswr_dev_show_stack(epswr_figure_t *eps, int32_t code)
  { FILE *wr = eps->wr;
    affirm(wr != NULL, "no wr");
    fprintf(eps->wr, "\n%d stack pop\n", code);
    if (eps->verbose) { fprintf(stderr, "epswr_show_stack(%d)\n", code); }
    fflush(wr);
  }

void epswr_dev_flush (epswr_figure_t *eps)
  { FILE *wr = eps->wr;
    affirm(wr != NULL, "no wr");
    fflush(wr);
  }

void epswr_dev_compute_joining_arc_radius
  ( double ptx0, double pty0, 
    double ptx1, double pty1, 
    double ptx2, double pty2, 
    double rad,
    double *arad
  )
  {
    /* Compute vector {ux,uy} and distance {du} from point 1 to point 0: */
    double ux = ptx0 - ptx1, uy = pty0 - pty1;
    double du = hypot(ux,uy);
    /* Compute vector {vx,vy} and distance {dv} from point 1 to point 2: */
    double vx = ptx2 - ptx1, vy = pty2 - pty1;
    double dv = hypot(vx,vy);
    if ((du < 1.0e-8) || (dv < 1.0e-8))
      { /* Side with zero length, no rounding: */
        (*arad) = 0.0;
      }
    else
      { /* Sides with nonzero length. */
        /* Normalize {ux,uy} and {vx,vy}: */
        ux /= du; uy /= du;
        vx /= dv; vy /= dv;
        /* Compute unit vector {wx,wy} of bisector: */
        double wx = ux + vx;
        double wy = uy + vy;
        double dw = hypot(wx,wy);
        if (dw == 0)
          { /* Angle is 180 degrees, return arc with zero length: */
            (*arad) = 0.0;
          }
        else
          { /* Angle is not 180. */
            /* Normalize {wx,wy}: */
            wx /= dw; wy /= dw;
            /* Get cosine {ct} and sine {st} of half-angle: */
            double ct = wx*ux + wy*uy;
            double st = fabs(wx*uy - wy*ux);
            if (st < 1.0e-8)
              { /* Angle is essentially zero or 360. */
                (*arad) = 0.0;
              }
            else
              { /* Get distance {ds} from point 1 to ends of arc: */
                double dm = rad/st;
                double ds = dm * ct;
                double dsmax = fmin(du,dv);
                if (ds > 1.01*dsmax)
                  { /* Length removed by rounding is excessive: */
                    (*arad) = 0.0;
                  }
                else if (ds < 0.99*dsmax)
                  { /* Seems OK: */
                    (*arad) = rad;
                  }
                else 
                  { /* Reduce the radius slightly to ensure some bit of side remains: */
                    (*arad) = 0.99*rad*(dsmax/ds);
                  }
              }
          }
      }
  }
