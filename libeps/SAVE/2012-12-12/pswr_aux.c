/* See pswr.h */
/* Last edited on 2012-12-12 03:49:41 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

#include <affirm.h>
#include <bool.h>
#include <jsfile.h>
#include <jstime.h>
#include <jsstring.h>

#include <pswr.h>
#include <pswr_def.h>
#include <pswr_aux.h>
#include <pswr_vis.h>

/* STREAM STATE */

bool_t pswr_within_canvas(PSStream *ps)
  { return ps->curSlot != INT_MAX; }

bool_t pswr_next_slot_available(PSStream *ps)
  { return (ps->curSlot != INT_MAX) && (ps->curSlot+1 < ps->hPicCount*ps->vPicCount); }

bool_t pswr_canvas_dirty(PSStream *ps)
  { return (ps->curSlot != INT_MAX) && (ps->curSlot >= 0); }

void pswr_begin_canvas(PSStream *ps, const char *pageName)
  { 
    affirm(! pswr_within_canvas(ps), "bad state");
    ps->curCanv++;
    affirm(ps->curCanv > 0, "curCanv bug");
    char *pageNameLocal = (char *)pageName;
    if (ps->file == NULL)
      { /* Provide a default for {pageName}, if NULL or empty: */
        if ((pageName == NULL) || (*pageName == '\000'))
          { pageNameLocal = NULL;
            asprintf(&pageNameLocal, "%06d", ps->curCanv);
          }
      }
   if (ps->eps)
      { if(ps->file == NULL)
          { /* Open new file: */
            char *fileName = NULL;
            asprintf(&fileName, "%s%s.eps", ps->prefix, pageNameLocal);
            ps->file = open_write(fileName, TRUE);
            free(fileName);
          }
        pswr_begin_encapsulated_figure(ps);
      }
    else
      { affirm(ps->file != NULL, "file bug");
        pswr_begin_document_page(ps, pageNameLocal);
      }
    /* Write canvas preamble commands to file: */
    double *fc = ps->fillColor;
    fc[0] = fc[1] = fc[2] = 0.500;
    pswr_write_canvas_graphics_setup_cmds(ps->file, fc);
    /* Mark writer state as {within_canvas}: */
    ps->curSlot = -1;
    /* Set the plot window to the full canvas: */
    double hSize = ps->hCanvasSize;
    double vSize = ps->vCanvasSize;
    /* Set window in writer (note: clobbers {ps->curSlot}, resets cell grid): */
    pswr_set_window(ps, 0.0, hSize, 0.0, vSize, 0.0, hSize, 0.0, vSize);
    /* Ensure that next picture slot is the first one: */
    ps->curSlot = -1;
    /* Free {pageNameLocal} if locally geerated: */
    if (pageNameLocal != pageName) { free(pageNameLocal); }
  }

void pswr_end_canvas(PSStream *ps)
  { 
    affirm(pswr_within_canvas(ps), "bad state");
    affirm(ps->curCanv >= 0, "curCanv bug");
    affirm(ps->file != NULL, "no file");
    if (ps->eps)
      { /* We have an open file with stuff in it: */
        pswr_end_encapsulated_figure(ps);
        fclose(ps->file); 
        ps->file = NULL;
      }
    else
      { pswr_end_document_page(ps); }
    /* Indicate {between_canvas} state: */
    ps->curSlot = INT_MAX;
  }

/* FILE INTIALIZATION AND FINALIZATION */

void pswr_begin_encapsulated_figure(PSStream *ps)
  { affirm(ps->eps, "bad file type");
    
    /* Initialize the font list with the default fonts: */
    pswr_initialize_font_list(&(ps->nFonts), &(ps->fonts));

    /* Initialize local state: */
    double *fc = ps->fillColor;
    fc[0] = fc[1] = fc[2] = 0.500;

    FILE *file = ps->file;
    char *date = today();
    pswr_write_eps_file_header(file, date, 0.0, ps->hCanvasSize, 0.0, ps->vCanvasSize);
    if(ps->verbose) { fprintf(stderr, "date = %s\n", date); } 
    pswr_write_canvas_header(file, "1", 1);
    fflush(file);
    free(date);
  }
  
void pswr_end_encapsulated_figure(PSStream *ps)
  { FILE *file = ps->file;
    affirm(ps->eps, "use of pswr_end_figure in non-EPS file");
    pswr_write_canvas_trailer(file, TRUE);
    pswr_write_eps_file_trailer(file, ps->nFonts, ps->fonts);
    fflush(file);
    pswr_free_font_list(&(ps->nFonts), &(ps->fonts));
  }

void pswr_begin_document(PSStream *ps)
  { affirm(! ps->eps, "bad file type");

    /* Initialize the font list with the default fonts: */
    pswr_initialize_font_list(&(ps->nFonts), &(ps->fonts));

    /* Initialize local state: */
    FILE *file = ps->file;
    char *date = today();
    if (ps->verbose) { fprintf(stderr, "date = %s\n", date); }
    pswr_write_ps_file_header(file, ps->paperSize, date);
    fflush(file);
    free(date);
  }

void pswr_begin_document_page(PSStream *ps, const char *pageName)
  { affirm(! ps->eps, "use of pswr_begin_document_page in EPS file");
    affirm(! pswr_within_canvas(ps), "page begin/end bug");
    
    FILE *file = ps->file;
    char *date = today();
    if (ps->verbose) { fprintf(stderr, "page %s date = %s\n", pageName, date); }
    pswr_write_canvas_header(file, pageName, ps->curCanv);
    pswr_write_date_and_page_show_cmds(file, date, pageName);
    free(date);
  }

void pswr_end_document_page(PSStream *ps)
  { FILE *file = ps->file;
    affirm(! ps->eps, "use of pswr_end_document_page in EPS file");
    affirm(pswr_within_canvas(ps), "page begin/end bug");
    pswr_write_canvas_trailer(file, FALSE);
  }

void pswr_end_document(PSStream *ps)
  { affirm(! ps->eps, "use of pswr_end_document in EPS file");
    affirm(! pswr_within_canvas(ps), "page begin/end bug");
    /* In case the client forgot to call pswr_end_page: */
    FILE *file = ps->file;
    pswr_write_ps_file_trailer(file, ps->curCanv, ps->nFonts, ps->fonts);
    fflush(file);
    pswr_free_font_list(&(ps->nFonts), &(ps->fonts));
  }

/* PAGE LAYOUT TOOLS */

void pswr_fix_slot_size
  ( double usableSize, 
    double picMargin,
    double capSize,
    int *count,
    double *picSize, 
    double *slotSize
  )
  { 
    if ((*picSize) <= 0)
      { if ((*count) <= 0)
          { (*slotSize) = usableSize;
            (*count) = 1;
          }
        else 
          { (*slotSize) = usableSize / (*count);
          }
        (*picSize) = (*slotSize) - 2*picMargin - capSize;
      }
    else
      { (*slotSize) = (*picSize) + 2*picMargin + capSize; }
      
    if ((*count) <= 0)
      { (*count) = (int)floor((usableSize + 0.1)/(*slotSize));
        if ((*count) <= 0) { (*count) = 1; }
      }
  }

/* INTERNAL WINDOW SETUP */

void pswr_save_window_data
  ( PSStream *ps,
    double xMin, double xMax,
    double yMin, double yMax,
    double hMin, double hMax,
    double vMin, double vMax
  )
  { /* Save the actual plot window: */
    ps->hMin = hMin; ps->hMax = hMax;
    ps->vMin = vMin; ps->vMax = vMax;
    /* Save the client scales: */
    ps->xMin = xMin; ps->xMax = xMax;
    ps->yMin = yMin; ps->yMax = yMax;
    /* Compute client-to-file scales: */
    ps->xScale = (hMax - hMin)/(xMax - xMin);
    ps->yScale = (vMax - vMin)/(yMax - yMin);
    /* Check for equal scales, allowing for some rounding errors: */
    { double scalesum = fabs(ps->xScale) + fabs(ps->yScale);
      double scalediff = fabs(ps->xScale - ps->yScale);
      affirm(scalediff/scalesum < 0.0001, "unequal scales");
    }
  }

void pswr_save_grid_data(PSStream *ps, int xn, int yn)
  { /* Save client grid size: */
    ps->xGridCount = xn; 
    ps->yGridCount = yn;
  }

/* DOCUMENT FONT LIST */

void pswr_initialize_font_list(int *nFontsP, char ***fontsP)
  {
    affirm((*fontsP) == NULL, "non-NULL font list");
    /* "Courier" is default, for page numbers etc. */
    int nFonts = 1;
    char **fonts = (char **)malloc(sizeof(char *));
    fonts[0] = txtcat("Courier", ""); 
    /* Update caller's variables: */
    (*nFontsP) = nFonts;
    (*fontsP) = fonts;
  }

void pswr_add_font(const char *font, int *nFontsP, char ***fontsP)
  {
    int nFonts = (*nFontsP);
    char **fonts = (*fontsP);
    /* Check whether font already is in the table: */
    int i;
    for (i = 0; i < nFonts; i++)
      { if (strcmp(fonts[i], font) == 0) { return; } }
    /* Ensure that there is space for a new entry: */
    if (((nFonts - 1) & nFonts) == 0)
      { /* Table is full ({nFonts} is a power of 2), reallocate it. */
        fonts = realloc(fonts, 2*nFonts*sizeof(char*)); 
        affirm(fonts != NULL, "out of mem"); 
      }
    /* Store a copy of the name: */
    fonts[nFonts] = txtcat(font, "");
    nFonts++;
    /* Update caller's variables: */
    (*nFontsP) = nFonts;
    (*fontsP) = fonts;
  }

void pswr_free_font_list(int *nFontsP, char ***fontsP)
  {
    int nFonts = (*nFontsP);
    char **fonts = (*fontsP);
    int i;
    for (i = 0; i < nFonts; i++) { free(fonts[i]); }
    free(fonts);
    (*fontsP) = NULL;
    (*nFontsP) = 0;
  }

/* INTERNAL FILE WRITING PROCEDURES */

void pswr_write_eps_file_header
  ( FILE *file, 
    const char *date,
    double hMin, double hMax,
    double vMin, double vMax
  )
  { fprintf(file, 
      "%%!PS-Adobe-3.0 EPSF-2.0\n"
    );
    fprintf(file, 
      "%%%%CreationDate: %s\n", date
    );
    fprintf(file, 
      "%%%%BoundingBox: %.0f %.0f %.0f %.0f\n",
      floor(hMin), floor(vMin), ceil(hMax), ceil(vMax)
    );
    fprintf(file, 
      "%%%%Pages: 1\n"
    );
    fprintf(file, 
      "%%%%DocumentFonts: (atend)\n"
    );
    fprintf(file, 
      "%%%%EndComments\n"
    );
    
    fprintf(file, 
      "/pswr$dict 6400 dict def \n"
      "pswr$dict begin\n"
    );
    pswr_write_proc_defs(file);
    fprintf(file, 
      "end\n"
    );

    fprintf(file, 
      "%%%%EndProlog\n"
    );
  }

void pswr_write_eps_file_trailer(FILE *file, int nFonts, char **fonts)
  { fprintf(file, "%%%%Trailer\n" );
    pswr_write_font_list(file, nFonts, fonts);
  }

void pswr_write_ps_file_header(FILE *file, const char *paperSize, const char *date)
  { fprintf(file, 
      "%%!PS-Adobe-3.0\n"
      "%%%%Pages: (atend)\n"
    );
    fprintf(file, 
      "%%%%DocumentFonts: (atend)\n"
    );
    fprintf(file, 
      "%%%%DocumentPaperSizes: %s\n", paperSize
    );
    fprintf(file, 
      "%%%%EndComments\n"
    );
    
    fprintf(file, 
      "/pswr$dict 6400 dict def \n"
      "pswr$dict begin\n"
    );
    pswr_write_proc_defs(file);
    fprintf(file, 
      "end\n"
    );

    fprintf(file, 
      "%%%%EndProlog\n"
    );
  }

void pswr_write_canvas_header(FILE *file, const char *pageName, int pageSeq)
  { if ((pageName != NULL) && (*pageName != '\000'))
      { fprintf(file, "%%%%Page: %d %d\n\n", pageSeq, pageSeq); }
    else
      { fprintf(file, "%%%%Page: %s %d\n\n", pageName, pageSeq); }

    fprintf(file, 
      "%% Save the standard coord system:\n"
      "/savedstate save def\n"
      "\n"
      "pswr$dict begin\n"
      "\n"
    );
    fflush(file);
  }

void pswr_write_canvas_trailer(FILE *file, int eps)
  { if (! eps) { fprintf(file, "showpage\n"); }
    fprintf(file, 
      "end %% pswr$dict\n"
      "\n"
      "savedstate restore\n"
      "%% Now we are back to the initial Postscript state system.\n"
      "\n"
    );
    fprintf(file, 
      "%%%%EndPage\n"
      "\n"
    );
    fflush(file);
  }
  
void pswr_write_ps_file_trailer(FILE *file, int nCanv, int nFonts, char **fonts)
  { fprintf(file, 
      "%%%%Trailer\n"
    );
    pswr_write_font_list(file, nFonts, fonts);
    fprintf(file, 
      "%%%%Pages: %d\n", nCanv
    );
  }
  
void pswr_write_font_list(FILE *file, int nFonts, char **fonts)
  { int i;
    fprintf(file, "%%%%DocumentFonts:");
    for (i = 0; i < nFonts; i++)
      { fprintf(file, " %s", fonts[i]); }
    fprintf(file, "\n");
  }

void pswr_write_date_and_page_show_cmds(FILE *file, const char *date, const char *page)
  { 
    if (date == NULL) { date = ""; }
    if (page == NULL) { page = ""; }
    char *fill = (strlen(page) == 0 ? "" : "  page ");
    char buf[strlen(date) + strlen(fill) + strlen(page) + 10]; /* 10 for safety... */
    sprintf(buf, "%s%s%s", date, fill, page); 
    fprintf(file, 
      "%% Print date and page:\n"
      "gsave\n"
      "  /Courier findfont 10 scalefont setfont\n"
      "  72 40 moveto\n"
    );
    fprintf(file, "  ");
    pswr_write_ps_string(file, buf, "\\267");
    fprintf(file, " show\n");
    fprintf(file, 
      "grestore\n"
      "\n"
    );
  }

void pswr_write_canvas_graphics_setup_cmds(FILE *file, double *fc)
  { fprintf(file, 
      "%% Round joints and caps:\n"
      "1 setlinecap 1 setlinejoin\n"
      "\n"
      "%% Black thin lines:\n"
      "0 setlinewidth 0 setgray [ ] 0 setdash\n"
      "\n"
    );
    pswr_write_fill_color_set_cmd(file, fc);

    /* Just in case the client forgets to set it up: */
    pswr_write_label_font_setup_cmds(file, "Courier", 8.0);    
  }
  
void pswr_write_window_setup_cmds
  ( FILE *file,
    double hMin, double hMax,
    double vMin, double vMax
  )
  { /* In spite of their names, the Postscript variables
      {xmin,xmax,ymin,ymax} are not client coordinates, 
      but actual plotting coordinates (in pt), namely 
      {hMin,hMax,vMin,vMax}. */
      
    fprintf(file, "/xmin %f def   %% min plottable x\n", hMin);
    fprintf(file, "/xmax %f def   %% max plottable x\n", hMax);

    fprintf(file, "/ymin %f def   %% min plottable y\n", vMin);
    fprintf(file, "/ymax %f def   %% max plottable y\n", vMax);

    fprintf(file, 
      "%% Set clipping path to boundary of plot area:\n"
      "initclip\n"
      "newpath\n"
      "  xmin ymin moveto\n"
      "  xmax ymin lineto\n"
      "  xmax ymax lineto\n"
      "  xmin ymax lineto\n"
      "  xmin ymin lineto\n"
      "clip\n"
      "\n"
    );

    /* Caption cursor. */
    
    fprintf(file, 
      "%% Left margin cursor for caption text:\n"
      "/dytext 10 pt mul def\n"
      "/xtext xmin def\n"
      "/ytext ymin def\n"
      "\n"
    );

    fprintf(file, 
      "%% Caption font setup:\n"
      "/captionfont\n"
      "  /Courier findfont\n"
      "  dytext scalefont\n"
      "def\n"
      "\n"
    );

    fflush(file);
  }

void pswr_write_grid_setup_cmds(FILE *file, int xn, int yn)
  { fprintf(file, "/xn %d def     %% grid cells along x axis\n", xn);
    fprintf(file, "/yn %d def     %% grid cells along y axis\n", yn);

    fprintf(file, "/xstep xmax xmin sub xn div def  %% x-size of grid cell\n");
    fprintf(file, "/ystep ymax ymin sub yn div def  %% y-size of grid cell\n");

    fflush(file);
  }

void pswr_write_label_font_setup_cmds(FILE *file, const char *font, double size)
  { fprintf(file, 
      "%% Label font setup:\n"
      "/labelfont\n"
      "  /%s findfont %.3f pt mul scalefont\n"
      "def\n" 
      "\n",
      font, size
    );
  }  

void pswr_write_color(FILE *file, double *clr)
  { fprintf(file, "  %5.3f %5.3f %5.3f", clr[0], clr[1], clr[2]); }

void pswr_write_fill_color_set_cmd(FILE *file, double *fc)
  { pswr_write_color(file, fc);
    fprintf(file, " sfc\n");
  }

void pswr_write_ps_string(FILE *file, const char *text, const char *newline)
  { const char *p;
    putc('(', file);
    for (p = text; *p != 0; p++)
      { if (*p == '\n')
          { fprintf(file, "%s", newline); }
        else if (*p == '(')
          { putc('\\', file); putc('(', file); }
        else if (*p == ')')
          { putc('\\', file); putc(')', file); }
        else if (*p == '\t')
          { putc(' ', file); putc(' ', file); }
        else if (*p == '\\')
          { putc('\\', file); putc('\\', file); }
        else if ((*p < ' ') || (*p > '~'))
          { fprintf(file, "\\%03o", *p); }
        else
          { putc(*p, file); }
      }
    fprintf(file, ")");
  }

void pswr_write_proc_defs(FILE *file)
  { 
    /* Global constants and variables: */
    
    fprintf(file, 
      "%% Units of measure:\n"
      "/pt 1.0 def\n"
      "/in pt 72.0 mul def \n"
      "/mm pt 72.0 25.4 div mul def\n"
      "\n"
    );
  
    fprintf(file, 
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

    fprintf(file, 
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

    fprintf(file, 
      "%% Draw an X-value grid line:\n"
      "%%   {x} xgrd --> \n"
      "/xgrd\n"
      "{ gsave\n"
      "  newpath\n"
      "    dup ymin moveto\n"
      "    ymax lineto\n"
      "    stroke\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fprintf(file, 
      "%% Draw an Y-value grid line:\n"
      "%%   {y} ygrd --> \n"
      "/ygrd\n"
      "{ gsave\n"
      "  newpath\n"
      "    dup xmin exch moveto\n"
      "    xmax exch lineto\n"
      "    stroke\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fprintf(file, 
      "%% Draw all grid lines:\n"
      "%%   gridlines --> \n"
      "/gridlines\n"
      "{ gsave\n"
      "    initclip\n"
      "    0 1 xn {\n"
      "      xstep mul xmin add xgrd\n"
      "    } for\n"
      "    0 1 yn {\n"
      "      ystep mul ymin add ygrd\n"
      "    } for\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fprintf(file,
      "%% Draw a frame around the current plot window:\n"
      "%%   frame --> \n"
      "/wframe\n"
      "{ gsave\n"
      "    initclip\n"
      "    newpath\n"
      "    xmin ymin moveto\n"
      "    xmax ymin lineto\n"
      "    xmax ymax lineto\n"
      "    xmin ymax lineto\n"
      "    xmin ymin lineto\n"
      "    closepath stroke\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );
    
    /* Fill color: */
    
    fprintf(file, 
      "%% Set fill color operator:\n"
      "%%   {R} {G} {B} sfc --> \n"
      "/sfc\n"
      "{ /fB exch def\n"
      "  /fG exch def\n"
      "  /fR exch def\n"
      "} bind def\n"
      "\n"
    );

    /* Figures: */
    
    fprintf(file, 
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
      "        { gsave fR fG fB setrgbcolor fill grestore stroke pop pop }\n"
      "        { 1 eq { fR fG fB setrgbcolor fill } if\n"
      "          1 eq { stroke } if\n"
      "        }\n"
      "      ifelse\n"
      "    }\n"
      "  ifelse\n"
      "} bind def\n"
      "\n"
    );
    
    fprintf(file, 
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
      "        { gsave fR fG fB setrgbcolor eofill grestore stroke pop pop }\n"
      "        { 1 eq { fR fG fB setrgbcolor eofill } if\n"
      "          1 eq { stroke } if\n"
      "        }\n"
      "      ifelse\n"
      "    }\n"
      "  ifelse\n"
      "} bind def\n"
      "\n"
    );
    
    fprintf(file, 
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
    
    fprintf(file, 
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

    fprintf(file, 
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
    
    fprintf(file, 
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
    
    fprintf(file, 
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
    
    fprintf(file, 
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

    fprintf(file, 
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

    fprintf(file, 
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

    fprintf(file, 
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
    
    fprintf(file, 
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

    fprintf(file, 
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

    fprintf(file, 
      "%% Grid cell operator:\n"
      "%%   {draw} {fill} {xi} {yi} cel --> \n"
      "/cel\n"
      "{ exch dup \n"
      "  %% --- draw fill yi, xi, xi\n"
      "  xstep mul xmin add exch 1 add xstep mul xmin add\n"
      "  %% --- draw fill yi, xlo, xhi\n"
      "  3 2 roll dup\n"
      "  %% --- draw fill xlo, xhi, yi, yi\n"
      "  ystep mul ymin add exch 1 add ystep mul ymin add\n"
      "  %% --- draw fill xlo, xhi, ylo, yhi\n"
      "  rec\n"
      "} bind def\n"
      "\n"
    );

    /* Labels and captions: */
    
    fprintf(file, 
      "%% Label printing operator:\n"
      "%%   {str} {xa} {ya} {xc} {yc} lbsh --> \n"
      "/lbsh\n"
      "{ labelfont setfont\n"
      "  newpath moveto\n"
      "    %% --- str, xa, ya\n"
      "  gsave 2 index false charpath flattenpath pathbbox grestore\n"
      "    %% --- str, xa, ya, lox, loy, hix, hiy\n"
      "  3 index 3 index currentpoint \n"
      "    %% --- str, xa, ya, lox, loy, hix, hiy, lox, loy, cx, cy\n"
      "  exch 4 1 roll exch sub\n"
      "  3 1 roll sub exch\n"
      "    %% --- str, xa, ya, lox, loy, hix, hiy, cx-lox, cy-loy\n"
      "  rmoveto\n"
      "    %% --- str, xa, ya, lox, loy, hix, hiy\n"
      "  exch 4 1 roll exch sub \n"
      "  3 1 roll sub exch\n"
      "    %% --- str, xa, ya, dx, dy\n"
      "  exch 4 1 roll mul -1 mul\n"
      "  3 1 roll mul -1 mul exch\n"
      "    %% --- str, -dx*xa, -dy*ya\n"
      "  rmoveto\n"
      "    %% --- str\n"
      "  show\n"
      "} bind def\n"
      "\n"
    );

    fprintf(file, 
      "%% Operator to move caption line cursor to new line:\n"
      "%%   nl --> \n"
      "/nl\n"
      "{ /ytext ytext dytext sub def\n"
      "} bind def\n"
      "\n"
    );

    fprintf(file, 
      "%% Operator to print aligned caption line at CP without clipping:\n"
      "%%   {s} {xa} capsh --> \n"
      "/tmpstr 50 string def\n"
      "/capsh\n"
      "{ gsave\n"
      "    initclip\n"
      "    captionfont setfont\n"
      "    gsave\n"
      "      newpath 0 0 moveto\n"
      "      1 index false charpath flattenpath pathbbox\n"
      "    grestore\n"
      "      %% --- str, xa, lox, loy, hix, hiy\n"
      "      %% --- str, xa, lox, loy, hix, hiy\n"
      "    4 2 roll pop pop pop\n"
      "      %% --- str, xa, hix\n"
      "    xmax xmin sub\n"
      "      %% --- str, xa, hix, xmax-xmin\n"
      "    exch sub mul 0\n"
      "      %% --- str, xa*((xmax-xmin) - hix)), 0\n"
      "    newpath xtext ytext moveto\n"
      "      %% --- str, xa*((xmax-xmin) - hix)), 0\n"
      "    rmoveto\n"
      "      %% --- str\n"
      "    show\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fprintf(file, 
      "%% Operator to print string at CP without clipping:\n"
      "%%   {s} shw --> \n"
      "/shw\n"
      "{ gsave\n"
      "    initclip\n"
      "    captionfont setfont show\n"
      "  grestore\n"
      "} bind def\n"
      "\n"
    );

    fflush(file);
  }
