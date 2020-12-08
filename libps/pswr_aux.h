/* pswr_aux.h -- internal prototypes for pswr.c. */
/* Last edited on 2007-12-26 12:16:58 by stolfi */

#ifndef pswr_aux_H
#define pswr_aux_H

#include <pswr.h>
#include <bool.h>

#include <stdio.h>

/* STREAM STATUS */

/* A stream that has been created but not been closed can be in two
  states: 
  
    * {between_canvas} ({curCanv >= 0},{curSlot = INT_MAX}), or 
    
    * {within_canvas} ({curCanv > 0}, {curSlot \in -1..hPicCount*vPicCount-1}). 
  
  For a standalone document, these states mean:

    * {between_canvas}: the file {psFile} is open and the file preamble
      has been written out to it. The page preamble has not been
      written yet.

    * {within_canvas}: the file is open and both the file and page
      preambles have been written.

  For an encapsulated stream, the states mean:

    * {between_canvas}: the file {psFile} is NULL.

    * {within_canvas}: the encapsulated picture file is open and all
      necessary preambles have been written to it.
      
  A newly-created stream is in the {between_canvas} state. The
  transitions between these states are effected by the procedures
  {pswr_begin_page} and {pswr_end_page}.  Drawing commands may be
  issued only in the {within_canvas} state. */

bool_t pswr_within_canvas(PSStream *ps);
  /* TRUE between a {pswr_begin_canvas} and the next {pswr_end_canvas}. */

bool_t pswr_next_slot_available(PSStream *ps);
  /* TRUE if there is some unused picture slot in the current canvas,
    besides the current one. */

bool_t pswr_canvas_dirty(PSStream *ps);
  /* TRUE if {pswr_new_picture} or {pswr_set_window} were issued
    to the current canvas. */
    
/* CANVAS OPEN AND CLOSE */

void pswr_begin_canvas(PSStream *ps, const char *pageName);
  /* Starts a new canvas of {ps}, and sets the current
    picture slot to -1. Must be followed eventually by a call to
    {pswr_end_canvas}. */

void pswr_end_canvas(PSStream *ps);
  /* Terminates the current canvas. Must be followed by
    {pswr_begin_canvas} or {pswr_close_stream}. */

/* CANVAS LAYOUT */

void pswr_fix_slot_size
  ( double usableSize, 
    double picMargin,
    double capSize,
    int *count,
    double *picSize, 
    double *slotSize
  );
  /* Adjusts {*count}, {*picSize}, and {*slotSize}, if needed,
    given the parameters {usableSize}, {picMargin}, and {capSize}.
    
    The parameters are assumed to refer to a string of picture slots
    arranged in either horizontal or vertical direction. We will
    assume vertical, since horizontal is analogous except that
    {capSize} is probably zero. All parameters are taken to be in pt.
    
    The {usableSize} argument is the total heigh of the canvas minus the
    page margins. The {picMargin} is the amount of blank space that
    must be allowed on either side of each picture slot, except before
    the first slot and after the last one. The {capSize} is the extra
    space, for caption text, that must be allowed below each slot. The
    number of pictures in the string is {*count}. The size of a
    picture area, excluding all margins, is {*picSize}; and the
    distance between successive slots is {*slotSize}.
  
    If {count <= 0} and {picSize <= 0}, then count is set to 1 
    and {picSize} is set to {usableWidth}. */

/* ENCAPSULATED FIGURES */

void pswr_begin_encapsulated_figure(PSStreamPrivate *ps);
  /* Initializes the underlying stream as an Encapsulated PostScript file.
    
    Writes the file's preamble, which defines some auxiliary
    Postscript operators and constants, sets up the caption font,
    resets some graphics variables (line width and color, fill color,
    etc.), etc. The figure's bounding box (for the "%%BoundingBox"
    comment) is taken from the fields {ps->hPageSize,ps->vPageSize}.
    The plotting window is initialized to cover the whole figure, with
    client coordinates equal to the device coordinates. */
    
void pswr_end_encapsulated_figure(PSStreamPrivate *ps);
  /* Finalizes an Encapsulated PostScript figure. The file
    is flushed but is left open. */

/* NON-ENCAPSULATED MULTIPAGE DOCUMENTS */

void pswr_begin_document(PSStreamPrivate *ps);
  /* Initializes a file to contain a non-encapsulated, possibly
    multi-page Postscript document.
    
    The procedure writes the file preamble, which defines some
    auxiliary Postscript operators and constants, defines the 
    caption font, etc.
    
    The client must initialize each page with
    {pswr_begin_document_page}, before issuing any drawing commands,
    and finalize it with {pswr_end_document_page}. */

void pswr_begin_document_page(PSStreamPrivate *ps, const char *page);
  /* Starts a new page in a non-encapsulated document.  
  
    The {page} string is the logical page number; if not NULL and
    non-empty, it will be printed at the bottom of the page. It is is
    also used as the first argument of the "%%Page" structuring
    directive (where it defaults to the sequential page number,
    counting from 1). 
    
    At each new page, the plotting window is reset to the whole page,
    and the whole Postscript state (including line width and color,
    fill color, etc.) is restored to some default state. */

void pswr_end_document_page(PSStreamPrivate *ps);
  /* Finalizes a page in a non-encapsulated document. */

void pswr_end_document(PSStreamPrivate *ps);
  /* Finalizes a multipage document.  The file is flushed but left open. */
    
/* INTERNAL SCALING/CLIPPING PROCEDURES */

void pswr_save_window_data
  ( PSStream *ps,
    double xMin, double xMax,
    double yMin, double yMax,
    double hMin, double hMax,
    double vMin, double vMax
  );
  /* Saves the given data as the current plot rectangle
    ({ps->hMin,ps->hMax,ps->vMin,ps->vMax}) and the current client coordinate
    ranges ({ps->xMin,ps->xMax,ps->yMin,ps->yMax}). */
    
void pswr_save_grid_data(PSStream *ps, int xn, int yn);
  /* Saves the given data as the current cell grid counts
    ({ps->xGridCount,ps->yGridCount}). */
    
/* FILE PREMBLES AND POSTAMBLES */

/* The file preamble (produced by both {pswr_write_eps_file_header} and
  {pswr_write_ps_file_header}) creates a Postscript dictionary
  "ps$dict" that contains definitions for special Postscript
  operators and some state variables, such as the plotting rectangle
  and grid size. This dictionary is closed before the end of the
  preamble, and must be re-opened at each canvas. */

void pswr_write_eps_file_header
  ( FILE *psFile, 
    const char *date,
    double hSize,
    double vSize
  );
  /* Writes the file preamble for an encapsulated Postscript file. The
    figure's "%%BoundingBox" will be set to be from {0,0} to 
    {hSize,vSize}, rounded outwards to integers. The 
    The sizes are assmed to be in pt
    
    The preamble includes the creation of a Postscript dictionary where
    a number of operators and variables are defined, including
    Postscript variables {hSize} and {vSize}, which are assumed to
    be in pt. This procedure must be followed by {pswr_write_canvas_header}
    to complete the setup for plotting commands. */

void pswr_write_eps_file_trailer(FILE *psFile, int nfonts, char **fonts);
  /* Writes the file postamble for an encapsulated Postscript file,
    namely the "%%Trailer" comment followed by the "%%DocumentFonts:
    {F}" (where {F} is the given list of font names). This call must
    be preceded by {pswr_write_canvas_trailer}. */

void pswr_write_ps_file_header
  ( FILE *file, 
    const char *paperSize, 
    bool_t landscape, 
    double hSize,
    double vSize,
    const char *date
  );
  /* Writes the file preamble for a
    non-encapsulated Postscript document file, ending with "%%EndProlog".
    The {papersize} string is written as the structuring comment
    "%%DocumentPageSizes:". 
    
    The preamble includes the creation of a Postscript dictionary where
    a number of operators and variables are defined, including
    Postscript variables {hSize} and {vSize}, with the given values,
    which are assumed to be in pt. If the {landscape} argument is true,
    there will be commands in the preamble that rotate and translate the
    coordinate system accordingly. In this case, hSize and {vSize} must
    be page dimensions in that coordinate system. This procedure must be
    followed by {pswr_write_canvas_header} to complete the setup for
    plotting commands. */
    
void pswr_write_ps_file_trailer(FILE *psFile, int npages, int nfonts, char **fonts);
  /* Writes the trailer for a non-encapsulated Postscript document.
    The trailer consists of a "%%Trailer" line, plus the directives
    "%%Pages: {N}" (where {N} is the given page count) and
    "%%DocumentFonts: {F}" (where {F} is the given list of font
    names). This call must be preceded by {pswr_write_canvas_trailer}. */
    
/* DOCUMENT FONT LIST */

void pswr_initialize_font_list(int *nfontsP, char ***fontsP);
  /* Initializes the font list with the "Courier" font. */

void pswr_add_font(const char *font, int *nfontsP, char ***fontsP);
  /* Appends a copy of the given {font} string to the font list. */

void pswr_free_font_list(int *nfontsP, char ***fontsP);
  /* Frees all storage used by the font list and sets it to NULL. */

/* CANVAS PREAMBLES AND POSTAMBLES */

void pswr_write_canvas_header(FILE *psFile, const char *pageName, int pageSeq);
  /* Writes the preamble for a new canvas of an EPS or non-EPS file.
    May be called only once per EPS file, and once per page of a
    non-EPS file. The preamble includes the "%%Page" directive, if
    applicable, and commands that save the current Postscript state
    (with "save"), and open the "ps$dict" dictionary. */
    
void pswr_write_canvas_trailer(FILE *psFile, int eps);
  /* Writes the postamble for a canvas in an EPS or
    non-EPS file.  Must be matched to a preceding {pswr_write_canvas_header}
    
    The postamble includes a "showpage" command, if {eps} is false. In
    any case, the postamble includes commands that close the "ps$dict"
    dictionary, and restore the Postcript state that was saved when
    the canvas was started. */

/* MISCELLANEOUS WRITING PROCEDURES */

void pswr_write_date_and_page_show_cmds(FILE *file, const char *date, const char *page);
  /* Writes to {ps->psFile} some commands that print the 
    given date and canvas number near the bottom of the
    current canvas.  */

void pswr_write_canvas_graphics_setup_cmds(FILE *psFile, double *fc);
  /* Writes to {ps->psFile} some commands that setup some graphics
    variables, suitable for the beginning of a new canvas. In
    particular the stroke color is reset to black, and the fill color
    to {fc}. */
  
void pswr_write_window_setup_cmds
  ( FILE *psFile,
    double hMin, double hMax,
    double vMin, double vMax
  );
  /* Writes to {psFile} some commands that setup the Postscript
    variables {xMin,xMax,yMin,yMax} to the values
    {hMin,hMax,vMin,vMax}, and set the clip path to that rectangle.
    Those commands also set the Postscript variables {xtext,ytext},
    which define the start of the current caption line; and {dytext}
    that defines the caption interline spacing (in pt). Finally, they
    set the {captionfont} varriable to the caption font, properly
    scaled. */

void pswr_write_grid_setup_cmds(FILE *psFile, int xn, int yn);
  /* Writes to {psFile} some commands that setup the Postscript
    variables {xn,yn} that define the cell grid */

void pswr_write_label_font_setup_cmds(FILE *psFile, const char *font, double size);
  /* Writes to {psFile} commands that create an instance of the font
    named {font}, scaled by {size} pt, and save it in the Postscript
    variable {labelfont} (which is used by {pswr_label}). */

void pswr_write_ps_string(FILE *psFile, const char *text, const char *newline);
  /* Writes the given text out  to /psFile/, in Postscript form, properly
    quoting any special chars and parentheses. 
    
    Replaces any embedded '\n' by the given /newline/ string --- which
    could be, for example, ") show mynl (". */

void pswr_write_color(FILE *psFile, double *clr);
  /* Writes the color value {clr[0..2]} to the Postscript file,
    with 3 decimals, preceded by spaces.  */

void pswr_write_font_list(FILE *psFile, int nfonts, char **fonts);
  /* Writes the structured comment "%%DocumentFonts: {F}" 
    where {F} is the given list of font names. */

void pswr_write_fill_color_set_cmd(FILE *psFile, double *fc);
  /* Writes the Postoscript command that saves {fc[0..2]}
    as the RGB fill color for subsequent figure commands. */

#endif
