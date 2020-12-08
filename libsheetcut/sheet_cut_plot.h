#ifndef sheet_cut_plot_H
#define sheet_cut_plot_H
 
/* Plotting layouts of rectangular plates on sheet stock. */
/* Last edited on 2020-11-05 17:18:50 by jstolfi */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>

#include <r2.h>
#include <bool.h>
#include <epswr.h>
#include <jsmath.h>

#include <sheet_cut.h>

#define sheet_cut_plot_SCALE (72/25.4)
  /* Plot scale (plot points per real mm). */
  
#define sheet_cut_plot_FIG_MARGIN_PT (8.0)
  /* EPS figure margin and text line spacing in points. */

#define sheet_cut_plot_FONT_NAME "Courier"
#define sheet_cut_plot_FONT_SIZE_PT (18.0)
  /* Font name and nominal font height in points. */
    
epswr_figure_t *sheet_cut_plot_new_figure
  ( char *outPrefix, 
    char *sheet_tag,
    char *sheet_mat,
    double sheet_thk,
    r2_t sheet_dim,
    double scrap_mrg,
    double cur_px,
    double cur_py,
    double cur_ht
  );
  /* Starts a new EPS figure for another physical sheet of material in
    which pieces have been placed.
    
    The {sheet_tag} is a string that identifies the sheet among many sheets. Assumes that
    the sheet material, thickness, and physical size are {sheet_thk}, {sheet_mat},
    and {sheet_dim}. Assumes that the useful area is the physical sheet
    minus a border of width {scrap_mrg} all around.
    
    Will set the plot window to the usable area of the physical sheet,
    with the origin at its lower left corner. 
    
    If {cur_px},{cur_py}, and {cur_ht} are not {NAN}, 
    also draws in blue a step-like polygonal that starts 
    horizontally at the left margin with ordinate {cur_ht},
    drops to ordinate {cur_py} at abscissa {cur_px},
    and continues horizontally to the right margin. */

void sheet_cut_plot_end_figure(epswr_figure_t *eps);
  /* Ends the given EPS figure.  */
 
void sheet_cut_plot_all_nodes(epswr_figure_t *eps, sheet_cut_node_t *pc, r2_t org_pc);
  /* Draws the given node {pc} and all its descendant sub-nodes
    on the figure {eps}, assumed to have index {sheet_ix}.  
    Assumes that {sheet_cut_plot_new_figure} has been called on {eps}.  
  
    If {pc} is a plate, just draws it with {sheet_cut_plot_plate}.
    
    If {pc} is a block, draws its bounding box in red with
    {sheet_cut_plot_block_bbox}. Then recurively plots
    {pc}'s descendant nodes over that.
    
    In either case, the lower left corner of {pc} (if plate) or of its
    bounding box (if block) will be drawn at {org_pc + pc.pos}.
    
    As a side effect, all sub-nodes will be actually flipped when
    appropriate, and all fields {.sub_flip} will then be set to {FALSE},
    before the plotting itself. Assumes that {pc} itself does not need
    flipping.
    
    This procedure is usually called with {pc} being the root block of a
    whole sheet, and {org_pc} being the coordinates of the low corner of
    the usable area of the sheet. */
 
void sheet_cut_plot_block_bbox(epswr_figure_t *eps, r2_t pos, r2_t size);
  /* Draws the outline of a block's bounding box in red,
    assuming that its low corner is {pos} and its size is {size}. */

void sheet_cut_plot_plate(epswr_figure_t *eps, r2_t pos, r2_t size, char *tag);
  /* Draws a plate on the current figure, as a rectangle
    with lower left corner at {pos} and the given {size}.
    Assumes that {sheet_cut_plot_new_figure} has been called. 
    
    Also writes inside the plate, near the low and high corners,
    the {tag}, {size}, and the corner coordinates. Assumes that 
    the plate is big enough for that. */
   
void sheet_cut_plot_labels(epswr_figure_t *eps, r2_t pos, r2_t size, char *tag);
  /* Writes on the EPS figure the corner coordinates, size, and tag of a plate
    whose lower left corner is {pos} and whose size is {size}.
    
    Corner coordinates will be written near the corresponding corners.
    The tag and size will be written at the center of the plate.
    
    Assumes that the {size} of the plate is sufficient to contain 
    the labels, and nothing else will overlap with it. */

#endif
