/* Produceing EPS drawings of asexual descent trees. */
#ifndef drtree_plot_H
#define drtree_plot_H
/* Last edited on 2023-06-25 00:47:00 by stolfi */

#define drtree_plot_H_COPYRIGHT \
  "Duh?"

#define _GNU_SOURCE
#include <stdint.h>

#include <frgb.h>
#include <epswr.h>

#include <drtree.h>

/*
  The procedures in this interface generate a plot of an evolution history,
  described by an array of {drtree_node_t} records {dt[0..ni-1]}.
  The nodes must be in topological order, namely each non-null node must 
  occur after its parent, if any.
  
  GENERAL DIAGRAM LAYOUT

  The diagram is supposed to show the individual life spans and 
  parent-child relations in a given time range {tMin..tMax},
  which must contain the life span {q.tbr..q.tdt} of every non-null
  individual {q = dt[iq]} for {iq} in {0..ni-1}.

  The diagram's area is implicitly divided into an array of cells with
  {ncols = tMax-tMin-1} columns and {nrow} rows. Each column {j} in
  {0..ncols-1} corresponds to a time period {tMin + j}.

  Determining the layout of the diagram then means choosing a
  non-negative row {rdr(q)} to every non-null individual {q}.

  LIFE SPAN TRACES

  The graphical representation of the life span of a non-null individual
  {q} is assumed to occupy the cells in columns {q.tbr-tMin..q.tdt-tMin}
  of row {rdr(q)}. Those cells are the /trace/ or {q}.

  The row assignments must be such that the traces of any two visible
  individuals in the same row are separated by at least one clear
  (unoccupied) cell.

  BIRTH ARROWS

  If {q} is not null and not a root, a /parenting arrow/ will be drawn in the column
  {q.tbr-tMin} from the trace of {q}'s parent {p} to the cell
  representing the birth of {q}.  Recall that {p} must be non-null,
  and {q.tbr} must be in {p.tbr+1..p.tdt}, so the trace of {p} must include the
  column {q.tbr-tMin}.
  
  No parenting arrow will be drawn if {q} is null or is a root (has no parent).

  LAYOUT PLANNING GOAL
  
  Layout planning procedures generally try to minimize the number of
  rows {nrows} used, within the constraints above and possibly other
  constraints. In particular, every row will have the trace of at least one
  non-null individual. In the worst case, {nrows} could be equal to the number of
  non-null individuals; but it is often much less. */
 
void drtree_plot_individuals
  ( epswr_figure_t *eps, 
    int32_t tMin, 
    int32_t ncols, 
    int32_t nrows, 
    double Xstep, 
    double Ystep, 
    int32_t ni, 
    drtree_node_t dt[], 
    int32_t rdr[], 
    bool_t *fill,
    int32_t *chf
  );
  /* Plots onto {eps} the individuals described in {dt[0..ni-1]}.
  
    Each plot cell will have width {Xstep}, and each row has height {Ystep},
    both in millimeters.

    If {Xstep} and {Ystep} are positive, rows are numbered {0..nrows-1}
    from bottom up, and columns are numbered {0..ncols-1} from left to
    right. Negative values can be used to reverse the order.
    
    The birth of each individual is shown as a large dot. If {fill} is
    not {NULL}, it must be a vector {fill[0..ni-1]} specifies whether
    the dot is filled or hollow.  If {fill} is {NULL},
    all dots are filled.
    
    If {chf} is not {NULL}, it must be a vector {chf[0..ni-1]} that
    partitions the individuals into /families/. The value of {chf[iq]}
    must be either {-1} to mean 'no family' or the index of the family's
    /chief/, family member with smallest index. In particular, {q} is a
    chief if and only if {chf[iq]=iq}. The procedure will try to use a
    distinct color for each famaily, and some drab color for 'no family'
    individuals. If {chf} is {NULL}, it will asssume the same value {-1}
    for all individuals. */

epswr_figure_t *drtree_plot_create_eps_figure
  ( char *name,
    int32_t ncols, 
    int32_t nrows,
    double Xstep,
    double Ystep
  );
  /* Creates an EPS figure object that writes to the file "{name}.eps"
    
    The figure will assume a plot grid of {ncols} columns and {nrows} rows,
    with cells of width {Xstep} and height {Ystep} (millimeters), with a 
    small margin all around. */

void drtree_plot_time_line
  ( epswr_figure_t *eps, 
    int32_t tMin, 
    int32_t ncols, 
    int32_t nrows, 
    double Xstep, 
    double Ystep, 
    frgb_t *rgb,
    int32_t tRef
  );
  /* Plots a vertical line with color {rgb} across the drawing, at time {tRef}
    (that is, column {tRef-tMin}). If that column is outside {0..ncols-1}, does nothing.
    See {drtree_plot_individuals} for the meaning of the other parameters. */

#endif

