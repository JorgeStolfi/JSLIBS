/* pnmift_path_cost_fn.h - some path cost functions for PNM segmentation */
/* Last edited on 2016-04-01 01:29:10 by stolfilocal */

#ifndef pnmift_path_cost_fn_H
#define pnmift_path_cost_fn_H

#include <ift.h>

#include <pnmift_arc_cost_fn.h>

typedef ift_path_cost_t pnmift_path_cost_t;

typedef pnmift_path_cost_t pnmift_path_cost_fn_t(pnmift_path_cost_t sC, pnmift_arc_cost_t aC);
  /* Type of a function that returns the cost of a path that ends with an
    arc {(s,t)}, given the cost {sC} of the path up to {s} and the 
    cost {aC} of the arc {(s,t)}. */

pnmift_path_cost_fn_t *pnmift_path_cost_fn_from_name(char *name);
  /* Returns the function {pnmift_path_cost_fn_{name}} from the list below. */
  
#define pnmift_path_cost_fn_INFO \
  "        max\n" \
  "        " pnmift_path_cost_fn_max_INFO "\n" \
  "\n" \
  "        sum\n" \
  "        " pnmift_path_cost_fn_sum_INFO "\n" \
  "\n" \
  "        mono\n" \
  "        " pnmift_path_cost_fn_mono_INFO ""

pnmift_path_cost_t pnmift_path_cost_fn_max(pnmift_path_cost_t sC, pnmift_arc_cost_t aC);

#define pnmift_path_cost_fn_max_INFO \
  "  Maximum of the initial vertex cost and of all arc" \
  " costs in the path.  To obtain the watershed transform of a" \
  " monochrome image: (1) set the the trivial path cost of pixel {p} to the image" \
  " value at {p} if {p} is a local minimum, to {+oo} otherwise;" \
  " (2) define the arc cost {w(p,q)} as  the value of the image at" \
  " {q}; (3) and run the IFT with this path cost and 4- or 8-neighborhood."

pnmift_path_cost_t pnmift_path_cost_fn_sum(pnmift_path_cost_t sC, pnmift_arc_cost_t aC);

#define pnmift_path_cost_fn_sum_INFO \
  "  Sum of the initial vertex cost and of all arc costs" \
  " along path. To obtain the approximate geodesic Voronoi of a" \
  " rectangular cell grid: (1) set the trivial path cost of {p} to zero if {p} is a" \
  " seed cell, to {+oo} otherwise; (2) define the arc cost {w(p,q)} as the" \
  " effective length of the segment {pq}; (3) run the IFT with this path cost" \
  " and a large neighborhood."
  
pnmift_path_cost_t pnmift_path_cost_fn_mono(pnmift_path_cost_t sC, pnmift_arc_cost_t aC);

#define pnmift_path_cost_fn_mono_INFO \
  "  If the arc costs are non-decreasing along the path, the cost of the" \
  " last arc; else {+oo}.  To find all local minima of a monochrome image: (1) set" \
  " the trivial path cost of every pixel {p} to the value of the image at {p}; (2) define" \
  " the arc cost {w(p,q)} as the image value at {q}; (3) run the IFT with this" \
  " path cost and FIFO tie-breaking.  To find one representative out of every" \
  " local minimum basin, use LIFO tie-breaking instead."

#endif
