#ifndef haf_draw_H
#define haf_draw_H
/* Generating EPS drawing of half-edge meshes and data structures. */
/* Last edited on 2023-10-05 12:21:35 by stolfi */

#define half_draw_H_copyright \
  "Copyright (C) 2023 Jorge Stolfi, UNICAMP.\n\n" jslibs_copyright

#include <jslibs_copyright.h>
#include <bool.h>
#include <epswr.h>
#include <r2.h>

#include <haf.h>

/* DRAWING THE MESH ELEMENTS */
    
typedef struct haf_draw_data_t
  { int64_t ne;      /* Number of undirected edges. */
    haf_arc_t edge; /* Natural arc on each edge, indexed {0..ne-1}. */
    uint64_t eid0;   /* Lowest edge ID. */
    int64_t nv;      /* Number of vertices. */
    uint64_t *iorg;  /* IDs of origin vertex of each arc. indexed {0..2*ne-1}. */
    r3_t *vpos;      /* Coordinates of each vertex, indexed {0..nv-1}. */
    bool_t *vshow;   /* Whether each vertex is to be drawn, indexed {0..nv-1}. */
    int64_t nf;      /* Number of faces. */
    uint64_t *ileft; /* ID of left face of each arc, indexed {0..2*ne-1}. */
    r3_t *fctr;      /* Nominal center of each face, indexed {0..nf-1}. */
    bool_t *fshow;   /* Whether each face center or record is to be shown, indexed {0..nf-1}. */
  } haf_draw_data_t;
  /* Data needed to draw the mesh or the records and pointers of a half-edge structure. */

void haf_draw_edge(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a);
  /* Draws the edge represented by {a} and its opposite {b}. 
  
    ??
    If {laba} and {labb} are {NULL}, and {nodes} is false, draws the
    unoriented edge from {a.org} to {a.dst}, displaced laterally by
    {(a.bend - b.bend)/2}. Uses the current pen settings.
    
    If laba} and {labb} are {NULL}, and {nodes} is true, draws instead
    the boxes representing the arc records, with the {.sym}, {.org}, and
    {.left} pointers. The {sym} pointers are displaced laterally by
    {a.bend} and {b.bend}
    
    If {laba} or {labb} are not both {NULL}, the drawing of the edge and
    boxes is suppressed. Instead, if {laba} is not null, and {nodes} is
    false, draws an arrohead on the edge at about 1/3 of the way from
    {a.org} to {a.dst}, and writes the string {laba} nearby, to the left
    of the arc {a}, with the current label font and fill color. The
    label is placed near the arrowhead if {nodes} is false, or near the
    record box if {nodes} is true.
    
    Ditto for {labb}, if not {NULL}, but with respect to the opposite arc {b}.*/

void haf_draw_arc_arrow(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a);
  /* Draws an arrowhead on a Bézier curve, about 1/3 of the way. */

void haf_draw_vert(epswr_figure_t *eps, haf_draw_data_t *dd, uint64_t ixv);
  /* Draws the vertex {v} as a dot at {v.pos}.  The drawing is suppressed if {v.show} is false. */
      
void haf_draw_face(epswr_figure_t *eps, haf_draw_data_t *dd, uint64_t ixf);
  /* Draws the face {f}, as a dot at {f.ctr}.   The drawing is suppressed if {v.show} is false. */

/* DRAWING THE DATA STRUCTURE */

void haf_draw_arc_record(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a);
  /* Draws the box representing the {haf_arc_t} record {a} (but not {a.sym}).
    The drawing is suppressed if {a.org.show} is false. */

void haf_draw_vert_record(epswr_figure_t *eps, haf_draw_data_t *dd, uint64_t ixv);
  /* Draws a circle representing the vertex record {v} at the position {v.pos}.
    The drawing is suppressed if {v.show} is false. */

void haf_draw_face_record(epswr_figure_t *eps, haf_draw_data_t *dd, uint64_t ixf);
  /* Draws a circle representing the specified record {a} (but not {a.sym}).
    The drawing is suppressed if {a.org.show} is false. */

void haf_draw_arc_sym_link(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a);
  /* Draws the link {a.sym}. */
  
void haf_draw_arc_next_link(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a);
  /* Draws the link {a.next}. */

void haf_draw_arc_org_link(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a);
  /* Draws the link {a.org}. */
  
void haf_draw_arc_left_link(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a);
  /* Draws the link {a.left}. */

void haf_draw_vert_out_link(epswr_figure_t *eps, haf_draw_data_t *dd, uint64_t ixv);
  /* Draws the link {v.out}. */

void haf_draw_face_side_link(epswr_figure_t *eps, haf_draw_data_t *dd, uint64_t ixf);
  /* Draws the link {f.side}. */

/* LABELING THINGS */

void haf_draw_arc_label(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a, char *lab, double du, double dv);
  /* Draws the label {lab} of an arc {a},
    near the position of the box or arrowhead.  The 
    ref point of the text will be displaced by {(du,dv)} from that point.  The
    coordinates of {dsp} are interpreted in the system whose axes are the
    direction tangent to the edge and the left-pointing normal. */

void haf_draw_vert_label(epswr_figure_t *eps, haf_draw_data_t *dd, uint64_t ixv, char *lab, double dx, double dy);
  /* Draws the label {lab} of a vertex {v}, with the text's ref point 
    displaced by {(dx,dy)} from its position. */

void haf_draw_face_label(epswr_figure_t *eps, haf_draw_data_t *dd, uint64_t ixf, char *lab, double dx, double dy);
  /* Draws the label {lab} of a face {f}, with teh reference point displaced by {(dx,dy)} from 
    its nominal center. */

/* GENERIC DRAWING COMMANDS */

void haf_draw_label(epswr_figure_t *eps, haf_draw_data_t *dd, r2_t *p, char *lab, double hAlign, double vAlign);
  /* Draws the label {lab} with the relative point {hAlign,vAlign} of the text at the position {p}. */

void haf_draw_link(epswr_figure_t *eps, haf_draw_data_t *dd, r2_t *p, r2_t *q, double bend);
  /* Draws a link arrow from {p} to {q} that bends {bend} mm to the left of the
    line from {p} to {q}. */

#endif
