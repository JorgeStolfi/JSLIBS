/* See {haf_draw.h}. */
/* Last edited on 2023-10-05 19:36:55 by stolfi */

#define haf_draw_C_copyright \
  "Copyright © 2023 State University of Campinas (UNICAMP).\n\n" jslibs_copyright

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>

#include <jslibs_copyright.h>
#include <bool.h>
#include <vec.h>
#include <affirm.h>
#include <epswr.h>
#include <r2.h>
#include <r2_bezier.h>
#include <r3.h>

#include <haf.h>
#include <haf_draw.h>
  
#define haf_draw_arc_record_TIME 0.35
  /* Relative position of arc box along arc. */

/* Detail dimensions (all in mm): */

#define haf_draw_vert_RAD 1.0
  /* Radius of vertex dots in mesh view. */
  
#define haf_draw_face_RAD 0.0
  /* Radius of face dot in mesh view. */

#define draw_edge_arc_arrow_WIDTH 3.0
#define draw_edge_arc_arrow_LENGTH 4.0
  /* Dimensions of arc arrowhead on mesh view. */

#define haf_draw_arc_record_HU 3.0
#define haf_draw_arc_record_HV 2.5
  /* Half-sides of an arc record box in data structure view. */

#define haf_draw_vert_record_RAD 3.0
  /* Radius of vertex record in data structure view. */

#define haf_draw_face_record_RAD 3.0
  /* Radius of face record in data structure view. */

#define haf_draw_link_arrow_WIDTH 1.5
#define haf_draw_link_arrow_LENGTH 2.0
  /* Dimensions of link arrowhead on data structure view. */

#define haf_draw_link_dot_RAD 0.35
  /* Radius of link tail dot on data structure view. */

void haf_draw_choose_lab_align(r2_t *dp, double *hAlign_P, double *vAlign_P);
  /* Chooses the parameters {hAlign} and {vAlign} for {epswr_label}
    based on the direction of the displacement {dp} from the labebeld
    point. The {displacement {dp} must be in the world coordinate
    system, not the local tangent system. */

void haf_draw_get_edge_bezier_points(haf_arc_t a, r2_t *p0, r2_t *p1, r2_t *p2, r2_t *p3);
  /* Returns the Bézier control points {p0,p1,p2,p3} for the edge of the arc {a}. 
    The endpoints {p0} and {p1} will be {a.org.pos} and {a.dst.pos}.  The other 
    two points are computed based on the endpoints and the bending parameter {a.bend}. */

void haf_draw_get_arc_point(haf_arc_t a, double t, r2_t *p, r2_t *pu, r2_t *pv);
  /* Returns the point {*p} that is on the edge of arc {a}, {t} of the way from {a.org}
  to {a.dst}. If {pv} is not {NULL}, stores the unit tangent direction vector in {*pu},
  and, {pv} is not {NULL}, stores into it the unit normal direction pointing to the left.  */

r2_t haf_draw_get_arc_record_point(r2_t *p, r2_t *pu, r2_t *pv, double su, double sv);
  /* Returns the absolute {(x,y)} coordinates of the point of an arc record box that has 
    coordinates {(su*HU,sv*HV)} in the box's coordinate system, which has origin {p}, 
    axis direction vectors {pu} and {pv}; where {HU} is {haf_draw_arc_record_HU}, and
    similarly for {HV}. */

r2_t haf_draw_get_vert_record_point(r2_t *p, r2_t *q);
  /* Returns a point where the line from {p} to {q} crosses
    the edge of the box that represents the vertex record with center {q}. */

r2_t haf_draw_get_face_record_point(r2_t *p, r2_t *q);
  /* Returns a point where the line from {p} to {q} crosses
    the edge of the box that represents the face record with center {q}. */

/* DRAW MESH VIEW */

void haf_draw_edge(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a)
  { 
    haf_arc_t b = a->sym;
    haf_vert_t *org = a->org;
    haf_vert_t *dst = b->org;
    
    if ((! org->show) && (! dst->show)) { return; } 
    r2_t p0, p1, p2, p3;
    haf_draw_get_edge_bezier_points(a, &p0, &p1, &p2, &p3);
    if ((! org->show) || (! dst->show))
      { r2_t p01, p12, p23,  p012,  p123, p0123;
        r2_bezier_split
          ( 0, 1, &p0, &p1, &p2, &p3, 
            0.5, &p01, &p12, &p23, &p012, &p123, &p0123, NULL 
          );
        if (! dst->show)
          { p1 = p01; p2 = p012; p3 = p0123;}
        else if (! org->show)
          { p0 = p0123; p1 = p123; p2 = p23; }
        else
          { assert(FALSE); }
      }
      epswr_curve(eps, dd, p0.c[0], p0.c[1], p1.c[0], p1.c[1], p2.c[0], p2.c[1], p3.c[0], p3.c[1]);
  }
    
void haf_draw_arc_arrow(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a)
  { r2_t p, pu;
    double tm = haf_draw_arc_record_TIME;
    haf_draw_get_arc_point(a, tm, &p, &pu, NULL);
    double width = draw_edge_arc_arrow_WIDTH;
    double length = draw_edge_arc_arrow_LENGTH;
    r2_t q; r2_mix(1.0, &p, 2.0, &pu, &q);
    epswr_arrowhead(eps, dd, p.c[0], p.c[1], q.c[0], q.c[1], width, length, 1.0, TRUE, TRUE);
  }

void haf_draw_vert(epswr_figure_t *eps, haf_draw_data_t *dd, uint64_t kv)
  {
    if (v->show)
      { double rad = haf_draw_vert_RAD;
        epswr_dot(eps, dd, v->pos.c[0], v->pos.c[1], rad, TRUE, TRUE);
      }
  }

void haf_draw_face(epswr_figure_t *eps, haf_draw_data_t *dd, haf_face_t *f)
  {
    if (f->show)
      { double rad = haf_draw_face_RAD;
        if (rad > 0) { epswr_dot(eps, dd, f->ctr.c[0], f->ctr.c[1], rad, TRUE, TRUE); }
      }
  }

void haf_draw_get_edge_bezier_points(haf_arc_t a, r2_t *p0, r2_t *p1, r2_t *p2, r2_t *p3)
  {
    haf_arc_t b = a->sym;
    haf_vert_t *org = a->org;
    haf_vert_t *dst = b->org;

    (*p0) = (r2_t){{ org->pos.c[0], org->pos.c[1] }}; 
    (*p3) = (r2_t){{ dst->pos.c[0], dst->pos.c[1] }}; 
    double ebend = (a->bend - b->bend)/2;
    r2_bezier_from_bend(p0, p3, ebend, p1, p2);
  }

void haf_draw_get_arc_point(haf_arc_t a, double t, r2_t *p, r2_t *pu, r2_t *pv)
  { r2_t p0, p1, p2, p3;
    haf_draw_get_edge_bezier_points(a, &p0, &p1, &p2, &p3);
    r2_bezier_eval(0, 1, &p0, &p1, &p2, &p3, t, p, pu);
    (void)r2_dir(pu, pu);
    if ((pu != NULL) && (pv != NULL))
      { (*pv) = (r2_t){{ -pu->c[1], +pu->c[0] }}; }
  }

/* DATA STRUCTURE VIEW */

void haf_draw_arc_record(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a)
  { if (a->org->show)
      { r2_t p, pu, pv;
        double tm = haf_draw_arc_record_TIME;
        haf_draw_get_arc_point(a, tm, &p, &pu, &pv);
        double hu = haf_draw_arc_record_HU; 
        double hv = haf_draw_arc_record_HV;
        r2_t qu, qv;
        r2_scale(hu, &pu, &qu);
        r2_scale(hv, &pv, &qv);
        epswr_parallelogram
          ( eps,
            p.c[0], p.c[1],
            qu.c[0], qu.c[1],
            qv.c[0], qv.c[1],
            TRUE, TRUE
          );
      }
  }

void haf_draw_vert_record(epswr_figure_t *eps, haf_draw_data_t *dd, uint64_t kv)
  {
    if (v->show)
      { double rad = haf_draw_vert_record_RAD;
        epswr_dot(eps, dd, v->pos.c[0], v->pos.c[1], rad, TRUE, TRUE);
      }
  }

void haf_draw_face_record(epswr_figure_t *eps, haf_draw_data_t *dd, haf_face_t *f)
  {
    if (f->show)
      { double rad = haf_draw_face_record_RAD;
        epswr_dot(eps, dd, f->ctr.c[0], f->ctr.c[1], rad, TRUE, TRUE);
      }
  }

void haf_draw_arc_next_link(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a)
  { haf_arc_t b = a->next;
    if (a->org->show)
      { /* Link exits from top left quarter of {a}'s box: */
        r2_t p, pu, pv;
        double tm = haf_draw_arc_record_TIME;
        haf_draw_get_arc_point(a, tm, &p, &pu, &pv);
        r2_t pp = haf_draw_get_arc_record_point(&p, &pu, &pv, +0.5, +0.5);
        
        /* Link goes to left side of {b}'s box, near bottom: */
        r2_t q, qu, qv;
        haf_draw_get_arc_point(b, tm, &q, &qu, &qv);
        r2_t qq = haf_draw_get_arc_record_point(&q, &qu, &qv, -0.7, +1.0);
        
        double bend = -0.5;
        /* Cut in half if other end is invisible: */
        if (! b->org->show) { r2_mix(0.5, &pp, 0.5, &qq, &qq); bend = bend/4; }
        haf_draw_link(eps, dd, &pp, &qq, bend);
      }
  }

void haf_draw_arc_sym_link(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a)
  { haf_arc_t b = a->sym;
    if (a->org->show)
      { /* Link exits from top right quarter of {a}'s box: */
        r2_t p, pu, pv;
        double tm = haf_draw_arc_record_TIME;
        haf_draw_get_arc_point(a, tm, &p, &pu, &pv);
        r2_t pp = haf_draw_get_arc_record_point(&p, &pu, &pv, +0.5, -0.5);
        
        /* Link goes near bottom right corner of {b}'s box: */
        r2_t q, qu, qv;
        haf_draw_get_arc_point(b, tm, &q, &qu, &qv);
        r2_t qq = haf_draw_get_arc_record_point(&q, &qu, &qv, +1.0, 00.0);
        
        double bend = -0.5;
        haf_draw_link(eps, dd, &pp, &qq, bend);
      }
  }

void haf_draw_arc_org_link(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a)
  { uint64_t kv = a->org;
    if (v->show)
      { /* Link exits from bottom half of {a}'s box: */
        r2_t p, pu, pv;
        double tm = haf_draw_arc_record_TIME;
        haf_draw_get_arc_point(a, tm, &p, &pu, &pv);
        r2_t pp = haf_draw_get_arc_record_point(&p, &pu, &pv, -0.5, -0.5);
        
        /* Link goes to edge of {v}'s box: */
        r2_t q = (r2_t){{ v->pos.c[0], v->pos.c[1] }};
        r2_t qq = haf_draw_get_vert_record_point(&p, &q);
        
        double bend = 00.0;
        haf_draw_link(eps, dd, &pp, &qq, bend);
      }
  
  }

void haf_draw_arc_left_link(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a)
  {
    haf_face_t *f = a->left;
    if (a->org->show)
      { /* Make link starts at bot left quarter of {a}'s box: */
        r2_t p, pu, pv; 
        double tm = haf_draw_arc_record_TIME;
        haf_draw_get_arc_point(a, tm, &p, &pu, &pv);
        r2_t pp = haf_draw_get_arc_record_point(&p, &pu, &pv, -0.5, +0.5);
    
        /* Link ends on edge of {f}'s box: */
        r2_t q = (r2_t){{ f->ctr.c[0], f->ctr.c[1] }};
        r2_t qq = haf_draw_get_face_record_point(&p, &q);
        double bend = 00.0;
        /* Cut in half if other end is invisible: */
        if (! f->show) { r2_mix(0.5, &pp, 0.5, &qq, &qq); bend = bend/4; }
        haf_draw_link(eps, dd, &pp, &qq, bend);
      }
  }

void haf_draw_vert_out_link(epswr_figure_t *eps, haf_draw_data_t *dd, uint64_t kv)
  {
    haf_arc_t a = v->out;
    if (v->show)
      { /* Link ends on bootom right corner of {a}'s box: */
        r2_t q, qu, qv; 
        double tm = haf_draw_arc_record_TIME;
        haf_draw_get_arc_point(a, tm, &q, &qu, &qv);
        r2_t qq = haf_draw_get_arc_record_point(&q, &qu, &qv, -1.0, -1.0);
    
        /* Link starts at center of vertex's box: */
        r2_t p = (r2_t){{ v->pos.c[0], v->pos.c[1] }};
        double bend = -0.075*r2_dist(&p, &qq);
        haf_draw_link(eps, dd, &p, &qq, bend);
      }
  }

void haf_draw_face_side_link(epswr_figure_t *eps, haf_draw_data_t *dd, haf_face_t *f)
  {
    haf_arc_t a = f->side;
    if (f->show)
      { /* Link ends on top quarter of left side of {a}'s box: */
        r2_t q, qu, qv; 
        double tm = haf_draw_arc_record_TIME;
        haf_draw_get_arc_point(a, tm, &q, &qu, &qv);
        r2_t qq = haf_draw_get_arc_record_point(&q, &qu, &qv, +0.5, +1.00);
    
        /* Link starts at center of face's box: */
        r2_t p = (r2_t){{ f->ctr.c[0], f->ctr.c[1] }};
        double bend = +0.075*r2_dist(&p, &qq);
        /* Cut in half if arc box is invisible: */
        if (! a->org->show) { r2_mix(0.5, &p, 0.5, &qq, &qq); bend = bend/4; }
        haf_draw_link(eps, dd, &p, &qq, bend);
      }
  }

r2_t haf_draw_get_arc_record_point(r2_t *p, r2_t *pu, r2_t *pv, double u, double v)
  { 
    r2_t q; r2_mix(u*haf_draw_arc_record_HU, pu, v*haf_draw_arc_record_HV, pv, &q);
    r2_add(p, &q, &q);
    return q;
  }

r2_t haf_draw_get_vert_record_point(r2_t *p, r2_t *q)
  { double d = r2_dist(p, q);
    double rad = haf_draw_vert_record_RAD;
    double s = rad/d;
    r2_t qq;
    r2_mix(s, p, 1-s, q, &qq);
    return qq;
  }

r2_t haf_draw_get_face_record_point(r2_t *p, r2_t *q)
  { double d = r2_dist(p, q);
    double rad = haf_draw_face_record_RAD;
    double s = rad/d;
    r2_t qq;
    r2_mix(s, p, 1-s, q, &qq);
    return qq;
  }

void haf_draw_link(epswr_figure_t *eps, haf_draw_data_t *dd, r2_t *p, r2_t *q, double bend)
  { 
    r2_t p1, p2;
    r2_bezier_from_bend(p, q, bend, &p1, &p2);
    epswr_curve
      ( eps, 
        p->c[0], p->c[1], 
        p1.c[0], p1.c[1], 
        p2.c[0], p2.c[1], 
        q->c[0], q->c[1]
      );
    double rad = haf_draw_link_dot_RAD;
    epswr_dot(eps, dd, p->c[0], p->c[1], rad, TRUE, TRUE);
    double ahwid = haf_draw_link_arrow_WIDTH;
    double ahlen = haf_draw_link_arrow_LENGTH;
    epswr_arrowhead
      ( eps, p2.c[0], p2.c[1], q->c[0], q->c[1],
        ahwid, ahlen, 1.0, 
        TRUE, TRUE
      );
  }

/* LABELS */

void haf_draw_face_label(epswr_figure_t *eps, haf_draw_data_t *dd, haf_face_t *f, char *lab, double dx, double dy)
  { r2_t ctr = (r2_t){{ f->ctr.c[0], f->ctr.c[1] }};
    double hAlign, vAlign;
    r2_t dsp = (r2_t){{ dx, dy }};
    haf_draw_choose_lab_align(&dsp, &hAlign, &vAlign);
    r2_t q; r2_add(&dsp, &ctr, &q);
    epswr_label(eps, dd, lab, "Rg", q.c[0], q.c[1], 0, TRUE, hAlign, vAlign, TRUE, FALSE);
  }
  
void haf_draw_vert_label(epswr_figure_t *eps, haf_draw_data_t *dd, uint64_t kv, char *lab, double dx, double dy)
  { r2_t pos = (r2_t){{ v->pos.c[0], v->pos.c[1] }};
    double hAlign, vAlign;
    r2_t dsp = (r2_t){{ dx, dy }};
    haf_draw_choose_lab_align(&dsp, &hAlign, &vAlign);
    r2_t q; r2_add(&dsp, &pos, &q);
    epswr_label(eps, dd, lab, "Rg", q.c[0], q.c[1], 0, TRUE, hAlign, vAlign, TRUE, FALSE);
  }

void haf_draw_arc_label(epswr_figure_t *eps, haf_draw_data_t *dd, haf_arc_t a, char *lab, double du, double dv)
  { r2_t p, pu, pv;
    double tm = haf_draw_arc_record_TIME;
    haf_draw_get_arc_point(a, tm, &p, &pu, &pv);
    r2_t dp; r2_mix(du, &pu, dv, &pv, &dp);
    double hAlign, vAlign;
    haf_draw_choose_lab_align(&dp, &hAlign, &vAlign);
    r2_t q; r2_add(&dp, &p, &q);
    epswr_label(eps, dd, lab, "Rg", q.c[0], q.c[1], 0, TRUE, hAlign, vAlign, TRUE, FALSE);
  }

void haf_draw_choose_lab_align(r2_t *dp, double *hAlign_P, double *vAlign_P)
  { 
    (*hAlign_P) = (dp->c[0] > 0 ? 0.0 : (dp->c[0] < 0 ? 1.0 : 0.5));
    (*vAlign_P) = (dp->c[1] > 0 ? 0.0 : (dp->c[1] < 0 ? 1.0 : 0.5));
  }
