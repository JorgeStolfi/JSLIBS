/* See hr2.h */
/* Last edited on 2011-12-21 20:30:01 by stolfi */ 

/* Based on HR2.m3, created 1994-05-04 by J. Stolfi. */

#define _GNU_SOURCE
#include <hr2.h>

#include <r2.h>
#include <r3.h>
#include <rn.h>
#include <r3x3.h>
#include <affirm.h>
#include <sign.h>
#include <sign_get.h>

#include <math.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

hr2_point_t hr2_from_r2(r2_t *c)
  {
    return (hr2_point_t){{{1.0, c->c[0], c->c[1]}}};
  }

r2_t r2_from_hr2(hr2_point_t *p)
  {
    double w = p->c.c[0];
    demand(w != 0.0, "point at infinity");
    return (r2_t){{p->c.c[1]/w, p->c.c[2]/w}};
  }

double hr2_pt_pt_diff(hr2_point_t *p, hr2_point_t *q)
  {
    return r3_angle(&(p->c), &(q->c));
  }

sign_t hr2_side(hr2_point_t *p, hr2_line_t *L)
  {
    double dd = r3_dot(&(p->c), &(L->f));
    return sign_double(dd);
  }

sign_t hr2_orient(hr2_point_t *p, hr2_point_t *q, hr2_point_t *r)
  {
    double dd = r3_det(&(p->c), &(q->c), &(r->c));
    return sign_double(dd);
  }

hr2_line_t hr2_join(hr2_point_t *p, hr2_point_t *q)
  {
    r3_t f;
    r3_cross(&(p->c), &(q->c), &f);
    return (hr2_line_t){f};
  }

hr2_point_t hr2_meet(hr2_line_t *K, hr2_line_t *L)
  {
    r3_t c;
    r3_cross(&(K->f), &(L->f), &c);
    return (hr2_point_t){c};
  }

bool_t hr2_pmap_is_identity(hr2_pmap_t *m)
  {
    return r3x3_is_unif_scaling(&(m->dir), m->dir.c[0][0]); 
  }

hr2_point_t hr2_map_point(hr2_point_t *p, hr2_pmap_t *m)
  {
    r3_t q; 
    r3x3_map_row(&(p->c), &(m->dir), &q);
    return (hr2_point_t){q};
  }

hr2_line_t hr2_map_line(hr2_line_t *L, hr2_pmap_t *m)
  {
    r3_t f;
    r3x3_map_col(&(m->inv), &(L->f), &f);
    return (hr2_line_t){f};
  }

hr2_pmap_t hr2_comp_map(hr2_pmap_t *m, hr2_pmap_t *n)
  {
    r3x3_t dir, inv;
    r3x3_mul(&(m->dir), &(n->dir), &dir);
    r3x3_mul(&(n->inv), &(m->inv), &inv);
    return (hr2_pmap_t){dir, inv};
  }

hr2_pmap_t hr2_inv_map(hr2_pmap_t *m)
  {
    return (hr2_pmap_t){m->inv, m->dir};
  }

r2_t hr2_line_dir(hr2_line_t *L)
  {
    double dx = L->f.c[2];
    double dy = -L->f.c[1];
    double length = hypot(dx, dy);
    return (r2_t){{dx/length, dy/length}};
  }

r2_t hr2_line_normal(hr2_line_t *L)
  {
    double nx = L->f.c[1];
    double ny = L->f.c[2];
    double length = hypot(nx, ny);
    return (r2_t){{nx/length, ny/length}};
  }

r2_t hr2_point_point_dir(hr2_point_t *frm, hr2_point_t *tto)
  {
    double fw = frm->c.c[0];
    double tw = tto->c.c[0];
    
    double fx = frm->c.c[1];
    double tx = tto->c.c[1];
    double dx = fw * tx - tw * fx;
    
    double fy = frm->c.c[2];
    double ty = tto->c.c[2];
    double dy = fw * ty - tw * fy;
    
    double length = hypot(dx, dy);
    return (r2_t){{dx/length, dy/length}};
  }

hr2_pmap_t hr2_pmap_from_points(hr2_point_t *p, hr2_point_t *q, hr2_point_t *r, hr2_point_t *u)
  {
    int i;
    hr2_pmap_t m; /* The resulting map. */
    
    /* Compute weights {(a,b,c)=w.c[0..2]} for rows of {Q} that map {[1,1,1]} to {u}: */
    r3_t w;
    { /* Compute a matrix {Q} that maps the cardinal points to {p,q,r} as given: */
      r3x3_t Q;
      Q.c[0][0] = p->c.c[0]; Q.c[0][1] = p->c.c[1]; Q.c[0][2] = p->c.c[2];
      Q.c[1][0] = q->c.c[0]; Q.c[1][1] = q->c.c[1]; Q.c[1][2] = q->c.c[2];
      Q.c[2][0] = r->c.c[0]; Q.c[2][1] = r->c.c[1]; Q.c[2][2] = r->c.c[2];
      r3x3_inv(&Q, &Q);
      w.c[0] = u->c.c[0];  w.c[1] = u->c.c[1];  w.c[2] = u->c.c[2];
      r3x3_map_row(&w, &Q, &w);
    }
    
    /* Make the weights positive, so that {p,q,r} are strictly honored: */
    for (i = 0; i < NH; i++) { w.c[i] = fabs(w.c[i]); }

    /* Ensure that {m.dir} maps the cardinal points to {p,q,r} and some unit point to {u}: */
    r3x3_t *P = &(m.dir);
    P->c[0][0] = w.c[0]*p->c.c[0]; P->c[0][1] = w.c[0]*p->c.c[1]; P->c[0][2] = w.c[0]*p->c.c[2];
    P->c[1][0] = w.c[1]*q->c.c[0]; P->c[1][1] = w.c[1]*q->c.c[1]; P->c[1][2] = w.c[1]*q->c.c[2];
    P->c[2][0] = w.c[2]*r->c.c[0]; P->c[2][1] = w.c[2]*r->c.c[1]; P->c[2][2] = w.c[2]*r->c.c[2];

    /* Compute the inverse map: */
    r3x3_inv(&(m.dir), &(m.inv));

    return m;

  }
