/* See hr2.h */
/* Last edited on 2024-11-20 11:52:34 by stolfi */ 

/* Based on HR2.m3 created 1994-05-04 by J. Stolfi. */

#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include <sign_get.h>
#include <sign.h>
#include <affirm.h>

#include <r2x2.h>
#include <r3x3.h>
#include <rn.h>
#include <r3.h>
#include <r2.h>


#include <hr2.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

r2_t r2_from_hr2_nan(hr2_point_t *p)
  {
    double w = p->c.c[0];
    double m = fmax(fabs(p->c.c[1]), fabs(p->c.c[2]));
    r2_t c;
    if (fabs(w) <= m*1e-200) 
      { c = (r2_t){{ NAN, NAN }}; }
    else
      { c = (r2_t){{ p->c.c[1]/w, p->c.c[2]/w }}; }
    return c;
  }

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

hr2_point_t hr2_point_at_infinity(r2_t *dir)
  { hr2_point_t p;
    p.c.c[0] = 0.0;
    p.c.c[1] = dir->c[0];
    p.c.c[2] = dir->c[1];
    return p;
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
  
double hr2_dist(hr2_point_t *a, hr2_point_t *b)
  {
    double aw = 1.0/a->c.c[0];
    double bw = 1.0/b->c.c[0];
    
    double ax = a->c.c[1];
    double bx = b->c.c[1];
    double dx = ax*aw - bx*bw;
    
    double ay = a->c.c[2];
    double by = b->c.c[2];
    double dy = ay*aw - by*bw;
    return hypot(dx, dy);
  }
    
double hr2_dist_sqr(hr2_point_t *a, hr2_point_t *b)
  {
    double aw = 1.0/a->c.c[0];
    double bw = 1.0/b->c.c[0];
    
    double ax = a->c.c[1];
    double bx = b->c.c[1];
    double dx = ax*aw - bx*bw;
    
    double ay = a->c.c[2];
    double by = b->c.c[2];
    double dy = ay*aw - by*bw;
    return dx*dx + dy*dy;
  }

hr2_point_t hr2_point_mix(double pt, hr2_point_t *p, double qt, hr2_point_t *q)
  { hr2_point_t r;
    r.c.c[0] = pt * p->c.c[0] + qt * q->c.c[0];
    r.c.c[1] = pt * p->c.c[1] + qt * q->c.c[1];
    r.c.c[2] = pt * p->c.c[2] + qt * q->c.c[2];
    return r;
  }

void hr2_L_inf_normalize_point(hr2_point_t *p)
  { double m = r3_L_inf_norm(&(p->c));
    if (m != 1.0)
      {	p->c.c[0] /= m;
	p->c.c[1] /= m;
	p->c.c[2] /= m;
      }
  }
  
void hr2_L_inf_normalize_line(hr2_line_t *L)
  { double m = r3_L_inf_norm(&(L->f));
    if (m != 1.0)
      {	L->f.c[0] /= m;
	L->f.c[1] /= m;
	L->f.c[2] /= m;
      }
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

hr2_point_t hr2_point_throw(void)
  { hr2_point_t p;
    r3_throw_dir(&(p.c));
    return p;
  }

hr2_line_t hr2_line_throw(void)
  { hr2_line_t A;
    r3_throw_dir(&(A.f));
    return A;
  }

void hr2_point_print(FILE *wr, char *pre, hr2_point_t *p, char *fmt, char *suf)
  { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    fputs("[ ", wr);
    if (fmt == NULL) { fmt = "24.26e"; }
    for (uint32_t i = 0;  i < NH; i++)
      { if (i != 0) { fputc(' ', wr); }
        fprintf(wr, fmt, p->c.c[i]);
      }
    fputs(" ]", wr);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }

void hr2_line_print(FILE *wr, char *pre, hr2_line_t *L, char *fmt, char *suf)
  { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, wr); }
    fputs("< ", wr);
    if (fmt == NULL) { fmt = "24.26e"; }
    for (uint32_t i = 0;  i < NH; i++)
      { if (i != 0) { fputc(' ', wr); }
        fprintf(wr, fmt, L->f.c[i]);
      }
    fputs(" >", wr);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, wr); }
  }
