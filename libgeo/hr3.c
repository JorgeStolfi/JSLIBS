/* See hr3.h */
/* Last edited on 2023-10-09 21:23:51 by stolfi */ 

#define _GNU_SOURCE
#include <stdint.h>
#include <hr3.h>

#include <r3.h>
#include <r4.h>
#include <r6.h>
#include <r3x3.h>
#include <r4x4.h>
#include <affirm.h>
#include <sign.h>
#include <sign_get.h>

#include <math.h>

/* Based on HR3.m3, created 1993-05-18 by Marcos C. Carrard. */
/* Based on H3.pas by J. Stolfi. */

#define NH 4
  /* Number of homogeneous coordinates in a point or coefficients in a plane. */

#define NG 6
  /* Number of Grassman coordinates of a line. */

#define NC 3
  /* Number of Cartesian coordinates in a point. */

void hr3_point_to_r3_nan(hr3_point_t *p, r3_t *c);
  /* Converts a point from homogeneous coordinates {p} to Cartesian
    coordinates {c}.  If {p} is at infinity, returns {(NAN,NAN,NAN)}. */

void hr3_point_to_r3_nan(hr3_point_t *p, r3_t *c)
  {
    double w = p->c.c[0];
    double m = fmax(fmax(fabs(p->c.c[1]), fabs(p->c.c[2])), fabs(p->c.c[3]));
    if (fabs(w) <= m*1e-200) 
      { (*c) = (r3_t){{ NAN, NAN, NAN }}; }
    else
      { (*c) = (r3_t){{ p->c.c[1]/w, p->c.c[2]/w, p->c.c[3]/w }}; } 
  }

hr3_point_t hr3_from_r3(r3_t *c)
  {
    return (hr3_point_t){{{1.0, c->c[0], c->c[1], c->c[2]}}};
  }
    
r3_t r3_from_hr3(hr3_point_t *p)
  {
    double w = p->c.c[0];
    demand(w != 0.0, "point at infinity");
    return (r3_t){{ p->c.c[1]/w, p->c.c[2]/w, p->c.c[3]/w }};
  }

double hr3_pt_pt_diff(hr3_point_t *p, hr3_point_t *q)
  {
    return r4_angle(&(p->c), &(q->c));
  }

sign_t hr3_side(hr3_point_t *p, hr3_plane_t *Q)
  {
    double dd = r4_dot(&(p->c), &(Q->f));
    return sign_double(dd);
  }       

sign_t hr3_orient(hr3_point_t *p, hr3_point_t *q, hr3_point_t *r, hr3_point_t *s)
  {
    double dd = r4_det(&(p->c), &(q->c), &(r->c), &(s->c));
    return sign_double(dd);
  }
  
hr3_line_t hr3_line_from_two_points(hr3_point_t *p, hr3_point_t *q)
  {
    double a0 = p->c.c[0];
    double a1 = p->c.c[1];
    double a2 = p->c.c[2];
    double a3 = p->c.c[3];
                      
    double b0 = q->c.c[0];
    double b1 = q->c.c[1];
    double b2 = q->c.c[2];
    double b3 = q->c.c[3];
    return 
      (hr3_line_t)
        {{{ a0*b1 - a1*b0,
            a0*b2 - a2*b0,
            a1*b2 - a2*b1,
            a0*b3 - a3*b0,
            a1*b3 - a3*b1,
            a2*b3 - a3*b2
        }}};
  }
  
hr3_plane_t hr3_plane_from_three_points(hr3_point_t *p, hr3_point_t *q, hr3_point_t *r)
  {
    r4_t f;
    r4_cross(&(p->c), &(q->c), &(r->c), &f);
    return (hr3_plane_t){f};
  }   

hr3_plane_t hr3_plane_from_line_and_point(hr3_line_t *L, hr3_point_t *r)
  {
    double a01 = L->k.c[0];
    double a02 = L->k.c[1];
    double a12 = L->k.c[2];
    double a03 = L->k.c[3];
    double a13 = L->k.c[4];
    double a23 = L->k.c[5];
    
    double b0 = r->c.c[0];
    double b1 = r->c.c[1];
    double b2 = r->c.c[2];
    double b3 = r->c.c[3];
    return 
      (hr3_plane_t)
        {{{
            - a12*b3 + a13*b2 - a23*b1,
            + a02*b3 - a03*b2 + a23*b0,
            - a01*b3 + a03*b1 - a13*b0,
            + a01*b2 - a02*b1 + a12*b0
        }}};
  }

hr3_line_t hr3_line_from_two_planes(hr3_plane_t *P, hr3_plane_t *Q)
  {
    double a0 = P->f.c[0];
    double a1 = P->f.c[1];
    double a2 = P->f.c[2];
    double a3 = P->f.c[3];
                      
    double b0 = Q->f.c[0];
    double b1 = Q->f.c[1];
    double b2 = Q->f.c[2];
    double b3 = Q->f.c[3];
    return 
      (hr3_line_t)
        {{{
            + a2*b3 - a3*b2,
            - a1*b3 + a3*b1,
            + a0*b3 - a3*b0,
            + a1*b2 - a2*b1,
            - a0*b2 + a2*b0,
            + a0*b1 - a1*b0
        }}};
  }
  
hr3_point_t hr3_point_from_three_planes(hr3_plane_t *P, hr3_plane_t *Q, hr3_plane_t *R)
  {
    r4_t c;
    r4_cross(&(R->f), &(Q->f), &(P->f), &c);
    return (hr3_point_t){c};
  } 
  
hr3_point_t hr3_point_from_line_and_plane(hr3_line_t *L, hr3_plane_t *R)
  {
    double a01 = L->k.c[0];
    double a02 = L->k.c[1];
    double a12 = L->k.c[2];
    double a03 = L->k.c[3];
    double a13 = L->k.c[4];
    double a23 = L->k.c[5];
     
    double b123 = R->f.c[0];
    double b023 = R->f.c[1];
    double b013 = R->f.c[2];
    double b012 = R->f.c[3];
    return 
      (hr3_point_t)
        {{{
            + a01*b023 + a02*b013 + a03*b012,
            - a01*b123 + a12*b013 + a13*b012,
            - a02*b123 - a12*b023 + a23*b012,
            - a03*b123 - a13*b023 - a23*b013
        }}};
  }

bool_t hr3_pmap_is_identity(hr3_pmap_t *M)
  {
    return r4x4_is_unif_scaling(&(M->dir), M->dir.c[0][0]); 
  }

hr3_point_t hr3_pmap_point(hr3_point_t *p, hr3_pmap_t *M)
  {
    r4_t c;
    r4x4_map_row(&(p->c), &(M->dir), &c);
    return (hr3_point_t){c};
  }
        
hr3_point_t hr3_pmap_inv_point(hr3_point_t *p, hr3_pmap_t *M)
  {
    r4_t c;
    r4x4_map_row(&(p->c), &(M->inv), &c);
    return (hr3_point_t){c};
  }
        
hr3_plane_t hr3_pmap_plane(hr3_plane_t *P, hr3_pmap_t *M)
  {
    r4_t f;
    r4x4_map_col(&(M->inv), &(P->f), &f);
    return (hr3_plane_t){f};
  }
    
hr3_plane_t hr3_pmap_inv_plane(hr3_plane_t *P, hr3_pmap_t *M)
  {
    r4_t f;
    r4x4_map_col(&(M->dir), &(P->f), &f);
    return (hr3_plane_t){f};
  }

r3_t hr3_pmap_r3_point(r3_t *p, hr3_pmap_t *M)
  {
    hr3_point_t ph = (hr3_point_t){{{ 1.0, p->c[0], p->c[1], p->c[2] }}};
    hr3_point_t qh;
    r4x4_map_row(&(ph.c), &(M->dir), &(qh.c));
    r3_t qc; hr3_point_to_r3_nan(&qh, &qc);
    return qc;
  }

r3_t hr3_pmap_inv_r3_point(r3_t *p, hr3_pmap_t *M)
  {
    hr3_point_t ph = (hr3_point_t){{{ 1, p->c[0], p->c[1], p->c[2] }}};
    hr3_point_t qh;
    r4x4_map_row(&(ph.c), &(M->inv), &(qh.c));
    r3_t qc; hr3_point_to_r3_nan(&qh, &qc);
    return qc;
  }

hr3_pmap_t hr3_pmap_translation(r3_t *v)
  {
    hr3_pmap_t M; 
    r4x4_ident(&(M.dir));
    M.dir.c[0][1] = +v->c[0];
    M.dir.c[0][2] = +v->c[1];
    M.dir.c[0][3] = +v->c[2];
    r4x4_ident(&(M.inv));
    M.inv.c[0][1] = -v->c[0];
    M.inv.c[0][2] = -v->c[1];
    M.inv.c[0][3] = -v->c[2];
    return M;
  }

hr3_pmap_t hr3_pmap_scaling(r3_t *scale)
  {
    hr3_pmap_t M;
    for (int32_t i = 0; i < NH; i++)
      { for (int32_t j = 0; j < NH; j++)
          { if ((i == 0) && (j == 0))
              { M.dir.c[i][j] = M.inv.c[i][j] = 1.0; }
            else if (i == j)
              { M.dir.c[i][j] = scale->c[j-1];
                M.inv.c[i][j] = 1/scale->c[j-1];
              }
            else
              { M.dir.c[i][j] = M.inv.c[i][j] = 0.0; }
          }
      }
    return M;
  
  }

hr3_pmap_t hr3_pmap_u_v_rotation(r3_t *u, r3_t *v)
  { /* Compute the matrix {r} for {R^3} to {R^3}: */
    r3x3_t r;
    r3x3_u_v_rotation(u, v, &r);
    
    /* Convert to a projective map (note that inverse is just transpose): */
    hr3_pmap_t M;
    for (int32_t i = 0; i < NH; i++)
      { for (int32_t j = 0; j < NH; j++)
          { if ((i == 0) && (j == 0))
              { M.dir.c[i][j] = M.inv.c[i][j] = 1.0; }
            else if ((i == 0) || (j == 0))
              { M.dir.c[i][j] = M.inv.c[i][j] = 0.0; }
            else
              { M.dir.c[i][j] = r.c[i-1][j-1]; 
                M.inv.c[i][j] = r.c[j-1][i-1]; 
              }
          }
      }
    return M;
  }

hr3_pmap_t hr3_pmap_compose(hr3_pmap_t *M, hr3_pmap_t *N)
  {
    r4x4_t dir, inv; 
    r4x4_mul(&(M->dir), &(N->dir), &dir);
    r4x4_mul(&(N->inv), &(M->inv), &inv);
    return (hr3_pmap_t){dir, inv};
  }    
  
hr3_pmap_t hr3_pmap_inv(hr3_pmap_t *M)
  {
    return (hr3_pmap_t){M->inv, M->dir};
  }

hr3_pmap_t hr3_pmap_inv_compose(hr3_pmap_t *M, hr3_pmap_t *N)
  {
    r4x4_t dir, inv;
    r4x4_mul(&(M->inv), &(N->dir), &dir);
    r4x4_mul(&(N->inv), &(M->dir), &inv);
    return (hr3_pmap_t){dir, inv};
  }

r3_t hr3_line_dir(hr3_line_t *L)
  {
    double a01 = L->k.c[0];
    double a02 = L->k.c[1];
    double a03 = L->k.c[3];
    r3_t d = (r3_t){{a01, a02, a03}};
    r3_dir(&d, &d);
    return d;
  }

r3_t hr3_plane_normal(hr3_plane_t *P)
  {
    double nx = P->f.c[1];
    double ny = P->f.c[2];
    double nz = P->f.c[3];
    double length = hypot(hypot(nx, ny), nz);
    return (r3_t){{nx/length, ny/length, nz/length}};
  }    

r3_t hr3_point_point_dir(hr3_point_t *frm, hr3_point_t *tto)
  {
    double fw = frm->c.c[0];
    double tw = tto->c.c[0];
    
    double fx = frm->c.c[1];
    double tx = tto->c.c[1];
    double dx = fw * tx - tw * fx;
    
    double fy = frm->c.c[2];
    double ty = tto->c.c[2];
    double dy = fw * ty - tw * fy;
    
    double fz = frm->c.c[3];
    double tz = tto->c.c[3];
    double dz = fw * tz - tw * fz;
    
    double length = hypot(hypot(dx, dy), dz);
    
    return (r3_t){{dx/length, dy/length, dz/length}};
  }  
  
double hr3_dist(hr3_point_t *a, hr3_point_t *b)
  {
    double aw = 1.0/a->c.c[0];
    double bw = 1.0/b->c.c[0];
    
    double ax = a->c.c[1];
    double bx = b->c.c[1];
    double dx = ax*aw - bx*bw;
    
    double ay = a->c.c[2];
    double by = b->c.c[2];
    double dy = ay*aw - by*bw;
    
    double az = a->c.c[3];
    double bz = b->c.c[3];
    double dz = az*aw - bz*bw;
    return hypot(hypot(dx, dy), dz);
  }
    
double hr3_dist_sqr(hr3_point_t *a, hr3_point_t *b)
  {
    double aw = 1.0/a->c.c[0];
    double bw = 1.0/b->c.c[0];
    
    double ax = a->c.c[1];
    double bx = b->c.c[1];
    double dx = ax*aw - bx*bw;
    
    double ay = a->c.c[2];
    double by = b->c.c[2];
    double dy = ay*aw - by*bw;
    
    double az = a->c.c[3];
    double bz = b->c.c[3];
    double dz = az*aw - bz*bw;
    return dx * dx + dy * dy + dz * dz;
  }

hr3_point_t hr3_point_at_infinity(r3_t *dir)
  { hr3_point_t p;
    p.c.c[0] = 0.0;
    p.c.c[1] = dir->c[0];
    p.c.c[2] = dir->c[1];
    p.c.c[3] = dir->c[2];
    return p;
  }

hr3_point_t hr3_point_mix(double pt, hr3_point_t *p, double qt, hr3_point_t *q)
  { hr3_point_t r;
    r.c.c[0] = pt * p->c.c[0] + qt * q->c.c[0];
    r.c.c[1] = pt * p->c.c[1] + qt * q->c.c[1];
    r.c.c[2] = pt * p->c.c[2] + qt * q->c.c[2];
    r.c.c[3] = pt * p->c.c[3] + qt * q->c.c[3];
    return r;
  }

void hr3_L_inf_normalize_point(hr3_point_t *p)
  { double m = r4_L_inf_norm(&(p->c));
    if (m != 1.0)
      {	p->c.c[0] /= m;
	p->c.c[1] /= m;
	p->c.c[2] /= m;
	p->c.c[3] /= m;
      }
  }
  
void hr3_L_inf_normalize_plane(hr3_plane_t *P)
  { double m = r4_L_inf_norm(&(P->f));
    if (m != 1.0)
      {	P->f.c[0] /= m;
	P->f.c[1] /= m;
	P->f.c[2] /= m;
	P->f.c[3] /= m;
      }
  }

void hr3_L_inf_normalize_line(hr3_line_t *L)
  { double m = r6_L_inf_norm(&(L->k));
    if (m != 1.0)
      {	L->k.c[0] /= m;
	L->k.c[1] /= m;
	L->k.c[2] /= m;
	L->k.c[3] /= m;
	L->k.c[4] /= m;
	L->k.c[5] /= m;
      }
  }

hr3_pmap_t hr3_pmap_aff_from_mat_and_disp(r3x3_t *E, r3_t *d)

  {
    hr3_pmap_t M;

    r3x3_t F; r3x3_inv(E, &(F));
    r3_t p;
    r3x3_map_row(d, &(F), &p);

    for (int32_t i = 0; i < NH; i++)
      { for (int32_t j = 0; j < NH; j++)
          { double Pij, Qij; /* Elements of {M.dir} and {M.inv}. */
            if (i == 0)
              { if (j == 0)
                  { Pij = 1.0; Qij = 1.0; }
                else
                  { Pij = d->c[j-1]; Qij = p.c[j-1]; }
              }
            else
              { if (j == 0)
                  { Pij = 0.0; Qij = 0.0; }
                else
                  { Pij = E->c[i-1][j-1]; Qij = F.c[i-1][j-1]; }
              }
            M.dir.c[i][j] = Pij; 
            M.inv.c[i][j] = Qij;
          }
      }
    
    return M;
  }
  
hr3_pmap_t hr3_pmap_aff_from_four_points(r3_t *o, r3_t *p, r3_t *q, r3_t *r)
  {
    hr3_pmap_t M; 

    for (int32_t j = 0; j < NH; j++)
      { M.dir.c[0][j] = (j == 0 ? 1.0 : o->c[j-1]);
        M.dir.c[1][j] = (j == 0 ? 0.0 : p->c[j-1] - o->c[j-1]);
        M.dir.c[2][j] = (j == 0 ? 0.0 : q->c[j-1] - o->c[j-1]);
        M.dir.c[3][j] = (j == 0 ? 0.0 : r->c[j-1] - o->c[j-1]);
      }

    r4x4_inv(&(M.dir), &(M.inv));
    
    return M;
  }
  
hr3_pmap_t hr3_pmap_from_five_points(hr3_point_t *p, hr3_point_t *q, hr3_point_t *r, hr3_point_t *s, hr3_point_t *u)
  {
    hr3_pmap_t M; /* The resulting map. */
    
    /* Compute weights {(a,b,c)=w.c[0..2]} for rows of {Q} that map {[1,1,1]} to {u}: */
    r4_t w;
    { /* Compute a matrix {Q} that maps the cardinal points to {p,q,r,s} as given: */
      r4x4_t Q;
      for (int32_t j = 0; j < NH; j++)
        { M.dir.c[0][j] = p->c.c[j];
          M.dir.c[1][j] = q->c.c[j];
          M.dir.c[2][j] = r->c.c[j];
          M.dir.c[3][j] = s->c.c[j];
        }
      /* Map {u} by the inverse of {Q}: */
      r4x4_inv(&(M.dir), &Q);
      w = u->c;
      r4x4_map_row(&w, &Q, &w);
    }
    
    /* Make the weights positive, so that {p,q,r,s} are strictly honored: */
    for (int32_t i = 0; i < NH; i++) { w.c[i] = fabs(w.c[i]); }

    /* Ensure that {M.dir} maps the cardinal points to {p,q,r,s} and some unit point to {u}: */
    for (int32_t i = 0; i < NH; i++)
      { for (int32_t j = 0; j < NH; j++)
          { M.dir.c[i][j] *= w.c[i];  }
      }

    /* Compute the inverse map: */
    r4x4_inv(&(M.dir), &(M.inv));

    return M;
  }

hr3_pmap_t hr3_pmap_persp(hr3_point_t *obs, hr3_point_t *foc, double rad, hr3_point_t *upp)
  { hr3_pmap_t M;
    r4x4_t Mt, Mr, Mc, Mrt;
    
    demand(foc->c.c[0] > 0.0, "focus must be finite and hither");
    demand(obs->c.c[0] >= 0.0, "observer must be hither or infinite");
    demand(upp->c.c[0] >= 0.0, "zenith must be hither or infinite"); 

    /* Start with a translation from {foc} to the origin: */
    for (int32_t i = 0; i < NH; i++){
      for (int32_t j = 0; j < NH; j++){
        if (i == j)
          { Mt.c[i][j] = foc->c.c[0]; }
        else if (i == 0)
          { Mt.c[i][j] = -(foc->c.c[j]); }
        else
          { Mt.c[i][j] = 0.0; }
      }
    }
    
    /* Compute the image frame vectors {r,s,t}: */
    r3_t t = hr3_point_point_dir(foc, obs); /* Vector {t} points out of image towards {obs}. */
    r3_t u = hr3_point_point_dir(foc, upp); /* Direction from {foc} towards zenith. */
    r3_t v, s;
    r3_decomp(&u, &t, &v, &s);
    /* Zenith reference point {upp} must not be on imagesys Z axis: */
    demand(r3_norm_sqr(&s) >= 1.0e-12, "bad zenith"); 
    r3_dir(&s, &s); /* Vector {s} is the image's vertical axis. */
    r3_t r;
    r3_cross(&s, &t, &r);
    r3_dir(&r, &r); /* Vector {r} is the image's horizontal axis. */

    /* Append the rotation matrix that moves {r,s,t} to X,Y,Z: */
    Mr.c[0][0] = 1.0;
    for (int32_t i = 1; i < NH; i++)
      { Mr.c[0][i] = 0.0;
        Mr.c[i][0] = 0.0;
        Mr.c[i][1] = r.c[i-1];
        Mr.c[i][2] = s.c[i-1];
        Mr.c[i][3] = t.c[i-1];
      }
    r4x4_mul(&Mt, &Mr, &Mrt);

    /* Do we need a conical projection step? */
    if (obs->c.c[0] == 0.0)
      { /* Observer is at infinity; cilindrical projection. */
        M.dir = Mrt;
      }
    else
      { /* Observer is finite; add conical projection step. */
        double d = hr3_dist(foc, obs);
        double uno = -1.0;
        if (fabs(d) > 1.0) { uno = -1.0/d; d = 1.0; }
        for (int32_t i = 0; i < NH; i++)
          for (int32_t j = 0; j < NH; j++)
            { if (i == j)
                { Mc.c[i][j] = d; }
              else
                { Mc.c[i][j] = 0.0; }
            }
        Mc.c[3][0] = uno;
        r4x4_mul(&Mrt, &Mc, &(M.dir));
      }

    /* Asked for scaling? */
    if (rad > 0.0)
      { /* Combine {M} with a uniform scale of {1/rad}: */
        for (int32_t i = 0; i < NH; i++)
          for (int32_t j = 1; j < NH; j++)
            { M.dir.c[i][j] /= rad; }
      }

    /* Compute inverse matrix: */
    r4x4_inv(&M.dir, &M.inv);
    return M;
  }
  
hr3_point_t hr3_point_throw(void)
  { hr3_point_t p;
    r4_throw_dir(&(p.c));
    return p;
  }
   
hr3_plane_t hr3_plane_throw(void)
  { hr3_plane_t A;
    r4_throw_dir(&(A.f));
    return A;
  }

void hr3_point_print (FILE *f, char *pre, hr3_point_t *p, char *fmt, char *suf)
  { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, f); }
    fputc('[', f);
    if (fmt == NULL) { fmt = "24.26e"; }
    for (int32_t i = 0; i < NH; i++)
      { if (i != 0) { fputc(' ', f); }
        fprintf(f, fmt, p->c.c[i]);
      }
    fputc(']', f);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, f); }
  }

void hr3_plane_print (FILE *f, char *pre, hr3_plane_t *P, char *fmt, char *suf)
  { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, f); }
    fputc('[', f);
    if (fmt == NULL) { fmt = "24.26e"; }
    for (int32_t i = 0; i < NH; i++)
      { if (i != 0) { fputc(' ', f); }
        fprintf(f, fmt, P->f.c[i]);
      }
    fputc(']', f);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, f); }
  }

void hr3_line_print (FILE *f, char *pre, hr3_line_t *L, char *fmt, char *suf)
  { if ((pre != NULL) && ((*pre) != 0)) { fputs(pre, f); }
    fputc('[', f);
    if (fmt == NULL) { fmt = "24.26e"; }
    for (int32_t i = 0; i < NG; i++)
      { if (i != 0) { fputc(' ', f); }
        fprintf(f, fmt, L->k.c[i]);
      }
    fputc(']', f);
    if ((suf != NULL) && ((*suf) != 0)) { fputs(suf, f); }
  }
  
void hr3_pmap_print (FILE *wr, hr3_pmap_t *M, char *pref, char *suff)
  { 
    hr3_pmap_gen_print(wr, M, "%12.6f", pref, "  ", "  ", "\n", "[ ", " ", " ]", suff);
    fflush(wr);
  }

void hr3_pmap_gen_print
  ( FILE *wr, hr3_pmap_t *M,
    char *fmt, char *pref,                /* Overall prefix. */
    char *rpref, char *rsep, char *rsuff, /* Row prefix, matrix separator, and suffix. */
    char *elp, char *esep, char *erp,     /* Delimiters for each row. */
    char *suff                            /* Overall sufffix. */
  )
  {
    /* Defaults: */
    if (fmt == NULL) { fmt = "%12.6f"; }
    if (elp == NULL) { elp = "[ "; }
    if (esep == NULL) { esep = " "; }
    if (erp == NULL) { erp = " ]"; }
    if (rpref == NULL) { rpref = "  "; }
    if (rsep == NULL) { rsep = " "; }
    if (rsuff == NULL) { rsuff = "\n"; }
    
    if (pref != NULL) { fputs(pref, wr); }
    for (int32_t i = 0; i < NH; i++)
      { fputs(rpref, wr);
        for (int32_t k = 0; k < 2; k++)
          { if (k != 0) { fputs(rsep, wr); }
            fputs(elp, wr);
            for (int32_t j = 0; j < NH; j++)
              { if (j != 0) { fputs(esep, wr); }
                fprintf(wr, fmt, (k == 0 ? M->dir : M->inv).c[i][j]);
              }
            fputs(erp, wr);
          }
        fputs(rsuff, wr);
      }
    if (suff != NULL) { fputs(suff, wr); }
  }
