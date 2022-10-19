/* See hr2.h */
/* Last edited on 2022-10-19 11:42:11 by stolfi */ 

/* Based on HR2.m3 created 1994-05-04 by J. Stolfi. */

#define _GNU_SOURCE
#include <stdint.h>
#include <math.h>
#include <assert.h>

#include <r2.h>
#include <r3.h>
#include <rn.h>
#include <r3x3.h>
#include <r2x2.h>
#include <affirm.h>
#include <sign.h>
#include <sign_get.h>

#include <hr2.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

/* INTERNAL PROTOPTYPES */

void hr2_pmap_normalize(r3x3_t *A);
  /* Divides the matrix {A} by the absolute value of element {A[0][0]}. 
    If that element is zero, the matrix is filled with {NAN}. */

/* IMPLEMENTATIONS */

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

bool_t hr2_pmap_is_identity(hr2_pmap_t *M)
  {
    return r3x3_is_unif_scaling(&(M->dir), M->dir.c[0][0]); 
  }

hr2_point_t hr2_pmap_point(hr2_point_t *p, hr2_pmap_t *M)
  {
    r3_t q; 
    r3x3_map_row(&(p->c), &(M->dir), &q);
    return (hr2_point_t){q};
  }

hr2_line_t hr2_pmap_line(hr2_line_t *L, hr2_pmap_t *M)
  {
    r3_t f;
    r3x3_map_col(&(M->inv), &(L->f), &f);
    return (hr2_line_t){f};
  }

r2_t hr2_pmap_r2_point(r2_t *p, hr2_pmap_t *M)
  {
    r3_t ph = (r3_t){{ 1, p->c[0], p->c[1] }};
    r3_t qh;
    r3x3_map_row(&ph, &(M->dir), &qh);
    double w = qh.c[0];
    double m = fmax(fabs(qh.c[1]), fabs(qh.c[2]));
    if (w <= m*1e-200) 
      { return (r2_t){{ NAN, NAN }}; }
    else
      { return (r2_t){{ qh.c[1]/w, qh.c[2]/w }}; } 
  }

hr2_pmap_t hr2_pmap_identity(void)
  { hr2_pmap_t M;
    r3x3_ident(&(M.dir));
    M.inv = M.dir;
    return M;
  }

hr2_pmap_t hr2_pmap_translation(r2_t *vec)
  { hr2_pmap_t M;
    r3x3_ident(&(M.dir));
    r3x3_ident(&(M.inv));

    M.dir.c[0][1] = +vec->c[0];
    M.dir.c[0][2] = +vec->c[1];

    M.inv.c[0][1] = -vec->c[0];
    M.inv.c[0][2] = -vec->c[1];
    return M;
  }

hr2_pmap_t hr2_pmap_rotation(double ang)
  {
    double c = cos(ang);
    double s = sin(ang);

    hr2_pmap_t M;
    r3x3_ident(&(M.dir));
    r3x3_ident(&(M.inv));

    M.dir.c[1][1] = +c;
    M.dir.c[1][2] = +s;
    M.dir.c[2][1] = -s;
    M.dir.c[2][2] = +c;

    M.inv.c[1][1] = +c;
    M.inv.c[1][2] = -s;
    M.inv.c[2][1] = +s;
    M.inv.c[2][2] = +c;
    return M;
  }

hr2_pmap_t hr2_pmap_scaling(r2_t *scale)
  {
    for (int32_t j = 0; j < 2; j++)
      { demand(scale->c[j] != 0.0, "cannot scale by zero"); }
    
    hr2_pmap_t M;
    r3x3_ident(&(M.dir));
    r3x3_ident(&(M.inv));

    M.dir.c[1][1] = scale->c[0];
    M.dir.c[2][2] = scale->c[1];

    M.inv.c[1][1] = 1/scale->c[0];
    M.inv.c[2][2] = 1/scale->c[1];
    return M;
  }

hr2_pmap_t hr2_pmap_rotation_and_scaling(double ang, double scale)
  {
    double c = cos(ang);
    double s = sin(ang);

    hr2_pmap_t M;
    M.dir.c[0][0] = 1.0;       M.inv.c[0][0] = 1.0;     
    M.dir.c[0][1] = 0.0;       M.inv.c[0][1] = 0.0;     
    M.dir.c[0][2] = 0.0;       M.inv.c[0][2] = 0.0;     

    M.dir.c[1][0] = 0.0;       M.inv.c[1][0] = 0.0;     
    M.dir.c[1][1] = +c*scale;  M.inv.c[1][1] = +c/scale;
    M.dir.c[1][2] = +s*scale;  M.inv.c[1][2] = -s/scale;

    M.dir.c[2][0] = 0.0;       M.inv.c[2][0] = 0.0;     
    M.dir.c[2][1] = -s*scale;  M.inv.c[2][1] = +s/scale;
    M.dir.c[2][2] = +c*scale;  M.inv.c[2][2] = +c/scale;

    return M;
  }

hr2_pmap_t hr2_pmap_compose(hr2_pmap_t *M, hr2_pmap_t *N)
  {
    r3x3_t dir, inv;
    r3x3_mul(&(M->dir), &(N->dir), &dir);
    r3x3_mul(&(N->inv), &(M->inv), &inv);
    return (hr2_pmap_t){dir, inv};
  }

hr2_pmap_t hr2_pmap_inv(hr2_pmap_t *M)
  {
    return (hr2_pmap_t){M->inv, M->dir};
  }

hr2_pmap_t hr2_pmap_inv_comp(hr2_pmap_t *M, hr2_pmap_t *N)
  {
    r3x3_t dir, inv;
    r3x3_mul(&(M->inv), &(N->dir), &dir);
    r3x3_mul(&(N->inv), &(M->dir), &inv);
    return (hr2_pmap_t){dir, inv};
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

hr2_pmap_t hr2_pmap_congruence_from_point_and_dir(r2_t *p, r2_t *u, bool_t flip)
  { 
    /* Normaluze {u} to unit length to get the image of the {(1,0)} vector: */
    r2_t du;  
    double mu = r2_dir(u, &du);
    demand(mu != 0, "invalid direction {u}");

    /* Get the unit vector {dv} that is to be the image of the {(0,1)} vector: */
    r2_t dv;
    if (flip)
      { dv = (r2_t){{ +du.c[1], -du.c[0] }}; }
    else
      { dv = (r2_t){{ -du.c[1], +du.c[0] }}; }
      
    /* Assemble the matrix: */
    hr2_pmap_t M; /* The resulting map. */
    r2_t q = (r2_t){{ -r2_dot(p,&du), -r2_dot(p, &dv) }};
    M.dir.c[0][0] = 1.0;        M.inv.c[0][0] = 1.0;    
    M.dir.c[0][1] = p->c[0];    M.inv.c[0][1] = q.c[0];
    M.dir.c[0][2] = p->c[1];    M.inv.c[0][2] = q.c[1];

    M.dir.c[1][0] = 0.0;        M.inv.c[1][0] = 0.0;    
    M.dir.c[1][1] = +du.c[0];   M.inv.c[1][1] = du.c[0];
    M.dir.c[1][2] = +du.c[1];   M.inv.c[1][2] = dv.c[0];

    M.dir.c[2][0] = 0.0;        M.inv.c[2][0] = 0.0;    
    M.dir.c[2][1] = +dv.c[0];   M.inv.c[2][1] = du.c[1];
    M.dir.c[2][2] = +dv.c[1];   M.inv.c[2][2] = dv.c[1];

    return M;
  }

hr2_pmap_t hr2_pmap_similarity_from_two_points(r2_t *p, r2_t *q, bool_t flip)
  {
    /* Compute the vector {u} that is the image of the {(1,0)} vector: */
    r2_t u; r2_sub(q, p, &u);

    /* Get the squared length {d2} of {u}: */
    double d2 = r2_norm_sqr(&u);
    demand(d2 != 0, "points {p,q} coincide");

    /* Get the vector {v} that is to be the image of the {(0,1)} vector: */
    r2_t v;
    if (flip)
      { v = (r2_t){{ +u.c[1], -u.c[0] }}; }
    else
      { v = (r2_t){{ -u.c[1], +u.c[0] }}; }
      
    /* Assemble the matrix: */
    hr2_pmap_t M; /* The resulting map. */
    r2_t r = (r2_t){{ -r2_dot(p,&u), -r2_dot(p, &v) }};
    M.dir.c[0][0] = 1.0;        M.inv.c[0][0] = d2;    
    M.dir.c[0][1] = p->c[0];    M.inv.c[0][1] = r.c[0];
    M.dir.c[0][2] = p->c[1];    M.inv.c[0][2] = r.c[1];

    M.dir.c[1][0] = 0.0;        M.inv.c[1][0] = 0.0;    
    M.dir.c[1][1] = +u.c[0];    M.inv.c[1][1] = u.c[0];
    M.dir.c[1][2] = +u.c[1];    M.inv.c[1][2] = v.c[0];

    M.dir.c[2][0] = 0.0;        M.inv.c[2][0] = 0.0;    
    M.dir.c[2][1] = +v.c[0];    M.inv.c[2][1] = u.c[1];
    M.dir.c[2][2] = +v.c[1];    M.inv.c[2][2] = v.c[1];

    return M;
  }

hr2_pmap_t hr2_pmap_from_four_points(hr2_point_t *p, hr2_point_t *q, hr2_point_t *r, hr2_point_t *u)
  {
    hr2_pmap_t M; /* The resulting map. */
    
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
    for (int32_t i = 0; i < NH; i++) { w.c[i] = fabs(w.c[i]); }

    /* Ensure that {M.dir} maps the cardinal points to {p,q,r} and some unit point to {u}: */
    r3x3_t *P = &(M.dir);
    P->c[0][0] = w.c[0]*p->c.c[0]; P->c[0][1] = w.c[0]*p->c.c[1]; P->c[0][2] = w.c[0]*p->c.c[2];
    P->c[1][0] = w.c[1]*q->c.c[0]; P->c[1][1] = w.c[1]*q->c.c[1]; P->c[1][2] = w.c[1]*q->c.c[2];
    P->c[2][0] = w.c[2]*r->c.c[0]; P->c[2][1] = w.c[2]*r->c.c[1]; P->c[2][2] = w.c[2]*r->c.c[2];

    /* Compute the inverse map: */
    r3x3_inv(&(M.dir), &(M.inv));

    return M;

  }

double hr2_pmap_mismatch_sqr(hr2_pmap_t *M, int32_t np, r2_t p1[], r2_t p2[])
  {
    double sum2 = 0.0;
    hr2_pmap_t N = hr2_pmap_inv(M);
    for (int32_t k = 0; k < np; k++)
      { r2_t *p1k = &(p1[k]);
        r2_t q1k = hr2_pmap_r2_point(p1k, M);
        r2_t *p2k = &(p2[k]);
        r2_t q2k = hr2_pmap_r2_point(p2k, &N);
        double d2 = r2_dist_sqr(&q1k, &q2k);
        sum2 += d2;
      }
    return sum2/np;
  }

double hr2_pmap_deform_sqr(r2_t ph[], hr2_pmap_t *M)
  {
    int32_t nk = 4;           /* Number of corners of the quadrilateral. */
    r2_t qh[nk];
    for (int32_t k = 0; k < nk; k++)
      { qh[k] = hr2_pmap_r2_point(&(ph[k]), M); }
    
    int32_t nd = nk*(nk-1)/2; /* Number of distances to probe. */
    assert(nd == 6);
    double logr[nd];
    int32_t kd = 0;
    for (int32_t ik = 1; ik < nk; ik++)
      { for (int32_t jk = 0; jk < ik; jk++)
          { double dp2 = r2_dist_sqr(&(ph[ik]), &(ph[jk]));
            double dq2 = r2_dist_sqr(&(qh[ik]), &(qh[jk]));
            assert(kd < nk);
            logr[kd] = log(dq2/dp2)/2;
            kd++;
          }
      }
    assert(kd == nd);
    
    /* Compute the variance of the logs: */
    double sum = 0;
    for (int32_t kd = 0; kd < nd; kd++) { sum += logr[kd]; }
    double avg = sum/nd;
    double sum2 = 0;
    for (int32_t kd = 0; kd < nd; kd++) { double dk = logr[kd] - avg; sum2 += dk*dk; }
    double var = sum2/(nd-1);
    return var;
  }

bool_t hr2_pmap_is_affine(hr2_pmap_t *M)
  { r3x3_t *A = &(M->dir);
    if (A->c[0][0] <= 0) { return FALSE; }
    if ((A->c[1][0] != 0.0) || (A->c[2][0] != 0.0)) { return FALSE; }
    return TRUE;
  }

void hr2_pmap_normalize(r3x3_t *A)
  { 
    double w = fabs(A->c[0][0]);
    if (w != 0)
      { for (int32_t i = 0; i < 3; i++)
          { for (int32_t j = 0; j < 3; j++)
             { A->c[i][j] /= w; }
          }
        assert(fabs(A->c[0][0]) == 1.0);
      }
    else
      { for (int32_t i = 0; i < 3; i++)
          { for (int32_t j = 0; j < 3; j++)
             { A->c[i][j] = NAN; }
          }
      }
  }

void hr2_pmap_print (FILE *wr, hr2_pmap_t *M, char *pref, char *suff)
  { 
    hr2_pmap_gen_print(wr, M, "%12.6f", pref, "[ ", " ", " ]", "  ", "  ", "\n", suff);
    fflush(wr);
  }

void hr2_pmap_gen_print
  ( FILE *wr, hr2_pmap_t *M,
    char *fmt, char *pref,                           /* Overall prefix. */
    char *elp, char *esep, char *erp,     /* Delimiters for each row. */
    char *rpref, char *rsep, char *rsuff, /* Row prefix, matrix separator, and suffix. */
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
                fprintf(wr, fmt, M->dir.c[i][j]);
              }
            fputs(erp, wr);
          }
        fputs(rsuff, wr);
      }
    if (suff != NULL) { fputs(suff, wr); }
  }

hr2_pmap_t hr2_pmap_aff_from_mat_and_disp(r2x2_t *E, r2_t *d)
  {
    hr2_pmap_t M;

    M.dir.c[0][0] = 1.0;
    M.dir.c[0][1] = d->c[0];
    M.dir.c[0][2] = d->c[1];
    
    M.dir.c[1][0] = 0.0;
    M.dir.c[1][1] = E->c[0][0];
    M.dir.c[1][2] = E->c[0][1];
    
    M.dir.c[2][0] = 0.0;
    M.dir.c[2][1] = E->c[1][0];
    M.dir.c[2][2] = E->c[1][1];
    
    r2x2_t F;
    r2x2_inv(E, &(F));
    r2_t psid;
    r2x2_map_row(d, &(F), &psid);
    
    M.inv.c[0][0] = 1.0;
    M.inv.c[0][1] = - psid.c[0];
    M.inv.c[0][2] = - psid.c[1];
    
    M.inv.c[1][0] = 0.0;
    M.inv.c[1][1] = F.c[0][0];
    M.inv.c[1][2] = F.c[0][1];
    
    M.inv.c[2][0] = 0.0;
    M.inv.c[2][1] = F.c[1][0];
    M.inv.c[2][2] = F.c[1][1];
    
    return M;
  }

double hr2_pmap_aff_mismatch_sqr(hr2_pmap_t *M, hr2_pmap_t *N)
  {
    r3x3_t A = M->dir;
    hr2_pmap_normalize(&A);
    r3x3_t B = N->dir;
    hr2_pmap_normalize(&B);
    /* Hope the math is right: */
    r3x3_t E, H;
    r3x3_sub(&A, &B, &E);
    r3x3_mul_tr(&E, &E, &H);
    double h2 = (H.c[1][1] + H.c[2][2])/2;
    double d2 = H.c[0][0];
    return h2 + d2;
  }

hr2_pmap_t hr2_pmap_aff_from_points(r2_t *o, r2_t *p, r2_t *q)
  {
    hr2_pmap_t M; 
    /* The resulting map. */
    
    M.dir.c[0][0] = 1.0;
    M.dir.c[0][1] = o->c[0];
    M.dir.c[0][2] = o->c[1];

    M.dir.c[1][0] = 0.0;
    M.dir.c[1][1] = p->c[0] - o->c[0];
    M.dir.c[1][2] = p->c[1] - o->c[1];

    M.dir.c[2][0] = 0.0;
    M.dir.c[2][1] = q->c[0] - o->c[0];
    M.dir.c[2][2] = q->c[1] - o->c[1];

    r3x3_inv(&(M.dir), &(M.inv));
    
    return M;
  }

hr2_pmap_t hr2_pmap_aff_from_point_pairs(int32_t np, r2_t p1[], r2_t p2[], double w[])
  {
    double debug = FALSE;
    
    if (debug) { fprintf(stderr, "--- computing the affine matrix ---\n"); }
    
    r3x3_t A;
    r3x3_ident(&A);
    if (np > 0)
      { /* Compute the barycenters of {p1} and {p2}: */
        r2_t bar1; r2_barycenter(np, p1, w, &bar1);
        if (debug) { r2_gen_print(stderr, &bar1, "%10.4f", "  bar1 = ( ", " ", " )\n"); }
        demand(r2_is_finite(&bar1), "barycenter of {p1} is undefined");

        r2_t bar2; r2_barycenter(np, p2, w, &bar2);
        if (debug) { r2_gen_print(stderr, &bar2, "%10.4f", "  bar2 = ( ", " ", " )\n"); }
        demand(r2_is_finite(&bar2), "barycenter of {p2} is undefined");

        r2_t d; /* Displacement vector */
        
        if (np == 1)
          { /* Just translate {bar1} to {bar2}: */
            r2_sub(&bar2, &bar1, &d);
          }
        else
          { /* Determine the linear map matrix {L = A.c[1..2][1..2]}: */
            r2x2_t L;
            if (np == 2)
              { /* Compute {L} by composing rotation & scale matrices: */
                r2_t q1; r2_sub(&(p1[0]), &bar1, &q1);
                r2x2_t R1; r2x2_rot_and_scale(&q1, &R1);

                r2_t q2; r2_sub(&(p2[0]), &bar2, &q2);
                r2x2_t R2; r2x2_rot_and_scale(&q2, &R2);

                r2x2_inv(&R1, &R1); r2x2_mul(&R1, &R2, &L);
              }
            else
              { /* Compute {L} by least squares: */

                r2x2_t E; r2x2_zero(&E); /* Moment matrix. */
                r2x2_t P; r2x2_zero(&P); /* Projection matrix. */
                for (int32_t k = 0; k < np; k++)
                  { double wk = (w == NULL ? 1.0 : w[k]);
                    /* Reduce points relative to barycenter: */
                    r2_t q1k, q2k;
                    r2_sub(&(p1[k]), &bar1, &q1k);
                    r2_sub(&(p2[k]), &bar2, &q2k);
                    /* Accumulate moments and projections: */
                    for (int32_t i = 0; i < 2; i ++)
                      { for (int32_t j = 0; j < 2; j++)
                          { E.c[i][j] += wk*q1k.c[i]*q1k.c[j];
                            P.c[i][j] += wk*q1k.c[i]*q2k.c[j];
                          }
                      }
                  }
                r2x2_t Z; r2x2_inv(&E, &Z);
                r2x2_mul(&Z, &P, &L);
              }
              
            /* Store linear part into {A}: */
            A.c[1][1] = L.c[0][0];
            A.c[1][2] = L.c[0][1];
            A.c[2][1] = L.c[1][0];
            A.c[2][2] = L.c[1][1];
            
            /* Compute the displacement taking {L} into account: */
            r2_t v1; r2x2_map_row(&bar1, &L, &v1);
            r2_sub(&bar2, &v1, &d);
          }

        /* Store the displacement vector {d}: */
        A.c[0][1] = d.c[0];
        A.c[0][2] = d.c[1];
      }

    if (debug) 
      { fprintf(stderr, "  matrix:\n");
        r3x3_gen_print(stderr, &(A), "%13.6e", "", "\n", "\n", "    [ ", " ", " ]");
      }
      
    hr2_pmap_t M;
    M.dir = A;
    r3x3_inv(&A, &(M.inv));
    return M;
  }
