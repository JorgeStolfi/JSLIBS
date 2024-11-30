/* See hr2_pmap.h */
/* Last edited on 2024-11-24 08:44:21 by stolfi */ 

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


#include <hr2_pmap.h>

#define NH 3
  /* Number of homogeneous coordinates in a point. */

#define NC 2
  /* Number of Cartesian coordinates in a point. */

/* INTERNAL FUNCTIONS */
    
void hr2_pmap_fill_rot_scale_matrix(r3x3_t *A, double ux, double uy);
  /* Fills {A} with the direct matrix of a linear projective
    map that performs a rotation and uniform scaling about the 
    origin that sends the Cartesian vector {(1,0)} to the vector
    {(ux,uy)} (which had better be nonzero). Specifically, sets the rows of {A}
    to {1,0,0}, {0,+ux,+uy}, and {0,-uy,+ux}. */

/* IMPLEMENTATIONS */

hr2_point_t hr2_pmap_point(hr2_point_t *p, hr2_pmap_t *M)
  {
    r3_t q; 
    r3x3_map_row(&(p->c), &(M->dir), &q);
    return (hr2_point_t){q};
  }

hr2_point_t hr2_pmap_inv_point(hr2_point_t *p, hr2_pmap_t *M)
  {
    r3_t q; 
    r3x3_map_row(&(p->c), &(M->inv), &q);
    return (hr2_point_t){q};
  }

hr2_line_t hr2_pmap_line(hr2_line_t *L, hr2_pmap_t *M)
  {
    r3_t f;
    r3x3_map_col(&(M->inv), &(L->f), &f);
    return (hr2_line_t){f};
  }

hr2_line_t hr2_pmap_inv_line(hr2_line_t *L, hr2_pmap_t *M)
  {
    r3_t f;
    r3x3_map_col(&(M->dir), &(L->f), &f);
    return (hr2_line_t){f};
  }

r2_t hr2_pmap_r2_point(r2_t *p, hr2_pmap_t *M)
  {
    hr2_point_t ph = (hr2_point_t){{{ 1.0, p->c[0], p->c[1] }}};
    hr2_point_t qh;
    r3x3_map_row(&(ph.c), &(M->dir), &(qh.c));
    r2_t qc = r2_from_hr2_nan(&qh);
    return qc;
  }

r2_t hr2_pmap_inv_r2_point(r2_t *p, hr2_pmap_t *M)
  {
    hr2_point_t ph = (hr2_point_t){{{ 1, p->c[0], p->c[1] }}};
    hr2_point_t qh;
    r3x3_map_row(&(ph.c), &(M->inv), &(qh.c));
    r2_t qc = r2_from_hr2_nan(&qh);
    return qc;
  }

hr2_pmap_t hr2_pmap_identity(void)
  { hr2_pmap_t M;
    r3x3_ident(&(M.dir));
    M.inv = M.dir;
    return M;
  }

hr2_pmap_t hr2_pmap_mirror(sign_t sgnx, sign_t sgny)
  { demand(sgnx != 0, "{sgnx} must be {-1} or {+1}");
    demand(sgny != 0, "{sgny} must be {-1} or {+1}");
    hr2_pmap_t M;
    r3x3_ident(&(M.dir));
    M.dir.c[1][1] = sgnx;
    M.dir.c[2][2] = sgny;
    M.inv = M.dir;
    return M;
  }

hr2_pmap_t hr2_pmap_xy_swap(void)
  { hr2_pmap_t M;
    r3x3_ident(&(M.dir));
    M.dir.c[1][1] = 0; M.dir.c[1][2] = 1;
    M.dir.c[2][1] = 1; M.dir.c[2][2] = 0;
    M.inv = M.dir;
    return M;
  }

hr2_pmap_t hr2_pmap_rotation(double ang)
  {
    double c = cos(ang);
    double s = sin(ang);

    hr2_pmap_t M;
    hr2_pmap_fill_rot_scale_matrix(&(M.dir), +c, +s);
    hr2_pmap_fill_rot_scale_matrix(&(M.dir), +c, -s);

    return M;
  }

hr2_pmap_t hr2_pmap_scaling(r2_t *scale)
  {
    for (uint32_t j = 0;  j < 2; j++)
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

hr2_pmap_t hr2_pmap_inv_compose(hr2_pmap_t *M, hr2_pmap_t *N)
  {
    r3x3_t dir, inv;
    r3x3_mul(&(M->inv), &(N->dir), &dir);
    r3x3_mul(&(N->inv), &(M->dir), &inv);
    return (hr2_pmap_t){dir, inv};
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
      /* Map {u} by the inverse of {Q}: */
      r3x3_inv(&Q, &Q);
      r3x3_map_row(&(u->c), &Q, &w);
    }
    
    /* Make the weights positive, so that {p,q,r} are strictly honored: */
    for (uint32_t i = 0;  i < NH; i++) { w.c[i] = fabs(w.c[i]); }

    /* Ensure that {M.dir} maps the cardinal points to {p,q,r} and some unit point to {u}: */
    r3x3_t *P = &(M.dir);
    P->c[0][0] = w.c[0]*p->c.c[0]; P->c[0][1] = w.c[0]*p->c.c[1]; P->c[0][2] = w.c[0]*p->c.c[2];
    P->c[1][0] = w.c[1]*q->c.c[0]; P->c[1][1] = w.c[1]*q->c.c[1]; P->c[1][2] = w.c[1]*q->c.c[2];
    P->c[2][0] = w.c[2]*r->c.c[0]; P->c[2][1] = w.c[2]*r->c.c[1]; P->c[2][2] = w.c[2]*r->c.c[2];

    /* Compute the inverse map: */
    r3x3_inv(&(M.dir), &(M.inv));

    return M;
  }

hr2_pmap_t hr2_pmap_from_four_r2_points(r2_t *p, r2_t *q, r2_t *r, r2_t *u)
  {
    hr2_point_t hp = hr2_from_r2(p);
    hr2_point_t hq = hr2_from_r2(q);
    hr2_point_t hr = hr2_from_r2(r);
    hr2_point_t hu = hr2_from_r2(u);
    hr2_pmap_t M = hr2_pmap_from_four_points(&hp, &hq, &hr, &hu);
    for (uint32_t j = 0;  j < NH; j++)
      { M.dir.c[1][j] -= M.dir.c[0][j];  
        M.dir.c[2][j] -= M.dir.c[0][j];
        M.inv.c[j][0] += M.inv.c[j][1] + M.inv.c[j][2];
      }
    return M;
  }

hr2_pmap_t hr2_pmap_r2_from_sign_class(uint32_t class)
  {
    hr2_pmap_t M; 
    
    M.dir.c[0][0] = 1;
    M.dir.c[0][1] = 0;
    M.dir.c[0][2] = 0;

    int32_t b = ((class & 2) == 0 ? +1 : -1);
    M.dir.c[1][0] = b - 1;
    M.dir.c[1][1] = b;
    M.dir.c[1][2] = 0;

    int32_t c = ((class & 1) == 0 ? +1 : -1);
    M.dir.c[2][0] = c - 1;
    M.dir.c[2][1] = 0;
    M.dir.c[2][2] = c;

    r3x3_inv(&(M.dir), &(M.inv));
    
    return M;
  }


hr2_pmap_t hr2_pmap_from_parms
  ( r2_t *persp,
    double shear,
    double skew,
    double scale,
    double ang,
    r2_t *disp
  )
  { 
    hr2_pmap_t M;
    r3x3_t *A = &(M.dir); 
    r3x3_ident(A);
    if ((persp != NULL) && ((persp->c[0] != 0) || (persp->c[1] != 0)))
      { /* Create a perspective map: */
        r3x3_t P;
        r3x3_ident(&P);
        P.c[1][0] = persp->c[0];
        P.c[2][0] = persp->c[1];
        /* Compose it before {A}: */
        (*A) = P;
      }
 
    if ((shear != 0) || (skew != 0))
      { /* Make matrix {S} of a shear linear map with given shear and skew amounts: */
        demand(isfinite(shear) && (fabs(shear) < 1.0), "invalid {shear}");
        demand(isfinite(skew) && (skew >= 1.0e-100) && (skew <= 1.0e+100), "invalid {skew}");
        double rsk = sqrt(skew);
        r3x3_t S;
        for (uint32_t i = 1;  i < NH; i++)
          { double mag = (i == 1 ? rsk : 1/rsk);
            for (uint32_t j = 1;  j < NH; j++)
              { S.c[i][j] = (i == j ? mag : shear); }
          }
        /* Compose it after {A}: */
        r3x3_mul(A, &S, A);
      }
    
    if ((ang != 0) || (scale != 1.0))
      { /* Generate matrix of rotation linear map with given angle and scale: */
        demand(isfinite(scale) && (scale >= 1.0e-100) && (scale <= 1.0e+100), "invalid {scale}");
        r3x3_t R;
        double c = cos(ang), s = sin(ang);
        hr2_pmap_fill_rot_scale_matrix(&R, c*scale, s*scale);
        /* Compose it after {A}: */
        r3x3_mul(A, &R, A);
      }
    
    if ((disp != NULL) && ((disp->c[0] != 0) || (disp->c[1] != 0)))
      { /* Modify {A} with translation by {disp}: */
        demand(isfinite(disp->c[0]) && isfinite(disp->c[1]), "invalid {disp}");
        A->c[0][1] = disp->c[0];
        A->c[0][2] = disp->c[1];
      }
   
    /* Build the projective map: */
    r3x3_inv(&(M.dir), &(M.inv));
    M.inv.c[1][0] = 0; /* Paranoia. */
    M.inv.c[2][0] = 0; /* Paranoia. */
    /* Ensure {M.inv[0,0] = 1}: */
    double B00 = M.inv.c[0][0];
    if (B00 != 1.0)
      { /* Scale matrix homogeneously by {1/B00}: */
        assert(isfinite(B00) && (B00 > 1.0e-200) && (B00 < 2.0e200));
        for (uint32_t i = 0;  i < NH; i++)
          { for (uint32_t j = 1;  j < NH; j++)
              { M.inv.c[i][j] /= B00; }
          }
        M.inv.c[0][0] = 1.0;
      }
    return M;
  }

void hr2_pmap_throw(hr2_pmap_t *M)
  {
    double det;
    do
      { for (uint32_t i = 0;  i < NH; i++)
          { r3_t a;
            r3_throw_cube(&a);
            for (uint32_t j = 0;  j < NH; j++) 
              { M->dir.c[i][j] = a.c[j]; }
          }
        det = r3x3_det(&(M->dir));
      }
    while ((fabs(det) < 1.0e-6) || (fabs(det) > 1.0e6));
    r3x3_inv(&(M->dir), &(M->inv));
  }

double hr2_pmap_diff_sqr(hr2_pmap_t *M, hr2_pmap_t *N)
  { double sum_d2 = 0;
    for (uint32_t sense = 0;  sense < 2; sense++)
      { r3x3_t *A = (sense == 0 ? &(M->dir) : &(M->inv));
        double Am = r3x3_norm(A) + 1.0e-200;
        r3x3_t *B = (sense == 0 ? &(N->dir) : &(N->inv));
        double Bm = r3x3_norm(B) + 1.0e-200;
        for (uint32_t i = 0;  i < NH; i++)
          { for (uint32_t j = 0;  j < NH; j++)
             { double Aij = A->c[i][j]/Am;
               double Bij = B->c[i][j]/Bm;
               double dij = Aij - Bij;
               sum_d2 += dij*dij;
             }
          }
      }
    return sum_d2;
  }

double hr2_pmap_mismatch_sqr(hr2_pmap_t *M, uint32_t np, r2_t p1[], r2_t p2[], double w[])
  {
    double sum_wD2 = 0.0;
    double sum_w = 1.0e-200; /* In case there are no points with positive weight. */
    for (uint32_t k = 0;  k < np; k++)
      { r2_t *p1k = &(p1[k]);
        r2_t *p2k = &(p2[k]);
        double wk = (w == NULL ? 1.0 : w[k]);
        r2_t q1k = hr2_pmap_r2_point(p1k, M);
        r2_t q2k = hr2_pmap_inv_r2_point(p2k, M);
        double D2k_dir = r2_dist_sqr(&q1k, p2k);
        double D2k_inv = r2_dist_sqr(p1k, &q2k);
            
        sum_wD2 += wk*(D2k_dir + D2k_inv);
        sum_w += wk;
      }
    return 0.5*sum_wD2/sum_w;
  }

double hr2_pmap_max_mismatch(hr2_pmap_t *M, uint32_t np, r2_t p1[], r2_t p2[])
  {
    /* Compute max distance between paired points mapped  by the initial guess: */
    double maxDist = 0.0;
    for (uint32_t kp = 0;  kp < np; kp++)
      { r2_t *p1k = &(p1[kp]);
        r2_t *p2k = &(p2[kp]);
        r2_t q1 = hr2_pmap_r2_point(p1k, M);
        double dk12 = r2_dist(&q1, p2k);
        if (dk12 > maxDist) { maxDist = dk12; }
        r2_t q2 = hr2_pmap_inv_r2_point(p2k, M);
        double dk21 = r2_dist(&q2, p1k);
        if (dk21 > maxDist) { maxDist = dk21; }
      }
    return maxDist;
  }

void hr2_pmap_show_point_mismatch(hr2_pmap_t *M, uint32_t np, r2_t p1[], r2_t p2[], double w[])
  { 
    for (uint32_t k = 0;  k < np; k++)
      { r2_t *p1k = &(p1[k]);
        r2_t *p2k = &(p2[k]);
        r2_t q1k = hr2_pmap_r2_point(p1k, M);
        r2_t q2k = hr2_pmap_inv_r2_point(p2k, M);
        double Dk_dir = r2_dist(&q1k, p2k);
        double Dk_inv = r2_dist(p1k, &q2k);
        
        fprintf(stderr, "        %3d", k);
        r2_gen_print(stderr, p1k, "%+12.8f", " p1 = ( ", " ", " )");
        r2_gen_print(stderr, &(q1k), "%+12.8f", " --> ( ", " ", " )");
        fprintf(stderr, " d = %12.8f\n", Dk_dir);
        fprintf(stderr, "           ");
        r2_gen_print(stderr, p2k, "%+12.8f", " p2 = ( ", " ", " )");
        r2_gen_print(stderr, &(q2k), "%+12.8f", " --> ( ", " ", " )");
        fprintf(stderr, " d = %12.8f\n", Dk_inv);
      }
  }

double hr2_pmap_deform_sqr(r2_t ph[], hr2_pmap_t *M)
  {
    uint32_t nk = (1 << NC); /* Number of corners of the quadrilateral. */
    r2_t qh[nk];
    for (uint32_t k = 0;  k < nk; k++)
      { qh[k] = hr2_pmap_r2_point(&(ph[k]), M); }
    
    uint32_t nd = 4 + 2; /* Number of distances to probe. */
    double logr[nd];
    uint32_t kd = 0;
    for (uint32_t ik = 1;  ik < nk; ik++)
      { for (uint32_t jk = 0;  jk < ik; jk++)
          { double dp2 = r2_dist_sqr(&(ph[ik]), &(ph[jk]));
            double dq2 = r2_dist_sqr(&(qh[ik]), &(qh[jk]));
            assert(kd < nd);
            logr[kd] = log(dq2/dp2)/2;
            kd++;
          }
      }
    assert(kd == nd);
    
    /* Compute the variance of the logs: */
    double sum = 0;
    for (uint32_t kd = 0;  kd < nd; kd++) { sum += logr[kd]; }
    double avg = sum/nd;
    double sum2 = 0;
    for (uint32_t kd = 0;  kd < nd; kd++) { double dk = logr[kd] - avg; sum2 += dk*dk; }
    double var = sum2/(nd-1);
    return var;
  }

bool_t hr2_pmap_is_valid(hr2_pmap_t *M, double tol)
  { bool_t debug = FALSE; /* Print diagnostics if NOT valid. */
    if (debug) { fprintf(stderr, "  > --- %s ---\n", __FUNCTION__); }
    
    /* Grab {M.dir} and {M.inv}, check nonzero norm, and normalize: */
    bool_t ok = TRUE;
    r3x3_t A[2]; /* Matrices of {M}, eventually normalized. */
    double det[2]; /* Their determinants. */
    char *which[2] = { "dir", "inv" };
    for (uint32_t dir = 0;  dir <= 1; dir++)
      { A[dir] = (dir == 0 ? M->dir : M->inv);
        double enorm = r3x3_norm(&(A[dir]))/3; /* RMS elem value of {A[dir]}. */
        if ((! isfinite(enorm)) || (enorm < 1.0e-100)) 
          { if (debug) { fprintf(stderr, "    norm(M.%s)/3 = %24.16e, bad or too small\n", which[dir], enorm); }
            ok = FALSE;
          }
        else
          { /* Normalize {A[dir]} to unit RMS elem value: */
            r3x3_scale(1.0/enorm, &(A[dir]), &(A[dir]));
          }
      }
    if (ok)
      { /* Both matrices are sufficiently far from zero. Check that both have
          definitely nonzero determinant.
      
          Tests show that if a singular {3×3} matrix with unit element
          RMS is perturbed by adding to each element a random number in
          the range {[-eps _ +eps]}, for small {eps}, the determinant of
          the perturbed matrix will have zero mean and deviation
          {~2*eps}. Thus if {det(A)} is on the order of
          {2*tol}, there is a significant chance that {A} is 
          a singular matrix with each element perturbed by {~tol}. */
        for (uint32_t dir = 0;  dir <= 1; dir++)
          { det[dir] = r3x3_det(&(A[dir]));
            if ((! isfinite(det[dir])) || (fabs(det[dir]) < 6*tol)) 
              { if (debug)
                  { fprintf(stderr, "    det(nrmz(M.%s)) = %24.16e", which[dir], det[dir]);
                    fprintf(stderr, " bad or too small  (min = %24.16e)\n", 6*tol); 
                  }
                ok =  FALSE;
              }
          }
      }
      
    if (ok)
      { /* Both matrices are sufficiently far form zero and have
          sufficiently nonzero det. Check that {M.dir} and {M.inv} are
          inverses apart from scaling factor and roundoff: */
        r3x3_t P; r3x3_adj(&(A[0]), &P); 
        double pnorm = r3x3_norm(&P)/3; /* RMS value of P elements. */
        assert(isfinite(pnorm) && (pnorm >= 1.0e-200));
        for (uint32_t i = 0;  (i < NH) && ok; i++)
          { for (uint32_t j = 0;  (j < NH) && ok; j++)
             { double Pij = P.c[i][j]/pnorm * (det[0] < 0 ? -1 : +1);
               double Bij = A[1].c[i][j];
               if ((! isfinite(Pij)) || (fabs(Pij-Bij) > 10*tol))
                 { if (debug)
                     { fprintf(stderr, "    (nrmz(inv(nrmz(M.dir)))[%d,%d] = %24.16e", i, j, Pij);
                       fprintf(stderr, "  (nrmz(M.inv))[%d,%d] = %24.16e", i, j, Bij);
                       fprintf(stderr, " too different\n"); 
                     }
                   ok =  FALSE;
                 }
             }
          }
      }
    if (debug && (! ok))
      { hr2_pmap_gen_print
          ( stderr, M, "%+24.16e", 
            "    input M = \n",
            "      ", "  ", "\n", 
            NULL,NULL,NULL,
            "\n"
          ); 
        hr2_pmap_t N;
        N.dir = A[0]; N.inv = A[1];
        hr2_pmap_gen_print
          ( stderr, &N, "%+24.16e", 
            "    normalized M = \n",
            "      ", "  ", "\n", 
            NULL,NULL,NULL,
            "\n"
          );
        fprintf(stderr, "    map is NOT valid\n"); 
      }
    if (debug && ok) { fprintf(stderr, "    map seems valid\n"); }
    if (debug) { fprintf(stderr, "  < --- %s ---\n", __FUNCTION__); }
    return ok;
  }

sign_t hr2_pmap_sign(hr2_pmap_t *M)
  { double det = r3x3_det(&(M->dir));
    if (det < 0)
      { return -1; }
    else if (det > 0)
      { return +1; }
    else
      { return 0; }
  }

void hr2_pmap_print (FILE *wr, hr2_pmap_t *M, char *pref, char *suff)
  { 
    hr2_pmap_gen_print(wr, M, "%12.6f", pref, "  ", "  ", "\n", "[ ", " ", " ]", suff);
    fflush(wr);
  }

void hr2_pmap_gen_print
  ( FILE *wr, hr2_pmap_t *M,
    char *fmt, 
    char *pref,                           /* Overall prefix. */
    char *rpref, char *rsep, char *rsuff, /* Row prefix, matrix separator, and suffix. */
    char *elp, char *esep, char *erp,     /* Delimiters for each matrix in each row. */
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
    for (uint32_t i = 0;  i < NH; i++)
      { fputs(rpref, wr);
        for (uint32_t k = 0;  k < 2; k++)
          { if (k != 0) { fputs(rsep, wr); }
            fputs(elp, wr);
            for (uint32_t j = 0;  j < NH; j++)
              { if (j != 0) { fputs(esep, wr); }
                fprintf(wr, fmt, (k == 0 ? M->dir : M->inv).c[i][j]);
              }
            fputs(erp, wr);
          }
        fputs(rsuff, wr);
      }
    if (suff != NULL) { fputs(suff, wr); }
  }

bool_t hr2_pmap_type_from_string(char *tname, hr2_pmap_type_t *type_P)
  { if ((strcmp(tname, "IDENTITY") == 0) || (strcmp(tname, "identity") == 0))
      { (*type_P) = hr2_pmap_type_IDENTITY; }
    else if ((strcmp(tname, "IDENTITY") == 0) || (strcmp(tname, "identity") == 0))
      { (*type_P) = hr2_pmap_type_IDENTITY; }
    else if ((strcmp(tname, "TRANSLATION") == 0) || (strcmp(tname, "translation") == 0))
      { (*type_P) = hr2_pmap_type_TRANSLATION; }
    else if ((strcmp(tname, "CONGRUENCE") == 0) || (strcmp(tname, "congruence") == 0))
      { (*type_P) = hr2_pmap_type_CONGRUENCE; }
    else if ((strcmp(tname, "SIMILARITY") == 0) || (strcmp(tname, "similarity") == 0))
      { (*type_P) = hr2_pmap_type_SIMILARITY; }
    else if ((strcmp(tname, "AFFINE") == 0) || (strcmp(tname, "affine") == 0))
      { (*type_P) = hr2_pmap_type_AFFINE; }
    else if ((strcmp(tname, "GENERIC") == 0) || (strcmp(tname, "projective") == 0))
      { (*type_P) = hr2_pmap_type_GENERIC; }
    else if ((strcmp(tname, "NONE") == 0) || (strcmp(tname, "none") == 0))
      { (*type_P) = hr2_pmap_type_NONE; }
    else 
      { return FALSE; }
    return TRUE;
  }

char *hr2_pmap_type_to_string(hr2_pmap_type_t type)
  { 
    switch(type)
      { case hr2_pmap_type_IDENTITY:
          return "IDENTITY";
        case hr2_pmap_type_TRANSLATION:
          return "TRANSLATION";
        case hr2_pmap_type_CONGRUENCE:
          return "CONGRUENCE";
        case hr2_pmap_type_SIMILARITY:
          return "SIMILARITY";
        case hr2_pmap_type_AFFINE:
          return "AFFINE";
        case hr2_pmap_type_NONE:
          return "NONE";
        case hr2_pmap_type_GENERIC:
          return "GENERIC";
      }
    assert(FALSE);
  }

void hr2_pmap_invert_sign(hr2_pmap_t *M)
  {
    r3x3_t *Md = &(M->dir), *Mi = &(M->inv);
    for (uint32_t j = 0;  j < NH; j++)
      { double td = Md->c[1][j]; Md->c[1][j] = Md->c[2][j]; Md->c[2][j] = td;
        double ti = Mi->c[j][1]; Mi->c[j][1] = Mi->c[j][2]; Mi->c[j][2] = ti;
      }
  }

void hr2_pmap_set_sign(hr2_pmap_t *M, sign_t sgn)
  { demand((sgn == -1) || (sgn == +1), "invalid {sgn}"); 
    double det = r3x3_det(&(M->dir));
    demand(fabs(det) > 1.0e-250, "invalid map - determinant is zero");
    if (sgn*det < 0) { hr2_pmap_invert_sign(M); }
  }

bool_t hr2_pmap_is_identity(hr2_pmap_t *M, double tol)
  { 
    bool_t debug = FALSE;
    if (debug) { fprintf(stderr, "      > --- %s ---\n", __FUNCTION__); }

    bool_t ok = TRUE;
    for (uint32_t dir = 0;  (dir <= 1) && ok; dir++)
      { r3x3_t *A = (dir == 0 ? &(M->dir) : &(M->inv));
        bool_t Aok = TRUE;
        double A00 = A->c[0][0];
        if (A00 < 1.0e-200) { Aok = FALSE; }
        
        if (Aok)
          { double Atol = A00*tol;

            if (fabs(A->c[0][1]) > Atol) { Aok = FALSE; }
            if (fabs(A->c[0][2]) > Atol) { Aok = FALSE; }

            if (fabs(A->c[1][0]) > Atol) { Aok = FALSE; }
            if (fabs(A->c[1][1] - A00) > Atol) { Aok = FALSE; }
            if (fabs(A->c[1][2]) > Atol) { Aok = FALSE; }

            if (fabs(A->c[2][0]) > Atol) { Aok = FALSE; }
            if (fabs(A->c[2][1]) > Atol) { Aok = FALSE; }
            if (fabs(A->c[2][2] - A00) > Atol) { Aok = FALSE; }
          }
        ok = (ok && Aok);
      }
    if (debug) 
      { hr2_pmap_gen_print
          ( stderr, M, "%+24.16e",
            "        M = \n",
            "        ", "  ", "\n", 
            NULL,NULL,NULL,
            "\n"
          );
        fprintf(stderr, "        result: map is%s the identity\n", (ok ? "" : " NOT")); 
      }
    if (debug) { fprintf(stderr, "      < --- %s ---\n", __FUNCTION__); }
    return ok;
  }

bool_t hr2_pmap_is_translation(hr2_pmap_t *M, double tol)
  {
    bool_t debug = FALSE; /* Prints diagnostics anyway. */
    if (debug) { fprintf(stderr, "      > --- %s ---\n", __FUNCTION__); }
    bool_t ok = TRUE;
    bool_t is_affine = hr2_pmap_is_affine(M, tol);
    if (debug) { fprintf(stderr, "        {hr2_pmap_is_affine(M,tol)} = %c\n", "FT"[is_affine]); }
    ok = (ok && is_affine);
    
    auto bool_t check_matrix(r3x3_t *A);
      /* Checks whether the matrix {A} (either {M.dir} or {M.inv} is 
        the {3x3} matrix of a translation {hr2_pmap_t},
        within the specified tolerance. */
    
    for (uint32_t dir = 0;  (dir <= 1) && ok; dir++)
      { r3x3_t *A = (dir == 0 ? &(M->dir) : &(M->inv));
        if (debug) { fprintf(stderr, "          checking {M.%s}\n", (dir == 0 ? "dir" : "inv")); }
        ok = check_matrix(A);
      }
    if (debug) { fprintf(stderr, "        result: map is%s a translation\n", (ok ? "" : " NOT")); }
    if (debug) { fprintf(stderr, "      < --- %s ---\n", __FUNCTION__); }
    return ok;

    bool_t check_matrix(r3x3_t *A)
      { double A00 = A->c[0][0];
        if (A00 < 1.0e-200) { return FALSE; }
        double Atol = A00*tol;
        if (debug) { fprintf(stderr, "        Atol = %24.26e\n", Atol); }

        for (uint32_t i = 1;  i < NH; i++)
          { for (uint32_t j = 1;  j < NH; j++)
              { double Aij = A->c[i][j];
                double Eij = (i == j ? A00 : 0.0);
                if (fabs(Aij - Eij) > Atol) 
                  { if (debug) 
                      { fprintf(stderr, "          A%d%d = %24.16e =", i, j, Aij);
                        if (Eij != 0) { fprintf(stderr, " %24.16e +", Eij); }
                        fprintf(stderr, " Atol * %24.16e\n", (Aij - Eij)/Atol);
                      }
                    return FALSE; 
                  }
              }
          }
        return TRUE;
      }
   }
 
bool_t hr2_pmap_is_congruence(hr2_pmap_t *M, double tol)
  {
    return hr2_pmap_is_similarity(M, 1.0, tol);
  }
   
bool_t hr2_pmap_is_similarity(hr2_pmap_t *M, double scale, double tol)
  {
    bool_t debug = TRUE; /* Prints diagnostics anyway. */
    if (debug) { fprintf(stderr, "      > --- %s ---\n", __FUNCTION__); }
    if (debug) { fprintf(stderr, "        tol = %24.16e\n", tol); }
    
    demand(isnan(scale) || (scale > 1.0e-100), "invalid {scale}");
    
    if (debug) { hr2_pmap_gen_print(stderr, M, "%+24.16e", "        M =\n",  "        ",NULL,NULL, NULL,NULL,NULL, "\n"); }
    bool_t ok = TRUE;
    bool_t is_affine = hr2_pmap_is_affine(M, tol);
    if (debug) { fprintf(stderr, "        {hr2_pmap_is_affine(M,tol)} = %c\n", "FT"[is_affine]); }
    ok = (ok && is_affine);
    
    auto bool_t check_matrix(r3x3_t *A, double sc);
      /* Checks whether the matrix {A} (either {M.dir} or {M.inv} is the
        {3x3} matrix of a similarity {hr2_pmap_t}, within the specified
        tolerance. If {sc} is not {NAN}, requires the scaling factor to
        be {sc}. */
    
    for (uint32_t dir = 0;  (dir <= 1) && ok; dir++)
      { if (debug) { fprintf(stderr, "          checking {M.%s}\n", (dir == 0 ? "dir" : "inv")); }
        r3x3_t *A = (dir == 0 ? &(M->dir) : &(M->inv));
        double sc = (isnan(scale) ? NAN : (dir == 0 ? scale : 1.0/scale));
            
        ok = check_matrix(A, sc);
      }

    if (debug) 
      { fprintf(stderr, "        result: map is%s a similarity", (ok ? "" : " NOT"));
        if (! isnan(scale)) { fprintf(stderr, " with scale %24.16e", scale); }
        fprintf(stderr, "\n");
      }
    if (debug) { fprintf(stderr, "      < --- %s ---\n", __FUNCTION__); }
    return ok;
      
    bool_t check_matrix(r3x3_t *A, double sc)
      { double A00 = A->c[0][0];
        if (debug) { fprintf(stderr, "          A00 = %24.16e\n", A00); }

        double ulen = sqrt(A->c[1][1]*A->c[1][1] + A->c[1][2]*A->c[1][2])/A00;
        double vlen = sqrt(A->c[2][1]*A->c[2][1] + A->c[2][2]*A->c[2][2])/A00;;
        double uvcos = (A->c[1][1]*A->c[2][1] + A->c[1][2]*A->c[2][2])/(A00*A00)/(ulen*vlen);
        if (debug) { fprintf(stderr, "          |u| = %24.16e  |v| = %24.16e  uvcos = %24.16e\n", ulen, vlen, uvcos); } 

        /* The {u,v} vectors must be orthogonal: */
        if (fabs(uvcos) > 2*tol) { return FALSE; }

        double uvdif = ulen - vlen;
        if (debug) { fprintf(stderr, "          |u| - |v| = %24.16e = %24.16e * tol\n", uvdif, uvdif/tol); }

        /* Their length must be equal: */
        if (fabs(uvdif) > tol) { return FALSE; }

        if (! isnan(sc))
          { if (debug) 
              { fprintf(stderr, "          sc = %24.16e\n", sc); 
                fprintf(stderr, "          |u|/sc = %24.16e = %24.16e + 1\n", ulen/sc, ulen/sc-1); 
                fprintf(stderr, "          |v|/sc = %24.16e = %24.16e + 1\n", vlen/sc, vlen/sc-1); 
              } 

            double det = (A->c[1][1]*A->c[2][2] - A->c[1][2]*A->c[2][1])/(A00*A00);
            /* The determinant should be the same as {|u|*|v|}, so: */
            assert(fabs(sqrt(fabs(det)/(ulen*vlen)) - 1) < 10*tol);
            double erel = sqrt(fabs(det))/sc;
            if (debug) 
              { fprintf(stderr, "          det = %24.16e = %24.16e * sc^2\n",  det, det/(sc*sc)); 
                fprintf(stderr, "          sqrt(|det|)/sc  = %24.16e = %24.16e + 1\n", erel, erel-1); 
              } 
            if (fabs(erel - 1) > 4*tol) { return FALSE; }
          }
        return TRUE;
      }
  }

bool_t hr2_pmap_is_affine(hr2_pmap_t *M, double tol)
  {
    bool_t debug = TRUE; /* Prints diagnostics anyway. */
    if (debug) { fprintf(stderr, "      > --- %s ---\n", __FUNCTION__); }
    if (debug) { fprintf(stderr, "        tol = %24.16e\n", tol); }
    
    bool_t ok = TRUE;

    bool_t is_valid = hr2_pmap_is_valid(M, tol);
    if (debug) { fprintf(stderr, "        {hr2_pmap_is_valid(M)} = %c\n", "FT"[is_valid]); }
    ok = (ok && is_valid);
    
    auto bool_t check_matrix(r3x3_t *A);
      /* Checks whether the matrix {A} (either {M.dir} or {M.inv} is the
        {3x3} matrix of an affine {hr2_pmap_t} of definitely nonzero 
        determinant with sign {sgn}, within the specified tolerance. */

    for (uint32_t dir = 0;  (dir <= 1) && ok; dir++)
      { if (debug) { fprintf(stderr, "          checking {A = M.%s}\n", (dir == 0 ? "dir" : "inv")); }
        r3x3_t *A = (dir == 0 ? &(M->dir) : &(M->inv));
        ok = check_matrix(A);
        if (debug) { fprintf(stderr, "          %s.\n", (ok ? "passed" : "failed")); }
        if (debug) { fprintf(stderr, "\n"); }
      }
    if (debug) { fprintf(stderr, "        result: map is%s affine\n", (ok ? "" : " NOT")); }
    if (debug) { fprintf(stderr, "      < --- %s ---\n", __FUNCTION__); }
    return ok;

    bool_t check_matrix(r3x3_t *A)
      { double A00 = A->c[0][0];
        double Atol = A00*tol;
        if (debug) { fprintf(stderr, "          A00 =    %24.16e        Atol = A00 * tol = %24.16e\n", A00, Atol); }
        /* A00 must be definitely positive: */
        if (A00 < 1.0e-200) { return FALSE; }
        
        /* Other column 0 elements should be essentially zero: */
        for (uint32_t i = 1;  i < NH; i++)
          { double Ai0 = A->c[i][0]/A00;
            if (debug) { fprintf(stderr, "          A%d0 =    %24.16e * A00 = %24.16e * Atol\n", i, Ai0, Ai0/tol); }
            if (fabs(Ai0) > Atol) { return FALSE; }
          }
          
        /* Tests show that if a singular {2×2} matrix 
          is perturbed by adding to each element a random number in
          the range {[-eps _ +eps]}, for small {eps}, the determinant of
          the perturbed matrix will have zero mean and deviation
          {~2*eps*An} where {An} is the RMS element value. Thus if {det(A)} is on the order of
          {2*tol*An}, there is a significant chance that {A} is 
          a singular matrix with each element perturbed by 
          a relative amount of {~tol}.
        */
        double A11 = A->c[1][1]/A00;
        double A12 = A->c[1][2]/A00;
        double A21 = A->c[2][1]/A00;
        double A22 = A->c[2][2]/A00;
        if (debug) 
          { fprintf(stderr, "          A[1,1] = %24.16e * A00\n", A11); 
            fprintf(stderr, "          A[1,2] = %24.16e * A00\n", A12); 
            fprintf(stderr, "          A[2,1] = %24.16e * A00\n", A21);
            fprintf(stderr, "          A[2,2] = %24.16e * A00\n", A22);
          }
        if (debug) { }
        if (! (isfinite(A11) && isfinite(A12) && isfinite(A21) && isfinite(A22))) { return FALSE; }
        
        /* Compute the RMS elem value of the {2×2} linear submatrix: */
        double An = hypot(hypot(A11,A12), hypot(A21,A22))/2; 
        
        /* Compute the determinant of the {2×2} linear submatrix (which is that of the map): */ 
        double det = A11*A22 - A12*A21;
        if (debug) { fprintf(stderr, "          det = %24.16e = %24.16e * An * tol\n", det, det/An/tol); }

        if (det < 2*An*tol) { return FALSE; }
        return TRUE;
      }
  }

bool_t hr2_pmap_is_type(hr2_pmap_t *M, hr2_pmap_type_t type, sign_t sgn, double tol)
  { 
    double det = r3x3_det(&(M->dir));
    
    hr2_pmap_t N = (*M);
    if (det < 0)
      { if (sgn == +1)
          { return FALSE; } 
        else
          { /* Invert the sign of {N} and test for given type with {sgn=+1}: */
            hr2_pmap_invert_sign(&N);
          }
      }
    /* Now the determinant may be too small, or the sign is right: */
    
    switch (type)
      { case hr2_pmap_type_IDENTITY:
          return hr2_pmap_is_identity(&N, tol);
        case hr2_pmap_type_TRANSLATION:
          return hr2_pmap_is_translation(&N, tol);
        case hr2_pmap_type_CONGRUENCE:
          return hr2_pmap_is_congruence(&N, tol);
        case hr2_pmap_type_SIMILARITY:
          return hr2_pmap_is_similarity(&N, NAN, tol);
        case hr2_pmap_type_AFFINE:
          return hr2_pmap_is_affine(&N, tol);
        case hr2_pmap_type_GENERIC:
          return hr2_pmap_is_valid(&N, tol);
        default:
          demand(FALSE, "unrecognized or invalid map type");
      }
  }
  
void hr2_pmap_fill_rot_scale_matrix(r3x3_t *A, double ux, double uy)
  { r3x3_ident(A);
    A->c[1][1] = +ux;
    A->c[1][2] = +uy;
    A->c[2][1] = -uy;
    A->c[2][2] = +ux;
  }
