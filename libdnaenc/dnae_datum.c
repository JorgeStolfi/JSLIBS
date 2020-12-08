/* See {dnae_datum.h}. */
/* Last edited on 2014-08-27 21:05:06 by stolfilocal */

#define dnae_datum_C_COPYRIGHT \
  "Copyright © 2006  by the State University of Campinas (UNICAMP)"

#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>

#include <vec.h>
#include <bool.h>
#include <affirm.h>

#include <dnae_nucleic.h>
#include <dnae_sample.h>
#include <dnae_datum.h>
#include <dnae_seq.h>

/* IMPLEMENTATIONS: */

/* INTERPOLATION */

dnae_datum_t dnae_datum_mix
  ( double sx, 
    dnae_datum_t *fx, 
    dnae_datum_scale_t *xscale,
    double sy, 
    dnae_datum_t *fy, 
    dnae_datum_scale_t *yscale,
    dnae_datum_scale_t *rscale
  )
  { dnae_datum_t fz; 
    int c;
    for (c = 0; c < dnae_CHANNELS; c++)
      { double xc = dnae_sample_decode(fx->c[c], xscale->f[c]);
        double yc = dnae_sample_decode(fy->c[c], yscale->f[c]);
        double zc = sx*xc + sy*yc;
        fz.c[c] = dnae_sample_encode(zc, rscale->f[c]);
      }
    return fz;
  }

double dnae_datum_euc_distsq(dnae_datum_t *fx, dnae_datum_scale_t *xscale, dnae_datum_t *fy, dnae_datum_scale_t *yscale)
  { double d2 = 0;
    int c;
    for (c = 0; c < dnae_CHANNELS; c++)
      { double xc = dnae_sample_decode(fx->c[c], xscale->f[c]);
        double yc = dnae_sample_decode(fy->c[c], yscale->f[c]);
        double dc = xc - yc;
        d2 += dc*dc;
      }
    return d2;
  }

/* QUEST OF THE CORRECT DATUM DISTANCE ESTIMATORS */

/* RIGOROUS ANALYSIS

  Assume given {N} /sample smoothing weights/ {ws[0..N-1]}
  and {N} /distance smoothing weights/ {wd[0..N-1]}.
  
  Suppose we have two datums {X,Y}, which are the weighted averages of
  raw datums strings {x[0..N-1]} and {y[0..N-1]}, respectively:
  
    {X = H(x) = SUM{ ws[i]*x[i] : i = 0..N-1}}
    {Y = H(y) = SUM{ ws[i]*y[i] : i = 0..N-1}}
  
  We are interested in the /local squared difference/ {S2(x,y)} of the
  strings {x,y}, defined as the weighted average of the squared
  differences of corresponding datums, namely
  
    {S2(x,y) = SUM{ wd[i]*|x[i] - y[i]|^2 : i = 0..N-1 }}
  
  We need the expected value of the random variable {S2(x,y)} assuming
  that the datums {x[0..N-1]} and {y[0..N-1]} are taken from some
  fixed distribution, given only the datums {X} and {Y}. Namely, we
  wish to compute
  
    {D2(X,Y) = E(S2(x,y) @ H(x)~X & H(y)~Y)}
    
  where {E(A @ B)} denotes the expectation of random variable {A}
  given condition {B}.  By the expectation-of-sum theorem have
  
    {D2(X,Y) = SUM{ wd[i] * F2[i](X,Y): i = 0..N-1 }}
    
  where
    
    {F2[i](X,Y) = E(|x[i] - y[i]|^2 @ H(x)~X & H(y)~Y)}
    
  In our case, the vectors {x[i]} and {y[i]} are taken from the set
  {V} of the vertices of a regular tetrahedron whose edge has length
  {sqrt(8)}. So {|x[i] - y[i]|^2} is 0 if {x[i] = y[i]}, and 8 if
  {x[i] != y[i]}. So
  
    {F2[i](X,Y) = 8 * (1 - Q[i](X,Y))}
    
  where
  
    {Q[i](X,Y) = Pr(x[i] = y[i] @ H(x)~X & H(y)~Y)}
  
  Assuming that, /a priori/, the choice of
  {x[i]} is indepedent of {y[0..N-1]}, and vice-versa,
  we have
  
    {Q[i](X,Y) = SUM{ Pr(x[i]=u @ H(x)~X)*Pr(y[i]=u @ H(y)~Y) : u in V }}
  
  So we need to estimate 
  
    {G[i](X,u) = Pr(x[i]=u @ H(x)~X)}
  
  Good question. 
  
*/

/* QUADRATIC INTERPOLATION OF SPECIAL CASES
  
  Another approach is to look for a second-degree approximation
  {T2(X,Y)} to the function {D2(X,Y)}, that has the right symmetries
  and right values for special cases of {X,Y}.
  
  The function {T2} should be symmetric under arbitrary channel
  permutations and exchange of {X} with {Y}. Therefore it must be a
  linear combination of the homogeneous functions of degree {0,1,2}
  with those symmetries:
  
    {L[0](X,Y) = 1}
    {L[1](X,Y) = X0*Y0 + X1*Y1 + X2*Y2 = <X|Y>}
    {L[2](X,Y) = X0^2 + X1^2 + X2^2 + Y0^2 + Y1^2 + Y2^2 = |X|^2 + |Y|^2}
    {L[3](X,Y) = X0 + X1 + X2 + Y0 + Y1 + Y2}
    {L[4](X,Y) = X0*X1 + X1*X2 + X2*X0 + Y0*Y1 + Y1*Y2 + Y2*Y0}
    {L[5](X,Y) = X0*Y1 + X1*Y2 + X2*Y0 + X1*Y0 + X2*Y1 + X0*y2}
    
  where {<X|Y>} is the dot product of the two vectors.
  Moreover, the function {T2} should be unchanged by simultaneous 
  negation of {X0,X1,Y0,Y1}.  That excludes {L[3],L[4],L[5]}, leaving only
  {L[0],L[2],L[1]}; namely, {1,|X|^2+|Y|^2,<X|Y>}. 
    
  We should also require that {T2(X,Y)} is exact when {X} and {Y}
  are `pure' datums (corners of the datum simplex).  The three 
  basis functions have the following values for such points:
  
    If {X==Y}:   {L[0] = +1}  {L[1] = +3}   {L[2] = +6}  {D2 = 0}
    If {X!=Y}:   {L[0] = +1}  {L[1] = -1}   {L[2] = +6}  {D2 = 8}
  
  Another known value of {D2(X,Y)} is when {X} and {Y} are the origin
  {O=(0,0,0)} (the center of the datum simplex). By symmetry, each raw
  datum {x[i]} is equally likely to be any of the four simplex corners
  'A', 'T', 'C', 'G'; and ditto for {y[i]}. In that case, we have
  
    {G[i](X,u) = Pr(x[i]=u @ H(x)~X) = 1/4}
    
    {Q[i](X,Y) = SUM{ G[i](X,u)*G[i](Y,u) : u in V } = 4/16 = 1/4}
  
    {F2[i](X,Y) = 8 * (1 - Q[i](X,Y)) = 8*(3/4) = 6} 
    
    {D2(X,Y) = SUM{ wd[i] * F2[i](X,Y) } = 6 * SUM{ wd[i] } = 6 }
    
  The data point is therefore
  
    If {X=Y=O}:  {L[0] = +1}  {L[1] = 00}   {L[2] = 00}  {D2 = 6}
    
  We can also consider the case when {X} and {Y} are edge midpoints,
  namely two of the six points {(±1,0,0), (0,±1,0), (0,0,±1)}.
  
  Then {x[i]} is equally likely to be either one of the two endpoints 
  of the edge, and ditto for {y[i]}.  We have three cases to consider,
  (a) the two edges coincide (6 cases in 36), (b) the two edges 
  share a vertex (24 in 36) and (c) the two edges are opposite (6 in 36).
  
  In all three cases we have {G[i](X,u) = 1/2} when {u} is an endpoint of the 
  edge that contains the point {X}.
  
  In case (a), {Q[i](X,Y) = 1/4 + 1/4 + 0 + 0 = 1/2}, {F2[i](X,Y) = 4}, {D2(X,Y) = 4}.
  
  In case (b), {Q[i](X,Y) = 1/4 + 0 + 0 + 0 = 1/4}, {F2[i](X,Y) = 6}, {D2(X,Y) = 6}.
  
  In case (c), {Q[i](X,Y) = 0 + 0 + 0 + 0 = 0}, {F2[i](X,Y) = 8}, {D2(X,Y) = 8}.
  
  The edge-related data points are therefore

    Case (a):    {L[0] = +1}  {L[1] = +1}   {L[2] = +2}  {D2 = 4}
    Case (b):    {L[0] = +1}  {L[1] = 00}   {L[2] = +2}  {D2 = 6}
    Case (c):    {L[0] = +1}  {L[1] = -1}   {L[2] = +2}  {D2 = 8}
  
  By considering the symmetry requirements and all data points above,
  we conclude that the only solution is {T2(X,Y) = 6*L[0](X,Y) -
  2*L[1](X,Y) = 6 - 2*<X|Y>}. */
  
/* QUADRATIC LEAST-SQUARES FITTING OF RANDOM CASES
  
  The empirical approach is as follows: choose weights {ws[0..N-1]} and
  {wd[0..N-1]}. Generate a large number of test cases {x[0..N-1]} and
  {y[0..N-1]}. For each case, compute {X=H(x)}, {Y=h(y)}, and
  {Z=S2(x,y)}. Now find a degree-2 function {T2(X,Y)} using the basis
  functions {L[0..5]} that provides the best fit for {Z} over all
  those data points.
  
  These numerical experiments indicate that the optimum quadratic distance
  estimator is {T2(X,Y) = 6 - K*<X|Y>} where {K} is exactly 2 for uniform 
  weights ({ws[i] = 1/N} for all {i}) and somewhat higher than 2
  for other weights.  For binomial weights {ws[i] = choose(N-1,i)/2^N}
  the optimal {K} is given by the table
   
     { N =     1     2     3     4     5     6     7  } 
     {     ----- ----- ----- ----- ----- ----- -----  }
     { K = 2.000 2.000 2.222 2.242 2.255 2.269 2.277  }
  
  In these tests we used {wd[i] = ws[i]}. The above results may mean
  that {D2} is more predictable when the {wd[i]} are different
  from the {ws[i]}.  Which ones?
  
  */

double dnae_datum_diffsq(dnae_datum_t *fx, dnae_datum_scale_t *xscale, dnae_datum_t *fy, dnae_datum_scale_t *yscale)
  { /* We use the formula {T2(X,Y) = 6 - 2*<X|Y>} to estimate {S2(X,Y)}. */
    double p = 0;  /* Dot product {<fx|fy>} */
    int c;
    for (c = 0; c < dnae_CHANNELS; c++)
      { double xc = dnae_sample_decode(fx->c[c], xscale->f[c]);
        double yc = dnae_sample_decode(fy->c[c], yscale->f[c]);
        p += xc*yc;
      }
    double T2 = 6 - 2*p;
    /* Just in case that we get datums outside the simplex: */
    if (T2 < 0) { T2 = 0; }
    /* Return {T2} normalized so that the maximum abs value is 1.
      Assuming that {X} and {Y} are in the datum simplex,
      the range of {<X|Y>} is {[-1 _ +3]}, so the maximum value 
      of {T2(X,Y)} is {6 - 2(-1) = 8}. */
    return T2/8;
  }

double dnae_datum_step_diffsq
  ( dnae_datum_t *fx0, 
    dnae_datum_t *fy0, 
    dnae_datum_t *fx1, 
    dnae_datum_t *fy1, 
    dnae_datum_scale_t *xscale, 
    dnae_datum_scale_t *yscale
  )
  { /* Compute the mean value {mp} of {<X(t)|Y(t)>} when
      {X(t)} interpolates linearly between {fx0} and {fx1}, and
      {Y(t)} interpolates between {fy0} and {fy1},
      as {t} ranges from 0 to 1: */
    double mp = 0;  /* Mean value of {<X(t)|Y(t)>} */
    int c;
    for (c = 0; c < dnae_CHANNELS; c++)
      { double x0c = dnae_sample_decode(fx0->c[c], xscale->f[c]);
        double y0c = dnae_sample_decode(fy0->c[c], yscale->f[c]);
        double x1c = dnae_sample_decode(fx1->c[c], xscale->f[c]);
        double y1c = dnae_sample_decode(fy1->c[c], yscale->f[c]);
        mp += (2*x0c*y0c + x0c*y1c + x1c*y0c + 2*x1c*y1c)/6;
      }
    /* Compute the mean value of {T2(X(t),Y(t))}, as in {dnae_datum_diffsq}: */
    double mT2 = 6 - 2*mp;
    if (mT2 < 0) { mT2 = 0; }
    /* Return {mT2} normalized, as in {dnae_datum_diffsq}: */
    return mT2/8;
  }

double dnae_datum_half_step_diffsq
  ( dnae_datum_t *fx0, 
    dnae_datum_t *fy0, 
    dnae_datum_t *fx1, 
    dnae_datum_t *fy1, 
    dnae_datum_scale_t *xscale, 
    dnae_datum_scale_t *yscale
  )
  { 
    /* Compute the weighted integral {iwp} of {w(t)*<X(t)|Y(t)>} when {X(t)} interpolates
      between {fx0} and {fx1}, and {Y(t)} interpolates between {fy0} and {fy1},
      as {t} ranges from 0 to 1, where {w} is the half-tent weight function,
      {w(t)=1-t}: */
    double iwp = 0; /* Integral of {w(t)*<X(t)|Y(t)>} over [0_1]. */
    int c;
    for (c = 0; c < dnae_CHANNELS; c++)
      { double x0c = dnae_sample_decode(fx0->c[c], xscale->f[c]);
        double y0c = dnae_sample_decode(fy0->c[c], yscale->f[c]);
        double x1c = dnae_sample_decode(fx1->c[c], xscale->f[c]);
        double y1c = dnae_sample_decode(fy1->c[c], yscale->f[c]);
        iwp += (3*x0c*y0c + x0c*y1c + x1c*y0c + x1c*y1c)/12;
      }
    double iw = 0.5; /* Integral of {w(t)} over [0_1]. */
    /* Compute the mean value of {w(t)*T2(X(t),Y(t))}, as in {dnae_datum_diffsq}: */
    double iwT2 = 6*iw - 2*iwp;
    if (iwT2 < 0) { iwT2 = 0; }
    /* Return {iwT2} normalized just as in {dnae_datum_diffsq}: */
    return iwT2/8;
  }

void dnae_datum_to_nucleic_densities(dnae_datum_t *d, dnae_datum_scale_t *dscale, double *A, double *T, double *C, double *G)
  { double d0 = dnae_sample_decode(d->c[0], dscale->f[0]);
    double d1 = dnae_sample_decode(d->c[1], dscale->f[1]);
    double d2 = dnae_sample_decode(d->c[2], dscale->f[2]);
    (*A) = ((+d0) + (+d1) + (+d2) + 1)/4;
    (*T) = ((-d0) + (-d1) + (+d2) + 1)/4;
    (*C) = ((+d0) + (-d1) + (-d2) + 1)/4;
    (*G) = ((-d0) + (+d1) + (-d2) + 1)/4;
  }

void dnae_datum_decoded_from_nucleic_char(char b, int *d)
  { int A, T, C, G;
    dnae_nucleic_value(b, &A, &T, &C, &G);
    d[0] = (A-T)+(C-G);
    d[1] = (G-C)+(A-T);
    d[2] = (A+T)-(C+G);
  }

dnae_datum_t dnae_datum_encoded_from_nucleic_char(char b)
  { int d[3];
    dnae_datum_decoded_from_nucleic_char(b, d);
    dnae_datum_t r;
    double scale = dnae_NUCLEIC_RAW_SCALE;
    r.c[0] = dnae_sample_encode((double)d[0], scale);
    r.c[1] = dnae_sample_encode((double)d[1], scale);
    r.c[2] = dnae_sample_encode((double)d[2], scale);
    return r;
  }

dnae_datum_vec_t dnae_datum_vec_from_nucleic_string(char *s)
  { /* Get count {nbas} of nucleotide characters: */
    int nbas = 0; 
    char *p = s;
    while ((*p) != 0) { if (is_dna_basis(*p)) { nbas++; } p++; }
    /* Allocate datum vector: */
    dnae_datum_vec_t dv = dnae_datum_vec_new(nbas);
    /* Convert to numeric and encode: */
    int j = 0;
    p = s;
    while ((*p) != 0)
      { char c = *p;
        if (is_dna_basis(c)) 
          { dv.e[j] = dnae_datum_encoded_from_nucleic_char(c);
            j++;
          }
        p++;
      }
    assert(j == nbas);
    assert(dv.ne == nbas);
    return dv;
  }

vec_typeimpl(dnae_datum_vec_t,dnae_datum_vec,dnae_datum_t);

void dnae_datum_encoded_write(FILE *wr, dnae_datum_t *d, char *lp, char *sep, char *rp)
{ if (lp != NULL) { fputs(lp, wr); }
    int c;
    for (c = 0; c < dnae_CHANNELS; c++)
      { if ((c != 0) && (sep != NULL)) { fputs(sep, wr); }
        fprintf(wr, "%+6d", d->c[c]); 
      } 
    if (rp != NULL) { fputs(rp, wr); }
  }
    
void dnae_datum_decoded_write(FILE *wr, dnae_datum_t *d, dnae_datum_scale_t *dscale, char *lp, char *sep, char *rp)
  { if (lp != NULL) { fputs(lp, wr); }
    int c;
    for (c = 0; c < dnae_CHANNELS; c++)
      { if ((c != 0) && (sep != NULL)) { fputs(sep, wr); }
        fprintf(wr, "%+11.7f", dnae_sample_decode(d->c[c], dscale->f[c])); 
      } 
    if (rp != NULL) { fputs(rp, wr); }
  }

