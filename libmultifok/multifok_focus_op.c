/* See {multifok_focus_op.h}. */
/* Last edited on 2023-01-09 10:45:37 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <vec.h>
#include <i2.h>
#include <wt_table.h>
#include <affirm.h>
#include <bool.h>

#include <multifok_focus_op.h>

void multifok_focus_op_set_samples_3x3
  ( double x[], 
    double scale,
    double x00, double x01, double x02, 
    double x10, double x11, double x12, 
    double x20, double x21, double x22 
  );
  /* Specific for {NW = 3} ({NS = 9}): sets the elements {x[0..NS-1]} to the 
    values {x00,x01,x02,x10,...,x22} multiplied by {scale}. */

char *multifok_focus_op_mop_code(int32_t i);
  /* Returns a newly allocated string with the value of {i} encoded as
    "o" for {i=0}, "m{k}" for {i = -k}, "p{i}" for {i = +k}, with {k} in
    {1..NW/2}; except that {k} is supressed if it is 1. */

char *multifok_focus_op_sample_name(char *tag, int32_t ix, int32_t iy);
  /* Returns a newly allocated string with the index pair {ix,iy}
    encoded as "{tage}{X}{Y}", where {X} and {Y} are the indices {ix}
    and {iy} encoded with {multifok_focus_op_mop_code}. */

void multifok_focus_op_sample_names(int32_t NW, char *tag, char *name[]);
  /* Sets {name[0..NS-1]} to {NS = NW*NW} newly allocated strings that are scrutable names of
    the samples in the window. The names have the form "{tag}{X}{Y}" where  {X} and {Y} are the indices {ix}
    and {iy} of the sample relative to the center sample, encoded with {multifok_focus_op_mop_code}. */

/* NON-NORMALIZED BASES 

  The procedures below generate bases of various types, but ignoring the sample weights
  (so the elements are not orthonormal).

*/

void multifok_focus_op_basis_canonical(int32_t NB, int32_t NW, double *phi[], double wc[], char *name[]);
  /* Requires {NW} odd and {NB=NW*NW}. Sets {phi[0..NB-1]} to the canonical basis element 
    that is 1 at sample {kb} and 0 everywhere else.  All weights {wc} will be 1.
    The element names will be "P{X}{Y}" where {X} and {Y} are as in {multifok_focus_op_sample_names}. */

void multifok_focus_op_basis_deriv_ops(int32_t NB, int32_t NW, double *phi[], double wc[], char *name[]);
  /* Requires {NW=3} and {NB=NW*NW}. Sets {phi[0..NB-1]} with the sample coefficients for computing the derivative operators up to
    order 2, plus 3-saddles and checker.  The weight {wc[0]} of the window mean will be 0, the others will be 1. 
    The names will be "F" (window sum), "DX", "DY", "DXX", etc. */
  
void multifok_focus_op_basis_hartley(int32_t NB, int32_t NW, double *phi[], double wc[], char *name[]);
  /* Requires {NW} odd and {NB=NW*NW}. Sets {phi[0..NB-1]} to the Hartley waves with absolute {X} and {Y} frequencies 
    up tp {NW/2}.   Each coefficient {wc[0..NB-1]} will be the square of the frequency
    vector modulus. The element names will be "H{X}{Y}" where {X} and {Y} are as in 
    {multifok_focus_op_sample_names}, only referring to frequencies instead of sample indices. */

/* BASIS NORMALIZATION */

int32_t multifok_focus_op_orthize(int32_t NW, double ws[], int32_t NB, double *phi[], double wc[], char *name[]);
  /* Makes the elements {phi[0..NB-1]} orthonormal with respect to the dot product with weights {ws[0..NS-1]} where
    {NS=NW*NW}.
    
    The procedure may conclude that some basis elements are nearly
    dependent on others, and exclude them from the normalized basis. In
    that case, the number of non-discarded elements {NK} is returned as
    result, and those elements are rearranged to be {phi[0..NK-1]}.
    
    The procedure also adjusts the coefficient weights {wc[0..NB-1]} and
    names {name[0..NB-1]} to account for rearrangements and the
    normalization of the respective elements. On output the coefficient
    {wc[k]} and {name[k]} still corresponds to {phi[k]} for {k} in
    {0..NK-1}. */
    
void multifok_focus_op_normalize_prod_range(int32_t NW, double ws[], int32_t NB, double *phi[], double wc[]);
  /* Scales each basis element {phi[0..NB-1]} so that its dot product with any window in {[0_1]^NS} is 
    in the range {[-1 _ +1]}. */

/* IMPLEMENTATIONS */
  
int32_t multifok_focus_op_num_samples(int32_t NW)
  { demand(((NW % 2) == 1) && (NW >= 3), "window size must be odd");
    return NW*NW;
  }
 
double *multifok_focus_op_sample_weights(int32_t NW)
  { /* Create a unidimensional weight table: */
    double u[NW]; /* Unidimensional weights. */
    bool_t normalize = TRUE;
    wt_table_fill_binomial(NW, u, normalize);
    
    /* Now fill the bidimensional table: */
    int32_t NS = multifok_focus_op_num_samples(NW);
    double *ws = notnull(malloc(NS*sizeof(double)), "no mem");
    for (int32_t iy = 0; iy < NW; iy++)
      { for (int32_t ix = 0; ix < NW; ix++)
          { double wxy = u[iy]*u[ix];
            ws[iy*NW + ix] = wxy;
          }
      }

    return ws;
  }

void multifok_focus_op_normalize_samples(int32_t NW, double x[], double ws[], double noise)
  {
    int32_t NS = NW*NW;
    double sum_wx = 0;
    double sum_w = 0;
    for (int32_t ks = 0; ks < NS; ks++) 
      { sum_wx += ws[ks]*x[ks]; 
        sum_w += ws[ks];
      }
    double avg = sum_wx/sum_w;
    double ns2 = noise*noise;
    double sum_wd2 = 1.0e-200; /* In case the actual dev and {noise} are zero. */
    for (int32_t ks = 0; ks < NS; ks++) 
      { double d = x[ks] - avg;
        sum_wd2 += ws[ks]*(d*d + ns2);
      }
    double dev = sqrt(sum_wd2/sum_w);
    for (int32_t ks = 0; ks < NS; ks++) 
      { x[ks] = (x[ks] - avg)/dev; }
  
  }

double multifok_focus_op_prod(int32_t NW, double x[], double y[], double ws[])
  { int32_t NS = multifok_focus_op_num_samples(NW);
    double prod = 0.0;
    for (int32_t k = 0; k < NS; k++) 
      { double xk = x[k];
        double yk = y[k];
        double wk = ws[k];
        prod += wk*xk*yk;
      }
    return prod;
  }

double multifok_focus_op_dist_sqr(int32_t NW, double x[], double y[], double ws[])
  { int32_t NS = multifok_focus_op_num_samples(NW);
    double d2 = 0.0;
    for (int32_t k = 0; k < NS; k++) 
      { double dk = x[k] - y[k];
        double wk = ws[k];
        d2 += wk*dk*dk;
      }
    return d2;
  }

double multifok_focus_op_score_from_basis(int32_t NW, double x[], double noise, double ws[], int32_t NB, double *phi[], double wc[])
  { 
    int32_t NS = multifok_focus_op_num_samples(NW);
    double c[NS];
    multifok_focus_op_normalize_samples(NW, x, ws, noise);
    for (int32_t kb = 0; kb < NB; kb++) { c[kb] = multifok_focus_op_prod(NW, x, phi[kb], ws); }
    double score = 0.0;
    for (int32_t kb = 0; kb < NB; kb++) { score += wc[kb]*c[kb]*c[kb]; }
    return score;
  }

void multifok_focus_op_basis_make
  ( int32_t NW, 
    double ws[],
    multifok_focus_op_basis_type_t bType,
    bool_t ortho,
    int32_t *NB_P, 
    double ***phi_P, 
    double **wc_P,
    char ***name_P
  )
  { int32_t NS = multifok_focus_op_num_samples(NW);
    
    int32_t NB = NS;  /* May be reduced later. */
    double *wc = notnull(malloc(NB*sizeof(double)), "no mem"); 
    double **phi = notnull(malloc(NB*sizeof(double*)), "no mem"); 
    for (int32_t i = 0; i < NS; i++) { phi[i] = notnull(malloc(NS*sizeof(double)), "no mem"); }
    char **name = notnull(malloc(NB*sizeof(char*)), "no mem"); 
    
    if (bType == multifok_focus_op_basis_type_NONE)
      { /* Canonical basis: */
        multifok_focus_op_basis_canonical(NB, NW, phi, wc, name);
      }
    else if (bType == multifok_focus_op_basis_type_DIFF)
      { /* Differential operator basis: */
        demand(NW == 3, "diff basis requires {NW=3}");
        multifok_focus_op_basis_deriv_ops(NB, NW, phi, wc, name);
      }
    else if (bType == multifok_focus_op_basis_type_HART)
      { /* Hartley wave basis: */
        multifok_focus_op_basis_hartley(NB, NW, phi, wc, name);
      }
    else
      { assert(FALSE); }

    if (ortho)
      { /* Orthonormalize basis. */
        int32_t NK = multifok_focus_op_orthize(NW, ws, NB, phi, wc, name);
        /* Free unused space: */
        if (NK < NB)
          { /* Free unused elements: */
            for (int32_t k = NK; k < NB; k++) { free(phi[k]); }
            phi = realloc(phi, NK*sizeof(double*));
            wc = realloc(wc, NK*sizeof(double));
          }
        NB = NK;
      }
    else
      { /* Scale each element so that the dot product with any window in {[0_1]^NS} is always in {[-1 _ +1]}: */
        multifok_focus_op_normalize_prod_range(NW, ws, NB, phi, wc);
     }
    
    (*NB_P) = NB;
    (*phi_P) = phi;
    (*wc_P) = wc;
    (*name_P) = name;
  }

void multifok_focus_op_basis_free(int32_t NB, double **phi, double *wc, char **name)
  {
    for (int32_t kb = 0; kb < NB; kb++) { free(phi[kb]); free(name[kb]); }
    free(phi);
    free(wc);
    free(name);
  }

void multifok_focus_op_basis_canonical(int32_t NB, int32_t NW, double *phi[], double wc[], char *name[])
  {
    int32_t NS = NW*NW;
    assert(NB == NS);

    for (int32_t kb = 0; kb < NB; kb++)
      { for (int32_t ks = 0; ks < NS; ks++)
          { phi[kb][ks] = (kb == ks ? 1.0 : 0.0); }
        wc[kb] = 1.0;
      }
    multifok_focus_op_sample_names(NW, "P", name);
  }
  
void multifok_focus_op_basis_deriv_ops(int32_t NB, int32_t NW, double *phi[], double wc[], char *name[])
  {
    assert(NB == 9);
    assert(NW == 3);
    
    auto void setbas
      ( int32_t ib, 
        double wc_ib,
        char *name_ib,
        double x00, double x01, double x02, 
        double x10, double x11, double x12, 
        double x20, double x21, double x22
      );
      /* Sets the name {name[ib]} to {nm}, the coeff weight {wc[ib]} to {w}, and the basis element 
        samples {phi[ib][0..8]} to {x00, ..., x22}. */

    setbas(0, 00.000, "F",  +1.0, +1.0, +1.0,  +1.0, +1.0, +1.0,  +1.0, +1.0, +1.0f );  /* Mean. */
    setbas(1, -0.500, "DX", -1.0, 00.0, +1.0,  -1.0, 00.0, +1.0,  -1.0, 00.0, +1.0f );  /* d/dX. */
    setbas(2, -0.500, "DY", -1.0, -1.0, -1.0,  00.0, 00.0, 00.0,  +1.0, +1.0, +1.0f );  /* d/dY. */
    setbas(3, -0.800, "DXX",-1.0, +1.0, -1.0,  -1.0, +1.0, -1.0,  -1.0, +1.0, -1.0f );  /* d^2/dX^2. */
    setbas(4, -0.800, "DYY",-1.0, -1.0, -1.0,  +1.0, +1.0, +1.0,  -1.0, -1.0, -1.0f );  /* d^2/dY^2. */
    setbas(5, -1.000, "DXY",-1.0, 00.0, +1.0,  00.0, 00.0, 00.0,  +1.0, 00.0, -1.0f );  /* d^2/dXdY. */
    setbas(6, -0.640, "S3", -1.0, +1.0, -1.0,  00.0, 00.0, 00.0,  +1.0, -1.0, +1.0f );  /* Saddle3. */
    setbas(7, -0.640, "C3", +1.0, 00.0, -1.0,  -1.0, 00.0, +1.0,  +1.0, 00.0, -1.0f );  /* Cosaddle3. */
    setbas(8, -0.333, "Q",  +1.0, -1.0, +1.0,  -1.0, +1.0, -1.0,  +1.0, -1.0, +1.0f );  /* Checker. */

    return;
    
    void setbas
      ( int32_t ib, 
        double wc_ib,
        char *name_ib,
        double x00, double x01, double x02, 
        double x10, double x11, double x12, 
        double x20, double x21, double x22
      )
      {
        wc[ib] = wc_ib;
        name[ib] = NULL; asprintf(&(name[ib]), "%s", name_ib); /* Make sure {name[ib]} is a heap-allocated string. */
        double *pib = phi[ib];
        pib[0] = x00;
        pib[1] = x01;
        pib[2] = x02;
        pib[3] = x10;
        pib[4] = x11;
        pib[5] = x12;
        pib[6] = x20;
        pib[7] = x21;
        pib[8] = x22;
      }
  }

void multifok_focus_op_basis_hartley(int32_t NB, int32_t NW, double *phi[], double wc[], char *name[])
  {
    assert((NW % 2 ) == 1);
    int32_t HW = NW/2;

    int32_t NS = NW*NW;
    assert(NB == NS);
    
    /* Generate tables {fre[0..NB-1]} of frequency vectors: */
    i2_t freq[NB];
    int32_t kb = 0; /* Basis element index. */
    for (int32_t fy = -HW; fy <= +HW; fy++)
      { for (int32_t fx = -HW; fx <= +HW; fx++)
          { freq[kb] = (i2_t){{ fx, fy}}; kb++; }
      }
    assert(kb == NB);
    
    /* Sort the table by increasing freq modulus: */
    for (int32_t ib = 0; ib < NB; ib++)
      { /* The {kb} smallest freqs are {freq[0..kb-1]} */
        int32_t jbMin = ib;
        for (int32_t jb = ib+1; jb < NB; jb++)
          { if (i2_norm_sqr(&(freq[jb])) < i2_norm_sqr(&(freq[jbMin])))
              { jbMin = jb; }
          }
        if (jbMin != ib)
          { i2_t ft = freq[ib]; freq[ib] = freq[jbMin]; freq[jbMin] = ft; }
      }
    
    /* Now fill the basis: */
    for (int32_t kb = 0; kb < NB; kb++)
      { 
        int32_t fx = freq[kb].c[0];
        int32_t fy = freq[kb].c[1];
        for (int32_t iy = -HW; iy <= +HW; iy++)
          { for (int32_t ix = -HW; ix <= +HW; ix++)
              { int32_t ks = (iy+HW)*NW + (ix+HW);
                double z = 0.125 + ((double)(fx*ix + fy*iy))/((double)NW); /* Hartley phase */
                double f = sin(2*M_PI*z);
                phi[kb][ks] = f;
              }
          }
        wc[kb] = hypot(fx,fy);
        name[kb] = multifok_focus_op_sample_name("H", fx, fy);
      }
  }
  
char *multifok_focus_op_mop_code(int32_t i)
  { char *u = NULL;
    char s = (i < 0 ? 'm' : (i > 0 ? 'p' : 'o'));
    if (abs(i) <= 1) { asprintf(&u, "%c", s); } else { asprintf(&u, "%c%d", s, abs(i)); }
    return u;
  }
  
char *multifok_focus_op_sample_name(char *tag, int32_t ix, int32_t iy)
  {
    char *ux = multifok_focus_op_mop_code(ix);
    char *uy = multifok_focus_op_mop_code(iy);
    char *uxy = NULL;
    asprintf(&uxy, "%s%s%s", tag, ux, uy);
    free(ux); 
    free(uy);
    return uxy;
  }    

void multifok_focus_op_sample_names(int32_t NW, char *tag, char *name[])
  {
    assert((NW % 2 ) == 1);
    int32_t HW = NW/2;
    int32_t NS = NW*NW;
    
    /* Set {u[i+HW]} to text versions of coordinate {i}, for {i} in {-HW..+HW}: */
    char *u[NW];
    for (int32_t i = -HW; i <= +HW; i++)
      { u[i + HW] = multifok_focus_op_mop_code(i); }
    
    /* Assemble sample names: */
    int32_t ks = 0; /* Windos sample index. */
    for (int32_t iy = -HW; iy <= +HW; iy++)
      { for (int32_t ix = -HW; ix <= +HW; ix++)
          { char *uxy = NULL;
            asprintf(&uxy, "%s%s%s", tag, u[ix+HW], u[iy+HW]);
            name[ks] = uxy;
            ks++;
         }
      }
    assert(ks == NS);
    
    /* Release the auxiliary coordinate names: */
    for (int32_t i = -HW; i <= +HW; i++){ free(u[i + HW]); }
  }

void multifok_focus_op_normalize_prod_range(int32_t NW, double ws[], int32_t NB, double *phi[], double wc[])
  {
    int32_t NS = NW*NW;
    for (int32_t kb = 0; kb < NB; kb++)
      { /* Determne the range {[sMin _ sMax]} of the dot product of {phi[kb]} and: */
        double sMax = 0.0; /* Sum of positive basis coordinates. */
        double sMin = 0.0; /* Sum of negative basis coordinates. */
        for (int32_t ks = 0; ks < NS; ks++) 
          { double wpk = ws[ks]*phi[kb][ks];
            if (wpk > 0.0) { sMax += wpk; }
            if (wpk < 0.0) { sMin += wpk; }
          }
        assert((sMin <= 0) && (sMax >= 0));
        /* Rescale basis samples to ensure the dot product is in the range {[-1 _ +1]}: */
        double scale = fmax(-sMin, sMax);
        assert(scale > 0.0);
        for (int32_t ks = 0; ks < NS; ks++) { phi[kb][ks] /= scale; }
      }
  }
 
int32_t multifok_focus_op_orthize(int32_t NW, double ws[], int32_t NB, double *phi[], double wc[], char *name[])
  { 
    bool_t verbose = TRUE;
    int32_t NS = multifok_focus_op_num_samples(NW);
    int32_t NK = 0; /* Nymber of basis elements that were kept. */
    for (int32_t kb = 0; kb < NB; kb++)
      { 
        /* Original norm of element {phi[kb]}: */
        double d_old = sqrt(multifok_focus_op_prod(NW, phi[kb], phi[kb], ws));
        
        /* Make element {phi[kb]} orthogonal to the previous elements:*/
        for (int32_t rb = 0; rb < NK; rb++)
          { double ci = multifok_focus_op_prod(NW, phi[kb], phi[rb], ws);
            for (int32_t js = 0; js < NS; js++)
              { phi[kb][js] = phi[kb][js] - ci*phi[rb][js]; }
          }
        /* Norm of residual element {phi[kb]} after removing previous elements: */
        double d_new = sqrt(multifok_focus_op_prod(NW, phi[kb], phi[kb], ws));
        if (verbose) { fprintf(stderr, "d old = %16.12f new = %16.12f\n", d_old, d_new); }
        if (d_new > 1.0e-8*d_old)
          { /* Residual seems significant: */
            for (int32_t js = 0; js < NS; js++)
              { phi[NK][js] = phi[kb][js]/d_new; }
            wc[NK] = wc[kb]*d_new;
            name[NK] = name[kb];
            NK++;
          }
        else
          { /* Element was mostly linear comb of previous ones, discard: */
            if (verbose) { fprintf(stderr, "!! eliminated element {phi[%d]}\n", kb); }
          }
      }
    
    return NK;
  }

void multifok_focus_op_set_samples_3x3
  ( double x[], 
    double scale,
    double x00, double x01, double x02, 
    double x10, double x11, double x12, 
    double x20, double x21, double x22 
  )
  { x[0] = scale * x00;
    x[1] = scale * x01;
    x[2] = scale * x02;
    x[3] = scale * x10;
    x[4] = scale * x11;
    x[5] = scale * x12;
    x[6] = scale * x20;
    x[7] = scale * x21;
    x[8] = scale * x22;
  }

void multifok_focus_op_remap_samples(int32_t NW, double x[], double ws[], int32_t NB, double *phi[], double c[])
  { int32_t NS = multifok_focus_op_num_samples(NW);
    for (int32_t kb = 0; kb  < NS; kb++) 
      { c[kb] = multifok_focus_op_prod(NW, x, phi[kb], ws); }
  }

void multifok_focus_op_basis_print
  ( FILE *wr,
    int32_t NW, 
    int32_t NB, 
    double *phi[],
    double wc[],
    char *name[]
  )
  { int32_t NS = multifok_focus_op_num_samples(NW);
    
    /* Print basis: */
    fprintf(wr, "--- basis ---------------------------------------------\n");
    fprintf(wr, "%3s %-12s %12s ", "", "", "");
    for (int32_t j = 0; j < NS; j++) { fprintf(wr, " %10d", j); }
    fprintf(wr, "\n");
    for (int32_t ib = 0; ib < NB; ib++)
      { fprintf(wr, "%3d %-12s %12.6f ", ib, name[ib], wc[ib]);
        for (int32_t js = 0; js < NS; js++)
          { double pij = phi[ib][js];
            fprintf(wr, " %+10.6f", pij);
          }
        fprintf(wr, "\n");
      }
    fprintf(wr, "-------------------------------------------------------\n");
  }

void multifok_focus_op_basis_ortho_check
  ( FILE *wr,
    int32_t NW, 
    double ws[], 
    int32_t NB, 
    double *phi[]
  )
  { /* Print basis products: */
    fprintf(wr, "--- momentum matrix -----------------------------------\n");
    fprintf(wr, "%3s", "");
    for (int32_t jb = 0; jb < NB; jb++) { fprintf(wr, " %10d", jb); }
    fprintf(wr, "\n");
    for (int32_t ib = 0; ib < NB; ib++)
      { fprintf(wr, "%3d", ib);
        for (int32_t jb = 0; jb <= ib; jb++)
          { double pij = multifok_focus_op_prod(NW, phi[ib], phi[jb], ws);
            fprintf(wr, " %+10.6f", pij);
          }
        fprintf(wr, "\n");
      }
    fprintf(wr, "-------------------------------------------------------\n");
  }
    
void multifok_focus_op_check(int32_t NW, multifok_focus_op_basis_type_t bType, bool_t ortho)
  { double *ws = multifok_focus_op_sample_weights(NW);
    
    int32_t NB;
    double *wc = NULL;
    double **phi = NULL;
    char **name = NULL;
    multifok_focus_op_basis_make(NW, ws, bType, ortho, &NB, &phi, &wc, &name);
    multifok_focus_op_basis_print(stderr, NW, NB, phi, wc, name);
    if (ortho) { multifok_focus_op_basis_ortho_check(stderr, NW, ws, NB, phi); }
    free(ws);
    multifok_focus_op_basis_free(NB, phi, wc, name);
  } 

#define multifok_focus_op_C_COPYRIGHT \
    "© 2017 by the State University of Campinas (UNICAMP)"

