/* See {multifok_focus_op.h}. */
/* Last edited on 2023-01-25 18:08:26 by stolfi */

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

void multifok_focus_op_sample_names(int32_t NW, char *tag, char *sname[]);
  /* Sets {sname[0..NS-1]} to {NS = NW*NW} newly allocated strings that are scrutable names of
    the samples in the window. The names have the form "{tag}{X}{Y}" where  {X} and {Y} are the indices {ix}
    and {iy} of the sample relative to the center sample, encoded with {multifok_focus_op_mop_code}. */

/* NON-NORMALIZED BASES 

  The procedures below generate bases {bas[0..NB-1][0..NS-1]} of various
  types, but ignoring the sample weights (so the elements are not
  orthonormal). Each procedure defines the number {NB} of elements in
  the basis, and fills {bas[0..NB-1][0..NS-1]} and {belName[0..NB-1]}. 
  
  Each procedure assumes that arrays {bas[]} and {belName[]} have been
  allocated with {NS=NW*NW} elements each, and returns as result the
  number {NB} of elements actually used.
*/

int32_t multifok_focus_op_basis_fill_CANC(int32_t NW, double *bas[], char *belName[]);
  /* The window size {NW} must be odd. Sets {*NB_P} to and {NB=NS}, and
    sets {bas[0..NB-1]} to the canonical basis elements. That is, sets
    {bas[kb][ks]} to 1 if {kb==ks}, and 0 otherwise. The element names
    will be "P{X}{Y}" where {X} and {Y} are as in
    {multifok_focus_op_sample_names}. */

int32_t multifok_focus_op_basis_fill_LAPL(int32_t NW, double *bas[], char *belName[]);
  /* Requires {NW=3}. Sets {*NB_P} to {NB=5} and {bas[0..NB-1][0..NS-1]}
    to the sample coefficients used to compute the derivative operators
    "DX", "DY", "DXX", "DXY", and "DYY". */

int32_t multifok_focus_op_basis_fill_CUBE(int32_t NW, double *bas[], char *belName[]);
  /* Requires {NW=3}. Sets {*NB_P} to {NB=7} and
    {bas[0..NB-1][0..NS-1]} to the sample coefficients used to compute
    the derivative operators "DX", "DY", "DXX", "DXY", "DYY", "C3", and
    "S3" (3-saddles). */

int32_t multifok_focus_op_basis_fill_DIFF(int32_t NW, double *bas[], char *belName[]);
  /* Requires {NW=3}. Sets {*NB_P} to {NB=8} and
    {bas[0..NB-1][0..NS-1]} to the sample coefficients used to compute
    the derivative operators "DX", "DY", "DXX", "DXY", "DYY", "C3",
    "S3" (3-saddles) and "Q" (corner checker). */
  
int32_t multifok_focus_op_basis_fill_HART(int32_t NW, double *bas[], char *belName[]);
  /* Requires {NW} odd. Sets {*NB_P} to {NB=NS}, and sets
    {bas[0..NB-1][0..NS-1]} to the Hartley waves with absolute {X} and
    {Y} frequencies up tp {NW/2}. The element names will be "H{X}{Y}"
    where {X} and {Y} are as in {multifok_focus_op_sample_names}, only
    referring to frequencies instead of sample indices. */

void multifok_focus_op_set_basis_elem_and_name_3x3
  ( double *bas[], 
    char *belName[],
    int32_t ib, 
    char *na,
    double x00, double x01, double x02, 
    double x10, double x11, double x12, 
    double x20, double x21, double x22
  );
  /* Assumes the window is 3x3. Sets {bas[ib][0..NS-1]} to {{x00,x01,.. x22}}, and {belName[ib]} to 
    a fresh heap-allocated copy of {na}. */ 
  
/* BASIS NORMALIZATION */

int32_t multifok_focus_op_orthize(int32_t NW, double ws[], int32_t NB, double *bas[], char *belName[]);
  /* Makes the elements {bas[0..NB-1]} orthonormal with respect to the dot product with weights {ws[0..NS-1]} where
    {NS=NW*NW}.
    
    The procedure may conclude that some basis elements are nearly
    dependent on others, and exclude them from the normalized basis. In
    that case, the number of non-discarded elements {NK} is returned as
    result, and those elements are rearranged to be {bas[0..NK-1]}.
    
    The procedure also adjusts the names {belName[0..NB-1]} to account for
    rearrangements of the respective elements, but not for the
    orthonormalization. On output the name {belName[k]} still corresponds
    to {bas[k]} for {k} in {0..NK-1}. */
    
void multifok_focus_op_normalize_prod_range(int32_t NW, double ws[], int32_t NB, double *bas[]);
  /* Scales each basis element {bas[0..NB-1]} so that its dot product with any window in {[0_1]^NS} is 
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

void multifok_focus_op_basis_make
  ( int32_t NW, 
    double ws[],
    multifok_focus_op_basis_type_t bType,
    bool_t ortho,
    int32_t *NB_P, 
    double ***bas_P,
    char ***belName_P
  )
  { int32_t NS = multifok_focus_op_num_samples(NW);
    
    int32_t NB_max = NS;  /* May be reduced later. */
    double **bas = notnull(malloc(NB_max*sizeof(double*)), "no mem"); 
    for (int32_t i = 0; i < NB_max; i++) { bas[i] = notnull(malloc(NS*sizeof(double)), "no mem"); }
    char **belName = notnull(malloc(NB_max*sizeof(char*)), "no mem"); 
    
    int32_t NB; /* Number of basis elements actually created. */
    
    if (bType == multifok_focus_op_basis_type_CANC)
      { NB = multifok_focus_op_basis_fill_CANC(NW, bas, belName); }
    else if (bType == multifok_focus_op_basis_type_LAPL)
      { NB = multifok_focus_op_basis_fill_LAPL(NW, bas, belName); }
    else if (bType == multifok_focus_op_basis_type_CUBE)
      { NB = multifok_focus_op_basis_fill_CUBE(NW, bas, belName); }
    else if (bType == multifok_focus_op_basis_type_DIFF)
      { NB = multifok_focus_op_basis_fill_DIFF(NW, bas, belName); }
    else if (bType == multifok_focus_op_basis_type_HART)
      { NB = multifok_focus_op_basis_fill_HART(NW, bas, belName); }
    else
      { assert(FALSE); }

    if (ortho)
      { /* Orthonormalize basis. */
        NB = multifok_focus_op_orthize(NW, ws, NB, bas, belName);
        /* Free unused space: */
        if (NB < NB_max)
          { /* Free unused elements: */
            for (int32_t k = NB; k < NB_max; k++) { free(bas[k]); }
            bas = realloc(bas, NB*sizeof(double*));
            belName = realloc(belName, NB*sizeof(char*));
          }
      }
    else
      { /* Scale each element so that the dot product with any window in {[0_1]^NS} is always in {[-1 _ +1]}: */
        multifok_focus_op_normalize_prod_range(NW, ws, NB, bas);
      }
    
    (*NB_P) = NB;
    (*bas_P) = bas;
    (*belName_P) = belName;
  }
  
#define setbas multifok_focus_op_set_basis_elem_and_name_3x3

int32_t multifok_focus_op_basis_fill_CANC(int32_t NW, double *bas[], char *belName[])
  {
    int32_t NS = NW*NW;
    int32_t NB = NS; /* Actual basis size. */

    for (int32_t kb = 0; kb < NB; kb++)
      { for (int32_t ks = 0; ks < NS; ks++)
          { bas[kb][ks] = (kb == ks ? 1.0 : 0.0); }
      }
    multifok_focus_op_sample_names(NW, "P", belName);
    
    return NB;
  }

int32_t multifok_focus_op_basis_fill_LAPL(int32_t NW, double *bas[], char *belName[])
  {
    demand(NW == 3, "{LAPL} basis requires {NW=3}");
    int32_t NB = 5; /* Actual basis size. */
    assert(NB <= NW*NW);

    int32_t kb = 0;

    setbas(bas,belName, kb, "DX", -1.0, 00.0, +1.0,  -1.0, 00.0, +1.0,  -1.0, 00.0, +1.0f ); kb++;  /* d/dX. */
    setbas(bas,belName, kb, "DY", -1.0, -1.0, -1.0,  00.0, 00.0, 00.0,  +1.0, +1.0, +1.0f ); kb++;  /* d/dY. */
    setbas(bas,belName, kb, "DXX",-1.0, +1.0, -1.0,  -1.0, +1.0, -1.0,  -1.0, +1.0, -1.0f ); kb++;  /* d^2/dX^2. */
    setbas(bas,belName, kb, "DYY",-1.0, -1.0, -1.0,  +1.0, +1.0, +1.0,  -1.0, -1.0, -1.0f ); kb++;  /* d^2/dY^2. */
    setbas(bas,belName, kb, "DXY",-1.0, 00.0, +1.0,  00.0, 00.0, 00.0,  +1.0, 00.0, -1.0f ); kb++;  /* d^2/dXdY. */

    assert(kb == NB);
    
    return NB;
  }
   
int32_t multifok_focus_op_basis_fill_CUBE(int32_t NW, double *bas[], char *belName[])
  {
    demand(NW == 3, "{CUBE} basis requires {NW=3}");
    int32_t NB = 7;
    assert(NB <= NW*NW);

    int32_t kb = 0;
    
    setbas(bas,belName, kb, "DX", -1.0, 00.0, +1.0,  -1.0, 00.0, +1.0,  -1.0, 00.0, +1.0f ); kb++; /* d/dX. */
    setbas(bas,belName, kb, "DY", -1.0, -1.0, -1.0,  00.0, 00.0, 00.0,  +1.0, +1.0, +1.0f ); kb++; /* d/dY. */
    setbas(bas,belName, kb, "DXX",-1.0, +1.0, -1.0,  -1.0, +1.0, -1.0,  -1.0, +1.0, -1.0f ); kb++; /* d^2/dX^2. */
    setbas(bas,belName, kb, "DYY",-1.0, -1.0, -1.0,  +1.0, +1.0, +1.0,  -1.0, -1.0, -1.0f ); kb++; /* d^2/dY^2. */
    setbas(bas,belName, kb, "DXY",-1.0, 00.0, +1.0,  00.0, 00.0, 00.0,  +1.0, 00.0, -1.0f ); kb++; /* d^2/dXdY. */
    setbas(bas,belName, kb, "S3", -1.0, +1.0, -1.0,  00.0, 00.0, 00.0,  +1.0, -1.0, +1.0f ); kb++; /* Saddle3. */
    setbas(bas,belName, kb, "C3", +1.0, 00.0, -1.0,  -1.0, 00.0, +1.0,  +1.0, 00.0, -1.0f ); kb++; /* Cosaddle3. */
    
    assert(kb == NB);
    
    return NB;
  }
 
int32_t multifok_focus_op_basis_fill_DIFF(int32_t NW, double *bas[], char *belName[])
  {
    demand(NW == 3, "{DIFF} basis requires {NW=3}");
    int32_t NB = 8;
    assert(NB <= NW*NW);

    int32_t kb = 0;
    
    setbas(bas,belName, kb, "DX", -1.0, 00.0, +1.0,  -1.0, 00.0, +1.0,  -1.0, 00.0, +1.0f ); kb++; /* d/dX. */
    setbas(bas,belName, kb, "DY", -1.0, -1.0, -1.0,  00.0, 00.0, 00.0,  +1.0, +1.0, +1.0f ); kb++; /* d/dY. */
    setbas(bas,belName, kb, "DXX",-1.0, +1.0, -1.0,  -1.0, +1.0, -1.0,  -1.0, +1.0, -1.0f ); kb++; /* d^2/dX^2. */
    setbas(bas,belName, kb, "DYY",-1.0, -1.0, -1.0,  +1.0, +1.0, +1.0,  -1.0, -1.0, -1.0f ); kb++; /* d^2/dY^2. */
    setbas(bas,belName, kb, "DXY",-1.0, 00.0, +1.0,  00.0, 00.0, 00.0,  +1.0, 00.0, -1.0f ); kb++; /* d^2/dXdY. */
    setbas(bas,belName, kb, "S3", -1.0, +1.0, -1.0,  00.0, 00.0, 00.0,  +1.0, -1.0, +1.0f ); kb++; /* Saddle3. */
    setbas(bas,belName, kb, "C3", +1.0, 00.0, -1.0,  -1.0, 00.0, +1.0,  +1.0, 00.0, -1.0f ); kb++; /* Cosaddle3. */
    setbas(bas,belName, kb, "Q",  +1.0, -1.0, +1.0,  -1.0, +1.0, -1.0,  +1.0, -1.0, +1.0f ); kb++; /* Checker. */
    
    assert(kb == NB);
    
    return NB;
  }

int32_t multifok_focus_op_basis_fill_HART(int32_t NW, double *bas[], char *belName[])
  {
    assert((NW % 2 ) == 1);
    int32_t HW = NW/2;
    int32_t NS = NW*NW;
    int32_t NB = NS;
    
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
                bas[kb][ks] = f;
              }
          }
        belName[kb] = multifok_focus_op_sample_name("H", fx, fy);
      }
      
    return NB;
  }

void multifok_focus_op_basis_free(int32_t NB, double **bas, char **belName)
  {
    for (int32_t kb = 0; kb < NB; kb++) { free(bas[kb]); free(belName[kb]); }
    free(bas);
    free(belName);
  }

void multifok_focus_op_term_indices_from_names(int32_t NB, char *belName[], int32_t NT, char *tname[], i2_t tix[])
  {
    auto int32_t find_belName(char *beg, char *end);
      /* Returns the index {ib} such that {belName[ib]} is the string between 
        {*beg} (inclusive) and (*end) (exclusive). */

    for (int32_t kt = 0; kt < NT; kt++)
      {
        char *tnk = tname[kt];

        int32_t ib, jb;
        if (strcmp(tnk, "1") == 0)
          { /* Constant term: */
            ib = -1; jb = -1;
          }
        else 
          { char *ast = strchrnul(tnk, '*');
            ib = find_belName(tnk, ast);
            if ((*ast) == 0)
              { /* Linear term: */
                jb = -1;
              }
            else
              { char *end = strchr(tnk, 0);
                jb = find_belName(ast+1, end);
              }
          }
        tix[kt] = (i2_t){{ ib, jb }};
      }

    return;

    int32_t find_belName(char *beg, char *end)
      { int32_t n = (int32_t)(end - beg);
        int32_t ib = -1;
        for (int32_t kb = 0; kb < NB; kb++)
          { char *bk = belName[kb];
            int32_t m = (int32_t)strlen(bk);
            if (n != m) { continue; }
            if (strncmp(beg, bk, n) == 0) 
              { ib = kb; break; }
          }
        demand(ib >= 0, "basis name not recognized");
        return ib;
      }
  }


void multifok_focus_op_compute_basis_coeffs(int32_t NW, double x[], double ws[], int32_t NB, double *bas[], double coeff[])
  { for (int32_t kb = 0; kb  < NB; kb++) 
      { coeff[kb] = multifok_focus_op_prod(NW, x, bas[kb], ws); }
  }

double multifok_focus_op_score_from_basis_coeffs
  ( int32_t NB, 
    double coeff[],
    int32_t NT,
    i2_t tix[],
    double wt[],
    double term[], /* OUT */
    bool_t squared
  )
  { 
    double score = 0.0;
    if (NT > 0)
      { /* Compute and combine term values: */
        for (int32_t kt = 0; kt < NT; kt++) 
          { int32_t ib = tix[kt].c[0]; demand(ib < NB, "bad basis index {ib}");
            int32_t jb = tix[kt].c[1]; demand(jb < NB, "bad basis index {jb}");
            double tvk = NAN;
            if (ib < 0)
              { demand(jb < 0, "inconsitent {tix}");
                tvk = 1.0;
              }
            else if (jb < 0)
              { tvk = coeff[ib]; }
            else
              { tvk = coeff[ib]*coeff[jb]; }
            term[kt] = tvk;
            score += wt[kt]*tvk;
          }
      }
      
    /* Convert from squared to linear scale: */
    if (squared) 
      { score = fmax(0.0, score);
        score = sqrt(score);
      }
    
    return score;
  }

void multifok_focus_op_normalize_window_samples
  ( int32_t NW, 
    double x[], 
    double ws[], 
    double noise, 
    double *avg_P,
    double *dev_P
  )
  {
    int32_t NS = NW*NW;
    double sum_wx = 0;
    double sum_w = 0;
    for (int32_t ks = 0; ks < NS; ks++) 
      { sum_wx += ws[ks]*x[ks]; 
        sum_w += ws[ks];
      }
    assert(sum_w > 0);
    double avg = sum_wx/sum_w;
    double sum_wd2 = 0;
    for (int32_t ks = 0; ks < NS; ks++) 
      { double d = x[ks] - avg;
        sum_wd2 += ws[ks]*d*d;
      }
    double dev = sqrt(sum_wd2/sum_w);
    noise = fmax(1.0e-200, noise); /* To avoid division of zero by zero. */
    double mag = hypot(dev, noise);
    for (int32_t ks = 0; ks < NS; ks++) 
      { x[ks] = (x[ks] - avg) / mag; }
    (*avg_P) = avg;
    (*dev_P) = dev;
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

void multifok_focus_op_set_basis_elem_and_name_3x3
  ( double *bas[], 
    char *belName[],
    int32_t ib, 
    char *na,
    double x00, double x01, double x02, 
    double x10, double x11, double x12, 
    double x20, double x21, double x22
  )
  {
    belName[ib] = NULL; asprintf(&(belName[ib]), "%s", na); /* Make sure {belName[ib]} is a heap-allocated string. */
    double *pib = bas[ib];
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

void multifok_focus_op_sample_names(int32_t NW, char *tag, char *sname[])
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
            sname[ks] = uxy;
            ks++;
         }
      }
    assert(ks == NS);
    
    /* Release the auxiliary coordinate names: */
    for (int32_t i = -HW; i <= +HW; i++){ free(u[i + HW]); }
  }

void multifok_focus_op_normalize_prod_range(int32_t NW, double ws[], int32_t NB, double *bas[])
  {
    int32_t NS = NW*NW;
    for (int32_t kb = 0; kb < NB; kb++)
      { /* Determne the range {[sMin _ sMax]} of the dot product of {bas[kb]} and: */
        double sMax = 0.0; /* Sum of positive basis coordinates. */
        double sMin = 0.0; /* Sum of negative basis coordinates. */
        for (int32_t ks = 0; ks < NS; ks++) 
          { double wpk = ws[ks]*bas[kb][ks];
            if (wpk > 0.0) { sMax += wpk; }
            if (wpk < 0.0) { sMin += wpk; }
          }
        assert((sMin <= 0) && (sMax >= 0));
        /* Rescale basis samples to ensure the dot product is in the range {[-1 _ +1]}: */
        double scale = fmax(-sMin, sMax);
        assert(scale > 0.0);
        for (int32_t ks = 0; ks < NS; ks++) { bas[kb][ks] /= scale; }
      }
  }
 
int32_t multifok_focus_op_orthize(int32_t NW, double ws[], int32_t NB, double *bas[], char *belName[])
  { 
    bool_t verbose = TRUE;
    int32_t NS = multifok_focus_op_num_samples(NW);
    int32_t NK = 0; /* Nymber of basis elements that were kept. */
    for (int32_t kb = 0; kb < NB; kb++)
      { 
        /* Original norm of element {bas[kb]}: */
        double d_old = sqrt(multifok_focus_op_prod(NW, bas[kb], bas[kb], ws));
        
        /* Make element {bas[kb]} orthogonal to the previous elements:*/
        for (int32_t rb = 0; rb < NK; rb++)
          { double ci = multifok_focus_op_prod(NW, bas[kb], bas[rb], ws);
            for (int32_t js = 0; js < NS; js++)
              { bas[kb][js] = bas[kb][js] - ci*bas[rb][js]; }
          }

        /* Norm of residual element {bas[kb]} after removing previous elements: */
        double d_new = sqrt(multifok_focus_op_prod(NW, bas[kb], bas[kb], ws));
        if (verbose) { fprintf(stderr, "d old = %16.12f new = %16.12f\n", d_old, d_new); }
        if (d_new > 1.0e-8*d_old)
          { /* Residual seems significant: */
            for (int32_t js = 0; js < NS; js++)
              { bas[NK][js] = bas[kb][js]/d_new; }
            belName[NK] = belName[kb];
            NK++;
          }
        else
          { /* Element was mostly linear comb of previous ones, discard: */
            if (verbose) { fprintf(stderr, "!! eliminated element {bas[%d]}\n", kb); }
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

void multifok_focus_op_basis_print
  ( FILE *wr,
    int32_t NW, 
    int32_t NB, 
    double *bas[],
    char *belName[]
  )
  { int32_t NS = multifok_focus_op_num_samples(NW);
    
    /* Print basis: */
    fprintf(wr, "--- basis ---------------------------------------------\n");
    fprintf(wr, "%3s %-12s ", "", "");
    for (int32_t j = 0; j < NS; j++) { fprintf(wr, " %10d", j); }
    fprintf(wr, "\n");
    for (int32_t ib = 0; ib < NB; ib++)
      { fprintf(wr, "%3d %-12s ", ib, belName[ib]);
        for (int32_t js = 0; js < NS; js++)
          { double pij = bas[ib][js];
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
    double *bas[]
  )
  { /* Print basis products: */
    fprintf(wr, "--- momentum matrix -----------------------------------\n");
    fprintf(wr, "%3s", "");
    for (int32_t jb = 0; jb < NB; jb++) { fprintf(wr, " %10d", jb); }
    fprintf(wr, "\n");
    for (int32_t ib = 0; ib < NB; ib++)
      { fprintf(wr, "%3d", ib);
        for (int32_t jb = 0; jb <= ib; jb++)
          { double pij = multifok_focus_op_prod(NW, bas[ib], bas[jb], ws);
            fprintf(wr, " %+10.6f", pij);
          }
        fprintf(wr, "\n");
      }
    fprintf(wr, "-------------------------------------------------------\n");
  }
    
void multifok_focus_op_check(int32_t NW, multifok_focus_op_basis_type_t bType, bool_t ortho)
  { double *ws = multifok_focus_op_sample_weights(NW);
    
    int32_t NB;
    double **bas = NULL;
    char **belName = NULL;
    multifok_focus_op_basis_make(NW, ws, bType, ortho, &NB, &bas, &belName);
    multifok_focus_op_basis_print(stderr, NW, NB, bas, belName);
    if (ortho) { multifok_focus_op_basis_ortho_check(stderr, NW, ws, NB, bas); }
    free(ws);
    multifok_focus_op_basis_free(NB, bas, belName);
  } 

char *multifok_focus_op_basis_type_to_text(multifok_focus_op_basis_type_t bType)
  { switch(bType)
      { case multifok_focus_op_basis_type_LAPL: return "LAPL";
        case multifok_focus_op_basis_type_CUBE: return "CUBE";
        case multifok_focus_op_basis_type_DIFF: return "DIFF";
        case multifok_focus_op_basis_type_HART: return "HART";
        case multifok_focus_op_basis_type_CANC: return "CANC";
        default: assert(FALSE);
      }
  }
    
multifok_focus_op_basis_type_t multifok_focus_op_basis_type_from_text(char *name, bool_t fail)
  { if (strcmp(name, "LAPL") == 0) { return multifok_focus_op_basis_type_LAPL; }
    if (strcmp(name, "CUBE") == 0) { return multifok_focus_op_basis_type_CUBE; }
    if (strcmp(name, "DIFF") == 0) { return multifok_focus_op_basis_type_DIFF; }
    if (strcmp(name, "HART") == 0) { return multifok_focus_op_basis_type_HART; }
    if (strcmp(name, "CANC") == 0) { return multifok_focus_op_basis_type_CANC; }
    if (fail) { demand(FALSE, "invalid basis name"); } else { return -1; }
  }

#define multifok_focus_op_C_COPYRIGHT \
    "© 2022 by the State University of Campinas (UNICAMP)"

