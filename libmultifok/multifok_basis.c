/* See {multifok_basis.h}. */
/* Last edited on 2023-01-30 13:00:03 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <vec.h>
#include <jsmath.h>
#include <wt_table.h>
#include <affirm.h>
#include <bool.h>
#include <i2.h>

#include <multifok_window.h>

#include <multifok_basis.h>

/* NON-NORMALIZED BASES 

  The procedures below generate raw bases {bas[0..NB-1][0..NS-1]} of various
  types. The linops  Each procedure defines the number {NB} of elements in
  the basis, and fills {bas[0..NB-1][0..NS-1]} and {belName[0..NB-1]}. 
  
  Each procedure assumes that arrays {bas[]} and {belName[]} have been
  allocated with {NS=NW*NW} elements each, and returns as result the
  number {NB} of elements actually used.
*/

int32_t multifok_basis_fill_CANC(int32_t NW, double *bas[], char *belName[]);
  /* The window size {NW} must be odd.  Sets the basis size {NB} to {NS}, and
    {bas[0..NB-1]} to the canonical basis elements. That is, sets
    {bas[kb][ks]} to 1 if {kb==ks}, and 0 otherwise. 
    
    The element names  will be "S{X}{Y}" where {X} and {Y} are as in
    {multifok_window_sample_names}. */

int32_t multifok_basis_fill_LAPL_CUBE_DIFF(int32_t NW, multifok_basis_type_t bType, double *bas[], char *belName[]);
  /* Requires {NW=3}. 
  
    If {bType} is {multifok_basis_type_DIFF}, sets {NB} to 8, and
    {bas[0..NB-1][0..NS-1]} to the sample coefficients used to compute
    the derivative operators "DX", "DY", "DXX", "DXY", "DYY",
    "S3" (3-saddle), "C3" (3-cosaddle), and "Q" (checkerboard).
  
    If {bType} is {multifok_basis_type_CUBE}, sets {NB} to 7,
    and fills the basis with the same elements above, but
    omitting the "Q" (checkerboard) element.
  
    If {bType} is {multifok_basis_type_LAPL}, sets {NB} to 5,
    omitting the elements "S3", "C3", and "Q". */
  
int32_t multifok_basis_fill_HART(int32_t NW, double *bas[], char *belName[]);
  /* Requires {NW} odd. Sets {*NB_P} to {NB=NS}, and sets
    {bas[0..NB-1][0..NS-1]} to the Hartley waves with absolute {X} and
    {Y} frequencies up tp {NW/2}. The element names will be "H{X}{Y}"
    where {X} and {Y} are as in {multifok_window_sample_names}, only
    referring to frequencies instead of sample indices. */

void multifok_basis_set_elem_and_name_3x3
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

int32_t multifok_basis_orthize(int32_t NW, int32_t NB, double *bas[], char *belName[]);
  /* Makes the elements {bas[0..NB-1]} orthonormal with respect to the dot product
    {mulifok_window_prod}.
    
    The procedure may conclude that some basis elements are nearly
    dependent on others, and exclude them from the normalized basis. In
    that case, the number of non-discarded elements {NK} is returned as
    result, and those elements are rearranged to be {bas[0..NK-1]}.
    
    The procedure also adjusts the names {belName[0..NB-1]} to account for
    rearrangements of the respective elements, but not for the
    orthonormalization. On output the name {belName[kb]} still corresponds
    to {bas[kb]} for {kb} in {0..NK-1}. */
    
void multifok_basis_normalize_prod_range(int32_t NW, int32_t NB, double *bas[]);
  /* Scales each basis element {bas[0..NB-1]} so that its dot product with any array {s[0..NS-1]}
    of window samples is always in the range {[-1 _ +1]}, assuming that the window samples
    are in that range. */

void multifok_basis_make
  ( multifok_basis_type_t bType,
    int32_t NW, 
    double ws[],
    bool_t ortho,
    int32_t *NB_P, 
    double ***bas_P,
    char ***belName_P
  )
  { int32_t NS = multifok_window_num_samples(NW);
    
    int32_t NB_max = NS;  /* May be reduced later. */
    double **bas = notnull(malloc(NB_max*sizeof(double*)), "no mem"); 
    for (int32_t i = 0; i < NB_max; i++) { bas[i] = notnull(malloc(NS*sizeof(double)), "no mem"); }
    char **belName = notnull(malloc(NB_max*sizeof(char*)), "no mem"); 
    
    int32_t NB; /* Number of basis elements actually created. */
    
    if (bType == multifok_basis_type_CANC)
      { NB = multifok_basis_fill_CANC(NW, bas, belName); }
    else if (bType == multifok_basis_type_LAPL)
      { NB = multifok_basis_fill_LAPL_CUBE_DIFF(NW, bType, bas, belName); }
    else if (bType == multifok_basis_type_CUBE)
      { NB = multifok_basis_fill_LAPL_CUBE_DIFF(NW, bType, bas, belName); }
    else if (bType == multifok_basis_type_DIFF)
      { NB = multifok_basis_fill_LAPL_CUBE_DIFF(NW, bType, bas, belName); }
    else if (bType == multifok_basis_type_HART)
      { NB = multifok_basis_fill_HART(NW, bas, belName); }
    else
      { assert(FALSE); }

    if (ws != NULL)
      { /* Apodize basis elements: */
        for (int32_t kb = 0; kb < NB; kb++) 
          { double *bask = bas[kb];
            for (int32_t ks = 0; ks < NS; ks++) 
              { bask[ks] *= ws[ks]; }
          }
      }

    if (ortho)
      { /* Orthonormalize basis. */
        NB = multifok_basis_orthize(NW, NB, bas, belName);
        /* Free unused space: */
        if (NB < NB_max)
          { /* Free unused elements: */
            for (int32_t kb = NB; kb < NB_max; kb++) { free(bas[kb]); }
            bas = realloc(bas, NB*sizeof(double*));
            belName = realloc(belName, NB*sizeof(char*));
          }
      }
    
    (*NB_P) = NB;
    (*bas_P) = bas;
    (*belName_P) = belName;
  }
  
#define setbas multifok_basis_set_elem_and_name_3x3

int32_t multifok_basis_fill_CANC(int32_t NW, double *bas[], char *belName[])
  {
    demand(NW == 3, "not implemented for this {NW}");
    
    int32_t NS = NW*NW;

    int32_t NB = 0; /* Actual basis size. */

    setbas(bas,belName, NB, "Soo", 00.0, 00.0, 00.0,  00.0, +1.0, 00.0,  00.0, 00.0, 00.0 ); NB++; 

    setbas(bas,belName, NB, "Smo", 00.0, 00.0, 00.0,  +1.0, 00.0, 00.0,  00.0, 00.0, 00.0 ); NB++; 
    setbas(bas,belName, NB, "Spo", 00.0, 00.0, 00.0,  00.0, 00.0, +1.0,  00.0, 00.0, 00.0 ); NB++; 
    setbas(bas,belName, NB, "Som", 00.0, +1.0, 00.0,  00.0, 00.0, 00.0,  00.0, 00.0, 00.0 ); NB++; 
    setbas(bas,belName, NB, "Sop", 00.0, 00.0, 00.0,  00.0, 00.0, 00.0,  00.0, +1.0, 00.0 ); NB++; 

    setbas(bas,belName, NB, "Smm", +1.0, 00.0, 00.0,  00.0, 00.0, 00.0,  00.0, 00.0, 00.0 ); NB++; 
    setbas(bas,belName, NB, "Spp", 00.0, 00.0, 00.0,  00.0, 00.0, 00.0,  00.0, 00.0, +1.0 ); NB++; 
    setbas(bas,belName, NB, "Smp", 00.0, 00.0, 00.0,  00.0, 00.0, 00.0,  +1.0, 00.0, 00.0 ); NB++; 
    setbas(bas,belName, NB, "Spm", 00.0, 00.0, +1.0,  00.0, 00.0, 00.0,  00.0, 00.0, 00.0 ); NB++; 
    
    assert(NB == NS);
    
    return NB;
  }

int32_t multifok_basis_fill_LAPL_CUBE_DIFF(int32_t NW, multifok_basis_type_t bType, double *bas[], char *belName[])
  {
    demand(NW == 3, "{DIFF} basis requires {NW=3}");

    int32_t NB = 0;
    
    setbas(bas,belName, NB, "DX", -1.0, 00.0, +1.0,  -1.0, 00.0, +1.0,  -1.0, 00.0, +1.0 ); NB++; /* d/dX. */
    setbas(bas,belName, NB, "DY", -1.0, -1.0, -1.0,  00.0, 00.0, 00.0,  +1.0, +1.0, +1.0 ); NB++; /* d/dY. */
    setbas(bas,belName, NB, "DXX",-1.0, +1.0, -1.0,  -1.0, +1.0, -1.0,  -1.0, +1.0, -1.0 ); NB++; /* d^2/dX^2. */
    setbas(bas,belName, NB, "DYY",-1.0, -1.0, -1.0,  +1.0, +1.0, +1.0,  -1.0, -1.0, -1.0 ); NB++; /* d^2/dY^2. */
    setbas(bas,belName, NB, "DXY",-1.0, 00.0, +1.0,  00.0, 00.0, 00.0,  +1.0, 00.0, -1.0 ); NB++; /* d^2/dXdY. */
    
    if (bType == multifok_basis_type_LAPL) { assert(NB == 5); return NB; }
    
    setbas(bas,belName, NB, "S3", -1.0, +1.0, -1.0,  00.0, 00.0, 00.0,  +1.0, -1.0, +1.0 ); NB++; /* Saddle3. */
    setbas(bas,belName, NB, "C3", +1.0, 00.0, -1.0,  -1.0, 00.0, +1.0,  +1.0, 00.0, -1.0 ); NB++; /* Cosaddle3. */
    
    if (bType == multifok_basis_type_CUBE) { assert(NB == 7); return NB; }
    
    setbas(bas,belName, NB, "Q",  +1.0, -1.0, +1.0,  -1.0, +1.0, -1.0,  +1.0, -1.0, +1.0 ); NB++; /* Checker. */
    
    if (bType == multifok_basis_type_CUBE) { assert(NB == 8); return NB; }
    
    demand(FALSE, "invalid {bType}");
    return -1;
  }

int32_t multifok_basis_fill_HART(int32_t NW, double *bas[], char *belName[])
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
        belName[kb] = multifok_window_sample_name("H", fx, fy);
      }
      
    return NB;
  }

void multifok_basis_free(int32_t NB, double **bas, char **belName)
  {
    for (int32_t kb = 0; kb < NB; kb++) { free(bas[kb]); free(belName[kb]); }
    free(bas);
    free(belName);
  }

void multifok_basis_compute_coeffs(int32_t NW, double x[], int32_t NB, double *bas[], double coeff[])
  { for (int32_t kb = 0; kb  < NB; kb++) 
      { coeff[kb] = multifok_window_prod(NW, x, bas[kb]); }
  }

void multifok_basis_set_elem_and_name_3x3
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

int32_t multifok_basis_orthize(int32_t NW, int32_t NB, double *bas[], char *belName[])
  { 
    bool_t verbose = TRUE;
    double tiny = 1.0e-8; /* Insignificant relative length. */
    int32_t NS = multifok_window_num_samples(NW);
    int32_t NK = 0; /* Nymber of basis elements that were kept. */
    for (int32_t kb = 0; kb < NB; kb++)
      { 
        /* Original norm of element {bas[kb]}: */
        double d_old = sqrt(multifok_window_prod(NW, bas[kb], bas[kb]));
        
        /* Make element {bas[kb]} orthogonal to the previous elements:*/
        for (int32_t rb = 0; rb < NK; rb++)
          { double ci = multifok_window_prod(NW, bas[kb], bas[rb]);
            if (verbose && (ci > tiny*d_old))
              { fprintf(stderr, "  subtracting {%+12.8f*bas[%d]} from {bas[%d]}\n", ci, rb, kb); }
            for (int32_t js = 0; js < NS; js++)
              { bas[kb][js] = bas[kb][js] - ci*bas[rb][js]; }
          }

        /* Norm of residual element {bas[kb]} after removing previous elements: */
        double d_new = sqrt(multifok_window_prod(NW, bas[kb], bas[kb]));
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

void multifok_basis_print
  ( FILE *wr,
    int32_t NW, 
    int32_t NB, 
    double *bas[],
    char *belName[]
  )
  { int32_t NS = multifok_window_num_samples(NW);
    int32_t HW = (NW-1)/2;
    
    fprintf(wr, "--- basis ---------------------------------------------\n");
    
    /* Column numbers (sample indices): */
    fprintf(wr, "%3s %-12s   ", "", "");
    for (int32_t j = 0; j < NS; j++) { fprintf(wr, " %10d", j); }
    fprintf(wr, "\n");
    
    /* Column names (sample names): */
    fprintf(wr, "%3s %-12s   ", "", "");
    for (int32_t j = 0; j < NS; j++)
      { int32_t ix = (j % NW) - HW;
        int32_t iy = (j / NW) - HW;
        char *sname = multifok_window_sample_name("S", ix, iy);
        fprintf(wr, " %10s", sname);
      }
    fprintf(wr, "\n");
    
    for (int32_t ib = 0; ib < NB; ib++)
      { /* Find the smallest significant abs value in the row: */
        double bi_min = +INF;
        for (int32_t js = 0; js < NS; js++) 
          { double bij = fabs(bas[ib][js]);
            if ((bij > 0.002) && (bij < bi_min)) { bi_min = bij; }
          }
        /* Row (element) index and name: */
        fprintf(wr, "%3d %-12s [ ", ib, belName[ib]);
        for (int32_t js = 0; js < NS; js++)
          { fprintf(wr, " %+10.6f", bas[ib][js]/bi_min); }
        /* Scale factor for row: */
        fprintf(wr, " ] * %12.8f", bi_min);
        fprintf(wr, "\n");
      }
    fprintf(wr, "-------------------------------------------------------\n");
  }

void multifok_basis_ortho_check(FILE *wr, int32_t NW, int32_t NB, double *bas[])
  { /* Print basis products: */
    fprintf(wr, "--- momentum matrix -----------------------------------\n");
    fprintf(wr, "%3s", "");
    for (int32_t jb = 0; jb < NB; jb++) { fprintf(wr, " %10d", jb); }
    fprintf(wr, "\n");
    for (int32_t ib = 0; ib < NB; ib++)
      { fprintf(wr, "%3d", ib);
        for (int32_t jb = 0; jb <= ib; jb++)
          { double pij = multifok_window_prod(NW, bas[ib], bas[jb]);
            fprintf(wr, " %+10.6f", pij);
          }
        fprintf(wr, "\n");
      }
    fprintf(wr, "-------------------------------------------------------\n");
  }
    
void multifok_basis_check(int32_t NW, multifok_basis_type_t bType, bool_t ortho)
  { double *ws = multifok_window_sample_weights(NW);
    
    int32_t NB;
    double **bas = NULL;
    char **belName = NULL;
    multifok_basis_make(bType, NW, ws, ortho, &NB, &bas, &belName);
    multifok_basis_print(stderr, NW, NB, bas, belName);
    if (ortho) { multifok_basis_ortho_check(stderr, NW, NB, bas); }
    free(ws);
    multifok_basis_free(NB, bas, belName);
  } 

char *multifok_basis_type_to_text(multifok_basis_type_t bType)
  { switch(bType)
      { case multifok_basis_type_LAPL: return "LAPL";
        case multifok_basis_type_CUBE: return "CUBE";
        case multifok_basis_type_DIFF: return "DIFF";
        case multifok_basis_type_HART: return "HART";
        case multifok_basis_type_CANC: return "CANC";
        default: assert(FALSE);
      }
  }
    
multifok_basis_type_t multifok_basis_type_from_text(char *name, bool_t fail)
  { if (strcmp(name, "LAPL") == 0) { return multifok_basis_type_LAPL; }
    if (strcmp(name, "CUBE") == 0) { return multifok_basis_type_CUBE; }
    if (strcmp(name, "DIFF") == 0) { return multifok_basis_type_DIFF; }
    if (strcmp(name, "HART") == 0) { return multifok_basis_type_HART; }
    if (strcmp(name, "CANC") == 0) { return multifok_basis_type_CANC; }
    if (fail) { demand(FALSE, "invalid basis name"); } else { return -1; }
  }

void multifok_basis_write_elem_names(FILE *wr, int32_t NB, char *belName[])  
  { for (int32_t kb = 0; kb < NB; kb++)
      { fprintf(wr, "%s\n", belName[kb]); }
    fflush(wr);
  }

#define multifok_basis_C_COPYRIGHT \
    "© 2022 by the State University of Campinas (UNICAMP)"

