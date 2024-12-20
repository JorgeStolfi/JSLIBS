/* See {multifok_basis.h}. */
/* Last edited on 2024-12-05 14:28:15 by stolfi */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <vec.h>
#include <jsmath.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <bool.h>
#include <i2.h>

#include <multifok_window.h>

#include <multifok_basis.h>

/* NON-NORMALIZED BASES 

  The procedures below assume that arrays {basis->bas[][]} and
  {basis->belName[]} have been allocated with {NS=NW*NW} elements
  each.  Each procedure uses only {NB <= NS} of those
  elements.  It sets {basis->NB}, {basis->bas[0..NB-1][0..NS-1]},
  and {basis->belName[0..NB-1]} as appropriate for specific basis types. */

void multifok_basis_fill_CANC(multifok_basis_t *basis);
  /* The window size {NW} must be odd.  Sets the basis size {NB} to {NS}, and
    {bas[0..NB-1]} to the canonical basis elements. That is, sets
    {bas[kb][ks]} to 1 if {kb==ks}, and 0 otherwise.  */

void multifok_basis_fill_LAPL_CUBE_DIFF(multifok_basis_type_t type, multifok_basis_t *basis);
  /* Requires {NW=3}.  Sets {basis->NB} and {basis->bas} for the "LAPL", "CUBE", or "DIFF" basis,
    according to the {type}. */
  
void multifok_basis_fill_HART(multifok_basis_t *basis);
  /* Requires {NW} odd. Sets {basis->NB} and {basis->bas} for the 
    for the Hartley basis. */

void multifok_basis_add_elem_and_name_3x3
  ( multifok_basis_t *basis, 
    char *na,
    double x00, double x01, double x02, 
    double x10, double x11, double x12, 
    double x20, double x21, double x22
  );
  /* Assumes the window is 3x3. Sets {bas[ib][0..NS-1]} to {{x00,x01,.. x22}}, and {belName[ib]} to 
    a fresh heap-allocated copy of {na}. */ 
  
/* BASIS NORMALIZATION */

void multifok_basis_orthize(multifok_basis_t *basis, bool_t verbose);
  /* Makes the elements {basis->bas[0..basis->NB-1]} orthonormal with
    respect to the dot product {multifok_window_prod} and window sample
    weights {basis->ws[]}.
    
    The procedure may conclude that some basis elements are nearly
    dependent on others, and exclude them from the normalized basis. In
    that case, the number {basis->NB} is reduced and those elements are
    compacted to be {basis->bas[0..basis->NB-1]}.
    
    The procedure also adjusts the names {belName[0..NB-1]} to account for
    rearrangements of the respective elements, but not for the
    orthonormalization. On output the name {belName[kb]} still corresponds
    to {bas[kb]} for {kb} in {0..NK-1}. */
    
void multifok_basis_normalize_prod_range(multifok_basis_t *basis);
  /* Scales each basis element {bas[0..NB-1]} so that its dot product
    with any array {s[0..NS-1]} of window samples is always in the range
    {[-1 _ +1]}, assuming that the window samples are in that range. */

multifok_basis_t *multifok_basis_make
  ( multifok_basis_type_t type,
    uint32_t NW, 
    double ws[],
    bool_t ortho,
    bool_t verbose
  )
  { demand((NW > 0) && ((NW % 2) == 1), "window size {NW} must be odd");
    
    uint32_t NS = multifok_window_num_samples(NW);
  
    multifok_basis_t *basis = talloc(1, multifok_basis_t);
    uint32_t NB_max = NS;  /* May be reduced later. */
    basis->bas = talloc(NB_max, double*);
    basis->belName = talloc(NB_max, char*); 
    for (uint32_t i = 0;  i < NB_max; i++) 
      { basis->bas[i] = talloc(NS, double);
        basis->belName [i] = NULL;
      }
    basis->NW = NW;
    basis->ws = ws;
   
    basis->NB = 0; /* Actual basis size. */
    
    if (type == multifok_basis_type_CANC)
      { multifok_basis_fill_CANC(basis);  }
    else if (type == multifok_basis_type_LAPL)
      { multifok_basis_fill_LAPL_CUBE_DIFF(type, basis); }
    else if (type == multifok_basis_type_CUBE)
      { multifok_basis_fill_LAPL_CUBE_DIFF(type, basis); }
    else if (type == multifok_basis_type_DIFF)
      { multifok_basis_fill_LAPL_CUBE_DIFF(type, basis); }
    else if (type == multifok_basis_type_HART)
      { multifok_basis_fill_HART(basis); }
    else
      { assert(FALSE); }

    basis->ortho = ortho;
    if (ortho) { multifok_basis_orthize(basis, verbose); }

    /* Free unused space: */
    if (basis->NB < NB_max)
      { /* Free unused elements: */
        for (uint32_t kb = basis->NB; kb < NB_max; kb++)
          { free(basis->bas[kb]); 
            if (basis->belName[kb] != NULL) { free(basis->belName[kb]); }            
          }
        basis->bas = retalloc(basis->bas, basis->NB, double*);
        basis->belName = retalloc(basis->belName, basis->NB, char*);
      }

    return basis;
  }
  
#define addbas multifok_basis_add_elem_and_name_3x3

void multifok_basis_fill_CANC(multifok_basis_t *basis)
  {
    demand(basis->NW == 3, "{CANC} basis requires {NW=3}");
    uint32_t NS = multifok_window_num_samples(basis->NW);
    basis->NB = 0;

    addbas(basis, "Soo", 00.0, 00.0, 00.0,  00.0, +1.0, 00.0,  00.0, 00.0, 00.0 ); 

    addbas(basis, "Smo", 00.0, 00.0, 00.0,  +1.0, 00.0, 00.0,  00.0, 00.0, 00.0 ); 
    addbas(basis, "Spo", 00.0, 00.0, 00.0,  00.0, 00.0, +1.0,  00.0, 00.0, 00.0 ); 
    addbas(basis, "Som", 00.0, +1.0, 00.0,  00.0, 00.0, 00.0,  00.0, 00.0, 00.0 ); 
    addbas(basis, "Sop", 00.0, 00.0, 00.0,  00.0, 00.0, 00.0,  00.0, +1.0, 00.0 ); 

    addbas(basis, "Smm", +1.0, 00.0, 00.0,  00.0, 00.0, 00.0,  00.0, 00.0, 00.0 ); 
    addbas(basis, "Spp", 00.0, 00.0, 00.0,  00.0, 00.0, 00.0,  00.0, 00.0, +1.0 ); 
    addbas(basis, "Smp", 00.0, 00.0, 00.0,  00.0, 00.0, 00.0,  +1.0, 00.0, 00.0 ); 
    addbas(basis, "Spm", 00.0, 00.0, +1.0,  00.0, 00.0, 00.0,  00.0, 00.0, 00.0 ); 
    
    assert(basis->NB == NS);
  }

void multifok_basis_fill_LAPL_CUBE_DIFF(multifok_basis_type_t type, multifok_basis_t *basis)
  {
    uint32_t NW = basis->NW;
    demand(NW == 3, "{LAPL}, {DIFF}, and {CUBE} bases requires {NW=3}");
    basis->NB = 0;
    
    addbas(basis, "F",  +1.0, +1.0, +1.0,  +1.0, +1.0, +1.0,  +1.0, +1.0, +1.0 ); /* F. */

    addbas(basis, "FX", -1.0, 00.0, +1.0,  -1.0, 00.0, +1.0,  -1.0, 00.0, +1.0 ); /* dF/dX. */
    addbas(basis, "FY", -1.0, -1.0, -1.0,  00.0, 00.0, 00.0,  +1.0, +1.0, +1.0 ); /* dF/dY. */
    addbas(basis, "FXX",+1.0, -2.0, +1.0,  +1.0, -2.0, +1.0,  +1.0, -2.0, +1.0 ); /* d^2F/dX^2. */
    addbas(basis, "FYY",+1.0, +1.0, +1.0,  -2.0, -2.0, -2.0,  +1.0, +1.0, +1.0 ); /* d^2F/dY^2. */
    addbas(basis, "FXY",-1.0, 00.0, +1.0,  00.0, 00.0, 00.0,  +1.0, 00.0, -1.0 ); /* d^2F/dXdY. */
    
    if (type == multifok_basis_type_LAPL) { assert(basis->NB == 6); return; }
    
    addbas(basis, "S3", -1.0, +2.0, -1.0,  00.0, 00.0, 00.0,  +1.0, -2.0, +1.0 ); /* Saddle3. */
    addbas(basis, "C3", +1.0, 00.0, -1.0,  -2.0, 00.0, +2.0,  +1.0, 00.0, -1.0 ); /* Cosaddle3. */
    
    if (type == multifok_basis_type_CUBE) { assert(basis->NB == 8); return; }
    
    addbas(basis, "Q",  +1.0, -2.0, +1.0,  -2.0, +4.0, -2.0,  +1.0, -2.0, +1.0 ); /* Checker. */
    
    if (type == multifok_basis_type_DIFF) { assert(basis->NB == 9); return; }
    
    demand(FALSE, "invalid {type}");
  }

void multifok_basis_fill_HART(multifok_basis_t *basis)
  {
    uint32_t NW = basis->NW;
    demand((NW % 2) == 1, "window size {NW} must be odd");
    uint32_t HW = (NW-1)/2;
    uint32_t NS = multifok_window_num_samples(NW);

    uint32_t NB = NS;
    
    /* Generate tables {fre[0..NB-1]} of frequency vectors: */
    i2_t freq[NB];
    uint32_t kb = 0; /* Basis element index. */
    for (int32_t fy = -(int32_t)HW; fy <= +(int32_t)HW; fy++)
      { for (int32_t fx = -(int32_t)HW; fx <= +(int32_t)HW; fx++)
          { freq[kb] = (i2_t){{ fx, fy}}; kb++; }
      }
    assert(kb == NB);
    
    /* Sort the table by increasing freq modulus: */
    for (uint32_t ib = 0;  ib < NB; ib++)
      { /* The {kb} smallest freqs are {freq[0..kb-1]} */
        uint32_t jbMin = ib;
        for (uint32_t jb = ib+1; jb < NB; jb++)
          { if (i2_norm_sqr(&(freq[jb])) < i2_norm_sqr(&(freq[jbMin])))
              { jbMin = jb; }
          }
        if (jbMin != ib)
          { i2_t ft = freq[ib]; freq[ib] = freq[jbMin]; freq[jbMin] = ft; }
      }
    
    /* Now fill the basis: */
    for (uint32_t kb = 0;  kb < NB; kb++)
      { 
        int32_t fx = freq[kb].c[0];
        int32_t fy = freq[kb].c[1];
        for (int32_t iy = -(int32_t)HW; iy <= +(int32_t)HW; iy++)
          { for (int32_t ix = -(int32_t)HW; ix <= +(int32_t)HW; ix++)
              { int32_t ks = (iy+(int32_t)(HW*NW)) + (ix+(int32_t)HW);
                double z = 0.125 + ((double)(fx*ix + fy*iy))/((double)NW); /* Hartley phase */
                double f = sin(2*M_PI*z);
                basis->bas[kb][ks] = f;
              }
          }
        basis->belName[kb] = multifok_window_sample_name("H", fx, fy);
      }
    basis->NB = NB;
    return;
  }

void multifok_basis_free(multifok_basis_t *basis)
  {
    for (uint32_t kb = 0;  kb < basis->NB; kb++) 
      { free(basis->bas[kb]); 
        if (basis->belName[kb] != NULL) { free(basis->belName[kb]); }
      }
    free(basis->bas);
    free(basis->belName);
    free(basis);
  }

void multifok_basis_compute_coeffs(double x[], multifok_basis_t *basis, double coeff[])
  { for (uint32_t kb = 0;  kb  < basis->NB; kb++) 
      { coeff[kb] = multifok_window_prod(basis->NW, x, basis->bas[kb]); }
  }

void multifok_basis_add_elem_and_name_3x3
  ( multifok_basis_t *basis,
    char *na,
    double x00, double x01, double x02, 
    double x10, double x11, double x12, 
    double x20, double x21, double x22
  )
  {
    uint32_t NW = basis->NW;
    demand(NW == 3, "window size should be 3x3");
    uint32_t NB_max = multifok_window_num_samples(NW);
    uint32_t ib = basis->NB;
    assert(ib < NB_max); /* Too many elemets. */
    /* Make sure {belName[ib]} is a heap-allocated copy of {na}: */
    basis->belName[ib] = jsprintf("%s", na); 
    double *pib = basis->bas[ib];
    pib[0] = x00;
    pib[1] = x01;
    pib[2] = x02;
    pib[3] = x10;
    pib[4] = x11;
    pib[5] = x12;
    pib[6] = x20;
    pib[7] = x21;
    pib[8] = x22;
    basis->NB++;
  }

void multifok_basis_orthize(multifok_basis_t *basis, bool_t verbose)
  { 
    uint32_t NW = basis->NW;
    uint32_t NS = multifok_window_num_samples(NW);
    double tiny = 1.0e-8; /* Insignificant relative length. */
    uint32_t NK = 0; /* Number of basis elements that were kept. */
    for (uint32_t kb = 0;  kb < basis->NB; kb++)
      { 
        /* Original norm of element {basis->bas[kb]}: */
        double d_old = sqrt(multifok_window_prod(NW, basis->bas[kb], basis->bas[kb]));
        
        /* Make element {basis->bas[kb]} orthogonal to the previous elements:*/
        for (uint32_t rb = 0;  rb < NK; rb++)
          { double ci = multifok_window_prod(NW, basis->bas[kb], basis->bas[rb]);
            if (verbose && (ci > tiny*d_old))
              { fprintf(stderr, "  subtracting {%+12.8f*bas[%d]} from {bas[%d]}\n", ci, rb, kb); }
            for (uint32_t js = 0;  js < NS; js++)
              { basis->bas[kb][js] = basis->bas[kb][js] - ci*basis->bas[rb][js]; }
          }

        /* Norm of residual element {bas[kb]} after removing previous elements: */
        double d_new = sqrt(multifok_window_prod(NW, basis->bas[kb], basis->bas[kb]));
        if (verbose) { fprintf(stderr, "d old = %16.12f new = %16.12f\n", d_old, d_new); }
        if (d_new > tiny*d_old)
          { /* Residual seems significant: */
            for (uint32_t js = 0;  js < NS; js++)
              { basis->bas[NK][js] = basis->bas[kb][js]/d_new; }
            basis->belName[NK] = basis->belName[kb];
            NK++;
          }
        else
          { /* Element was mostly linear comb of previous ones, discard: */
            if (verbose) 
              { fprintf(stderr, "!! eliminated element {bas[%d] = %s}\n", kb, basis->belName[kb]); }
          }
      }
    basis->NB = NK;
  }

void multifok_basis_print(FILE *wr, multifok_basis_t *basis)
  { uint32_t NW = basis->NW;
    uint32_t NS = multifok_window_num_samples(NW);
    uint32_t HW = (NW-1)/2;
    
    fprintf(wr, "--- basis ---------------------------------------------\n");
    
    /* Column numbers (sample indices): */
    fprintf(wr, "%3s %-12s   ", "", "");
    for (uint32_t j = 0;  j < NS; j++) { fprintf(wr, " %10d", j); }
    fprintf(wr, "\n");
    
    /* Column names (sample names): */
    fprintf(wr, "%3s %-12s   ", "", "");
    for (int32_t j = 0;  j < NS; j++)
      { int32_t ix = (j % (int32_t)NW) - (int32_t)HW;
        int32_t iy = (j / (int32_t)NW) - (int32_t)HW;
        char *sname = multifok_window_sample_name("S", ix, iy);
        fprintf(wr, " %10s", sname);
      }
    fprintf(wr, "\n");
    
    for (uint32_t ib = 0;  ib < basis->NB; ib++)
      { /* Find the smallest significant abs value in the row: */
        double bi_min = +INF;
        for (uint32_t js = 0;  js < NS; js++) 
          { double bij = fabs(basis->bas[ib][js]);
            if ((bij > 0.002) && (bij < bi_min)) { bi_min = bij; }
          }
        /* Row (element) index and name: */
        fprintf(wr, "%3d %-12s [ ", ib, basis->belName[ib]);
        for (uint32_t js = 0;  js < NS; js++)
          { fprintf(wr, " %+10.6f", basis->bas[ib][js]/bi_min); }
        /* Scale factor for row: */
        fprintf(wr, " ] * %12.8f", bi_min);
        fprintf(wr, "\n");
      }
    fprintf(wr, "-------------------------------------------------------\n");
  }

void multifok_basis_ortho_check(FILE *wr, multifok_basis_t *basis)
  { fprintf(wr, "--- moment matrix -----------------------------------\n");
    fprintf(wr, "%4s", "");
    for (uint32_t jb = 0;  jb < basis->NB; jb++) { fprintf(wr, " %10d", jb); }
    fprintf(wr, "\n");
    double orth_tol = 1.0e-7;  /* Tolerance for deviation from orthogonality. */
    double norm_tol = 1.0e-7;  /* Tolerance for deviation from normality. */
    double max_orth_err = 0.0; /* Max deviation from orthogonality. */
    double max_norm_err = 0.0; /* Max deviation from normality. */
    for (uint32_t ib = 0;  ib < basis->NB; ib++)
      { fprintf(wr, "%4d", ib);
        for (uint32_t jb = 0;  jb <= ib; jb++)
          { double pij = multifok_window_prod(basis->NW, basis->bas[ib], basis->bas[jb]);
            fprintf(wr, " %+10.6f", pij);
            if (ib == jb)
              { double norm_err = fabs(pij - 1.0);
                if (norm_err > max_norm_err) { max_norm_err = norm_err; }
              }
            else
              { double orth_err = fabs(pij);
                if (orth_err > max_orth_err) { max_orth_err = orth_err; }
              }
          }
        fprintf(wr, "\n");
      }
    char *msg_orth = (max_orth_err > orth_tol ? "(NOT orthogonal!)" : "(orthogonal)");
    char *msg_norm = (max_norm_err > norm_tol ? "(NOT normalized!)" : "(normalized)");
    fprintf(wr, "  max orthogonality error = %10.8f %s\n", max_orth_err, msg_orth);
    fprintf(wr, "  max normalization error = %10.8f %s\n", max_norm_err, msg_norm);
    fprintf(wr, "-------------------------------------------------------\n");
  }
    
void multifok_basis_module_check(uint32_t NW, multifok_basis_type_t type, bool_t ortho, bool_t verbose)
  { double *ws = multifok_window_weights_binomial(NW);
    
    multifok_basis_t *basis = multifok_basis_make(type, NW, ws, ortho, verbose);
    multifok_basis_print(stderr, basis);
    if (ortho) { multifok_basis_ortho_check(stderr, basis); }
    free(ws);
    multifok_basis_free(basis);
  } 

char *multifok_basis_type_to_text(multifok_basis_type_t type)
  { switch(type)
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

void multifok_basis_elem_names_write(FILE *wr, uint32_t NB, char *belName[])  
  { for (uint32_t kb = 0;  kb < NB; kb++)
      { fprintf(wr, "%s\n", belName[kb]); }
    fflush(wr);
  }

void multifok_basis_elem_names_write_named(char *outPrefix, multifok_basis_t *basis)  
  { char *fname = jsprintf("%s-bnames.txt", outPrefix);
    FILE *wr = open_write(fname, TRUE);
    multifok_basis_elem_names_write(wr, basis->NB, basis->belName);
    fclose(wr);
  }

#define multifok_basis_C_COPYRIGHT \
    "© 2022 by the State University of Campinas (UNICAMP)"

