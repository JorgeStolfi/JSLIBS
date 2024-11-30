/* See {multifok_term.h}. */
/* Last edited on 2024-10-11 21:11:51 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <vec.h>
#include <affirm.h>
#include <bool.h>
#include <fget.h>
#include <jsfile.h>

#include <multifok_window.h>
#include <multifok_basis.h>

#include <multifok_term.h>

void multifok_term_read_set_weights_and_names
  ( FILE *rd,
    bool_t weights,
    int32_t *NT_P,
    char ***termName_P,
    double *wt_P[],
    bool_t verbose
  );
  /* Reads from {rd} the weights {wt[0..NT-1]} (if {weights} is true)
    and the term names {termName[0..NT-1]}, as per {multifok_term_read_set}.
    
    Deduces {NT} from the number of data lines in the file, and returns
    it in {*NT_P}. The tables {wt} and {termName} are allocated by the
    procedure and returned in {*termName_P} and {*wt_P}. */
  
void multifok_term_indices_from_names
  ( multifok_basis_t *basis, 
    int32_t NT, 
    char *termName[], 
    int32_t *NP_P, 
    multifok_term_prod_t **prix_P, 
    bool_t verbose
  );
  /* The parameter {termName[0..NT-1]} must be a list of formulas of
    quadratic operator terms. Each formula must be a sum of products of
    pairs of coefficients, like "FX*FY+FXX*FYY". Each product in that
    sum must be "{el1}*{el2}" where {el1=belName[jb1]} and
    {el2=belName[jb2]} for some jb1,jb2} with {0 <= jb1,jb2 < NB}. As
    a special case, a product may be just "1".
  
    The procedure determines the total number {NP} of products appearing
    in the {termName[0..NT-1]}, such as "FX*FY" or "Soo*Sop", and
    creates a table {prix[0..NP-1]} that describes the products and the
    terms they belong to. See {multifok_term_values_from_basis_coeffs}
    for the meaning of {prix}.
    
    In each product, the two factors are internally swapped if needed so
    that {jb1 <= jb2}. After this adjustment, each product should appear
    only once. Thus {NP} is at most {NB*(NP+1)/2}. Also each term must
    have at least one product, so {NT} is at most equal to {NP}.
    
    The product count {NP} and the table {prix} are returned in {*NP_P}
    and {*prix_P}. The {pname} fields in the {prix} table will be newly
    allocated strings. */

void multifok_term_parse_product
  ( char *pbeg, 
    int32_t nc, 
    int32_t NB, 
    char *belName[], 
    int32_t *jb1_P, 
    int32_t *jb2_P,
    char ** pr_P
  );
  /* Parses the string between {pbeg} and {pbeg+nc} as "1" or the product "{el1}*{el2}"
    of two basis coeffs, where {el1 = belName[jb1]} and {el2 = belName[jb2]},
    for some {jb1,jb2} such that {0 <= jb1 <= jb2 < NB}. 
    
    Fails if the string does not have that form.  If it succeeds, returns {jb1,jb2} and a 
    new copy {pr} of that string in {*jb1_P,*jb2_P,*pr_P}.  */

void multifok_term_values_from_basis_coeffs
  ( double coeff[],
    multifok_term_set_t *tset,
    double term[]
  )
  {
    int32_t NT = tset->NT;
    int32_t NP = tset->NP;
    int32_t NB = tset->basis->NB;
    
    for (uint32_t kt = 0;  kt < NT; kt++) { term[kt] = 0; }
    for (uint32_t ip = 0;  ip < NP; ip++) 
      { int32_t jb1 = tset->prix[ip].jb1; 
        demand((jb1 >= 0) && (jb1 < NB), "bad basis index {jb1}");
        int32_t jb2 = tset->prix[ip].jb2; 
        demand((jb2 >= 0) && (jb2 < NB), "bad basis index {jb2}");
        double tvk;
        if (jb1 == -1)
          { demand(jb2 == -1, "inconsitent {jb1,jb1} for unit term");
            tvk = 1.0;
          }
        else
          { demand((jb1 >= 0) && (jb1 < NB), "invalid basis coeff index {jb1}"); 
            demand((jb2 >= 0) && (jb2 < NB), "invalid basis coeff index {jb2}"); 
            tvk = coeff[jb1]*coeff[jb2];
          }
        int32_t kt = tset->prix[ip].kt;
        term[kt] += tvk;
      }
  }

multifok_term_set_t *multifok_term_set_read
  ( FILE *rd,
    bool_t weights,
    multifok_basis_t *basis,
    double *wt_P[],
    bool_t verbose
  )
  {
    multifok_term_set_t *tset = talloc(1, multifok_term_set_t);
    
    tset->basis = basis;

    if (verbose) { fprintf(stderr, "reading %sterms...\n", (weights ? "weights and " : "")); }
    multifok_term_read_set_weights_and_names(rd, weights, &(tset->NT), &(tset->termName), wt_P, verbose);
    if (verbose) { fprintf(stderr, "parsing term names into product index table...\n"); }
    multifok_term_indices_from_names(basis, tset->NT, tset->termName, &(tset->NP), &(tset->prix), verbose);
    
    return tset;
  }
  
multifok_term_set_t *multifok_term_set_read_named
  ( char *fname,
    bool_t weights,
    multifok_basis_t *basis,
    double **wt_P,
    bool_t verbose
  )
  {
    FILE *rd = open_read(fname, TRUE);
    multifok_term_set_t *tset = multifok_term_set_read(rd, weights, basis, wt_P, verbose);
    fclose(rd);
    return tset;
  }
      
void multifok_term_read_set_weights_and_names
  ( FILE *rd,
    bool_t weights,
    int32_t *NT_P,
    char ***termName_P,
    double *wt_P[],
    bool_t verbose
  )
  {
    int32_t NT = 0; /* Number of terms seen in file. */
    double_vec_t wt = double_vec_new(50);
    string_vec_t termName = string_vec_new(50); /* Term formulas. */
    while (TRUE)
      { bool_t ok = fget_test_comment_or_eol(rd, '#', NULL);
        if (ok) { continue; }
        if (fget_test_eof(rd)) { break; }

        /* There is something there. Read the term index: */
        int32_t kt = fget_int32(rd);
        demand(kt == NT, "unexpected term index");

        /* Read the term weight, if any: */
        double wtk = 1.0;
        if (weights) { wtk = fget_double(rd); }
        double_vec_expand(&wt, NT);
        wt.e[kt] = wtk;
        
        /* Read the term name, if any: */
        char *tnamek = fget_string(rd);
        string_vec_expand(&termName, NT);
        termName.e[NT] = tnamek;
        
        if (verbose) { fprintf(stderr, "  %3d %12.8f %s\n", NT, wtk, (tnamek == NULL ? "" : tnamek)); }
        fget_comment_or_eol(rd, '#', NULL);
        NT++;
      }
    string_vec_trim(&termName, NT);
    double_vec_trim(&wt, NT);
    
    (*NT_P) = NT;
    (*termName_P) = termName.e;
    if (wt_P == NULL) { free(wt.e); } else { (*wt_P) = wt.e; }
  }

void multifok_term_indices_from_names
  ( multifok_basis_t *basis,
    int32_t NT, 
    char *termName[], 
    int32_t *NP_P, 
    multifok_term_prod_t **prix_P, 
    bool_t verbose
  )
  {
    int32_t NB = basis->NB;
    int32_t NP_max = NB*(NB + 1)/2 + 1; /* Allocated size of {prodName,prix}. */
    multifok_term_prod_t *prix = talloc(NP_max, multifok_term_prod_t);
    
    int32_t NP = 0; /* Total number of products in all terms. */
    
    for (uint32_t kt = 0;  kt < NT; kt++)
      { char *tnk = termName[kt];
        if (verbose) { fprintf(stderr, "  parsing term %d = \"%s\"\n", kt, tnk); }
        char *pbeg = tnk; /* Netx char to parse in {tnk}. */
        bool_t plus_next = FALSE; /* True is the next char may be '+'. */
        while ((*pbeg) != 0)
          { /* if (verbose) { fprintf(stderr, "    looking at \"%s\"\n", pbeg); } */
            if ((*pbeg) == '+') 
              { demand(plus_next, "spurious '+' char in term");
                pbeg++;
                plus_next = FALSE;
                continue;
              }
            else
              { demand(! plus_next, "missing '+' char in term"); }
            /* There seems to be another product */
            char *pend = pbeg + 1; /* Past the end of the product. */
            /* if (verbose) { fprintf(stderr, "      (*pbeg) = %c\n", (*pbeg)); } */
            if ((*pbeg) != '1')
              { demand(((*pbeg) >= 'A') && ((*pbeg) <= 'Z'), "invalid basis name");  }
            while 
              ( ( ((*pend) >= '0') && ((*pend) <= '9') ) ||
                ( ((*pend) >= 'A') && ((*pend) <= 'Z') ) ||
                ( ((*pend) >= 'a') && ((*pend) <= 'z') ) ||
                ( ((*pend) == '*') )
              )
              { pend++; }
            int32_t nc = (int32_t)(pend - pbeg);
            char *pr;
            int32_t jb1, jb2;
            multifok_term_parse_product(pbeg, nc, NB, basis->belName, &jb1, &jb2, &pr);

            /* Check for repeated products: */
            for (uint32_t rp = 0;  rp < NP; rp++)
              { demand(strcmp(pr, prix[rp].pname) != 0, "repeated product in term names"); }

            /* Store product in tables: */
            demand (NP < NP_max, "too many products");
            int32_t ip = NP;
            prix[ip] = (multifok_term_prod_t){ jb1, jb2, kt, pr };
            NP++;
    
            if (verbose) { fprintf(stderr, "   product %d = %d %d %.*s +-> term %d\n", ip, jb1, jb2, nc, pr, kt); }

            plus_next = TRUE;
            pbeg = pend;
          }
      }
    if (NP < NP_max)
      { prix = (multifok_term_prod_t*)notnull(realloc(prix, NP*sizeof(multifok_term_prod_t)), "no mem"); }
    (*NP_P) = NP;
    (*prix_P) = prix;
  
    return;
  }

void multifok_term_parse_product
  ( char *pbeg, 
    int32_t nc, 
    int32_t NB, 
    char *belName[], 
    int32_t *jb1_P, 
    int32_t *jb2_P,
    char ** pr_P
  )
  { 
    char *pend = pbeg+nc;
   
    auto int32_t find_belName(char *beg, char *end);
      /* Returns the index {ib} such that {belName[ib]} is the string between 
        {*beg} (inclusive) and (*end) (exclusive). */
  
    int32_t jb1, jb2;
    if ((nc == 1) && (strncmp(pbeg, "1", nc) == 0))
      { /* Constant term: */
        jb1 = -1; jb2 = -1;
      }
    else 
      { char *past = strchr(pbeg, '*');
        demand((past != NULL) && (past < pend), "linear terms not allowed");
        jb1 = find_belName(pbeg, past);
        jb2 = find_belName(past+1, pend);
        /* Ensure factors are in canonical order: */
        if (jb1 > jb2) { int32_t tmp = jb1; jb1 = jb2; jb2 = tmp; }
      }
    /* Create a fresh copy {pr} of product name, with factors sorted: */
    if (jb1 == -1)
      { char *pr = jsprintf("1"); }
    else
      { char *pr = jsprintf("%s*%s", belName[jb1], belName[jb2]); }
    
    (*jb1_P) = jb1;
    (*jb2_P) = jb2;
    (*pr_P) = pr;
    
    return;

    int32_t find_belName(char *beg, char *end)
      { int32_t n = (int32_t)(end - beg);
        int32_t ib = -1;
        for (uint32_t kb = 0;  kb < NB; kb++)
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

void multifok_term_set_names_write
  ( FILE *wr, 
    int32_t NT,
    char *termName[]
  )
  { for (uint32_t kt = 0;  kt < NT; kt++)
      { fprintf(wr, "%s\n", termName[kt]); }
    fflush(wr);
  }
 
void multifok_term_set_names_write_named
  ( char *outPrefix, 
    int32_t NT, 
    char *termName[]
  )
  { char *fname = NULL;
    char *fname = jsprintf("%s-tnames.txt", outPrefix);
    FILE *wr = open_write(fname, TRUE);
    multifok_term_set_names_write(wr, NT, termName);
    fclose(wr);
  }

void multifok_term_set_write
  ( FILE *wr, 
    multifok_term_set_t *tset,
    bool_t weights,
    double wt[]
  )
  { int32_t NT = tset->NT;
    for (uint32_t kt = 0;  kt < NT; kt++)
      { if (weights)
          { double wtk = (wt == NULL ? 1.0 : wt[kt]);
            fprintf(wr, "%4d %12.8f %s\n", kt, wtk, tset->termName[kt]); 
          }
        else
          { fprintf(wr, "%4d 1 %s\n", kt, tset-> termName[kt]); }
      }
    fflush(wr);
  }

void multifok_term_set_write_named
  ( char *outPrefix,
    multifok_term_set_t *tset,
    bool_t weights,
    double wt[]
  )
  { char *fname = NULL;
    char *fname = jsprintf("%s-twts.txt", outPrefix);
    FILE *wr = open_write(fname, TRUE);
    multifok_term_set_write(wr, tset, weights, wt);
    fclose(wr);
  }

void multifok_term_set_product_table_write
  ( FILE *wr, 
    int32_t NP,
    multifok_term_prod_t *prix,
    bool_t verbose
  )
  {
    demand(NP >= 1, "invalid {NP}");
    for (uint32_t ip = 0;  ip < NP; ip++)
      { multifok_term_prod_t *pri = &(prix[ip]);
        fprintf(wr, "%4d %4d %4d", ip, pri->jb1, pri->jb2);
        fprintf(wr, "  %4d %s\n", pri->kt, pri->pname);
      }
    fflush(wr);
  }

void multifok_term_set_product_table_write_named
  ( char *outPrefix, 
    int32_t NP, 
    multifok_term_prod_t prix[],
    bool_t verbose
  )
  {
    char *fname = jsprintf("%s-prix.txt", outPrefix);
    FILE *wr = open_write(fname, TRUE);
    multifok_term_set_product_table_write(wr, NP, prix, verbose);
    fclose(wr);
    free(fname);
  }

#define multifok_term_C_COPYRIGHT \
    "© 2022 by the State University of Campinas (UNICAMP)"

