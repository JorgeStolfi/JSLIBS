/* See {multifok_term.h}. */
/* Last edited on 2023-01-30 19:30:40 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <vec.h>
#include <wt_table.h>
#include <affirm.h>
#include <bool.h>
#include <fget.h>

#include <multifok_window.h>
#include <multifok_basis.h>
#include <multifok_term.h>

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
  ( int32_t NB, 
    double coeff[],
    int32_t NP,
    multifok_term_prod_t prix[],
    int32_t NT,
    double term[]
  )
  {
    for (int32_t kt = 0; kt < NT; kt++) { term[kt] = 0; }
    for (int32_t ip = 0; ip < NP; ip++) 
      { int32_t jb1 = prix[ip].jb1; demand((jb1 >= 0) && (jb1 < NB), "bad basis index {jb1}");
        int32_t jb2 = prix[ip].jb2; demand((jb2 >= 0) && (jb2 < NB), "bad basis index {jb2}");
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
        int32_t kt = prix[ip].kt;
        term[kt] += tvk;
      }
  }

void multifok_term_indices_from_names
  ( int32_t NB, 
    char *belName[],
    int32_t NT, 
    char *termName[], 
    int32_t *NP_P, 
    multifok_term_prod_t **prix_P, 
    bool_t verbose
  )
  {
    int32_t NP_max = NB*(NB + 1)/2 + 1; /* Allocated size of {prodName,prix}. */
    multifok_term_prod_t *prix = (multifok_term_prod_t*)notnull(malloc(NP_max*sizeof(multifok_term_prod_t)), "no mem");
    
    int32_t NP = 0; /* Total number of products in all terms. */
    
    for (int32_t kt = 0; kt < NT; kt++)
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
              { demand(((*pbeg) >= 'A') && ((*pbeg) <= 'Z'), "invalid term name");  }
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
            multifok_term_parse_product(pbeg, nc, NB, belName, &jb1, &jb2, &pr);

            /* Check for repeated products: */
            for (int32_t rp = 0; rp < NP; rp++)
              { demand(strcmp(pr, prix[rp].name) != 0, "repeated product in term names"); }

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
        demand(jb1 <= jb2, "product is not canonical ({jb1 > jb2})");
      }
    /* Create a fresh copy {pr} of product name: */
    char *pr = NULL;
    asprintf(&pr, "%*s", nc, pbeg);
    
    (*jb1_P) = jb1;
    (*jb2_P) = jb2;
    (*pr_P) = pr;
    
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

void multifok_term_read_index_table
  ( FILE *rd,
    int32_t NB,
    char *belName[],
    int32_t *NP_P,
    multifok_term_prod_t **prix_P,
    int32_t *NT_P,
    char ***termName_P,
    bool_t verbose
  )
  {
    int32_t NP_max = NB*(NB+1)/2 + NB + 1; /* Number of coeffs, canonical coeff pairs, and the unit term. */

    multifok_term_prod_t *prix = (multifok_term_prod_t *)notnull(malloc(NP_max*sizeof(multifok_term_prod_t)), "no mem");
    for (int32_t ip = 0; ip < NP_max; ip++) 
      { prix[ip] = (multifok_term_prod_t){ -1,-1, -1, NULL }; }
    
    int32_t NJJ = NB*NB; /* Number of possible basis element pairs {jb1,jb2} in any order. */
    bool_t jjseen[NJJ]; /* {jjseen[jb1*NB+jb2]} is true if the pair {jb1,jb2} already occurred. */
    for (int32_t jjb = 0; jjb < NJJ; jjb++) { jjseen[jjb] = FALSE; }

    int32_t NT_max = NP_max;
    int32_t nr_kt[NT_max]; /* Number of products assigned to each term. */
    for (int32_t kt = 0; kt < NT_max; kt++) { nr_kt[kt] = -1; }
    
    int32_t NT = 0; /* Number of terms seen in file. */
    int32_t NP = 0; /* Number of products read from file. */
   
    auto void read_data_line(int32_t *jb1_P, int32_t *jb2_P, int32_t *kt_P, char **pr_P);
      /* Reads a non-blank data line and returns the data in {*jb1_P,*jb2_P,*kt_P,*pr_P}. */

    if (verbose) { fprintf(stderr, "reading index triples...\n"); }

    while (TRUE)
      { fget_skip_spaces(rd);
        int32_t c = fgetc(rd);
        if (c == EOF) 
          { break; }
        else if (c == '\n') 
          { /* Blank line, ignore. */ }
        else if (c == '#') 
          { /* Comment line, skip rest: */ 
            fget_skip_to_eol(rd); 
          }
        else if ((c < '0') || (c > '9'))
          { /* Bad line: */
            demand(FALSE, "bad line format, expected digit.");
          }
        else
          { /* Looks like a data line: */
            ungetc(c, rd);
            int32_t jb1, jb2, kt;
            char *pr;
            read_data_line(&jb1, &jb2, &kt, &pr);
            if (kt == NT)
              { /* Should be a new term: */
                nr_kt[kt] = 1;
                NT++;
              }
            else
              { /* Previously started term: */
                assert(nr_kt[kt] > 0);
                nr_kt[kt]++;
              }
            prix[NP] = (multifok_term_prod_t){ jb1, jb2, kt, pr };
            if (verbose) { fprintf(stderr, "  %3d  %3d %3d  %3d %s\n", NP, jb1, jb2, kt, pr); }
            NP++;
          }
      }
      
   if (NP < NP_max) { prix = (multifok_term_prod_t*)notnull(realloc(prix, NP*sizeof(multifok_term_prod_t)), "no mem"); }
   if (verbose) { fprintf(stderr, "building term names...\n"); }
   char **termName = (char **)notnull(malloc(NT*sizeof(char*)), "no mem");
   for (int32_t kt = 0; kt < NT; kt++) { termName[kt] = NULL; }
   for (int32_t ip = 0; ip < NP; ip++)
     { char *pr = prix[ip].name;
       int32_t kt = prix[ip].kt;
       char *otm = termName[kt];
       char *tm = NULL;
       if (otm == NULL)
         { asprintf(&tm, "%s", pr); }
       else
         { asprintf(&tm, "%s+%s", otm, pr);
           free(otm);
         }
       termName[kt] = tm;
     }
     if (verbose) 
       { for (int32_t kt = 0; kt < NT; kt++) 
           { fprintf(stderr, "  %3d  %s\n", kt, termName[kt]); }
       }
       
    (*NP_P) = NP;
    (*prix_P) = prix;
    (*NT_P) = NT;
    (*termName_P) = termName;

    return;
    
    /* INTERNAL IMPLEMENTATIONS */
    
    void read_data_line(int32_t *jb1_P, int32_t *jb2_P, int32_t *kt_P, char **pr_P)
      { int32_t ip = fget_int32(rd);
        demand(ip == NP, "unexpected product index");

        int32_t jb1 = fget_int32(rd);
        demand((jb1 >= 0) && (jb1 < NB), "invalid basis coeff index {jb1}");

        int32_t jb2 = fget_int32(rd);
        demand((jb2 >= 0) && (jb2 < NB), "invalid basis coeff index {jb2}");
        demand(jb1 <= jb2, "non-canonical produc {jb1>jb2}");
        
        int32_t jjb = jb1*NB + jb2;
        demand(! jjseen[jjb], "repeated produc {jb1,jb2}");
        jjseen[jjb] = TRUE;
        
        int32_t kt = fget_int32(rd);
        demand(kt > 0, "invalid {kt}");
        demand(kt <= NT, "{kt} skipped term");

        char *pr_read = fget_string(rd);
        /* Validate {pr_read} with {jb1,jb2}: */
        char *pr = NULL;
        asprintf(&pr, "%s*%s", belName[jb1], belName[jb2]);
        demand(strcmp(pr_read, pr) == 0, "product name does not match {jb1,jb2}");
        free(pr_read);
        
        fget_comment_or_eol(rd, '#');
        demand(kt <= NT, "{kt} skipped term");
        
        (*jb1_P) = jb1;
        (*jb2_P) = jb2;
        (*kt_P) = kt;
        (*pr_P) = pr;
      }
  }


void multifok_term_write_names(FILE *wr, int32_t NT, char *termName[])  
  { for (int32_t kt = 0; kt < NT; kt++)
      { fprintf(wr, "%s\n", termName[kt]); }
    fflush(wr);
  }


#define multifok_term_C_COPYRIGHT \
    "© 2022 by the State University of Campinas (UNICAMP)"

