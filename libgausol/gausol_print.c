/* See gausol_print.h */
/* Last edited on 2024-11-29 18:05:17 by stolfi */

#include <stdint.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include <assert.h>
#include <affirm.h>
#include <jsprintf.h>
#include <bool.h>
#include <jsmath.h>

#include <gausol_print.h>

/* IMPLEMENTATIONS */

void gausol_print_array
  ( FILE *wr,
    uint32_t indent,
    char *fmt,
    char *head,
    uint32_t m, uint32_t prow[], uint32_t rh,
    uint32_t n, uint32_t pcol[], uint32_t rv,
    char *Mname, double M[], 
    char *foot
  )
  { 
    gausol_print_system
      ( wr, indent, fmt, head, 
        m,prow,rh,
        n,pcol,rv, Mname,M, 
        0,NULL,NULL,
        0,NULL,NULL,
        foot
      );
  }

void gausol_print_system
  ( FILE *wr, 
    uint32_t indent,
    char *fmt, 
    char *head, 
    uint32_t m, uint32_t prow[], uint32_t rh,
    uint32_t n, uint32_t pcol[], uint32_t rv, char *Aname, double A[], 
    uint32_t p, char *Bname, double B[], 
    uint32_t q, char *Cname, double C[],
    char *foot
  )
  { if (Aname == NULL) { Aname = ""; }
    char *Aeq = (strlen(Aname) > 0 ? " = " : "");
    demand(((m == 0) || (n == 0)) == (A == NULL), "{A} should be {NULL} iff size is zero");
    
    if (Bname == NULL) { Bname = ""; }
    char *Beq = (strlen(Bname) > 0 ? " = " : "");
    demand(((m == 0) || (p == 0)) == (B == NULL), "{B} should be {NULL} iff size is zero");
    
    if (Cname == NULL) { Cname = ""; }
    char *Ceq = (strlen(Cname) > 0 ? " = " : "");
    demand(((m == 0) || (q == 0)) == (C == NULL), "{C} should be {NULL} iff size is zero");
    
    auto void print_data_row
      ( uint32_t i1, char *name, char *eq, bool_t pneq,
        uint32_t m1, uint32_t n1, double M[], 
        uint32_t rv1, uint32_t pcol1[]
      );
      /* Prints on the current line of {wr} row {i1} of the matrix {M},
        surrounded by "[ " and " ]", without a final newline. The matrix
        assumed to have {m1} rows and {n1} columns.
        
        Specifically, prints two spaces. Then, if {pneq} is true, prints
        {name} and {eq}, otherwise prints as many spaces as the lengths
        of those two strings. Then prints "[ ". Then prints the {n1}
        elements {M[i1,j]} for {j} in {0..n1-1}, formatted by {fmt} and
        separated by single spaces. If {pcol1} is not
        {NULL}, prints {M[i1,pcol1[j]]} instead of {M[i1,j]}.
        
        If {rv1} is in {1..n1-1}, prints a vertical bar separating the 
        first {rv1} elements from the last {n1-rv1}.
        
        The matrix {M} must be non-{NULL} and {m1,n1} must be nonzero,
        with {i1<m1}. */
        
    auto void print_dash_row
      ( char *name, char *eq, bool_t pneq, uint32_t elwd,
        uint32_t n1, uint32_t rv1
      );
      /* Prints {n1} strings of {elwd} dashes.
        
        Specifically, prints two spaces. Then, if {pneq} is true, prints
        {name} and {eq}, otherwise prints as many spaces as the lengths
        of those two strings. Then prints "[ ". Then prints the {n1}
        strings of dashes, separated by single spaces. Then prints " ]".
        
        If {rv1} is in {1..n1-1}, prints a "+" between the first {rv1} strings
        and the last {n1-rv1}. */
        
    /* Hack to get the normal width {elwd} of one element: */
    char *xel = jsprintf(fmt, M_PI);
    uint32_t elwd = (uint32_t)strlen(xel);
    free(xel);
        
    bool_t print_perm = ((prow != NULL) || (pcol != NULL));
    /* Determine if the row-permuted matrice have a row of dashes: */
    bool_t has_dashes = ((rh > 0) && (rh < m));
    for (uint32_t vers = 0; vers <= 1; vers++)
      { /* {vers=0} unpermuted, {vers=1} permuted. */
        char *xunp = (vers == 0 ? " (unpermuted)" : " (permuted)");
        if (! print_perm) { xunp = ""; }
        if ((vers == 0) || print_perm)
          { if (vers == 1) { fprintf(wr, "\n"); }
            if (head != NULL) { fprintf(wr, "%*s%s%s\n", indent, "", head, xunp); }
            for (uint32_t i = 0; i < m; i++)
              { uint32_t i1 = ((prow == NULL) || (vers == 0) ? i : prow[i]); assert(i1 < m);
                uint32_t *pcol1 = (vers == 0 ? NULL : pcol);
                uint32_t rh1 = (vers == 0 ? 0 : rh);
                uint32_t rv1 = (vers == 0 ? 0 : rv);
                bool_t pneq_dash, pneq_data;
                gausol_print_row_has_name_eq(i, m, rv1, &pneq_dash, &pneq_data);
                if ((vers == 1) && has_dashes && (i == rh1))
                  { /* Print the rows of dashes: */
                    if (indent > 0) { fprintf(wr, "%*s", indent, ""); }
                    if (A != NULL) { print_dash_row(Aname, Aeq, pneq_dash, elwd, n, rv1); }
                    if (B != NULL) { print_dash_row(Bname, Beq, pneq_dash, elwd, p,   0); }
                    if (C != NULL) { print_dash_row(Cname, Ceq, pneq_dash, elwd, q,   0); }
                    fprintf(wr, "\n");
                  }
                if (indent > 0) { fprintf(wr, "%*s", indent, ""); }
                if (A != NULL) { print_data_row(i1, Aname, Aeq, pneq_data, m, n, A, rv1, pcol1); }
                if (B != NULL) { print_data_row(i1, Bname, Beq, pneq_data, m, p, B,   0,  NULL); }
                if (C != NULL) { print_data_row(i1, Cname, Ceq, pneq_data, m, q, C,   0,  NULL); }
                fprintf(wr, "\n");
              }
            if (foot != NULL) { fprintf(wr, "%*s%s\n", indent, "", foot); }
          }
      }
      
    return;
      
    auto void print_left(bool_t pneq, char *name, char *eq, uint32_t skip);
      /* Prints two spaces to {wr}.  If {pnew} is true, prints {name} and {eq}, otherwise prints {skip} 
        blanks.  Then prints "[". */
       
    void print_dash_row
      ( char *name, char *eq, bool_t pneq, uint32_t elwd,
        uint32_t n1, uint32_t rv1
      )
      { uint32_t skip = (uint32_t)(strlen(name) + strlen(eq));
        /* Determine which row (dashes or data), if any, the {name,eq} should go: */
        print_left(pneq, name, eq, skip);
        for (uint32_t j = 0;  j < n1; j++)
          { fputc(' ', wr); 
            if ((rv1 >= 1) && (rv1 <= n1-1) && (j == rv1)) { fputs("+ ", wr); }
            for (uint32_t k = 0; k < elwd; k++) { fputc('-', wr); }
          }
        fprintf(wr, " ]");
      }

    void print_data_row
      ( uint32_t i1, char *name, char *eq, bool_t pneq,
        uint32_t m1, uint32_t n1, double M[], 
        uint32_t rv1, uint32_t pcol1[]
      )
      { assert(i1 < m1);
        uint32_t skip = (uint32_t)(strlen(name) + strlen(eq));
        print_left(pneq, name, eq, skip);
        for (uint32_t j = 0; j < n1; j++)
          { fputc(' ', wr); 
            uint32_t j1 = (pcol1 == NULL ? j : pcol[j]);
            if ((rv1 >= 1) && (rv1 <= n1-1) && (j == rv1)) { fputs("| ", wr); }
            fprintf(wr, fmt, M[i1*n1 + j1]);
          }
        fprintf(wr, " ]");
      }
        
    void print_left(bool_t pneq, char *name, char *eq, uint32_t skip)
      { if (pneq)
          { fprintf(wr, "  %s%s[", name, eq); }
        else
          { fprintf(wr, "  %*s[", skip, ""); }
      }
  }

void gausol_print_row_has_name_eq
  ( uint32_t i,
    uint32_t m,
    uint32_t rank,
    bool_t *at_dash_P,
    bool_t *at_data_P
  )
  { 
    demand(i < m, "invalid row index {i}");
    assert(m > 0); /* So that {m-1} is OK. */
    bool_t has_dashes = ((rank >= 1) && (rank < m-1));
    bool_t at_dash, at_data;
    if (! has_dashes)
      { /* There is no dash line: */
        at_data = (i == (m-1)/2);
        at_dash = FALSE; 
      }
    else
      { /* There is a dash line. */
        /* Determine the row index, counting the dash line as a row: */
        uint32_t imid = m/2;
        if (rank == imid)
          { /* The {name,eq} go at the dash line: */
            at_data = FALSE;
            at_dash = (i == rank);
          }
        else
          { /* The {name,eq} go on a data line: */
            at_dash = FALSE;
            at_data = (rank > imid ? (i == imid) : (i+1 == imid));
          }
      }
    (*at_dash_P) = at_dash;
    (*at_data_P) = at_data;
  }

