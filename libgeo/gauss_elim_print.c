/* See gauss_elim.h */
/* Last edited on 2024-11-25 02:00:53 by stolfi */

#include <stdint.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include <assert.h>
#include <affirm.h>
#include <bool.h>
#include <rmxn.h>
#include <jsmath.h>

#include <gauss_elim.h>
#include <gauss_elim_triangularize.h>
#include <gauss_elim_diagonalize.h>
#include <gauss_elim_normalize.h>

/* IMPLEMENTATIONS */

void gauss_elim_print_array
  ( FILE *wr,
    uint32_t indent,
    char *fmt,
    char *head,
    uint32_t m,
    uint32_t n,
    char *Mname,
    double M[],
    char *foot
  )
  { 
    gauss_elim_print_system(wr, indent, fmt, head, m, n,Mname,M, 0,NULL,NULL, 0,NULL,NULL, foot);
  }

void gauss_elim_print_system
  ( FILE *wr, 
    uint32_t indent,
    char *fmt, 
    char *head, 
    uint32_t m, 
    uint32_t n, 
    char *Aname,
    double A[], 
    uint32_t p, 
    char *Bname,
    double B[], 
    uint32_t q,
    char *Cname,
    double C[], 
    char *foot
  )
  { if (Aname == NULL) { Aname = ""; }
    char *Aeq = (strlen(Aname) > 0 ? " = " : "");
    
    if (Bname == NULL) { Bname = ""; }
    char *Beq = (strlen(Bname) > 0 ? " = " : "");
    
    if (Cname == NULL) { Cname = ""; }
    char *Ceq = (strlen(Cname) > 0 ? " = " : "");
    
    auto void print_row(uint32_t i, char *name, char *eq, uint32_t r, uint32_t s, double M[]);
      /* Prints row {i} of the matrix {M}, assumed to have {r} rows and {s} columns, 
        on the current line, without newline. */
    
    if (head != NULL) { fprintf(wr, "%*s%s\n", indent, "", head); }
    for (uint32_t i = 0;  i < m; i++)
      { if (indent > 0) { fprintf(wr, "%*s", indent, ""); }
        if (A != NULL) { print_row(i, Aname, Aeq, m, n, A); }
        if (B != NULL) { print_row(i, Bname, Beq, m, p, B); }
        if (C != NULL) { print_row(i, Cname, Ceq, m, q, C); }
        fprintf(wr, "\n");
      }
    if (foot != NULL) { fprintf(wr, "%*s%s\n", indent, "", foot); }
    return;
    
    void print_row(uint32_t i, char *name, char *eq, uint32_t r, uint32_t s, double M[])
      { int32_t skip = (int32_t)(strlen(name) + strlen(eq));
        if (i == (int32_t)(r-1)/2)
          { fprintf(wr, "  %s%s[ ", name, eq); }
        else
          { fprintf(wr, "  %*s[ ", skip, ""); }
        for (uint32_t k = 0;  k < s; k++)
          { fprintf(wr, " "); fprintf(wr, fmt, M[i*s + k]); }
        fprintf(wr, " ]");
      }
  }
