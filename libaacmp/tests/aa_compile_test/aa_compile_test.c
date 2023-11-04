/* Last edited on 2023-10-15 03:40:51 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <pcode.h>
#include <affirm.h>
#include <ia.h>
#include <jsstring.h>

#include <aa_compile.h>

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  { char *ename;
    FILE *p_file;
    FILE *c_file;
    pcode_proc_t p;
    affirm (argc == 2, "usage: aa_compile_test exprname");
    ename = argv[1];
    p_file = fopen(txtcat3("out/", ename, ".pcode"), "r");
    affirm (p_file != NULL, "can't open input file");
    c_file = fopen(txtcat3("out/", ename, ".c"), "w");
    affirm (c_file != NULL, "can't open output file");
    p = pcode_parse (p_file);
    aa_compile(&p, NULL, NULL, 0, ename, c_file);
    return(0);
  }
