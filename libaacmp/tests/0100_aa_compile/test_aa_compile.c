/* Last edited on 2024-12-05 10:19:15 by stolfi */

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
    affirm (argc == 2, "usage: test_aa_compile exprname");
    ename = argv[1];
    p_file = fopen(txtcat3("out/", ename, ".pcode"), "r");
    affirm (p_file != NULL, "can't open input file");
    c_file = fopen(txtcat3("out/", ename, ".c"), "w");
    affirm (c_file != NULL, "can't open output file");
    p = pcode_parse (p_file);
    aa_compile(&p, NULL, NULL, 0, ename, c_file);
    return(0);
  }
