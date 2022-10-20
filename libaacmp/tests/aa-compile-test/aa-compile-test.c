#include <aa_compile.h>
#include <stdint.h>
#include <ia.h>
#include <affirm.h>
#include <pcode.h>
#include <stdlib.h>
#include <stdio.h>

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
    p_file = fopen(txtcat(ename, ".pcode"), "r");
    affirm (p_file != NULL, "can't open input file");
    c_file = fopen(txtcat(ename, ".c"), "w");
    affirm (c_file != NULL, "can't open output file");
    p = pcode_parse (p_file);
    aa_compile(&p, NULL, NULL, 0, "foo", c_file, NULL);
    return(0);
  }
