/* See pcode.h */
/* Last edited on 2016-04-01 01:03:45 by stolfilocal */

#include "pcode.h"

#include <affirm.h>
#include <jsfile.h>

#include <stdio.h>
#include <malloc.h>
#include <string.h>
#include <stdlib.h>

/*** INSTRUCTION INFO TABLE ***/

typedef struct {
    pcode_op_t op;   /* The internal opcode */
    char *name;      /* The external opcode */
    int has_arg;     /* TRUE if takes an argument */
    int arg_is_reg;  /* TRUE if takes a register argument */
    int stack_used;  /* Number of operands expected on stack */
    int stack_delta; /* Number of entries added/removed from stack */
  } pcode_op_info_t;
  
#define MAXOPCODE pcode_op_return

pcode_op_info_t pcode_op_info[MAXOPCODE+1] =
  /* Must be listed here in order of opcode: */
  { 
    {pcode_op_given,  "GIVEN",  1, 0,  0,  0},
    {pcode_op_const,  "CONST",  1, 0,  0,  1},
    {pcode_op_load,   "LOAD",   1, 1,  0,  1},
    {pcode_op_store,  "STORE",  1, 1,  1, -1},
    {pcode_op_add,    "+",      0, 0,  2, -1},
    {pcode_op_sub,    "-",      0, 0,  2, -1},
    {pcode_op_neg,    "NEG",    0, 0,  1,  0},
    {pcode_op_inv,    "INV",    0, 0,  1,  0},
    {pcode_op_mul,    "*",      0, 0,  2, -1},
    {pcode_op_div,    "/",      0, 0,  2, -1},
    {pcode_op_sqr,    "SQR",    0, 0,  1,  0},
    {pcode_op_sqrt,   "SQRT",   0, 0,  1,  0},
    {pcode_op_abs,    "ABS",    0, 0,  1,  0},
    {pcode_op_max,    "MAX",    0, 0,  2, -1},
    {pcode_op_min,    "MIN",    0, 0,  2, -1},
    {pcode_op_return, "RETURN", 1, 0,  0,  0}
  };

/*** INTERNAL PROTOTYPES ***/
  
char *pcode_parse_title(FILE *f, char **next_line);

void pcode_parse_instr(char *s, pcode_op_info_t **opinfo, int *arg);

/*** IMPLEMENTATIONS ***/

pcode_proc_t pcode_parse (FILE *f)
  {
    pcode_proc_t proc;
    char *s;

    /* Read title and number of arguments */
    
    proc.title = pcode_parse_title(f, &s);
    proc.nstack = 0;
    proc.nin = 0;    /* Just in case (GCC pacifier). */
    proc.nout = 0;   /* Just in case (GCC pacifier). */

    /* Parse instructions until RETURN: */
    { pcode_instr_t *pi;   /* The instructions read so far */
      int ni = 0;          /* Current instruction count */
      int mi = 10;         /* Current allocated size of $*pi$ */

      int top = -1;        /* Simulated top-of-stack */
      int max_top = -1;    /* Simulated top-of-stack */
      
      int max_reg_arg = -1;  /* Largest register number ever mentioned */
      
      pcode_op_info_t *opinfo;
      int arg;             

      /* should check whether any reg is used before defined. */

      pi = (pcode_instr_t *) malloc(mi * sizeof(pcode_instr_t));
      do
	{
	  affirm(s != NULL, "pcode_parse: missing RETURN");

          pcode_parse_instr(s, &opinfo, &arg);
          
          if (opinfo != NULL)
            { 
	      /* Process GIVEN directive: */
	      if (opinfo->op == pcode_op_given)
		{ affirm(ni == 0, "pcode_parse: GIVEN out of place");
		  affirm(arg >= 0, "pcode_parse: invalid num of inputs");
		  proc.nin = arg;
		  max_reg_arg = arg+1;
		}
	      else 
		{ affirm(ni > 0, "pcode_parse: missing GIVEN directive"); }

	      /* Process RETURN directive:*/
	      if (opinfo->op == pcode_op_return) 
		{ affirm(arg >= 0, "pcode_parse: invalid num of outputs");
		  affirm(arg <= top+1, "pcode_parse: return of undefined register");
		  proc.nout = arg;
		}

	      /* Update number of registers used: */
	      if (opinfo->arg_is_reg)
		{ affirm(arg >= 0, "pcode_parse: invalid register argument");
		  if (arg >= max_reg_arg) max_reg_arg = arg;
		}

	      /* simulate stack motion: */
	      affirm (top >= opinfo->stack_used - 1, "pcode_parse: eval stack underflow");
	      top += opinfo->stack_delta;
	      affirm (top >= -1, "pcode_parse: eval stack underflow");
	      if (top > max_top) max_top = top;

	      /* store in array */
	      if (ni >= mi)
		{ mi *= 2; 
		  pi = (pcode_instr_t *) realloc ((void*) pi, mi*sizeof(pcode_instr_t)); 
		}
	      affirm(pi != NULL, "pcode_parse: alloc failed");
	      pi[ni].op = opinfo->op;
	      pi[ni].arg = arg;
	      ni++;
            }

          free ((void *) s);
          s = read_line(f);
	}
      while (opinfo->op != pcode_op_return);
      proc.nregs = max_reg_arg+1;
      proc.nstack = max_top+1;
      proc.nops = ni;
      proc.code = (pcode_instr_t *) realloc ((void *) pi, ni*sizeof(pcode_instr_t));
      affirm(proc.code != NULL, "pcode_parse: alloc failed");
    }
    return (proc);
  }

char *pcode_parse_title(FILE *f, char **next_line)
  /* 
    Parses the title, returns it as a single null-terminated string,
    with lines separated by '\n'.  Also stores in *next_line the first input line
    after the title. */
  { 
    char **t;      /* The title lines */
    int nt = 0;    /* Number of title lines read */
    int nc = 0;    /* Number of title characters read */
    int maxt = 10; /* current allocated size of $t$ */
    int i;
    char *p, *q, *s, *title;

    t = (char **)malloc(sizeof(char*)*maxt);

    s = read_line(f);
    while ((s != NULL) && (s[0] == '#'))
      { if (nt >= maxt)
	  { maxt *= 2; t = (char**) realloc((void*)t, sizeof(char*)*maxt); }
	affirm(t != NULL, "pcode_parse_title: alloc failed");
	t[nt++] = s;
	nc += (int)strlen(s);
	s = read_line(f);
      }
    *next_line = s;

    title = (char *) malloc(sizeof(char) * (nc + (nt - 1) + 1));
    affirm(title != NULL, "pcode_parse_title: alloc failed");
    p = title; 
    for (i=0; i<nt; i++)
      { q = t[i];
        if (i>0) *(p++) = '\n';
	while (*q != '\000') *(p++) = *(q++); 
	*p = '\000';
	free((void*) t[i]);
      }
    free ((void*) t);
    return(title);
  }
  
void pcode_parse_instr(char *s, pcode_op_info_t **opinfo, int *arg)
  {
    pcode_op_info_t *info;
    int iarg;
    char *tok;
    
    /* get and decode the opcode: */
    tok = strtok(s, " \n\t");
    if (tok == NULL)
      { *opinfo = NULL; *arg = 0; return; }
    
    /* decode the opcode */
    { int j;
      for (j=0; (j <= MAXOPCODE) && (strcmp(tok, pcode_op_info[j].name) != 0); j++) { }
      affirm(j <= MAXOPCODE, "pcode_parse_instr: invalid op code");
      info = &(pcode_op_info[j]);
    }

    if (info->has_arg)
      { /* Get and decode the argument */
	char *a = strtok(NULL, " \n\t");
	affirm(a != NULL, "pcode_parse_instr: missing argument");
	iarg = atoi(a);
      }
    else
      { iarg = 0; }
    
    /* Check for bogus data: */
    { char *a = strtok(NULL, " \n\t");
      affirm (a == NULL, "pcode_parse_instr: extra characters in line");
    }
    
    /* Return results */
    *opinfo = info;
    *arg = iarg;
  }

void pcode_print (FILE *f, pcode_proc_t proc)
  { pcode_op_info_t *opinfo;
    pcode_op_t op;
    int arg; 
    int i = 0;
    fprintf(f, "%s\n", proc.title);
    do
      { op = proc.code[i].op;
        arg = proc.code[i].arg;
        opinfo = &pcode_op_info[op];
        affirm(opinfo->op == op, "pcode_print: bad table");
          
        fprintf(f, "%s", pcode_op_info[op].name);
        if (opinfo->has_arg) fprintf(f, " %d", arg);
        fprintf(f, "\n");
        if (op == pcode_op_given) 
          affirm(arg = proc.nin, "pcode_print: inconsistent input arity");
        if (op == pcode_op_return) 
          affirm(arg = proc.nout, "pcode_print: inconsistent output arity");
        i++;
      }
    while (op != pcode_op_return);
    affirm(i == proc.nops, "pcode_print: inconsistent code length");
    fflush(f);
  }
   
