/* See pcode.h */
/* Last edited on 2024-11-20 08:06:28 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

#include <jsfile.h>
#include <affirm.h>

#include <pcode.h>

/*** INSTRUCTION INFO TABLE ***/

typedef struct {
    pcode_op_t op;   /* The internal opcode */
    char *name;      /* The external opcode */
    int32_t has_arg;     /* TRUE if takes an argument */
    int32_t arg_is_reg;  /* TRUE if takes a register argument */
    int32_t stack_used;  /* Number of operands expected on stack */
    int32_t stack_delta; /* Number of entries added/removed from stack */
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
  
char *pcode_parse_title(FILE *f, char **next_line_P);

void pcode_parse_instr(char *s, pcode_op_info_t **opinfo, int32_t *arg);

/*** IMPLEMENTATIONS ***/

void pcode_proc_null(pcode_proc_t *pc)
  { pc->title  = NULL;
    pc->nin    = 0;
    pc->nout   = 0;
    pc->nregs  = 0;
    pc->nstack = 0;
    pc->nops   = 0;
    pc->code   = NULL;
  }
  
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
      int32_t ni = 0;          /* Current instruction count */
      int32_t mi = 10;         /* Current allocated size of $*pi$ */

      int32_t top = -1;        /* Simulated top-of-stack */
      int32_t max_top = -1;    /* Simulated top-of-stack */
      
      int32_t max_reg_arg = -1;  /* Largest register number ever mentioned */
      
      pcode_op_info_t *opinfo;
      int32_t arg;             

      /* should check whether any reg is used before defined. */

      pi = talloc(mi, pcode_instr_t);
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
		  pi = talloc(pi, pcode_instr_t); 
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
      proc.code = retalloc (pi, ni, pcode_instr_t);
      affirm(proc.code != NULL, "pcode_parse: alloc failed");
    }
    return (proc);
  }

char *pcode_parse_title(FILE *f, char **next_line_P)
  /* 
    Parses the title, returns it as a single null-terminated string,
    with lines separated by '\n'.  Also stores in {*next_line_P} the first input line
    after the title. */
  { 
    uint32_t nt = 0;    /* Number of title lines read */
    uint32_t nc = 0;    /* Number of title characters read */
    uint32_t maxt = 10; /* Current allocated size of {t} */

    char **t = talloc(maxt, char*);  /* The title lines */

    char *s = read_line(f);
    while ((s != NULL) && (s[0] == '#'))
      { if (nt >= maxt) { maxt *= 2; t = retalloc(t, maxt, char*); }
	t[nt++] = s;
        nc += (uint32_t)strlen(s);
	s = read_line(f);
      }
    *next_line_P = s;

    char *title = talloc(nc + (nt - 1) + 1, char);
    char *p = title; 
    for (uint32_t i = 0;  i < nt; i++)
      { char *q = t[i];
        if (i>0) *(p++) = '\n';
	while (*q != '\000') *(p++) = *(q++); 
	*p = '\000';
	free(t[i]);
      }
    free (t);
    return(title);
  }
  
void pcode_parse_instr(char *s, pcode_op_info_t **opinfo, int32_t *arg)
  {
    pcode_op_info_t *info = NULL;
    int32_t iarg = 0;
    
    /* Get the opcode: */
    char *tok = strtok(s, " \n\t");
    if (tok != NULL) 
      { /* Decode the opcode */
        { int32_t j = 0;
          while ((j <= MAXOPCODE) && (strcmp(tok, pcode_op_info[j].name) != 0)) { j++; }
          affirm(j <= MAXOPCODE, "pcode_parse_instr: invalid op code");
          info = &(pcode_op_info[j]);
        }

        if (info->has_arg)
          { /* Get and decode the argument */
            char *a = strtok(NULL, " \n\t");
            affirm(a != NULL, "pcode_parse_instr: missing argument");
            iarg = atoi(a);
          }

        /* Check for bogus data: */
        { char *a = strtok(NULL, " \n\t");
          affirm (a == NULL, "pcode_parse_instr: extra characters in line");
        }
      }
    
    /* Return results */
    (*opinfo) = info;
    (*arg) = iarg;
  }

void pcode_print (FILE *f, pcode_proc_t proc)
  { pcode_op_info_t *opinfo;
    pcode_op_t op;
    int32_t arg; 
    int32_t i = 0;
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
   
