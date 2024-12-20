/* See aa_compile.h */
/* Last edited on 2024-12-05 10:19:08 by stolfi */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <ia.h>
#include <affirm.h>
#include <pcode.h>

#include <aa_compile.h>

/* INTERNAL PROTOTYPES */

void aa_compute_value_scopes
  ( pcode_proc_t *p, 
    int32_t maxv,
    /* out */
    int32_t *defined, 
    int32_t *last_used,
    int32_t *result,
    int32_t *nv
  );
  /* Examines the p-code $p$ and sets the tables 
    
      defined[iv]   = index of operation that defines value number $iv$
      
      last_used[iv] = index of last operation that uses value number $iv$
      
      result[ir]    = number of value that is to be returned as result number $ir$
      
    Also returns the number $nv$ of distinct values computed by $p$.
      
    Assumes the tables have been allocated with $maxv$ entries each.
  */
  
void aa_compute_epsilon_tables
  ( pcode_proc_t *p, 
    int32_t *defined, 
    int32_t *last_used, 
    int32_t nv, 
    bool_t **arg_depends,
    bool_t **res_depends,
    int32_t shared_ne,
    int32_t maxe,
    /* out */
    int32_t *erreps,
    bool_t **depends, 
    int32_t **reduce, 
    int32_t *ne
  );
  /* Examines the p-code $p$, and the tables $arg_depends$ and $res_depends$,
    filling the tables
    
      erreps[iv]      = number $ie$ of the epsilon-term that subsumes the 
                        rounding, truncation, and eps-reduction errors of 
                        the operation that created value $iv$; or -1 if none.
    
      depends[iv][ie] = TRUE if value number $iv$ depends on epsilon number $ie$
      
      reduce[iv][je]  = number $ie$ (usually equal to erreps[iv]) of the 
                        epsilon-term that should absorb term number
                        $je$ of value number $iv$; or -1 if that term
                        is to be preserved.
                        
   The procedure assumes that these tables have been allocated with at
   least $nv$ rows and $maxe$ columns.  The procedure also returns the 
   number $ne$ of epsilon-variables actually used in the computation.  
 */

void aa_generate_code
  ( pcode_proc_t *p, 
    char *proc_name, 
    FILE *c_file, 
    bool_t **arg_depends, 
    bool_t **res_depends, 
    int32_t shared_ne,
    int32_t *defined, 
    int32_t *last_used, 
    int32_t *result, 
    int32_t nv,
    int32_t *erreps,
    bool_t **depends, 
    int32_t **reduce, 
    int32_t ne
  );
  /* Writes to $c_file$ a C procedure for evaluating the p-code $p$
    using affine arithmetic.  The arguments are the tables and counts
    given to $aa_compile$ or built by by $aa_compute_value_scopes$
    and $aa_compute_epsilon_tables$.
  */

void aa_update_dep_count(
    bool_t *var_depends,  /* One row of the $depends$ table */
    int32_t ne, 
    int32_t *dep_count, 
    int32_t increment
  );
  /*
    Used internally by aa_compute_epsilon_tables. 
  */

char **aa_compute_value_names(int32_t nv);
  /* Returns a table with the synthetic name assigned to value number $iv$,
    for $iv$ in $[0..nv-1]$. */

void aa_generate_proc_header
  ( pcode_proc_t *p,
    char *proc_name,
    FILE *c_file, 
    int32_t *result,
    bool_t **arg_depends, 
    bool_t **res_depends, 
    int32_t shared_ne,
    char **v_name
  );
  /* Called by aa_generate_code.
    Writes to $c_file$ the C procedure header.
  */

void aa_generate_local_decls
  ( pcode_proc_t *p,
    char *proc_name,
    FILE *c_file, 
    int32_t nv,
    bool_t **depends, 
    int32_t ne,
    char **v_name
  );
  /* Called by aa_generate_code.
    Writes to $c_file$ the internal declarations of the 
    generated procedure 
  */
  
void aa_generate_body
  ( pcode_proc_t *p,
    char *proc_name,
    FILE *c_file, 
    int32_t *defined, 
    int32_t *last_used, 
    int32_t nv, 
    int32_t *erreps,
    bool_t **depends,
    int32_t **reduce,
    int32_t ne,
    char **v_name
  );
  /* Called by aa_generate_code.
    Writes to $c_file$ the body of the generated procedure.
  */
  
void aa_generate_proc_trailer
  ( pcode_proc_t *p,
    char *proc_name,
    FILE *c_file,
    int32_t nv,
    int32_t *result,
    bool_t **res_depends, 
    int32_t shared_ne,
    bool_t **depends,
    int32_t ne,
    char **v_name
  );
  /*
    Writes the procedure trailer to $c_file$.
  */
  
void aa_generate_op_const
  ( FILE *c_file, 
    char *name_v,   /* Name of result. */
    int32_t arg,        /* The integer to convert. */
    int32_t ie          /* Id of error epsilon, or -1 if exact. */
  );
  /* Called by aa_generate_body.
    Generates code to compute the AA representation of 
    an integer constant called $*name_v$, with value $arg$.
    If $arg$ is not representable exactly, then the affine form
    will ascribe the rounding error to epsilon number $ie$ 
    (which in that case should be non-negative).
  */

void aa_generate_unary_op
  ( FILE *c_file, 
    char *name_v,       /* Name of result. */
    pcode_op_t op,
    char *name_x,       /* Name of argument. */
    int32_t ie,             /* Id of epsilon of error term, or -1 if exact. */
    bool_t *depends_x,    /* Epsilons appearing in x. */
    int32_t *reduce_v,      /* Epsilon reduction table for result. */
    int32_t ne              /* Number of epsilons. */
  );
  /* Generates code for a unary operation $v = op(x)$.

    The variable names are actually $*name_v$ and $*name_x$.

    The argument $x$ is assumed to depend on the epsilons $je$ for which
    $depends_x[je] = TRUE$.  

    By default, the result $v$ will depend on these same
    epsilons, plus the special epsilon $ie$ which subsumes all rounding
    or truncation errors of $op$.  
    
    However, if $ke = reduce_v[je]$ is
    not -1 for some epsilon number $je$, the corresponding epsilon-term
    $je$ will be omitted from the expansion of $v$, and added onto the
    term $ke$.  (Usually, $ke$ is the same as $ie$.)

    The index $ie$ may be -1 if the operation is always exact and 
    $reduce_v$ is all -1.
  */
  
void aa_generate_binary_op
  ( FILE *c_file, 
    char *name_v,   /* Name of result. */
    pcode_op_t op,
    char *name_x,       /* Name of first argument. */
    char *name_y,       /* Name of second argument. */
    int32_t ie,             /* Id of epsilon of error term, or -1 if exact. */
    bool_t *depends_x,    /* Epsilons appearing in x. */
    bool_t *depends_y,    /* Epsilons appearing in y. */
    int32_t *reduce_v,      /* Epsilon reduction table for result. */
    int32_t ne              /* Number of epsilons. */
  );
  /* Generates code for a binary operation $v = op(x, y)$.

    The variable names are actually $*name_v$, $*name_x$, and $*name_y$.

    The argument $x$ is assumed to depend on the epsilons $je$ for which
    $depends_x[je] = TRUE$; and similarly for $y$.  

    By default, the result $v$ will depend on the union of these
    epsilons, plus the special epsilon $ie$ which subsumes all rounding
    or truncation errors of $op$.  However, th $reduce_v$ table may
    specify that some of the epsilon-terms may be lumped, as
    in aa_generate_unary_op.

    The index $ie$ may be -1 if the operation is always exact and 
    $reduce_v$ is all -1.
  */
  
void aa_generate_op_neg
  ( FILE *c_file, 
    char *name_v,       /* Name of result. */
    char *name_x,       /* Name of argument. */
    int32_t ie,             /* Id of epsilon of error term, or -1 if exact. */
    bool_t *depends_x,    /* Epsilons appearing in x. */
    int32_t *reduce_v,      /* Epsilon reduction table for result. */
    int32_t ne              /* Number of epsilons. */
  );
  /*
    Called by aa_generate_body.  Generates code to evaluate
    $v = -x$.  The operation is exact if there is no reduction.
  */

void aa_generate_op_sqrt
  ( FILE *c_file, 
    char *name_v,       /* Name of result. */
    char *name_x,       /* Name of argument x. */
    int32_t ie,             /* Id of epsilon of error term. */
    bool_t *depends_x,    /* Epsilons appearing in the x. */
    int32_t *reduce_v,      /* Epsilon reduction table for result. */
    int32_t ne              /* Number of epsilons. */
  );
  /*
    Called by aa_generate_body.  Generates code to evaluate
    $v = sqrt(x)$.  Not exact.
  */
  
void aa_generate_op_add
  ( FILE *c_file, 
    char *name_v,       /* Name of result. */
    char *name_x,       /* Name of argument x. */
    char *name_y,       /* Name of argument y. */
    int32_t ie,             /* Id of epsilon of error term. */
    bool_t *depends_x,    /* Epsilons appearing in the x. */
    bool_t *depends_y,    /* Epsilons appearing in the y. */
    int32_t *reduce_v,      /* Epsilon reduction table for result. */
    int32_t ne              /* Number of epsilons. */
  );
  /*
    Called by aa_generate_body.  Generates code to evaluate
    $v = sqrt(x)$.  Not exact (due to rounding errors).
  */
  
void aa_generate_op_mul
  ( FILE *c_file, 
    char *name_v,       /* Name of result. */
    char *name_x,       /* Name of argument x. */
    char *name_y,       /* Name of argument y. */
    int32_t ie,             /* Id of epsilon of error term. */
    bool_t *depends_x,    /* Epsilons appearing in the x. */
    bool_t *depends_y,    /* Epsilons appearing in the y. */
    int32_t *reduce_v,      /* Epsilon reduction table for result. */
    int32_t ne              /* Number of epsilons. */
  );
  /*
    Called by aa_generate_body.  Generates code to compute 
    $v = x * y$.  Not exact.
  */
  
/* IMPLEMENTATIONS */

void aa_compile
  ( pcode_proc_t *p,
    bool_t **arg_depends,
    bool_t **res_depends,
    int32_t shared_ne,
    char *proc_name,
    FILE *c_file
  )
  {
    int32_t maxv, nv;   /* Max and actual num of distinct values computed */
    int32_t maxe, ne;   /* Max and actual num of distinct epsilons appearing in values. */

    int32_t *defined;   /* defined[iv] == iop  iff value number $iv$ 
                       is created by instruction p->code[iop]. */

    int32_t *last_used; /* last_used[iv] == iop  iff p->code[iop] is the 
                       last operation that uses value $iv$. */
                       
    int32_t *erreps;    /* erreps[iv] is the number of the epsilon-term that
                       subsumes the rounding and truncation errors of 
                       value $iv$. */

    bool_t **depends; /* depends[iv][ie] is TRUE  iff value [iv] 
                       depends formally on epsilon[ie]. */

    int32_t **reduce;   /* reduce[iv][je] is ie  iff the term epsilon[je] that 
                       should appear on value [iv] is to be lumped 
                       into epsilon[ie]. */
                        
    int32_t *result;    /* result[ir] == iv means that the $ir$th result 
                       of $p$ is value [iv]. */

    /* Allocate and initialize tables: */
    
    maxv = p->nops;
    maxe = shared_ne + p->nin + p->nops;
    
    /* Allocate internal tables: */
    defined = (int32_t *) malloc(maxv * sizeof (int32_t));
    last_used = (int32_t *) malloc(maxv * sizeof (int32_t));
    depends = (bool_t **) malloc(maxv * sizeof(bool_t*));
    erreps = (int32_t *) malloc (maxv * sizeof(int32_t));
    reduce = (int32_t **) malloc(maxv * sizeof(int32_t*));
    result = (int32_t *) malloc(p->nout * sizeof(int32_t));
    { int32_t iv;
      for (iv = 0; iv < maxv; iv++)
        { depends[iv] = (bool_t *) malloc(maxe * sizeof(bool_t));
          reduce[iv] = (int32_t *) malloc(maxe * sizeof(int32_t));
        }
    }
    
    aa_compute_value_scopes
      ( p, maxv, 
        /* out */
        defined, last_used, result, &nv
      );
    aa_compute_epsilon_tables
      ( p, defined, last_used, nv, 
        arg_depends, res_depends, shared_ne, maxe,
        /* out */
        erreps, depends, reduce, &ne
      );
    aa_generate_code
      ( p, proc_name, c_file, 
        arg_depends, res_depends, shared_ne,
        defined, last_used, result, nv,
        erreps, depends, reduce, ne
      );
  }
          
void aa_compute_value_scopes
  ( pcode_proc_t *p, 
    int32_t maxv,
    /* out */
    int32_t *defined, 
    int32_t *last_used,
    int32_t *result,
    int32_t *nv
  )
  { 
    /* Initialize tables: */
    (*nv) = 0;
    { int32_t iv;
      for (iv = 0; iv < maxv; iv++) 
        { defined[iv] = -1; last_used[iv] = -1; }
    }

    /* Process pcode: */
    {
      int32_t *stack = (int32_t *) malloc(p->nstack * sizeof(int32_t));
      int32_t *reg = (int32_t *) malloc(p->nregs * sizeof(int32_t));

      int32_t top = -1;
      int32_t iop = 1;
      pcode_op_t op;
      int32_t arg1, arg2;
      int32_t *res1;

      /* Fill stack and regs with garbage: */
      { int32_t i;
	for (i = 0; i < p->nstack; i++) { stack[i] = -1; }
	for (i = 0; i < p->nregs; i++) { reg[i] = -1; }
      }

      /* Process the "given" instruction: */
      affirm(
	p->code[0].op == pcode_op_given, 
	"pcode does not start with GIVEN"
      );
      { int32_t ig;
	for (ig = 0; ig < p->nin; ig++) 
	  { defined[ig] = 0; reg[ig] = ig; }
	(*nv) = p->nin;
      }

      /* Process the body instructions: */
      for (iop = 1; iop < p->nops - 1; iop++)
	{ 
	  op = p->code[iop].op;
	  switch (op)
	    {
	      case pcode_op_given:
		affirm(FALSE, "GIVEN inside body");
		break;
	      case pcode_op_load:
		arg1 = reg[p->code[iop].arg];
                res1 = &(stack[top+1]);
                *res1 = arg1;
		top++;
		break;
	      case pcode_op_store:
		arg1 = stack[top];
                reg[p->code[iop].arg] = arg1; 
		top--;
		break;
	      case pcode_op_const:
                res1 = &(stack[top+1]);
		*res1 = *nv; defined[*nv] = iop; (*nv)++;
		top++;
		break;
	      case pcode_op_return:
		affirm(FALSE, "RETURN inside body");
		break;
	      case pcode_op_inv:
	      case pcode_op_neg:
	      case pcode_op_sqr:
	      case pcode_op_sqrt:
	      case pcode_op_abs:
		arg1 = stack[top];
                res1 = &(stack[top]);
                affirm(arg1 != -1, "undef operand");
		last_used[arg1] = iop;
		*res1 = *nv; defined[*nv] = iop; (*nv)++;
		break;
	      case pcode_op_add:
	      case pcode_op_sub:
	      case pcode_op_mul:
	      case pcode_op_div:
	      case pcode_op_max:
	      case pcode_op_min:
                arg1 = stack[top-1]; 
                arg2 = stack[top];
                res1 = &(stack[top-1]);
		affirm(arg1 != -1,"undef operand");
		affirm(arg2 != -1, "undef operand");
		last_used[arg1] = iop;
		last_used[arg2] = iop;
		*res1 = *nv; defined[*nv] = iop; (*nv)++;
		top--;
		break;
	      default:
		affirm(FALSE, "bad op code");
	    }
	}

      /* Process the RETURN instruction: */
      iop = p->nops -1;
      affirm(
	p->code[iop].op == pcode_op_return, 
	"pcode does not end with RETURN"
      );
      { int32_t ir;
	for (ir = 0; ir < p->nout; ir++) 
	  { int32_t iv = stack[ir];
	    affirm(iv != -1,"undef returned value");
	    result[ir] = iv;
	    last_used[iv] = iop; 
	  }
      }
    }
  }
  
void aa_compute_epsilon_tables
  ( pcode_proc_t *p, 
    int32_t *defined, 
    int32_t *last_used, 
    int32_t nv, 
    bool_t **arg_depends,
    bool_t **res_depends,
    int32_t shared_ne,
    int32_t maxe, 
    /* out */
    int32_t *erreps,
    bool_t **depends, 
    int32_t **reduce, 
    int32_t *ne
  )
  {
    int32_t *dep_count; /* dep_count[ie] is the number of live values 
                       that currently depend on epsilon[ie]
                       (recomputed for each instruction). */
    
    /* Initialize tables: */
    (*ne) = 0;
    dep_count = (int32_t *) malloc(maxe * sizeof(int32_t));
    { int32_t iv, ie;
      for (iv = 0; iv < nv; iv++)
        { erreps[iv] = -1;
          for (ie = 0; ie < maxe; ie++)
            { depends[iv][ie] = FALSE; reduce[iv][ie] = -1; }
        }
      for (ie = 0; ie < maxe; ie++) { dep_count[ie] = 0; }
    }
    
    { 
      int32_t *stack = (int32_t *) malloc(p->nstack * sizeof(int32_t));
      int32_t *reg = (int32_t *) malloc(p->nregs * sizeof(int32_t));

      int32_t iv = p->nin;
      int32_t iop, je;
      int32_t top = -1;
      pcode_op_t op;
      int32_t arg1, arg2;
      int32_t *res1;

      /* Fill stack and regs with garbage: */
      { int32_t i;
	for (i = 0; i < p->nstack; i++) { stack[i] = -1; }
	for (i = 0; i < p->nregs; i++) { reg[i] = -1; }
      }

      /* Process the "given" instruction: */
      affirm(
	p->code[0].op == pcode_op_given, 
	"pcode does not start with GIVEN"
      );
      { int32_t ig, ie;
	for (ig = 0; ig < p->nin; ig++) 
	  { iv = ig;
	    for (ie = 0; ie < shared_ne; ie++)
	      { depends[iv][ie] = arg_depends[iv][ie]; 
		if (depends[iv][ie]) { dep_count[ie]++; }
	      }
	    /* Include "specific" epsilons: */
	    ie = shared_ne + iv;
	    depends[iv][ie] = TRUE;
            erreps[iv] = ie;
	    dep_count[ie]++;
	    reg[ig] = iv; 
	  }
	(*ne) = shared_ne + p->nin;
	iv = p->nin;
      }

      /* Process the body instructions: */
      for (iop = 1; iop < p->nops - 1; iop++)
	{ 
	  op = p->code[iop].op;
	  switch (op)
	    {
	      case pcode_op_given:
		affirm(FALSE, "GIVEN inside body");
		break;
	      case pcode_op_load:
		{ arg1 = reg[p->code[iop].arg];
                  res1 = &(stack[top+1]);
                  *res1 = arg1;
		  top++;
                }
		break;
	      case pcode_op_store:
		{ arg1 = stack[top];
                  reg[p->code[iop].arg] = arg1; 
		  top--;
                }
		break;
	      case pcode_op_const:
		{ res1 = &(stack[top+1]);
		  affirm(defined[iv] == iop, "defined[iv] != iop");
                  { Interval ia = ia_int_const(p->code[iop].arg);
                    if (ia.lo != ia.hi)
                      { int32_t ie = (*ne);
                        depends[iv][ie] = TRUE; erreps[iv] = ie; (*ne)++;
                        dep_count[ie]++;
                      }
                  }
                  *res1 = iv; 
                  iv++;
                  top++;
                }
		break;
	      case pcode_op_return:
		affirm(FALSE, "RETURN inside body");
		break;
	      case pcode_op_neg:
                { int32_t ie; bool_t neweps = FALSE;
                  arg1 = stack[top];
                  res1 = &(stack[top]);
                  affirm(arg1 != -1, "undef operand");
                  affirm(last_used[arg1] >= iop, "last_used bug");
                  affirm(defined[iv] == iop, "defined[iv] != iop"); 
                  if (last_used[arg1] == iop) 
                    aa_update_dep_count(depends[arg1], *ne, dep_count, -1);
                  /* Compute epsilon reduction and dependency data for value $iv$: */
                  /* Note that NEG is exact, so there is no new epsilon unless */
                  /*   we are reducing: */
                  ie = (*ne);
                  for (je = shared_ne; je < (*ne); je++) 
                    { if (depends[arg1][je])
                        { if (dep_count[je] == 0)
                            { reduce[iv][je] = ie; neweps = TRUE; }
                          else
                            { depends[iv][je] = TRUE; }
                        }
                    }
                  if (neweps) 
                    { depends[iv][ie] = TRUE; erreps[iv] = ie; (*ne)++; }
                  if (last_used[iv] > iop) 
                    aa_update_dep_count(depends[iv], (*ne), dep_count, +1);
                  *res1 = iv;
                  iv++;
                }
		break;
	      case pcode_op_inv:
	      case pcode_op_sqr:
	      case pcode_op_sqrt:
	      case pcode_op_abs:
		{ int32_t ie;
                  arg1 = stack[top];
                  res1 = &(stack[top]);
                  affirm(arg1 != -1, "undef operand");
                  affirm(last_used[arg1] >= iop, "last_used bug");
                  if (last_used[arg1] == iop) 
                    aa_update_dep_count(depends[arg1], (*ne), dep_count, -1);
                  affirm(defined[iv] == iop, "defined[iv] != iop"); 
                  /* Compute epsilon reduction and dependency data for value $iv$: */
                  ie = (*ne);
                  for (je = shared_ne; je < (*ne); je++) 
                    { if (depends[arg1][je])
                        { if (dep_count[je] == 0)
                            { reduce[iv][je] = ie; }
                          else
                            { depends[iv][je] = TRUE; }
                        }
                    }
                  /* These operations always generate a new epsilon: */
                  depends[iv][ie] = TRUE; erreps[iv] = ie; (*ne)++;
                  if (last_used[iv] > iop) 
                    aa_update_dep_count(depends[iv], (*ne), dep_count, +1);
                  *res1 = iv; 
                  iv++; 
                }
		break;
	      case pcode_op_add:
	      case pcode_op_sub:
	      case pcode_op_mul:
	      case pcode_op_div:
	      case pcode_op_max:
	      case pcode_op_min:
		{ int32_t ie;
                  arg1 = stack[top-1]; 
                  arg2 = stack[top];
                  res1 = &(stack[top-1]);
                  affirm(arg1 != -1,"undef operand");
                  affirm(last_used[arg1] >= iop, "last_used bug");
                  if (last_used[arg1] == iop) 
                    aa_update_dep_count(depends[arg1], (*ne), dep_count, -1);
                  affirm(arg2 != -1, "undef operand");
                  affirm(last_used[arg2] >= iop, "last_used bug");
                  if (last_used[arg2] == iop) 
                    aa_update_dep_count(depends[arg2], (*ne), dep_count, -1);
                  affirm(defined[iv] == iop, "defined[iv] != iop"); 
                  /* Compute epsilon reduction and dependency data for value $iv$: */
                  ie = (*ne);
                  for (je = shared_ne; je < (*ne); je++) 
                    { if (depends[arg1][je] | depends[arg2][je])
                        { if (dep_count[je] == 0)
                            { reduce[iv][je] = ie; }
                          else
                            { depends[iv][je] = TRUE; }
                        }
                    }
                  depends[iv][ie] = TRUE; erreps[iv] = ie; (*ne)++;
                  if (last_used[iv] > iop) 
                    aa_update_dep_count(depends[iv], (*ne), dep_count, +1);
                  *res1 = iv; 
                  iv++;
                  top--;
                }
		break;
	      default:
		affirm(FALSE, "bad op code");
	    }
	}

      /* Process the RETURN instruction: */
      iop = p->nops -1;
      affirm(
	p->code[iop].op == pcode_op_return, 
	"pcode does not end with RETURN"
      );
      { int32_t ir;
	for (ir = 0; ir < p->nout; ir++) 
	  { int32_t iv = stack[ir];
	    affirm(iv != -1, "undef returned value");
	    affirm(last_used[iv] == iop, "last_used bug");
	  }
      }
    }

  }
  
void aa_update_dep_count
  ( bool_t *var_depends,  /* One row of the $depends$ table */
    int32_t ne, 
    int32_t *dep_count, 
    int32_t increment
  )
  {
    int32_t ie;
    for (ie = 0; ie < ne; ie++)
      { if (! (var_depends[ie]))
          { dep_count[ie] += increment; }
      }
  }
      
void aa_generate_code
  ( pcode_proc_t *p, 
    char *proc_name, 
    FILE *c_file, 
    bool_t **arg_depends, 
    bool_t **res_depends, 
    int32_t shared_ne,
    int32_t *defined, 
    int32_t *last_used, 
    int32_t *result, 
    int32_t nv,
    int32_t *erreps,
    bool_t **depends, 
    int32_t **reduce, 
    int32_t ne
  )
  {
    char **v_name = aa_compute_value_names(nv);
    
    aa_generate_proc_header
      ( p, proc_name, c_file, 
        result, arg_depends, res_depends, shared_ne,
        v_name
      );
    aa_generate_local_decls
      ( p, proc_name, c_file, 
        nv,
        depends, ne,
        v_name
      );
    aa_generate_body
      ( p, proc_name, c_file,  
        defined, last_used, nv,
        erreps, depends, reduce, ne,
        v_name
      );
    aa_generate_proc_trailer
      ( p, proc_name, c_file, 
        nv, result, res_depends, shared_ne, 
        depends, ne,
        v_name
      );
  }

char **aa_compute_value_names(int32_t nv)
  {
    char **v_name = (char**) malloc(nv * sizeof(char*));
    int32_t iv;
    for (iv = 0; iv < nv; iv++)
      { char *p = ((char*) malloc(4 * sizeof(char))) + 3;
	int32_t k = iv;
	*p = '\000';
	p--; *p = (char)('a' + (k % 26)); k = k / 26;
	if (k > 0) { p--; *p = (char)('a' + (k % 26) - 1); k = k / 26; }
	if (k > 0) { p--; *p = (char)('a' + (k % 26) - 1); k = k / 26; }
	affirm(k == 0, "too many variables!");
	v_name[iv] = p;
      }
    return(v_name);
  }

void aa_generate_proc_header
  ( pcode_proc_t *p,
    char *proc_name,
    FILE *c_file, 
    int32_t *result,
    bool_t **arg_depends, 
    bool_t **res_depends, 
    int32_t shared_ne,
    char **v_name
  )
  { int32_t iv, ir, ie;
    
    /* Input parameters: */
    fprintf(c_file, "void %s (\n", proc_name);
    for (iv = 0; iv < p->nin; iv++)
      { fprintf(c_file, "    Float %s_ctr,\n", v_name[iv]); 
        for (ie = 0; ie < shared_ne; ie++)
          { if ( arg_depends[iv][ie] )
              { fprintf(c_file, "    Float %s_%d,\n", v_name[iv], ie); }
          }
        fprintf(c_file, "    Float %s_%d,\n", v_name[iv], shared_ne + iv);
      }
    
    /* Output parameters: */
    for (ir = 0; ir < p->nout; ir++)
      { iv = result[ir];
        fprintf(c_file, "    Float *r%s_ctr,\n", v_name[ir]); 
        for (ie = 0; ie < shared_ne; ie++)
          { if ( res_depends[iv][ie] )
              { fprintf(c_file, "    Float *r%s_%d,\n", v_name[ir], ie); }
          }
        fprintf(c_file, "    Float *r%s_%d,\n", v_name[ir], shared_ne + p->nin + ir);
      }
    
    fprintf(c_file, "    int32_t *nops\n");

    fprintf(c_file, "  )\n");
    
    fprintf(c_file, "  {\n");
    
  }
  
void aa_generate_local_decls
  ( pcode_proc_t *p,
    char *proc_name,
    FILE *c_file, 
    int32_t nv,
    bool_t **depends, 
    int32_t ne,
    char **v_name
  )
  { int32_t iv;
    int32_t ie;
    fprintf(c_file, "    Float alpha, beta, zeta, gamma, delta, tmp;\n");
    for (iv = p->nin; iv < nv; iv++)
      { fprintf(c_file, "    Float %s_ctr, %s_rho;", v_name[iv], v_name[iv]); 
        fprintf(c_file, "    Interval %s_range;\n", v_name[iv]); 
        for (ie = 0; ie < ne; ie++)
          { if ( depends[iv][ie] )
              { fprintf(c_file, "    Float %s_%d;\n", v_name[iv], ie); }
          }
      }
  }

void aa_generate_body
  ( pcode_proc_t *p,
    char *proc_name,
    FILE *c_file, 
    int32_t *defined, 
    int32_t *last_used, 
    int32_t nv, 
    int32_t *erreps,
    bool_t **depends,
    int32_t **reduce,
    int32_t ne,
    char **v_name
  )
  { int32_t *stack = (int32_t *) malloc(p->nstack * sizeof(int32_t));
    int32_t *reg = (int32_t *) malloc(p->nregs * sizeof(int32_t));
    int32_t top = -1;
    int32_t iop;
    int32_t iv, ie;
    pcode_op_t op;
    int32_t arg1, arg2;
    int32_t *res1;

    /* Fill stack and regs with garbage: */
    { int32_t i;
      for (i = 0; i < p->nstack; i++) { stack[i] = -1; }
      for (i = 0; i < p->nregs; i++) { reg[i] = -1; }
    }

    /* Process the "given" instruction: */
    affirm
      ( p->code[0].op == pcode_op_given, 
        "pcode does not start with GIVEN"
      );
    { int32_t ig;
      for (ig = 0; ig < p->nin; ig++) 
	{ int32_t iv = ig;
	  reg[ig] = iv;
	}
      iv = p->nin;
    }

    /* Process the body instructions: */
    for (iop = 1; iop < p->nops - 1; iop++)
      { op = p->code[iop].op;
	switch (op)
	  {
	    case pcode_op_given:
	      affirm(FALSE, "GIVEN inside body");
	      break;
	    case pcode_op_load:
	      arg1 = reg[p->code[iop].arg];
	      res1 = &(stack[top+1]);
	      *res1 = arg1;
	      top++;
	      break;
	    case pcode_op_store:
	      arg1 = stack[top];
	      reg[p->code[iop].arg] = arg1; 
	      top--;
	      break;
	    case pcode_op_const:
	      res1 = &(stack[top+1]);
	      affirm(defined[iv] == iop, "defined[iv] != iop");
              ie = erreps[iv];
	      if (ie != -1) 
                { affirm(depends[iv][ie], "depends bug"); }
	      aa_generate_op_const(c_file, v_name[iv], p->code[iop].arg, ie);
              *res1 = iv; 
	      iv++;
	      top++;
	      break;
	    case pcode_op_return:
	      affirm(FALSE, "RETURN inside body");
	      break;
	    case pcode_op_neg:
            case pcode_op_inv:
            case pcode_op_sqr:
            case pcode_op_sqrt:
            case pcode_op_abs:
	      arg1 = stack[top];
	      res1 = &(stack[top]);
	      affirm(arg1 != -1, "undef operand");
	      affirm(last_used[arg1] >= iop, "last_used bug");
	      affirm(defined[iv] == iop, "defined[iv] != iop");
              ie = erreps[iv];
	      if (ie != -1) 
                { affirm(depends[iv][ie], "depends bug"); }
	      aa_generate_unary_op
                ( c_file, 
                  v_name[iv], op, v_name[arg1], ie, 
                  depends[arg1], reduce[iv], ne
                );
	      *res1 = iv; 
	      iv++;
	      break;
	    case pcode_op_add:
            case pcode_op_sub:
            case pcode_op_mul:
            case pcode_op_div:
            case pcode_op_max:
            case pcode_op_min:
	      arg1 = stack[top-1]; 
	      arg2 = stack[top];
	      res1 = &(stack[top-1]);
	      affirm(arg1 != -1, "undef operand");
	      affirm(last_used[arg1] >= iop, "last_used bug");
	      affirm(arg2 != -1, "undef operand");
	      affirm(last_used[arg2] >= iop, "last_used bug");
	      affirm(defined[iv] == iop, "defined[iv] != iop"); 
	      ie = erreps[iv];
	      if (ie != -1) 
                { affirm(depends[iv][ie], "depends bug"); }
	      aa_generate_binary_op
                (c_file, 
                  v_name[iv], op, v_name[arg1], v_name[arg2],
                  ie,
                  depends[arg1], depends[arg2], reduce[iv], ne
                );
	      *res1 = iv;
	      iv++;
	      top--;
	      break;
	    default:
	      affirm(FALSE, "bad op code");
	  }
      }

    /* Process the RETURN instruction: */
    iop = p->nops - 1;
    affirm
    ( p->code[iop].op == pcode_op_return, 
      "pcode does not end with RETURN"
    );
  }
    
void aa_generate_proc_trailer
  ( pcode_proc_t *p,
    char *proc_name,
    FILE *c_file,
    int32_t nv,
    int32_t *result,
    bool_t **res_depends, 
    int32_t shared_ne,
    bool_t **depends,
    int32_t ne,
    char **v_name
  )
  { int32_t ir, ie, je;
    ie = shared_ne + nv;
    fprintf(c_file, "    ROUND_UP;\n");
    for (ir = 0; ir < p->nout; ir++) 
      { int32_t iv = result[ir];
        affirm(iv != -1, "undef returned value");
        fprintf(c_file, "    *r%s_ctr = %s_ctr;\n", v_name[ir], v_name[iv]);
        fprintf(c_file, "    *r%s_%d = Zero;\n", v_name[ir], ie);
        for (je = 0; je < ne; je++)
          if (depends[iv][je])
            { if (je < shared_ne && res_depends[ir][je])
                fprintf
                  ( c_file, "    *r%s_%d = %s_%d;\n", 
                    v_name[ir], je, v_name[iv], je
                  );
              else
                fprintf
                  ( c_file, "    *r%s_%d += ABS(%s_%d);\n",
                    v_name[ir], ie, v_name[iv], je
                  );
            }
        ie++;
      }
    fprintf(c_file, "  }\n\n");
  }

void aa_generate_op_const
  ( FILE *c_file, 
    char *name_v,   /* Name of result. */
    int32_t arg,        /* The integer to convert. */
    int32_t ie          /* Id of error epsilon, or -1 if exact. */
  )
  { if (ie != -1)
      { /* Large integer, we have to include a rounding term. */
        fprintf(c_file, "    %s_range = ia_from_int(%d);\n", name_v, arg);
        /* The code below should neither oferflow nor underflow: */
        fprintf(c_file, "    ROUND_NEAR;\n");
        fprintf(c_file, "    %s_ctr = (Half * %s_range.lo) + (Half * %s_range.hi);\n",
           name_v, name_v, name_v
        );
        /* At this point, %s_ctr should be inside %s_range. */
        /* The code below should neither oferflow nor underflow: */
        fprintf(c_file, "    ROUND_UP;\n");
        fprintf(c_file, "    { Float r1 = %s_range.hi - %s_ctr;\n", name_v, name_v);
        fprintf(c_file, "      Float r2 = %s_ctr - %s_range.lo;\n", name_v, name_v);
        fprintf(c_file, "      %s_%d = FMAX(r1, r2);\n", name_v, ie);
        fprintf(c_file, "    }\n");
        fprintf(c_file, "    %s_rho = %s_%d;\n", name_v, name_v, ie);
      }
    else
      { /* Small integer, no epsilons to worry about: */
        fprintf(c_file, "    %s_ctr = (Float) %d;\n", name_v, arg);
        fprintf(c_file, "    %s_rho = Zero;\n", name_v);
        fprintf(c_file, "    %s_range = (Interval) {%s_ctr, %s_ctr};\n", 
          name_v, name_v, name_v
        );
      }
  }

void aa_generate_unary_op
  ( FILE *c_file, 
    char *name_v,   /* Name of result. */
    pcode_op_t op,
    char *name_x,       /* Name of argument. */
    int32_t ie,             /* Id of epsilon of error term, or -1 if exact. */
    bool_t *depends_x,    /* Epsilons appearing in x. */
    int32_t *reduce_v,      /* Epsilon reduction table for result. */
    int32_t ne              /* Number of epsilons. */
  )
  { switch(op)
      { case pcode_op_neg: 
          aa_generate_op_neg(c_file, name_v, name_x, ie, depends_x, reduce_v, ne);
          break;
        case pcode_op_sqrt: 
          aa_generate_op_sqrt(c_file, name_v, name_x, ie, depends_x, reduce_v, ne);
          break;
        case pcode_op_abs: 
        case pcode_op_sqr: 
        case pcode_op_inv: 
          affirm(FALSE, "unimplemented op code");
        default: 
          affirm(FALSE, "bad op code");
      }
  }

void aa_generate_binary_op
  ( FILE *c_file, 
    char *name_v,   /* Name of result. */
    pcode_op_t op,
    char *name_x,       /* Name of first argument. */
    char *name_y,       /* Name of second argument. */
    int32_t ie,             /* Id of epsilon of error term, or -1 if exact. */
    bool_t *depends_x,    /* Epsilons appearing in x. */
    bool_t *depends_y,    /* Epsilons appearing in y. */
    int32_t *reduce_v,      /* Epsilon reduction table for result. */
    int32_t ne              /* Number of epsilons. */
  )
  { 
    switch(op)
      { case pcode_op_add: 
          aa_generate_op_add
            ( c_file, 
              name_v, name_x, name_y, ie, depends_x, depends_y, reduce_v, ne
            );
          break;
        case pcode_op_mul: 
          aa_generate_op_mul
            ( c_file, 
              name_v, name_x, name_y, ie, depends_x, depends_y, reduce_v, ne
            );
          break;
        case pcode_op_max: 
        case pcode_op_min: 
        case pcode_op_div: 
          affirm(FALSE, "unimplemented op code");
        default: 
          affirm(FALSE, "bad op code");
      }
  }

void aa_generate_op_neg
  ( FILE *c_file, 
    char *name_v,       /* Name of result. */
    char *name_x,       /* Name of argument. */
    int32_t ie,             /* Id of epsilon of error term, or -1 if exact. */
    bool_t *depends_x,    /* Epsilons appearing in x. */
    int32_t *reduce_v,      /* Epsilon reduction table for result. */
    int32_t ne              /* Number of epsilons. */
  )
  { int32_t je;
    fprintf(c_file, "    %s_ctr = - %s_ctr;\n", name_v, name_x);
    if (ie != -1)
      { /* Result is inexact, must be because of reductions: */
	fprintf(c_file, "    %s_%d = Zero;\n",  name_v, ie);
	fprintf(c_file, "    ROUND_UP;\n");
	for (je = 0; je < ne; je++)
	  { if (depends_x[je]) 
	      { if (reduce_v[je] == -1)
		  fprintf(c_file, "    %s_%d = - %s_%d;\n", name_v, je, name_x, je);
		else
		  { affirm(reduce_v[je] == ie, "reduce inconsistent!");
		    fprintf(c_file, "    %s_%d += FABS(%s_%d);\n", 
		      name_v, ie, name_x, je
		    );
		  }
	      }
	  }
      }
    else
      { /* No reductions, result is exact: */
	for (je = 0; je < ne; je++)
	  { if (depends_x[je]) 
	      { affirm(reduce_v[je] == -1, "reduce inconsistent!");
		fprintf(c_file, "    %s_%d = - %s_%d;\n", name_v, je, name_x, je);
	      }
	  }
      }
    fprintf(c_file, "    %s_range = ia_neg(%s_range);\n", name_v, name_x);
  }

void aa_generate_op_sqrt
  ( FILE *c_file, 
    char *name_v,       /* Name of result. */
    char *name_x,       /* Name of argument x. */
    int32_t ie,             /* Id of epsilon of error term. */
    bool_t *depends_x,    /* Epsilons appearing in the x. */
    int32_t *reduce_v,      /* Epsilon reduction table for result. */
    int32_t ne              /* Number of epsilons. */
  )
  { int32_t je;
    fprintf(c_file, "    cheb_sqrt(%s_range, &alpha, &zeta, &gamma, &delta);\n", name_x);
    fprintf(c_file, "    %s_%d = delta;\n",  name_v, ie);
    fprintf(c_file, "    flt_affine(%s_ctr, alpha, zeta, gamma, &%s_ctr, &%s_%d);\n", 
      name_x, name_v, name_v, ie
    );
    for (je = 0; je < ne; je++)
      { if (depends_x[je]) 
	  { if (reduce_v[je] == -1)
	      { fprintf(c_file, "    flt_scale(%s_%d, alpha, zeta, &%s_%d, &%s_%d);\n",
		  name_x, je, name_v, je, name_v, ie);
	      }
	    else
	      { affirm(reduce_v[je] == ie, "reduce inconsistent!");
		fprintf(c_file, "    flt_scale(%s_%d, alpha, zeta, &tmp, &%s_%d);\n",
		  name_x, je, name_v, ie
                );
		fprintf(c_file, "    %s_%d += FABS(tmp);\n", name_v, ie);
	      }
	  }
      }
    fprintf(c_file, "    %s_range = ia_sqrt(%s_range);\n", name_v, name_x);
  }

void aa_generate_op_add
  ( FILE *c_file, 
    char *name_v,       /* Name of result. */
    char *name_x,       /* Name of argument x. */
    char *name_y,       /* Name of argument y. */
    int32_t ie,             /* Id of epsilon of error term. */
    bool_t *depends_x,    /* Epsilons appearing in the x. */
    bool_t *depends_y,    /* Epsilons appearing in the y. */
    int32_t *reduce_v,      /* Epsilon reduction table for result. */
    int32_t ne              /* Number of epsilons. */
  )
  { int32_t je;
    fprintf(c_file, "    %s_%d = Zero;\n",  name_v, ie);
    fprintf(c_file, "    %s_rho = Zero;\n",  name_v);
    fprintf(c_file, "    flt_add(%s_ctr, %s_ctr, &%s_ctr, &%s_%d);\n", 
      name_x, name_y, name_v, name_v, ie
    );
    for (je = 0; je < ne; je++)
      { if (depends_x[je] | depends_y[je]) 
	  { if (reduce_v[je] == -1)
	      { fprintf(c_file, "    flt_add(%s_%d, %s_%d, &%s_%d, &%s_%d);\n",
		  name_x, je, name_y, je, name_v, je, name_v, ie
                );
		fprintf(c_file, "    ROUND_UP; %s_rho += FABS(%s_%d);\n", name_v, name_v, je);
	      }
	    else
	      { affirm(reduce_v[je] == ie, "reduce inconsistent!");
		fprintf(c_file, "    flt_add(%s_%d, %s_%d, &tmp, &%s_%d);\n",
		  name_x, je, name_y, je, name_v, ie
                );
		fprintf(c_file, "    ROUND_UP; %s_%d += FABS(tmp);\n", name_v, ie);
	      }
	  }
      }
    fprintf(c_file, "    %s_rho += FABS(%s_%d);\n", name_v, name_v, ie);
    fprintf(c_file, "    %s_range = ia_meet(ia_add(%s_range, %s_range), ia_const(%s_ctr, %s_rho))\n",
      name_v, name_x, name_y, name_v, name_v
    );
  }

void aa_generate_op_mul
  ( FILE *c_file, 
    char *name_v,       /* Name of result. */
    char *name_x,       /* Name of argument x. */
    char *name_y,       /* Name of argument y. */
    int32_t ie,             /* Id of epsilon of error term. */
    bool_t *depends_x,    /* Epsilons appearing in the x. */
    bool_t *depends_y,    /* Epsilons appearing in the y. */
    int32_t *reduce_v,      /* Epsilon reduction table for result. */
    int32_t ne              /* Number of epsilons. */
  )
  { fprintf(c_file, "    %s_%d = delta;\n",  name_v, ie);
    fprintf(c_file, "    %s_rho = Zero;\n",  name_v);
    fprintf(c_file, "    flt_mul(%s_ctr, %s_ctr, &%s_ctr, &%s_%d);\n", 
      name_x, name_y, name_v, name_v, ie
    );
    
    affirm(FALSE, "unimplemented");
    
  }
