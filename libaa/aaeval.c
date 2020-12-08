/* See "aaeval.h" */
/* Last edited on 2002-12-29 16:04:05 by stolfi */

#include <aaeval.h>
#include <aa.h>
#include <pcode.h>
#include <flt.h>
#include <ia.h>
#include <affirm.h>

void aa_eval (AAP reg[], AAP stack[], pcode_instr_t fcode[])
  { int top = -1; 
    AAP *x, *y, *z;
    int iop = 0;
    pcode_op_t op;

    while ((op = fcode[iop].op)  != pcode_op_return)
      {
        x = &(stack[top-1]); y = x + 1; z = y + 1;
        switch (op)
	  {
	    case pcode_op_given:
              break;
            case pcode_op_load:
	      *z = reg[fcode[iop].arg];
	      top++;
	      break;
	    case pcode_op_const:
	      *z = aa_int_const(fcode[iop].arg);
	      top++;
	      break;
	    case pcode_op_store:
	      reg[fcode[iop].arg] = *y; 
	      top--;
	      break;
	    case pcode_op_return:
	      fatalerror("aa_eval: unexpected return!");
	      break;
	    case pcode_op_add:
	      *x  = aa_add(*x, *y);
	      top--;
	      break;
	    case pcode_op_sub:
	      *x  = aa_sub(*x, *y);
	      top--;
	      break;
	    case pcode_op_neg:
	      *y  = aa_neg(*y);
	      break;
	    case pcode_op_inv:
	      *y  = aa_inv(*y);
	      break;
	    case pcode_op_mul:
	      *x  = aa_mul(*x, *y);
	      top--;
	      break;
	    case pcode_op_div:
	      *x  = aa_div(*x, *y);
	      top--; break;
	    case pcode_op_sqr:
	      *y  = aa_sqr(*y);
	      break;
	    case pcode_op_sqrt:
	      *y  = aa_sqrt(*y);
	      break;
	    case pcode_op_abs:
	      *y  = aa_abs(*y);
	      break;
	    case pcode_op_max:
	      *x  = aa_max((*x), (*y));
	      top--;
	      break;
	    case pcode_op_min:
	      *x  = aa_min((*x), (*y));
	      top--;
	      break;
            default:
              fatalerror("aa_eval: bad op code");
	  }
        iop++;
      }
  }

void aa_eval_diff (DiffAAP reg[], DiffAAP stack[], pcode_instr_t fcode[])
  {
    int top = -1; 
    DiffAAP *x, *y, *z;
    int iop = 0;
    pcode_op_t op;

    while ((op = fcode[iop].op)  != pcode_op_return)
      {
        x = &(stack[top-1]); y = x + 1; z = y + 1;
        switch (op)
	  {
	    case pcode_op_given:
              break;
            case pcode_op_load:
	      *z = reg[fcode[iop].arg];
	      top++;
	      break;
	    case pcode_op_const:
	      z->f  = aa_int_const(fcode[iop].arg);
              z->df = aa_zero();
	      top++;
	      break;
	    case pcode_op_store:
	      reg[fcode[iop].arg] = *y; 
	      top--;
	      break;
	    case pcode_op_return:
	      fatalerror("aa_eval: unexpected return!");
	      break;
	    case pcode_op_add:
	      x->df = aa_add(x->df,y->df);
	      x->f  = aa_add(x->f, y->f);
	      top--;
	      break;
	    case pcode_op_sub:
	      x->df = aa_sub(x->df,y->df);
	      x->f  = aa_sub(x->f, y->f);
	      top--;
	      break;
	    case pcode_op_neg:
	      y->df = aa_neg(y->df);
	      y->f  = aa_neg(y->f);
	      break;
	    case pcode_op_inv:
	      y->df = aa_neg(aa_div(y->df, aa_sqr(x->f)));
              y->f  = aa_inv(y->f);
	      break;
	    case pcode_op_mul:
	      x->df = aa_add(aa_mul(x->f, y->df), aa_mul(x->df,y->f));
	      x->f  = aa_mul(x->f, y->f);
	      top--;
	      break;
	    case pcode_op_div:
	      x->df = aa_div(aa_sub(aa_mul(x->df, y->f), aa_mul(x->f, y->df)),
			     aa_sqr(y->f));
	      x->f  = aa_div(x->f, y->f);
	      top--; break;
	    case pcode_op_sqr:
	      y->df = aa_mul(aa_scale(y->f,2,1), y->df);
	      y->f  = aa_sqr(y->f);
	      break;
	    case pcode_op_sqrt:
	      { AAP r = aa_sqrt(y->f);
	        y->df = aa_div(y->df, aa_scale(r,1,2));
                y->f  = r;
              }
	      break;
	    case pcode_op_abs:
	      { Interval yr = aa_range(y->f);
		if ((yr.lo <= Zero) && (yr.hi >= Zero))
		  { y->df = aa_join(y->df, aa_neg(y->df)); }
		else
		  { y->df = (yr.lo > Zero ? y->df : aa_neg(y->df)); }
              }
              y->f  = aa_abs(y->f);
	      break;
	    case pcode_op_max:
	      { Interval xr = aa_range(x->f);
		Interval yr = aa_range(y->f);
		if ((yr.lo <= xr.hi) && (yr.hi >= xr.lo))
		  { x->df = aa_join(x->df, y->df); }
		else
		  { x->df = (xr.lo > yr.lo ? x->df : y->df); }
              }
              x->f  = aa_max(x->f, y->f);
	      top--;
	      break;
	    case pcode_op_min:
	      { Interval xr = aa_range(x->f);
		Interval yr = aa_range(y->f);
		if ((yr.lo <= xr.hi) && (yr.hi >= xr.lo))
		  { x->df = aa_join(x->df, y->df); }
		else
		  { x->df = (xr.lo > yr.lo ? y->df : x->df); }
              }
              x->f  = aa_min(x->f, y->f);
	      top--;
	      break;
            default:
              fatalerror("aa_eval_diff: bad op code");
	  }
        iop++;
      }
  }
