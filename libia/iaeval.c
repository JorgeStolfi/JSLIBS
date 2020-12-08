/* See "iaeval.h" */

#include <iaeval.h>
#include <ia.h>
#include <pcode.h>
#include <flt.h>
#include <affirm.h>

void ia_eval (Interval reg[], Interval stack[], pcode_instr_t fcode[])
  {
    int top = -1; 
    Interval *x, *y, *z;
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
	      *z = ia_int_const(fcode[iop].arg);
	      top++;
	      break;
	    case pcode_op_store:
	      reg[fcode[iop].arg] = *y; 
	      top--;
	      break;
	    case pcode_op_return:
	      fatalerror("ia_eval: unexpected return!");
	      break;
	    case pcode_op_add:
	      *x  = ia_add(*x, *y);
	      top--;
	      break;
	    case pcode_op_sub:
	      *x  = ia_sub(*x, *y);
	      top--;
	      break;
	    case pcode_op_neg:
	      *y  = ia_neg(*y);
	      break;
	    case pcode_op_inv:
	      *y  = ia_inv(*y);
	      break;
	    case pcode_op_mul:
	      *x  = ia_mul(*x, *y);
	      top--;
	      break;
	    case pcode_op_div:
	      *x  = ia_div(*x, *y);
	      top--; break;
	    case pcode_op_sqr:
	      *y  = ia_sqr(*y);
	      break;
	    case pcode_op_sqrt:
	      *y  = ia_sqrt(*y);
	      break;
	    case pcode_op_abs:
	      *y  = ia_abs(*y);
	      break;
	    case pcode_op_max:
	      *x  = ia_max((*x), (*y));
	      top--;
	      break;
	    case pcode_op_min:
	      *x  = ia_min((*x), (*y));
	      top--;
	      break;
            default:
              fatalerror("ia_eval: bad op code");
	  }
        iop++;
      }
  }

void ia_eval_diff (IntervalDiff reg[], IntervalDiff stack[], pcode_instr_t fcode[])
  {
    int top = -1; 
    IntervalDiff *x, *y, *z;
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
	      z->f = ia_int_const(fcode[iop].arg);
	      z->df = (Interval) {Zero, Zero};
	      top++;
	      break;
	    case pcode_op_store:
	      reg[fcode[iop].arg] = *y; 
	      top--;
	      break;
	    case pcode_op_return:
	      fatalerror("ia_eval: unexpected return!");
	      break;
	    case pcode_op_add:
	      x->df = ia_add(x->df,y->df);
	      x->f  = ia_add(x->f, y->f);
	      top--;
	      break;
	    case pcode_op_sub:
	      x->df = ia_sub(x->df,y->df);
	      x->f  = ia_sub(x->f, y->f);
	      top--;
	      break;
	    case pcode_op_neg:
	      y->df = ia_neg(y->df);
	      y->f  = ia_neg(y->f);
	      break;
	    case pcode_op_inv:
	      y->df = ia_neg(ia_div(y->df, ia_sqr(x->f)));
              y->f  = ia_inv(y->f);
	      break;
	    case pcode_op_mul:
	      x->df = ia_add(ia_mul(x->f, y->df), ia_mul(x->df,y->f));
	      x->f  = ia_mul(x->f, y->f);
	      top--;
	      break;
	    case pcode_op_div:
	      x->df = ia_div(ia_sub(ia_mul(x->df, y->f), ia_mul(x->f, y->df)),
			     ia_sqr(y->f));
	      x->f  = ia_div(x->f, y->f);
	      top--; break;
	    case pcode_op_sqr:
	      y->df = ia_mul(ia_scale(y->f,2,1), y->df);
	      y->f  = ia_sqr(y->f);
	      break;
	    case pcode_op_sqrt:
	      { Interval r = ia_sqrt(y->f);
	        y->df = ia_div(y->df, ia_scale(r,1,2));
                y->f  = r;
              }
	      break;
	    case pcode_op_abs:
	      if ((y->f.lo <= Zero) && (y->f.hi >= Zero))
                { y->df = ia_join(y->df, ia_neg(y->df)); }
              else
                { y->df = (y->f.lo > Zero ? y->df : ia_neg(y->df)); }
              y->f  = ia_abs(y->f);
	      break;
	    case pcode_op_max:
	      if ((x->f.hi >= y->f.lo) && (y->f.hi >= x->f.lo))
                { x->df = ia_join(x->df, y->df); }
              else
                { x->df = (x->f.lo > y->f.lo ? x->df : y->df); }
              x->f  = ia_max((x->f), (y->f));
	      top--;
	      break;
	    case pcode_op_min:
	      if ((x->f.hi >= y->f.lo) && (y->f.hi >= x->f.lo))
                { x->df = ia_join(x->df, y->df); }
              else
                { x->df = (x->f.lo > y->f.lo ? y->df : x->df); }
              x->f  = ia_min((x->f), (y->f));
	      top--;
	      break;
            default:
              fatalerror("ia_eval_diff: bad op code");
	  }
        iop++;
      }
  }
