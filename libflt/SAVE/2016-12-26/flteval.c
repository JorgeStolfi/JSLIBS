/* See "flteval.h" */
/* Last edited on 2016-04-01 02:43:53 by stolfilocal */

#include "flteval.h"
#include <flt.h>
#include <pcode.h>
#include <affirm.h>
#include <math.h>

void flt_eval (Float reg[], Float stack[], pcode_instr_t fcode[])
  {
    int top = -1; 
    Float *x, *y, *z;
    int iop = 0;
    pcode_op_t op;

    ROUND_NEAR;
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
	      *z = (Float) fcode[iop].arg;
	      top++;
	      break;
	    case pcode_op_store:
	      reg[fcode[iop].arg] = *y; 
	      top--;
	      break;
	    case pcode_op_return:
	      fatalerror("flt_eval: unexpected return!");
	      break;
	    case pcode_op_add:
	      *x  = (Float)((*x) + (*y));
	      top--;
	      break;
	    case pcode_op_sub:
	      *x  = (Float)((*x) - (*y));
	      top--;
	      break;
	    case pcode_op_neg:
	      *y  = (Float)(- (*y));
	      break;
	    case pcode_op_inv:
	      *y  = (Float)(1.0 / (*y));
	      break;
	    case pcode_op_mul:
	      *x  = (Float)((*x) * (*y));
	      top--;
	      break;
	    case pcode_op_div:
	      *x  = (Float)((*x) / (*y));
	      top--; 
              break;
	    case pcode_op_sqr:
	      *y  = (Float)((*y) * (*y));
	      break;
	    case pcode_op_sqrt:
	      *y  = (Float)sqrt((double)(*y));
	      break;
	    case pcode_op_abs:
	      *y  = (Float)(FABS(*y));
	      break;
	    case pcode_op_max:
	      *x  = (Float)(FMAX((*x), (*y)));
	      top--;
	      break;
	    case pcode_op_min:
	      *x  = (Float)(FMIN((*x), (*y)));
	      top--;
	      break;
            default:
              fatalerror("flt_eval: bad op code");
	  }
        iop++;
      }
  }

void flt_eval_diff (FloatDiff reg[], FloatDiff stack[], pcode_instr_t fcode[])
  {
    int top = -1; 
    FloatDiff *x, *y, *z;
    int iop = 0;
    pcode_op_t op;

    ROUND_NEAR;
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
	      z->f = (Float) fcode[iop].arg;
	      z->df = Zero;
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
	      x->df = (Float)((x->df)+(y->df));
	      x->f  = (Float)((x->f)+(y->f));
	      top--;
	      break;
	    case pcode_op_sub:
	      x->df = (Float)((x->df)-(y->df));
	      x->f  = (Float)((x->f)-(y->f));
	      top--;
	      break;
	    case pcode_op_neg:
	      y->df = (Float)(-(y->df));
	      y->f  = (Float)(-(y->f));
	      break;
	    case pcode_op_inv:
	      y->df = (Float)(- (y->df) / ((x->f) * (x->f)));
              y->f  = (Float)(1.0 / (y->f));
	      break;
	    case pcode_op_mul:
	      x->df = (Float)((x->f) * (y->df) + (x->df) * (y->f));
	      x->f  = (Float)((x->f) * (y->f));
	      top--;
	      break;
	    case pcode_op_div:
	      x->df = (Float)((x->df) / (y->f) - (x->f)*(y->df)/((y->f)*(y->f)));
	      x->f  = (Float)((x->f) / (y->f));
	      top--; break;
	    case pcode_op_sqr:
	      y->df = (Float)(2.0 * (y->f) * (y->df));
	      y->f  = (Float)((y->f) * (y->f));
	      break;
	    case pcode_op_sqrt:
	      { double r = sqrt((double)(y->f));
	        y->df = (Float)(0.5 * (y->df) / r);
                y->f  = (Float)(r);
              }
	      break;
	    case pcode_op_abs:
	      y->df = (Float)(((y->f) > Zero ? (y->df) : -(y->df)));
              y->f  = (Float)(FABS(y->f));
	      break;
	    case pcode_op_max:
	      x->df = (Float)(((x->f) >= (y->f) ? (x->df) : (y->df)));
              x->f  = (Float)(FMAX((x->f), (y->f)));
	      top--;
	      break;
	    case pcode_op_min:
	      x->df = (Float)(((x->f) >= (y->f) ? (y->df) : (x->df)));
              x->f  = (Float)(FMIN((x->f), (y->f)));
	      top--;
	      break;
            default:
              fatalerror("flt_eval_diff: bad op code");
	  }
        iop++;
      }
  }

void flt_eval_diff4 (FloatDiff4 reg[], FloatDiff4 stack[], pcode_instr_t fcode[])
  {
    int top = -1; 
    FloatDiff4 *x, *y, *z;
    int iop = 0;
    pcode_op_t op;
    int i;
    
#   define FORI for(i=0; i<4; i++)

    ROUND_NEAR;
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
	      z->f = (Float)(fcode[iop].arg);
	      FORI { z->df[i] = Zero; }
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
	      FORI { x->df[i] = (Float)((x->df[i])+(y->df[i])); }
	      x->f  = (Float)((x->f)+(y->f));
	      top--;
	      break;
	    case pcode_op_sub:
	      FORI { x->df[i] = (Float)((x->df[i])-(y->df[i])); }
	      x->f  = (Float)((x->f)-(y->f));
	      top--;
	      break;
	    case pcode_op_neg:
	      FORI { y->df[i] = (Float)(-(y->df[i])); }
	      y->f  = (Float)(-(y->f));
	      break;
	    case pcode_op_inv:
	      FORI { y->df[i] = (Float)(- (y->df[i]) / ((x->f) * (x->f))); }
              y->f  = (Float)(1.0 / (y->f));
	      break;
	    case pcode_op_mul:
	      FORI { x->df[i] = (Float)((x->f) * (y->df[i]) + (x->df[i]) * (y->f)); }
	      x->f  = (Float)((x->f) * (y->f));
	      top--;
	      break;
	    case pcode_op_div:
	      FORI { x->df[i] = (Float)((x->df[i]) / (y->f) - (x->f)*(y->df[i])/((y->f)*(y->f))); }
	      x->f  = (Float)((x->f) / (y->f));
	      top--; break;
	    case pcode_op_sqr:
	      FORI { y->df[i] = (Float)(2.0 * (y->f) * (y->df[i])); }
	      y->f  = (Float)((y->f) * (y->f));
	      break;
	    case pcode_op_sqrt:
	      { double r = sqrt((double)(y->f));
	        FORI { y->df[i] = (Float)(0.5 * (y->df[i]) / r); }
                y->f  = (Float)(r);
              }
	      break;
	    case pcode_op_abs:
	      FORI { y->df[i] = (Float)(((y->f) > Zero ? (y->df[i]) : -(y->df[i]))); }
              y->f  = (Float)(FABS(y->f));
	      break;
	    case pcode_op_max:
	      FORI { x->df[i] = (Float)(((x->f) >= (y->f) ? (x->df[i]) : (y->df[i]))); }
              x->f  = (Float)(FMAX((x->f), (y->f)));
	      top--;
	      break;
	    case pcode_op_min:
	      FORI { x->df[i] = (Float)(((x->f) >= (y->f) ? (y->df[i]) : (x->df[i]))); }
              x->f  = (Float)(FMIN((x->f), (y->f)));
	      top--;
	      break;
            default:
              fatalerror("flt_eval_diff4: bad op code");
	  }
        iop++;
      }
      
#   undef FORI
  }
