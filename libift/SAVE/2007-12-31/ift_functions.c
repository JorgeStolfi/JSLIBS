/* ift_functions.c - implementation of ift_functions.h */
/* Last edited on 2008-05-24 12:05:32 by stolfi */

#include "ift_functions.h"
#include "ift.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
  
/* FUNCTION SELECTION BY NAME */

PathCostFn *get_cost_function(char *name)
  {
    if (strcmp(name, "maxval") == 0)
      { return &pf_maxval; }
    else if (strcmp(name, "sumval") == 0)
      { return &pf_sumval; }
    else if (strcmp(name, "maxdiff") == 0)
      { return &pf_maxdiff; }
    else if (strcmp(name, "sumdiff") == 0)
      { return &pf_maxdiff; }
    else if (strcmp(name, "monoroot") == 0)
      { return &pf_monoroot; }
    else if (strcmp(name, "fuzzconn") == 0)
      { return &pf_fuzzconn; }
    else 
      { IFT_ERROR("bad function name"); }
    return NULL;
  }

/* PATH COST FUNCTIONS */

PathCost pf_maxdiff(PixelNode *s, RelArc *a, PixelNode *t, ImageGraph *G, int stage)
  {
    PathCost tC_new = 0;
    if ((s == NULL) != (a == NULL)) { IFT_ERROR("s/a inconstency"); }
    switch(stage)
      { 
        case 0: 
          break;
        case 1: 
          if (t == NULL) { IFT_ERROR("t is null"); }
          if ((s != NULL) && (a != NULL))
            { double diff = abs_pixel_diff(&(s->y), &(t->y), G->channels);
              double ndiff = diff/((double)G->channels);  /* Scaled to `[0_1]' */
              PathCost aC = arc_cost_quantize(ndiff/a->len);
              tC_new = (aC > s->C ? aC : s->C);
            }
          break;
        case 2: 
          break;
        default:
          IFT_ERROR("bad stage");
      }
    return tC_new;
  }

PathCost pf_sumdiff(PixelNode *s, RelArc *a, PixelNode *t, ImageGraph *G, int stage)
  {
    PathCost tC_new = 0;
    if ((s == NULL) != (a == NULL)) { IFT_ERROR("s/a inconstency"); }
    switch(stage)
      { 
        case 0: 
          break;
        case 1: 
          if (t == NULL) { IFT_ERROR("t is null"); }
          if ((s != NULL) && (a != NULL))
            { double diff = abs_pixel_diff(&(s->y), &(t->y), G->channels);
              double ndiff = diff/((double)G->channels);  /* Scaled to `[0_1]' */
              PathCost aC = arc_cost_quantize(ndiff/a->len);
              tC_new = s->C + aC;
            }
          break;
        case 2: 
          break;
        default:
          IFT_ERROR("bad stage");
      }
    return tC_new;
  }

PathCost pf_maxval(PixelNode *s, RelArc *a, PixelNode *t, ImageGraph *G, int stage)
  {
    PathCost tC_new = 0;
    if ((s == NULL) != (a == NULL)) { IFT_ERROR("s/a inconstency"); }
    switch(stage)
      { 
        case 0: 
          break;
        case 1: 
          if (t == NULL) { IFT_ERROR("t is null"); }
          { double tval = abs_pixel_value(&(t->y), G->channels);
            double ntval = tval/((double)G->channels);  /* Scaled to `[0_1]' */
            tC_new = arc_cost_quantize(ntval);
            if ((s != NULL) && (tC_new < s->C)) { tC_new = s->C; }
          } 
          break;
        case 2: 
          break;
        default:
          IFT_ERROR("bad stage");
      }
    return tC_new;
  }

PathCost pf_sumval(PixelNode *s, RelArc *a, PixelNode *t, ImageGraph *G, int stage)
  {
    PathCost tC_new = 0;
    if ((s == NULL) != (a == NULL)) { IFT_ERROR("s/a inconstency"); }
    switch(stage)
      { 
        case 0: 
          break;
        case 1: 
          if (t == NULL) { IFT_ERROR("t is null"); }
          { double tval = abs_pixel_value(&(t->y), G->channels);
            double ntval = tval/((double)G->channels);  /* Scaled to `[0_1]' */
            tC_new = arc_cost_quantize(ntval);
            if (s != NULL) { tC_new += s->C; }
          } 
          break;
        case 2: 
          break;
        default:
          IFT_ERROR("bad stage");
      }
    return tC_new;
  }

PathCost pf_monoroot(PixelNode *s, RelArc *a, PixelNode *t, ImageGraph *G, int stage)
  {
    PathCost tC_new = 0;
    if ((s == NULL) != (a == NULL)) { IFT_ERROR("s/a inconstency"); }
    switch(stage)
      { 
        case 0: 
          break;
        case 1: 
          if (t == NULL) { IFT_ERROR("t is null"); }
          { double tval = abs_pixel_value(&(t->y), G->channels);
            double ntval = tval/((double)G->channels);  /* Scaled to `[0_1]' */
            PathCost tint = arc_cost_quantize(ntval);
            if (s == NULL) 
              { tC_new = tint; }
            else 
              { PathCost sint = s->C;
                if (tint < sint)
                  { tC_new = INFINITE_PATH_COST; }
                else
                  { tC_new = tint; }
              }
          } 
          break;
        case 2: 
          break;
        default:
          IFT_ERROR("bad stage");
      }
    return tC_new;
  }

PathCost pf_fuzzconn(PixelNode *s, RelArc *a, PixelNode *t, ImageGraph *G, int stage)
  { PathCost tC_new = 0;
    IFT_ERROR("not implemented yet");
    return tC_new;
  }

/* AUXILIARY FUNCTIONS */  

PathCost arc_cost_quantize(double x)
  {
    PathCost y;
    if (x == INFINITY) { return INFINITE_ARC_COST; }
    if (x > 1.0) { IFT_ERROR("raw arc cost > 1.0"); }
    if (x < 0.0) { IFT_ERROR("raw arc cost < 0.0"); }
    y = floor(x * ((double)MAX_FINITE_ARC_COST) + 0.9999);
    return y;
  }

double abs_pixel_diff(PixelValue *p, PixelValue *q, int channels)
  {
    int k;
    double s = 0;
    for(k = 0; k < channels; k++) { s += fabs(p->c[k] - q->c[k]); }
    return s;
  }

double abs_pixel_value(PixelValue *p, int channels)
  {
    int k;
    double s = 0;
    for(k = 0; k < channels; k++) { s += fabs(p->c[k]); }
    return s;
  }

