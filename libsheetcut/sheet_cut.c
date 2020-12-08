/* See {sheet_cut.h} */
/* Last edited on 2019-12-08 13:28:40 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <r2.h>
#include <vec.h>
#include <bool.h>
#include <jsmath.h>
#include <jsstring.h>
#include <affirm.h>

#include <sheet_cut.h>
    
vec_typeimpl(sheet_cut_node_vec_t, sheet_cut_node_vec, sheet_cut_node_ref_t);

bool_t sheet_cut_same_mat_thk(sheet_cut_node_t *pca, sheet_cut_node_t *pcb)
  { 
    return ((strcmp(pca->mat, pcb->mat) == 0) && (pca->thk == pcb->thk));
  }

double sheet_cut_area_m2(r2_t size)
  { 
    double area = (size.c[0] * size.c[1])/1.0e6;
    return area;
  }
  
void sheet_cut_print_node(FILE *wr, int32_t ind, sheet_cut_node_t *pc, bool_t sub)
  { 
    assert(pc->mat != NULL);
    assert(pc->nsub >= 0);
    fprintf(wr, "%*s%s", ind, "", (pc->nsub == 0 ? "plate" : "block"));
    fprintf(wr, " tag = %s", (pc->tag != NULL ? pc->tag : "NULL"));
    fprintf(wr, " mat = %s thk = %.1f\n", pc->mat, pc->thk);
    
    fprintf(wr, "%*ssize = (%.1f,%.1f)", ind, "", pc->size.c[0], pc->size.c[1]);
    if (pc->nsub > 0)
      { fprintf(wr, " %d children%s\n", pc->nsub, (pc->sub_flip ? " (flipped)" : ""));
        if (sub)
          { for (int32_t ich = 0; ich < pc->nsub; ich++)
              { sheet_cut_node_t *chi = pc->sub.e[ich];
                sheet_cut_print_node(wr, ind+2, chi, sub);
              }
          }
      }
    fflush(wr);
  }

sheet_cut_node_t* sheet_cut_new_plate(char *mat, double thk, r2_t size, char *tag)
  { demand(size.c[0] > 0, "invalid plate X size"); 
    demand(size.c[1] > 0, "invalid plate Y size"); 
    demand(thk > 0.0, "invalid plate thickness");    
  
    sheet_cut_node_t *pc = notnull(malloc(sizeof(sheet_cut_node_t)), "no mem");
    
    /* General attributes: */
    pc->mat = mat;
    pc->thk = thk;
    pc->pos = (r2_t){{ 0.0, 0.0 }};
    pc->size = size;
    pc->tag = tag;
    
    /* No children: */
    pc->sub = sheet_cut_node_vec_new(0);
    pc->nsub = 0;
    pc->sub_flip = FALSE;
    
    return pc;
  }
  
sheet_cut_node_t* sheet_cut_new_block(sheet_cut_node_t *pc)
  { sheet_cut_node_t *blk = notnull(malloc(sizeof(sheet_cut_node_t)), "no mem");

    /* General attributes: */
    blk->mat = pc->mat;
    blk->thk = pc->thk;
    blk->pos = (r2_t){{ 0.0, 0.0 }};
    r2_add(&(pc->pos), &(pc->size), &(blk->size));
    blk->tag = txtcat(pc->tag,".P");
    
    /* Only child: */
    blk->sub_flip = FALSE;
    blk->sub = sheet_cut_node_vec_new(10);
    blk->sub.e[0] = pc;
    blk->nsub = 1;
    
    return blk;
  }

void sheet_cut_flip_node(sheet_cut_node_t *pc)
  { 
    { double tmp = pc->pos.c[0];  pc->pos.c[0] = pc->pos.c[1];  pc->pos.c[1] = tmp; }
    { double tmp = pc->size.c[0]; pc->size.c[0] = pc->size.c[1]; pc->size.c[1] = tmp; }
    pc->sub_flip = (! pc->sub_flip);
  }
  
sheet_cut_node_t* sheet_cut_add_child(sheet_cut_node_t *blk,  sheet_cut_node_t *pc)
  { 
    if (pc != NULL)
      { 
        demand(pc->mat != NULL, "undefined child material");
        demand(! isnan(pc->thk), "undefined child thickness");
        demand((pc->pos.c[0] >= 0.0) && (pc->pos.c[1] >= 0.0), "invalid child position");
        demand((pc->size.c[0] > 0.0) && (pc->size.c[1] > 0.0), "invalid child size");

        if (blk == NULL)
          { /* Create a new block record, with {pc} as the only child: */
            blk = sheet_cut_new_block(pc); 
          }
        else
          { /* Add {pc} to the block's children. */
            demand(strcmp(blk->mat, pc->mat) == 0, "incompatible child and block materials");
            demand(blk->thk == pc->thk, "incompatible child and block thicknesses");
            demand((blk->pos.c[0] >= 0.0) && (blk->pos.c[1] >= 0.0), "invalid block position");
            demand((blk->size.c[0] > 0.0) && (blk->size.c[1] > 0.0), "invalid block size");

            if (blk->sub_flip) 
              { /* Flip {pc} so that it will be un-flipped when used: */
                sheet_cut_flip_node(pc); 
              }

            /* Append {pc} to the list of children: */
            sheet_cut_node_vec_expand(&(blk->sub), blk->nsub);
            blk->sub.e[blk->nsub] = pc;
            blk->nsub++;

            /* Expand the block's bounding box: */
            r2_t pc_hi; r2_add(&(pc->pos), &(pc->size), &pc_hi); /* Upper corner of child rel to block. */
            blk->size.c[0] = fmax(blk->size.c[0], pc_hi.c[0]);
            blk->size.c[1] = fmax(blk->size.c[1], pc_hi.c[1]);
          }
      }
    assert((pc == NULL) || (blk != NULL));
    return blk;
  }

void sheet_cut_enum(sheet_cut_node_t *pc, r2_t org_pc, sheet_cut_visit_proc_t *visit)
  {
    if (pc->nsub == 0)
      { /* A plate; just visit it once: */
        visit(pc, org_pc, FALSE);
      }
    else
      { /* A block; visit it, enum descendants, revisit it. */
        
        if (pc->sub_flip)
          { /* Hard-flip all children and reset {pc->sub_flip}: */
            for (int32_t i = 0; i < pc->nsub; i++)
              { sheet_cut_flip_node(pc->sub.e[i]); }
            pc->sub_flip = FALSE;
          }
        
        /* Visit the block (pre-order): */
        visit(pc, org_pc, FALSE);
        
        /* Enumerate the children  nodes: */
        r2_t pos_pc; r2_add(&org_pc, &(pc->pos), &pos_pc); /* Coords of corner of {pc}. */
        for (int32_t i = 0; i < pc->nsub; i++) 
          { sheet_cut_enum(pc->sub.e[i], pos_pc, visit); }
        
        /* Visit the block again (post-order): */
        visit(pc, org_pc, TRUE);
      }
  }

