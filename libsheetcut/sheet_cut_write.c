/* See {sheet_cut_write.h} */
/* Last edited on 2024-12-05 10:40:13 by stolfi */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <math.h>

#include <r2.h>
#include <bool.h>
#include <affirm.h>
#include <jsmath.h>

#include <sheet_cut.h>
#include <sheet_cut_write.h>

/* IMPLEMENTATIONS */  

void sheet_cut_write_all_plates
  ( FILE *wr, 
    int32_t sheet_ix, 
    int32_t *plate_ctP, 
    sheet_cut_node_t *pc,
    r2_t org_pc
  )
  { 
    auto void visit(sheet_cut_node_t *pr, r2_t org_pr, bool_t post);
      /* Visiting procedure for {sheet_cut_enum}. 
      Will be called for each plate or block {pr} reachable from {pc}.
      
      If {pr} is a plate, prints its info with 
      {sheet_cut_write_plate} and increment {*plate_ctP}.
      Otherwise does nothing. */
      
    sheet_cut_enum(pc, org_pc, visit);
    return;
    
    /* Inpternal implementations: */
      
    void visit(sheet_cut_node_t *pr, r2_t org_pr, bool_t post)
      { 
        assert(sheet_cut_same_mat_thk(pr, pc));
        if (pr->nsub == 0)
          { /* Plate; print it: */
            r2_t low_pr; r2_add(&org_pr, &(pr->pos), &low_pr); /* Plot coords of {pr}'s low corner. */
            sheet_cut_write_plate
              ( wr, sheet_ix, pr->mat, pr->thk,
                (*plate_ctP), low_pr, pr->size, pr->tag
              );
            (*plate_ctP)++;
          }
      }
  
  }

void sheet_cut_write_plate
  ( FILE *wr, 
    int32_t sheet_ix, 
    char *mat, 
    double thk, 
    int32_t plate_ix, 
    r2_t pos, 
    r2_t size, 
    char* tag
  )
  { 
    fprintf
      ( wr, 
        "%03d %-10s %6.1f  %05d  %6.1f %6.1f  %6.1f %6.1f %s\n", 
        sheet_ix, mat, thk, plate_ix, pos.c[0], pos.c[1], size.c[0], size.c[1], tag
      );
  }
