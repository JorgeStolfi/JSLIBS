#ifndef sheet_cut_write_H
#define sheet_cut_write_H
 
/* Writing out layouts of rectangular plates on sheet stock. */
/* Last edited on 2024-12-05 10:40:15 by stolfi */

#include <stdint.h>
#include <stdio.h>

#include <r2.h>
#include <bool.h>
#include <jsmath.h>

#include <sheet_cut.h>

void sheet_cut_write_all_plates
  ( FILE *wr, 
    int32_t sheet_ix, 
    int32_t *plate_ctP, 
    sheet_cut_node_t *pc,
    r2_t org_pc
  );
  /* Writes to {wr} a one-line descrition of every plate 
    (leaf node) that descends from {pc}, using {sheet_cut_write_plate}.
    
    Assumes that the low corner of {pc} is at {org_pc + pc->pos}.
    
    Also increments {plate_ct} for each plate found. The caller
    must initialize that counter. The value of that counter will be used
    as the plate index for {sheet_cut_write_plate}.
    
    As a side effect, all nodes that descend from {pc} will be actually
    flipped when appropriate, and all fields {.sub_flip} will then be
    set to {FALSE}, before writing any plates out. Assumes that {pc}
    itself does not need flipping.
    
    This procedure is usually called with {pc} being the root block of a
    stock sheet, and {org_pc} being the coordinates of the low corner of
    the usable area of the sheet. */

void sheet_cut_write_plate
  ( FILE *wr, 
    int32_t sheet_ix, 
    char *mat, 
    double thk, 
    int32_t plate_ix, 
    r2_t pos, 
    r2_t size, 
    char* tag
  );
  /* Writes to {wr} the information about a placed plate, in the format 
    described by {sheet_cut_write_plate_INFO} below. */
    
#define sheet_cut_write_plate_INFO \
  "The description of a placed plate contains 9 fields\n" \
  "\n" \
  "    \"{SHIX} {MAT} {THK}  {PLIX}  {PX} {PY}  {DX} {DY} {TAG}\"\n" \
  "\n" \
  "  where {SHIX} is a sequential sheet index, {MAT,THK} are the" \
  " material type (a string) and thickness (a number), {PLIX} is" \
  " the sequential of the plate in the sheet, {PX,PY} are the" \
  " coordinates of the low corner of the plate (relatve to" \
  " the usable area of the sheet), {DX,DY} are the plate's dimensions, and" \
  " {TAG} is the plate's identifying tag." 

#endif
