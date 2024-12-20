#ifndef sheet_cut_H
#define sheet_cut_H
 
/* Basic types for packing rectangular plates on sheet stock. */
/* Last edited on 2024-12-05 10:40:06 by stolfi */

#include <stdint.h>
#include <stdio.h>

#include <r2.h>
#include <vec.h>
#include <bool.h>
#include <jsmath.h>

typedef struct sheet_cut_node_t * sheet_cut_node_ref_t;
   
vec_typedef(sheet_cut_node_vec_t, sheet_cut_node_vec, sheet_cut_node_ref_t);
  /* A {sheet_cut_node_vec_t} is an extensible, self-sized vector of
    pointers to {sheet_cut_node_t} records. */

typedef struct sheet_cut_node_t
  { char *mat;                 /* Sheet material. */
    double thk;                /* Sheet thickness. */
    r2_t pos;                  /* Lower corner of node rel to lower corner of parent block.  */
    r2_t size;                 /* Dimensions (width and height) of the plate or block. */
    /* For plates, possibly also for blocks: */
    char *tag;                 /* A label that identifies the plate globally. */
    /* For blocks: */
    sheet_cut_node_vec_t sub; /* The children nodes. */
    int32_t nsub;              /* Number of children. */
    bool_t sub_flip;           /* If {TRUE}, children positions and sizes are to be XY-flipped. */
  } sheet_cut_node_t;
  /* 
    A record {pc} of this type describes a /plate/, a rectangular piece of sheet
    material to be cut from some sheet stock; or a /block/ of such plates,
    that are to be cut out of the same sheet stock with fixed positions
    relative to each other, and should be treated as a single rectangular object 
    for placement purposes.
    
    The node is a plate if {n = pc.sub.ne} is zero.  Otherwise, 
    it is a block, and its immediate /children/ nodes are {s[0..n-1] =
    pc.sub.e[0..n-1]}.  The children can be themselves blocks, and so
    on recursively.
    
    The children links connect the nodes as a tree, rather than a DAG or
    general graph. That is, each node {pc} has at most one /parent
    block/ that includes {pc} among its children; and any path through
    child links eventually reaches a leaf node (a plate).
    
    In any case, the stock sheet out of which {pc} should
    be cut out has type {pc.mat} (e. g., "plywood", "MDF", "steel") 
    and thickness {pc.thk}. All the nodes that descend from {pc}  
    must have the same {.mat} and {.thk} fields.
    
    The field {pc.pos} is the low (lower left) corner this plate or
    block relative to the low corner of the parent block. Its
    coordinates cannot be negative. If {pc} has no parent, {pc.pos} is
    normally {(0,0)}. The position {s[i].pos} of each child node, if any,
    is relative to {pc}'s low corner.
    
    If {pc} is a plate, {pc.size} is simply its size. If {pc} is a block, 
    {pc.size} is the size of a rectangle with lower left corner
    {(0,0)} that encloses the box {s[i].size} of every child node, each shifted 
    by the corresponding position {s[i].pos}.  That is, the rectangle
    from {(0,0)} to {pc.size} contains the relative low corner {s[i].pos}
    and the relative high corner {s[i].pos + s[i].size} of every child node
    {s[i]}.
    
    However, if {.sub_flip} is {TRUE}, the relative position {s[i].pos} and the
    size {s[i].size} of each child must be assumed to be XY-flipped with
    respect to what is actually stored in those fields. Namely, the child
    {s[i]} should be placed at {pc.pos + FLIP(s[i].pos)} and its size
    should be assumed to be {FLIP(s[i].size); where {FLIP((X,Y)) =
    (Y,X)}. If {s[i]} is itself a block, its contents should be
    XY-flipped too. */  
    
sheet_cut_node_t* sheet_cut_new_plate(char *mat, double thk, r2_t size, char *tag);
  /* Creates a new plate node {pc}, with the specified properties.
    The position {pc.pos} will be set to {(0,0)}. */
    
sheet_cut_node_t* sheet_cut_new_block(sheet_cut_node_t *pc);
  /* Creates a block node {blk}, initially with {pc} as the only child.
    The block will have {blk.pos = (0,0)}, {blk.sub_flip = FALSE}, and
    the appropriate {blk.size}. The tag of the new block will be that of
    {pc} with ".P" appended.*/

sheet_cut_node_t* sheet_cut_add_child(sheet_cut_node_t *blk,  sheet_cut_node_t *pc);
  /* Adds the node {pc} as a child of {blk}, which must be
    a block node.  Returns the expanded block {blk}.
    
    However, if the node {pc} is {NULL}, the operation is a no-op: the 
    block {blk} is returned unchanged (even if it is {NULL}).  
    
    Otherwise, {pc} must have non-negative relative position coordinates
    {pc.pos.c[0..1]}; and positive dimensions {pc.size.c[0..1]}.
     
    If {pc} is not {NULL} but {blk} is NULL, the procedure will create
    and return a new block record with {pc} as its only
    child, with {sheet_cut_new_block(pc)}. Note that the new block will
    have {blk.pos = (0,0)} and {blk.sub_flip = FALSE}.
    
    If neither {pc} nor {blk} are null, {pc} is appended to the list of
    children {s[0..n-1]} of {blk}. In this case, {pc-mat} and {pc-thk}
    must match those of all the other children of {blk}, namely
    {blk->mat} and {blk->thk}. The new child will be flipped if
    {blk.sub_flip} is true.

    In any case, if the resulting block {blk} is not {NULL},
    its bounding box size {blk.size} will be udated aas necessary
    to include all the rectanges {s[i].size} shifted by {s[i].pos};
    including that of the new child {pc} if not {NULL}. */

void sheet_cut_flip_node(sheet_cut_node_t *pc);
  /* Swaps all X and Y coordinates of the node {pc} and all its contents. Namely,
    swaps X and Y in {pc.pos} and {pc.size}, and inverts {pc.sub_flip}. */
 
bool_t sheet_cut_same_mat_thk(sheet_cut_node_t *pca, sheet_cut_node_t *pcb);
  /* Returns {TRUE} iff the nodes {pca} and {pcb} have the same
    material and thickness. */
    
double sheet_cut_area_m2(r2_t size);
  /* Returns the area of a rectangle, in square meters, given its usable 
    dimensions {size} in millimeters. */
    
void sheet_cut_print_node(FILE *wr, int32_t ind, sheet_cut_node_t *pc, bool_t sub);
  /* Writes the data of the node {pc} to {wr}, indented by {ind}, in readable format.
    If {sub} is true and {pc} is a block, also prints all descendant, recursively. */
    
typedef void sheet_cut_visit_proc_t (sheet_cut_node_t *pr, r2_t org_pr, bool_t post); 
  /* Type of the node visit procedure expected by {sheet_cut_enum}. */

void sheet_cut_enum(sheet_cut_node_t *pc, r2_t org_pc, sheet_cut_visit_proc_t *visit);
  /* Enumerates all nodes (blocks and plates) in the 
    tree whose root is the node {pc}, in double order 
    (preorder and postorder).
    
    For each such node {pr}, calls {visit(pr, org_pr, FALSE)} on first
    encounter, before visiting all its descendants. If {pr} is a block,
    calls again {visit(pr, org_pr, TRUE)} after visiting all of {pr}'s
    descendants.
    
    In these calls, {org_pr} is the assumed coords of the low corner of
    {pr}'s parent, assuming that the low corner of {pc}'s parent is
    {org_pc}. So, the low corner of the plate or block {pr} is assumed to be at
    coordinates {org_pr + pr->pos}. In particular, the low corner of
    {pc} is assumed to be at {org_pc + pc->pos}.
    
    Assumes that {pc} itself does not require any flipping. Before
    calling {visit} on the children of any block {blk} reachable
    from {pc}, including possibly {pc}, checks if {blk.sub_flip} is
    {TRUE}. If so, applies {sheet_cut_flip_node} to all those children,
    and sets {blk.sub_flip} to {FALSE}. Thus, after the first call to
    {sheet_cut_enum}, all nodes reachable from {pc} will have
    {.sub_flip} set to {FALSE}, and their coordinates and sizes will no
    longer require any flipping. */

#endif
