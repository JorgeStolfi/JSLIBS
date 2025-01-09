/* imgsys - linear system solving for image-related problems. */

#ifndef pst_imgsys_H
#define pst_imgsys_H

/* Created on 2005-12-04 by Jorge Stolfi, unicamp, <stolfi@ic.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2025-01-09 06:08:55 by stolfi */
/* See the copyright and authorship notice at the end of this file. */

#include <stdint.h>
#include <stdio.h>

#include <bool.h>
#include <float_image.h>

#define pst_imgsys_MAX_COEFFS 5
  /* Max coefficients in an equation. */

#define MAX_COEFFS pst_imgsys_MAX_COEFFS
  /* Shorter local name. */

typedef struct pst_imgsys_equation_t
  { uint32_t nt;             /* {nt} is the number of terms in the equation. */
    uint32_t uid[MAX_COEFFS];   /* {uid[k]} is the index of some unknown in the equation. */
    double cf[MAX_COEFFS];     /* Coefficient of that unknown. */
    double rhs;              /* Right-hand side of equation. */
    double wtot;             /* Weight of the equation and of its main variable. */
  } pst_imgsys_equation_t;
  /* An {pst_imgsys_equation_t} record {eq} represents a linear equation 
    with at most {MAX_COEFFS} nonzero terms, namely
      {SUM {eq.cf[j]*h[eq.uid[j]] : j = 0..eq.nt-1} == rhs},
    where {h[j]} is the unknown with index {j}.
    The variable with index {eq.uid[0]} is the main variable of the equation. */

typedef struct pst_imgsys_t
  { uint32_t N;                 /* Number of unknowns and equations. */
    pst_imgsys_equation_t* eq;  /* The equations are {eq[0..N-1]}. */
    /* Debugging info: */
    int32_t NX, NY;             /* Number of cols and rows of the assumed grid. */
    int32_t *col;               /* Maps index of an unknown/equation to a {x} index, or -1. */
    int32_t *row;               /* Maps index of an unknown/equation to a {y} index, pr -1. */
    int32_t *uid;                /* Maps {x + y*NX_Z} to index of unknown, or -1. */
  } pst_imgsys_t;
  /* A {pst_imgsys_t} record {S} represents {N} linear equations
    {eq[0..N-1]} on {N} unknowns {h[0..N-1}. In a proper system, the
    main variable of equation {k} has index {k} --- that is,
    {eq[k].uid[0] == k} --- for all {k} in {0..N-1}.
    
    The tables {S.col,S.row,S.uid} are used only when printing the
    system. Specifically, {S.col[0..N-1]} and {S.row[0..N-1]} define a
    partial one-to-one mapping from unknown/equation index {k} in
    {0..N-1} to a point {(x,y)} of a subset of a 2D integer grid with
    {x} in {0..NX-1} and {y} in {0..NY-1}. If the unknown {k} does not
    correspond to a grd point, then {S.col[k]=S.row[k]=-1}.
    
    The table {S.uid[0..NX*NY-1]} is the inverse mapping, from the col
    and row indices {(x,y)} to the corresponding unknown index {k =
    S.uid[x+y*NX]}; or {-1} if there is no unknown corresponding to that
    col/row pair. */
    
pst_imgsys_t *pst_imgsys_new_grid(int32_t NX, int32_t NY);
  /* Creates a new linear system {S} with {N=NX*NY} equations on {N} unknowns.
    Each equation {eqk=S.eq[k]} intially has zero righ-hand side {eqk.rhs} and
    total weight {eqk.wtot}, and one term ({eqk.nt=1}) with unknown 
    {eqk.uid[0]=k} and coefficient {eqk.cf[0]=0}.   Namely, the equation is
    indeterminate, {0*Z[k]=0}.
    
    The tables {S.col,S.row,S.uid} are filled assuming that the unknowns
    (and therefore equations) correspond to a grid {NX} columns and {NY}
    rows, with indices assigned row by row; that is, {S.col[k]= k mod
    NX}, {S.row[k] = k/NX}, and {S.uid[k] = k}, for {k} in {0..N-1}. */

pst_imgsys_t *pst_imgsys_from_eqs
  ( uint32_t N,
    pst_imgsys_equation_t *eq,
    int32_t NX,
    int32_t NY,
    int32_t *col,
    int32_t *row,
    int32_t *uid
  );
  /* Like {pst_imgsys_new} but uses the given vectors {eq,col,row,uid}
    instead of allocating a new one. The arrays {eq,col,row} must have
    {N} elements. The array {uid} must have {NX*NY} elements, and must be
    the inverse of {col,row} as explained under {pst_imgsys_t}. */

void pst_imgsys_check_valid(pst_imgsys_t *S, uint32_t N, int32_t NX, int32_t NY);
  /* Checks the basic consistency of the equation system. */

void pst_imgsys_free(pst_imgsys_t *S);
  /* Deallocates all storage used by {S}, including the index tables. */
      
bool_t pst_imgsys_equation_is_null(uint32_t uid, pst_imgsys_equation_t *eq, uint32_t N);
  /* Returns true iff the equation {eqk} is indeterminate,
    meaning that it is essentially {0*Z[uidk]=0} for its
    main unknown {uidk}.  Also does some consistency checks. */

void pst_imgsys_fill_holes(pst_imgsys_t *S);
  /*  Modifies every equation {eqk=S.eq[k]} that is indeterminate
    {0*Z[x,y]=0}, which should have with zero {eqk.wtot},
    so that it will set the corresponding unknown height {k} to the average of
    its neighors, with equal weights and an arbitrary positive {wtot}.
    
    The equations of those neighbors should not depend on {Z[x,y]}, and
    are not modified -- unless they too are indeterminate. The effect of
    this fudging is to fill in any indeterminate regions of the height
    map with a smooth surface that is attached at the edges to the
    determinate heights. */

void pst_imgsys_remove_holes(pst_imgsys_t *S);
  /* Removes from the system {S} any unknown {Z[k]=Z[x,y} whose equation {S->eq[k]}
   is indeterminate, namely {0*Z[k]=0}.  This will reduce the number {S->N}
   (but not {S->NX} and {S->NY}), renumber the equations and unknowns
   to the new range {0..S->N-1}, and set {S->uid[x + y*NX]} to {-1}. */

void pst_imgsys_copy_image_to_sol_vec(pst_imgsys_t *S, float_image_t *Z, double h[], double vdef);
  /* Sets the height value vector {h[0..S.N-1]} of system {S} to the samples of the
    height map {Z}. Specifically, sets {h[k]} to {Z[0,x,y]}, where {x} and {y} are {S->col[k]}
    and {S->row[k]}.
    
    However, if {x} and {y} are {-1} (meaning that the {h[k]} is
    not associated to an heigh map pixel), {h[k]} is set to {vdef}. */
    
void pst_imgsys_copy_sol_vec_to_image(pst_imgsys_t *S, double h[], float_image_t *Z, float vdef);
  /* Stores the height value vector {h[0..S.N-1]} of system {S} into the height map {Z}.
    Specifically, fills channel 0 of {Z} with {vdef}, then copies each height value {h[k}} 
    into {Z[0,x,y]}, where {x} and {y} are {S->col[k]} and {S->row[k]}.
    
    However, ignores {h[k]} if {x} and {y} are {-1} (meaning that the
    height value {h[k]} is not associated to an heigh map pixel). Note that
    samples of {Z} that are not associated with any height value height
    value of {S} will be left {vdef}. */
          
void pst_imgsys_extract_system_weight_image(pst_imgsys_t *S, float_image_t *W);
  /* Stores into channel 0 of {W} the equation weights {.wtot} of the system {S}.
    The map {W} must have {S.NX} cols and {S.NY} rows. */ 

/* DEBUGGING */
    
typedef void pst_imgsys_report_sys_proc_t(uint32_t level, pst_imgsys_t *S); 
  /* Type of a client-given procedure that may be called
    by recursive integrators to report the system used at each scale.
    Uses {col} and {row} to map indices of unknowns to pixel indices. */   

/* I/O */

#define pst_imgsys_FILE_TYPE "pst_imgsys_t"
#define pst_imgsys_FILE_VERSION "2024-01-07"
    
void pst_imgsys_write(FILE *wr, pst_imgsys_t *S);
  /* Writes the system {S} to stream {wr}.  
    
    The whole system is preceded by a line "begin {TYPE} (version of {VERSION})"
    and terminated by a line "end {TYPE}".
    
    The next three lines have "N = {S.N}", "NX = {S.NX}", and "NY = {S.NY}".
    
    Then follows one line for each unknown and equation, with the format 
    
      "{k} = {xk} {yk} w = {eqk.wtot} rhs = {eqk.rhs} nt = {eqk.nt} {TERMS}"
      
    where {k} is the unknow/equation index in {0..N-1}, {xk} is {S.row[k]},
    {yk} is {S.col[k]}, {eqk} us {S.eq[k]},
    and {TERMS} is a  list of {eqk.nt} pairs "{eqk.uid[j]} {eqk.cf[j]}"
    where {j} varies in {0..eqk.nt-1}.  The term's 
    unknown index {uidj = eqk.uid[j]} wil be followed by " = {x},{y}" if 
    {col[uidj]} and {row[uidj]} are not {-1}.  */

void pst_imgsys_write_report(pst_imgsys_t *S, char *filePrefix, int32_t level, char *tag, uint32_t indent);
  /* Writes the system {S} to a file called
    "{filePrefix}-{level}-{tag}.sys". If {tag} is null or empty the
    "-{tag}" is omitted. Uses {pst_imgsys_write}. Diagnostic messages
    are indented by {indent} spaces. */

#undef MAX_COEFFS
  /* Clients please use full name. */

#endif

/*
**
** Copyright (C) Jorge Stolfi, Unicamp.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appear in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty. Neither the author nor its employers are liable to
** any damages which may result from its use.
*/
