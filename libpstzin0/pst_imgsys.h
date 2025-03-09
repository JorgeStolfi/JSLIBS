/* imgsys - linear system solving for image-related problems. */

#ifndef pst_imgsys_H
#define pst_imgsys_H

/* Created on 2005-12-04 by Jorge Stolfi, unicamp, <stolfi@ic.unicamp.br> */
/* Based on the work of Rafael Saracchini, U.F.Fluminense. */
/* Last edited on 2025-03-03 03:47:56 by stolfi */
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
  { uint32_t nt;               /* {nt} is the number of terms in the equation. */
    uint32_t uid[MAX_COEFFS];  /* {uid[k]} is the index of some variable in the equation. */
    double cf[MAX_COEFFS];     /* Coefficient of that variable. */
    double rhs;                /* Right-hand side of equation. */
    double wtot;               /* Weight of the equation and of its main variable. */
  } pst_imgsys_equation_t;
  /* An {pst_imgsys_equation_t} record {eq} represents a linear equation 
    with at most {MAX_COEFFS} nonzero terms, namely
      {SUM {eq.cf[j]*z[eq.uid[j]] : j = 0..eq.nt-1} == rhs},
    where {z[j]} is the variable with index {j}.
    The variable with index {eq.uid[0]} is the main variable of the equation. */

typedef struct pst_imgsys_t
  { uint32_t N;                 /* Number of variables and equations. */
    pst_imgsys_equation_t* eq;  /* The equations are {eq[0..N-1]}. */
    /* Debugging info: */
    int32_t NX, NY;             /* Number of cols and rows of the assumed grid. */
    int32_t *col;               /* Maps index of a variable/equation to a {x} index, or -1. */
    int32_t *row;               /* Maps index of a variable/equation to a {y} index, pr -1. */
  } pst_imgsys_t;
  /* A {pst_imgsys_t} record {S} represents {N} linear equations
    {eq[0..N-1]} on {N} variables {z[0..N-1}. In a proper system, the
    main variable of equation {k} has index {k} --- that is,
    {eq[k].uid[0] == k} --- for all {k} in {0..N-1}.
    
    The tables {S.col,S.row} are used only when printing the
    system. Specifically, {S.col[0..N-1]} and {S.row[0..N-1]} define a
    partial one-to-one mapping from variable/equation index {k} in
    {0..N-1} to a point {(x,y)} of a subset of a 2D integer grid with
    {x} in {0..NX-1} and {y} in {0..NY-1}. If the variable {k} does not
    correspond to a grd point, then {S.col[k]=S.row[k]=-1}. */
    
pst_imgsys_t *pst_imgsys_new_grid(int32_t NX, int32_t NY);
  /* Creates a new linear system {S} with {N=NX*NY} equations on {N} variables.
    Each equation {eqk=S.eq[k]} intially has zero righ-hand side {eqk.rhs} and
    zero total weight {eqk.wtot}, and one term ({eqk.nt=1}) with variable 
    {eqk.uid[0]=k} and coefficient {eqk.cf[0]=0}.   Namely, the equation is
    indeterminate, {0*Z[k]=0}.
    
    The tables {S.col,S.row,S.uid} are filled assuming that the variables
    (and therefore equations) correspond to a grid {NX} columns and {NY}
    rows, with indices assigned row by row; that is, {S.col[k]= k mod
    NX}, {S.row[k] = k/NX}, and {S.uid[k] = k}, for {k} in {0..N-1}. */

pst_imgsys_t *pst_imgsys_from_eqs
  ( uint32_t N,
    pst_imgsys_equation_t *eq,
    int32_t *col,
    int32_t *row
  );
  /* Like {pst_imgsys_new} but uses the given vectors {eq,col,row,uid}
    instead of allocating a new one. The arrays {eq,col,row} must have
    {N} elements. */

void pst_imgsys_check_valid(pst_imgsys_t *S, uint32_t N, int32_t NX, int32_t NY);
  /* Checks the basic consistency of the equation system.
    The parameters {NX,NY} may be {-1}; if not, the {S.col} and {S.row}
    entries must be in {0..NX-1]} and {0..NY-1}, respectively. */

void pst_imgsys_free(pst_imgsys_t *S);
  /* Deallocates all storage used by {S}, including the index tables. */
      
bool_t pst_imgsys_equation_is_null(uint32_t uid, pst_imgsys_equation_t *eq, uint32_t N);
  /* Returns true iff the equation {eqk} is indeterminate,
    meaning that it is essentially {0*Z[uidk]=0} for its
    main variable {uidk}.  Also does some consistency checks. */
    
float_image_t* pst_imgsys_make_weight_image(pst_imgsys_t *S);
  /* Returns a single-channel image {SW} with {S.NX} cols and {S.NY} rows where pixel 
    {SW[X,Y]} is the reliability weight of the equation that defines {Z[X,Y]},
    as implied by the {S.col} and {S.row} tables. */

void pst_imgsys_copy_image_to_sol_vec(pst_imgsys_t *S, float_image_t *Z, double z[], double vdef);
  /* Sets the height value vector {z[0..S.N-1]} of system {S} to the samples of channel 0 the
    height map {Z}. Specifically, sets {z[k]} to {Z[0,x,y]}, where {x} and {y} are {S->col[k]}
    and {S->row[k]}.  Any other channels of {Z} are ignored.
    
    However, if {x} and {y} are {-1} (meaning that the {z[k]} is
    not associated to an heigh map pixel), {z[k]} is set to {vdef}. */
    
void pst_imgsys_copy_sol_vec_to_image(pst_imgsys_t *S, double z[], float_image_t *Z, float vdef);
  /* Stores the height value vector {z[0..S.N-1]} of system {S} into the height map {Z}.
    Specifically, fills channel 0 of {Z} with {vdef}, then copies each height value {z[k}} 
    into {Z[0,x,y]}, where {x} and {y} are {S->col[k]} and {S->row[k]}.  
    
    If {Z} has two or more chanels, also initialized channel 1 with zeros, 
    then copies the total weight of each equation {S.eq[k]} to {Z[1,x,y]}.
    
    However, skips the assignment if {x} and {y} are {-1} (meaning that the
    height value {z[k]} is not associated to an heigh map pixel). */

/* DEBUGGING */
    
typedef void pst_imgsys_report_sys_proc_t(int32_t level, pst_imgsys_t *S); 
  /* Type of a client-given procedure that may be called
    by recursive integrators to report the system used at each scale. */   

/* I/O */

#define pst_imgsys_FILE_TYPE "pst_imgsys_t"
#define pst_imgsys_FILE_VERSION "2025-03-04"
    
void pst_imgsys_write(FILE *wr, pst_imgsys_t *S, char *fmt);
  /* Writes the system {S} to stream {wr}.  
    
    The whole system is preceded by a line "begin {TYPE} (version of {VERSION})"
    and terminated by a line "end {TYPE}".
    
    The next three lines have "N = {S.N}", "NX = {S.NX}", and "NY = {S.NY}".
    
    Then follows one line for each variable and equation, with the format 
    
      "{k} = {xk} {yk} w = {eqk.wtot} rhs = {eqk.rhs} nt = {ntk} {TERMS}"
      
    where {k} is the unknow/equation index in {0..N-1}, {xk} is
    {S.row[k]}, {yk} is {S.col[k]}, {eqk} us {S.eq[k]}, {ntk=eqk.nt} is
    the number of left-hand terms, and {TERMS} is a list of {ntk}
    entries "{eqk.cf[j]}*Z[{eqk.uid[j]}]" where {j} varies in
    {0..eqk.nt-1}. The term's variable index {uidj = eqk.uid[j]} will be
    followed by "={x},{y}" if {x=col[uidj]} and {y=row[uidj]} are not
    {-1}.  
    
    The string {fmt} must be a valid float format spec with forced sign,
    such as "%+10.8f" or "%+24.16e", and is used to print the right-hand
    side {eqk.rhs} and each coefficient {eqk.cf[j]}. */

void pst_imgsys_write_named(char *fileName, pst_imgsys_t *S, char *fmt, int32_t indent);
  /* Writes the system {S} to a the file "{fileName}".
    Messages to {stderr} are indented by {indent} spaces. */

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
