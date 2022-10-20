#ifndef dna_seq_view_tools_H
#define dna_seq_view_tools_H

/* Tools for DNA sequence visualization. */
/* Last edited on 2022-10-20 10:51:09 by stolfi */

#define _GNU_SOURCE
#include <GL/gl.h>
#include <GL/glu.h>
#include <stdint.h>

#include <bool.h>
#include <r3.h>
#include <frgb.h>

#include <msm_rung.h>

#include <dnae_seq.h>

void dna_seq_view_tools_draw_ball(r3_t* p, double rad, GLUquadric* quad);
  /* Draws a ball at {p} with radius {rad}. */

void dna_seq_view_tools_draw_stick(r3_t* p, bool_t ptrim, r3_t* q, bool_t qtrim, double rad, double trimrad, GLUquadric* quad);
  /* Draws a stick from {p} to {q} given radius {rad}.  If {ptrim} is true, trims the stick at the {p} end
    by {trimrad}.  Ditto if {trimq}, at the {q} end. */

void dna_seq_view_tools_draw_tetrahedron (frgb_t *color);
  /* Paints the tetrahedron. */

void dna_seq_view_tools_draw_sequence
  ( dnae_seq_t *z, 
    double magnify, 
    r3_t *pert, 
    double radius, 
    int32_t step,
    frgb_t *color, 
    int32_t ini, 
    int32_t fin
  );
  /* Paints the sequence {z}, scaling every point by {magnify}
    and adding the {pert} vector to it. Each datum whose index is a multiple of {step} is painted
    as a sphere with given {radius}.  The spheres and the connecting lines are
    painted with the given {color}, unless they are among the first {ini} or last {fin} datums. */

void dna_seq_view_tools_draw_paired
  ( msm_rung_vec_t *gv,
    dnae_seq_t *x,
    int32_t inix,
    int32_t finx,
    dnae_seq_t *y,
    int32_t iniy,
    int32_t finy,
    double magnify,
    bool_t perturb,
    double radius,
    int32_t step,
    frgb_t *color_x,
    frgb_t *color_y, 
    frgb_t *color_p,
    double dif_cutoff
  );
  /* Paints sequences {x} and {y} with colors {color_x} and {color_y}, and
    connected pairs informed by the pairing {gv} with sticks colored with {color_p}.
    the paramters {magnify}, {radius}, and {step} are those
    of {dna_seq_view_tools_draw_sequence}.
    The parameters {inix,finx} are {ini,fin} for the {x} sequence,
    and {iniy,finy} are the same for the {y} sequence. 
    If {perturb} is TRUE, displaces the sequences by slight
    perturbations, proportional to {radius}.  

    If {dif_cutoff} is non-negative, highlights the rungs that connect datums
    whose distance (before magnification and perturbation) is less than or
    equal to {dif_cutoff}. */

void dna_seq_view_tools_draw_rungs
  ( dnae_seq_t *x, 
    dnae_seq_t *y, 
    msm_rung_vec_t *gv, 
    double magnify,
    r3_t *x_pert,
    r3_t *y_pert,
    double radius, 
    frgb_t* color,
    double dif_cutoff
  );
  /* Draws the rungs {gv.e[0..gv.ne-1]} as rods of the given {color} 
    connecting the corresponding elements of sequences {x} and {y}.  Does 
    not draw the sequences themselves.  Assumes that the datums of {x} and {y}
    are drawn as spheres of given {radius}, after being scaled by 
    {magnify} and displaced by {x_pert} and {y_pert}, respectively.
    
    The {dif_cutoff} has the same meaning as in {dna_seq_view_tools_draw_paired}. */

#endif
