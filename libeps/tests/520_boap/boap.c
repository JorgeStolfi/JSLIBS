/* Tech drawing for the front gate of the storage area of Boaretto da Silva 113 */
/* Last edited on 2024-12-05 10:15:26 by stolfi */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>

#include <r2.h>
#include <r3.h>
#include <bool.h>
#include <sign.h>
#include <affirm.h>
#include <jsfile.h>

#include <epswr.h>
#include <epswr_dim.h>

#include <boap_dim.h>
#include <boap_bar.h>
#include <boap_elems.h>
#include <boap_model.h>

typedef struct boap_iron_dims_t
  { /* Big 'L' irons: */
    double L_wd;       /* Width in both cross-sesction axes. */
    double L_th;       /* Thickness in both cross-sesction axes. */
    frgb_t L_color;
    /* Vitral 'T' irons: */
    double T_wd0;      /* Width of the "arms" of the cross-section. */
    double T_th0;      /* Thickness of the "arms". */
    double T_th1;      /* Thickness of the "leg" including arms thickness. */
    double T_wd1;      /* Width of the "leg" including arms thickness. */
    frgb_t T_color;
    /* Vitral flat irons: */
    double F_wd;       /* Width of the cross-section. */
    double F_th;       /* Thicknes of cross-section. */
    frgb_t F_color;
    /* Big metalons: */
    double Mbig_wd;    /* Width in both cross-sesction axes. */
    double Mbig_th;    /* Thickness in both cross-sesction axes. */
    frgb_t Mbig_color;
    /* Small metalons: */
    double Msma_wd;    /* Width in both cross-sesction axes. */
    double Msma_th;    /* Thickness in both cross-sesction axes. */
    frgb_t Msma_color;
  }  boap_iron_dims_t;
  /* Dimensions and colors of stock iron bar types. */
  
boap_iron_dims_t *boap_iron_dims_define(void);
  /* Defines the dimensions of iron bars to be used. */

void boap_draw_iron_cross_sections(boap_iron_dims_t *ird);
  /* Writes an EPS file "out/porch_irons.eps" with the cross-sections
    of the irons used in the porch. */

void boap_draw_front_wall(boap_iron_dims_t *ird, char *pname, bool_t both);
  /* Writes EPS files called "out/porch_{pname}_{subname}.eps" with the drawings 
    of the front wall and gate of the porch. The origin is at the left 
    bottom corner of the gate, at the back of the main tube frame.
    
    If {both} is {TRUE}, shows both together without dimensions,
    with {subname} equal to "both".  If {both} is false,
    else shows them in separate figures, with {subname} equal to "skel" or
    "pans", with dimensions.  */

void boap_draw_rear_wall(boap_iron_dims_t *ird, char *pname, bool_t both);
  /* Writes EPS files called "out/porch_{pname}_{n}.eps" with the drawings
    of the rear wall and gate of the porch. The World coordinate axes are: {X}
    horizontal towards the street, {Y} horizontal along the street (left
    to right seen from {+X}), {Z} up. The origin is at the left 
    bottom corner of the gate, at the back of the main tube frame.  */

void boap_draw_side_wall(boap_iron_dims_t *ird, char *pname, bool_t both);
  /* Writes EPS files called "out/porch_{pname}_{n}.eps" with the drawings
   window of the porch. The World coordinate axes are: {X}
   x.horizontal towards the house, {Y} horizontal towards the street
   (left to right seen from {+X}), {Z} up. The origin is at the left
   bottom corner of the window, at the back of the main tube frame. */

void boap_draw_roof(boap_iron_dims_t *ird, char *pname, bool_t both);
  /* Writes EPS files called "out/porch_{pname}_{n}.eps" with the drawings
    of the roof of the porch. The World coordinate axes are: {X}
    horizontal towards the street, {Y} horizontal along the street (left
    to right seen from {+X}), {Z} up. The origin is at the left 
    bottom corner of the roof frame, at the back of the main tube frame.  */

boap_model_t *boap_make_front_skeleton
  ( boap_iron_dims_t *ird, 
    double p_Xmin, double p_Xmax, 
    double p_Ymin, double p_Ymax,
    double p_Zmin, double p_Zmax,
    double s_Ymin, double s_Ymax,
    double s_Zmax,
    bool_t dims
  );
  /* Returns a geometric model with the "skeleton" of the front wall,
    consisting mostly of 'L' irons creating six gaps where the door and
    panels are to be inserted.
    
    The skeleton spans the box {[p_Xmin _ p_Xmax] × [p_Ymin _ p_Ymax] × [p_Zmin _ p_Zmax]}.
    
    The main horz 'L' bars are snug below {s_Zmax} and {p_Zmax}.
    
    There are vertical 'L' bars snug right of {p_Ymin} and {s_Ymin}, and
    snug left of {s_Ymax} and {p_Ymax}, stretching the whole height of
    the skeleton.
    
    The 'L' bars have width {L_XYZwd} in both transversaldirections, and thickness {L_XYZth}. */

boap_model_t *boap_make_front_panels
  ( boap_iron_dims_t *ird, 
    double pans_Xmin, double pans_Xmax,
  
    int8_t lpan_Vbar_N, double lpan_Ymin, double lpan_Ymax,
    int8_t door_Vbar_N, double door_Ymin, double door_Ymax,
    int8_t rpan_Vbar_N, double rpan_Ymin, double rpan_Ymax,
    double glass_Ysize,
    
    int8_t bpan_Hbar_N, double bpan_Zmin, double bpan_Zmax,
    int8_t tpan_Hbar_N, double tpan_Zmin, double tpan_Zmax,
    double glass_Zsize, double hole_Zsize,
    
    bool_t dims
  );
  /* Returns a model for the panels of the front wall.
  
  The parameters {glass_Ysize,glass_Zsize} are the dimensions of a
  typical glass pane in the bottom panel vitrals. The paramter
  {hole_Zsize} is the height of a typical hole ofthe top panel open
  grilles.
  
  The parameters {lpan_Vbar_N,door_Vbar_N,rpan_Vbar_N} are the counts of
  vertical bars in the left, center, and right panels. The parameters
  {bpan_Hbar_N,tpan_Hbar_N} are the counts of horizontal bars in the
  even-indexed columns of the bottom and top vitrals/grilles.
  
  The parameters {lpan_Ymin,lpan_Ymax,...tpan_Zmax} are the coord ranges
  of the panels, including their metalon frames but excluding the gaps
  between panel and skeleton.
  
  The other parameters are as in {boap_add_front_panel}. */

void boap_add_front_panel
  ( boap_model_t *mod, 
    int8_t ixp, 
    int8_t iyp,
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    int8_t NV, double Vbar_Ystep,
    int8_t NH, double Hbar_Zstep,
    bool_t dims
  );
  /* Adds to {mod} one of the six parts of the front wall that are insterted
    into the skeleton, given its coordinate ranges.

    The index {ixp} can be 0 for the left window and the grille above it,
    1 for the door and the grille above it, 2 for the right window
    and the grille above it.
    
    The index {iyp} is 0 for the door and side windows, and 1 for the 
    grilles above them. 
    
    The panel will span the box {[Xmin _ Xmax] × [Ymin _ Ymax] × [Zmin _
    Zmax]}. The panel will have a frame made of metalon with {Mwd_big ×
    Mwd_big} square cross section. Inside it there may be a pleated "skirt"
    and a grille, depending on {iyp}.
    
    If {iyp} is 0, there will be a skirt, and the grille will be
    designed to be the support of a vitral. It will consist of 'T' irons
    with an extra framelet of flat irons all around, inside the main
    frame. The parameters {T_YZwd,F_YZwd,T_YXth,FT_Xth} are like in
    {boap_add_front_panel_vitral}. The parameter {Mwd_sma} will be
    ignored.
    
    If{iyp} is 1, the grille will be open, consisting of 
    metalons with cross-section {Mwd_sma × Mwd_sma}.
    The parameters {T_YZwd,F_YZwd,T_YZth,FT_Xth} will be ignored.
    
    In any case, there will be {NV} vertical bars with their centers
    spaced {Vbar_Ystep} apart. There will be {NH} horizontal bars in
    even columns, with centers spaced {Hbar_Zstep} apart. The bars in
    odd columns will be vertically offset by {Hbar_Zstep/2}. */

void boap_add_front_panel_frame
  ( boap_model_t *mod, 
    int8_t which,
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    double Zdisp,
    bool_t dims
  );
  /* Adds to {mod} the frame of a panel on the front wall, on the left
    ({which=0}) or right ({which = 1}) side of the door. The position
    and dimensions are of the frame. The frame consists of metalons
    all around, and a horz metalon displaced {Zdist} up from the bottom. */

void boap_add_front_panel_skirt
  ( boap_model_t *mod, 
    int8_t ixp,
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    bool_t dims
  );
  /* Adds to {mod} the "pleated skirt" in the bottom half of a panel on the front wall.
    The parameter {ixp} is 0 for the left window, 1 for the door, and 2 for the 
    right window. The position and dimensions EXCLUDE the frame.  */

void boap_add_front_panel_vitral
  ( boap_model_t *mod, 
    int8_t ixp,
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    int8_t NV, double Vbar_Ystep, 
    int8_t NH, double Hbar_Zstep,
    bool_t dims
  );
  /* Adds to {mod} the "vitral" grille in the top half of a side panel or of the
    door. The parameter {ixp} is 0 for the left window, 1 for the door,
    and 2 for the right window. The position and dimensions EXCLUDE the
    frame.  

    The vitral is a grille of 'T' bars with narrow flat bars all
    around, meant to receive glass panes.
    
    The width of the 'T' irons and 'F' irons in the {Y} and {Z}
    directions will be {T_YZwd} and {F_YZwd}, respectively.  The thickness 
    of the  'T' irons in the {Y} and {Z} directions will be {T_YZth},
    and the common thickness of all irons in the {X} direction will be {FT_Xth}.
    
    There will be {NV} vertical bars with their centers 
    spaced {Vbar_Ystep} apart. There  will be {NH} horizontal bars
    in even columns, with centers spaced {Hbar_Zstep} apart.  The bars in odd columns
    will be vertically offset by {Hbar_Zstep/2}. */

void boap_add_front_panel_open_grille 
  ( boap_model_t *mod, 
    int8_t ixp,
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    int8_t NV, double Vbar_Ystep, 
    int8_t NH, double Hbar_Zstep,
    bool_t dims
  );
  /* Adds to {mod} an open grille in a panel. The parameter {ixp} is 0 for the 
    grille over the left window, 1 for the one over the door,
    and 2 for the one over the right window. The position and dimensions EXCLUDE the
    frame. 
    
    The grille is made of metalons with square {Mwd × Mwd}
    cross-section. The difference {Xmax-Xmin} must be equal to {Mwd}.
    The difference {Xmax-Xmin} must be equal to {Mwd}.
    
    There will be {NV} vertical bars with their centers 
    spaced {Vbar_Ystep} apart. There  will be {NH} horizontal bars
    in even columns, with centers spaced {Hbar_Zstep} apart.  The bars in odd columns
    will be vertically offset by {Hbar_Zstep/2}. */
  
void boap_add_grille
  ( boap_model_t *mod, 
    char type,
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    int8_t NV, double Vbar_Ystep, 
    int8_t NH, double Hbar_Zstep,
    bool_t dims
  );
  /* Adds to {mod} a grille with staggered bars,
    The grille will fit in the specified space. 
  
    If {type} is 'M', the grille will consist of metalon tubes with
    square cross-section {Mwd × Mwd}, and {T_YZwd,F_YZwd,T_YZth,FT_Xth} will be ignored.
    
    If {type} is 'T', the grille will consist of 'T'-bars and there will
    be an "inner frame" of flat irons. The parameters
    {T_YZwd,F_YZwd,T_YXth,FT_Xth} are like in {boap_add_front_panel_vitral}.
    The parameter {Mwd} will be ignored.
    
    There will be {NV} vertical bars with centers spaced {Vbar_Ystep}
    apart. The first and last gap will be equal but may have some other
    width.
    
    There will be {NV} vertical bars with their centers 
    spaced {Vbar_Ystep} apart. There  will be {NH} horizontal bars
    in even columns, with centers spaced {Hbar_Zstep} apart.  The bars in odd columns
    will be vertically offset by {Hbar_Zstep/2}. */

void boap_add_rect_tube_frame
  ( boap_model_t *mod, 
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    bool_t dims,
    frgb_t *color
  );
  /* Adds to {mod} the frame of a panel that fits 
    snugly inside the box {[Xmin _ Xmax]}  by {[Ymin _ Ymax]}
    by {[Zmin _ Zmax]}.  The bounds of each coordinate can be given 
    in any order.  
    
    Assumes that the tube cross sections are rectangular, measuring
    {M_Xwd} in the {X} direction and {M_YZwd} in the {Y} or {Z}
    driections. */

void boap_add_flat_frame
  ( boap_model_t *mod, 
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    bool_t dims
  );
  /* Adds to {mod} a rectangular frame in the given box. The flame 
    consists of flat irons of width {F_YZwd} in the {Y} and {Z}
    directions and thickness {F_Xth}  in the {X} direction. */
    
void boap_add_rect_tube
  ( boap_model_t *mod, 
    int8_t axC,
    double Xctr, double Yctr, double Zctr, 
    double Xsize, double Ysize, double Zsize,
    frgb_t *color
  );
  /* Adds to {mod} a tube with rectangular cross-section given its center
    and extent in each coordinate direction.  The tube's main axis line is
    assumed to be parallel to coordinate axis {axC} 
    ({0=X}, {1=Y}, {2=Z}). */

void boap_add_bar
  ( boap_model_t *mod, 
    char type,
    double Xctr, double Yctr, double Zctr, 
    int8_t axC, double Csize, 
    int8_t axA, double Awd, double Athk,
    int8_t axB, double Bwd, double Bthk,
    frgb_t *color
  );
  /* Adds to {mod} an 'L' or 'T' bar, as specified by the character {type}.
    
    The "center" of the bar is at {(Xctr,Yctr,Zctr)}. For an 'L' bar,
    the center is on the edge that corresponds to the outer corner of
    the 'L'. For a 'T' bar, the center is on the midline of the flat
    face that corresponds to the the top edge of the 'T'.
    
    Either way, the bar onsists of two flat plates that extend
    parallel to axis {axC} by a total length {Csize}.  
    
    One of the plates corresponds to the horizontal side of the 'L'
    or 'T' and is {Awd} wide in the direction {axA} and {Bthk}
    thick in the direction {axB}. 
    
    The other plate  corresponds to the {B} side of the 'L' or 'T' 
    and is {Bwd} wide in the direction {axB} and {Athk} thick in 
    the direction {axA}.
    
    The parameters {Awd,Athk,Bwd,Bthk} can be negative to indicate
    direction opposite to the respective axis. However, {Awd} and {Athk}
    must have the same sign, and {Bwd} and {Bthk} must have the same
    sign. */

void boap_add_bar_for_cross_section
  ( boap_model_t *mod,
    char type,
    double Yctr, double Zctr,
    double wd0, double th0,
    double wd1, double th1,
    frgb_t *color
  );
  /* Adds the bar to the cross-section sampler. */

void boap_add_dim
  ( boap_model_t *mod, 
    double xa, double ya,
    double xb, double yb,
    double agap, double alen,
    double elen,
    bool_t inner
  );
  /* Adds to {mod} a dimension element that shows the distance between 
    points {a=(xa,ya)} and {b=(xb,yb)} as a label with extension and dimension lines.  
    
    The parameters {agap,alen} define the start and length (mm) of the
    extension line for the point {a}. The parameters {bgap,blen}
    have the same meaning for {b}. If {inner}, the dimension line
    is drawn between the extension lines, otherwise outside them. */

frgb_t *boap_make_color(double R, double G, double B);
  /* Allocates a new {frgb_t} record and fills it with the given
    RGB color. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    boap_iron_dims_t *ird = boap_iron_dims_define(); 
    boap_draw_iron_cross_sections(ird);
    boap_draw_front_wall(ird, "front", TRUE);
    boap_draw_front_wall(ird, "front", FALSE);
    /* boap_draw_rear_wall(ird, "rear", TRUE); */
    /* boap_draw_roof(ird, "roof", TRUE); */
    return 0;
  }
  
boap_iron_dims_t *boap_iron_dims_define(void)
  { 
    boap_iron_dims_t *ird = notnull(malloc(sizeof(boap_iron_dims_t)), "no mem");

    ird->Msma_wd = 20.0;   /* Width of small metalon. */
    ird->Msma_th =  1.2;   /* Thickness of small metalon. */
    ird->Msma_color = (frgb_t){{ 0.600f, 0.300f, 1.000f }};

    ird->Mbig_wd = 30.0;   /* Width of big metalon. */
    ird->Mbig_th =  1.2;   /* Thickness of small metalon. */
    ird->Mbig_color = (frgb_t){{ 1.000f, 0.300f, 0.600f }};

    ird->L_wd = 32.0;      /* Width of skeleton 'L' bars. */
    ird->L_th = 3.0;       /* Thickness of 'L' bars. */
    ird->L_color = (frgb_t){{ 1.000f, 0.750f, 0.200f }};

    ird->T_wd0 = 23.0;       /* Width of "arms" of 'T' bars. */
    ird->T_th0 = 2.0;        /* Thickness of "arms" of 'T'. */
    ird->T_wd1 = 23.0;       /* Length of "leg" of 'T' bars. !!! Check! !!! */
    ird->T_th1 = 3.0;        /* Thickness of "leg" part of 'T' irons. */
    ird->T_color = (frgb_t){{ 0.600f, 8.000f, 0.200f }};

    ird->F_wd =  9.0;       /* Width of flat ('F') bars. */
    ird->F_th = 2.0;        /* Thickness of 'T' and 'F' bars in the {X} direction. */
    ird->F_color = (frgb_t){{ 0.400f, 0.800f, 1.000f }};
    
    return ird;
  }

void boap_draw_iron_cross_sections(boap_iron_dims_t *ird)
  {
    boap_model_t *mod = boap_model_new();
    
    double Yc, Zc;  /* Bar center. */
    double wd0, th0;    /* Width and thickness of part 0 of cross-section. */
    double wd1, th1;    /* Width and thickness of part 1 of cross-section. */

    /* Big metalon: */
    Yc = 0.0; Zc = 0.0;
    wd0 = wd1 = ird->Mbig_wd; th0 = th1 = ird->Mbig_th;
    boap_add_bar_for_cross_section(mod, 'M', Yc,Zc, wd0, th0, wd1, th1, &(ird->Mbig_color));

    /* Small metalon: */
    Yc = 100.0; Zc = 0.0;
    wd0 = wd1 = ird->Msma_wd; th0 = th1 = ird->Msma_th;
    boap_add_bar_for_cross_section(mod, 'M', Yc,Zc, wd0, th0, wd1, th1, &(ird->Msma_color));

    /* 'L' iron: */
    Yc = 0.0; Zc = 100.0;
    wd0 = wd1 = ird->L_wd; th0 = th1 = ird->L_th;
    boap_add_bar_for_cross_section(mod, 'L', Yc,Zc, wd0, th0, wd1, th1, &(ird->L_color));

    /* 'T' iron: */
    Yc = 100.0; Zc = 100.0;
    wd0 = ird->T_wd0; th0 = ird->T_th0; wd1 = ird->T_wd1; th1 = ird->T_th1;
    boap_add_bar_for_cross_section(mod, 'T', Yc,Zc, wd0, th0, wd1, th1, &(ird->T_color));

    /* 'F' iron: */
    Yc = 0.0; Zc = 200.0;
    wd0 = ird->F_wd; th0 = ird->F_th; wd1 = NAN; th1 = NAN;
    boap_add_bar_for_cross_section(mod, 'F', Yc,Zc, wd0, th0, wd1, th1, &(ird->F_color));
    
    double tot_Ymin = -25, tot_Ymax = +125;
    double tot_Zmin = -45, tot_Zmax = +255;

    epswr_figure_t *epsf = boap_new_figure("irons", NULL, tot_Ymin, tot_Ymax, tot_Zmin, tot_Zmax, FALSE);
    boap_model_draw(epsf, mod, 0, 0.0);
    epswr_end_figure(epsf);
  } 
  
void boap_add_bar_for_cross_section
  ( boap_model_t *mod,
    char type,
    double Yctr, double Zctr,
    double wd0, double th0,
    double wd1, double th1,
    frgb_t *color
  )
  {
    double Xctr = 0.0; 
    double Xsize = 100;
    boap_add_bar
      ( mod, type,
        Xctr, Yctr, Zctr,
        0, Xsize,
        1, wd0, th1,
        2, wd1, th0,
        color
      );
    
    /* Show dimensions: */
    double Xa = Xctr, Xb = Xctr;
    double Ya_wd, Za_wd, Zb_wd, Yb_wd;
    double Ya_th, Za_th, Zb_th, Yb_th;
    double agap_wd, bgap_wd, elen_wd;
    double agap_th, bgap_th, elen_th;

    { /* Show dimensions of part 0 of cross-section: */
      agap_wd = bgap_wd = -1.0; elen_wd = -7.0;
      agap_th = bgap_th = +1.0; elen_th = +7.0;
      if (type == 'M')
        { Ya_wd = Yctr - wd0/2; Yb_wd = Yctr + wd0/2;
          Za_wd = Zctr - wd1/2; Zb_wd = Za_wd;
          Ya_th = Yctr + wd0/2; Yb_th = Ya_th;
          Za_th = Zctr - wd1/2 + th0; Zb_th = Zctr - wd1/2;
        } 
      else if (type == 'L')
        { Ya_wd = Yctr; Yb_wd = Yctr + wd0;
          Za_wd = Zctr; Zb_wd = Za_wd;
          Ya_th = Yctr + wd0; Yb_th = Ya_th;
          Za_th = Zctr + th0; Zb_th = Zctr;
        } 
      else if (type == 'T')
        { Ya_wd = Yctr - wd0/2; Yb_wd = Yctr + wd0/2;
          Za_wd = Zctr; Zb_wd = Za_wd;
          Ya_th = Yctr + wd0/2; Yb_th = Ya_th;
          Za_th = Zctr + th0; Zb_th = Zctr;
        }
      else if (type == 'F')
        { Ya_wd = Yctr - wd0/2; Yb_wd = Yctr + wd0/2;
          Za_wd = Zctr - th0/2; Zb_wd = Za_wd;
          Ya_th = Yctr + wd0/2; Yb_th = Ya_th;
          Za_th = Zctr + th0/2; Zb_th = Zctr - th0/2;
        }
      else
        { assert(FALSE); }
      boap_model_add_dim(mod, boap_dim_new(Xa, Ya_wd, Za_wd, Xb, Yb_wd, Zb_wd, agap_wd, bgap_wd, elen_wd, TRUE, 0));
      boap_model_add_dim(mod, boap_dim_new(Xa, Ya_th, Za_th, Xb, Yb_th, Zb_th, agap_th, bgap_th, elen_th, FALSE, 1));
    }
    
    if (type != 'F')
      { /* Show dimensions of part 1 of cross-section: */
        agap_wd = bgap_wd = -1.0; elen_wd = -7.0;
        agap_th = bgap_th = +1.0; elen_th = +7.0;
        if (type == 'M')
          { Ya_wd = Yctr - wd0/2; Yb_wd = Ya_wd;
            Za_wd = Zctr + wd1/2; Zb_wd = Zctr - wd1/2;
            Ya_th = Yctr - wd0/2; Yb_th = Ya_th + th1;
            Za_th = Zctr + wd1/2; Zb_th = Za_th;
          } 
        else if (type == 'L')
          { Ya_wd = Yctr; Yb_wd = Ya_wd;
            Za_wd = Zctr + wd1; Zb_wd = Zctr;
            Ya_th = Yctr; Yb_th = Ya_th + th1;
            Za_th = Zctr + wd1; Zb_th = Za_th;
          } 
        else if (type == 'T')
          { Ya_wd = Yctr - th1/2; Yb_wd = Ya_wd;
            Za_wd = Zctr + wd1; Zb_wd = Zctr;
            agap_wd = -1.0;  bgap_wd = -1.0 - wd0/2; 
            Ya_th = Yctr - th1/2; Yb_th = Ya_th + th1;
            Za_th = Zctr + wd1; Zb_th = Za_th;
          }
        else
          { assert(FALSE); }
        boap_model_add_dim(mod, boap_dim_new(Xa, Ya_wd, Za_wd, Xb, Yb_wd, Zb_wd, agap_wd, bgap_wd, elen_wd, TRUE, 0));
        boap_model_add_dim(mod, boap_dim_new(Xa, Ya_th, Za_th, Xb, Yb_th, Zb_th, agap_th, bgap_th, elen_th, FALSE, 1));
      }
  }

void boap_draw_front_wall(boap_iron_dims_t *ird, char *pname, bool_t both)
  {
    /* Wall position and dimensions: */

    double tot_Xpos = 0.0; 
    double tot_Ypos = 0.0, tot_Ysize = 1771.0; /* Space is 1710; leave 2cm gap. */
    double tot_Zpos = 0.0, tot_Zsize = 3000.0;  

    double tot_Xmin = tot_Xpos;
    double tot_Ymin = tot_Ypos, tot_Ymax = tot_Ymin + tot_Ysize;
    double tot_Zmin = tot_Zpos, tot_Zmax = tot_Zmin + tot_Zsize;
    
    /* Skeleton: */
   
    double pan_Xgap = 1.0;  /* {X} gap between panel and frame. */
    double pan_Ygap = 5.0;  /* {Y} gap between panel and frame. */
    double pan_Zgap = 4.0;  /* {Z} gap between panel and frame. */
    
    /* Number of vertical bars in panels: */
    int8_t lpan_Vbar_N = 0; /* Left panels. */
    int8_t door_Vbar_N = 6; /* Door panel and its open grille. */
    int8_t rpan_Vbar_N = 4; /* Right panels. */
    
    /* Compute width {glass_Ysize} of glass panes (tight): */
    /* Width of left panel:  {LPwd = 2*ird->Mbig_wd + (lpan_Vbar_N + 1)*u + lpan_Vbar_N*ird->T_th1}. */
    /* Width of door leaf:   {DPwd = 2*ird->Mbig_wd + (door_Vbar_N + 1)*u + door_Vbar_N*ird->T_th1}. */
    /* Width of right panel: {RPwd = 2*ird->Mbig_wd + (rpan_Vbar_N + 1)*u + rpan_Vbar_N*ird->T_th1}. */
    /* where {u} is the width (tight) of a glass pane. */
    /* Total width: {tot_Ysize = LPwd + DPwd + RPwd + 4*ird->L_th + 6*pan_Ygap}. */
    /* Hence: */
    
    double LPA = lpan_Vbar_N + 1, LPB = 2*ird->Mbig_wd + lpan_Vbar_N*ird->T_th1;
    double DPA = door_Vbar_N + 1, DPB = 2*ird->Mbig_wd + door_Vbar_N*ird->T_th1;
    double RPA = rpan_Vbar_N + 1, RPB = 2*ird->Mbig_wd + rpan_Vbar_N*ird->T_th1;
    
    double glass_Ysize = (tot_Ysize - 4*ird->L_th - 6*pan_Ygap - (LPB + DPB + RPB))/(LPA + DPA + RPA);
    double glass_Zsize = 295;   /* Height of glass panes (tight fit) in vitrals  of bot panels. */
    double hole_Zsize = 246;    /* Height of holes in open grille of top panels. */
    
    fprintf(stderr, "TYPICAL GLASS PANE: %.2f x %.2f mm\n", glass_Ysize, glass_Zsize);

    /* ------------------------------------------------------------ */
    /* {Y} sizes of panels: */
    
    double lpan_Ysize = 2*ird->Mbig_wd + (lpan_Vbar_N+1)*glass_Ysize + lpan_Vbar_N*ird->T_th1;
    double door_Ysize = 2*ird->Mbig_wd + (door_Vbar_N+1)*glass_Ysize + door_Vbar_N*ird->T_th1;
    double rpan_Ysize = 2*ird->Mbig_wd + (rpan_Vbar_N+1)*glass_Ysize + rpan_Vbar_N*ird->T_th1;

    /* {Y} coord ranges of panels: */
    double lpan_Ymin = tot_Ymin + ird->L_th + pan_Ygap;
    double lpan_Ymax = lpan_Ymin + lpan_Ysize;
    
    double door_Ymin = lpan_Ymax + pan_Ygap + ird->L_th + pan_Ygap;   
    double door_Ymax = door_Ymin + door_Ysize;

    double rpan_Ymin = door_Ymax + pan_Ygap + ird->L_th + pan_Ygap;
    double rpan_Ymax = rpan_Ymin + rpan_Ysize;
    
    assert(fabs((rpan_Ymax + pan_Ygap + ird->L_th) - tot_Ymax) < 0.00001);
    
    /* ------------------------------------------------------------ */
    /* {Z} sizes and ranges of bottom and top panels: */
    
    double door_Zsize = 2180.0; /* Height of door leaf. */

    double bpan_Zmin = tot_Zmin;
    double bpan_Zsize = door_Zsize;
    double bpan_Zmax = bpan_Zmin + bpan_Zsize;
    
    double tpan_Zmin = bpan_Zmax + pan_Zgap + ird->L_th + pan_Zgap;
    double tpan_Zmax = tot_Zmax - ird->L_th - pan_Zgap;
    /* double tpan_Zsize = tpan_Zmax - tpan_Zmin; */

    /* Number of horizontal bars in panels: */
    int8_t bpan_Hbar_N = 3; /* Bottom panels. */
    int8_t tpan_Hbar_N = 2; /* Top panels. */

    /* {Y} span of door jambs: */
    double jamb_Ymin = door_Ymin - (pan_Ygap + ird->L_th);
    double jamb_Ymax = door_Ymax + (pan_Ygap + ird->L_th);
    double jamb_Zmax = bpan_Zmax + (pan_Zgap + ird->L_th);
    
    /* Master frame: */
    double skel_Xmin = tot_Xmin - (pan_Xgap + ird->L_th);
    double skel_Xmax = skel_Xmin + ird->L_wd; 
    double skel_Ymin = tot_Ymin;
    double skel_Ymax = tot_Ymax;  
    double skel_Zmin = tot_Zmin;
    double skel_Zmax = tot_Zmax;  
    
    bool_t dims = (! both);
    
    boap_model_t *skel_mod = boap_make_front_skeleton
      ( ird,
        skel_Xmin, skel_Xmax,
        skel_Ymin, skel_Ymax,
        skel_Zmin, skel_Zmax,
        jamb_Ymin, jamb_Ymax, 
        jamb_Zmax,
        dims
      );

    double pans_Xmin = skel_Xmin + ird->L_th + pan_Xgap;
    double pans_Xmax = pans_Xmin + ird->Mbig_wd;
    
    boap_model_t *pans_mod = boap_make_front_panels
      ( ird, 
        pans_Xmin, pans_Xmax,
        
        lpan_Vbar_N, lpan_Ymin, lpan_Ymax,
        door_Vbar_N, door_Ymin, door_Ymax,
        rpan_Vbar_N, rpan_Ymin, rpan_Ymax,
        glass_Ysize,

        bpan_Hbar_N, bpan_Zmin, bpan_Zmax,
        tpan_Hbar_N, tpan_Zmin, tpan_Zmax,
        glass_Zsize, hole_Zsize,
        
        dims
      );

    /* Plot window ranges: */
    if (both)
      { /* Draws the two models on the same figure: */
        /* !!! Should merge the models. !!! */
        assert(! dims);
        double fig_Ymin = tot_Ymin - 150;
        double fig_Ymax = tot_Ymax + 150;
        double fig_Zmin = tot_Zmin - 150;
        double fig_Zmax = tot_Zmax + 150;
        epswr_figure_t *epsf = boap_new_figure(pname, "both", fig_Ymin, fig_Ymax, fig_Zmin, fig_Zmax, FALSE);
        boap_model_draw(epsf, skel_mod, 0, NAN);
        boap_model_draw(epsf, pans_mod, 0, NAN);
        epswr_end_figure(epsf);
      }
    else
      { /* Draws the two models as separate figures: */
        assert(dims);
        { double fig_Ymin = tot_Ymin - 250;
          double fig_Ymax = tot_Ymax + 150;
          double fig_Zmin = tot_Zmin - 250;
          double fig_Zmax = tot_Zmax + 450;
          epswr_figure_t *epsf = boap_new_figure(pname, "skel", fig_Ymin, fig_Ymax, fig_Zmin, fig_Zmax, FALSE);
          boap_model_draw(epsf, skel_mod, 0, NAN);
          epswr_end_figure(epsf);
        }
        { double fig_Ymin = tot_Ymin - 350;
          double fig_Ymax = tot_Ymax + 850;
          double fig_Zmin = tot_Zmin - 350;
          double fig_Zmax = tot_Zmax + 750;
          epswr_figure_t *epsf = boap_new_figure(pname, "pans", fig_Ymin, fig_Ymax, fig_Zmin, fig_Zmax, FALSE);
          boap_model_draw(epsf, pans_mod, 0, NAN);
          epswr_end_figure(epsf); epsf = NULL;
        }
      }
  }
  
boap_model_t *boap_make_front_panels
  ( boap_iron_dims_t *ird, 
    double pans_Xmin, double pans_Xmax,
  
    int8_t lpan_Vbar_N, double lpan_Ymin, double lpan_Ymax,
    int8_t door_Vbar_N, double door_Ymin, double door_Ymax,
    int8_t rpan_Vbar_N, double rpan_Ymin, double rpan_Ymax,
    double glass_Ysize,
    
    int8_t bpan_Hbar_N, double bpan_Zmin, double bpan_Zmax,
    int8_t tpan_Hbar_N, double tpan_Zmin, double tpan_Zmax,
    double glass_Zsize, double hole_Zsize,
    bool_t dims
  )
  {
    double Yplode = (dims ? 350 : 0);  /* {Y} explosion amount */
    double Zplode = (dims ? 300 : 0);  /* {Z} explosion amount */

    boap_model_t *mod = boap_model_new();
    
    double pan_Xmin = pans_Xmin, pan_Xmax = pans_Xmax;
    int8_t nxp = 3; /* Num of columns of panels. */
    int8_t nyp = 2; /* Numb of rows of panels. */
    for (int8_t ixp = 0; ixp < nxp; ixp++)
      { 
        /* Panel/door {Y} coord range: */
        double pan_Ymin, pan_Ymax;
        if (ixp == 0)
          { pan_Ymin = lpan_Ymin; pan_Ymax = lpan_Ymax; }
        else if (ixp == 1)
          { pan_Ymin = door_Ymin; pan_Ymax = door_Ymax; }
        else if (ixp == 2)
          { pan_Ymin = rpan_Ymin; pan_Ymax = rpan_Ymax; }
          
        /* Apply explosion shifts: */
        pan_Ymin += ixp*Yplode;
        pan_Ymax += ixp*Yplode;

        for (int8_t iyp = 0; iyp < nyp; iyp++)
          { 
            /* Panel {Z} coord range:*/
            double pan_Zmin, pan_Zmax;
            if (iyp == 0)
              { pan_Zmin = bpan_Zmin; pan_Zmax = bpan_Zmax; }
            else if (iyp == 1)
              { pan_Zmin = tpan_Zmin; pan_Zmax = tpan_Zmax; }

            /* Apply explosion shifts: */
            pan_Zmin += iyp*Zplode;
            pan_Zmax += iyp*Zplode;

            double Vbar_Ystep = glass_Ysize + ird->T_th1; 
            double Hbar_Zstep = (iyp == 0 ? glass_Zsize + ird->T_th1 : hole_Zsize + ird->Msma_wd);


            /* Number of vert and horz bars in panel's grille: */
            int8_t NV = (ixp == 0 ? lpan_Vbar_N : (ixp == 1 ? door_Vbar_N : rpan_Vbar_N)); 
            int8_t NH = (iyp == 0 ? bpan_Hbar_N : tpan_Hbar_N) ; /* Number of horz bars in even columns. */

            boap_add_front_panel
              ( mod, ixp, iyp, ird,
                pan_Xmin, pan_Xmax,
                pan_Ymin, pan_Ymax,
                pan_Zmin, pan_Zmax,
                NV, Vbar_Ystep, 
                NH, Hbar_Zstep, 
                dims
              );
            
            if (dims)
              { if (ixp == 0)
                  { /* Show panel {Z} dimensions: */
                    double xDim0 = pan_Ymin, yDim0 = pan_Zmin;
                    double xDim1 = xDim0,     yDim1 = pan_Zmax;
                    boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1,  +1.0,  +1.0, +18.0, TRUE);
                  }
                if (iyp == 0)
                  { /* Show panel {Y} dimensions: */
                    double xDim0 = pan_Ymin, yDim0 = pan_Zmin;
                    double xDim1 = pan_Ymax, yDim1 = yDim0;
                    boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1,  -1.0,  -1.0, -4.0, TRUE);
                  }
              }
          }
      }

    return mod;
  }

boap_model_t *boap_make_front_skeleton
  ( boap_iron_dims_t *ird, 
    double p_Xmin, double p_Xmax, 
    double p_Ymin, double p_Ymax,
    double p_Zmin, double p_Zmax,
    double s_Ymin, double s_Ymax,
    double s_Zmax,
    bool_t dims
  )
  {
    frgb_t *L_color = &(ird->L_color);
    
    boap_model_t *mod = boap_model_new();
    
    /* Top 'L' bar: */
    double tp_Xctr = p_Xmin;
    double tp_Yctr = (p_Ymin + p_Ymax)/2;
    double tp_Ysize = p_Ymax - p_Ymin;
    double tp_Zctr = p_Zmax;
    boap_add_bar
      ( mod, 'L',
        tp_Xctr, tp_Yctr, tp_Zctr,
        1, tp_Ysize,  
        0, +ird->L_wd, +ird->L_th,
        2, -ird->L_wd, -ird->L_th,
        L_color
      );
    
    /* Lintel 'L' bar: */
    double sb_Xctr = p_Xmin;
    double sb_Yctr = (p_Ymin + p_Ymax)/2;
    double sb_Ysize = p_Ymax - p_Ymin;
    double sb_Zctr = s_Zmax;
    boap_add_bar
      ( mod, 'L',
        sb_Xctr, sb_Yctr, sb_Zctr,
        1, sb_Ysize,  
        0, +ird->L_wd, +ird->L_th,
        2, -ird->L_wd, -ird->L_th,
        L_color
      );
    
    /* Outer left post: */
    double post0_Xctr = p_Xmin;
    double post0_Yctr = p_Ymin;
    double post0_Zctr = (p_Zmin + p_Zmax)/2;
    double post0_Zsize = p_Zmax - p_Zmin;
    boap_add_bar
      ( mod, 'L',
        post0_Xctr, post0_Yctr, post0_Zctr,
        2, post0_Zsize, 
        0, +ird->L_wd, +ird->L_th,
        1, +ird->L_wd, +ird->L_th,
        L_color
      );

    /* Left jamb: */
    double jamb0_Xctr = p_Xmin;
    double jamb0_Yctr = s_Ymin;
    double jamb0_Zctr = (p_Zmin + p_Zmax)/2;
    double jamb0_Zsize = p_Zmax - p_Zmin;
    boap_add_bar
      ( mod, 'L',
        jamb0_Xctr, jamb0_Yctr, jamb0_Zctr,
        2, jamb0_Zsize, 
        0, +ird->L_wd, +ird->L_th,
        1, +ird->L_wd, +ird->L_th,
        L_color
      );

    /* Right jamb: */
    double jamb1_Xctr = p_Xmin ;
    double jamb1_Yctr = s_Ymax;
    double jamb1_Zctr = (p_Zmin +p_Zmax)/2;
    double jamb1_Zsize = p_Zmax - p_Zmin;
    boap_add_bar
      ( mod, 'L',
        jamb1_Xctr, jamb1_Yctr, jamb1_Zctr,
        2, jamb1_Zsize, 
        0, +ird->L_wd, +ird->L_th,
        1, -ird->L_wd, -ird->L_th,
        L_color
      );

    /* Right post: */
    double post1_Xctr = p_Xmin ;
    double post1_Yctr = p_Ymax;
    double post1_Zctr = (p_Zmin +p_Zmax)/2;
    double post1_Zsize = p_Zmax - p_Zmin;
    boap_add_bar
      ( mod, 'L',
        post1_Xctr, post1_Yctr, post1_Zctr,
        2, post1_Zsize, 
        0, +ird->L_wd, +ird->L_th,
        1, -ird->L_wd, -ird->L_th,
        L_color
      );

    if (dims)
      { /* Show wall and jamb dimensions: */
        { /* Total width: */
          double xDim0 = p_Ymin, yDim0 = p_Zmax;
          double xDim1 = p_Ymax, yDim1 = yDim0;
          boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1, +1.0, +1.0, +12.0, TRUE);
        }      
        { /* Total height: */
          double xDim0 = p_Ymin, yDim0 = p_Zmin;
          double xDim1 = xDim0,     yDim1 = p_Zmax;
          boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1,  +1.0, +1.0, +12.0, TRUE);
        }
        { /* Distance from left edge of wall and left door jamb: */
          double xDim0 = p_Ymin,  yDim0 = p_Zmax;
          double xDim1 = s_Ymin,  yDim1 = yDim0;
          boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1, +1.0, +1.0, +6.0, TRUE);
        }
        { /* Door width including jambs: */
          double xDim0 = s_Ymin, yDim0 = p_Zmax;
          double xDim1 = s_Ymax, yDim1 = yDim0;
          boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1, +1.0, +1.0, +6.0, TRUE);
        }
        { /* Distance from right jamb to right edge of wall: */
          double xDim0 = s_Ymax, yDim0 = p_Zmax;
          double xDim1 = p_Ymax, yDim1 = yDim0;
          boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1, +1.0, +1.0, +6.0, TRUE);
        }
        { /* Lintel height: */
          double xDim0 = p_Ymin, yDim0 = p_Zmin;
          double xDim1 = xDim0,     yDim1 = s_Zmax;
          boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1,  +1.0,  +1.0, +6.0, TRUE);
        }
        { /* Height above lintel: */
          double xDim0 = p_Ymin, yDim0 = s_Zmax;
          double xDim1 = xDim0,     yDim1 = p_Zmax;
          boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1,  +1.0,  +1.0, +6.0, TRUE);
        }
      }    
    return mod;
  }

void boap_add_front_panel
  ( boap_model_t *mod, 
    int8_t ixp, 
    int8_t iyp,
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    int8_t NV, double Vbar_Ystep,
    int8_t NH, double Hbar_Zstep,
    bool_t dims
  )
  {
    frgb_t *Mbig_color = &(ird->Mbig_color);
    
    double Pthk = 2.0;   /* Plate thickness includin pleats. */
    
    fprintf(stderr, "--- panel[%d,%d] --------------------------\n", ixp,iyp);

    boap_add_rect_tube_frame
      ( mod, ird,
        Xmin, Xmax,
        Ymin, Ymax,
        Zmin, Zmax,
        dims,
        Mbig_color
      );
    
    if (iyp == 0)
      { /* Panel has two sections: */
        double Zdisp = 960; /* Distance from bot of bottom panel to top of transverse bar. */
        
        /* Draw middle metalon: */
        double md_Zctr = Zmin + Zdisp - 0.5*ird->Mbig_wd;
        double md_Xctr = (Xmin + Xmax)/2;
        double md_Yctr = (Ymin + Ymax)/2;
        double md_Ysize = (Ymax - Ymin) - 2*ird->Mbig_wd;
        boap_add_rect_tube
          ( mod, 1,
            md_Xctr, md_Yctr, md_Zctr, 
            ird->Mbig_wd, md_Ysize, ird->Mbig_wd,
            Mbig_color
          );
        if (dims)
          { /* Show midbar {Z} dimensions: */
            double xDim0 = Ymin;
            double xDim1 = xDim0;
            double yDim0 = md_Zctr - ird->Mbig_wd/2;
            double yDim1 = md_Zctr + ird->Mbig_wd/2;
            boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1,  +1.0,  +1.0, +13.0, FALSE);
          }

        /* Draw skirt. */
        double sk_Xmin = Xmin + (Xmax - Xmin - Pthk)/2, sk_Xmax = sk_Xmin + Pthk;
        double sk_Ymin = Ymin + ird->Mbig_wd, sk_Ymax = Ymax - ird->Mbig_wd;
        double sk_Zmin = Zmin + ird->Mbig_wd, sk_Zmax = Zmin + Zdisp - ird->Mbig_wd;
        boap_add_front_panel_skirt
          ( mod, ixp, ird,
            sk_Xmin, sk_Xmax,
            sk_Ymin, sk_Ymax,
            sk_Zmin, sk_Zmax,
            dims
          );

        /* Draw the vitral grille: */
        double vt_Xmin = Xmin, vt_Xmax = vt_Xmin + ird->T_wd1;
        double vt_Ymin = Ymin + ird->Mbig_wd, vt_Ymax = Ymax - ird->Mbig_wd;
        double vt_Zmin = Zmin + Zdisp, vt_Zmax = Zmax - ird->Mbig_wd;
        boap_add_front_panel_vitral
          ( mod, ixp, ird,
            vt_Xmin, vt_Xmax, 
            vt_Ymin, vt_Ymax, 
            vt_Zmin, vt_Zmax,
            NV, Vbar_Ystep, 
            NH, Hbar_Zstep,
            dims
          );
      }
    else
      { /* Panel is an open grille:*/
        double og_Xmin = Xmin + (Xmax - Xmin - ird->Msma_wd)/2, og_Xmax = og_Xmin + ird->Msma_wd;
        double og_Ymin = Ymin + ird->Mbig_wd, og_Ymax = Ymax - ird->Mbig_wd;
        double og_Zmin = Zmin + ird->Mbig_wd, og_Zmax = Zmax - ird->Mbig_wd;
        boap_add_front_panel_open_grille
          ( mod, ixp, ird,
            og_Xmin, og_Xmax, 
            og_Ymin, og_Ymax, 
            og_Zmin, og_Zmax,
            NV, Vbar_Ystep, 
            NH, Hbar_Zstep,
            dims
          );
      }
    fprintf(stderr, "-----------------------------------------\n");
  }

void boap_add_front_panel_skirt
  ( boap_model_t *mod, 
    int8_t ixp,
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    bool_t dims
  )
  {
    frgb_t *P_color = boap_make_color(0.200, 0.700, 0.000);
    
    double P_th = 1.0; /* Thickness of skirt. */
    
    /* Number of half-pleats rounded up to odd, and half-pleat step: */ 
    double Ystep_ideal = 91.0; /* Ideal width of 1 half period of the pleated "skirt". */
    int8_t NP = (int8_t)floor((Ymax - Ymin - P_th)/Ystep_ideal + 0.10);
    if ((NP % 2) == 0) { NP++; }
    double Ystep = (Ymax - Ymin - P_th)/NP;
    
    double Xctr = Xmin, Xsize = (Xmax - Xmin)/2;       
    double Yctr = (Ymin + Ymax)/2, Ysize = Ystep + P_th;
    double Zctr = (Zmin + Zmax)/2, Zsize = (Zmax - Zmin);

    /* Loop until filling space: */
    double Xdir = -1; /* {X} sense of concavity of next half-pleat(s). */
    int8_t kP = 0;
    while (TRUE)
      { double Yctr1 = Yctr + kP*Ystep;
        double Yctr0 = Yctr - kP*Ystep;
        if (Yctr1 > Ymax) { break; }
        
        boap_add_bar
          ( mod, 'U',
            Xctr, Yctr1, Zctr,
            2, Zsize, 
            1, Ysize, P_th, 
            0, Xsize, P_th,
            P_color
          );
        if (kP != 0)
          {  boap_add_bar
              ( mod, 'U',
                Xctr, Yctr0, Zctr,
                2, Zsize, 
                1, Ysize, P_th, 
                0, Xsize, P_th,
                P_color
              );
          }
        /* Advance to next pair of pleats: */
        Xdir = -Xdir;
        Xctr = Xmin + (Xmax - Xmin) - (Xctr - Xmin);
        kP++;
    }
    if (dims)
      { /* Show dimensions: */
        { /* Show skirt {Z} dimensions: */
          double xDim0 = Ymin;
          double xDim1 = xDim0;
          double yDim0 = Zmin;
          double yDim1 = Zmax;
          boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1,  +1.0,  +1.0, +12.0, TRUE);
        }
      }
  }

void boap_add_front_panel_vitral
  ( boap_model_t *mod, 
    int8_t ixp,
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    int8_t NV, double Vbar_Ystep, 
    int8_t NH, double Hbar_Zstep,
    bool_t dims
  )
  {
    assert(fabs((Xmax - Xmin) - ird->T_wd0) < 0.001);
    boap_add_grille
      ( mod, 'T', ird,
        Xmin, Xmax,  Ymin, Ymax, Zmin, Zmax,
        NV, Vbar_Ystep,
        NH, Hbar_Zstep,
        dims
      );
  }

void boap_add_front_panel_open_grille 
  ( boap_model_t *mod, 
    int8_t ixp,
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    int8_t NV, double Vbar_Ystep, 
    int8_t NH, double Hbar_Zstep,
    bool_t dims
  )
  {
    assert(fabs((Xmax - Xmin) - ird->Msma_wd) < 0.001);
    boap_add_grille
      ( mod, 'M', ird,
        Xmin, Xmax,  Ymin, Ymax, Zmin, Zmax,
        NV, Vbar_Ystep, 
        NH, Hbar_Zstep,
        dims
      );
  }

void boap_add_grille
  ( boap_model_t *mod, 
    char type,
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    int8_t NV, double Vbar_Ystep, 
    int8_t NH, double Hbar_Zstep,
    bool_t dims
  )
  { 
    /* {X} center and extent of all tubes: */
    double bar_Xctr;  /* {X} of all bar centers. */
    double bar_Xsize; /* {X} dimension of 'T'-irons or metalons. */
    double Vbar_wd;   /* {Y} dimension of 'T'-irons or metalons in vert bars. */
    double Hbar_wd;   /* {Z} dimension of 'T'-irons or metalons in horz bars. */
    double Vbar_nwd, Hbar_nwd; /* Nominal bar widths for gap size computation. */
    frgb_t *color;

    if (type == 'T')
      { /* Grille is made of 'T' irons, to support glass panes. */
        bar_Xctr = Xmin; /* The 'T'-iron axis is on the top edge of the 'T'. */
        bar_Xsize = Xmax - Xmin; /* {X} dimension of bars. */
        Vbar_wd = ird->T_wd0; /* {Y} dimension of bars. */
        Hbar_wd = ird->T_wd0; /* {Z} dimension of bars. */

        Vbar_nwd = ird->T_th1; /* Nominal bar width for gap size compt. */
        Hbar_nwd = ird->T_th1; /* Nominal bar width for gap size compt. */
        
        boap_add_flat_frame
          ( mod, ird,
            Xmin, Xmin + ird->F_th, 
            Ymin, Ymax,
            Zmin, Zmax,
            dims
          );
          
        color = &(ird->T_color);
      }
    else if (type == 'M')
      { /* Grille is open, made of metalons. */
        bar_Xsize = ird->Msma_wd; /* {X} dimension of bars. */
        bar_Xctr = (Xmax + Xmin)/2; /* The metalon axis is in the center. */
        Vbar_wd = ird->Msma_wd; /* {Y} dimension of bars. */
        Hbar_wd = ird->Msma_wd; /* {Z} dimension of bars. */

        Vbar_nwd = ird->Msma_wd; /* Nominal bar width for gap size compt. */
        Hbar_nwd = ird->Msma_wd; /* Nominal bar width for gap size compt. */

        color = &(ird->Msma_color);
      }
    else
      { assert(FALSE); }


    /* {Z} center and {Z} extent of vertical grille beams: */
    double Vbar_Zsize = Zmax - Zmin;       /* {Z} extent of vert beams. */
    double Vbar_Zctr = (Zmin + Zmax)/2;    /* {Z} coord of center of vert beams. */

   /* {Y} coord of first vert beam, or of right frame if none: */
    double Vbar_DY0; 
    if (NV > 0)
      { /* Split space equally among first and last gaps: */
        Vbar_DY0 = ((Ymax - Ymin) - (NV - 1)*Vbar_Ystep)/2;
      }
    else
      { Vbar_DY0 = Ymax - Ymin; }
    fprintf(stderr, "  %2d vert beams at Ymin + %7.1f + k*%7.1f mm\n", NV, Vbar_DY0, Vbar_Ystep);

     /* {Z} coord of first vert beam in even columns, or of top frame if none: */
    double Hbar_DZ0; 
    if (NH > 0)
      { /* Split space equally among first and last gaps: */
        Hbar_DZ0 = ((Zmax - Zmin) - (NH - 1)*Hbar_Zstep)/2;
      }
    else
      { Hbar_DZ0 = Zmax - Zmin; }

    double Vini_Yctr = Ymin - Vbar_nwd/2;   /* Center {Y} of imaginary vert bar at left. */
    double Vfin_Yctr = Ymax + Vbar_nwd/2;   /* Center {Y} of imaginary vert bar at right. */
    double Vmin_Yctr = Ymin + Vbar_wd/2; /* Min center {Y} of real vert bar. */
    double Vmax_Yctr = Ymax - Vbar_wd/2; /* Max center {Y} of real vert bar. */

    double Vprev_Yctr = Vini_Yctr;
    double Vnext_Yctr = Ymin + Vbar_DY0;
    double Hbar_Zshift = 0.0; /*  {Z} shift for horz bars between {Vprev_Yctr} and {Vnext_Yctr}. */
    for (int8_t kV = 0; kV <= NV; kV++)
      { 
        /* Correction for last gap: */
        if (kV >= NV) { Vnext_Yctr = Vfin_Yctr; }

        if ((Vprev_Yctr >= Vmin_Yctr) && (Vprev_Yctr <= Vmax_Yctr))
          { /* Draw vert beam at {Vprev_Yctr}: */
            if (type == 'M')
              { boap_add_rect_tube
                  ( mod, 2, 
                    bar_Xctr, Vprev_Yctr, Vbar_Zctr, 
                    bar_Xsize, Vbar_wd, Vbar_Zsize,
                    color
                  );
              }
            else
              { boap_bar_t *bar = boap_bar_new
                  ( 'T', 
                    Xmin, Vprev_Yctr, Vbar_Zctr, 
                    1, Vbar_wd,   ird->T_th1,
                    0, bar_Xsize, ird->T_th0,
                    2, Vbar_Zsize, 
                    color
                  );
                boap_bar_trim(bar, 0, -1, ird->F_wd);
                boap_bar_trim(bar, 0, +1, ird->F_wd);
                boap_model_add_bar(mod,bar);
              } 
          }

        /* Draw horz beams between {Vprev_Yctr} and {Vnext_Yctr}: */

        /* {Y} coord ranges of the horz bars to the left of bar {kV}: */
        double Hbar_Ymin = Vprev_Yctr + Vbar_nwd/2;
        double Hbar_Ymax = Vnext_Yctr - Vbar_nwd/2;
        double Hbar_Yctr = (Hbar_Ymin + Hbar_Ymax)/2;  /* {Y} coord of center of horz beams. */
        double Hbar_Ysize = Hbar_Ymax - Hbar_Ymin;

        /* Loop on horizontal bars: */
        double Hini_Zctr = Zmin - Hbar_nwd/2;    /* Center {Z} of imaginary horz bar at bottom. */
        double Hfin_Zctr = Zmax + Hbar_nwd/2;    /* Center {Z} of imaginary horz bar at top. */
        double Hmin_Zctr = Zmin + Hbar_wd/2;  /* Min center {Z} of real horz bar. */
        double Hmax_Zctr = Zmax - Hbar_wd/2;  /* Max center {Z} of real horz bar. */

        double Hprev_Zctr = Hini_Zctr;
        double Hnext_Zctr = Zmin + Hbar_DZ0;
        double H_Ytrim = (ird->T_wd0 - ird->T_th1)/2; /* How much to trim from typical horz T-bars */
        /* Adjust start by even-odd shift: */
        Hnext_Zctr = Hnext_Zctr - Hbar_Zshift;
        if (Hnext_Zctr <= Hprev_Zctr + Hbar_wd) { Hnext_Zctr += Hbar_Zstep; }

        for (int8_t kH = 0; kH <= NH+1; kH++)
          { 
            /* Correction for last gap: */
            if (Hnext_Zctr > Hmax_Zctr + 0.001) { Hnext_Zctr = Hfin_Zctr; }

            if ((Hprev_Zctr >= Hmin_Zctr) && (Hprev_Zctr <= Hmax_Zctr))
              { /* Draw bar at {Hprev_Zctr}: */
                if (type == 'M')
                  { boap_add_rect_tube
                      ( mod, 1, 
                        bar_Xctr, Hbar_Yctr, Hprev_Zctr, 
                        bar_Xsize, Hbar_Ysize, Hbar_wd,
                        color
                      );
                  }
                else
                  { boap_bar_t *bar = boap_bar_new
                      ( 'T', 
                        Xmin, Hbar_Yctr, Hprev_Zctr, 
                        2, Vbar_wd,   ird->T_th1,
                        0, bar_Xsize, ird->T_th0,
                        1, Hbar_Ysize, 
                        color
                      );
                    boap_bar_trim(bar, 0, -1, (Vprev_Yctr < Ymin ? ird->F_wd : H_Ytrim));
                    boap_bar_trim(bar, 0, +1, (Vnext_Yctr > Ymax ? ird->F_wd : H_Ytrim));
                    boap_model_add_bar(mod,bar);
                  } 
              }
              
            if (dims && (Hprev_Zctr < Zmax) && ((kV == 0) || (kV == NV-1)))
              { /* Show {Z} spacing between horz bars: */
                double xDim0 = (kV == 0 ? Ymin : Vnext_Yctr);
                double xDim1 = xDim0;
                double yDim0 = (Hprev_Zctr < Zmin ? Zmin : Hprev_Zctr);
                double yDim1 = (Hnext_Zctr > Zmax ? Zmax : Hnext_Zctr);
                double abgap = (kV == 0 ? +1.0 : -1.5);
                double elen = (kV == 0 ? +6.0 : -14.0);
                boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1, abgap, abgap, elen, TRUE);
              }
             
             /* Advance to next gap: */
             Hprev_Zctr = Hnext_Zctr;
             Hnext_Zctr += Hbar_Zstep;
          }

        if (dims)
          { /* Show width of gaps between vert beams: */
            double xDim0 = (Vprev_Yctr < Ymin ? Ymin : Vprev_Yctr);
            double xDim1 = (Vnext_Yctr > Ymax ? Ymax : Vnext_Yctr);
            double yDim0 = Zmax;
            double yDim1 = yDim0;
            boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1, +3.0, +3.0, +6.0, TRUE);
          }

        /* Advance to next column: */
        Vprev_Yctr = Vnext_Yctr;
        Vnext_Yctr += Vbar_Ystep;
        Hbar_Zshift = (Hbar_Zshift == 0 ? 0.5*Hbar_Zstep : 0);
      }
    if (dims)
      { 
        { /* Show grille {Y} dimensions: */
          double xDim0 = Ymin;
          double xDim1 = Ymax;
          double yDim0 = Zmax;
          double yDim1 = yDim0;
          boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1,  +1.0, +1.0, +12.0, TRUE);
        }
        { /* Show grille {Z} dimensions: */
          double xDim0 = Ymin;
          double xDim1 = xDim0;
          double yDim0 = Zmin;
          double yDim1 = Zmax;
          boap_add_dim(mod, xDim0, yDim0, xDim1, yDim1,  +1.0, +1.0, +12.0, TRUE);
        }
      }
  }

void boap_add_flat_frame
  ( boap_model_t *mod, 
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    bool_t dims
  )
  {
    frgb_t *F_color = &(ird->F_color);
    
    double bar_Xctr, bar_Yctr, bar_Zctr;     /* Coords of center of a bar. */
    double bar_Xsize, bar_Ysize, bar_Zsize;  /* Dimensions of a bar. */
    
    bar_Xctr = Xmin + ird->F_th/2;
    bar_Xsize = ird->F_th;
    assert(fabs(bar_Xsize - (Xmax -Xmin)) < 0.00001);
    
    /* Draw bars parallel to {Y} axis: */
    bar_Yctr = (Ymin + Ymax)/2;
    bar_Ysize = Ymax - Ymin;
    bar_Zsize = ird->F_wd;
    for (int8_t side = 0; side <= 1; side++) 
      { /* Draw bar on low ({side=0})or high ({side=1}) side: */
        bar_Zctr = (side == 0 ? Zmin + ird->F_wd/2 : Zmax - ird->F_wd/2);
        boap_add_bar
          ( mod, 'F',
            bar_Xctr,  bar_Yctr,  bar_Zctr, 
            1, bar_Ysize, 
            2, bar_Zsize, NAN,
            0,   NAN, bar_Xsize,
            F_color
          );
      }

    /* Draw bars parallel to {Z} axis: */
    bar_Zctr = (Zmin + Zmax)/2;
    bar_Zsize = Zmax - Zmin - 2*ird->F_wd;
    bar_Ysize = ird->F_wd;
    for (int8_t side = 0; side <= 1; side++) 
      { /* Draw bar on low ({side=0})or high ({side=1}) side: */
        bar_Yctr = (side == 0 ? Ymin + ird->F_wd/2 : Ymax - ird->F_wd/2);
        boap_add_bar
          ( mod, 'F',
            bar_Xctr,  bar_Yctr,  bar_Zctr, 
            2, bar_Zsize, 
            1, bar_Ysize, NAN,
            0,   NAN, bar_Xsize,
            F_color
          );
      }

  }

void boap_add_rect_tube_frame
  ( boap_model_t *mod, 
    boap_iron_dims_t *ird, 
    double Xmin, double Xmax,
    double Ymin, double Ymax,
    double Zmin, double Zmax,
    bool_t dims,
    frgb_t *color
  )
  {
    /* Sort interval bounds: */
    if (Xmin > Xmax) { double tmp = Xmin; Xmin = Xmax; Xmax = tmp; }
    if (Ymin > Ymax) { double tmp = Ymin; Ymin = Ymax; Ymax = tmp; }
    if (Zmin > Zmax) { double tmp = Zmin; Zmin = Zmax; Zmax = tmp; }
    
    demand(fabs((Xmax - Xmin) - ird->Mbig_wd) < 0.0001, "inconsistent metalon {X} width {ird->Mbig_wd}");
    
    /* {X} center and extent of all tubes: */
    double Xctr = Xmin + ird->Mbig_wd/2; /* {X} of all tubes. */

    /* Parameters of horz frame tubes: */
    double H_Yctr = (Ymin + Ymax)/2;      /* {Y} of tubes. */
    double H_Ysize = Ymax - Ymin;         /* {Y} extent of tubes. */
    double H0_Zctr = Zmin + ird->Mbig_wd/2; /* {Z} of bot tube. */
    double H1_Zctr = Zmax - ird->Mbig_wd/2; /* {Z} of top tube. */

    boap_add_rect_tube(mod, 1, Xctr, H_Yctr, H0_Zctr, ird->Mbig_wd, H_Ysize, ird->Mbig_wd, color);
    boap_add_rect_tube(mod, 1, Xctr, H_Yctr, H1_Zctr, ird->Mbig_wd, H_Ysize, ird->Mbig_wd, color);

    /* Parameters of vert frame tubes: */
    double V0_Yctr = Ymin + ird->Mbig_wd/2; /* {Y} of left tube. */
    double V1_Yctr = Ymax - ird->Mbig_wd/2; /* {Y} of right tube. */
    double V_Zctr = (Zmin + Zmax)/2;      /* {Z} of tubes. */
    double V_Zsize = Zmax - Zmin - 2*ird->Mbig_wd;  /* {Z} extent of tubes. */
    
    boap_add_rect_tube(mod, 2, Xctr, V0_Yctr, V_Zctr, ird->Mbig_wd, ird->Mbig_wd, V_Zsize, color);
    boap_add_rect_tube(mod, 2, Xctr, V1_Yctr, V_Zctr, ird->Mbig_wd, ird->Mbig_wd, V_Zsize, color);
  }

void boap_add_dim
  ( boap_model_t *mod,
    double xa, double ya,
    double xb, double yb,
    double agap, 
    double bgap, 
    double elen,
    bool_t inner
  )
  { 
    boap_dim_t *dim = boap_dim_new(0, xa, ya, 0, xb, yb, agap, bgap, elen, inner, 0);
    boap_model_add_dim(mod, dim);
  }

void boap_add_rect_tube
  ( boap_model_t *mod, 
    int8_t axC,
    double Xctr, double Yctr, double Zctr, 
    double Xsize, double Ysize, double Zsize,
    frgb_t *color
  )
  {
    double M_ABth = 1.0; /* Thickness of walls. */
    
    double size[3] = { Xsize, Ysize, Zsize }; /* Extents of bar in each world coordinate. */
    int8_t axA = (int8_t)((axC + 1) % 3);
    int8_t axB = (int8_t)((axA + 1) % 3);
    boap_bar_t *bar = boap_bar_new
      ( 'M',
        Xctr, Yctr, Zctr,
        axA, size[axA], M_ABth,
        axB, size[axB], M_ABth,
        axC, size[axC],
        color
      );
    boap_model_add_bar(mod, bar);
  }
  
void boap_add_bar
  ( boap_model_t *mod, 
    char type,
    double Xctr, double Yctr, double Zctr, 
    int8_t axC, double Csize, 
    int8_t axA, double Awd, double Athk,
    int8_t axB, double Bwd, double Bthk,
    frgb_t *color
  )
  {
    boap_bar_t *bar = boap_bar_new
      ( type,
        Xctr, Yctr, Zctr,
        axA, Awd, Athk,
        axB, Bwd, Bthk,
        axC, Csize,
        color
      );
    boap_model_add_bar(mod, bar);
  }

frgb_t *boap_make_color(double R, double G, double B)
  {
    frgb_t *color = notnull(malloc(sizeof(frgb_t)), "no mem");
    (*color) = (frgb_t){{ (float)R, (float)G, (float)B }};
    return color;
  }
