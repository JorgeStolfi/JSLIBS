/* Tech drawing for the new center gate of Boaretto da Silva 113 */
/* Last edited on 2024-12-05 10:15:23 by stolfi */

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include <r2.h>
#include <r3.h>
#include <bool.h>
#include <sign.h>
#include <affirm.h>
#include <jsfile.h>

#include <epswr.h>
#include <epswr_dim.h>

/* All dimensions are in mm. The World coordinate axes are: {X}
  horizontal towards the street, {Y} horizontal along the street (left
  to right seen from {+X}), {Z} up. The Plot coordinate axes are: {x}
  horizontal, left to right, {y} vertical, up. */
  
#define D (1.0)
  /* Scaling factor for dimension marking placement. */

void boag_draw_gate(char *fname, bool_t cut);
  /* Writes an EPS file called {fname} with a drawing of the gate.
    If {cut} is {FALSE}, draws the front view, from the {+X} direction. 
    Otherwise draws an horizontal section.
    
    The origin is at the left bottom corner of the gate, at the
    back of the main tube frame. Note that pillars and 
    some details may extend into negative coordinates.  */

void boag_draw_leaf_and_post
  ( epswr_figure_t *epsf, 
    bool_t cut,
    sign_t which, 
    double Xpos, double Ypos, double Zpos, 
    double Xsize, double Ysize, double Zsize,
    double Ypgap,
    double XYpsize,
    double Zpdepth,
    double *YminP, 
    double *YmaxP
  );
  /* Draws on {epsf} one leaf of the gate, including its supporting
    post, hinges, latches, etc.
    
    For the left leaf, the {which} parameter should be {0} and {Ypos}
    should be the {Y} coordinate of the left edge of the grille. For the
    right leaf, {which} should be {1}, and {Ypos} should be the {Y}
    coordinate of the right edge of the grille.
    
    In either case, {Xpos} is the {X} coordinate of the back plane of
    the grille, {Zpos} is the {Z} coordinate of the bottom edge of the
    grille and {Xsize}, {Ysize}, and {Zsize} should be the dimensions of
    the grille.  
    
    Parameter {Ypgap} is the {Y} distance between the leaf and its
    supporting post. {XYpsize} are the {X} and {Y} dimensions of the
    post. Parameter {Zpdepth} is how much the post extends below the
    bottom edge of the gate leaf.
    
    Returns in {*YminP,*YmaxP} the total range of {Y} coordinates,
    INCLUDING the side posts. */

void boag_draw_leaf
  ( epswr_figure_t *epsf, 
    bool_t cut,
    sign_t which, 
    double Xpos, double Ypos, double Zpos, 
    double Xsize, double Ysize, double Zsize,
    double Zhinge0, double Zhinge1,
    double Ypgap
  );
  /* Draws on {epsf} one leaf of the gate, including latches and 
    other trimmings, WITHOUT its supporting post and the post-side half 
    of the hinges.
    
    Parameters {Zhinge0} and {Zhinge1} are the {Z} coordinates of the
    middle of the hinges. The other parameters are
    the same as those of {boag_draw_leaf_and_post}. */
  
void boag_draw_post
  ( epswr_figure_t *epsf, 
    bool_t cut,
    sign_t which, 
    double Xpos, double Ypos, double Zpos,
    double Xsize, double Zsize,
    double Zhinge0, double Zhinge1,
    double Ypgap,
    double XYpsize,
    double Zpdepth
  );
  /* Draws on {epsf} the supporting post and hinges of one 
    leaf of the gate. 
    
    The other parameters are the same as 
    those of {boag_draw_leaf_and_post}. */

void boag_draw_leaf_frame
  ( epswr_figure_t *epsf, 
    bool_t cut,
    double Xmin, double Xmax,
    double Ymin, double Ymax, double tube_Ysize,
    double Zmin, double Zmax, double tube_Zsize
  );
  /* Draws into {epsf} the frame of a leaf that fits 
    snugly inside the box {[Xmin _ Xmax]}  by {[Ymin _ Ymax]}
    by {[Zmin _ Zmax]}.  The bounds of each coordinate can be given 
    in any order.  
    
    Assumes that the tube cross sections are {|Xmax-Xmin|} in {X},
    {tube_Ysize} in {Y}, and {tube_Zsize} in {Z}. */

void boag_draw_leaf_beams
  ( epswr_figure_t *epsf, 
    bool_t cut,
    double Xmin, double Xmax,
    double Ymin, double Ymax, double tube_Ysize,
    double Zmin, double Zmax, double tube_Zsize
  );
  /* Draws into {epsf} the vertical and horizontal beams of
    a grille.  The beams will fit snugly inside the box 
    {[Xmin _ Xmax]}  by {[Ymin _ Ymax]} by {[Zmin _ Zmax]}.
    
    The bounds of each coordinate can be given in any order,
    except that the horizontal beams next to the original {Ymin}
    will be in the "even" position.
    
    Assumes that the tube cross sections are {|Xmax-Xmin|} in {X},
    {tube_Ysize} in {Y}, and {tube_Zsize} in {Z}. */

void boag_draw_half_hinge
  ( epswr_figure_t *epsf, 
    bool_t cut,
    int32_t which,
    double Xctr, double Yctr, double Zctr,
    double Yhdist
  );
  /* Draws one half of a gate hinge.  
  
    If {which} is 0, draws the bottom half of the hinge, attached to the
    post; if {which} is 1, draws the top half, attached to the leaf.
    Assumes that the center of the whole hinge's axis is
    {(Xctr,Yctr,Zctr)}. Assumes that the leaf's edge lies at
    {Yctr+Yhdist}, where {Yhdist} may be negative. */

void boag_draw_plate
  ( epswr_figure_t *epsf, 
    bool_t cut,
    double Xctr, double Yctr, double Zctr, 
    double Xsize, double Ysize, double Zsize
  );
  /* Draws an axis-aligned box wth center at {(Xctr,Yctr,Zctr)}
    and extents {Xsize {Ysize}, and {Zsize} along the three axes. */

void boag_draw_tube
  ( epswr_figure_t *epsf, 
    bool_t cut,
    int32_t axis,
    double Xctr, double Yctr, double Zctr, 
    double Xsize, double Ysize, double Zsize
  );
  /* Draws into {epsf} a tube with rectangular cross-section given its center
    and extent in each coordinate direction.  The tube's axis is
    assumed to be parallel to coordinate axis {axis} 
    ({0=X}, {1=Y}, {2=Z}). */

void boag_draw_vert_cylinder
  ( epswr_figure_t *epsf, 
    bool_t cut,
    double Xctr, double Yctr, 
    double Zmin, double Zmax,
    double Hrad
  );
  /* Draws a cylinder of radius {Hrad}, whose axis is parallel to the {Z} axis,
    with axis at {X=Xctr} and {Y=Yctr}, spanninn in {Z}
    from {Zmin} to {Zmax}. */

void boag_show_dim
  ( epswr_figure_t *epsf, 
    double xa, double ya,
    double xb, double yb,
    double agap,
    double bgap, 
    double elen,
    int32_t nfrac,
    bool_t inner
  );
  /* Shows the distance between points {a=(xa,ya)} and {b=(xb,yb)}
    as a label with extension and dimension lines.  
    
    The parameters {agap,alen} define the start and length (mm) of the
    extension line for the point {a}. The parameters {bgap,blen}
    have the same meaning for {b}. The dimension will be printed with {nfrac}
    decimals after point; if 0, point is omitted.  If {inner}, the dimension line
    is drawn between the extension lines, otherwise outside them. */

/* IMPLEMENTATIONS */

int32_t main (int32_t argc, char **argv)
  {
    boag_draw_gate("out/gate_plan.eps", FALSE);
    /* boag_draw_gate("out/gate_hcut.eps", TRUE); */

    return 0;
  }

void boag_draw_gate(char *fname, bool_t cut)
  {
    /* Device coordinate extent (pt; aspect ratio 13:10): */
    double hSize = 780.0;
    double vSize = 600.0;
    double hvMarg = 6.0;

    FILE *wr = open_write(fname, TRUE);
    bool_t verbose = FALSE;
    epswr_figure_t *epsf = epswr_new_figure(wr, hSize, vSize, hvMarg, hvMarg, hvMarg, hvMarg, verbose);
    epswr_set_fill_color(epsf, 1.000,0.800,0.400);
    epswr_set_label_font(epsf, "Courier-Bold", 10.0);
    epswr_set_pen(epsf, 0.000,0.000,0.000, 0.25, 0.0,0.0);
    
    /* Nominal dimensions of gate leaves: */
    double gate_Xsize = 25.0;    /* Thickness of leaf grille. */
    double gate0_Ysize = 870.0;  /* Width of left leaf. */
    double gate_Ygap = 5.0;      /* Gap between leaves. */
    double gate1_Ysize = 1353.0; /* Width of right leaf. */
    double gate_Ysize = gate0_Ysize + gate_Ygap + gate1_Ysize; /* Width of both leaves and gap. */
    double gate_Zsize = 2075.0;  /* Height of leaves. */
    
    double XYpsize = 80; /* {X} and {Y} dimensions of post. */
    double Ypgap = 37; /* {Y} gap between post and leaf. */
    double Zpdepth = 41; /* Extent of post below bottom edge of the leaf. */

    /* Nominal positions of gate leaves: */
    double gate_Xpos = 0.0;         /* Back of leaf grilles. */
    double gate_Ypos = 0.0;         /* Left edge of both grilles. */
    double gate0_Ypos = gate_Ypos;  /* Left edge of leaf grille. */
    double gate1_Ypos = gate0_Ypos + gate_Ysize;  /* Right edge of right grille. */
    double gate_Zpos = 0.0;         /* Bottom edge of leaf grilles. */

    /* Define Client coords of the plot window (mm; aspect ratio 13:10):): */
    double xSize = 3640;
    double ySize = 2800;
    double xCtr, yCtr;  /* Approx coordinates of center of gate: */
    if (cut)
      { /* Plot coordinates are {x=Y}, {y=-X}: */
        xCtr = gate_Ypos + gate_Ysize/2; 
        yCtr = gate_Xpos - gate_Xsize/2;
      }
    else
      { /* Plot coordinates are {x=Y}, {y=Z}: */
        xCtr = gate_Ypos + gate_Ysize/2; 
        yCtr = gate_Zpos + gate_Zsize/2;
      }
    double xMin = xCtr - xSize/2, xMax = xMin + xSize;
    double yMin = yCtr - ySize/2, yMax = yMin + ySize;
    epswr_set_client_window(epsf, xMin, xMax, yMin, yMax);
    
    double gate0_Ymin, gate0_Ymax;
    boag_draw_leaf_and_post
      ( epsf, cut, 0, 
        gate_Xpos, gate0_Ypos, gate_Zpos, 
        gate_Xsize, gate0_Ysize, gate_Zsize,
        Ypgap, XYpsize, Zpdepth,
        &gate0_Ymin, &gate0_Ymax
      );
      
    double gate1_Ymin, gate1_Ymax;
    boag_draw_leaf_and_post
      ( epsf, cut, 1, 
        gate_Xpos, gate1_Ypos, gate_Zpos, 
        gate_Xsize, gate1_Ysize, gate_Zsize,
        Ypgap, XYpsize, Zpdepth,
        &gate1_Ymin, &gate1_Ymax
      );
      
    if (! cut)
      { /* Show gate dimensions: */
        { /* Total width INCLUDING posts: */
          double xDim0 = gate0_Ymin, yDim0 = gate_Zpos + gate_Zsize;
          double xDim1 = gate1_Ymax, yDim1 = yDim0;
          boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, +8, +8, +16.0, 0, TRUE);
        }
        { /* Total width BETWEEN posts: */
          double xDim0 = gate0_Ymin + XYpsize, yDim0 = gate_Zpos - Zpdepth;
          double xDim1 = gate1_Ymax - XYpsize, yDim1 = yDim0;
          boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, -8, -8, -14.0, 0, TRUE);
        }
        { /* Gap between leaves: */
          double xDim0 = gate0_Ypos + gate0_Ysize, yDim0 = gate_Zpos + gate_Zsize;
          double xDim1 = gate1_Ypos - gate1_Ysize,               yDim1 = gate_Zpos + gate_Zsize;
          boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, +8, +8,  +10.0, 0, FALSE); /*@*/
        }
        { /* Total width OF GATE LEAVES plus gap between them: */
          double xDim0 = gate_Ypos,              yDim0 = gate_Zpos;
          double xDim1 = gate_Ypos + gate_Ysize, yDim1 = gate_Zpos;
          boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, -8, -8, -10.0, 0, TRUE);
        }
        { /* Total height EXCLUDING posts: */
          double xDim0 = gate0_Ypos, yDim0 = gate_Zpos;
          double xDim1 = gate0_Ypos, yDim1 = gate_Zpos + gate_Zsize;
          boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, +8, +8, +13.0, 0, TRUE); /*@*/
        }
      }
    epswr_end_figure(epsf);
  }

void boag_draw_leaf_and_post
  ( epswr_figure_t *epsf, 
    bool_t cut,
    sign_t which, 
    double Xpos, double Ypos, double Zpos, 
    double Xsize, double Ysize, double Zsize,
    double Ypgap,
    double XYpsize,
    double Zpdepth,
    double *YminP, 
    double *YmaxP
  )
  {
    double Zhinge0 = 200 - Zpdepth; /* {Z} coordinate of bottom hinge. */
    double Zhinge1 = Zsize - 200; /* {Z} coordinate of top hinge. */
    
    boag_draw_leaf(epsf, cut, which, Xpos, Ypos, Zpos, Xsize, Ysize, Zsize, Zhinge0, Zhinge1, Ypgap);
    boag_draw_post(epsf, cut, which, Xpos, Ypos, Zpos, Xsize, Zsize, Zhinge0, Zhinge1, Ypgap, XYpsize, Zpdepth);
    if (which == 0)
      { (*YminP) = Ypos - Ypgap - XYpsize;
        (*YmaxP) = Ypos + Ysize;
      }
    else
      { (*YminP) = Ypos - Ysize;
        (*YmaxP) = Ypos + Ypgap + XYpsize;
      }
  }

void boag_draw_leaf
  ( epswr_figure_t *epsf, 
    bool_t cut,
    sign_t which, 
    double Xpos, double Ypos, double Zpos, 
    double Xsize, double Ysize, double Zsize,
    double Zhinge0, double Zhinge1,
    double Ypgap
  )
  {
    sign_t dir = 1 - 2*which; /* Y direction from hinges to free edge. */
    
    /* Cross-section of grille tubing: */
    double tube_Xsize = Xsize;
    double tube_YZsize = Xsize;
    
    fprintf(stderr, "leaf %d:\n", which); 
    
    boag_draw_leaf_frame
      ( epsf, cut, 
        Xpos, Xpos + tube_Xsize, 
        Ypos, Ypos + dir*Ysize, tube_YZsize,
        Zpos, Zpos + Zsize, tube_YZsize
      );
    
    boag_draw_leaf_beams
      ( epsf, cut, 
        Xpos, Xpos + tube_Xsize, 
        Ypos + dir*tube_YZsize, Ypos + dir*(Ysize - tube_YZsize), tube_YZsize,
        Zpos + tube_YZsize, Zpos + (Zsize - tube_YZsize), tube_YZsize
      );
      
    /* Draw hinges (top halves): */
    double Xhinges = Xpos + tube_Xsize/2;
    double Yhinges = Ypos - dir*Ypgap/2;
    double Yhdist = dir*Ypgap/2; /* Singned {Y} istance from hinge axis to leaf edge. */

    boag_draw_half_hinge(epsf, cut, 1, Xhinges, Yhinges, Zhinge0, Yhdist);
    boag_draw_half_hinge(epsf, cut, 1, Xhinges, Yhinges, Zhinge1, Yhdist);
    
    if (which == 0)
      { /* Draw letter box, key lock and handle, latch stop, gap cover strip: */
      }
    else
      { /* Draw vertical latch: */
      }
      
    if (! cut)
      { /* Show leaf dimensions: */
        { /* Leaf width: */
          double xDim0_raw = Ypos,             yDim0 = Zpos + Zsize;
          double xDim1_raw = Ypos + dir*Ysize, yDim1 = Zpos + Zsize;
          double xDim0 = fmin(xDim0_raw, xDim1_raw);
          double xDim1 = fmax(xDim0_raw, xDim1_raw);
          boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, +8, +8, +4.0, 0, TRUE);
        }
      }
  }
    
void boag_draw_leaf_frame
  ( epswr_figure_t *epsf, 
    bool_t cut,
    double Xmin, double Xmax,
    double Ymin, double Ymax, double tube_Ysize,
    double Zmin, double Zmax, double tube_Zsize
  )
  {
    /* Sort interval bounds: */
    if (Xmin > Xmax) { double tmp = Xmin; Xmin = Xmax; Xmax = tmp; }
    if (Ymin > Ymax) { double tmp = Ymin; Ymin = Ymax; Ymax = tmp; }
    if (Zmin > Zmax) { double tmp = Zmin; Zmin = Zmax; Zmax = tmp; }
    
    /* {X} center and extent of all tubes: */
    double tube_Xsize = Xmax - Xmin; /* {X} extent of all tubes. */
    double Xctr = Xmin + tube_Xsize/2; /* {X} of all tubes. */

    /* Parameters of horz frame tubes: */
    double H_Yctr = (Ymin + Ymax)/2;      /* {Y} of tubes. */
    double H_Ysize = Ymax - Ymin;         /* {Y} extent of tubes. */
    double H0_Zctr = Zmin + tube_Zsize/2; /* {Z} of bot tube. */
    double H1_Zctr = Zmax - tube_Zsize/2; /* {Z} of top tube. */

    boag_draw_tube(epsf, cut, 1, Xctr, H_Yctr, H0_Zctr, tube_Xsize, H_Ysize, tube_Zsize);
    boag_draw_tube(epsf, cut, 1, Xctr, H_Yctr, H1_Zctr, tube_Xsize, H_Ysize, tube_Zsize);

    /* Parameters of vert frame tubes: */
    double V0_Yctr = Ymin + tube_Ysize/2; /* {Y} of left tube. */
    double V1_Yctr = Ymax - tube_Ysize/2; /* {Y} of right tube. */
    double V_Zctr = (Zmin + Zmax)/2;      /* {Z} of tubes. */
    double V_Zsize = Zmax - Zmin - 2*tube_Zsize;  /* {Z} extent of tubes. */
    
    boag_draw_tube(epsf, cut, 2, Xctr, V0_Yctr, V_Zctr, tube_Xsize, tube_Ysize, V_Zsize);
    boag_draw_tube(epsf, cut, 2, Xctr, V1_Yctr, V_Zctr, tube_Xsize, tube_Ysize, V_Zsize);
  }

void boag_draw_leaf_beams
  ( epswr_figure_t *epsf, 
    bool_t cut,
    double Xmin, double Xmax,
    double Ymin, double Ymax, double tube_Ysize,
    double Zmin, double Zmax, double tube_Zsize
  )
  {
    /* Start and direction of counting:*/
    sign_t dir = (Ymax > Ymin ? +1 : -1);
    double Yorg = Ymin;
    
    /* Sort interval bounds: */
    if (Xmin > Xmax) { double tmp = Xmin; Xmin = Xmax; Xmax = tmp; }
    if (Ymin > Ymax) { double tmp = Ymin; Ymin = Ymax; Ymax = tmp; }
    if (Zmin > Zmax) { double tmp = Zmin; Zmin = Zmax; Zmax = tmp; }
    
    /* {X} center and extent of all tubes: */
    double tube_Xsize = Xmax - Xmin; /* {X} extent of all tubes. */
    double Xctr = Xmin + tube_Xsize/2; /* {X} of all tubes. */

    /* Number, {Y} spacing, {Z} center, and {Z} extent of vertical grille beams: */
    double ideal_Ygap = 96.0; /* Ideal gap between vert beams. */
    int32_t YN = (int32_t)floor((Ymax - Ymin + tube_Ysize)/(ideal_Ygap + tube_Ysize) + 0.5) - 1;
    double Ygap = (Ymax - Ymin + tube_Ysize)/(YN + 1) - tube_Ysize; /* Actual gap between vert beams. */
    double V_Zctr = (Zmin + Zmax)/2; /* {Z} coord of center of vert beams. */
    double V_Zsize = Zmax - Zmin; /* {Z} extent of vert beams. */

    /* Max number and spacing of horizontal grille beams. */
    double ideal_Zgap = 315.0; /* Ideal gap between horz beams: */
    double Zex = 20; /* Extra space to trim at bottom. */
    int32_t ZN = (int32_t)floor((Zmax - Zmin + Zex + tube_Zsize)/(ideal_Zgap + tube_Zsize) + 0.5);
    double Zgap = (Zmax - Zmin + Zex + tube_Zsize)/ZN - tube_Zsize; /* Actual gap between horz beams. */
    
    fprintf(stderr, "  %2d vert beams, Ygap = %7.3f mm\n", YN, Ygap);
    fprintf(stderr, "  %2d horz beams, Zgap = %7.3f mm\n", ZN, Zgap);
    
    assert(YN %2 == 0);  /* For best looks. */

    for (int32_t kY = 0;  kY <= YN; kY++)
      { 
        if (kY > 0)
          { /* Draw ver beam {kY}: */
            double V_Ymax = Yorg + kY*dir*(Ygap + tube_Ysize);
            double V_Ymin = V_Ymax - dir*tube_Ysize;
            if (V_Ymin > V_Ymax) { double tmp = V_Ymin; V_Ymin = V_Ymax; V_Ymax = tmp; }
            double V_Yctr = (V_Ymin + V_Ymax)/2;
            boag_draw_tube(epsf, cut, 2, Xctr, V_Yctr, V_Zctr, tube_Xsize, tube_Ysize, V_Zsize);
          }

        /* Draw horz beams to the right of vert beam {kY}: */
        
        /* {Y} center of horz beams in gap: */
        double H_Ymin = Yorg + kY*dir*(Ygap + tube_Ysize);
        double H_Ymax = H_Ymin + dir*Ygap;
        if (H_Ymin > H_Ymax) { double tmp = H_Ymin; H_Ymin = H_Ymax; H_Ymax = tmp; }

        double H_Yctr = (H_Ymin + H_Ymax)/2; /* {Y} coord of center of horz beams. */
        double H_Ysize = Ygap; /* {Y} extent of horz beams. */
        int32_t H_phase = ((kY + 1) % 2); /* 0 if first vert gap is whole, 1 if half-height. */
             
        for (int32_t kZ = 0;  kZ < ZN; kZ++)
           { if ((H_phase == 1) || (kZ > 0))
               { /* Draw beam {kZ}: */
                 double H_Zctr = Zmin - Zex - tube_Zsize/2 + (kZ + 0.5*H_phase)*(Zgap + tube_Zsize);
                 boag_draw_tube(epsf, cut, 1, Xctr, H_Yctr, H_Zctr, tube_Xsize, H_Ysize, tube_Zsize);
               }
           }
      }

    if (dir == -1)
      { /* Show dimensions: */
        { /* Gap between vertical beams: */
          double xDim0 = Ymax - Ygap, yDim0 = Zmin;
          double xDim1 = Ymax,        yDim1 = Zmin;
          boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, -8, -8, -6.0, 1, FALSE);
        }
        
        { /* Vertical measurements in column 0: */
          double xDim0 = Ymax, xDim1 = xDim0;
          { /* Gap between frame and first horizonal beam: */
            double yDim0 = Zmin;
            double yDim1 = yDim0 - Zex + Zgap - 0.5*(Zgap + tube_Zsize);
            boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, +12, +12, 00.0, 0, TRUE);
          }
          { /* Gap between first and second horizonal beams: */
            double yDim0 = Zmin - Zex + 0.5*(Zgap + tube_Zsize);
            double yDim1 = yDim0 + Zgap;
            boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, +12, +12, 00.0, 0, TRUE);
          }
          { /* Gap between last two horizonal beams in column 0: */
            double yDim0 = Zmin - Zex + 4.5*(Zgap + tube_Zsize);
            double yDim1 = yDim0 + Zgap;
            boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, +12, +12, 00.0, 0, TRUE);
          }
          { /* Gap between last horiz beam and frame: */
            double yDim0 = Zmin - Zex + 5.5*(Zgap + tube_Zsize);
            double yDim1 = yDim0 + Zgap - 0.5*(Zgap + tube_Zsize);
            boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, +12, +12, 00.0, 0, TRUE);
          }
        }
       
        { /* Measurements in column 1: */
          double xDim0 = Ymax - (Ygap + tube_Ysize), xDim1 = xDim0;
          { /* Gap between frame and first horizonal beam in column 1: */
            double yDim0 = Zmin;
            double yDim1 = yDim0 - Zex + Zgap;
            boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, +12, +12, 00.0, 0, TRUE);
          }
          { /* Gap between first and second horizonal beams in column 1: */
            double yDim0 = Zmin - Zex + (Zgap + tube_Zsize);
            double yDim1 = yDim0 + Zgap;
            boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, +12, +12, 00.0, 0, TRUE);
          }
          { /* Gap between first and second horizonal beams in column 1: */
            double yDim0 = Zmin - Zex + 5*(Zgap + tube_Zsize);
            double yDim1 = yDim0 + Zgap;
            boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, +12, +12, 00.0, 0, TRUE);
          }
        }
      }
  }

void boag_draw_post
  ( epswr_figure_t *epsf, 
    bool_t cut,
    sign_t which, 
    double Xpos, double Ypos, double Zpos,
    double Xsize, double Zsize,
    double Zhinge0, double Zhinge1,
    double Ypgap,
    double XYpsize,
    double Zpdepth
  )
  {
    sign_t dir = 1 - 2*which;  /* Direction of leaf. */
    
    /* Undershoot of post: */
    double post_Zpdepth = Zpdepth; 

    /* Cross-section of supporting post: */
    double post_Xsize = XYpsize;             /* Cross-section {X} size. */
    double post_Ysize = XYpsize;             /* Cross-section {Y} size. */
    double post_Zsize = Zsize + post_Zpdepth; /* Length. */
    
    /* Post center: */
    double post_Xctr = Xpos + Xsize/2; /* Sic: {Xpos} is the back of the leaf grille. */
    double post_Yctr = Ypos - dir*(Ypgap + post_Ysize/2);
    double post_Zctr = Zpos - post_Zpdepth + post_Zsize/2;
    
    /* Draw the post: */
    boag_draw_tube
      ( epsf, cut, 2, 
        post_Xctr, post_Yctr, post_Zctr,
        post_Xsize, post_Ysize, post_Zsize
      );
      
    /* Draw hinges (bottom halves): */
    double Xhinges = Xpos + Xsize/2;
    double Yhinges = Ypos - dir*Ypgap/2;
    double Yhdist = dir*Ypgap/2;
    boag_draw_half_hinge(epsf, cut, 0, Xhinges, Yhinges, Zhinge0, Yhdist);
    boag_draw_half_hinge(epsf, cut, 0, Xhinges, Yhinges, Zhinge1, Yhdist);
      
    if (which == 0)
      { /* Show dimensions: */
        { /* Width of post: */
          double xDim0 = post_Yctr - post_Ysize/2, yDim0 = Zpos + Zsize;
          double xDim1 = post_Yctr + post_Ysize/2, yDim1 = yDim0;
          boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, +8, +8, +4.0, 0, TRUE);
        }
        { /* Gap between post and leaf: */
          double xDim0 = post_Yctr + post_Ysize/2, yDim0 = Zpos + Zsize;
          double xDim1 = Ypos,                     yDim1 = yDim0;
          boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, +8, +8, +10.0, 0, FALSE);
        }
        { /* Depth of post below grilles: */
          double xDim0 = Ypos, yDim0 = Zpos - post_Zpdepth;
          double xDim1 = Ypos, yDim1 = Zpos;
          boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, +130, +8, +7.0, 0, FALSE); /*@*/
        }
      }
    else
      { /* Dimensions on the right side. */
        { /* Distance from bottom of post to bottom hinge: */
          double xDim0 = Yhinges, yDim0 = Zpos - post_Zpdepth;
          double xDim1 = Yhinges, yDim1 = Zhinge0;
          boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, -108, -108, -6.0, 0, TRUE);
        }
        { /* Distance from top of post to top hinge: */
          double xDim0 = Yhinges, yDim0 = Zhinge1;
          double xDim1 = Yhinges, yDim1 = Zpos + Zsize;
          boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, -108, -108, -5.0, 0, TRUE);
        }
        { /* Total height of post: */
          double xDim0 = post_Yctr + post_Ysize/2, yDim0 = Zpos - post_Zpdepth;
          double xDim1 = post_Yctr + post_Ysize/2, yDim1 = Zpos + Zsize;
          boag_show_dim(epsf, xDim0, yDim0, xDim1, yDim1, -8, -8, -10.0, 0, TRUE);
        }
      }
  }

void boag_draw_half_hinge
  ( epswr_figure_t *epsf, 
    bool_t cut,
    int32_t which,
    double Xctr, double Yctr, double Zctr,
    double Yhdist
  )
  {
    sign_t dir = (Yhdist > 0 ? +1 : -1);  /* Direction of leaf. */

    double Hrad = 10.0;   /* Hinge radius. */
    double Sthk = 5.0;    /* Thickness of support plates. */
    double Slen = 23.0;   /* Length of transversal support plate. */
    double Zsize = 30;    /* Height of a half-hinge. */
    
    /* Find {Z} range of half-hinge and direction of support: */
    double Zmin, Zmax;  /* {Z} range of hinge and support. */
    sign_t hdir; /* {Y} direction from hinge axis to support. */ 

    if (which == 0)
      { /* Bottom half-hinge: */
        Zmin = Zctr - Zsize;
        Zmax = Zctr;
        hdir = -dir;
      }
    else
      { /* Top half_hinge: */
        Zmin = Zctr;
        Zmax = Zctr + Zsize;
        hdir = +dir;
      }
   
    boag_draw_vert_cylinder(epsf, cut, Xctr, Yctr, Zmin, Zmax, Hrad);
    
    /* Support plates: */
    double sup_Zctr = (Zmin + Zmax)/2; 
    double sup_Zsize = Zsize;

    /* Support plate parallel to {YZ} plane: */
    double sup0_Ysize = fabs(Yhdist) - Hrad - Sthk;
    if (sup0_Ysize > 0)
      { double sup0_Xctr = Xctr, sup0_Xsize = Sthk;
        double sup0_Yctr = Yctr + hdir*(Hrad + sup0_Ysize/2);
        boag_draw_plate(epsf, cut, sup0_Xctr, sup0_Yctr, sup_Zctr,  sup0_Xsize, sup0_Ysize, sup_Zsize);
      }

    /* Support plate parallel to {XZ} plane: */
    double sup1_Xctr = Xctr, sup1_Xsize = Slen;
    double sup1_Ya = Yctr + hdir*fabs(Yhdist), sup1_Yb = sup1_Ya - hdir*Sthk;
    double sup1_Yctr = (sup1_Ya + sup1_Yb)/2, sup1_Ysize = fabs(sup1_Ya - sup1_Yb);
    boag_draw_plate(epsf, cut, sup1_Xctr, sup1_Yctr, sup_Zctr,  sup1_Xsize, sup1_Ysize, sup_Zsize);
  }
  
void boag_draw_tube
  ( epswr_figure_t *epsf, 
    bool_t cut,
    int32_t axis,
    double Xctr, double Yctr, double Zctr, 
    double Xsize, double Ysize, double Zsize
  )
  {
    if (cut)
      { /* Draw horizontal cross-section: */
        affirm(FALSE, "not implemented");
      }
    else
      { /* Cylindrical view from {+X}: */
        demand(axis != 0, "Tube along {X} axis not implemented");
        epswr_set_fill_color(epsf, 1.000,0.800,0.400);
        epswr_set_pen(epsf, 0.000,0.000,0.000, 0.25, 0.0,0.0);
        epswr_centered_rectangle(epsf, Yctr,Zctr, Ysize,Zsize, 0, TRUE, TRUE);
      }
  }
  
void boag_draw_plate
  ( epswr_figure_t *epsf, 
    bool_t cut,
    double Xctr, double Yctr, double Zctr, 
    double Xsize, double Ysize, double Zsize
  )
  {
    if (cut)
      { /* Draw horizontal cross-section: */
        affirm(FALSE, "not implemented");
      }
    else
      { /* Cylindrical view from {+X}: */
        epswr_set_fill_color(epsf, 1.000,0.800,0.400);
        epswr_set_pen(epsf, 0.000,0.000,0.000, 0.25, 0.0,0.0);
        epswr_centered_rectangle(epsf, Yctr,Zctr, Ysize,Zsize, 0, TRUE, TRUE);
      }
  }
  
void boag_draw_vert_cylinder
  ( epswr_figure_t *epsf, 
    bool_t cut,
    double Xctr, double Yctr, 
    double Zmin, double Zmax,
    double Hrad
  )
  {
    if (cut)
      { /* Draw horizontal cross-section: */
        affirm(FALSE, "not implemented");
      }
    else
      { /* Cylindrical view from {+X}: */
        double xMin = Yctr - Hrad, xMax = Yctr + Hrad;
        double yMin = Zmin, yMax = Zmax;
        double xc = Yctr, yc = (Zmin + Zmax)/2;
        double wd = 2*Hrad, ht = Zmax - Zmin;
        epswr_set_fill_color(epsf, 1.000,0.800,0.400);
        epswr_set_pen(epsf, 0.000,0.000,0.000, 0.25, 0.0,0.0);
        epswr_centered_rectangle(epsf, xc, yc, wd, ht, 0, TRUE, TRUE);
      }
  }
  
void boag_show_dim
  ( epswr_figure_t *epsf, 
    double xa, double ya,
    double xb, double yb,
    double agap,
    double bgap, 
    double elen,
    int32_t nfrac,
    bool_t inner
  )
  { sign_t dir = ( elen < 0 ? -1 : +1); /* Direction of extension lines rel {a-->b}. */
    double dpos = 1.0;
    double dlen = (inner ? 100000 : 5.0);
    double hoff = 0.0;
    double voff = (inner ? 2.5 : dir*3.5);
    double dab, xr, yr, rot;
    epswr_set_fill_color(epsf, 0.000, 0.000, 0.000);
    epswr_set_pen(epsf, 0.000,0.000,0.000, 0.25, 0.0,0.0);
    epswr_set_verbose(epsf, TRUE);
    epswr_dim_linear
      ( epsf, xa, ya, xb, yb, &dab,
        agap, bgap, elen, 
        inner, dpos, dlen, 
        hoff, voff, &xr, &yr, &rot
      );
    char *label = NULL;
    char *label = jsprintf("%.*f", nfrac, dab);
    epswr_label(epsf, label, "0", xr, yr, rot, FALSE, 0.5, 0.5, TRUE, FALSE);
    epswr_set_verbose(epsf, FALSE);
    free(label);
  }
