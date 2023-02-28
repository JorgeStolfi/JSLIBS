/* Floorlans of IC-1+2 building */
/* Last edited on 2023-02-27 21:04:37 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include <r2.h>
#include <vec.h>
#include <epswr.h>
#include <affirm.h>
#include <jsfile.h>
#include <frgb.h>

#include <archdraw.h>
#include <archdraw_ic.h>

#define OUT_PREFIX "ic12"

adrw_point_vec_t adrw_ic12_define_points(void);
int32_t adrw_ic12_get_office_number(int32_t wing, int32_t side, int32_t slot);
int32_t adrw_ic12_get_first_corner_of_office(int32_t noff);
void adrw_ic12_get_office_corners(int32_t wing, int32_t side, int32_t slot, r2_t *p00, r2_t *p01, r2_t *p10, r2_t *p11);
void adrw_ic12_get_office_type_and_span(int32_t noff, adrw_ic_space_type_t *toffP, int32_t *xspanP);
char **adrw_ic12_get_office_descriptions(void);

void adrw_ic12_append_building_outline(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]);
void adrw_ic12_append_pillars(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]);
void adrw_ic12_append_hallways(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]);
void adrw_ic12_append_offices(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[], char *descr[]);
void adrw_ic12_append_normal_office
  ( adrw_building_t *B, 
    adrw_point_vec_t *P, 
    int32_t wing, /* Wing (0 = south, 1 = north). */
    int32_t side, /* Side of corridor (0 = south, 1 = north). */
    int32_t slot, /* Office slot along one side of a wing (0 = easternmost). */
    adrw_unit_style_t *style[], /* Plot styles, indexed by type. */
    char *descr[]              /* Office descriptions, indexed by {noff}. */
  );
void adrw_ic12_append_coord_courses_office(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[], char *descr[]);
void adrw_ic12_append_secr_courses(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[], char *descr[]);
void adrw_ic12_append_central_cabinet(adrw_building_t *B, adrw_point_vec_t *P, int32_t wing,adrw_unit_style_t *style[], char *descr[]);
void adrw_ic12_append_auditorium_seating(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]);
void adrw_ic12_append_entrance_arrow(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]);

void adrw_ic12_plot_all(adrw_building_t *B, int32_t nx, int32_t ny, bool_t show_dots);

#define ADDPOLY(PTS,KFL,LABEL,DESC,NMODS,TYPE,...)      \
  do \
    { int32_t v[] = { __VA_ARGS__ , -1 }; \
      adrw_unit_t *rm = adrw_make_poly((LABEL),(DESC),(NMODS),(PTS),v,0.0,(KFL),(TYPE),style[(TYPE)]); \
      adrw_append_unit(B, rm); \
    } \
  while(0)

#define ADDSEATS(PTS,KFL,TYPE,PTA,PTB,SZX,SZY)       \
  do \
    { adrw_append_seats(B,(PTS),(PTA),(PTB),(SZX),(SZY),(KFL),(TYPE),style[(TYPE)]); } \
  while(0)

#define ADDBOX(PTS,KFL,LABEL,DESC,TYPE,CTR,WDX,WDY)      \
  do \
    { adrw_unit_t *rm = adrw_make_box((LABEL),(DESC),(PTS),(CTR),0.0,0.0,0.0,(WDX),(WDY),0.0,(KFL),(TYPE),style[(TYPE)]); \
      adrw_append_unit(B, rm); \
    } \
  while(0)

#define ADDPIPE(PTS,KFL,LABEL,DESC,TYPE,...) \
  do \
    { int32_t v[] = { __VA_ARGS__ , -1 }; \
      adrw_unit_t *rm = adrw_make_poly((LABEL),(DESC),0.0,(PTS),v,0.0,(KFL),(TYPE),style[(TYPE)]); \
      adrw_append_unit(B, rm); \
    } \
  while(0)

int32_t main (int32_t argc, char **argv)
  {
    int32_t ntypes = adrw_ic_space_type_MAX+1;
    adrw_unit_style_t **style = adrw_ic_define_type_styles();
    char **type_tag = adrw_ic_define_type_tags();
    bool_t *type_is_printable = adrw_ic_define_printable_types();
    bool_t *type_is_movable = adrw_ic_define_movable_types();
    int32_t *color_key_types = adrw_ic_define_color_key_types();
    
    adrw_point_vec_t P = adrw_ic12_define_points();
    
    char **descr = adrw_ic12_get_office_descriptions();
    adrw_building_t *B = adrw_make_building();
    
    adrw_ic12_append_building_outline(B, &P, style);
    adrw_ic12_append_hallways(B, &P, style);
    adrw_ic12_append_offices(B, &P, style, descr);
    adrw_ic12_append_pillars(B, &P, style);
    adrw_ic12_append_entrance_arrow(B, &P, style);
    
    /* Print the point table: */
    { FILE *wr = open_write("out/" OUT_PREFIX "_P.txt", TRUE);
      adrw_print_points(wr, &P);
      fclose(wr);
    }

    /* Plot the floorplan etc.: */
    { /* Plot the floorplan: */
      adrw_ic12_plot_all(B, 1, 1, FALSE);
      /* Plot the legend: */
      adrw_plot_type_legend("out", OUT_PREFIX, "T", color_key_types, ntypes, 5, type_tag, style);
      /* Plot the unit area bars: */
      adrw_ic_plot_histogram_bars("out", OUT_PREFIX, B, type_is_movable, type_tag, style);
      /* adrw_ic12_plot_all(epsf, B, 1, 1, TRUE); */
    }

    /* Print the units table: */
    { FILE *wr = open_write("out/" OUT_PREFIX "_U.tex", TRUE);
      bool_t TeX = TRUE;
      adrw_print_building(wr, B, type_is_printable, type_tag, TRUE, TRUE, FALSE, TeX);
      fclose(wr);
    }
    
    return 0;
  }

void adrw_ic12_plot_all(adrw_building_t *B, int32_t nx, int32_t ny, bool_t show_dots)
  {
    /* Building dimensions: */
    double xwid = 7200;
    double ywid = 4000;
    /* Plot domain in cm, for all three floors. */
    double xmrg = 300;
    double ymrg = 300;
    double xmin = 00 - xmrg, xmax = xwid + xmrg;
    double ymin = 00 - ymrg, ymax = ywid + ymrg;

    for (int32_t ox = 0; ox < nx; ox++)
      { for (int32_t oy = 0; oy < ny; oy++)
          {
            fprintf(stderr, "=== PLOTTING PAGE [%d,%d] OF [%d,%d] ===\n", ox,oy,nx,ny);
            epswr_figure_t *epsf = adrw_new_figure
              ( "out", OUT_PREFIX "_A", "P", xmin, xmax, ymin, ymax,
                ox, nx, oy, ny, "Atual - Predio IC-1+2"
              );
            adrw_plot_building(epsf, B, show_dots);
            epswr_end_figure(epsf);
          }
      }
  }
   
int32_t adrw_ic12_get_office_number(int32_t wing, int32_t side, int32_t slot)
  { return wing*50 + slot*2 + side + 1; }

int32_t adrw_ic12_get_first_corner_of_office(int32_t noff)
  { return 400 + 4*noff; }
 
adrw_point_vec_t adrw_ic12_define_points(void)
  {
    adrw_point_vec_t P = adrw_point_vec_new(100);
    int32_t np = 0;
    
    auto void s(int32_t j, double dX, double dY, int32_t i);
      /* Defines point {P[j]+(dX,dY)} as {P[i]}.  If {P[i]} was defined
        previously, checks whether the definitions agree. */
       
    void s(int32_t j, double dX, double dY, int32_t i)
      { char *lab = NULL;
        asprintf(&lab, "P%03d", i);
        adrw_append_point(lab, i, j, j, j, dX, dY, 0.0, &P, &np); 
      } 
    
    /* All dimensions in centimeters. */
    
    double Xor = 00; /* X of W edge of IC1. */
    double Yor = 00; /* Y of S edge of IC1. */

    fprintf(stderr, "=== DEFINING POINTS ===\n");
    
    fprintf(stderr, "--- main corners ---\n");
    s( -1,   Xor,   Yor, 101);
    s(101, +7220,    00, 102);
    s(101,    00, +1520, 106);
    s(102,    00, +1520, 103);
    
    s(101,    00, +4020, 201); 
    s(201, +7220,    00, 202);
    s(201,    00, -1520, 206);
    s(202,    00, -1520, 203);
    
    fprintf(stderr, "--- pillars ---\n");
    s(101,   +10,   +10, 160);
    s(160,  +600,    00, 161);
    s(161,  +600,    00, 162);
    s(162,  +600,    00, 163);
    s(163,  +600,    00, 164);
    s(164,  +600,    00, 165);
    s(165,  +600,    00, 166);
    s(166,  +600,    00, 167);
    s(167,  +600,    00, 168);
    s(168,  +600,    00, 169);
    s(169,  +600,    00, 170);
    s(170,  +600,    00, 171);
    s(171,  +600,    00, 172);
    s(160, +7200,    00, 172); /* Redundant check. */
     
    s(106,   +10,   -10, 180);
    s(180,  +600,    00, 181);
    s(181,  +600,    00, 182);
    s(182,  +600,    00, 183);
    s(183,  +600,    00, 184);
    s(184,  +600,    00, 185);
    s(185,  +600,    00, 186);
    s(186,  +600,    00, 187);
    s(187,  +600,    00, 188);
    s(188,  +600,    00, 189);
    s(189,  +600,    00, 190);
    s(190,  +600,    00, 191);
    s(191,  +600,    00, 192);
    s(180, +7200,    00, 192); /* Redundant check. */
    s(160, +7200, +1500, 192); /* Redundant check. */

    s(201,   +10,   -10, 260);
    s(260,  +600,    00, 261);
    s(261,  +600,    00, 262);
    s(262,  +600,    00, 263);
    s(263,  +600,    00, 264);
    s(264,  +600,    00, 265);
    s(265,  +600,    00, 266);
    s(266,  +600,    00, 267);
    s(267,  +600,    00, 268);
    s(268,  +600,    00, 269);
    s(269,  +600,    00, 270);
    s(270,  +600,    00, 271);
    s(271,  +600,    00, 272);
    s(260, +7200,    00, 272); /* Redundant check. */
     
    s(206,   +10,   +10, 280);
    s(280,  +600,    00, 281);
    s(281,  +600,    00, 282);
    s(282,  +600,    00, 283);
    s(283,  +600,    00, 284);
    s(284,  +600,    00, 285);
    s(285,  +600,    00, 286);
    s(286,  +600,    00, 287);
    s(287,  +600,    00, 288);
    s(288,  +600,    00, 289);
    s(289,  +600,    00, 290);
    s(290,  +600,    00, 291);
    s(291,  +600,    00, 292);
    s(280, +7200,    00, 292); /* Redundant check. */
    s(260, +7200, -1500, 292); /* Redundant check. */

    fprintf(stderr, "--- building outlines ---\n");
    s(103, -2700,    00, 104);
    s(188,  -290,   +10, 104); /* Redundant check. */
    
    s(203, -2700,    00, 204);
    s(288,  -290,   -10, 204); /* Redundant check. */
    
    s(186,  +395,   +10, 105);
    s(286,  +395,   -10, 205);
    
    s(187,   +10,  -600, 115);
    s(287,   +10,  +600, 215);
   
    s(105,    00,  +165, 107);
    s(107,  -425,    00, 108);
    s(108,    00,  +150, 109);
    s(109,  -615,    00, 110);

    s(205,    00,  -165, 207);
    s(207,  -425,    00, 208);
    s(208,    00,  -150, 209);
    s(209,  -615,    00, 210);
    s(210,    00,  -350, 110); /* Redundant check. */
    
    /* Assumes that exterior walls are 20 cm thick. */
    /* Assumes that interior brick walls are 16 cm thick. */
    /* Assumes that interior board walls are 6 cm thick. */

    fprintf(stderr, "--- hallways ---\n");
    s(101,   +20,  +610, 111);
    s(102,   -20,  +610, 112);
    s(106,   +20,  -610, 113);
    s(103,   -20,  -610, 114);
    s(104,   -20,  -610, 116);
    
    s(201,   +20,  -610, 211);
    s(211,  +287,    00, 221);
    s(221,    00,   -40, 222);
    s(222, +1811,    00, 217);
    s(217,    00,   +40, 218);
    s(202,   -20,  -610, 212);
    s(206,   +20,  +650, 213);
    s(213, +2098,    00, 223);
    s(223,    00,   -30, 224);
    s(203,   -20,  +610, 214);
    s(204,   -20,  +610, 216);
    
    fprintf(stderr, "--- office corners for main wings ---\n");
    int32_t wing; /* Wing (0 = south, 1 = north). */
    int32_t side; /* Side of corridor (0 = south, 1 = north). */
    int32_t slot; /* Office slot along one side of a wing (0 = easternmost). */
    for (wing = 0; wing < 2; wing++)
      { for (side = 0; side < 2; side++)
          { for (slot = 0; slot < 24; slot++)
              { int32_t noff = adrw_ic12_get_office_number(wing,side,slot); /* Office number (regularized). */
                /* Note that some numbers are changed here because {xspan} is always westward: */
                /*   LIS   is room 95 (labeled 93), spanning 95 and 97. */
                /*   LSD   is room 91 (labeled 91), spanning 91 and 93. */
                /*   LRC   is room 96 (labeled 96), spanning 96 and 98. */
                /*   LAS   is room 82 (labeled 84), spanning 82 and 84. */
                /*   LIV   is room 72 (labeled 74), spanning 72 and 74. */
                /*   LSC   is room 66 (labeled 68), spanning 66 and 68. */
                /*   Dr0   is room 71 (labeled 71), spanning 71 and 73. */
                /*   Dr1   is room 86 (labeled ??), spanning 86 and 88. */
                /*   Dr2   is room 90 (labeled 90). */
                /*   Dr3   is room 92 (labeled 92), spanning 92 and 94. */
                /*   Copa  is room 46 (labeled 46), spanning 46 and 48. */
                /*   Aud   is room 85 (labeled 85), spanning 85, 87, and 89. */
                /*   These glitches should be fixed by allowing negative span. */
                int32_t fco = adrw_ic12_get_first_corner_of_office(noff); /* Number of easternmost outer corner. */
                fprintf(stderr, "  office %d - corners %d..%d\n", noff, fco, fco + 3);
                
                /* Get X and Y cords of office (excluding walls) rel to floor origin: */
                r2_t p00;  /* Outer east corner of office: */ 
                r2_t p01;  /* Inner east corner of office: */ 
                r2_t p10;  /* Outer west corner of office: */   
                r2_t p11;  /* Inner west corner of office: */   
                adrw_ic12_get_office_corners(wing, side, slot, &p00, &p01, &p10, &p11);

                /* Set the corners: */
                s(101, p00.c[0], p00.c[1], fco+0);
                s(101, p01.c[0], p01.c[1], fco+1);
                s(101, p10.c[0], p10.c[1], fco+2);
                s(101, p11.c[0], p11.c[1], fco+3);
              }
          }
      }

    fprintf(stderr, "--- auditorium seating ---\n");
    int32_t auc = adrw_ic12_get_first_corner_of_office(85);
    s(auc+0,   00,  +12, 130);
    s(130,     00, +240, 131);
    s(131,   -640,   00, 132);
    s(132,     00, -240, 133);
    s(130,   -640,   00, 133); /* Redundant check. */
    
    s(auc+1,   00,  -12, 134);
    s(134,     00, -240, 135);
    s(135,   -640,   00, 136);
    s(136,     00, +240, 137);
    s(134,   -640,   00, 137); /* Redundant check. */
    
    s(135,    -80,   00, 138); /* Lone seats in back. */
    
    fprintf(stderr, "--- key corners in central block ---\n");
    s(187,  -10,  +10, 800);
    s(287,  -10,  -10, 900);
    
    s(800,   00, +165, 810);
    s(900,   00, -165, 910);
    
    fprintf(stderr, "--- office corners coord courses ---\n");
    s(110,  +20,  +20, 801);
    s(210,  +20,  -20, 901);
    s(109,   00,  +20, 802);
    s(209,   00,  -20, 902);
    
    fprintf(stderr, "--- office corners xerox room ---\n");
    s(800,   00,   00, 805);
    s(810,   00, +145, 806);
    s(105,  +20,   00, 807);
    s(107,  +20, +145, 808);
    
    fprintf(stderr, "--- office corners deposit room ---\n");
    s(900,   00,  00, 905);
    s(910,   00,  -5, 906);
    s(205,  +20,  00, 907);
    s(207,  +20,  -5, 908);
    
    fprintf(stderr, "--- office corners secr. courses ---\n");
    s(108,  +20,  +20, 804);
    s(107,   +5,  +20, 812);
    s(808,  -15,  +15, 813);
    s(806,   00,  +15, 803);
    s(910,   00,  -20, 903);
    s(208,  +20,  -20, 904);
    
    fprintf(stderr, "--- entrance arrow ---\n");
    s(104, -160, +490, 850);
    s(850,  +50,  +50, 851);
    s(851,   00,  -30, 852);
    s(852, +180,   00, 853);
    s(853,   00,  -40, 854);
    s(854, -180,   00, 855);
    s(855,   00,  -30, 856);
    s(850,  +50,  -50, 856); /* Redundant check. */

    show_holes(&P, 0,np-1);

    adrw_point_vec_trim(&P, np);
    
    return P;
  }

void adrw_ic12_get_office_corners(int32_t wing, int32_t side, int32_t slot, r2_t *p00, r2_t *p01, r2_t *p10, r2_t *p11)
  {
    int32_t noff = adrw_ic12_get_office_number(wing,side,slot); /* Office number (regularized). */

    /* X and Y cords of office (excluding walls): */
    double xe, xw;    /* East and west X coords of office: */
    double yow;  /* Actual Y coord of outer west corner of office: */
    double yiw;  /* Actual Y coord of inner west corner of office: */
    double yoe;  /* Actual Y coord of outer east corner of office: */
    double yie;  /* Actual Y coord of inner east corner of office: */

    /* Compute coords for typical office: */
    xe = 10 + 300*(24-slot) - 8;
    xw = xe - 284;
    if (side == 0)
      { yow = yoe = 10 + 2500*wing + 10;
        yiw = yie = yow + 574;
      }
    else
      { yow = yoe = 10 + 1500 + 2500*wing - 10;
        yiw = yie = yow - 574;
      }

    /* Account for thicker walls at ends and adjacent to central hallway: */
    if (slot == 0)
      { xe += -2; }
    else if (slot == 23)
      { xw += +2; }
    else if ((slot == 9) && (side != wing))
      { xw += +2; }
    else if ((slot == 11) && (side != wing))
      { xe += -2; }
    
    /* Account for displaced walls: */
    if (noff == 17)
      { xw += 150; }
    else if (noff == 19)
      { xw += 200; xe += 150; }
    else if (noff == 21)
      { xe += 200; }
    
    
    /* Account for displaced walls in PhD rooms: */
    if (noff == 86)
      { xw += +50; }
    else if (noff == 88)
      { xw += +100; xe += +50; }
    else if (noff == 90)
      { xw += -100; xe += +100; }
    else if (noff == 92)
      { xw += -50; xe += -100; }
    if (noff == 94)
      { xe += -50; }
    
    /* Account for thinner board walls (6 cm instead of 16 cm): */
    if ((noff == 85) || (noff == 86))
      { xw += -5; }
    else if ((noff >= 87) && (noff <= 96))
      { xw += -5; xe += +5;}
    else if ((noff == 97) || (noff == 98)) 
      { xe += +5; }
      
    /* Apply office-specific Y adjustments to corners: */
    if (wing == 1)
      { if (side == 1)
          { if ((slot >= 17) && (slot <= 22)) 
              { yiw += -50; yie += -50; }
            else if (slot == 23)
              { yiw += -10; yie += -50; }
          }
        else
          { if ((slot >= 17) && (slot <= 23)) 
              { yiw += +50; yie += +50; }
          }
      }

    fprintf(stderr, "  x =  [ %5.0f _ %5.0f ]\n", xw, xe);
    fprintf(stderr, "  yw = [ %5.0f _ %5.0f ]\n", yow, yiw);
    fprintf(stderr, "  ye = [ %5.0f _ %5.0f ]\n", yoe, yie);
    
    /* Return results (east to west, outer then inner): */
    (*p00) = (r2_t){{ xe, yoe }};
    (*p01) = (r2_t){{ xe, yie }};
    (*p10) = (r2_t){{ xw, yow }};
    (*p11) = (r2_t){{ xw, yiw }};
  }

void adrw_ic12_append_pillars(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[])
  {
    fprintf(stderr, "--- pillars ---\n");
    int32_t wing; /* Wing (0 = south, 1 = north). */
    int32_t side; /* Side of wing (0 = south, 1 = north). */
    int32_t pill; /* Pillar index (0 = westernmost). */
    for (wing = 0; wing < 2; wing++)
      { for (side = 0; side < 2; side++)
          { for (pill = 0; pill <= 12; pill++)
              { 
                int32_t ctr = 160 + 100*wing + 20*side + pill;
                int32_t npil = 26*wing + 13*side + pill;
                char *lpil = NULL;
                asprintf(&lpil, "P%02d", npil);
                char *dpil = NULL;
                asprintf(&dpil, "Pilar %02d", npil);
                ADDBOX
                  ( P, 0, lpil, dpil, adrw_ic_space_type_PIL,
                    ctr, 16, 20
                  );
              }
          }
      }
  }

void adrw_ic12_append_building_outline(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[])
  {
    fprintf(stderr, "--- building outline ---\n");
    ADDPOLY
      ( P, 0, "Ext", "Building outline", 0.0, adrw_ic_space_type_ETC,
        101, 102, 103, 104, 204,
        203, 202, 201, 206, 205,
        207, 208, 209, 210, 110,
        109, 108, 107, 105, 106,
        101 
      );
  }

void adrw_ic12_append_hallways(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[])
  {
    fprintf(stderr, "--- hallways ---\n");
    ADDPOLY
      ( P, 0, "Cor", "Corredores", 0.0, adrw_ic_space_type_HAL,
        111, 112, 114, 116, 216,
        214, 212, 218, 217, 222, 
        211, 213, 223, 224, 215, 
        115, 113, 111
      );
  }
  
void adrw_ic12_append_offices(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[], char *descr[])
  {
    int32_t wing; /* Wing (0 = south, 1 = north). */
    int32_t side; /* Side of corridor (0 = south, 1 = north). */
    int32_t slot; /* Office slot along one side of a wing (0 = easternmost). */
    for (wing = 0; wing < 2; wing++)
      { for (slot = 0; slot < 24; slot++)
          { for (side = 0; side < 2; side++)
              { adrw_ic12_append_normal_office(B, P, wing, side, slot, style, descr); }
          }
      }
    adrw_ic12_append_auditorium_seating(B, P, style);
    adrw_ic12_append_coord_courses_office(B, P, style, descr);
    adrw_ic12_append_secr_courses(B, P, style, descr);
    adrw_ic12_append_central_cabinet(B, P, 0, style, descr);
    adrw_ic12_append_central_cabinet(B, P, 1, style, descr);
  }

void adrw_ic12_append_normal_office
  ( adrw_building_t *B, 
    adrw_point_vec_t *P, 
    int32_t wing, /* Wing (0 = south, 1 = north). */
    int32_t side, /* Side of corridor (0 = south, 1 = north). */
    int32_t slot, /* Office slot along one side of a wing (0 = easternmost). */
    adrw_unit_style_t *style[], /* Plot styles, indexed by type. */
    char *descr[]              /* Office descriptions, indexed by {noff}. */
  )
  {              
    int32_t noff = adrw_ic12_get_office_number(wing, side, slot); /* Office number. */
    char *loff = NULL;
    asprintf(&loff, "%02d", noff);
    char *doff = descr[noff];
    /*  Decide office type {toff} and span {xspan,yspan}: */
    adrw_ic_space_type_t toff; /* Office type. */
    int32_t xspan; /* How many modules it spans in X direction (1 or more). */
    adrw_ic12_get_office_type_and_span(noff, &toff, &xspan);
    double modules = xspan;
    if (noff == 17)
      { modules = 0.5; }
    else if (noff == 21)
      { modules = 1.5; }
    /* Plot the office: */                      
    if (toff == adrw_ic_space_type_NEX)
      { /* Skip it. */ }
    else if (xspan == 1)
      { /* Get the numbers of the four corners, decreasing X. */
        int32_t fco = adrw_ic12_get_first_corner_of_office(noff); /* Number of easternmost outer corner. */
        int32_t c00 = fco + 0;  /* X0, Yot0. */
        int32_t c01 = c00 + 1;  /* X0, Yin0. */
        int32_t c10 = c00 + 2;  /* X1, Yot1. */
        int32_t c11 = c01 + 2;  /* X1, Yin1. */
        fprintf(stderr, "--- office %d - corners %d,%d,%d,%d ---\n", noff, c00,c01,c10,c11);
        ADDPOLY
          ( P, 0, loff, doff, modules, toff,
            c01, c11, c10, c00, c01
          ); 
      }
    else if (xspan == 2)
      { /* Get the numbers of the eight corners, decreasing X. */
        int32_t fco = adrw_ic12_get_first_corner_of_office(noff); /* Number of easternmost outer corner. */
        int32_t c00 = fco + 0;  /* X0, Yot0. */
        int32_t c01 = c00 + 1;  /* X0, Yin0. */
        int32_t c10 = c00 + 2;  /* X1, Yot1. */
        int32_t c11 = c01 + 2;  /* X1, Yin1. */
        int32_t c20 = c00 + 8;  /* X3, Yot3. */
        int32_t c21 = c01 + 8;  /* X3, Yin3. */
        int32_t c30 = c20 + 2;  /* X2, Yot2. */
        int32_t c31 = c21 + 2;  /* X2, Yin2. */
        fprintf
          ( stderr, "--- office %d - corners %d %d  %d %d  %d %d  %d %d ---\n", 
            noff, 
            c00, c01, c10, c11, 
            c20, c21, c30, c31
          );
        ADDPOLY
          ( P, 0, loff, doff, modules, toff,
            c01, c11, c21, c31,
            c30, c20, c10, c00,
            c01 
          );
      }
    else if (xspan == 3)
      { /* Get the numbers of the 12 corners, decreasing X. */
        int32_t fco = adrw_ic12_get_first_corner_of_office(noff); /* Number of easternmost outer corner. */
        int32_t c00 = fco + 0;  /* X0, Yot0. */
        int32_t c01 = c00 + 1;  /* X0, Yin0. */
        int32_t c10 = c00 + 2;  /* X1, Yot1. */
        int32_t c11 = c01 + 2;  /* X1, Yin1. */
        int32_t c20 = c00 + 8;  /* X2, Yot2. */
        int32_t c21 = c01 + 8;  /* X2, Yin2. */
        int32_t c30 = c20 + 2;  /* X3, Yot3. */
        int32_t c31 = c21 + 2;  /* X3, Yin3. */
        int32_t c40 = c20 + 8;  /* X4, Yot4. */
        int32_t c41 = c21 + 8;  /* X4, Yin4. */
        int32_t c50 = c40 + 2;  /* X5, Yot5. */
        int32_t c51 = c41 + 2;  /* X5, Yin5. */
        fprintf
          ( stderr, "--- office %d - corners %d %d  %d %d  %d %d  %d %d  %d %d  %d %d ---\n", 
            noff, 
            c00, c01, c10, c11, 
            c20, c21, c30, c31, 
            c40, c41, c50, c51 
          );
        ADDPOLY
          ( P, 0, loff, doff, modules, toff,
            c01, c11, c21, c31, c41, c51, 
            c50, c40, c30, c20, c10, c00,
            c01
          );
      }
    else
      { affirm(FALSE, "invalid xspan/yspan combination"); }
  }
  
void adrw_ic12_append_coord_courses_office(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[], char *descr[]) 
  {
    int32_t noff = 103; /* Office number, actually 20C */
    /* Get the numbers of the four corners, decreasing X. */
    int32_t c10 = 801;  /* X0, Ylo0. */
    int32_t c11 = 901;  /* X0, Yhi0. */
    int32_t c00 = 802;  /* X1, Ylo1. */
    int32_t c01 = 902;  /* X1, Yhi1. */
    fprintf(stderr, "--- office %d - coord courses - corners %d %d  %d %d\n", noff, c00, c01, c10, c11);
    ADDPOLY
      ( P, 0, "20C", descr[noff], 1.0, adrw_ic_space_type_ADM,
        c01, c11, c10, c00, c01
      ); 
  }
  
void adrw_ic12_append_secr_courses(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[], char *descr[]) 
  {
    int32_t noff = 100; /* Office number, actually 20 */
    /* Get the numbers of the six corners, ccw order. */
    int32_t c0 = 804;  
    int32_t c1 = 812;
    int32_t c2 = 813;
    int32_t c3 = 803;
    int32_t c4 = 903;
    int32_t c5 = 904;
    fprintf(stderr, "--- office %d - secr courses corners %d %d %d %d %d %d\n", noff, c0, c1, c2, c3, c4, c5);
    ADDPOLY
      ( P, 0, "20", descr[noff], 2.0, adrw_ic_space_type_ADM,
        c0, c1, c2, c3, c4, c5, c0
      ); 
  }
  
void adrw_ic12_append_central_cabinet(adrw_building_t *B, adrw_point_vec_t *P, int32_t wing,adrw_unit_style_t *style[], char *descr[])
  {
    int32_t noff = 101 + wing; /* Office number, actually 20A or 20B */
    char *loff = (wing == 0 ? "20A" : "20B");
    /* Get the numbers of the four corners, decreasing X. */
    int32_t c10 = 805 + 100*wing;  /* Xhi, Ylo0. */
    int32_t c11 = 806 + 100*wing;  /* Xhi, Yhi0. */
    int32_t c00 = 807 + 100*wing;  /* Xlo, Ylo1. */
    int32_t c01 = 808 + 100*wing;  /* Xlo, Yhi1. */
    fprintf(stderr, "--- office %d - central cabinet - wing %d - corners %d %d %d %d\n", noff, wing, c00, c01, c10, c11);
    adrw_ic_space_type_t toff = (wing == 0 ? adrw_ic_space_type_SRV : adrw_ic_space_type_DEP);
    ADDPOLY
      ( P, 0, loff, descr[noff], 0.5, toff,
        c01, c11, c10, c00, c01
      ); 
  }
  
void adrw_ic12_append_auditorium_seating(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]) 
  {
    int32_t noff = 85; /* Office number, spans 85, 87, 89. */
    int32_t bank; /* Bank of seats (0 or 1). */
    for (bank = 0; bank < 2; bank++)
      { /* Get the numbers of the four corners in cyclic order. */
        int32_t fco = 130 + 4*bank; /* Index of back corner. */
        int32_t c0 = fco + 0; /* Index of a back corner. */
        int32_t c2 = fco + 2; /* Index of opposite front corner. */
        fprintf(stderr, "--- auditorium %d seats - bank %d - corners %d %d\n", noff, bank, c0, c2);
        ADDSEATS
          ( P, 0, adrw_ic_space_type_SIT,
            c0, c2, 80, 60
          ); 
      }
    { int32_t c0 = 131;
      int32_t c2 = 138;
      fprintf(stderr, "--- auditorium %d seats lone seat in back - corners %d %d\n", noff, c0, c2);
      ADDSEATS
        ( P, 0, adrw_ic_space_type_SIT,
          c0, c2, 80, 60
        );
    }
  }

void adrw_ic12_append_entrance_arrow(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[])
  {
    /* Get the numbers of the 7 corners in cyclic order. */
    int32_t fco = 850; /* Index of arrow tip. */
    int32_t c0 = fco + 0;
    int32_t c1 = fco + 1;
    int32_t c2 = fco + 2;
    int32_t c3 = fco + 3;
    int32_t c4 = fco + 4;
    int32_t c5 = fco + 5;
    int32_t c6 = fco + 6;
    fprintf(stderr, "--- entrance arrow - corners %d %d %d %d %d %d %d\n", c0, c1, c2, c3, c4, c5, c6);
    ADDPOLY
      ( P, 0, "Ent", "Entrada", 0.0, adrw_ic_space_type_ARR,
        c0, c1, c2, c3, c4, c5, c6, c0
      ); 
  }

void adrw_ic12_get_office_type_and_span(int32_t noff, adrw_ic_space_type_t *toffP, int32_t *xspanP)
  {
    adrw_ic_space_type_t toff = adrw_ic_space_type_ETC; /* Office type. */
    int32_t xspan = 1; /* How many modules it spans in -X direction (1 or more). */
    switch(noff)
      { 
        case   1 : toff = adrw_ic_space_type_DOC; break;
        case   2 : toff = adrw_ic_space_type_DOC; break;
        case   3 : toff = adrw_ic_space_type_DOC; break;
        case   4 : toff = adrw_ic_space_type_DOC; break;
        case   5 : toff = adrw_ic_space_type_DOC; break;
        case   6 : toff = adrw_ic_space_type_DOC; break;
        case   7 : toff = adrw_ic_space_type_DOC; break;
        case   8 : toff = adrw_ic_space_type_DOC; break;
        case   9 : toff = adrw_ic_space_type_DOC; break;
        case  10 : toff = adrw_ic_space_type_DOC; break;
        case  11 : toff = adrw_ic_space_type_DOC; break;
        case  12 : toff = adrw_ic_space_type_DOC; break;
        case  13 : toff = adrw_ic_space_type_DOC; break;
        case  14 : toff = adrw_ic_space_type_DOC; break;
        case  15 : toff = adrw_ic_space_type_DOC; break;
        case  16 : toff = adrw_ic_space_type_DOC; break;
        case  17 : toff = adrw_ic_space_type_DEP; break;
        case  18 : toff = adrw_ic_space_type_DOC; break;
        case  19 : toff = adrw_ic_space_type_REU; break;
        case  20 : toff = adrw_ic_space_type_NEX; break;
        case  21 : toff = adrw_ic_space_type_SRV; break;
        case  22 : toff = adrw_ic_space_type_BAN; break;
        case  23 : toff = adrw_ic_space_type_DOC; break;
        case  24 : toff = adrw_ic_space_type_BAN; break;
        case  25 : toff = adrw_ic_space_type_DOC; break;
        case  26 : toff = adrw_ic_space_type_DOC; break;
        case  27 : toff = adrw_ic_space_type_DOC; break;
        case  28 : toff = adrw_ic_space_type_DOC; break;
        case  29 : toff = adrw_ic_space_type_DOC; break;
        case  30 : toff = adrw_ic_space_type_DOC; break;
        case  31 : toff = adrw_ic_space_type_DOC; break;
        case  32 : toff = adrw_ic_space_type_DOC; break;
        case  33 : toff = adrw_ic_space_type_DOC; break;
        case  34 : toff = adrw_ic_space_type_DOC; break;
        case  35 : toff = adrw_ic_space_type_DOC; break;
        case  36 : toff = adrw_ic_space_type_DOC; break;
        case  37 : toff = adrw_ic_space_type_DOC; break;
        case  38 : toff = adrw_ic_space_type_DOC; break;
        case  39 : toff = adrw_ic_space_type_DOC; break;
        case  40 : toff = adrw_ic_space_type_DOC; break;
        case  41 : toff = adrw_ic_space_type_DOC; break;
        case  42 : toff = adrw_ic_space_type_DOC; break;
        case  43 : toff = adrw_ic_space_type_DOC; break;
        case  44 : toff = adrw_ic_space_type_DOC; break;
        case  45 : toff = adrw_ic_space_type_DOC; break;
        case  46 : toff = adrw_ic_space_type_REU; xspan = 2; break;
        case  47 : toff = adrw_ic_space_type_DOC; break;
        case  48 : toff = adrw_ic_space_type_NEX; break;
        case  51 : toff = adrw_ic_space_type_REU; xspan = 2; break;
        case  52 : toff = adrw_ic_space_type_ADM; xspan = 2; break;
        case  53 : toff = adrw_ic_space_type_NEX; break;
        case  54 : toff = adrw_ic_space_type_NEX; break;
        case  55 : toff = adrw_ic_space_type_ADM; break;
        case  56 : toff = adrw_ic_space_type_ADM; break;
        case  57 : toff = adrw_ic_space_type_ADM; break;
        case  58 : toff = adrw_ic_space_type_ADM; break;
        case  59 : toff = adrw_ic_space_type_ADM; break;
        case  60 : toff = adrw_ic_space_type_POS; break;
        case  61 : toff = adrw_ic_space_type_INF; break;
        case  62 : toff = adrw_ic_space_type_POS; break;
        case  63 : toff = adrw_ic_space_type_INF; break;
        case  64 : toff = adrw_ic_space_type_POS; break;
        case  65 : toff = adrw_ic_space_type_ADM; break;
        case  66 : toff = adrw_ic_space_type_POS; xspan = 2; break;
        case  67 : toff = adrw_ic_space_type_ADM; break;
        case  68 : toff = adrw_ic_space_type_NEX; break;
        case  69 : toff = adrw_ic_space_type_NEX; break;
        case  70 : toff = adrw_ic_space_type_POS; break;
        case  71 : toff = adrw_ic_space_type_POS; xspan = 2; break;
        case  72 : toff = adrw_ic_space_type_POS; xspan = 2; break;
        case  73 : toff = adrw_ic_space_type_NEX; break;
        case  74 : toff = adrw_ic_space_type_NEX; break;
        case  75 : toff = adrw_ic_space_type_BAN; break;
        case  76 : toff = adrw_ic_space_type_SRV; break;
        case  77 : toff = adrw_ic_space_type_BAN; break;
        case  78 : toff = adrw_ic_space_type_ADM; break;
        case  79 : toff = adrw_ic_space_type_POS; break;
        case  80 : toff = adrw_ic_space_type_POS; break;
        case  81 : toff = adrw_ic_space_type_INF; xspan = 2; break;
        case  82 : toff = adrw_ic_space_type_POS; xspan = 2; break;
        case  83 : toff = adrw_ic_space_type_NEX; break;
        case  84 : toff = adrw_ic_space_type_NEX; break;
        case  85 : toff = adrw_ic_space_type_AUD; xspan = 3; break;
        case  86 : toff = adrw_ic_space_type_POS; xspan = 2; break;
        case  87 : toff = adrw_ic_space_type_NEX; break;
        case  88 : toff = adrw_ic_space_type_NEX; break;
        case  89 : toff = adrw_ic_space_type_NEX; break;
        case  90 : toff = adrw_ic_space_type_POS; break;
        case  91 : toff = adrw_ic_space_type_POS; xspan = 2; break;
        case  92 : toff = adrw_ic_space_type_POS; xspan = 2; break;
        case  93 : toff = adrw_ic_space_type_NEX; break;
        case  94 : toff = adrw_ic_space_type_NEX; break;
        case  95 : toff = adrw_ic_space_type_POS; xspan = 2; break;
        case  96 : toff = adrw_ic_space_type_POS; xspan = 2; break;
        case  97 : toff = adrw_ic_space_type_NEX; break;
        case  98 : toff = adrw_ic_space_type_NEX; break;
        case  99 : toff = adrw_ic_space_type_NEX; break;
        case 100 : toff = adrw_ic_space_type_ADM; break;
        case 101 : toff = adrw_ic_space_type_SRV; break;
        case 102 : toff = adrw_ic_space_type_DEP; break;
        default :  demand(FALSE, "invalid office number");
      }
        
    /* Return results: */
    (*toffP) = toff;
    (*xspanP) = xspan;
  }

char **adrw_ic12_get_office_descriptions(void)
  { 
    int32_t max_noff = 103;
    int32_t no = max_noff + 1;
    char **descr = (char **)notnull(malloc(no*sizeof(char *)), "no mem");
    for (int32_t i = 0; i < no; i++) { descr[i] = "Inexistente"; }
    
    descr[  1] = "Escr.Rtorres";
    descr[  2] = "Escr.Helio";
    descr[  3] = "Escr.Rodolfo";
    descr[  4] = "Escr.Lee";
    descr[  5] = "Escr.Visitas";
    descr[  6] = "Escr.Beatriz";
    descr[  7] = "Escr.Ariadne";
    descr[  8] = "Escr.Cid";
    descr[  9] = "Escr.Heloisa";
    descr[ 10] = "Escr.Eduardo";
    descr[ 11] = "Escr.Ducatte";
    descr[ 12] = "Escr.Siome";
    descr[ 13] = "Escr.Cmrubira";
    descr[ 14] = "Escr.Hans";
    descr[ 15] = "Escr.Anamaria";
    descr[ 16] = "Escr.Buzato";
    descr[ 17] = "Qd.Eletr./Dep.";
    descr[ 18] = "Escr.Islene";
    descr[ 19] = "Reun.Limp./Dep.";
    /* 20: inexistente (passagem) */
    /* A Secr.Cursos tem n�mero oficial 20, mas n�o neste programa. */
    descr[ 21] = "Serv.Inform.";
    descr[ 22] = "Banh.Masc.Doc.";
    descr[ 23] = "Escr.Zanoni";
    descr[ 24] = "Banh.Fem.Doc.";
    descr[ 25] = "Escr.Cmbm";
    descr[ 26] = "Escr.Paulo";
    descr[ 27] = "Escr.Cecilia";
    descr[ 28] = "Escr.Rdahab";
    descr[ 29] = "Escr.Rezende";
    descr[ 30] = "Escr.Fkm";
    descr[ 31] = "Escr.Celia";
    descr[ 32] = "Escr.Nfonseca";
    descr[ 33] = "Escr.Stolfi";
    descr[ 34] = "Escr.Julioher";
    descr[ 35] = "Escr.Colabs.";
    descr[ 36] = "Escr.Eliane";
    descr[ 37] = "Escr.Ranido";
    descr[ 38] = "Escr.Wainer";
    descr[ 39] = "Escr.Edmundo";
    descr[ 40] = "Escr.RTC/RTP";
    descr[ 41] = "Escr.Arnaldo";
    descr[ 42] = "Escr.Srigo";
    descr[ 43] = "Escr.Neucimar";
    descr[ 44] = "Escr.Afalcao";
    descr[ 45] = "Escr.Guido";
    descr[ 46] = "Copa(46,48)"; /* N�mero oficial: 46. */
    descr[ 47] = "Escr.RTP/RTC";
    /* 48 inexistente (juntada com 46). */
    /* 49 inexistente. */
    /* 50 inexistente. */
    descr[ 51] = "Reun.Ger.(51,53)";   /* N�mero oficial: 53. */
    descr[ 52] = "Fin./Patr.(52,54)"; /* N�mero oficial: 54. */
    /* 53 inexistente (juntada com 51). */
    /* 54 inexistente (juntada com 52). */
    descr[ 55] = "Secr.Dir.";
    descr[ 56] = "RH/Exped.";
    descr[ 57] = "Escr.ATU";
    descr[ 58] = "Diretor";
    descr[ 59] = "Eng.Manut.";
    descr[ 60] = "LOCO[2]";
    descr[ 61] = "Analistas[1]";
    descr[ 62] = "LAD";
    descr[ 63] = "Analistas[2]";
    descr[ 64] = "LOCO[1]";
    descr[ 65] = "Secr.Deptos.";
    descr[ 66] = "LSC(66,68)"; /* N�mero oficial: 68. */
    descr[ 67] = "Secr.Ext.";
    /* 68 inexistente (juntada com 66). */
    /* 69 inexistente (passagem). */
    descr[ 70] = "Brazil-IP";
    descr[ 71] = "Escr.Dr.[0](71,73)"; /* N�mero oficial: 71. */
    descr[ 72] = "LIV(72,74)";       /* N�mero oficial: 74. */
    /* 73 inexistente (juntada com 71). */
    /* 74 inexistente (juntada com 72). */
    descr[ 75] = "Banh.Fem.Alu.";
    descr[ 76] = "Xerox/Impr.";
    descr[ 77] = "Banh.Masc.Alu.";
    descr[ 78] = "Of.Manut.";
    descr[ 79] = "Atend.PED";
    descr[ 80] = "Reun.Alu.";
    descr[ 81] = "Of.Eletron.(81,83)"; /* N�mero oficial: 83. */
    descr[ 82] = "LAS/LCA(82,84)";     /* N�mero oficial: 84. */
    /* 83 inexistente (juntada com 81). */
    /* 84 inexistente (juntada com 82). */
    descr[ 85] = "Audit.(85,87,89)"; /* N�mero oficial: 85. */
    descr[ 86] = "Escr.Dr.[1](86,88)"; /* N�mero oficial: 86. */
    /* 87 inexistente (juntada com 85). */
    /* 88 inexistente (juntada com 86). */
    /* 89 inexistente (juntada com 85). */
    descr[ 90] = "Escr.Dr[2]";
    descr[ 91] = "LSD(91,93)";       /* N�mero oficial: 91. */
    descr[ 92] = "Escr.Dr.[3](92,94)"; /* N�mero oficial: 92. */
    /* 93 inexistente (juntada com 91). */
    /* 94 inexistente (juntada com 92). */
    descr[ 95] = "LIS(95,97)";       /* N�mero oficial: 93. */
    descr[ 96] = "LRC(96,98)";       /* N�mero oficial: 96. */
    /* 97 inexistente (juntada com 95). */
    /* 98 inexistente (juntada com 96). */
    /* 99 inexistente */
    descr[100] = "Secr.Cursos";       /* N�mero oficial: 20. */
    descr[101] = "Xerox/Impr.";       /* N�mero oficial: 20A. */
    descr[102] = "Almoxarifado";      /* N�mero oficial: 20B. */
    descr[103] = "Coord.Cursos";      /* N�mero oficial: 20C. */
    
    return descr;
  }
 
