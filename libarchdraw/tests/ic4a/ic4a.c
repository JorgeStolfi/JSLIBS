/* Floorlans of IC-4a building */
/* Last edited on 2017-02-26 02:33:36 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <r2.h>
#include <vec.h>
#include <epswr.h>
#include <affirm.h>
#include <jsfile.h>
#include <frgb.h>

#include <archdraw.h>
#include <archdraw_ic.h>

adrw_point_vec_t adrw_ic4a_define_points(void);
void adrw_ic4a_define_points_in_floor(adrw_point_vec_t *P, int kfl, int *np);
int adrw_ic4a_get_office_number_in_floor(int wing, int side, int slot);
int adrw_ic4a_get_first_corner_of_office_in_floor(int noff);
void adrw_ic4a_get_office_corners(int kfl, int wing, int side, int slot, r2_t *p00, r2_t *p01, r2_t *p10, r2_t *p11);
void adrw_ic4a_get_office_type_and_span(int kfl, int noff, adrw_space_type_t *toffP, int *xspanP, int *yspanP);
char **adrw_ic4a_get_office_descriptions(void);

void adrw_ic4a_plot_all(epswr_figure_t *epsf, adrw_building_t *B, int nx, int ny, bool_t show_dots);

void adrw_ic4a_append_building_outline(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]);
void adrw_ic4a_append_pillars(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]);
void adrw_ic4a_append_pillar(adrw_building_t *B, adrw_point_vec_t *P, int kfl, int npil, double wdx, double wdy, adrw_unit_style_t *style[]);
void adrw_ic4a_append_hallways(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]);
void adrw_ic4a_append_offices(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[], char *descr[]);
void adrw_ic4a_append_normal_office
  ( adrw_building_t *B, 
    adrw_point_vec_t *P, 
    int kfl,  /* Floor (0, 1, 2). */
    int wing, /* Wing in floor (0 = west, 1 = east). */
    int side, /* Side of corridor (0 = south, 1 = north). */
    int slot, /* Office slot along one side of a wing (0 = next to stairs). */
    adrw_unit_style_t *style[],
    char *descr[]
  );
void adrw_ic4a_append_power_cabinet(adrw_building_t *B, adrw_point_vec_t *P, int kfl, adrw_unit_style_t *style[], char *descr[]);
void adrw_ic4a_append_stairblock_storage(adrw_building_t *B, adrw_point_vec_t *P, int kfl, adrw_unit_style_t *style[], char *descr[]);
void adrw_ic4a_append_bathroom(adrw_building_t *B, adrw_point_vec_t *P, int kfl, int wing, adrw_unit_style_t *style[], char *descr[]);
void adrw_ic4a_append_auditorium_seating(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]);
void adrw_ic4a_append_entrance_arrow(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]);

#define ADDPOLY(PTS,KFL,LABEL,DESC,NMODS,TYPE,...) \
  do \
    { int v[] = { __VA_ARGS__ , -1 }; \
      int bc = 1000*((KFL)+1); /* Base of corner numbering for this floor. */ \
      int i = 0; while (v[i] >= 0) { v[i] += bc; i++; } \
      adrw_unit_t *rm = adrw_make_poly((LABEL),(DESC),(NMODS),(PTS),v,0.0,(KFL),(TYPE),style[(TYPE)]); \
      adrw_append_unit(B, rm); \
    } \
  while(0)

#define ADDSEATS(PTS,KFL,TYPE,PTA,PTB,SZX,SZY)       \
  do \
    { int bc = 1000*((KFL)+1); /* Base of corner numbering for this floor. */ \
      adrw_append_seats(B,(PTS),bc+(PTA),bc+(PTB),(SZX),(SZY),(KFL),(TYPE),style[(TYPE)]); \
    } \
  while(0)

#define ADDBOX(PTS,KFL,LABEL,DESC,TYPE,CTR,WDX,WDY)      \
  do \
    { int bc = 1000*((KFL)+1); /* Base of corner numbering for this floor. */ \
      adrw_unit_t *rm = adrw_make_box((LABEL),(DESC),(PTS),bc+(CTR),0.0,0.0,0.0,(WDX),(WDY),0.0,(KFL),(TYPE),style[(TYPE)]); \
      adrw_append_unit(B, rm); \
    } \
  while(0)

#define ADDPIPE(PTS,KFL,LABEL,DESC,TYPE,...) \
  do \
    { int v[] = { __VA_ARGS__ , -1 }; \
      int bc = 1000*(KFL+1); /* Base of corner numbering for this floor. */ \
      int i = 0; while (v[i] >= 0) { v[i] += bc; i++; } \
      adrw_unit_t *rm = adrw_make_poly((LABEL),(DESC),(PTS),v,0.0,(KFL),(TYPE),style[(TYPE)]); \
      adrw_append_unit(B, rm); \
    } \
  while(0)

int main (int argc, char **argv)
  {
    int ntypes = adrw_ic_space_type_MAX+1;
    adrw_unit_style_t **style = adrw_ic_define_type_styles();
    char **type_tag = adrw_ic_define_type_tags();
    bool_t *type_is_printable = adrw_ic_define_printable_types();
    bool_t *type_is_movable = adrw_ic_define_movable_types();
    int *color_key_types = adrw_ic_define_color_key_types();
    
    char **descr = adrw_ic4a_get_office_descriptions();
    adrw_point_vec_t P = adrw_ic4a_define_points();
    adrw_building_t *B = adrw_make_building();
    
    adrw_ic4a_append_building_outline(B, &P, style);
    adrw_ic4a_append_hallways(B, &P, style);
    adrw_ic4a_append_offices(B, &P, style, descr);
    adrw_ic4a_append_auditorium_seating(B, &P, style);
    adrw_ic4a_append_pillars(B, &P, style);
    adrw_ic4a_append_entrance_arrow(B, &P, style);
    
    /* Print the point table: */
    { FILE *wr = open_write("out/ic4a-p.txt", TRUE);
      adrw_print_points(wr, &P);
      fclose(wr);
    }

    /* Plot the floorplan: */
    { bool_t eps_fmt = TRUE;
      epswr_figure_t *epsf = pswr_new_stream("out/ic4a-", NULL, eps_fmt, "d", "letter", FALSE, 788.0, 508.0);
      adrw_ic4a_plot_all(epsf, B, 1, 1, FALSE);
      adrw_plot_type_legend(epsf, "tkey", color_key_types, ntypes, 5, type_tag, style);
      adrw_ic_plot_histogram_bars(epsf, B, type_is_movable, type_tag, style);
      /* adrw_ic4a_plot_all(epsf, B, 1, 1, TRUE); */
      pswr_close_stream(epsf);
    }
    
    /* Print the units table: */
    { FILE *wr = open_write("out/ic4a-u.tex", TRUE);
      bool_t TeX = TRUE;
      adrw_print_building(wr, B, type_is_printable, type_tag, TRUE, TRUE, FALSE, TeX);
      fclose(wr);
    }

    return 0;
  }

void adrw_ic4a_plot_all(epswr_figure_t *epsf, adrw_building_t *B, int nx, int ny, bool_t show_dots)
  {
    /* Plot domain in cm, for all three floors together. */
    /* Use nominal width = 7200 cm to get same scale as IC-1+2 plots. */
    double xwid = 6490;
    double ywid = 1190;
    double xmrg = 300 + (7200 - xwid)/2;
    double ymrg = 300;
    double xmin = 00 - xmrg, xmax =          xwid + xmrg;
    double ymin = 00 - ymrg, ymax = 2*1500 + ywid + ymrg;

    int ox, oy;
    for (ox = 0; ox < nx; ox++)
      for (oy = 0; oy < ny; oy++)
        {
          fprintf(stderr, "=== PLOTTING PAGE [%d,%d] OF [%d,%d] ===\n", ox,oy,nx,ny);
          adrw_start_page(epsf, xmin, xmax, ymin, ymax, ox, nx, oy, ny, "Proposta - IC-4a");
          pswr_set_label_font(epsf, "Courier", 8.0);
          adrw_plot_building(epsf, B, show_dots);
        }
  }
  
int adrw_ic4a_get_office_number_in_floor(int wing, int side, int slot)
  { return wing*16 + slot*2 + side; }

int adrw_ic4a_get_first_corner_of_office_in_floor(int noff)
  { return 300 + 4*noff; }
  
adrw_point_vec_t adrw_ic4a_define_points(void)
  {
    adrw_point_vec_t P = adrw_point_vec_new(100);
    int np = 0;
    
    fprintf(stderr, "=== DEFINING POINTS ===\n");
    
    adrw_ic4a_define_points_in_floor(&P, 0, &np);
    adrw_ic4a_define_points_in_floor(&P, 1, &np);
    adrw_ic4a_define_points_in_floor(&P, 2, &np);
    
    adrw_point_vec_trim(&P, np);
    
    return P;
  }
  
void adrw_ic4a_define_points_in_floor(adrw_point_vec_t *P, int kfl, int *np)
  {
    /* All dimensions in centimeters. */

    /* Floor origin: */
    double Xor =   00;            /* X of outer edge of west-most wall. */
    double Yor = +150 + 1500*kfl; /* Y of outer edge of south-most wall. */
    double Zor = 350*kfl;         /* Z of floor. */

    auto void s(int dj, double dX, double dY, int di);
      /* Defines {P[dj]+(dX,dY,350*kfl)} as the coordinates of
        point {P[di]}.  If {P[di]} was defined
        previously, checks whether the definitions agree.
        The point numbers {di,dj} are relative to the floor. */
       
    void s(int dj, double dX, double dY, int di)
      { char *lab = NULL;
        int i = 1000*(kfl+1) + di;
        int j = (dj < 0 ? dj : 1000*(kfl+1) + dj);
        double dZ = (dj < 0 ? Zor : 00);
        asprintf(&lab, "P%04d", i);
        adrw_append_point(lab, i, j, j, j, dX, dY, dZ, P, np); 
      } 
    
    
    /* Assuming exterior N/W and interior walls are 10 cm thick. */
    /* Assuming exterior E/W walls, stairwell W wall and eleator E wall are 20 cm thick. */
    /* Assuming stairwell E wall and other elevator walls are 16 cm thick. */
    /* Assuming interior space is 10.20 m in Y direction, wall to wall. */
    /* Assuming spacing between pillars is 7.20 m, center to center. */

    fprintf(stderr, "--- building outlines ---\n");
    s( -1,   Xor,   Yor, 101);
    s(101, +2880,    00, 102);
    s(102,    00,  -120, 103);
    s(103,  +740,    00, 104);
    s(104,    00,  +120, 105);
    s(105, +2880,    00, 106);
    s(106,    00, +1040, 107);
    s(107, -2880,    00, 108);
    s(108,    00,   +20, 109);
    s(109,  -740,    00, 110);
    s(110,    00,   -20, 111);
    s(111, -2880,    00, 112);
    s(101,    00, +1040, 112); /* Redundant check. */ 
    
    fprintf(stderr, "--- pillar centers ---\n");
    { int fpi = 260;
      int wing; /* Wing (0 = west, 1 = east). */
      int side; /* Side of wing (0 = south, 1 = north). */
      int pill; /* Pillar index (0 = innermost). */
      int kfl;
      for (kfl = 0; kfl < 3; kfl++)
        { for (wing = 0; wing < 2; wing++)
            { for (side = 0; side < 2; side++)
                { for (pill = 0; pill <= 4; pill++)
                    { 
                      int ctr = fpi + 5*wing + 10*side + pill;
                      double xpill, ypill;  /* Pillar center rel to SW wall corner. */
                      if (wing == 0)
                        { xpill = 10 + 720*(4 - pill); }
                      else
                        { xpill = 10 + 720*(5 + pill); }
                      if (side == 0)
                        { ypill = (pill == 0 ? -130 : -20); }
                      else
                        { ypill = +1060; }
                      s(101, xpill, ypill, ctr);
                    }
                }
            }
          s(fpi + 10, +360,  +10, fpi + 20); /* Center pillar, N side. */
          s(fpi + 20,   +5, -500, fpi + 21); /* Internal pillar, center. */
          s(fpi + 10,   00, -500, fpi + 22); /* Internal pillar, W wing. */
          s(fpi + 15,   00, -500, fpi + 23); /* Internal pillar, E wing. */
        } 
    }
    
    fprintf(stderr, "--- hallways ---\n");
    s(112,   +20,  -420, 125);
    s(125, +2515,    00, 126);
    s(126,  +345,    00, 124);
    s(124,    00,   -80, 123);
    s(123,  +740,    00, 122);
    s(122,    00,   +80, 121);
    s(121,  +345,    00, 127);
    if (kfl == 0)
      { 
        s(127, +1800,    00, 120);
        s(107,  -735,  -420, 120); /* Redundant check. */
      }
    else
      { 
        s(127, +2515,    00, 120);
        s(107,   -20,  -420, 120); /* Redundant check. */
      }
    
    s(101,   +20,  +420, 113);
    s(113, +2515,    00, 114);
    s(114,    00,   -85, 115);
    s(115,  +365,    00, 128);
    s(128,  +700,    00, 129);
    s(129,  +365,    00, 117);
    s(117,    00,   +85, 118);
    if (kfl == 0)
      { 
        s(118, +1800,    00, 119);
        s(106,  -735,  +420, 119); /* Redundant check. */ 
      }
    else
      { 
        s(118, +2515,    00, 119);
        s(106,   -20,  +420, 119); /* Redundant check. */ 
      }
    
    fprintf(stderr, "--- stairwell, elevator, power ---\n");
    s(123,   +20,    00, 130);
    s(122,   -20,    00, 147);
    
    s(110,   +20,   -10, 131);
    s(109,   -20,   -10, 136);
    
    s(130,  +700,    00, 147); /* Redundant check. */
    s(136,    00,  -510, 147); /* Redundant check. */
    s(130,    00,  +510, 131); /* Redundant check. */
    
    s(131,  +342,    00, 132);
    s(132,   +16,    00, 135);
    s(135,  +342,    00, 136); /* Redundant check. */
    
    s(132,    00,  -510, 133);
    s(135,    00,  -510, 134);
    s(130,  +342,    00, 133); /* Redundant check. */
    
    s(136,    00,  -115, 137);
    s(137,  -172,    00, 138);
    s(138,    00,  -262, 144);
    
    if (kfl ==0)
      { 
        s(138,   -16,    00, 139);
        s(144,   -16,    00, 143);
      }
    else
      { 
        s(138,   -90,    00, 139);
        s(144,   -90,    00, 143);

        s(139,    00,   -16, 140); /* Front of storage cabinet. */
        s(143,    00,   +16, 142); /* Front of storage cabinet. */

        s(138,   -16,   -16, 161); /* Back of storage cabinet. */
        s(144,   -16,   +16, 141); /* Back of storage cabinet. */
      }
    
    s(130,  +166,    00, 191); /* Stair separator wall. */
    s(191,   +10,    00, 190);
    s(191,    00,  +344, 196);
    s(190,    00,  +344, 193);

    if (kfl ==0)
      { 
        
        s(191,  -166,   +10, 201);
        s(191,    00,   +10, 192);
        s(196,  -166,    00, 198);
        s(193,  +166,    00, 200);
        
        s(198,    00,   +10, 197);
        s(196,    00,   +10, 195);
        s(193,    00,   +10, 194);
        s(200,    00,   +10, 199);
        
      }
    else if (kfl == 2)
      {
        s(190,    00,   +10, 202);
        s(190,  +166,   +10, 203);
        s(190,  +166,    00, 133); /* Redundant check. */
      }
      
    s(138,    00,   -40, 162); /* Back of power/pipes cabinet */
    s(137,    00,   -40, 163); /* Back of power/pipes cabinet */
    
    s(162,    00,   -16, 145); /* Back of elevator. */
    s(163,    00,   -16, 146); /* Back of elevator. */

    s(145,    00,  -194, 164); /* Front of elevator. */
    s(146,    00,  -194, 165); /* Front of elevator. */
    
    s(164,    00,   -10, 166); /* Front of elevator doorwall. */
    s(165,    00,   -10, 167); /* Front of elevator doorwall. */

    fprintf(stderr, "--- office corners ---\n");
    int wing; /* Wing in floor (0 = west, 1 = east). */
    int side; /* Side of corridor (0 = south, 1 = north). */
    int slot; /* Office slot along one side of a wing (0 = next to stairs). */
    for (side = 0; side < 2; side++)
      { /* Compute basic Y cords of office: */
        double yot, yin; /* Inner and outer Y of office, excl. walls. */
        if (side == 0)
          { yot =   +10; yin = yot + 400; }
        else
          { yot = +1030; yin = yot - 400; }
        for (wing = 0; wing < 2; wing++)
          { for (slot = 0; slot < 8; slot++)
              { int noff = adrw_ic4a_get_office_number_in_floor(wing, side, slot);
                int fco = adrw_ic4a_get_first_corner_of_office_in_floor(noff);
                fprintf(stderr, "office %d corners %d..%d\n", 100*kfl + noff, 1000*kfl + fco, 1000*kfl + fco + 3);
                
                /* Get X and Y cords of office (excluding walls) rel to floor origin: */
                r2_t p00;  /* Outer proximal corner of office: */ 
                r2_t p01;  /* Inner proximal corner of office: */ 
                r2_t p10;  /* Outer distal corner of office: */   
                r2_t p11;  /* Inner distal corner of office: */   
                adrw_ic4a_get_office_corners(kfl, wing, side, slot, &p00, &p01, &p10, &p11);

                /* Set the corners: */
                s(101, p00.c[0], p00.c[1], fco+0);
                s(101, p01.c[0], p01.c[1], fco+1);
                s(101, p10.c[0], p10.c[1], fco+2);
                s(101, p11.c[0], p11.c[1], fco+3);
              }
          }
      }

    fprintf(stderr, "--- west bathroom corners ---\n");
    s(128,   00,  -10, 151);
    s(151, +345,   00, 150);
    s(150,   00, -435, 152);
    s(151,   00, -435, 153);
    s(103,  +20,  +10, 153); /* Redundant check. */
    
    fprintf(stderr, "--- east bathroom corners ---\n");
    s(129,   00,  -10, 156);
    s(156, -345,   00, 155);
    s(155,   00, -435, 157);
    s(156,   00, -435, 158);
    s(104,  -20,  +10, 158); /* Redundant check. */
    
    if (kfl ==0)
      { 
        fprintf(stderr, "--- auditorium seating ---\n");
        int auc = adrw_ic4a_get_first_corner_of_office_in_floor(30); /* Outer proximal corner */
        s(auc+2,  -10,   00, 181);
        s(181,   -600,   00, 180);
        s(181,     00, +720, 182);
        s(182,   -600,   00, 183);
        s(180,     00, +720, 183); /* Redundant check. */
        
        s(180,    -60,  +80, 184); /* Lone seat in back. */
      }
              
    show_holes(P, 1000*(kfl+1),1000*(kfl+2)-1);
    
    if (kfl ==0)
      { 
        fprintf(stderr, "--- entrance arrow ---\n");
        s(135,   +77,  -170, 220);
        s(220,   +50,   +50, 221);
        s(221,   -30,    00, 222);
        s(222,    00,  +210, 223);
        s(223,   -40,    00, 224);
        s(224,    00,  -210, 225);
        s(225,   -30,    00, 226);
        s(220,   -50,   +50, 226); /* Redundant check. */
      }
  }

void adrw_ic4a_get_office_corners(int kfl, int wing, int side, int slot, r2_t *p00, r2_t *p01, r2_t *p10, r2_t *p11)
  {
    int noff = adrw_ic4a_get_office_number_in_floor(wing, side, slot);
                
    /* X and Y cords of office (excluding walls) rel to floor origin: */
    double xp;   /* Proximal X coord of office: */
    double xd;   /* Distal X coord of office: */
    double yop;  /* Y coord of outer proximal corner of office: */
    double yip;  /* Y coord of inner proximal corner of office: */
    double yod;  /* Y coord of outer distal corner of office: */
    double yid;  /* Y coord of inner distal corner of office: */

    /* Compute coords for typical office, rel to SW corner: */
    if (wing == 0)
      { xp = 10 + 360*(8 - slot) - 5; 
        xd = xp - 350;
      }
    else
      { xp = 10 + 360*(10 + slot) + 5; 
        xd = xp + 350;
      }
    if (side == 0)
      { yod = yop = +10;
        yid = yip = yod + 400;
      }
    else
      { yod = yop = +1030;
        yid = yip = yod - 400;
      }

    /* Apply adjustments for end offices (thicker walls, bent corridor): */
    if (noff == 0)
      { yip += -83; xp += +5; }
    else if (noff == 1)
      { yip += -80; xp += -5; }
    else if (noff == 16)
      { yip += -83; xp += -5; }
    else if (noff == 17)
      { yip += -80; xp += +5; }
    else if ((noff == 14) || (noff == 15))
      { xd += +5; }
    else if ((noff == 30) || (noff == 31))
      { xd += -5; }

    /* Return results: */
    (*p00) = (r2_t){{ xp, yop }};
    (*p01) = (r2_t){{ xp, yip }};
    (*p10) = (r2_t){{ xd, yod }};
    (*p11) = (r2_t){{ xd, yid }};
  }

void adrw_ic4a_append_building_outline(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[])
  {
    int kfl;
    for (kfl = 0; kfl < 3; kfl++)
      { 
        fprintf(stderr, "--- building outline - floor %d ---\n", kfl);
        ADDPOLY
          ( P, kfl, "Ext", "Building outline", 0.0, adrw_ic_space_type_ETC,
            101, 102, 103, 104, 105, 
            106, 107, 108, 109, 110, 
            111, 112, 101
          );
      }
  }
  
void adrw_ic4a_append_pillars(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[])
  {
    fprintf(stderr, "--- pillars ---\n");
    int wing; /* Wing (0 = south, 1 = north). */
    int side; /* Side of wing (0 = south, 1 = north). */
    int pill; /* Pillar index (0 = westernmost). */
    int kfl;
    for (kfl = 0; kfl < 3; kfl++)
      { for (wing = 0; wing < 2; wing++)
          { for (side = 0; side < 2; side++)
              { for (pill = 0; pill <= 4; pill++)
                  { 
                    int npil = 5*wing + 10*side + pill;
                    adrw_ic4a_append_pillar(B, P, kfl, npil, 40, 60, style);
                  }
              }
          }
        adrw_ic4a_append_pillar(B, P, kfl, 20, 40, 40, style);
        adrw_ic4a_append_pillar(B, P, kfl, 21, 30, 50, style);
        adrw_ic4a_append_pillar(B, P, kfl, 22, 40, 40, style);
        adrw_ic4a_append_pillar(B, P, kfl, 23, 40, 40, style);
      }
  }

void adrw_ic4a_append_pillar(adrw_building_t *B, adrw_point_vec_t *P, int kfl, int npil, double wdx, double wdy, adrw_unit_style_t *style[])
  {
    int fpi = 260; /* Index of first pillar center point in floor. */
    int ctr = fpi + npil;
    char *lpil = NULL;
    asprintf(&lpil, "P%02d", npil);
    char *dpil = NULL;
    asprintf(&dpil, "Pilar %02d", npil);
    ADDBOX
      ( P, kfl, lpil, dpil, adrw_ic_space_type_PIL,
        ctr, wdx, wdy
      );
  }

void adrw_ic4a_append_hallways(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[])
  {
    int kfl;
    
    fprintf(stderr, "--- hallways and stairs - floor 0 ---\n");
    ADDPOLY
      ( P, 0, "Cor", "Corr./Escad.", 0.0, adrw_ic_space_type_HAL,
        113, 114, 128, 129, 118, 
        119, 120, 127, 122, 
        147, 167, 166, 144, 143,
        139, 138, 137, 136,
        135, 134, 130,
        123, 126, 125, 113
      );
    
    fprintf(stderr, "--- hallways and stairs - floor 1 ---\n");
    ADDPOLY
      ( P, 1, "Cor", "Corr./Escad.", 0.0, adrw_ic_space_type_HAL,
        113, 114, 128, 129, 118, 
        119, 120, 127, 122, 
        147, 167, 166, 144, 143,
        139, 138, 137, 136,
        135, 134, 133, 
        191, 196, 193, 190, 
        133, 132, 131, 130,
        123, 126, 125, 113
      );
    
    fprintf(stderr, "--- hallways and stairs - floor 2 ---\n");
    ADDPOLY
      ( P, 2, "Cor", "Corr./Escad.", 0.0, adrw_ic_space_type_HAL,
        113, 114, 128, 129, 118, 
        119, 120, 127, 122, 
        147, 167, 166, 144, 143,
        139, 138, 137, 136,
        135, 134, 133, 
        191, 196, 193, 202, 203, 
        132, 131, 130,
        123, 126, 125, 113
      );
    
    for (kfl = 0; kfl < 3; kfl++)
      { 
        int bc = 1000*(kfl+1);
        fprintf(stderr, "--- elevator - floor %d ---\n", kfl);
        /* Get the numbers of the four corners. */
        int c00 = 145;  /* Xp, Yot. */
        int c01 = 164;  /* Xp, Yin. */
        int c10 = 146;  /* Xd, Yot. */
        int c11 = 165;  /* Xd, Yin. */
        fprintf
          ( stderr, "  corners %d %d %d %d\n", 
            bc+c00, bc+c01, bc+c10, bc+c11
          );
        ADDPOLY
          ( P, kfl, "Ele", "Elevador", 0.0, adrw_ic_space_type_HAL,
            c00, c01, c11, c10, c00
          ); 
      }
  }

void adrw_ic4a_append_offices(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[], char *descr[])
  {
    int kfl;
    for (kfl = 0; kfl < 3; kfl++)
      { 
        int wing; /* Wing in floor (0 = west, 1 = east). */
        int side; /* Side of corridor (0 = south, 1 = north). */
        int slot; /* Office slot along one side of a wing (0 = next to stairs). */
        for (side = 0; side < 2; side++)
          { for (wing = 0; wing < 2; wing++)
              { for (slot = 0; slot < 8; slot++)
                  { adrw_ic4a_append_normal_office(B, P, kfl, wing, side, slot, style, descr); }
              }
          }
        adrw_ic4a_append_power_cabinet(B, P, kfl, style, descr);
        adrw_ic4a_append_stairblock_storage(B, P, kfl, style, descr);
        adrw_ic4a_append_bathroom(B, P, kfl, 0, style, descr);
        adrw_ic4a_append_bathroom(B, P, kfl, 1, style, descr);
      }
  }

void adrw_ic4a_append_normal_office
  ( adrw_building_t *B, 
    adrw_point_vec_t *P, 
    int kfl,  /* Floor (0, 1, 2). */
    int wing, /* Wing in floor (0 = west, 1 = east). */
    int side, /* Side of corridor (0 = south, 1 = north). */
    int slot, /* Office slot along one side of a wing (0 = next to stairs). */
    adrw_unit_style_t *style[],
    char *descr[]
  )
  {
    int noff = adrw_ic4a_get_office_number_in_floor(wing, side, slot); /* Office number in floor. */
    int boff = 100*(kfl + 1); /* Start of office numbering in this floor. */
    char *loff = NULL;
    asprintf(&loff, "%03d", boff+noff);
    char *doff = descr[boff+noff];
    int bc = 1000*(kfl + 1);  /* Start of corner numbering in this floor. */
    /*  Decide office type {toff} and span {xspan,yspan}: */
    adrw_space_type_t toff; /* Type of office. */
    int xspan; /* How many modules it spans in X direction (1 to 8). */
    int yspan; /* How many modules it spans in Y direction (1 or 2). */
    adrw_ic4a_get_office_type_and_span(kfl, noff, &toff, &xspan, &yspan);
    double modules = xspan*yspan;
    if (toff == adrw_ic_space_type_NEX)
      { /* Skip it. */
        fprintf(stderr, "  office %d does not exist - skipped\n", boff+noff);
      }
    else if (xspan == 1)
      { /* Get the numbers of the four corners. */
        int fco = adrw_ic4a_get_first_corner_of_office_in_floor(noff);
        int c00 = fco + 0;            /* Xin, Yot. */
        int c01 = c00 + (yspan == 1 ? 1 : 4);  /* Xin, Yin. */
        int c10 = c00 + 2;                     /* Xot, Yot. */
        int c11 = c01 + 2;                     /* Xot, Yin. */
        fprintf
          ( stderr, "  office %d - corners %d %d %d %d\n", 
            boff+noff, bc+c00, bc+c01, bc+c10, bc+c11
          );
        ADDPOLY
          ( P, kfl, loff, doff, modules, toff,
            c00, c01, c11, c10, c00
          ); 
      }
    else if (xspan == 2) 
      { /* Get the numbers of the eight corners. */
        int fco = adrw_ic4a_get_first_corner_of_office_in_floor(noff);
        int c00 = fco + 0;            /* Xin, Yot. */
        int c01 = c00 + (yspan == 1 ? 1 : 4);  /* Xin, Yin. */
        int c10 = c00 + 2;                     /* Xm1, Yo1. */
        int c11 = c01 + 2;                     /* Xm1, Yi1. */
        int c20 = c00 + 8;                     /* Xm2, Yo2. */
        int c21 = c01 + 8;                     /* Xm2, Yi2. */
        int c30 = c20 + 2;                     /* Xot, Yot. */
        int c31 = c21 + 2;                     /* Xot, Yin. */
        fprintf
          ( stderr, "  office %d - corners %d %d  %d %d  %d %d  %d %d\n", 
            boff+noff, 
            bc+c00, bc+c01, bc+c10, bc+c11, 
            bc+c20, bc+c21, bc+c30, bc+c31
          );
        ADDPOLY
          ( P, kfl, loff, doff, modules, toff,
            c01, c11, c21, c31,
            c30, c20, c10, c00,
            c01 
          );
      }
    else if (xspan == 3)
      { /* Get the numbers of the 12 corners. */
        int fco = adrw_ic4a_get_first_corner_of_office_in_floor(noff);
        int c00 = fco + 0;           /* Xin, Yot. */
        int c01 = c00 + (yspan == 1 ? 1 : 4); /* Xin, Yin. */
        int c10 = c00 + 2;                    /* Xm1, Yo1. */
        int c11 = c01 + 2;                    /* Xm1, Yi1. */
        int c20 = c00 + 8;                    /* Xm2, Yo2. */
        int c21 = c01 + 8;                    /* Xm2, Yi2. */
        int c30 = c20 + 2;                    /* Xm3, Yo3. */
        int c31 = c21 + 2;                    /* Xm3, Yi3. */
        int c40 = c20 + 8;                    /* Xm4, Yo4. */
        int c41 = c21 + 8;                    /* Xm4, Yi4. */
        int c50 = c40 + 2;                    /* Xot, Yot. */
        int c51 = c41 + 2;                    /* Xot, Yin. */
        fprintf
          ( stderr, "  office %d - corners %d %d  %d %d  %d %d  %d %d  %d %d  %d %d\n", 
            boff+noff, 
            bc+c00, bc+c01, bc+c10, bc+c11, 
            bc+c20, bc+c21, bc+c30, bc+c31, 
            bc+c40, bc+c41, bc+c50, bc+c51 
          );
        ADDPOLY
          ( P, kfl, loff, doff, modules, toff,
            c01, c11, c21, c31, c41, c51, 
            c50, c40, c30, c20, c10, c00,
            c01
          );
      }
    else
      { affirm(FALSE, "invalid xspan/yspan combination"); }
  }
          
void adrw_ic4a_append_stairblock_storage(adrw_building_t *B, adrw_point_vec_t *P, int kfl, adrw_unit_style_t *style[], char *descr[])
  { 
    if (kfl > 0)
      { fprintf(stderr, "--- storage cabinet next to elevator - floor %d ---\n", kfl);
        int noff = 32; /* Office number in floor. */
        int boff = 100*(kfl + 1); /* Start of office numbering in this floor. */
        char *loff = NULL;
        asprintf(&loff, "%03d", boff+noff);
        char *doff = descr[boff+noff];
        int bc = 1000*(kfl + 1);  /* Start of corner numbering in this floor. */
        double modules = 0.1;     /* Only {1.7 m^2} */
        /* Get the numbers of the four corners. */
        int c00 = 161;  /* Xd, Yot. */
        int c01 = 141;  /* Xd, Yin. */
        int c10 = 140;  /* Xp, Yot. */
        int c11 = 142;  /* Xp, Yin. */
        fprintf
          ( stderr, "  corners %d %d %d %d\n", 
            bc+c00, bc+c01, bc+c10, bc+c11
          );
        ADDPOLY
          ( P, kfl, loff, doff, modules, adrw_ic_space_type_DEP,
            c00, c01, c11, c10, c00
          ); 
      }
      
    if (kfl == 0)
      { fprintf(stderr, "--- outer storage cabinet under stairs - floor %d ---\n", kfl);
        int noff = 37; /* Office number in floor. */
        int boff = 100*(kfl + 1); /* Start of office numbering in this floor. */
        char *loff = NULL;
        asprintf(&loff, "%03d", boff+noff);
        char *doff = descr[boff+noff];
        int bc = 1000*(kfl + 1);  /* Start of corner numbering in this floor. */
        double modules = 0.4;     /* Only {5.3 m^2} */
        /* Get the numbers of the four corners. */
        int c00 = 131;  /* Xd, Yot. */
        int c01 = 197;  /* Xd, Yin. */
        int c10 = 132;  /* Xp, Yot. */
        int c11 = 199;  /* Xp, Yin. */
        fprintf
          ( stderr, "  corners %d %d %d %d\n", 
            bc+c00, bc+c01, bc+c10, bc+c11
          );
        ADDPOLY
          ( P, kfl, loff, doff, modules, adrw_ic_space_type_DEP,
            c00, c01, c11, c10, c00
          ); 
      }
      
    if (kfl == 0)
      { fprintf(stderr, "--- inner storage cabinet under stairs - floor %d ---\n", kfl);
        int noff = 36; /* Office number in floor. */
        int boff = 100*(kfl + 1); /* Start of office numbering in this floor. */
        char *loff = NULL;
        asprintf(&loff, "%03d", boff+noff);
        char *doff = descr[boff+noff];
        int bc = 1000*(kfl + 1);  /* Start of corner numbering in this floor. */
        double modules = 0.4;     /* Only {5.5 m^2} */
        /* Get the numbers of the four corners. */
        int c00 = 198;  /* Xd, Yot. */
        int c01 = 201;  /* Xd, Yin. */
        int c10 = 196;  /* Xp, Yot. */
        int c11 = 192;  /* Xp, Yin. */
        fprintf
          ( stderr, "  corners %d %d %d %d\n", 
            bc+c00, bc+c01, bc+c10, bc+c11
          );
        ADDPOLY
          ( P, kfl, loff, doff, modules, adrw_ic_space_type_DEP,
            c00, c01, c11, c10, c00
          ); 
      }
  }

void adrw_ic4a_append_power_cabinet(adrw_building_t *B, adrw_point_vec_t *P, int kfl, adrw_unit_style_t *style[], char *descr[])
  { fprintf(stderr, "--- power cabinet - floor %d ", kfl);
    int noff = 33; /* Office number in floor. */
    int boff = 100*(kfl + 1); /* Start of office numbering in this floor. */
    char *loff = NULL;
    asprintf(&loff, "%03d", boff+noff);
    char *doff = descr[boff+noff];
    int bc = 1000*(kfl + 1);  /* Start of corner numbering in this floor. */
    double modules = 0.0;     /* Only {0.6 m^2} */
    /* Get the numbers of the four corners. */
    int c00 = 138;  /* Xp, Yot. */
    int c01 = 162;  /* Xp, Yin. */
    int c10 = 137;  /* Xd, Yot. */
    int c11 = 163;  /* Xd, Yin. */
    fprintf
      ( stderr, " - corners %d %d %d %d---\n", 
        bc+c00, bc+c01, bc+c10, bc+c11
      );
    ADDPOLY
      ( P, kfl, loff, doff, modules, adrw_ic_space_type_SRV,
        c00, c01, c11, c10, c00
      ); 
  }
          
void adrw_ic4a_append_bathroom(adrw_building_t *B, adrw_point_vec_t *P, int kfl, int wing, adrw_unit_style_t *style[], char *descr[])
  { fprintf(stderr, "--- bathroom - floor %d - wing %d ---\n", kfl, wing);
    int noff = 34 + wing; /* Office number in floor. */
    int boff = 100*(kfl + 1); /* Start of office numbering in this floor. */
    char *loff = NULL;
    asprintf(&loff, "%03d", boff+noff);
    char *doff = descr[boff+noff];
    int bc = 1000*(kfl + 1);  /* Start of corner numbering in this floor. */
    double modules = 1.0;     /* Only {1.7 m^2} */
    /* Get the numbers of the four corners. */
    int fco = 150 + 5*wing;
    int c01 = fco + 0;  /* Xp, Yin. */
    int c11 = fco + 1;  /* Xd, Yin. */
    int c00 = fco + 2;  /* Xp, Yot. */
    int c10 = fco + 3;  /* Xd, Yot. */
    fprintf
      ( stderr, "  corners %d %d %d %d\n", 
        bc+c00, bc+c01, bc+c10, bc+c11
      );
    ADDPOLY
      ( P, kfl, loff, doff, modules, adrw_ic_space_type_BAN,
        c00, c01, c11, c10, c00
      ); 
  }

void adrw_ic4a_append_auditorium_seating(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[]) 
  {
    int kfl = 0;
    int noff = 28; /* Office number, spans 4 */
    int bank; /* Bank of seats (0 for now). */
    for (bank = 0; bank < 1; bank++)
      { /* Get the numbers of the four corners in cyclic order. */
        int fco = 180 + 4*bank; /* Index of first corner. */
        int c0 = fco + 0;
        int c2 = fco + 2;
        fprintf(stderr, "--- auditorium %d seats - bank %d - corners %d %d\n", noff, bank, c0, c2);
        ADDSEATS
          ( P, kfl, adrw_ic_space_type_SIT,
            c0, c2, 60, 80
          ); 
      }
    { int c0 = 180;
      int c2 = 184;
      fprintf(stderr, "--- auditorium %d seats lone seat in back - corners %d %d\n", noff, c0, c2);
      ADDSEATS
        ( P, kfl, adrw_ic_space_type_SIT,
          c0, c2, 60, 80
        ); 
    }
  }

void adrw_ic4a_append_entrance_arrow(adrw_building_t *B, adrw_point_vec_t *P, adrw_unit_style_t *style[])
  {
    /* Get the numbers of the 7 corners in cyclic order. */
    int fco = 220; /* Index of arrow tip. */
    int c0 = fco + 0;
    int c1 = fco + 1;
    int c2 = fco + 2;
    int c3 = fco + 3;
    int c4 = fco + 4;
    int c5 = fco + 5;
    int c6 = fco + 6;
    fprintf(stderr, "--- entrance arrow - corners %d %d %d %d %d %d %d\n", c0, c1, c2, c3, c4, c5, c6);
    ADDPOLY
      ( P, 0, "Ent", "Entrada", 0.0, adrw_ic_space_type_ARR,
        c0, c1, c2, c3, c4, c5, c6, c0
      ); 
  }

void adrw_ic4a_get_office_type_and_span(int kfl, int noff, adrw_space_type_t *toffP, int *xspanP, int *yspanP)
  {
    adrw_space_type_t toff = adrw_ic_space_type_ETC; /* Type of office. */
    int xspan = 1; /* How many modules it spans in X direction (1 to 8). */
    int yspan = 1; /* How many modules it spans in Y direction (1 or 2). */
    if (kfl == 0)
      { 
        switch(noff)
          { 
            case  0: toff = adrw_ic_space_type_DEP; break;
            case  1: toff = adrw_ic_space_type_ADM; break;
            case  2: toff = adrw_ic_space_type_REU; break;
            case  3: toff = adrw_ic_space_type_ADM; xspan = 2; break;
            case  4: toff = adrw_ic_space_type_SRV; break;
            case  5: toff = adrw_ic_space_type_NEX; break;
            case  6: toff = adrw_ic_space_type_DEP; break;
            case  7: toff = adrw_ic_space_type_ADM; xspan = 3; break;
            case  8: toff = adrw_ic_space_type_ADM; break;
            case  9: toff = adrw_ic_space_type_NEX; break;
            case 10: toff = adrw_ic_space_type_ADM; break;
            case 11: toff = adrw_ic_space_type_NEX; break;
            case 12: toff = adrw_ic_space_type_DEP; xspan = 2; break;
            case 13: toff = adrw_ic_space_type_INF; break;
            case 14: toff = adrw_ic_space_type_NEX; break;
            case 15: toff = adrw_ic_space_type_INF; break;
            
            case 16: toff = adrw_ic_space_type_SRV; xspan = 2; break;
            case 17: toff = adrw_ic_space_type_ADM; xspan = 2; break;
            case 18: toff = adrw_ic_space_type_NEX; break;
            case 19: toff = adrw_ic_space_type_NEX; break;
            case 20: toff = adrw_ic_space_type_DEP; break;
            case 21: toff = adrw_ic_space_type_ADM; break;
            case 22: toff = adrw_ic_space_type_ADM; break;
            case 23: toff = adrw_ic_space_type_ADM; break;
            case 24: toff = adrw_ic_space_type_ADM; break;
            case 25: toff = adrw_ic_space_type_ADM; break;
            case 26: toff = adrw_ic_space_type_ADM; break;
            case 27: toff = adrw_ic_space_type_ADM; break;
            case 28: toff = adrw_ic_space_type_AUD; xspan = 2; yspan = 2; break;
            case 29: toff = adrw_ic_space_type_NEX; break;
            case 30: toff = adrw_ic_space_type_NEX; break;
            case 31: toff = adrw_ic_space_type_NEX; break;
            default: affirm(FALSE, "invalid noff");
          }
      }
    else
      { if (noff == 0)
          { toff = adrw_ic_space_type_DEP; }
        else if (noff == 16)
          { toff = adrw_ic_space_type_SRV; }
        else if ((noff == 1) || (noff == 17))
          { toff = adrw_ic_space_type_REU; xspan = 2; }
        else if ((noff == 3) || (noff == 19))
          { toff = adrw_ic_space_type_NEX; }
        else
          { toff = adrw_ic_space_type_DOC; }
      }

    /* Return results: */
    (*toffP) = toff;
    (*xspanP) = xspan;
    (*yspanP) = yspan;
  }

char **adrw_ic4a_get_office_descriptions(void)
  { 
    int max_noff = 399;
    int no = max_noff + 1;
    char **descr = (char **)notnull(malloc(no*sizeof(char *)), "no mem");
    int i;
    for (i = 0; i < no; i++) { descr[i] = "Inexistente"; }
    
    /* TÉRREO */
    descr[100] = "Limpeza[0]";
    descr[101] = "Secr.Ext.";
    descr[102] = "Laz.Limp.";
    descr[103] = "RH/Exped(103,105)";
    descr[104] = "Xerox/Impr.[0]";
    /* 105 inexistente (anexado ao 103). */
    descr[106] = "Arq.Morto";
    descr[107] = "Fin./Patr(107,109,111)";
    descr[108] = "Eng.Manut";
    /* 109 inexistente (anexado ao 107). */
    descr[110] = "Of.Manut.";
    /* 111 inexistente (anexado ao 107). */
    descr[112] = "Almoxarif.";
    descr[113] = "Analistas";
    /* 114 inexistente (anexado ao 112). */
    descr[115] = "Of.Eletron.";
    
    descr[116] = "Serv.Inf.(116,118)";
    descr[117] = "Secr.Cursos(117,119)";
    /* 118 inexistente (anexado ao 116). */
    /* 119 inexistente (anexado ao 117). */
    descr[120] = "Arq.Cursos";
    descr[121] = "Coord.Cursos";
    descr[122] = "Secr.Deptos.";
    descr[123] = "Dir.Assoc.";
    descr[124] = "Secr.Diretor";
    descr[125] = "Reun.Diretor";
    descr[126] = "Escr.ATU";
    descr[127] = "Diretor";
    descr[128] = "Audit.(128--131)";
    /* 129 inexistente (anexado ao 128). */
    /* 130 inexistente (anexado ao 128). */
    /* 131 inexistente (anexado ao 128). */
    /* 132 inexistente. */
    descr[133] = "Qd.Eletr.[0]";
    descr[134] = "Banh.Masc.[0]";
    descr[135] = "Banh.Fem.[0]";
    descr[136] = "Limpeza[1]";
    descr[137] = "Depósito[0]";
    
    /* PRMEIRO ANDAR */
    descr[200] = "Limpeza[2]";
    descr[201] = "Reun.Doc.[1]";
    descr[202] = "Escr.Doc.";
    /* 203 inexistente (anexado ao 201). */
    for (i = 204; i <= 215; i++) { descr[i] = "Escr.Doc."; }
    
    descr[216] = "Impressora[1]";
    descr[217] = "Copa/Lazer(217,219)";
    descr[218] = "Escr.Doc.";
    /* 219 inexistente (anexado ao 217). */
    for (i = 220; i <= 231; i++) { descr[i] = "Escr.Doc."; }
    descr[232] = "Depósito[2]";
    descr[233] = "Qd.Eletr.[1]";
    descr[234] = "Banh.Masc.[1]";
    descr[235] = "Banh.Fem.[1]";
    
    /* SEGUNDO ANDAR */
    descr[300] = "Limpeza[3]";
    descr[301] = "Reun.Doc.[2]";
    descr[302] = "Escr.Doc.";
    /* 303 inexistente (anexado ao 301). */
    for (i = 304; i <= 315; i++) { descr[i] = "Escr.Doc."; }
    
    descr[316] = "Impressora[2]";
    descr[317] = "Reun.Doc.[3](317,319)";
    descr[318] = "Escr.Doc.";
    /* 319 inexistente (anexado ao 317). */
    for (i = 320; i <= 331; i++) { descr[i] = "Escr.Doc."; }
    descr[332] = "Deposito[2]";
    descr[333] = "Qd.Eletr.[2]";
    descr[334] = "Banh.Masc.[2]";
    descr[335] = "Banh.Fem.[2]";
    
    return descr;
  }
 
