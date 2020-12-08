/* Planta basica da casa na Boaretto da Silva 113 */
/* Last edited on 2012-12-07 20:54:43 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <r2.h>
#include <vec.h>
#include <epswr.h>
#include <affirm.h>
#include <jsfile.h>

#include <archdraw.h>

/* !!! First build the units then draw them. !!! */

adrw_point_vec_t define_points(void);

void plot_all(char *prefix, adrw_point_vec_t *P, int nx, int ny, bool_t show_dots);

void draw_perimeter_wall(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);
void draw_old_sidewalk(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);
void draw_sidewalk_cuts(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);
void draw_trimmed_sidewalk(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);
void draw_sidewalk_additions(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);
void draw_building(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);
void draw_drain_boxes(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);
void draw_drain_pipes(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots);

#define PSHAPE(LABEL,DESC,R,G,B,SHOW_DOTS,...)  \
  do \
    { int v[] = { __VA_ARGS__ , -1 }; \
      frgb_t color = (frgb_t){{(float)(R),(float)(G),(float)(B)}}; \
      frgb_t black = (frgb_t){{0,0,0}}; \
      adrw_unit_style_t *st = adrw_make_unit_style(&color,&black,0.10,NULL,0.0); \
      int kfloor = 0, type = 0; \
      adrw_unit_t *rm = adrw_make_poly((LABEL),(DESC),0.0,P,v,0.0,kfloor,type,st); \
      adrw_plot_unit(epsf, rm, (SHOW_DOTS)); \
      free(st); free(rm->pt.e); free(rm); \
    } \
  while(0)

#define DRBOX(LABEL,DESC,R,G,B,SHOW_DOTS,PT) \
  do \
    { frgb_t color = (frgb_t){{(float)(R),(float)(G),(float)(B)}}; \
      frgb_t black = (frgb_t){{0,0,0}}; \
      adrw_unit_style_t *st = adrw_make_unit_style(&color,&black,0.10,NULL,0.0); \
      int kfloor = 0, type = 0; \
      adrw_unit_t *rm = adrw_make_box((LABEL),(DESC),P,(PT),0.0,0.0,0.0,28,28,0.0,kfloor,type,st); \
      adrw_plot_unit(epsf, rm, (SHOW_DOTS)); \
      free(st); free(rm->pt.e); free(rm); \
    } \
  while(0)

#define DRPIPE(LABEL,DESC,R,G,B,SHOW_DOTS,PTA,PTB) \
  do \
    { int v[] = { (PTA), (PTB), -1 };        \
      frgb_t color = (frgb_t){{(float)(R),(float)(G),(float)(B)}}; \
      adrw_unit_style_t *st = adrw_make_unit_style(NULL,&color,0.50,NULL,0.0); \
      int kfloor = 0, type = 0; \
      adrw_unit_t *rm = adrw_make_poly((LABEL),(DESC),0.0,P,v,0.0,kfloor,type,st); \
      adrw_plot_unit(epsf, rm, (SHOW_DOTS)); \
      free(st); free(rm->pt.e); free(rm); \
    } \
  while(0)

int main (int argc, char **argv)
  {

    adrw_point_vec_t P = define_points();
    
    /* Print the point table: */
    FILE *wr = open_write("out/boar_p.txt", TRUE);
    adrw_print_points(wr, &P);
    fclose(wr);

    /* Plot the floorplan: */
    plot_all("out/boar_A", &P, 1, 1, FALSE);
    plot_all("out/boar_B", &P, 2, 3, FALSE);
    plot_all("out/boar_C", &P, 1, 1, TRUE);
    
    return 0;
  }

void plot_all(char *prefix, adrw_point_vec_t *P, int nx, int ny, bool_t show_dots)
  {
    double hSize = 460;
    double vSize = 640;
    
    double xmin = -100, xmax = +1700;
    double ymin = -100, ymax = +3300;
    
    int ox, oy;
    for (ox = 0; ox < nx; ox++)
      for (oy = 0; oy < ny; oy++)
        {
          char *fname = NULL;
          asprintf(&fname, "%s_%03d_%03d.eps", prefix, ox, oy);
          FILE *wr = open_write(fname, TRUE);
          bool_t verbose = TRUE;
          double hvMarg = 4.0;
          epswr_figure_t *epsf = epswr_new_figure(wr, hSize, vSize, hvMarg, hvMarg, hvMarg, hvMarg, verbose);

          adrw_start_page(epsf, xmin, xmax, ymin, ymax, ox, nx, oy, ny, "Situacao atual");
          draw_old_sidewalk(epsf, P, show_dots);
          draw_building(epsf, P, show_dots);
          draw_perimeter_wall(epsf, P, show_dots);

          adrw_start_page(epsf, xmin, xmax, ymin, ymax, ox, nx, oy, ny, "Remocao de ladrilhos");
          draw_old_sidewalk(epsf, P, show_dots);
          draw_sidewalk_cuts(epsf, P, show_dots);
          draw_building(epsf, P, show_dots);
          draw_perimeter_wall(epsf, P, show_dots);

          adrw_start_page(epsf, xmin, xmax, ymin, ymax, ox, nx, oy, ny, "Tubulacao e caixas de drenagem");
          draw_trimmed_sidewalk(epsf, P, show_dots);
          draw_drain_pipes(epsf, P, show_dots);
          draw_drain_boxes(epsf, P, show_dots);
          draw_building(epsf, P, show_dots);
          draw_perimeter_wall(epsf, P, show_dots);

          adrw_start_page(epsf, xmin, xmax, ymin, ymax, ox, nx, oy, ny, "Assentamento de ladrilhos");
          draw_trimmed_sidewalk(epsf, P, show_dots);
          draw_drain_boxes(epsf, P, show_dots);
          draw_sidewalk_additions(epsf, P, show_dots);
          draw_building(epsf, P, show_dots);
          draw_perimeter_wall(epsf, P, show_dots);

          epswr_end_figure(epsf);
          free(fname);
        }
  }

void draw_perimeter_wall(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    PSHAPE
      ( "pwall", "Perimeter wall", 1.000,0.750,0.250, show_dots,
         80,  81,  82,  83,  87,
         86,  85,  84,  80
      );
  }

void draw_old_sidewalk(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    PSHAPE
      ( "oswalk", "Old sidewalk", 0.800,0.800,0.950, show_dots,
         90,   3,   6,   7,   8,  
          9,  10,  11,  12,  13, 
         14,  16,  17,  19,  21, 
         24,  25,  27,  28,  28, 
         29,  30,  35,  36,  31,
         90
      );
  }

void draw_sidewalk_cuts(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    PSHAPE("cut_01", "Sidewalk cut-out", 1.000,0.250,0.000, show_dots,   3,   5,   4,   2,   3);
    PSHAPE("cut_02", "Sidewalk cut-out", 1.000,0.250,0.000, show_dots,  15,  14,  16,  17,  19,  21,  24,  23,  22,  20,  18,  15);
    PSHAPE("cut_03", "Sidewalk cut-out", 1.000,0.250,0.000, show_dots,  25,  27,  26,  25);
    PSHAPE("cut_04", "Sidewalk cut-out", 1.000,0.250,0.000, show_dots,  75,  76, 121, 123,  75);
    PSHAPE("cut_05", "Sidewalk cut-out", 1.000,0.250,0.000, show_dots,  77,  78,  91,  79,  77);
    PSHAPE("cut_06", "Sidewalk cut-out", 1.000,0.250,0.000, show_dots,  92,  93,  94,  95,  92);
    PSHAPE("cut_07", "Sidewalk cut-out", 1.000,0.250,0.000, show_dots, 100, 122,  99,  98,  48,  97,  96, 100);
    PSHAPE("cut_08", "Sidewalk cut-out", 1.000,0.250,0.000, show_dots, 101, 156, 142, 103, 104, 141, 155, 102, 101);
    PSHAPE("cut_09", "Sidewalk cut-out", 1.000,0.250,0.000, show_dots, 105, 107, 108, 144, 159, 106, 105);
    PSHAPE("cut_10", "Sidewalk cut-out", 1.000,0.250,0.000, show_dots, 109, 111, 112, 110, 109);
    PSHAPE("cut_11", "Sidewalk cut-out", 1.000,0.250,0.000, show_dots, 113, 115, 116, 148, 215, 114, 113);
    PSHAPE("cut_12", "Sidewalk cut-out", 1.000,0.250,0.000, show_dots, 117, 118, 120, 119, 117);
  }

void draw_trimmed_sidewalk(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    PSHAPE
      ( "tswalk", "Trimmed sidewalk", 0.800,0.800,0.950, show_dots,
         90,   2,   4,   5,   6,
          7,  75, 123, 121,  76,
          8,  77,  79,  91,  78,
          9,  92,  95,  94,  93,
         10,  11,  12,  96,
         97,  98,  14,  15, 103, 142, 156, 
        101, 102, 155, 141, 104,  18,  20,
        107, 105, 106, 159, 144, 108,  22,
         23, 111, 109, 110, 112, 
         26,  27, 115, 113, 114, 215, 148,
        116,  28, 120, 118, 117,
        119,  29,  30,  35,  36, 31,
         90
      );
    PSHAPE
      ( "tswalk", "Trimmed sidewalk", 0.800,0.800,0.950, show_dots,
        122, 100,  13,  99, 122
      );
  }
  
void draw_sidewalk_additions(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    PSHAPE("pav_00", "Sidewalk addition", 0.000,0.250,1.000, show_dots,  11,  37,  13, 100, 122, 139,  96,  12,  11);
    PSHAPE("pav_01", "Sidewalk addition", 0.000,0.250,1.000, show_dots,  32,   7,   6,   5,   4, 133, 132, 131,  32);
    PSHAPE("pav_02", "Sidewalk addition", 0.000,0.250,1.000, show_dots,  20,  18, 161,  33,  27,  26,  23,  22,  20);
    PSHAPE("pav_03", "Sidewalk addition", 0.000,0.250,1.000, show_dots,  34, 162, 160,  28,  34);
    PSHAPE("pav_04", "Sidewalk addition", 0.000,0.250,1.000, show_dots,  75,  76, 153, 152,  75);
    PSHAPE("pav_05", "Sidewalk addition", 0.000,0.250,1.000, show_dots,  77,  78, 136, 135,  77);
    PSHAPE("pav_06", "Sidewalk addition", 0.000,0.250,1.000, show_dots,  92,  93, 137, 138,  92);
    PSHAPE("pav_07", "Sidewalk addition", 0.000,0.250,1.000, show_dots, 122,  99,  98,  48, 122);
    PSHAPE("pav_08", "Sidewalk addition", 0.000,0.250,1.000, show_dots, 142, 103, 104, 141, 142);
    PSHAPE("pav_09", "Sidewalk addition", 0.000,0.250,1.000, show_dots, 143, 107, 108, 144, 143);
    PSHAPE("pav_10", "Sidewalk addition", 0.000,0.250,1.000, show_dots, 145, 163, 109, 111, 112, 110, 164, 146, 145);
    PSHAPE("pav_11", "Sidewalk addition", 0.000,0.250,1.000, show_dots, 147, 115, 116, 148, 147);
    PSHAPE("pav_12", "Sidewalk addition", 0.000,0.250,1.000, show_dots, 149, 150, 212, 118, 120, 119, 117, 211, 149);
  }
  
void draw_building(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    PSHAPE
      ( "house", "House", 1.000,0.750,0.500, show_dots,
         60,  63, 125, 126, 127,
        128,  38,  43,  44,  45,
         46,  47,  48,  49,  50,
         51,  52,  53,  54,  61,
         55,  56,  57,  58,  70,
         72,  73,  71,  65,  66, 
         67,  59,  60
      );
    PSHAPE
      ( "pil01", "Pillar", 1.000,0.750,0.500, show_dots,
         41, 89, 88, 39, 41
      );
  }
  
void draw_drain_boxes(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    DRBOX("box_01", "Drain box", 0.000,0.500,0.000, show_dots, 171);
    DRBOX("box_02", "Drain box", 0.000,0.500,0.000, show_dots, 172);
    DRBOX("box_03", "Drain box", 0.000,0.500,0.000, show_dots, 185);
    DRBOX("box_04", "Drain box", 0.000,0.500,0.000, show_dots, 173);
    DRBOX("box_05", "Drain box", 0.000,0.500,0.000, show_dots, 174);
    DRBOX("box_06", "Drain box", 0.000,0.500,0.000, show_dots, 175);
    DRBOX("box_07", "Drain box", 0.000,0.500,0.000, show_dots, 178);
    DRBOX("box_08", "Drain box", 0.000,0.500,0.000, show_dots, 176);
    DRBOX("box_09", "Drain box", 0.000,0.500,0.000, show_dots, 177);
    DRBOX("box_10", "Drain box", 0.000,0.500,0.000, show_dots, 179);
    DRBOX("box_11", "Drain box", 0.000,0.500,0.000, show_dots, 180);
    DRBOX("box_12", "Drain box", 0.000,0.500,0.000, show_dots, 181);
    DRBOX("box_13", "Drain box", 0.000,0.500,0.000, show_dots, 182);
    DRBOX("box_14", "Drain box", 0.000,0.500,0.000, show_dots, 183);
    DRBOX("box_15", "Drain box", 0.000,0.500,0.000, show_dots, 184);
  }
  
void draw_drain_pipes(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    DRPIPE("pipe_00", "Drain pipe", 0.000,0.500,0.750, show_dots, 185, 203);
    DRPIPE("pipe_01", "Drain pipe", 0.000,0.500,0.750, show_dots, 171, 192);
    DRPIPE("pipe_02", "Drain pipe", 0.000,0.500,0.750, show_dots, 192, 202);
    DRPIPE("pipe_03", "Drain pipe", 0.000,0.500,0.750, show_dots, 192, 172);
    DRPIPE("pipe_04", "Drain pipe", 0.000,0.500,0.750, show_dots, 172, 203);
    DRPIPE("pipe_05", "Drain pipe", 0.000,0.500,0.750, show_dots, 203, 193);
    DRPIPE("pipe_06", "Drain pipe", 0.000,0.500,0.750, show_dots, 193, 194);
    DRPIPE("pipe_07", "Drain pipe", 0.000,0.500,0.750, show_dots, 194, 173);
    DRPIPE("pipe_08", "Drain pipe", 0.000,0.500,0.750, show_dots, 194, 195);
    DRPIPE("pipe_09", "Drain pipe", 0.000,0.500,0.750, show_dots, 195, 174);
    DRPIPE("pipe_10", "Drain pipe", 0.000,0.500,0.750, show_dots, 175, 178);
    DRPIPE("pipe_11", "Drain pipe", 0.000,0.500,0.750, show_dots, 178, 196);
    DRPIPE("pipe_12", "Drain pipe", 0.000,0.500,0.750, show_dots, 196, 197);
    DRPIPE("pipe_13", "Drain pipe", 0.000,0.500,0.750, show_dots, 197, 177);
    DRPIPE("pipe_14", "Drain pipe", 0.000,0.500,0.750, show_dots, 177, 204);
    DRPIPE("pipe_15", "Drain pipe", 0.000,0.500,0.750, show_dots, 204, 176);
    DRPIPE("pipe_16", "Drain pipe", 0.000,0.500,0.750, show_dots, 204, 198);
    DRPIPE("pipe_17", "Drain pipe", 0.000,0.500,0.750, show_dots, 198, 179);
    DRPIPE("pipe_18", "Drain pipe", 0.000,0.500,0.750, show_dots, 198, 181);
    DRPIPE("pipe_19", "Drain pipe", 0.000,0.500,0.750, show_dots, 181, 180);
    DRPIPE("pipe_20", "Drain pipe", 0.000,0.500,0.750, show_dots, 181, 199);
    DRPIPE("pipe_21", "Drain pipe", 0.000,0.500,0.750, show_dots, 199, 182);
    DRPIPE("pipe_22", "Drain pipe", 0.000,0.500,0.750, show_dots, 199, 200);
    DRPIPE("pipe_23", "Drain pipe", 0.000,0.500,0.750, show_dots, 200, 201);
    DRPIPE("pipe_24", "Drain pipe", 0.000,0.500,0.750, show_dots, 201, 183);
    DRPIPE("pipe_25", "Drain pipe", 0.000,0.500,0.750, show_dots, 201, 190);
    DRPIPE("pipe_26", "Drain pipe", 0.000,0.500,0.750, show_dots, 190, 184);
  }
  
adrw_point_vec_t define_points(void)
  {
    adrw_point_vec_t P = adrw_point_vec_new(40);
    int np = 0;
    
    auto void s(int j, double dX, double dY, int i);
      /* Defines point {P[j]+(dX,dY)} as {P[i]}.  If {P[i]} was defined
        previously, checks whether the definitions agree. */
       
    void s(int j, double dX, double dY, int i)
      { char *lab = NULL;
        asprintf(&lab, "P%03d", i);
        adrw_append_point(lab, i, j, j, j, dX, dY, 0.0, &P, &np); 
      } 
    
    double Ybd = +500; /* Y of front of building. */

    /* Building points. */
    s( -1,    00,   Ybd,  60); /* dY = ??? */
    s( 60,    00,  +932,  63);
    s( 63,   +89,    00, 125);
    s(125,    00,   +10, 126);
    s(126,   +43,    00, 127);
    s(127,    00,   +19, 128);
    s(128,  +364,    00,  38);
    s( 38,    00,  +803,  43);
    s( 43,   +72,  +126,  44);
    s( 44,   +59,    00,  45);
    s( 45,    00,  +377,  46); 
    s( 46,  +303,    00,  47);
    s( 47,    00,  +175,  48);
    s( 48,  +433,    00,  49);
    s( 49,    00,  -467,  50);
    s( 50,   -23,    00,  51);
    s( 51,    00,  -557,  52);
    s( 52,   -50,   -85,  53);
    s( 53,    00,  -625,  54);
    s( 54,   +35,    00,  61);
    s( 61,    00,  -276,  55);
    s( 55,   +47,    00,  56);
    s( 56,    00,  -109,  57);
    s( 57,   -47,    00,  58);
    s( 61,    00,  -709,  59);
    s( 59,    00,   +26,  65);
    s( 65,    +5,    00,  66);
    s( 59,    +5,    00,  67);
    s( 59,  -445,    00,  62);
    s( 62,  -300,    00,  64); /* dX = ??? */
    
    /* Pillar in rear porch. */
    s( 63,  +163,  +332,  41);
    s( 41,    00,   -39,  39); /* dY = ??? */
    s( 41,   +39,    00,  89); /* dY = ??? */
    s( 89,    00,   -39,  88);
    
    /* Points of garden wall: */
    s( 55,    00,  -292,  70);
    s( 70,    00,   -22,  71);
    s( 70,  +176,    00,  72);
    s( 72,    00,   -22,  73);
    
    /* Points of perimeter wall: */
    s( 60,    00,  -Ybd,  80);
    s( 60,    00, +2646,  81); 
    s( 81, +1490,    00,  82);
    s( 82,    00, -2646,  83);
    s( 80,   -15,    00,  84);
    s( 81,   -15,   +15,  85);
    s( 82,   +15,   +15,  86);
    s( 83,   +15,    00,  87);
    
    /* Sidewalk points. */
    s( 60,    00,    00,  90);
    s( 64,    00,    00,  31);
    s( 63,    00,  +119,   3);
    s(  3,    00,  +267,  32);
    s(  3,  +140,    00,   6);
    s(  6,    00,  +267,   7);
    s(  7,  +267,    00,   8);
    s(  8,    00,  +475,   9);
    s(  9,  +131,  +238,  10);
    s( 10,    00,  +325,  11);
    s( 11,  +303,    00,  12);
    s( 12,    00,  +175,  13);
    s( 49,   +89,   +89, 158); /* Corner that was cut off. */
    s(158,   -38,    00,  14); 
    s(158,    00,   -28,  16);
    s( 50,   +89,   -80,  17);
    s( 51,   +89,  -109,  19);
    s( 52,   +89,   -26,  21);
    s( 21,   -44,   -89,  24);
    s( 54,   +95,  +120,  25);
    s( 25,   +29,   -29,  27);
    s( 59,   +89,   -89,  28);
    s( 62,    00,   -89,  29);
    s( 29,    00,  -280,  30);
    s( 29,  -300,    00,  36);
    s( 36,    00,  -280,  35);
    s( 43,   -89,   +29,   9); /* Redundant: */
    s( 46,   -89,   +89,  11); /* Redundant: */
    s(158,    00,  -636,  17); /* Redundant: */
    s( 19,    00,  -474,  21); /* Redundant: */
    s( 53,   +95,   -30,  24); /* Redundant: */
    s( 24,    00,  -475,  25); /* Redundant: */
    
    /* Points where garden wall crosses sidewalk: */
    s( 70,   +89,    00,  33); 
    s( 71,   +89,    00,  34);
    
    /* Points of cuts and repavs in sidewalk. */
    /* cut_01 = 3, 5, 4, 2, 3, */ 
    s(  3,    00,   -89,   2);
    s(  3,   +59,    00,   5);
    s(  2,   +59,    00,   4);
    
    /* cut_02 = 15, 14, 16, 17, 19, 21, 24, 23, 22, 20, 18, 15, */
    s( 14,   +29,   -29,  15);
    s( 17,    -9,    00,  18);
    s( 18,   -29,    00,  20);
    s( 21,   -15,   -26,  22);
    s( 22,   -29,    00,  23);
    /* Consistency: */
    /* s( 23,    00,   -29,  24); */

    /* cut_03 = 25, 27, 26, 25, */
    s( 25,    00,   -29,  26);
    
    /* cut_04 = 75, 76, 121, 123, 75, */
    s(  7,   +29,    00,  75);
    s( 75,   +29,    00,  76);
    s( 75,    00,   -54, 123);
    s( 76,    00,   -54, 121);
    
    /* cut_05 = 77, 78, 91, 79, 77, */
    s(  9,    00,   -30,  78);
    s( 78,    00,   -29,  77);
    s( 77,   +89,    00,  79);
    s( 78,   +89,    00,  91);
    
    /* cut_06 = 92, 93, 94, 95, 92, */
    s( 10,   +89,   -89,  94);
    s( 94,    00,   -29,  95);
    s( 94,  -138,    00,  93);
    s( 95,  -154,    00,  92);
    
    /* cut_07 = 100, 122, 99, 98, 48, 97, 96, 100, */
    s( 48,    00,   -29,  97);
    s( 97,   -89,    00,  96);
    s( 97,   -89,   +29, 100);
    s( 97,   -29,   +29, 122);
    s( 48,    00,   +89,  98);
    s( 48,   -29,   +89,  99);
    
    /* cut_08 = 101, 103, 104, 102, 101, */
    s( 49,    00,    +7, 101); /* dY = ??? */
    s(101,    00,   -29, 102);
    s( 49,   +29,    00, 142);
    s( 49,   +29,   -13, 141);
    s(101,   +29,    00, 156);
    s(102,   +29,    00, 155);
    s(142,   +51,    00, 103);
    s(141,   +51,    00, 104);
    
    /* cut_09 = 105, 107, 108, 144, 159, 106, 105, */
    s( 52,    00,    -6, 106);
    s( 52,   +29,    00, 144);
    s( 52,   +74,    00, 108);
    s(106,    00,   +29, 105);
    s(106,   +29,    00, 159);
    s(105,   +29,    00, 143);
    s(105,   +74,    00, 107);
    
    /* cut_10 = 109, 111, 112, 110, 109, */
    s( 24,    00,  -237, 111);
    s(111,    00,   -59, 112);
    s(111,   -95,    00, 109);
    s(112,   -95,    00, 110);
    s(110,    00,   +15, 164);
    s(164,    00,   +28, 163);
    s(164,   +29,    00, 146);
    s(163,   +29,    00, 145);
    
    /* cut_11 = 113, 115, 116, 148, 215, 114, 113, */
    s( 27,    00,   -59, 115);
    s(115,  -124,    00, 113);
    s(113,   +29,    00, 147);
    s(115,    00,   -29, 116);
    s( 54,    00,    00, 114);
    s(114,   +29,    00, 215);
    s(116,   -95,    00, 148);
    s(215,    00,    +3, 148); /* Redundant. */
    
    /* cut_12 = 117, 118, 120, 119, 117, */
    s( 28,   -59,    00, 120);
    s( 28,   -59,   +89, 118);
    s(118,   -58,    00, 117);
    s(120,   -58,    00, 119);
    
    /* pav_00 = 11, 37, 13, 100, 122, 139, 96, 12, 11 */
    s( 97,   -29,    00, 139);
    s( 11,    00,  +175,  37);

    /* pav_01 = 32, 7, 6, 5, 4, 133, 132, 131, 32 */
    s(  2,   +29,    00, 133);
    s(  2,   +29,   +29, 132);
    s(  2,    00,   +29, 131);

    /* pav_02 = 20, 18, 161, 33, 27, 26, 23, 22, 20 */
    s( 33,   +29,    00, 161);

    /* pav_03 = 34, 162, 160, 28, 34 */
    s( 34,   +29,    00, 162);
    s( 28,   +29,    00, 160);

    /* pav_04 = 75, 76, 153, 152, 75 */
    s(123,    00,   +29, 152);
    s(121,    00,   +29, 153);
    
    /* pav_05 = 77, 78, 136, 135, 77 */
    s( 79,   -29,    00, 135);
    s( 91,   -29,    00, 136);

    /* pav_06 = 92, 93, 137, 138, 92 */
    s( 94,   -29,    00, 137);
    s( 95,   -29,    00, 138);
    
    /* pav_07 = 122, 99, 98, 48, 122 */
    
    /* pav_08 = 142, 103, 104, 141, 142 */
    
    /* pav_09 = 143, 107, 108, 144, 143 */
    
    /* pav 10 = 145, 163, 109, 111, 112, 110, 164, 146, 145 */
    
    /* pav_11 = 147, 115, 116, 148, 147 */
    
    /* pav_12 = 149, 150, 120, 119, 149 */
    s( 59,   -16,    00, 211);
    s(211,   +28,    00, 212);
    s(211,    00,   -29, 149);
    s(212,    00,   -29, 150);

    /* Drain box centers: */
    s(  2,   +15,   +15, 171);  /* box_01 */
    s( 32,   +57,   +20, 172);  /* box_02 */
    s( 41,   +20,   +15, 185);  /* box_03 */
    s( 79,   -15,   +15, 173);  /* box_04 */ 
    s( 95,   -15,   +14, 174);  /* box_05 */
    s( 48,   -14,   -11, 175);  /* box_06 */
    s(175,    00,  +120, 178);  /* box_07 */
    s( 49,   +15,    -6, 176);  /* box_08 */
    s(158,   +11,   -30, 177);  /* box_09 */
    s( 52,   +15,    +8, 179);  /* box_10 */
    s( 53,   +15,  -294, 180);  /* box_11 */
    s(180,  +159,    00, 181);  /* box_12 */
    s(114,   +15,   +14, 182);  /* box_13 */
    s( 59,    -2,   -15, 183);  /* box_14 */
    s( 29,   +44,  -247, 184);  /* box_15 */ 
    
    /* Pipe ends and joints: */
    s(  3,   +57,    00, 191);
    s(  2,   +57,   +15, 192);
    s(  8,   -20,   +20, 193);
    s(173,   -94,    00, 194);
    s(174,  -225,    00, 195);
    s(158,   -40,   +20, 196);
    s(158,   +11,   -21, 197);
    s(179,  +106,    00, 198);
    s(182,  +156,    00, 199);
    s(160,   +20,   -20, 200);
    s(183,    00,   -94, 201);
    s(192,    00,   -30, 202);
    s( 41,   +20,   +74, 203);
    s(176,   +85,    00, 204);
    s(184,    00,  +228, 190);
    
    adrw_point_vec_trim(&P, np);
    show_holes(&P, 0, P.ne-1);
    
    return P;
  }
