/* See {boas_draw.h} */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <bool.h>
#include <vec.h>
#include <epswr.h>
#include <frgb.h>
#include <archdraw.h>

#include <boas_draw.h>

#define RGB_GRASS_IN   (frgb_t){{ 0.650f, 1.000f, 0.650f }}
#define RGB_GRASS_OUT  (frgb_t){{ 0.550f, 0.900f, 0.550f }}
#define RGB_CURB       (frgb_t){{ 0.650f, 0.450f, 0.650f }}

#define RGB_WALL       (frgb_t){{ 1.000f, 0.850f, 0.650f }}
#define RGB_WALL_CUT   (frgb_t){{ 1.000f, 0.250f, 0.150f }}
#define RGB_WALL_NEW   (frgb_t){{ 0.300f, 0.200f, 1.000f }}

#define RGB_FLOOR_WALK (frgb_t){{ 0.800f, 0.800f, 0.950f }} /* Walkway of blue-gray tiles. */
#define RGB_FLOOR_GAR  (frgb_t){{ 0.950f, 0.940f, 0.800f }} /* Garage floor,bege-gray tiles. */
#define RGB_FLOOR_RAW  (frgb_t){{ 0.700f, 0.700f, 0.750f }} /* Rough granite tiles on sidewalk and driveway. */
#define RGB_FLOOR_CUT  (frgb_t){{ 1.000f, 0.100f, 0.250f }} 
#define RGB_FLOOR_NEW  (frgb_t){{ 0.200f, 0.300f, 1.000f }}

#define RGB_DRAIN_BOX  (frgb_t){{ 0.400f, 0.600f, 0.400f }}
#define RGB_DRAIN_PIPE (frgb_t){{ 0.400f, 0.400f, 0.800f }}
#define RGB_BUILDING   (frgb_t){{ 1.000f, 0.850f, 0.650f }}

void plot_shape
  ( epswr_figure_t *epsf, 
    char *label, 
    char *desc,
    frgb_t color, 
    bool_t show_dots,
    adrw_point_vec_t *P,
    int v[]
  );

void plot_box
  ( epswr_figure_t *epsf, 
    char *label,
    char *desc, 
    frgb_t color,
    bool_t show_dots,
    adrw_point_vec_t *P,
    int ipt
  );

void plot_pipe
  ( epswr_figure_t *epsf, 
    char *label, 
    char *desc,
    frgb_t color, 
    bool_t show_dots,
    adrw_point_vec_t *P,
    int v[]
  );

#define PSHAPE(LABEL,DESC,COLOR,SHOW_DOTS,...)  \
  do \
    { int v[] = { __VA_ARGS__ , -1 }; \
      plot_shape(epsf,(LABEL),(DESC),(COLOR),(SHOW_DOTS),P,v); \
    } \
  while(0)

#define DRBOX(LABEL,DESC,COLOR,SHOW_DOTS,PT) \
  do \
    { \
      plot_box(epsf,(LABEL),(DESC),(COLOR),(SHOW_DOTS),P,(PT)); \
    } \
  while(0)

#define DRPIPE(LABEL,DESC,COLOR,SHOW_DOTS,PTA,PTB) \
  do \
    { int v[] = { (PTA), (PTB), -1 }; \
      plot_pipe(epsf,(LABEL),(DESC),(COLOR),(SHOW_DOTS),P,v); \
    } \
  while(0)
  
void boas_draw_land_plot(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    PSHAPE
      ( "landp", "Land plot", RGB_GRASS_IN, show_dots,
          80,  63,  81,  82,  83,  80
      );

    PSHAPE
      ( "swalk", "Owned sidewalk", RGB_GRASS_OUT, show_dots,
          4,   5,  83,  80,   4
      );

    PSHAPE
      ( "swalk", "Owned curb wall", RGB_CURB, show_dots,
           1,   2,  5,   4,   1
      );
  }

void boas_draw_cuts_walls(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    PSHAPE
      ( "swall", "Remove side door wall", RGB_WALL_CUT, show_dots,
         72,  70,  71,  73,  72
      );

    PSHAPE
      ( "fdoor1", "Remove front wall door jamb 1", RGB_WALL_CUT, show_dots,
        304, 351, 352, 362, 361, 344, 334, 304
      );
    PSHAPE
      ( "fdoor2", "Remove front wall door jamb 2", RGB_WALL_CUT, show_dots,
        305, 335, 345, 364, 363, 353, 354, 305
      );

  }

void boas_draw_trimmed_walls(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    PSHAPE
      ( "pwall", "Trimmed perimeter wall", RGB_WALL, show_dots,
         80,  63,  81,  82,  
         83,  87,  86,  85,  64,  84,  80
      );

    PSHAPE
      ( "fwall1", "Trimmed front wall section 1", RGB_WALL, show_dots,
         80, 301, 323, 322,  80
      );
    PSHAPE
      ( "fwall1", "Trimmed front wall section 2", RGB_WALL, show_dots,
         80, 301, 323, 322,  80
      );
    PSHAPE
      ( "fwall2", "Front wall section 3", RGB_WALL, show_dots,
        308, 309, 325, 324, 308
      );
    PSHAPE
      ( "fwall3", "Trimmed front wall section 3", RGB_WALL, show_dots,
        310, 311, 327, 326, 310
      );
    PSHAPE
      ( "fwall4", "Trimmed front wall section 4", RGB_WALL, show_dots,
        312,  83, 329, 328, 312
      );

    PSHAPE
      ( "fpill1", "Trimmed front wall pillar 1", RGB_WALL, show_dots,
        301, 302, 332, 331, 301
      );
    PSHAPE
      ( "fpill2", "Trimmed front wall pillar 2", RGB_WALL, show_dots,
        303, 304, 334, 333, 303
      );
    PSHAPE
      ( "fpill3", "Trimmed front wall pillar 3", RGB_WALL, show_dots,
        305, 306, 336, 335, 305
      );
    PSHAPE
      ( "fpill4", "Trimmed front wall pillar 4", RGB_WALL, show_dots,
        307, 308, 338, 337, 307
      );
    PSHAPE
      ( "fpill5", "Front wall pillar 5", RGB_WALL, show_dots,
        309, 310, 340, 339, 309
      );
  }

void boas_draw_new_walls(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    PSHAPE
      ( "swall", "New side door wall", RGB_WALL_NEW, show_dots,
        370, 371, 373, 372, 370
      );
    PSHAPE
      ( "rwall", "New side roof wall", RGB_WALL_NEW, show_dots,
        383, 382, 380, 381, 383
      );
  }

void boas_draw_cuts_pavement(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    // PSHAPE
    //   ( "cwalk", "Cut on tiled walkway", RGB_FLOOR_CUT, show_dots,
    //     370, 371, 375, 374, 370
    //   );

    PSHAPE
      ( "cyard1", "Cut on rough yard pavement 1", RGB_FLOOR_CUT, show_dots,
        481, 331, 332, 302, 512, 532, 481
      );
    PSHAPE
      ( "cyard2", "Cut on rough yard pavement 2", RGB_FLOOR_CUT, show_dots,
        533, 513, 303, 486, 533
      );
    PSHAPE
      ( "cyard3", "Cut on rough yard pavement 3", RGB_FLOOR_CUT, show_dots,
        408, 410, 460, 457, 408
      );
    PSHAPE
      ( "cyard5", "Cut on rough yard pavement 5", RGB_FLOOR_CUT, show_dots,
        491, 306, 516, 536, 491
      );
    PSHAPE
      ( "cyard6", "Cut on rough yard pavement 6", RGB_FLOOR_CUT, show_dots,
        537, 517, 307, 492, 537
      );

    PSHAPE
      ( "cside1", "Cut on rough sdewalk pavement 1", RGB_FLOOR_CUT, show_dots,
        302, 403, 522, 512, 302
      );
    PSHAPE
      ( "cside2", "Cut on rough sdewalk pavement 2", RGB_FLOOR_CUT, show_dots,
        513, 523, 406, 303, 513
      );
    PSHAPE
      ( "cside3", "Cut on rough sdewalk pavement 3", RGB_FLOOR_CUT, show_dots,
        407, 410, 409, 408, 407
      );
    PSHAPE
      ( "cside5", "Cut on rough sdewalk pavement 5", RGB_FLOOR_CUT, show_dots,
        306,  411, 526, 516, 306
      );
    PSHAPE
      ( "cside6", "Cut on rough sdewalk pavement 6", RGB_FLOOR_CUT, show_dots,
        517, 527, 414, 307, 517
      );

  }
  
void boas_draw_trimmed_pavement(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    PSHAPE
      ( "twalk1", "Trimmed walkway north", RGB_FLOOR_WALK, show_dots,
         63,   3,   6,   7,   8,  
          9,  10,  11,  12,  13, 
         14,  16,  74,  70,  
         58,  57,  56,  55,  54,  53,  52,  51,  50,  49,
         48,  47,  46,  45,  44,  43,  38,
        128, 127,  
         63
      );

    PSHAPE
      ( "twalk1", "Trimmed walkway middle", RGB_FLOOR_WALK, show_dots,
         75,  71, 370, 374,  
         75
      );

    PSHAPE
      ( "twalk2", "Trimmed walkway south", RGB_FLOOR_WALK, show_dots,
        375,  28,  29,  90,  67,  66,  65, 371,
        375
      );

    PSHAPE
      ( "tgaraf", "Trimmed garage floor", RGB_FLOOR_GAR, show_dots,
         62,  68,  69, 100,  98,  96,  91,  
         92,  93,  99, 103, 104,
         62
      );

    PSHAPE
      ( "tyardp", "Trimmed rough tiles in front yard", RGB_FLOOR_RAW, show_dots,
        450, 451, 481, 532, 512, 513, 533, 486, 456, 
        461, 491, 536, 516, 
        517, 537, 492,  29,  90,  60, 
        450
      );

    PSHAPE
      ( "tsidew", "Trimmed rough tiles on sidewalk", RGB_FLOOR_RAW, show_dots,
         401, 402, 522, 512, 513, 523, 526, 516, 
         517, 527, 415, 416, 
         401
      );
  }

void boas_draw_new_pavement(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    PSHAPE
      ( "swall", "Fill cut on walkway", RGB_FLOOR_NEW, show_dots,
         70,  71,  75,  74,  70
      );

    PSHAPE
      ( "swall", "Pave driveway 2", RGB_FLOOR_NEW, show_dots,
         524, 525, 535, 490, 491, 461, 
         456, 486, 534, 524
      );
  }
  
void boas_draw_drain_boxes(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    DRBOX("box_00", "Drain box", RGB_DRAIN_BOX, show_dots, 170);
    DRBOX("box_01", "Drain box", RGB_DRAIN_BOX, show_dots, 171);
    DRBOX("box_02", "Drain box", RGB_DRAIN_BOX, show_dots, 172);
    DRBOX("box_03", "Drain box", RGB_DRAIN_BOX, show_dots, 185);
    DRBOX("box_04", "Drain box", RGB_DRAIN_BOX, show_dots, 173);
    DRBOX("box_05", "Drain box", RGB_DRAIN_BOX, show_dots, 174);
    DRBOX("box_06", "Drain box", RGB_DRAIN_BOX, show_dots, 175);
    DRBOX("box_07", "Drain box", RGB_DRAIN_BOX, show_dots, 178);
    DRBOX("box_08", "Drain box", RGB_DRAIN_BOX, show_dots, 176);
    DRBOX("box_09", "Drain box", RGB_DRAIN_BOX, show_dots, 177);
    DRBOX("box_10", "Drain box", RGB_DRAIN_BOX, show_dots, 179);
    DRBOX("box_11", "Drain box", RGB_DRAIN_BOX, show_dots, 180);
    DRBOX("box_12", "Drain box", RGB_DRAIN_BOX, show_dots, 181);
    DRBOX("box_13", "Drain box", RGB_DRAIN_BOX, show_dots, 182);
    DRBOX("box_14", "Drain box", RGB_DRAIN_BOX, show_dots, 183);
    DRBOX("box_15", "Drain box", RGB_DRAIN_BOX, show_dots, 184);
    DRBOX("box_16", "Drain box", RGB_DRAIN_BOX, show_dots, 186);
  }
  
void boas_draw_drain_pipes(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    DRPIPE("pipe_00", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 185, 203);
    DRPIPE("pipe_01", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 171, 192);
    DRPIPE("pipe_02", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 192, 202);
    DRPIPE("pipe_03", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 192, 172);
    DRPIPE("pipe_04", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 172, 203);
    DRPIPE("pipe_05", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 203, 193);
    DRPIPE("pipe_06", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 193, 194);
    DRPIPE("pipe_07", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 194, 173);
    DRPIPE("pipe_08", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 194, 195);
    DRPIPE("pipe_09", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 195, 174);
    DRPIPE("pipe_10", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 175, 178);
    DRPIPE("pipe_11", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 178, 196);
    DRPIPE("pipe_12", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 196, 197);
    DRPIPE("pipe_13", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 197, 177);
    DRPIPE("pipe_14", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 177, 204);
    DRPIPE("pipe_15", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 204, 176);
    DRPIPE("pipe_16", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 204, 198);
    DRPIPE("pipe_17", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 198, 179);
    DRPIPE("pipe_18", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 198, 181);
    DRPIPE("pipe_19", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 181, 180);
    DRPIPE("pipe_20", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 181, 199);
    DRPIPE("pipe_21", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 199, 182);
    DRPIPE("pipe_22", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 199, 200);
    DRPIPE("pipe_23", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 200, 201);
    DRPIPE("pipe_24", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 201, 183);
    DRPIPE("pipe_25", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 201, 186);
    DRPIPE("pipe_26", "Drain pipe", RGB_DRAIN_PIPE, show_dots, 190, 184);
  }

void boas_draw_building(epswr_figure_t *epsf, adrw_point_vec_t *P, bool_t show_dots)
  {
    PSHAPE
      ( "house", "House", RGB_BUILDING, show_dots,
         63, 127, 128, 115, 116, 117, 118, 
         38,  43,  44,  45,
         46,  47,  48,  49,  50,
         51,  52,  53,  54,  61,
         55,  56,  57,  58,  65,  66,  67,
         62,  68,  69, 100,  98,  96,  91,  92,  93,
         63
      );

    PSHAPE
      ( "rpil0", "Rear porch pillar", RGB_BUILDING, show_dots,
         41, 89, 88, 39, 41
      );

    PSHAPE
      ( "fpil0", "Facade pillar 0", RGB_BUILDING, show_dots,
        114, 104, 103, 113, 114
      );

    PSHAPE
      ( "fpil1", "Facade pillar 1", RGB_BUILDING, show_dots,
        108, 109, 139, 138, 108
      );

    PSHAPE
      ( "fpil2", "Facade pillar 2", RGB_BUILDING, show_dots,
        110, 111, 141, 140, 110
      );
  }
  
void plot_shape
  ( epswr_figure_t *epsf, 
    char *label, 
    char *desc,
    frgb_t color, 
    bool_t show_dots,
    adrw_point_vec_t *P,
    int v[]
  )
  {
    frgb_t black = (frgb_t){{ 0.0f, 0.0f, 0.0f }}; 
    adrw_unit_style_t *st = adrw_make_unit_style(&color, &black, 0.10, NULL, 0.0); 
    int kfloor = 0, type = 0; 
    adrw_unit_t *rm = adrw_make_poly(label, desc, 0.0, P, v, 0.0, kfloor, type, st); 
    adrw_plot_unit(epsf, rm, show_dots); 
    free(st); free(rm->pt.e); free(rm); 
  }
      
void plot_box
  ( epswr_figure_t *epsf, 
    char *label,
    char *desc, 
    frgb_t color,
    bool_t show_dots,
    adrw_point_vec_t *P,
    int ipt
  )
  {
    frgb_t black = (frgb_t){{ 0.0f, 0.0f, 0.0f}}; 
    adrw_unit_style_t *st = adrw_make_unit_style(&color, &black, 0.10, NULL, 0.0); 
    int kfloor = 0, type = 0; 
    double bsize = 28; /* Box width and height. */
    double round = 0.0; /* Corner rounding radius. */
    adrw_unit_t *rm = adrw_make_box(label, desc, P, ipt, 0.0,0.0,0.0, bsize, bsize, round, kfloor, type, st); 
    adrw_plot_unit(epsf, rm, show_dots); 
    free(st); free(rm->pt.e); free(rm); 
  }

void plot_pipe
  ( epswr_figure_t *epsf, 
    char *label, 
    char *desc,
    frgb_t color, 
    bool_t show_dots,
    adrw_point_vec_t *P,
    int v[]
  )
  {
    adrw_unit_style_t *st = adrw_make_unit_style(NULL, &color, 0.50, NULL, 0.0); 
    int kfloor = 0, type = 0; 
    adrw_unit_t *rm = adrw_make_poly(label, desc, 0.0, P, v, 0.0, kfloor, type, st); 
    adrw_plot_unit(epsf, rm, show_dots); 
    free(st); free(rm->pt.e); free(rm); 
  }
