/* See {boas_define_points.h} */
/* Last edited on 2023-02-20 17:25:03 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <r2.h>
#include <bool.h>
#include <vec.h>
#include <epswr.h>
#include <affirm.h>
#include <jsfile.h>
#include <archdraw.h>

#include <boas_define_points.h>

adrw_point_vec_t boas_define_points(void)
  {
    adrw_point_vec_t P = adrw_point_vec_new(40);
    int np = 0;
    
    auto void s(int i, int jx, int jy, double dX, double dY);
      /* Defines point {P[i]} as {(P[jx].x,P[jy].y)+(dX,dY)}.  If {P[i]} was defined
        previously, checks whether the definitions agree. */
       
    void s(int i, int jx, int jy, double dX, double dY)
      { char *lab = NULL;
        char *lab = jsprintf("P%03d", i);
        adrw_append_point(lab, i, jx, jy, -1, dX, dY, 0.0, &P, &np); 
      } 
    
    double Xtf = +1509; /* X extent of land plot at street wall. */
    double Xtr = +1483; /* X extent of land plot at rear wall. */

    double Yfw =  +296; /* Y of front of street wall rel to curb edge. */
    /* double Ytt = +2990; */ /* Y extent of land plot (excl. sidewalk, curb, border walls). */
    
    double Ybd =  +398; /* Y of front of building rel to front of street wall. */
    double Ypo =  +961; /* Y of back of rear wall of living room rel front of building. */
    double Ysr =   +17; /* Y of back rear wall of living room rel back of rear wall of storage room. */
    double Ycs =   +15; /* Y width of curb wall. */

    double Ypq =  +119; /* Y of front of porchside flowerbed rel rear wall of living room. */
    double Ypp =  +358; /* Y of back edge of porch pavement rel rear wall of living room. */
    
    double Wot =    12; /* Thickness of outer walls. */
    double Wft =    22; /* Thickness of front walls. */
    double Wfp =    32; /* X and Y size of front wall pillars. */
    double Wgp =    35; /* X and Y size of garage and back porch pillars. */
    double Wcg =    24; /* W distance from front pillar to car way. */
    
    /* double Ygr =  +164; */ /* Y dist from curb edge to grassy areas on sidewalk */
   
    /* The "owned" land plot is the area delimited by the W, N, and E perimeter walls,
      excluding said walls, and the front wall, including it.  Plus the sidewalk
      and curb wall due south of that area. 
      
      Note that the rough pavement that we put on the sidewalk extends 
      a bit ouside our "owned" sidewalk on the W side, as if the W perimeter 
      wall belonged to us (which it doesn't). */
    
    /* Main reference points on left edge of claimed land plot: */
    s( 80,  -1,  -1,    00,    00); /* SW corn land plot = left front corn street wall. */
    s(  1,  80,  80,    00,  -Yfw); /* SW corn of owned sidewalk with curb wall. */
    s( 60,  80,  80,    00,  +Ybd); /* Left front corn building. */
    s( 40,  60,  60,    00,  +Ypo); /* Back left of rear wall of living room. */
    s( 63,  40,  40,    00,  -Ysr); /* Back left of rear wall of storage room. */
    s(  3,  40,  40,    00,  +Ypq); /* Front edge of porchside flowebed. */
    s( 42,  40,  40,    00,  +Ypp); /* Back edge of porch pavement. */

    /* Main corners of building: */
    s( 67,  60,  60, +1336,    00); /* SE corn, from SW corn building. */
    s( 49,  67,  67,   +35, +2436); /* NE corn, from SE corn building. */
    
    s( 38,  40,  40,  +500,    00); /* SE corn back porch walls. */

    /* Corners of land plot (excl. sidewalk and curb raised blocks): */
    s( 83,  80,  80,  +Xtf,    00); /* SE corn = right front corn street wall. */
    s( 82,  49,  49,  +129,  +156); /* NE corn. */
    s( 81,  82,  82,  -Xtr,    -4); /* NW corn. */

    /* Corners of owned sidewalk: */
    s(  2,   1,   1,  +Xtf,    00); /* SE corn of owned sidewalk with curb wall. */
    s(  4,  80,   1,    00,  +Ycs); /* SW corn of owned sidewalk without curb wall. */
    s(  5,   2,   4,    00,    00); /* SW corn of owned sidewalk without curb wall. */

    /* Outer points of perimeter wall: */
    s( 64,  63,  63,  -Wot,    00);
    s( 84,  80,  80,  -Wot,    00);
    s( 85,  81,  81,  -Wot,  +Wot);
    s( 86,  82,  82,  +Wot,  +Wot);
    s( 87,  83,  83,  +Wot,    00);
    
    /* Points of front wall: */
    s(301,  80,  80,   +22,    00);
    s(302, 301,  80,  +Wfp,    00); 
    s(303, 302,  80,  +252,    00); 
    s(304, 303,  80,  +Wfp,    00); 
    s(305, 304,  80,  +250,    00); 
    s(306, 305,  80,  +Wfp,    00); 
    s(307, 306,  80,  +256,    00); 
    s(308, 307,  80,  +Wfp,    00); 
    s(309, 308,  80,  +181,    00); 
    s(310, 309,  80,  +Wfp,    00); 
    s(311, 310,  80,  +181,    00); 
    s(312, 311,  80,    00,    00); /* Repeated point. */

    s(312,  83,  83,  -178,    00); /* !!! Check !!! Redundant. */

    s(321,  84,  84,    00,  +Wft); 
    s(322,  80, 321,    00,    00); 
    s(323, 301, 321,    00,    00); 
    s(324, 308, 321,    00,    00); 
    s(325, 309, 321,    00,    00); 
    s(326, 310, 321,    00,    00); 
    s(327, 311, 321,    00,    00); 
    s(328, 312, 321,    00,    00); 
    s(329,  83, 321,    00,    00); 
    s(330,  87, 321,    00,    00); 

    s(329, 322, 322,  +Xtf,    00); /* Redundant. */

    s(331, 301, 301,    00,  +Wfp);
    s(332, 302, 331,    00,    00); 
    s(333, 303, 331,    00,    00); 
    s(334, 304, 331,    00,    00); 
    s(335, 305, 331,    00,    00); 
    s(336, 306, 331,    00,    00); 
    s(337, 307, 331,    00,    00); 
    s(338, 308, 331,    00,    00); 
    s(339, 309, 331,    00,    00); 
    s(340, 310, 331,    00,    00); 

    s(340,  83, 331,   +50,    00); /* !!! Check !!! Redundant. */

    s(344, 333, 333,   +27,    00); 
    s(345, 336, 336,   -27,    00); 

    s(361, 304, 304,   +28,   +94);
    s(362, 361, 361,   +60,    00); 
    s(363, 362, 361,  +100,    00); 
    s(364, 363, 361,   +25,    00);

    s(364, 305, 361,   -50,    00); /* !!! Check !!! Redundant. */

    s(352, 362, 362,    00,  -Wft); 
    s(351, 352, 352,   -48,    00); 
    s(353, 363, 352,    00,    00); 
    s(354, 353, 352,   +14,    00); /* !!! Check !!! */

    /* Building points. */
    s(127,  63,  63,  +135,    00);

    s(128, 127, 127,    00,   +17);
    s(128,  40,  40,  +135,    00); /* Redundant. */
    
    s(115, 128, 128,  +102,    00);
    s(116, 115, 115,    00,   -30);
    s(117, 116, 116,  +180,    00);
    s(118, 117, 117,    00,   +30);
    
    s( 38, 118, 118,   +83,    00); /* Redundant, from living room back door */
    s( 38, 128, 128,  +365,    00); /* Redundant, from storage room recess. */

    s( 43,  38,  38,    00,  +802);
    s( 44,  43,  43,   +74,  +129);
    s( 45,  44,  44,   +58,    00); /* Corn outer walls of bathroom and office. */
    s( 46,  45,  45,    00,  +375); 
    s( 47,  46,  46,  +306,    00);
    s( 48,  47,  47,    00,  +176); /* NW outer corn main bedroom. */

    /* s( 49,  48,  48,  +433,    00); */ /* !!! NE corn building. Redundant, from back NW outer corn bedroom. */
    /* s( 49,  38,  38,  +868, +1482); */ /* !!! NE corn building. Redundant, from SE corn back porch. */
    
    s( 50,  49,  49,    00,  -466);
    s( 51,  50,  50,   -23,    00);
    s( 52,  51,  51,    00,  -560);
    s( 53,  52,  52,   -52,   -85);
    s( 54,  53,  53,    00,  -625);
    s( 61,  54,  54,   +35,    00);
    s( 55,  61,  61,    00,  -155);
    s( 56,  55,  55,   +47,    00);
    s( 57,  56,  56,    00,  -107);
    s( 58,  57,  57,   -47,    00);
    
    s( 59,  61,  61,    00,  -700);
    s( 65,  59,  59,    00,   +23);
    s( 66,  65,  65,    +5,    00);
    
    s( 67,  59,  59,    +5,    00); /* Redundant. */
    
    s( 62,  67,  67,  -456,    00);
    
    s( 68,  62,  62,    00,  +Wgp);
    s( 69,  68,  68,    +7,    00);
    s(100,  69,  69,    00,  +505);
    
    s( 98, 100, 100,  -283,    00);
    s( 96,  98,  98,    00,  -245);
    s( 91,  96,  96,  -296,    00);
    s( 92,  91,  91,    00,  +245);
    s( 93,  92,  92,  -308,    00);
    s( 93,  60,  60,    00,  +540);

    /* Facade pillars: */
    s( 99,  60,  68,    00,    00);
    s(103,  99,  99,    +7,    00);
    s(104, 103,  60,    00,    00);
    s(107, 104, 104,    00,    00); /* Replace by 104. */
    s(113,  99,  99,  -Wot,    00);
    s(114, 113, 104,    00,    00);

    s(108, 107,  60,  +291,    00);
    s(109, 108,  60,  +Wgp,    00);
    s(138, 108,  99,    00,    00);
    s(139, 109,  99,    00,    00);
    
    s(110, 109,  60,  +242,    00);
    s(111, 110,  60,  +Wgp,    00);
    s(140, 110,  99,    00,    00);
    s(141, 111,  99,    00,    00);

    /* Pillar in rear porch. */
    s( 39,  40,  40,  +164,  +268); 
    s( 41,  39,  39,    00,  +Wgp);
    s( 88,  39,  39,  +Wgp,    00);
    s( 89,  88,  41,    00,    00); 
 
    /* Points of smooth tiled walkway around the house. */
    s( 90, 307,  60,    00,    00);
    s( 29,  90,  90,    00,   -92);
    s( 28,  29,  29,  +567,    00);
    s( 28,  60,  60, +1449,   -92); /* Redundant. */
    
    s( 31,  64,  64,    00,    00);
    s( 32,   3,   3,    00,  +267);
    s(  6,   3,   3,   +30,    00);
    s(  7,   6,   6,    00,  +239);
    s(  7,  42,  42,   +30,    00); /* Redundant. */
    s(  8,   7,   7,  +267,    00);
    s(  9,   8,   8,    00,  +475);
    s( 10,   9,   9,  +131,  +238);
    s( 11,  10,  10,    00,  +325);
    s( 12,  11,  11,  +303,    00);
    s( 13,  12,  12,    00,  +175);

    s(158,  49,  49,   +89,   +89); /* Corner that was cut off. */
    s( 14, 158, 158,   -38,    00); 
    s( 16, 158, 158,    00,   -28);
    
    /* Points of east garden door wall: */
    s( 70,  58,  58,    00,  -182);
    s( 71,  70,  70,    00,  -Wft);
    s( 72,  83,  70,    00,    00);
    s( 73,  72,  71,    00,    00);
    
    s( 74,  28,  70,    00,    00);  /* On edge of walkway. */
    s( 75,  28,  71,    00,    00);  /* On edge of walkway. */
    
    s( 59,  71,  71,    00,  -233);  /* Redundant. */
    
    /* Points of displaced east garden door wall: */
    s(371,  65,  65,    00,   +20); /* !!! Check !!! */
    s(370, 371, 371,    00,  +Wft);
    s(372,  83, 370,    00,    00);
    s(373, 372, 371,    00,    00);
    
    s(374,  28, 370,    00,    00);  /* On edge of walkway. */
    s(375,  28, 371,    00,    00);  /* On edge of walkway. */
    
    s(383, 372, 372,    00,    00);
    s(382, 372, 372,    00,  +400); /* !!! Check !!! */
    s(380, 382, 382,  -Wot,    00);
    s(381, 380, 383,    00,    00);
  
    /* Points of rough pavement on sidewalk: */
    s(401,  84,   4,    00,    00);  
    s(402, 401,  84,    00,   -48);
    s(403, 302, 402,    00,    00);
    s(404, 403,  84,    00,    00);
    s(405, 303,  84,    00,    00);
    s(406, 405,  84,    00,  -130);
    s(407, 352, 406,    00,    00);
    s(408, 407,  84,    00,    00);
    s(409, 353,  84,    00,    00);
    s(410, 409, 406,    00,    00);
    s(411, 306, 406,    00,    00);
    s(412, 411,  84,    00,    00);
    s(413, 307,  84,    00,    00);
    s(414, 413, 406,    00,    00);
    s(415,  83, 406,    00,    00);
    s(416, 415, 401,    00,    00);

    /* Points of rough pavement on front yard: */
    s(450,  60,  60,    00,   -94); /* Westmost corner, S of flowerbed. */
    s(451, 450, 450,   +22,    00);
    s(456, 303, 450,    00,    00);
    s(457, 352, 450,    00,    00);
    s(458, 457,  84,    00,    00);
    s(459, 353,  84,    00,    00);
    s(460, 459, 450,    00,    00);
    s(461, 306, 450,    00,    00);
    
    /* Points on cuts/additions to rough pavement on front yard: */
    s(502, 302, 451,  +Wcg,    00);
    s(503, 303, 451,  -Wcg,    00);
    s(504, 304, 451,  +Wcg,    00);
    s(505, 305, 451,  -Wcg,    00);
    s(506, 306, 451,  +Wcg,    00);
    s(507, 307,  29,  -Wcg,    00);

    s(512, 502,  80,    00,    00);
    s(513, 503,  80,    00,    00);
    s(514, 504,  80,    00,    00);
    s(515, 505,  80,    00,    00);
    s(516, 506,  80,    00,    00);
    s(517, 507,  80,    00,    00);

    s(481, 451, 451,    00,   -59);
    s(486, 456, 481,    00,    00);
    s(487, 457, 481,    00,    00);
    s(490, 460, 481,    00,    00);
    s(491, 461, 481,    00,    00);
    s(492, 337, 481,    00,    00);

    s(532, 502, 481,    00,    00);
    s(533, 503, 481,    00,    00);
    s(534, 504, 481,    00,    00);
    s(535, 505, 481,    00,    00);
    s(536, 506, 481,    00,    00);
    s(537, 507, 481,    00,    00);

    /* Points on cuts/additions to rough pavement on sidewalk: */
    s(522, 502, 403,    00,    00);
    s(523, 503, 406,    00,    00);
    s(524, 504, 406,    00,    00);
    s(525, 505, 406,    00,    00);
    s(526, 506, 406,    00,    00);
    s(527, 507, 406,    00,    00);

    /* Drain box centers: */
    s(170,  40,  40,   +52,   +98);  /* box_00 */
    s(171,  40,  40,   +15,   +32);  /* box_01 */
    s(172,  41,  41,   -90,   +74);  /* box_02 */
    s(185,  41,  41,   +17,   +15);  /* box_03 */
    s(173,  43,  43,   -15,   -15);  /* box_04 */ 
    s(174,  45,  45,   -15,   +42);  /* box_05 */
    s(175,  47,  47,   -15,  +162);  /* box_06 */
    s(178,  48,  48,   -15,  +109);  /* box_07 */
    s(176,  49,  49,   +15,   -12);  /* box_08 */
    s(177,  82,  49,   -17,   +53);  /* box_09 */
    s(179,  52,  52,   +17,   +17);  /* box_10 */
    s(180,  54,  54,   +15,  +329);  /* box_11 */
    s(181, 177, 180,    00,    00);  /* box_12 */
    s(182,  54,  54,   +15,   +15);  /* box_13 */
    s(183,  59,  59,    -2,   -15);  /* box_14 */
    s(184, 337, 324,   +47,   +42);  /* box_15 */ 
    s(186,  29,  29,   +18,   -25);  /* box_16 */ 
    
    /* Pipe ends and joints: */
    s(191,  40,  40,   +57,  +119);
    s(192,  40,  40,   +52,   +75);

    s(195,  45, 174,  -247,    00);
    s(194, 195, 173,    00,    00);
    s(193, 195,  41,    00,   +74);

    /* s(193,  41,  41,  +113,   +74); */ /* Redundant. */

    s(196, 158, 158,   -40,   +20);
    s(197, 158, 158,   +11,   -21);
    s(198, 177, 179,    00,    00);
    s(199, 177, 182,    00,    00);
    s(200, 177, 186,    00,    00);
    s(201, 183, 186,    00,    00);
    s(202, 192, 192,    00,   -30);
    s(203,  41,  41,   +17,   +74);
    s(204, 176, 176,   +85,    00);
    s(190, 184, 186,    00,    00);                                                                                 
    
    adrw_point_vec_trim(&P, np);
    show_holes(&P, 0, P.ne-1);
    
    return P;
  }

