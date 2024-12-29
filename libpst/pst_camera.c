/* See pst_camera.h */
/* Last edited on 2024-12-28 06:36:51 by stolfi */ 

#include <math.h>
#include <values.h>
#include <assert.h>
#include <stdint.h>

#include <float_image_mscale.h>
#include <argparser.h>
#include <argparser_geo.h>
#include <r2.h>
#include <r3.h>
#include <hr3.h>
#include <affirm.h> 

#include <pst_camera.h>
#include <pst_basic.h>

#define X c[0]
#define Y c[1]
#define Z c[2]
  /* Cartesian coordinates of an {r2_t} or {r3_t}. */
  
#define hm c.c[0]
#define hx c.c[1]
#define hy c.c[2]
#define hz c.c[3]
  /* Homogeneous coordinates of an {hr2_point_t} or {hr3_point_t}. */

double pst_camera_focal_length(hr3_point_t *O)
  { assert(O->hm >= 0);
    assert(O->hz > 0);
    return O->hz/O->hm;
  }
  
double pst_camera_spread(hr3_point_t *O)  
  { assert(O->hm >= 0);
    assert(O->hz > 0);
    return O->hm/O->hz;
  }

hr2_point_t pst_camera_center(hr3_point_t *O)
  { assert(O->hm >= 0);
    return (hr2_point_t) {{{ O->hm, O->hx, O->hy }}}; 
  }

r2_t pst_camera_center_cartesian(hr3_point_t *O)
  { hr2_point_t Q = pst_camera_center(O);
    return r2_from_hr2(&Q);
  }

hr3_point_t pst_camera_viewpoint_from_center_spread(r2_t *Q, double G)
  { if (G == 0)
      { return (hr3_point_t) {{{ 0, 0, 0, 1 }}}; }
    else
      { return (hr3_point_t) {{{ 1, Q->X, Q->Y, 1/G }}}; }
  }

double pst_camera_min_focal_length(r2_t *Q, int32_t NX, int32_t NY)
  { /* Determine the max distance {dMax} from {(OX,OY)} to any image pixel: */
    double xMax = (Q == NULL ? NX : fmax(NX - Q->X, Q->X));
    double yMax = (Q == NULL ? NY : fmax(NY - Q->Y, Q->Y));
    double dMax = hypot(xMax, yMax);
    /* Choose a maximum angle {aMax} of any pixel from the optical axis: */
    double aMax = pst_camera_max_angle*(M_PI/180); /* In radians. */
    /* Return the corresponding focal length: */
    return dMax/tan(aMax); 
  }
 
pst_camera_t pst_camera_shrink(pst_camera_t *C, int32_t dx, int32_t dy, int32_t nw)
  { if (C->O.hm == 0)
      { /* Camera at infinity: */
        return (*C);
      }
    else
      { /* Must halve {F}, shrink and shift the center: */
        r2_t Q = pst_camera_center_cartesian(&(C->O));
        double F = pst_camera_focal_length(&(C->O));
        r2_t Q_r = float_image_mscale_point_shrink(&Q, dx, dy, nw); 
        hr3_point_t O_r = pst_camera_viewpoint_from_center_spread(&Q_r, F/2);
        return (pst_camera_t) { .O = O_r };
      }
  }

pst_camera_t pst_camera_expand(pst_camera_t *C, int32_t dx, int32_t dy, int32_t nw)
  { if (C->O.hm == 0)
      { /* Camera at infinity: */
        return (*C);
      }
    else
      { /* Must double {F}, expand and shift the center: */
        r2_t Q = pst_camera_center_cartesian(&(C->O));
        double F = pst_camera_focal_length(&(C->O));
        r2_t Q_e = float_image_mscale_point_expand(&(Q), dx, dy, nw); 
        hr3_point_t O_e = pst_camera_viewpoint_from_center_spread(&Q_e, F*2);
        return (pst_camera_t) { .O = O_e };
      }
  }

void pst_camera_args_parse
  ( argparser_t *pp,  /* The command line parsing state. */
    pst_camera_t *C,  /* (OUT) The camera parameters. */
    double *OAdj,   /* (OUT) The viewpoint adjustment amount, or NULL. */ 
    double *QAdj,   /* (OUT) The optical center adjustment amount, or NULL. */ 
    double *GAdj    /* (OUT) The spread adjustment amount, or NULL. */ 
  )
  { 
    /* Default is all {NAN}s: */
    C->O = (hr3_point_t){{{ NAN, NAN, NAN, NAN }}};
    
    /* Alternative model parsed from the command line: */
    r2_t Q = (r2_t){{ NAN, NAN }};
    double G = NAN;
    
    /* Provide zero defaults for the adjustment amounts: */
    if (OAdj != NULL) { *OAdj = 0.0; }
    if (QAdj != NULL) { *QAdj = 0.0; }
    if (GAdj != NULL) { *GAdj = 0.0; }
    
    /* To check for duplicates: */
    bool_t viewpoint_given = FALSE;
    bool_t center_given = FALSE;
    bool_t spread_given = FALSE;
    
    /* Max values allowed: */
    double hptMax = 1.0e+100; /* Homogeneous coordinate. */
    double cptMax = 1.0e+5;   /* Cartesian coordinate. */
    double sprMax = 1/pst_camera_min_focal_length(NULL, 1, 1);   /* Spread. */
    
    while (TRUE)
      { if (argparser_keyword_present_next(pp, "viewpoint"))
          { if (viewpoint_given) { argparser_error(pp, "duplicate \"viewpoint\""); }
            if (center_given || spread_given) 
              { argparser_error(pp, "\"viewpoint\" excludes \"center\",\"spread\""); }
            if (argparser_next_is_number(pp))
              { C->O.hm = argparser_get_next_double(pp, -hptMax, +hptMax);
                C->O.hx = argparser_get_next_double(pp, -hptMax, +hptMax);
                C->O.hy = argparser_get_next_double(pp, -hptMax, +hptMax);
                C->O.hz = argparser_get_next_double(pp, -hptMax, +hptMax);
                if (C->O.hm < 0) 
                  { argparser_error(pp, "negative viewpoint weight"); }
                if (C->O.hz <= 0)
                  { argparser_error(pp, "non-positive viewpoint Z"); }
                if (r4_norm(&(C->O.c)) == 0)
                  { argparser_error(pp, "indefinite viewpoint"); }
              }
            argparser_get_next_adjust(pp, OAdj, 0.0, hptMax );
            viewpoint_given = TRUE;
          }
        else if (argparser_keyword_present_next(pp, "center"))
          { if (center_given) { argparser_error(pp, "duplicate \"center\""); }
            if (viewpoint_given) 
              { argparser_error(pp, "\"viewpoint\" excludes \"center\",\"spread\""); }
            if (argparser_next_is_number(pp))
              { Q.X = argparser_get_next_double(pp, -cptMax, +cptMax);
                Q.Y = argparser_get_next_double(pp, -cptMax, +cptMax);
              }
            argparser_get_next_adjust(pp, QAdj, 0.0, cptMax );
            center_given = TRUE;
          }
        else if (argparser_keyword_present_next(pp, "spread"))
          { if (spread_given) { argparser_error(pp, "duplicate \"spread\""); }
            if (viewpoint_given) 
              { argparser_error(pp, "\"viewpoint\" excludes \"center\",\"spread\""); }
            if (argparser_next_is_number(pp))
              { G = argparser_get_next_double(pp, 0.0, sprMax); }
            argparser_get_next_adjust(pp, GAdj, 0.0, sprMax);
            spread_given = TRUE;
          }
        else
          { /* The next arg is not a parameter keyword -- assume end of specs: */
            break;
          }
      }
      
    /* If user opted for center/spread, pack them nto {C}: */
    if (center_given || spread_given)
      { assert(! viewpoint_given);
        /* Assemble the viewpoint from {Q,G}: */
        C->O = pst_camera_viewpoint_from_center_spread(&Q, G);
      }
  }

r3x3_t pst_camera_normal_correction_matrix(r3_t *dir)
  { r3x3_t M;
    double sR = hypot(dir->c[0],dir->c[1]);
    double cR = dir->c[2];
    if (fabs(sR) > 1.0e-14) 
      { r3x3_t A,R;
        double ux = dir->c[0]/sR;
        double uy = dir->c[1]/sR;
        A = (r3x3_t)
          { { { +ux,  0.0, -uy },
              { +uy,  0.0, +ux },
              { 0.0,  1.0, 0.0 } } };
        R = (r3x3_t)
          { { { +cR, -sR, 0.0 },
              { +sR, +cR, 0.0 },
              { 0.0, 0.0, 1.0 } } };
        r3x3_mul(&A, &R, &M);
        r3x3_mul_tr(&M, &A, &M);
      }else {
        r3x3_ident(&M);
      }
    return M;
  }

void pst_camera_print(FILE *wr, pst_camera_t *C, char *fmt)
  { 
    fprintf(wr, "{");
    for (uint32_t i = 0; i < 4; i++)
      { fprintf(wr, " ");
        fprintf(wr, fmt, C->O.c.c[i]);
      }
    fprintf(wr, " }");
    
    fflush(wr);
  }

void pst_camera_args_print(FILE *wr, pst_camera_t *C, char *fmtp, char *fmtr)
  { /* Choose the presentation format: */ 
    bool_t use_viewpoint = FALSE; /* Must use the "viewpoint" presentation. */
    if (isnan(C->O.hm)) { use_viewpoint = TRUE; }
    if ((C->O.hm == 0) && ((C->O.hx != 0) || (C->O.hy != 0))) { use_viewpoint = TRUE; }
    if (use_viewpoint)
      { pst_camera_args_adjust_print_viewpoint(wr, C, 0, fmtp, fmtr); }
    else
      { pst_camera_args_adjust_print_center_spread(wr, C, 0, 0, fmtp, fmtr); }
  }

void pst_camera_args_adjust_print_viewpoint
  ( FILE *wr, 
    pst_camera_t *C, 
    double OAdj, 
    char *fmtp,
    char *fmtd
  )
  { 
    fprintf(wr, "  viewpoint");
    double *Oc = C->O.c.c;
    for (uint32_t i = 0; i < 4; i++)
      { fprintf(wr, " ");
        fprintf(wr, (Oc[0] == 0 ? fmtd : fmtp), Oc[i]);
      }
    if ((! isnan(OAdj)) && (OAdj != 0))
      { fprintf(wr, " adjust ");
        fprintf(wr, (Oc[0] == 0 ? fmtd : fmtp), OAdj);
      }
    fprintf(wr, "\n");

    fflush(wr);
  }

void pst_camera_args_adjust_print_center_spread
  ( FILE *wr, 
    pst_camera_t *C, 
    double QAdj, 
    double GAdj, 
    char *fmtp,
    char *fmtr
  )
  { 
    r2_t Q = pst_camera_center_cartesian(&(C->O));
    double G = pst_camera_spread(&(C->O));

    fprintf(wr, "  center");
    if ((! isnan(Q.X)) || (! isnan(Q.Y)))
      { fprintf(wr, " ");
        fprintf(wr, fmtp, Q.X);
        fprintf(wr, " ");
        fprintf(wr, fmtp, Q.Y);
      }
    if ((! isnan(QAdj)) && (QAdj != 0))
      { fprintf(wr, " adjust ");
        fprintf(wr, fmtp, QAdj);
      }
    fprintf(wr, "\n");
    
    fprintf(wr, "  spread");
    if (! isnan(G)) 
      { fprintf(wr, " ");
        fprintf(wr, fmtr, G);
      }
    if ((! isnan(GAdj)) && (GAdj != 0))
      { fprintf(wr, " adjust ");
        fprintf(wr, fmtr, GAdj);
      }
    fprintf(wr, "\n");

    fflush(wr);
  }

