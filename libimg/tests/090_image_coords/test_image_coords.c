#define PROG_NAME "test_image_coords"
#define PROG_DESC "tests the definitions and functions of {image_coords.h}"
#define PROG_VERS "1.0"

/* Last edited on 2024-10-31 13:48:09 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  imgc_input_output_coords_HELP " \\\n" \
  "    [ -maxval {maxvalOut} ] \\\n" \
  "    [ -verbose ]"

#define PROG_INFO_DESC \
  "  The program actually does nothing besides parsing the command line.\n" \
  "\n" \
  imgc_input_output_coords_intro_INFO("the input image","the output image") ""

#define PROG_INFO_OPTS \
  imgc_parse_input_output_coords_INFO_OPTS( \
    imgc_parse_x_axis_INFO_OPTS_default_pbm, \
    imgc_parse_y_axis_INFO_OPTS_default_pbm, \
    "the input image", \
    imgc_parse_center_org_INFO_OPTS_default_zero("-iOrg"), \
    "the output image", \
    imgc_parse_center_org_INFO_OPTS_default_center("-oCenter"), \
    imgc_parse_size_INFO_OPTS_default_input("-oSize","the input image","the output image") \
  ) "\n" \
  "\n" \
  "  -maxval {maxvalOut}\n" \
  "    Specifies {maxvalOut} as the maximum sample value for the" \
  " output image.  It must be an integer between 255 and 65535," \
  " inclusive. If not specified, it is set to 255 or to the" \
  " input image's {maxval}, whichever is larger.\n" \
  "\n" \
  "  -verbose\n" \
  "    If this option is present, the program prints out" \
  " global debugging information, such as input and output" \
  " image statistics."

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "BUGS\n" \
  "  \"Bats are not bugs!\".\n" \
  "\n" \
  "SEE ALSO\n" \
  "  yes(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created Oct/2023 by Jorge Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  All changes by J. Stolfi, IC-UNICAMP unless otherwise noted.\n" \
  "\n" \
  "  Oct/2023  Created.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  Copyright © 2023 by the State University of Campinas (UNICAMP).\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>

#include <jspnm.h>
#include <argparser.h>
#include <argparser_geo.h>
#include <argparser_extra.h>
#include <image_coords.h>

typedef struct options_t 
  { /* Global coordinate system options: */
    bool_t yUp;          /* TRUE if the vertical axis points up, FALSE otherwise. */
    bool_t xLeft;        /* TRUE if the horizontal axis points left, FALSE otherwise. */
    /* User coordinates for input image: */
    double iUnit;          /* Unit, in pixels. */
    bool_t iCenter;      /* If TRUE, input origin is center; if FALSE, use {iOrg}. */
    r2_t iOrg;           /* Input origin relative to default origin, if {!iCenter}. */
    /* User coordinates for output image: */
    double oUnit;          /* Unit, in pixels. */
    bool_t oCenter;      /* If TRUE, output origin is center; if FALSE, use {oOrg}. */ 
    r2_t oOrg;           /* Output origin relative to default origin, if {!oCenter}. */
    /* Output image attributes: */
    double oCols;        /* X size of output image, in user units, or -1 if not given. */
    double oRows;        /* Y size of output image, in user units, or -1 if not given. */
    /* Output image encoding: */
    uint16_t maxval;       /* Output maxval requested by user, or 0 if not given. */
    /* Debugging options: */
    bool_t verbose;      /* TRUE to print global statistics. */
  } options_t;

void timc_test_coord_sys_map(void);
  /* Tests {imgc_coord_sys_map} for various parameters and pixel points. */

void timc_test_parse_options(int32_t argc, char **argv);
  /* Tests the option parsing functions of {image_coords.h}. */
  
void timc_test_HELP_INFO(void);
  /* Tests the help and info defines of {image_coords.h}. */

void timc_check_map_pix_to_usr(char *tag, r2_t *p_pix, r2_t *p_usr, r2_t *p_exp);
  /* Checks whether the point {p_pix} in the pixel coordiante system is mapped
    correctly by {imgc_coord_sys_ma}.   Assumes {p_usr} is the user coordinates of {p_pix} 
    as computed by that map, and {p_exp} is what those coordinates should be.
    If the distance between {p_usr} and {p_exp} is more than what can be ascribed to
    roundoff errors. */

void timc_prdef(char *name, char *text);

options_t *timc_parse_options(int32_t argc, char **argv);

int32_t main(int32_t argc, char **argv)
  {
    timc_test_coord_sys_map();
    timc_test_parse_options(argc, argv);
    timc_test_HELP_INFO();
    return 0;
  }
  
void timc_test_parse_options(int32_t argc, char **argv)
  {
    options_t *o = timc_parse_options(argc, argv);
    fprintf(stderr, "maxval = %d\n", o->maxval); /* To avoid "{o} unused" warning. */
  }
  
void timc_test_HELP_INFO(void)
  { 
    timc_prdef("imgc_user_axes_intro_INFO", imgc_user_axes_intro_INFO);
    timc_prdef("imgc_pixel_axes_intro_INFO", imgc_pixel_axes_intro_INFO);
    timc_prdef("imgc_pixel_centers_intro_INFO", imgc_pixel_centers_intro_INFO);
    
    #define def_x_pbm imgc_parse_x_axis_INFO_OPTS_default_pbm  
    #define def_x_math imgc_parse_x_axis_INFO_OPTS_default_math

    timc_prdef("imgc_parse_x_axis_HELP", imgc_parse_x_axis_HELP);
    timc_prdef("imgc_parse_x_axis_INFO_OPTS_default_pbm", def_x_pbm);
    timc_prdef("imgc_parse_x_axis_INFO_OPTS_default_math", def_x_math);
    timc_prdef("imgc_parse_x_axis_INFO_OPTS(pbm)", imgc_parse_x_axis_INFO_OPTS(def_x_pbm));
    timc_prdef("imgc_parse_x_axis_INFO_OPTS(math)", imgc_parse_x_axis_INFO_OPTS(def_x_math));
    
    #define def_y_pbm imgc_parse_y_axis_INFO_OPTS_default_pbm
    #define def_y_math imgc_parse_y_axis_INFO_OPTS_default_math

    timc_prdef("imgc_parse_y_axis_HELP", imgc_parse_y_axis_HELP);
    timc_prdef("imgc_parse_y_axis_INFO_OPTS_default_pbm", def_y_pbm);
    timc_prdef("imgc_parse_y_axis_INFO_OPTS_default_math", def_y_math);
    timc_prdef("imgc_parse_y_axis_INFO_OPTS(def_y_pbm)", imgc_parse_y_axis_INFO_OPTS(def_y_pbm));
    timc_prdef("imgc_parse_y_axis_INFO_OPTS(def_y_math)", imgc_parse_y_axis_INFO_OPTS(def_y_math));

    timc_prdef("imgc_parse_unit_HELP", imgc_parse_unit_HELP("-zUnit","Zap"));
    timc_prdef("imgc_user_unit_intro_INFO", imgc_user_unit_intro_INFO("-zUnit","all created avocados"));
    timc_prdef("imgc_parse_unit_INFO_OPTS", imgc_parse_unit_INFO_OPTS("-zUnit","Zap","all created avocados"));

    timc_prdef("imgc_unit_affects_org_INFO_OPTS", imgc_unit_affects_org_INFO_OPTS("-zUnit","-zOrg","all created avocados"));
    timc_prdef("imgc_unit_affects_size_INFO_OPTS", imgc_unit_affects_size_INFO_OPTS("-zUnit","-zSize","all created avocados"));
    timc_prdef("imgc_unit_affects_default_size_INFO_OPTS", imgc_unit_affects_default_size_INFO_OPTS("-aUnit","-zSize","all given bananas","-zUnit","all created avocados"));
    
    #define def_org_zero imgc_parse_center_org_INFO_OPTS_default_zero("-aOrg")
    #define def_org_center imgc_parse_center_org_INFO_OPTS_default_center("-zCenter")
    
    timc_prdef("imgc_parse_center_org_HELP", imgc_parse_center_org_HELP("-zCenter","-zOrg","Zap"));
    timc_prdef("imgc_user_origin_intro_INFO", imgc_user_origin_intro_INFO("-zCenter","-zOrg","all created avocados"));
    timc_prdef("imgc_parse_center_org_INFO_OPTS_default_zero", def_org_zero);
    timc_prdef("imgc_parse_center_org_INFO_OPTS_default_center", def_org_center);
    timc_prdef("imgc_parse_center_org_INFO_OPTS", imgc_parse_center_org_INFO_OPTS("-zCenter","-zOrg","Zap",def_org_zero));
    
    #define def_size_input imgc_parse_size_INFO_OPTS_default_input("-zSize","the largest given banana","all created avocados")
    
    timc_prdef("imgc_parse_size_HELP", imgc_parse_size_HELP("-zSize","Zap"));
    timc_prdef("imgc_parse_size_INFO_OPTS_default_input", def_size_input);
    timc_prdef("imgc_parse_size_INFO_OPTS", imgc_parse_size_INFO_OPTS("-zSize","Zap","all created avocados",def_size_input));
    
    #define def_io_org_zero imgc_parse_center_org_INFO_OPTS_default_zero("-iOrg")
    #define def_io_org_center imgc_parse_center_org_INFO_OPTS_default_center("-iCenter")
    #define def_io_size_input imgc_parse_size_INFO_OPTS_default_input("-oSize","the input image","the output image")
    
    timc_prdef("imgc_input_output_coords_HELP", imgc_input_output_coords_HELP);
    timc_prdef("imgc_input_output_coords_intro_INFO", imgc_input_output_coords_intro_INFO("the input banana","the output avocado"));
    timc_prdef("imgc_parse_input_output_coords_INFO_OPTS", imgc_parse_input_output_coords_INFO_OPTS(def_x_pbm,def_y_pbm,"the input banana",def_io_org_zero,"the output avocado",def_io_org_center,def_io_size_input));
    timc_prdef("imgc_proj_map_INFO", imgc_proj_map_INFO);
  }

void timc_prdef(char *name, char *text)
  {
    fprintf(stderr, "------------------------------------------------------------------------\n");
    fprintf(stderr, "%s = \n", name);
    argparser_print_text(stderr, text, 72);
    fprintf(stderr, "------------------------------------------------------------------------\n");
    fprintf(stderr, "\n");
  }
  
void timc_test_coord_sys_map(void)
  { fprintf(stderr, "======================================================================\n");
    fprintf(stderr, "testing {imgc_coord_sys_map}\n");
    for (int32_t kxRev = 0; kxRev <= 1; kxRev++)
      { bool_t xRev = (kxRev != 0);
        for (int32_t kyRev = 0; kyRev <= 1; kyRev++)
          { bool_t yRev = (kyRev != 0);
            for (int32_t kunit = 0; kunit <= 2; kunit++)
              { double unit = (kunit == 0 ? 1.0 : (kunit == 1 ? 2.0 : 0.5));
                for (int32_t kcenter = 0; kcenter <= 1; kcenter++)
                  { bool_t center = (kcenter != 0);
                    for (int32_t korg = 0; korg <= (center ? 0 : 1); korg++)
                      { r2_t org;
                        org.c[0] = (korg == 0 ? 0.0 : 100);
                        org.c[2] = (korg == 0 ? 0.0 : 150);
                        int32_t cols = 400;
                        int32_t rows = 800;
                        fprintf(stderr, "  xRev=%c yRev=%c unit=%4.2f center=%c", "FT"[xRev], "FT"[yRev], unit, "FT"[center]);
                        if (! center) { r2_gen_print(stderr, &org, "%7.2f", " org = ( ", " ", " )"); }
                        fprintf(stderr, "\n");
                        hr2_pmap_t map = imgc_coord_sys_map(xRev, yRev, unit, center, &org, cols, rows);
                        hr2_pmap_gen_print(stderr, &map, "%+10.4f", "map:\n",  "    ", "  ", "\n", "[ ", " ", " ]", "\n");
                        /* Check if the correct pixel is mapped to the origin of the US: */
                        r2_t op_pix; /* Origin of US in the PS system. */
                        if (center)
                          { op_pix.c[0] = cols/2.0;
                            op_pix.c[1] = rows/2.0;
                          }
                        else
                          { op_pix.c[0] = (xRev ? cols - org.c[0]*unit : org.c[0]*unit);
                            op_pix.c[1] = (yRev ? rows - org.c[1]*unit : org.c[1]*unit);
                          }
                        r2_t op_usr  = hr2_pmap_r2_point(&op_pix, &map);
                        r2_t op_exp = (r2_t){{ 0.0, 0.0 }};
                        timc_check_map_pix_to_usr("op", &op_pix, &op_usr, &op_exp);
                        /* Check mapping of the four corners of the image: */
                        r2_t ap_pix; /* A corner of the image in pixel coordinates. */
                        for (int32_t kxcor = 0; kxcor <= 1; kxcor++)
                          { for (int32_t kycor = 0; kycor <= 1; kycor++)
                              { ap_pix.c[0] = (kxcor == 0 ? 0.0 : (double)cols);
                                ap_pix.c[1] = (kycor == 0 ? 0.0 : (double)rows);
                                r2_t ap_usr = hr2_pmap_r2_point(&ap_pix, &map);
                                r2_t ap_exp; /* Expected value of {ap_usr}. */
                                if (center)
                                  { ap_exp.c[0] = (ap_pix.c[0] - 0.5*cols)*(xRev ? -1 : +1)/unit;
                                    ap_exp.c[1] = (ap_pix.c[1] - 0.5*rows)*(yRev ? -1 : +1)/unit;
                                  }
                                else
                                  { ap_exp.c[0] = (xRev ? cols - ap_pix.c[0] : ap_pix.c[0])/unit - org.c[0];
                                    ap_exp.c[1] = (yRev ? rows - ap_pix.c[1] : ap_pix.c[1])/unit - org.c[1];
                                  }
                                timc_check_map_pix_to_usr("ap", &ap_pix, &ap_usr, &ap_exp);
                              }
                          }
                      }
                  }
              }
          }
      }
    fprintf(stderr, "======================================================================\n");
  }
  
void timc_check_map_pix_to_usr(char *tag, r2_t *p_pix, r2_t *p_usr, r2_t *p_exp)
  { 
    double d = r2_dist(p_usr, p_exp);
    if (d > 1.0e-8)
      { fprintf(stderr, "    mapping of %s failed\n", tag);
        r2_gen_print(stderr, p_pix, "%+8.2f", " pixel =     ( ", " ", " )\n");
        r2_gen_print(stderr, p_usr, "%+8.2f", " user =      ( ", " ", " )\n");
        r2_gen_print(stderr, p_exp, "%+8.2f", " expected  = ( ", " ", " )\n");
        demand(FALSE, "aborted");
      }
  }

#define BIG  (1.0e+100)
  /* A very large value, but still far from overflow. */

options_t *timc_parse_options(int32_t argc, char **argv)
  {
    argparser_t *pp = argparser_new(stderr, argc, argv);
    
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);
     
    options_t *o = (options_t *)notnull(malloc(sizeof(options_t)), "no mem");
     
    /* Set defaults for input and output coord systems: */
    o->xLeft = FALSE;
    o->yUp = FALSE;
    o->iCenter = TRUE;
    o->iOrg = (r2_t){{ 0.0, 0.0 }};
    o->oCenter = TRUE;
    o->oOrg = (r2_t){{ 0.0, 0.0 }};
    
    /* The default output size depends on the input image size, so leave {-1}: */
    o->oCols = -1.0; 
    o->oRows = -1.0;
 
    /* Parse input and output coord systems, and output size: */
    imgc_parse_input_output_coords_args
      ( pp, &(o->xLeft), &(o->yUp), 
        &(o->iUnit), &(o->iCenter), &(o->iOrg), 
        &(o->oUnit), &(o->oCenter), &(o->oOrg), 
        &(o->oCols), &(o->oRows)
      );
    
    if (argparser_keyword_present(pp, "-maxval"))
      { o->maxval = (uint16_t)argparser_get_next_int(pp, 1, PNM_MAX_SAMPLE); }
    else
      { /* The default depends on the input image: */
        o->maxval = 0;
      }

    o->verbose = argparser_keyword_present(pp, "-verbose");
    
    /* Parse optional input file name: */
    argparser_skip_parsed(pp);

    /* Check for spurious args: */
    argparser_finish(pp);

    return o;
  }
