#define PROG_NAME "test_image_coords"
#define PROG_DESC "tests the definitions and functions of {image_coords.h}"
#define PROG_VERS "1.0"

/* Last edited on 2023-10-14 10:59:05 by stolfi */

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  imgc_input_output_coords_HELP " \\\n" \
  "    [ -maxval {MV_OUT} ] \\\n" \
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
  "  -maxval {MV_OUT}\n" \
  "    Specifies {MV_OUT} as the maximum sample value for the" \
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

void timc_prdef(char *name, char *text);
options_t *timc_parse_options(int32_t argc, char **argv);

int32_t main(int32_t argc, char **argv)
  {
    options_t *o = timc_parse_options(argc, argv);
    
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

    timc_prdef("imgc_parse_unit_HELP", imgc_parse_unit_HELP("-zUnit","_ZAP"));
    timc_prdef("imgc_user_unit_intro_INFO", imgc_user_unit_intro_INFO("-zUnit","all created avocados"));
    timc_prdef("imgc_parse_unit_INFO_OPTS", imgc_parse_unit_INFO_OPTS("-zUnit","_ZAP","all created avocados"));

    timc_prdef("imgc_unit_affects_org_INFO_OPTS", imgc_unit_affects_org_INFO_OPTS("-zUnit","-zOrg","all created avocados"));
    timc_prdef("imgc_unit_affects_size_INFO_OPTS", imgc_unit_affects_size_INFO_OPTS("-zUnit","-zSize","all created avocados"));
    timc_prdef("imgc_unit_affects_default_size_INFO_OPTS", imgc_unit_affects_default_size_INFO_OPTS("-aUnit","-zSize","all given bananas","-zUnit","all created avocados"));
    
    #define def_org_zero imgc_parse_center_org_INFO_OPTS_default_zero("-aOrg")
    #define def_org_center imgc_parse_center_org_INFO_OPTS_default_center("-zCenter")
    
    timc_prdef("imgc_parse_center_org_HELP", imgc_parse_center_org_HELP("-zCenter","-zOrg","_ZAP"));
    timc_prdef("imgc_user_origin_intro_INFO", imgc_user_origin_intro_INFO("-zCenter","-zOrg","all created avocados"));
    timc_prdef("imgc_parse_center_org_INFO_OPTS_default_zero", def_org_zero);
    timc_prdef("imgc_parse_center_org_INFO_OPTS_default_center", def_org_center);
    timc_prdef("imgc_parse_center_org_INFO_OPTS", imgc_parse_center_org_INFO_OPTS("-zCenter","-zOrg","_ZAP",def_org_zero));
    
    #define def_size_input imgc_parse_size_INFO_OPTS_default_input("-zSize","the largest given banana","all created avocados")
    
    timc_prdef("imgc_parse_size_HELP", imgc_parse_size_HELP("-zSize","_ZAP"));
    timc_prdef("imgc_parse_size_INFO_OPTS_default_input", def_size_input);
    timc_prdef("imgc_parse_size_INFO_OPTS", imgc_parse_size_INFO_OPTS("-zSize","_ZAP","all created avocados",def_size_input));
    
    #define def_io_org_zero imgc_parse_center_org_INFO_OPTS_default_zero("-iOrg")
    #define def_io_org_center imgc_parse_center_org_INFO_OPTS_default_center("-iCenter")
    #define def_io_size_input imgc_parse_size_INFO_OPTS_default_input("-oSize","the input image","the output image")
    
    timc_prdef("imgc_input_output_coords_HELP", imgc_input_output_coords_HELP);
    timc_prdef("imgc_input_output_coords_intro_INFO", imgc_input_output_coords_intro_INFO("the input banana","the output avocado"));
    timc_prdef("imgc_parse_input_output_coords_INFO_OPTS", imgc_parse_input_output_coords_INFO_OPTS(def_x_pbm,def_y_pbm,"the input banana",def_io_org_zero,"the output avocado",def_io_org_center,def_io_size_input));
    timc_prdef("imgc_proj_map_INFO", imgc_proj_map_INFO);
    
    fprintf(stderr, "maxval = %d\n", o->maxval); /* To avoid "{o} unused" warnng. */

    return 0;
  }

void timc_prdef(char *name, char *text)
  {
    fprintf(stderr, "------------------------------------------------------------------------\n");
    fprintf(stderr, "%s = \n", name);
    argparser_print_text(stderr, text, 72);
    fprintf(stderr, "------------------------------------------------------------------------\n");
    fprintf(stderr, "\n");
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
