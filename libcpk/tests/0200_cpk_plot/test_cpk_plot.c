#define PROG_NAME "test_cpk_plot"
#define PROG_DESC "tests the {libcpk} plotting routines"
#define PROG_VERSION "1.0"

/* Last edited on 2024-12-31 20:57:51 by stolfi */ 

#define PROG_C_COPYRIGHT "  Run \"" PROG_NAME " -info\" for details"

#define PROG_HELP \
#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    [ -mapName {MAP_NAME} ]"

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
  "AUTHOR\n" \
  "  Created 2024-12-31 by Jorge Stolfi, UNICAMP.\n" \
  "\n" \
  "MODIFICATION HISTORY\n" \
  "  2024-12-31 Created based on {test_cpk.c}. J. Stolfi.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " PROG_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS
  
#define PROG_INFO_DESC \
  "  The program reads a map description and draws it on an regular hexagonal grid."
  
#define PROG_INFO_OPTS \
  "  -mapName {MAP_NAME}\n" \
  "    This mandatory option specifies the name of the map, that will be read from \"in/{MAP_NAME}-urb.txt\"."

#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include <affirm.h>
#include <fget.h>
#include <r2.h>
#include <i2.h>
#include <interval.h>

#include <cpk_eps.h>
 
#include <cpk_basic.h>
#include <cpk_debug.h>
#include <cpk_io.h>
#include <cpk_coords.h>

typedef struct options_t
  { char *mapName;    /* Name of map to read. */
  } options_t;

/* INTERNAL PROTOTYPES */

int32_t main(int32_t argc, char **argv);
  
options_t* parse_options(int32_t argc, char **argv);
  /* Parses the command line options. */
   
r2_vec_t *tcpp__get_domain(char *mapName);
  /* Returns a {r2_vec_t} object with a set of polygons read from the
    file "in/{mapName}-urb.txt". Coordinates in the files should be
    Lon/Lat pairs (in fractional degrees). */

r2_vec_t tcpp__all_grid_points(r2_t org, double sx, double sy, uint32_t nx, uint32_t ny);
  /* Generates all points of an hexagonal grid with lower
    left pixel at {org}, {nx} points with step {sx} within each row,
    {ny} rows with step {sy} between rows. Odd-numbered rows
    are offset by {sx/2} in the X direction. */

/* IMPLEMENTATIONS */

int32_t main(int32_t argc, char **argv)
  {   
    /* Get/fake the miscellaneous options: */
    options_t *o = tcpp__parse_options(argc, argv);

    /* Read the geography data: */
    r2_vec_t C = tcpp__get_domain(o->mapName);

    /* Seed for random numbers: */
    uint32_t seed = 4635*417;
    
    cpk_eps_eps_fill_polygon(eps, &C);

    r2_vec_D = ???
    cpk_eps_eps_draw_polyline(eps, &D);

    r2_t ctr1 = (r2_t){{ ???, ??? }};
    double rad1 = ???;
    cpk_eps_eps_plot_crcle(eps, &ctr1, rad1, TRUE, TRUE)
    
    r2_t ctr2 = (r2_t){{ ???, ??? }};
    double rad2 = ???;
    cpk_eps_eps_plot_crcle(eps, &ctr2, rad2, TRUE, TRUE)
    
    r2_t a3 = ;
    r2_t b3 = ;
    double rad3 = 
    cpk_eps_eps_plot_rect_stroke(eps, a3, b3, rad3, TRUE, TRUE);
    
    /* Create the EPS writer: */
    epswr_figure_t *cpk_eps_eps_new_figure(interval_t B[], char *fname);
    
    return 0;
  }
  
r2_vec_t *tcpp__get_domain(char *mapName)
  {
    r2_vec_t *C;
    C = (r2_vec_t *)notnull(malloc(sizeof(r2_vec_t)), "no mem");
    /* Read main polygon: */
    C->Urb = tcpp__read_polys(mapName, "urb", TRUE);
    return C;
  }

r2_vec_t tcpp__read_polys(char *mapName, char *tag, bool_t closed)
  {
    char *path = jsprintf("in/%s-%s.txt", mapName, tag);
    FILE *rd = open_read(path, TRUE);
    free(path);
    
    /* Read file. */
    r2_vec_t P = r2_vec_new(100);
    uint32_t nP = 0;
    int32_t fP = -1; /* First point of current component, or -1 outside components. */
    while(1)
      { /* Skip blanks on same line:: */
        fget_skip_spaces(rd);
        /* Check next character: */
        int32_t c = fgetc(rd);
        if (c == EOF)
          { break; }
        else if (c == '\n')
          { /* Blank line, ignore: */ }
        else if (c == '\015')
          { /* CR: ignore, like LF, but eat next LF if any: */
            c = fgetc(rd);
            if (c == EOF) { break; }
            if (c != '\n') { ungetc(c, rd); }
          }
        else if 
          ( (c == '+') || (c == '-') || (c == '.') || 
            ((c >= '0') && (c <= '9'))
          )
          { /* Number, must be a point: */
            ungetc(c, rd);
            double Lon = fget_double(rd);
            double Lat = fget_double(rd);
            r2_t p = (r2_t){{Lon,Lat}};
            /* Rest of line must be empty: */
            fget_eol(rd);
            /* Check component status: */
            if (fP == -1)
              { /* First point of a component. */
                if (nP > 0)
                  { /* Not the first component, must insert {{#}-vertex: */
                    r2_vec_expand(&P, (vec_index_t)nP);
                    P.e[nP] = (r2_t){{INF,INF}}; 
                    nP++;
                  }
                /* Remember this point: */
                fP = (int32_t)nP;
              }
            else if (closed)
              { /* Check for component closure: */
                if (r2_eq(&p, &(P.e[fP])))
                  { /* End of component: */
                    fP = -1; 
                  }
              }
            /* Save point in list: */
            r2_vec_expand(&P, (vec_index_t)nP);
            P.e[nP] = p; 
            nP++;
          }
        else 
          { /* Non-number, non-blank: assume it is a component separator: */
            if ((fP != -1) && closed) 
              { demand(FALSE, "component is not closed"); }
            fP = -1;
            do { c = fgetc(rd); } while ((c != '\n') && (c != EOF));
            if (c == EOF) { break; }
          }
      }
    /* End of input, check for proper closure: */
    if ((fP != -1) && closed)
      { demand(FALSE, "component is not closed"); }
    /* Trim: */
    r2_vec_trim(&P, nP);
    return P;
  }
  
options_t* tcpp__parse_options(int32_t argc, char **argv)
  { 
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);
    argparser_process_help_info_options(pp);

    options_t *o = talloc(1, options_t);

    argparser_get_keyword(pp, "-mapName");
    o->mapName = argparser_get_next_non_keyword(pp);

    argparser_finish(pp);

    return o;
  }
