/* See ellipse_crs_args_parse.h */
/* Last edited on 2021-06-09 19:48:45 by jstolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <values.h>

#include <r2.h>
#include <affirm.h>
#include <argparser.h>
#include <argparser_geo.h>

#include <ellipse_crs.h>
#include <ellipse_crs_args.h>

void ellipse_crs_args_parse
  ( argparser_t *pp,   /* The command line parsing state. */
    ellipse_crs_t *EP, /* (OUT) The ellipse parameters. */
    double *ctrAdj,    /* (OUT) The max adjustment to {EP.ctr} coords, or NULL. */
    double *radAdj,    /* (OUT) The max adjustment to {EP.rad}, or NULL. */ 
    double *strAdj     /* (OUT) The max adjustment to {EP.str} coords, or NULL. */
  )
  { 
    /* Provide NAN defaults for geometric args: */
    EP->ctr = (r2_t){{ NAN, NAN}};
    EP->rad = NAN; 
    EP->str = (r2_t){{ NAN, NAN }};
    
    /* Provide zero defaults for the adjustment amounts: */
    if (ctrAdj != NULL) { *ctrAdj = 0.0; }
    if (radAdj != NULL) { *radAdj = 0.0; }
    if (strAdj != NULL) { *strAdj = 0.0; }
    
    /* To check for duplicates: */
    bool_t ctr_given = FALSE;
    bool_t rad_given = FALSE;
    bool_t str_given = FALSE;
    while (TRUE)
      { if (argparser_keyword_present_next(pp, "center"))
          { if (ctr_given) { argparser_error(pp, "duplicate sphere's \"center\""); }
            if (argparser_next_is_number(pp))
              { EP->ctr.c[0] = argparser_get_next_double(pp, -1.0e+5, 1.0e+5);
                EP->ctr.c[1] = argparser_get_next_double(pp, -1.0e+5, 1.0e+5);
              }
            argparser_get_next_adjust(pp, ctrAdj, 0.0, 1.0e+5);
            ctr_given = TRUE;
          }
        else if (argparser_keyword_present_next(pp, "radius"))
          { if (rad_given) { argparser_error(pp, "duplicate sphere's \"radius\""); }
            if (argparser_next_is_number(pp))
              { EP->rad = argparser_get_next_double(pp, 0.01, 1.0e+5); }
            argparser_get_next_adjust(pp, radAdj, 0.0, 1.0e+5);
            rad_given = TRUE;
          }
        else if (argparser_keyword_present_next(pp, "stretch"))
          { if (str_given) { argparser_error(pp, "duplicate sphere's \"stretch\""); }
            if (argparser_next_is_number(pp))
              { EP->str.c[0] = argparser_get_next_double(pp, -1.0e+10, +1.0e+10);
                EP->str.c[1] = argparser_get_next_double(pp, -1.0e+10, +1.0e+10);
              }
            argparser_get_next_adjust(pp, strAdj, 0.0, 1.0e+5);
            str_given = TRUE;
          }
        else
          { /* The next arg is not a parameter keyword -- assume end of specs: */
            break;
          }
      }
  }

void ellipse_crs_args_print(FILE *wr, ellipse_crs_t *E, char *fmt)
  { ellipse_crs_args_adjust_print(wr, E, 0, 0, 0, fmt); }

void ellipse_crs_args_adjust_print
  ( FILE *wr, 
    ellipse_crs_t *E, 
    double ctrAdj, 
    double radAdj, 
    double strAdj,
    char *fmt
  )
  { /* Try to determine the natural width of {fmt}-printed values: */
    char *tmp = NULL;
    asprintf(&tmp, fmt, 0);
    int32_t wd = (int32_t)strlen(tmp);
    free(tmp);
    
    fprintf(wr, "  center ");
    if ((! isnan(E->ctr.c[0])) || (! isnan(E->ctr.c[1])))
      { fprintf(wr, " ");
        fprintf(wr, fmt, E->ctr.c[0]);
        fprintf(wr, " ");
        fprintf(wr, fmt, E->ctr.c[1]);
      }
    if ((! isnan(ctrAdj)) && (ctrAdj != 0))
      { fprintf(wr, " adjust ");
        fprintf(wr, fmt, ctrAdj);
      }
    fprintf(wr, "\n");
    
    fprintf(wr, "  radius ");
    if (! isnan(E->rad)) 
      { fprintf(wr, " ");
        fprintf(wr, fmt, E->rad); 
      }
    if ((! isnan(radAdj)) && (radAdj != 0))
      { fprintf(wr, " ");
        fprintf(wr, "%*s", wd, "");
        fprintf(wr, " adjust "); 
        fprintf(wr, fmt, radAdj);
      }
    fprintf(wr, "\n");
    
    fprintf(wr, "  stretch");
    if ((! isnan(E->str.c[0])) || (! isnan(E->str.c[1])))
      { fprintf(wr, " ");
        fprintf(wr, fmt, E->str.c[0]);
        fprintf(wr, " ");
        fprintf(wr, fmt, E->str.c[1]); 
      }
    if ((! isnan(strAdj)) && (strAdj != 0))
      { fprintf(wr, " adjust ");
        fprintf(wr, fmt, strAdj);
      }
    fprintf(wr, "\n");

    fflush(wr);
  }
