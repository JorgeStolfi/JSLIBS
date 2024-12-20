/* See colorfield.h
** Last edited on 2008-05-25 03:22:48 by stolfi
**
** Copyright (C) 2003 by Jorge Stolfi, the University of Campinas, Brazil.
** See the rights and conditions notice at the end of this file.
*/

#include <colorfield.h>

#include <limits.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include <argparser.h>
#include <bool.h>

#include <frgb_ops.h>
#include <colorfield_unif.h>
#include <colorfield_wavy.h>
#include <colorfield_ramp.h>

cfld_params_t *cfld_compute_params
  ( cfld_args_t *fargs, 
    frgb_adjuster_t *adjust,
    bool_t logarithmic
  )
  { cfld_params_t *fp = malloc(sizeof(cfld_params_t));
    fp->kind = fargs->kind;
    fp->logarithmic = logarithmic;
    fp->params = NULL;
    switch (fargs->kind)
      { case cfld_UNIF:
          { cfld_unif_args_t *args = (cfld_unif_args_t *)(fargs->args);
            cfld_unif_params_t *params =  cfld_unif_compute_params
              (args, adjust, logarithmic);
            fp->params = (void *)params;
          }
          break;
        case cfld_RAMP:
          { cfld_ramp_args_t *args = (cfld_ramp_args_t *)(fargs->args);
            cfld_ramp_params_t *params = cfld_ramp_compute_params
              (args, adjust, logarithmic);
            fp->params = (void *)params;
          }
          break;
        case cfld_WAVY:
          { cfld_wavy_args_t *args = (cfld_wavy_args_t *)(fargs->args);
            cfld_wavy_params_t *params = cfld_wavy_compute_params
              (args, adjust, logarithmic);
            fp->params = (void *)params;
          }
          break;
      }
    return fp;
  }

void cfld_eval(cfld_params_t *fp, int32_t col, int32_t row, frgb_t *fv, int32_t chns)
  { switch(fp->kind)
      { case cfld_UNIF:
          { cfld_unif_params_t *params = (cfld_unif_params_t *)(fp->params);
            cfld_unif_eval(params, fp->logarithmic, col, row, fv, chns);
          }
          break;
        case cfld_RAMP:
          { cfld_ramp_params_t *params = (cfld_ramp_params_t *)(fp->params);
            cfld_ramp_eval(params, fp->logarithmic, col, row, fv, chns);
          }
          break;

        case cfld_WAVY:
          { cfld_wavy_params_t *params = (cfld_wavy_params_t *)(fp->params);
            cfld_wavy_eval(params, fp->logarithmic, col, row, fv, chns);
          }
          break;
      }
  }

cfld_args_t *cfld_parse(argparser_t *pp)
  { cfld_args_t *fargs = (cfld_args_t *)malloc(sizeof(cfld_args_t));
    fargs->args = NULL;
    if (argparser_keyword_present_next(pp, "uniform"))
      { fargs->kind = cfld_UNIF;
        fargs->args = cfld_unif_parse(pp);
      }
    else if (argparser_keyword_present_next(pp, "ramp"))
      { fargs->kind = cfld_RAMP;
        fargs->args = (void *)cfld_ramp_parse_general(pp);
      }
    else if (argparser_keyword_present_next(pp, "hWave"))
      { fargs->kind = cfld_WAVY;
        fargs->args = (void *)cfld_wavy_parse_simple(pp, 1, 0 );
      }
    else if (argparser_keyword_present_next(pp, "vWave"))
      { fargs->kind = cfld_WAVY;
        fargs->args = (void *)cfld_wavy_parse_simple(pp, 0, 1 ); 
      }
    else if (argparser_keyword_present_next(pp, "wave"))
      { fargs->kind = cfld_WAVY;
        fargs->args = (void *)cfld_wavy_parse_general(pp, 1 ); 
      }
    else if (argparser_keyword_present_next(pp, "wavePair"))
      { fargs->kind = cfld_WAVY;
        fargs->args = (void *)cfld_wavy_parse_general(pp, 2 ); }
    else
      { argparser_error(pp, "invalid or missing field kind"); }
    return fargs;
  }
  
cfld_args_t *cfld_make_args_uniform(frgb_t *color)
  {
    cfld_args_t *fargs = (cfld_args_t *)malloc(sizeof(cfld_args_t));
    fargs->kind = cfld_UNIF;
    fargs->args = cfld_unif_make_args(color);
    return fargs;
  }

/*
** AUTHORSHIP, RIGHTS AND CONDITIONS NOTICE.  This program was written
** in apr/2003 by Jorge Stolfi at the University of Campinas, Brazil.
** 
** Permission to use, copy, modify, and redistribute this software and
** its documentation for any purpose and without fee is hereby
** granted, provided that: (1) the copyright notice at the top of this
** file and this authorship, rights and conditions notice is retained
** in all derived source files and documentation; (2) no executable
** code derived from this file is published or distributed without the
** corresponding source code; and (3) these same rights are granted to
** any recipient of such code, under the same conditions.  This
** software is provided "as is" without express or implied warranty.
** Neither the author nor his employer are liable for any damages or
** losses that may result from its use. END OF NOTICE.
**
*/
 
