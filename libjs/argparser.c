/* See argparser.h. */
/* Last edited on 2024-11-22 02:17:45 by stolfi */

/* Copyright © 2003 Jorge Stolfi, Unicamp. See note at end of file. */
/* Based on Params.m3 by J.Stolfi, DEC-SRC, 1988.  */

#include <limits.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include <vec.h>
#include <affirm.h>
#include <jswsize.h>

#include <argparser.h>
#include <argparser_extra.h>

void argparser_check_all_parsed(argparser_t *pp, uint32_t num);
  /* Checks whether all arguments between {argv[1]} and {argv[num-1]}
    have been parsed.  Fails with a message if any wasn't. */

/* IMPLEMENTATIONS */

argparser_t *argparser_new(FILE *wr, int32_t argc, char **argv)
  { argparser_t *pp = (argparser_t *)notnull(malloc(sizeof(argparser_t)), "no mem");
    pp->arg.ne = (uint32_t)argc;
    pp->arg.e = argv;
    pp->wr = wr;
    pp->parsed = bool_vec_new((uint32_t)argc);
    pp->parsed.e[0] = TRUE;
    for (int32_t i = 1; i < argc; i++) { pp->parsed.e[i] = FALSE; }
    pp->next = 1;
    pp->nhelp = 0;
    pp->help = string_vec_new(10);
    pp->ninfo = 0;
    pp->info = string_vec_new(10);
    return pp;
  }

void argparser_set_help(argparser_t *pp, char *help)
  { string_vec_expand(&(pp->help), pp->nhelp);
    pp->help.e[pp->nhelp] = help;
    pp->nhelp++;
  }

void argparser_set_info(argparser_t *pp, char *info)
  { string_vec_expand(&(pp->info), pp->ninfo);
    pp->info.e[pp->ninfo] = info;
    pp->ninfo++;
  }
    
void argparser_process_help_info_options(argparser_t *pp)
  { int32_t info_wd = 90; /* Line width for the INFO part. */
    if (argparser_keyword_present(pp, "-info"))
      { argparser_print_info(pp, info_wd); exit(0); }
    if (argparser_keyword_present(pp, "-help"))
      { argparser_print_help_and_halt(pp, 0); }
  }

void argparser_error(argparser_t *pp, char *msg)
  { argparser_error_at(pp, msg, "after", pp->next-1); }

void argparser_print_help_and_halt(argparser_t *pp, int32_t status)
  { for (int32_t i = 0; i < pp->nhelp; i++)
      { fprintf(pp->wr, "%s", pp->help.e[i]); }
    fprintf(pp->wr, "\n");
    exit(status);
  }

bool_t argparser_keyword_present(argparser_t *pp, char *key)
  { demand((key != NULL) && ((*key) != 0), "invalid null/empty key");
    char **a = pp->arg.e;
    bool_t *p = pp->parsed.e;
    for (int32_t i = 0; i < pp->arg.ne; i++)
      { if ((! p[i]) && argparser_key_matches(key, a[i]))
          { pp->next = i + 1;
            p[i] = TRUE;
            return TRUE;
          }
      }
    return FALSE;
  }

void argparser_get_keyword(argparser_t *pp, char *key)
  { if (! argparser_keyword_present(pp, key))
      { argparser_arg_msg(pp, "", 0, "keyword \"%s\" not found.\n", key);
        argparser_print_help_and_halt(pp, 1);
      }
  }

char *argparser_get_next(argparser_t *pp)
  { char **a = pp->arg.e;
    bool_t *p = pp->parsed.e;
    if ((pp->next >= pp->arg.ne) || (p[pp->next]))
      { argparser_arg_msg(pp, "missing argument after ", pp->next-1, "%s.\n", ""); 
        argparser_print_help_and_halt(pp, 1); 
      }
    p[pp->next] = TRUE;
    pp->next++;
    return a[pp->next-1];
  }

char *argparser_get_next_keyword(argparser_t *pp)
  { if (! argparser_next_is_keyword(pp))
      { argparser_arg_msg(pp, "missing keyword after ", pp->next-1, "%s.\n", ""); 
        argparser_print_help_and_halt(pp, 1);
        return NULL;
      }
    else
      { return argparser_get_next(pp); }
  }

char *argparser_get_next_non_keyword(argparser_t *pp)
  { if (! argparser_next_is_non_keyword(pp))
      { argparser_arg_msg(pp, "missing argument after ", pp->next-1, "%s.\n", ""); 
        argparser_print_help_and_halt(pp, 1);
        return NULL;
      }
    else
      { return argparser_get_next(pp); }
  }

void argparser_get_keyword_next(argparser_t *pp, char *key)
  { bool_t *p = pp->parsed.e;
    if (argparser_is_next(pp, key))
      { p[pp->next] = TRUE;
        pp->next++;
      }
    else
      { argparser_arg_msg(pp, "", pp->next, " should be %s.\n", key); 
        argparser_print_help_and_halt(pp, 1);
      }
  }

int64_t argparser_get_next_int(argparser_t *pp, int64_t min, int64_t max) 
  { char *txt = argparser_get_next(pp);
    int64_t nn = argparser_parse_int_string(pp, pp->next-1, txt, min, max, NULL);
    return nn;
  }

uint64_t argparser_get_next_uint(argparser_t *pp, uint64_t min, uint64_t max) 
  { uint64_t nn;
    char *rest;
    char *txt = argparser_get_next(pp);
    nn = strtoull(txt, &rest, 10);
    if ((*rest) != '\000')
      { argparser_arg_msg(pp, "", pp->next-1, " should be an unsigned integer%s.\n", "");
         argparser_print_help_and_halt(pp, 1); 
      }
    if ((nn < min) || (nn > max))
      { argparser_arg_msg(pp, "", pp->next-1, " should be in %s", "");
        fprintf(pp->wr, "[%" uint64_u_fmt "..%" uint64_u_fmt "].\n", min, max);
        argparser_print_help_and_halt(pp, 1); 
      }
    return nn;
  }

double argparser_get_next_double(argparser_t *pp, double min, double max)
  { double x;
    char *txt = argparser_get_next(pp);
    char *rest;
    x = strtod(txt, &rest);
    if ((*rest) != '\000')
      { argparser_arg_msg(pp, "", pp->next-1, " should be a number%s.\n", "");
        argparser_print_help_and_halt(pp, 1); 
      }
    if ((x < min) || (x > max))
      { argparser_arg_msg(pp, "", pp->next-1, " should be in %s", "");
        fprintf(pp->wr, "[%g _ %g].\n", min, max);
        argparser_print_help_and_halt(pp, 1); 
      }
    return x;
  }

bool_t argparser_get_next_bool(argparser_t *pp)
  {
    if (argparser_next_is_number(pp))
      { double nv = argparser_get_next_double(pp, -DBL_MAX, +DBL_MAX);
        if (nv == 0) 
          { return FALSE; }
        else if (nv == 1)
          { return TRUE; }
        else
          { argparser_arg_msg(pp, "", pp->next-1, " should be 0 or 1%s.\n", "");
            argparser_print_help_and_halt(pp, 1); 
            return FALSE; /* To pacify the compiler. */
          }
      }
    else
      { char *txt = argparser_get_next(pp);
        if 
          ( (strcmp(txt, "T") == 0) || (strcmp(txt, "t") == 0) ||
            (strcmp(txt, "TRUE") == 0) || (strcmp(txt, "true") == 0) ||
            (strcmp(txt, "Y") == 0) || (strcmp(txt, "y") == 0) ||
            (strcmp(txt, "YES") == 0) || (strcmp(txt, "yes") == 0)
          )
          { return TRUE; }
        else if 
          ( (strcmp(txt, "F") == 0) || (strcmp(txt, "f") == 0) || 
            (strcmp(txt, "FALSE") == 0) || (strcmp(txt, "false") == 0) || 
            (strcmp(txt, "N") == 0) || (strcmp(txt, "n") == 0) || 
            (strcmp(txt, "NO") == 0) || (strcmp(txt, "no") == 0)
          )
          { return FALSE; }
        else 
          { argparser_arg_msg(pp, "", pp->next-1, " should be \"TRUE\", \"FALSE\", \"YES\", \"NO\", or initial%s.\n", "");
            argparser_print_help_and_halt(pp, 1); 
            return FALSE; /* To pacify the compiler. */
          }
      }
  }
 
char *argparser_next(argparser_t *pp)
  { char **a = pp->arg.e;
    bool_t *p = pp->parsed.e;
    if ((pp->next >= pp->arg.ne) || (p[pp->next]))
      { return NULL; }
    else
      { return a[pp->next]; }
  }

bool_t argparser_is_next(argparser_t *pp, char *key)
  { char **a = pp->arg.e;
    bool_t *p = pp->parsed.e;
    if (pp->next >= pp->arg.ne) { return FALSE; }
    return (! p[pp->next]) && (argparser_key_matches(key, a[pp->next])); 
  }

bool_t argparser_next_is_keyword(argparser_t *pp)
  { char **a = pp->arg.e;
    bool_t *p = pp->parsed.e;
    if ((pp->next >= pp->arg.ne) || (p[pp->next])) { return FALSE; }
    char *n = a[pp->next];
    return argparser_seems_keyword(n);
  }

bool_t argparser_next_is_non_keyword(argparser_t *pp)
  { char **a = pp->arg.e;
    bool_t *p = pp->parsed.e;
    if ((pp->next >= pp->arg.ne) || (p[pp->next])) { return FALSE; }
    char *n = a[pp->next];
    return (! argparser_seems_keyword(n));
  }

bool_t argparser_next_is_number(argparser_t *pp)
  { char **a = pp->arg.e;
    bool_t *p = pp->parsed.e;
    if ((pp->next >= pp->arg.ne) || (p[pp->next])) { return FALSE; }
    char *n = a[pp->next];
    if ((n[0] == '+') || (n[0] == '-')) { n++; }
    return ((n[0] >= '0') && (n[0] <= '9'));
  }

bool_t argparser_keyword_present_next(argparser_t *pp, char *key)
  { bool_t *p = pp->parsed.e;
    if (argparser_is_next(pp, key))
      { p[pp->next] = TRUE;
        pp->next++;
        return TRUE;
      }
    else
      { return FALSE; }
  }
  
#define argparser_show_bogus_max 5
  /* Max leftover args to print. */

void argparser_check_all_parsed(argparser_t *pp, uint32_t num)
  { int32_t bogus = 0;
    bool_t *p = pp->parsed.e;
    for (int32_t i = 1; i < num; i++)
      { if (! p[i])
          { bogus++;
            if (bogus <= argparser_show_bogus_max)
              { argparser_arg_msg(pp, "", i, " extraneous or misplaced.%s\n", ""); }
          }
      }
    if (bogus > argparser_show_bogus_max) 
      { fprintf(pp->wr, "(and %d more).\n", bogus - argparser_show_bogus_max); }
    if (bogus > 0) { argparser_print_help_and_halt(pp, 1); }
  }

void argparser_skip_parsed(argparser_t *pp)
  { bool_t *p = pp->parsed.e;
    pp->next = pp->arg.ne;
    while ((pp->next > 0) && (! p[pp->next-1])) { pp->next--; }
    /* Check for unparsed arguments: */
    argparser_check_all_parsed(pp, pp->next);
  }

void argparser_finish(argparser_t *pp)
  { argparser_check_all_parsed(pp, pp->arg.ne);
    free(pp->parsed.e);
    free(pp->help.e);
    free(pp);
  }

/* Copyright © 2003 by Jorge Stolfi.
**
** Permission to use, copy, modify, and distribute this software and its
** documentation for any purpose and without fee is hereby granted, provided
** that the above copyright notice appears in all copies and that both that
** copyright notice and this permission notice appear in supporting
** documentation.  This software is provided "as is" without express or
** implied warranty of any kind.
*/
 
