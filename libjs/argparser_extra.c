/* See argparser_extra.h. */
/* Last edited on 2023-10-10 12:27:42 by stolfi */

/* Copyright © 2003 Jorge Stolfi, Unicamp. See note at end of file. */
/* Based on Params.m3 by J.Stolfi, DEC-SRC, 1988.  */

#define _GNU_SOURCE
#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <errno.h>
#include <assert.h>

#include <affirm.h>
#include <jswsize.h>
#include <argparser.h>

#include <argparser_extra.h>

char *argparser_print_info_line(FILE *wr, char *lin, int32_t wd);
  /* Prints the string {lin} (terminated by '\000' or '\n') to {wr}.
     Tries to breaks it so that it does not exceed {wd} characters per
     line, replicating its initial indentation on each subsequent
     piece. Such breaks are introduced at internal blanks. Any blanks
     at the break points or at the end of {lin} are discarded. Returns
     a pointer to the line's terminating character, which is is not
     printed. */
  
bool_t argparser_seems_keyword(char *n)
  { if (n[0] != '-') { return FALSE; }
    if (n[1] == '-') { return TRUE; }
    return (((n[1] >= 'a') && (n[1] <= 'z')) || ((n[1] >= 'A') && (n[1] <= 'Z')));
  }
  
bool_t argparser_key_matches(char *a, char *b)
  { if (strcmp(a, b) == 0) 
      { return TRUE; }
    else if ((a[0] == '-') && (b[0] == '-'))
      { if ((a[1] == '-') && (strcmp(a+1,b) == 0)) { return TRUE; }
        if ((b[1] == '-') && (strcmp(a,b+1) == 0)) { return TRUE; }
      }
    return FALSE;
  }
    
int64_t argparser_parse_int_string(argparser_t *pp, int32_t index, char *arg, int64_t min, int64_t max, char **rest_P)
  { char *rest = NULL;
    errno = 0;
    int64_t v = (int64_t)strtoll(arg, &rest, 10);
    if (errno == EINVAL)
      { if (v == LONG_MAX)
          { argparser_error_at(pp, "integer overflow", "at", index); }
        else if (v == LONG_MIN)
          { argparser_error_at(pp, "integer underflow", "at", index); }
        else
          { assert(v == 0);
            argparser_error_at(pp, "invalid integer", "at", index);
          }
      }
    else if (rest == arg)
      { /* No digits parsed. Must be syntax error. */
        argparser_error_at(pp, "invalid integer", "at", index);
      }
    else if (((*rest) != '\000') && (rest_P == NULL))
      { /* No digits parsed. Must be syntax error. */
        argparser_error_at(pp, "invalid integer", "at", index);
      }
    else
      { if ((v > max) || (v < min)) 
         { argparser_arg_msg(pp, "", index, " should be in %s", "");
           fprintf(pp->wr, ("[%" int64_d_fmt "..%" int64_d_fmt "].\n"), min, max);
           argparser_print_help_and_halt(pp, 1);
         }
      }
    if (rest_P != NULL) { (*rest_P) = rest; }
    return v;
  }
 
void argparser_arg_msg(argparser_t *pp, char *msg1, int32_t index, char *msg2, char *val)
  { fprintf(pp->wr, "%s: %s", pp->arg.e[0], msg1); 
    if ((index > 0) && (index < pp->arg.ne))
      { fprintf(pp->wr, "parameter %d = \"%s\"", index, pp->arg.e[index]); }
    fprintf(pp->wr, msg2, val);
  }

void argparser_error_at(argparser_t *pp, char *msg, char* pos, int32_t index)
  { fprintf(pp->wr, "%s: ** %s\n", pp->arg.e[0], msg);
    if ((index > 0) && (index < pp->arg.ne))
      { fprintf(pp->wr, " %s argument %d = \"%s\"", pos, index, pp->arg.e[index]);
        fprintf(pp->wr, " (%s)", (pp->parsed.e[index] ? "parsed" : "unparsed"));
        fprintf(pp->wr, "\n");
      }
    argparser_print_help_and_halt(pp, 1);
  }

void argparser_print_info(argparser_t *pp, int32_t wd)
  { int32_t i;
    for (i = 0; i < pp->nhelp; i++)
      { argparser_print_text(pp->wr, pp->info.e[i], wd); }
    fprintf(pp->wr, "\n");
  }

void argparser_print_text(FILE *wr, char *info, int32_t wd)
  { char *lin = info;
    while ((*lin) != '\000')
      { /* Print line starting at {lin}, advance {lin} to its end: */
        lin = argparser_print_info_line(wr, lin, wd);
        /* Skip eol, if any: */
        if ((*lin) == '\n') { lin++; }
      }
  }

char *argparser_print_info_line(FILE *wr, char *lin, int32_t wd)
  { /* Compute the indentation of this input line: */
    int32_t indent = 0;
    while ((*lin) == ' ') { lin++; indent++; }
    do 
      { /* Get next output line, from {lin} to {end}: */
        char *end = lin;
        char *q = end;
        while(((*q) != '\000') && ((*q) != '\n'))
          { /* Advance {q} to end of current word: */
            while (((*q) != ' ') && ((*q) != '\000') && ((*q) != '\n')) { q++; }
            /* Can we add that word to the output line? */
            if ((end > lin) && ((q - lin) + indent > wd)) { break; }
            /* Yes, add it: */
            end = q;
            /* Advance {q} to beginning of next word: */
            while ((*q) == ' ') { q++; }
          }
        /* Print the current output line with fixed indent: */
        fprintf(wr, "%*s%.*s\n", indent, "", (int32_t)(end-lin), lin);
        fflush(wr);
        lin = end;
        /* Skip blanks: */
        while ((*lin) == ' ') { lin++; }
      }
    while (((*lin) != '\000') && ((*lin) != '\n'));
    return lin;
  }
