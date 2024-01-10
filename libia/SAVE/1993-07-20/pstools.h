/* Generic hacks for generating postscript files */

#ifndef PSTOOLS_H
#define PSTOOLS_H

#include <stdio.h>

void pst_begin_file(FILE *psfile);

void pst_end_file(FILE *psfile);

void pst_begin_page(FILE *psfile, int page);

void pst_end_page(FILE *psfile);

void pst_set_coords(
    FILE *psfile,
    char axis,            /* 'x' or 'y' */
    double ulo, double uhi, /* Interval in client coords */
    double wlo, double whi  /* Interval in current coords */
  );

void pst_put_text(FILE *psfile, char *text, char *newline);
  /* Writes a text string to /psfile/, in Postscript form, */
  /* handling special chars and parenthese.  */
  /* Replaces any embedded '\n' by the given /newline/ string. */

#endif
