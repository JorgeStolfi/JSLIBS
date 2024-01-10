
#include "pstools.h"
#include "foimisc.h"
#include "iomisc.h"
#include <stdio.h>

void pst_begin_file(FILE *psfile)
  {
    fprintf(psfile, "%%! Adobe-PostScript\n");
  }

void pst_end_file(FILE *psfile)
  {
    fprintf(psfile, "%%%%Trailer\n");
    fflush(psfile);
  }

void pst_begin_page(FILE *psfile, int page)
  { 
    char *date = today();
    fprintf(stderr, "date = %s\n", date);
    fprintf(psfile, "%%%%Page %d ?\n", page);
    fprintf(psfile, "\n");
    fprintf(psfile, "%% Print date:\n");
    fprintf(psfile, "gsave\n");
    fprintf(psfile, "  /Courier findfont\n");
    fprintf(psfile, "  10 scalefont setfont\n");
    fprintf(psfile, "  80 18 moveto\n");
    fprintf(psfile, "  (%s   page %d) show\n", date, page);
    fprintf(psfile, "grestore\n");
    fprintf(psfile, "\n");
    
  }

void pst_end_page(FILE *psfile)
  {
    fprintf(psfile, "showpage\n");
    fprintf(psfile, "\n");
    fflush(psfile);
  }

void pst_set_coords(
    FILE *psfile,
    char axis,            /* 'x' or 'y' */
    double ulo, double uhi, /* Interval in client coords */
    double wlo, double whi  /* Interval in current coords */
  )
  {
    double mag = (whi - wlo)/(uhi - ulo);
    if (axis == 'x')
      {
        fprintf(psfile, "  %f 0 translate\n", wlo);
        fprintf(psfile, "  %f 1 scale\n", mag);
        fprintf(psfile, "  %f 0 translate\n", -ulo);
      }
    else if (axis == 'y')
      {
        fprintf(psfile, "  0 %f translate\n", wlo);
        fprintf(psfile, "  1 %f scale\n", mag);
        fprintf(psfile, "  0 %f translate\n", -ulo);
      }
    else
      { error ("pst_set_coords: invalid axis"); }

    fflush(psfile);
  }


void pst_put_text(FILE *psfile, char *text, char *newline)
  {
    char *p;
    putc('(', psfile);
    for (p = text; *p != 0; p++)
      {
        if (*p == '\n')
          { fprintf(psfile, "%s", newline); }
        else if (*p == '(')
          { putc('\\', psfile); putc('(', psfile); }
        else if (*p == ')')
          { putc('\\', psfile); putc(')', psfile); }
        else if (*p == '\t')
          { putc(' ', psfile); putc(' ', psfile); }
        else if (*p == '\\')
          { putc('\\', psfile); putc('\\', psfile); }
        else if ((*p < ' ') || (*p > '~'))
          { fprintf(psfile, "\\%03o", *p); }
        else
          { putc(*p, psfile); }
      }
    fprintf(psfile, ")");
  }






