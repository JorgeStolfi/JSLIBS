/* Access to the PSWriter associated to a Postscript file */
/* Last edited on 2007-06-02 16:30:16 by stolfi */

#ifndef ps_compat_H
#define ps_compat_H

/* 
  This module is an extension of {ps.h}. It gives access to the
  {PSStream} associated to a Postscript file. Its purpose is to allow
  a more gradual upgrade of programs that used the old "ps.h"
  interface to the new "pswr.h" interface. */

#include <ps.h>
#include <pswr.h>
#include <stdio.h>

PSStream *ps_get_stream(FILE *psFile);
  /* Gets the current PS stream for file {psFile}. The file must have
    been initialized with {ps_begin_document} or ps_begin_figure. */

#endif
