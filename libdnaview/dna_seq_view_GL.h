#ifndef dna_seq_view_GL_H
#define dna_seq_view_GL_H

/* dna_seq_view_GL.h - model-independent graphics routines for dna_seq_view(1). */
/* Last edited on 2014-06-09 23:10:24 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <values.h>
#include <assert.h>

#include <bool.h>

#define dna_seq_view_GL_default_window_corner_H  64
#define dna_seq_view_GL_default_window_corner_V  64
  /* Default screen coordinates of window's top left corner. */

#define dna_seq_view_GL_default_window_HSize  960
#define dna_seq_view_GL_default_window_VSize  960
  /* Default window size. */

void dna_seq_view_GL_initialize_libraries(int *argc, char** argv);
  /* Intializes the GL library.  Also parses any X11/glut command line 
    arguments in {argv}, and deletes them.  Therefore, it should be 
    called BEFORE the program's own option parsing. */

void dna_seq_view_GL_initialize_window
  ( char *title,
    void display(void),
    void reshape(int width, int height),
    void keyboard(unsigned char key, int x, int y),
    void passivemouse(int x, int y),
    void activemouse(int x, int y),
    void special(int key, int x, int y)
  );
  /* Intializes the GL window attributes with the given title
    and methods.  Assumes that {glutInit} has been called. */

void dna_seq_view_GL_start_viewing(void);
  /* Sets some graphics state parameters and starts the GL main loop. */


#endif
