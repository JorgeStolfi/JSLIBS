#ifndef dna_seq_view_GL_H
#define dna_seq_view_GL_H

/* dna_seq_view_GL.h - model-independent graphics routines for dna_seq_view(1). */
/* Last edited on 2022-10-20 06:32:27 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
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

void dna_seq_view_GL_initialize_libraries(int32_t *argc, char** argv);
  /* Intializes the GL library.  Also parses any X11/glut command line 
    arguments in {argv}, and deletes them.  Therefore, it should be 
    called BEFORE the program's own option parsing. */

void dna_seq_view_GL_initialize_window
  ( char *title,
    void display(void),
    void reshape(int32_t width, int32_t height),
    void keyboard(unsigned char key, int32_t x, int32_t y),
    void passivemouse(int32_t x, int32_t y),
    void activemouse(int32_t x, int32_t y),
    void special(int32_t key, int32_t x, int32_t y)
  );
  /* Intializes the GL window attributes with the given title
    and methods.  Assumes that {glutInit} has been called. */

void dna_seq_view_GL_start_viewing(void);
  /* Sets some graphics state parameters and starts the GL main loop. */


#endif
