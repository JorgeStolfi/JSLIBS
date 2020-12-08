/* See dna_seq_view_GL.h */
/* Last edited on 2014-06-09 23:10:43 by stolfilocal */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <values.h>
#include <assert.h>

#include <GL/glu.h>
#include <GL/glut.h>

#include <bool.h>

#include <dna_seq_view_GL.h>

void dna_seq_view_GL_initialize_libraries(int *argc, char** argv)
  { 
    glutInitWindowPosition(dna_seq_view_GL_default_window_corner_H, dna_seq_view_GL_default_window_corner_V);
    glutInitWindowSize(dna_seq_view_GL_default_window_HSize, dna_seq_view_GL_default_window_VSize);
    glutInit(argc, argv);
  }
  
void dna_seq_view_GL_initialize_window
  ( char *title,
    void display(void),
    void reshape(int width, int height),
    void keyboard(unsigned char key, int x, int y),
    void passivemouse(int x, int y),
    void activemouse(int x, int y),
    void special(int key, int x, int y)
  )
  { 
    /* Setup window parameters: */
    glutInitDisplayMode(GLUT_RGB|GLUT_DEPTH|GLUT_DOUBLE);

    /* Create window: */
    bool_t success = glutCreateWindow(title);
    if (!success) 
      { fprintf(stderr, "Error opening main window\n"); exit(1); }

    /* Setup graphics modes: */
    glShadeModel(GL_FLAT);
    glEnable(GL_DEPTH_TEST);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE,1);

    /* Register event handling methods: */
    glutKeyboardFunc(keyboard);
    glutMotionFunc(activemouse);
    glutPassiveMotionFunc(passivemouse);
    glutSpecialFunc(special);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
  }
  
void dna_seq_view_GL_start_viewing(void)
  { 
    /* Start displaying: */
    glutMainLoop();
  }

