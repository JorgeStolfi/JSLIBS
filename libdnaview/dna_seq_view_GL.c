/* See dna_seq_view_GL.h */
/* Last edited on 2022-10-20 10:15:03 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <values.h>
#include <assert.h>

#define GLFW_INCLUDE_NONE
#include <GL/freeglut.h>
#include <GL/gl.h>
#include <GLFW/glfw3.h>

#include <bool.h>

#include <dna_seq_view_GL.h>

void dna_seq_view_GL_initialize_libraries(int32_t *argc, char** argv)
  { 
    glutInitWindowPosition(dna_seq_view_GL_default_window_corner_H, dna_seq_view_GL_default_window_corner_V);
    glutInitWindowSize(dna_seq_view_GL_default_window_HSize, dna_seq_view_GL_default_window_VSize);
    glutInit(argc, argv);
  }
  
void dna_seq_view_GL_initialize_window
  ( char *title,
    void display(void),
    void reshape(int32_t width, int32_t height),
    void keyboard(unsigned char key, int32_t x, int32_t y),
    void passivemouse(int32_t x, int32_t y),
    void activemouse(int32_t x, int32_t y),
    void special(int32_t key, int32_t x, int32_t y)
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

