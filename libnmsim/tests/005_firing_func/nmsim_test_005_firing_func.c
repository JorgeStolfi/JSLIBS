#define PROG_NAME "nmsim_test_005_firing_func"
#define PROG_DESC "basic tests of {limnmism} firing function procedures"
#define PROG_VERS "1.0"

/* Last edited on 2022-10-20 06:34:38 by stolfi */ 

#define PROG_COPYRIGHT \
  "Copyright © 2019  State University of Campinas (UNICAMP)"
  
#define PROG_AUTH \
  "Created by J. Stolfi on 2019-01-02"
  
#define PROG_HIST
  
#define PROG_HELP \
  "  " PROG_NAME ""

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <bool.h>
#include <argparser.h>
#include <jsfile.h>
#include <jsrandom.h>
#include <affirm.h>

#include <nmsim_firing_func.h>
#include <nmsim_test.h>

int32_t main(int32_t argc, char **argv);

/* IMPLEMENTATIONS: */

int32_t main(int32_t argc, char **argv)
  { 
    char *classes = "GLN"; /* Firing function classes to plot. */
    int32_t k = 0;
    while (classes[k] != '\0')
      { nmsim_firing_func_class_t class = classes[k];
        double V_M = -40;
        double V_D = 10;
        nmsim_firing_func_t Phi = nmsim_firing_func_make(class, V_M, V_D);
        double r = 8.0;
        int32_t ns = 1000;
        nmsim_test_firing_func(&Phi, r, ns);
        k++;
      }
    return 0;
  }
    
