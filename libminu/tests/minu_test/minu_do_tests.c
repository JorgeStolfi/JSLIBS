// Tests univariate minimizers.


#include <minu_gen.h>
#include <stdint.h>
#include <minu_js.h>
#include <minu_brent.h>
#include <minu_herm.h>
#include <minu_test.h>
#include <minu_test_problems.h>
#include <stdlib.h>
#include <math.h>
#include <pswr.h>
#include <values.h>
#include <stdio.h>

#define MaxFCalls (50)
#define Phi (1.61803398874989484821)   /* Golden ratio */

int32_t main(int32_t argc, char **argv)
{
  #define NPrb 6
  #define NOpt 3
  Problem prb[NPrb] =
    {
      (Problem){ 
        /*name*/     "quad",
        /*eval*/     QuadEval,
        /*error*/    QuadError,
        /*xRange*/   QuadXMin, QuadXMax, 
        /*yRange*/   QuadYMin, QuadYMax,
        /*start*/    QuadXStart,
        /*tol*/      1.0e-4,
        /*dist*/     QuadDist,
        /*maxCalls*/ MaxFCalls
      },
      (Problem){
        /*name*/ "bent",
        /*eval*/     BentEval,
        /*error*/    BentError,
        /*xRange*/   BentXMin, BentXMax,
        /*yRange*/   BentYMin, BentYMax,
        /*xStart*/   BentXStart,
        /*tol*/      1.0e-4,
        /*dist*/     BentDist,
        /*maxCalls*/ MaxFCalls
      },
      (Problem){
        /*name*/     "biqu",
        /*eval*/     BiquEval,
        /*error*/    BiquError,
        /*xRange*/   BiquXMin, BiquXMax,
        /*yRange*/   BiquYMin, BiquYMax,
        /*xStart*/   BiquXStart,
        /*tol*/      1.0e-4,
        /*dist*/     BiquDist,
        /*maxCalls*/ MaxFCalls
      },
      (Problem){
        /*name*/     "hole",
        /*eval*/     HoleEval,
        /*error*/    HoleError,
        /*xRange*/   HoleXMin, HoleXMax,
        /*yRange*/   HoleYMin, HoleYMax,
        /*xStart*/   HoleXStart,
        /*tol*/      1.0e-4,
        /*dist*/     HoleDist,
        /*maxCalls*/ MaxFCalls
      },
      (Problem){
        /*name*/    "wavy",
        /*eval*/     WavyEval,
        /*error*/    WavyError,
        /*xRange*/   WavyXMin, WavyXMax,
        /*yRange*/   WavyYMin, WavyYMax,
        /*xStart*/   WavyXStart,
        /*tol*/      1.0e-4,
        /*dist*/     WavyDist,
        /*maxCalls*/ MaxFCalls
      },
      (Problem){
        /*name*/     "corn",
        /*eval*/     CornEval,
        /*error*/    CornError,
        /*xRange*/   CornXMin, CornXMax,
        /*yRange*/   CornYMin, CornYMax,
        /*xStart*/   CornXStart,
        /*tol*/      1.0e-4,
        /*dist*/     CornDist,
        /*maxCalls*/ MaxFCalls
      }
    };
  Minimizer opt [NOpt] =
    { 
      (Minimizer){
        /*name*/     "brent",
        /*minimize*/  minu_brent_minimize,
        /*deriv*/     FALSE
      },
      (Minimizer){
        /*name*/     "js",
        /*minimize*/  minu_js_minimize,
        /*deriv*/     FALSE
      },
      (Minimizer){
        /*name*/     "herm",
        /*minimize*/  minu_herm_minimize,
        /*deriv*/     TRUE
      }
    };
  Performance perf[NOpt*NPrb];
  char psfname[80];
  PSStream *ps;
  int32_t i,j;
  
  for (i = 0; i < NOpt; i++)
    { sprintf(psfname, "aumt-%02d-sngl", i);
      ps = pswr_new_stream(psfname, NULL, FALSE, "doc", "letter", FALSE, 0,0);
      for (j = 0; j < NPrb; j++) 
        { int32_t ij = NPrb*i+j;
          minu_test_single(ps, j, &(opt[i]), &(prb[j]), TRUE);
          srandom(2567898753u);
          perf[ij] = minu_test_multiple(&(opt[i]), &(prb[j]), 100);
        }
      pswr_close_stream(ps);
    }

  /* Print summary: */
  fprintf(stderr, "\n");
  fprintf(stderr, "=============================================================================\n");
  fprintf(stderr, "* SUMMARY\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "%30s ", "");
  for (j = 0; j < NPrb; j++)
    { fprintf(stderr, "%8s ", prb[j].name); }
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
        
  for (i = 0; i < NOpt; i++)
    { fprintf(stderr, "%30s ", opt[i].name);
      for (j = 0; j < NPrb; j++)
        { int32_t ij = NPrb*i+j;
          fprintf(stderr, "%8.1f ", perf[ij].avgCalls);
        }
      fprintf(stderr, "\n");
      fprintf(stderr, "%30s ", "");
      for (j = 0; j < NPrb; j++)
        { int32_t ij = NPrb*i+j;
          fprintf(stderr,  "%8.1f ", perf[ij].avgError);
        }
      fprintf(stderr, "\n");
      fprintf(stderr, "%30s ", "");
      for (j = 0; j < NPrb; j++)
        { int32_t ij = NPrb*i+j;
          fprintf(stderr, "%8d ", perf[ij].nFailures);
        }
      fprintf(stderr, "\n");
      fprintf(stderr, "\n");
    }
  return 0;
}

