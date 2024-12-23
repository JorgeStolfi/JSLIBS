/* Basic tests univariate minimizers. */
/* Last edited on 2024-12-21 11:24:43 by stolfi  */

#include <stdio.h>
#include <values.h>
#include <stdint.h>
#include <math.h>
#include <stdlib.h>

#include <epswr.h>
#include <jsrandom.h>

#include <minu_herm.h>
#include <minu_brent.h>
#include <minu_js.h>
#include <minu_gen.h>

#include <test_minu_problems.h>
#include <test_minu_tools.h>

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
  for (uint32_t i_opt = 0;  i_opt < NOpt; i_opt++)
    { for (uint32_t i_prb = 0;  i_prb < NPrb; i_prb++) 
        { uint32_t ij = NPrb*i_opt+i_prb;
          test_minu_tools_single(i_opt, &(opt[i_opt]), i_prb, &(prb[i_prb]), TRUE);
          srandom(2567898753u);
          perf[ij] = test_minu_tools_multiple(i_opt, &(opt[i_opt]), i_prb, &(prb[i_prb]), 100);
        }
    }

  /* Print summary: */
  fprintf(stderr, "\n");
  fprintf(stderr, "=============================================================================\n");
  fprintf(stderr, "* SUMMARY - avgCalls/avgError/nFailures\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "%30s ", "");
  for (uint32_t i_prb = 0;  i_prb < NPrb; i_prb++)
    { fprintf(stderr, "%8s ", prb[i_prb].name); }
  fprintf(stderr, "\n");
  fprintf(stderr, "\n");
        
  for (uint32_t i_opt = 0;  i_opt < NOpt; i_opt++)
    { fprintf(stderr, "%30s ", opt[i_opt].name);
      for (uint32_t i_prb = 0;  i_prb < NPrb; i_prb++)
        { uint32_t ij = NPrb*i_opt+i_prb;
          fprintf(stderr, "%8.1f ", perf[ij].avgCalls);
        }
      fprintf(stderr, "\n");
      fprintf(stderr, "%30s ", "");
      for (uint32_t i_prb = 0;  i_prb < NPrb; i_prb++)
        { uint32_t ij = NPrb*i_opt+i_prb;
          fprintf(stderr,  "%8.1f ", perf[ij].avgError);
        }
      fprintf(stderr, "\n");
      fprintf(stderr, "%30s ", "");
      for (uint32_t i_prb = 0;  i_prb < NPrb; i_prb++)
        { uint32_t ij = NPrb*i_opt+i_prb;
          fprintf(stderr, "%8d ", perf[ij].nFailures);
        }
      fprintf(stderr, "\n");
      fprintf(stderr, "\n");
    }
  return 0;
}

