/* See cpk_main.h */
/* Last edited on 2024-12-31 16:15:54 by stolfi */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <bool.h>
#include <interval.h>
#include <r2.h>
#include <affirm.h>

#include <cpk_main.h>
#include <cpk_basic.h>
#include <cpk_build.h>
#include <cpk_greedy.h>
#include <cpk_grasp.h>
#include <cpk_graph.h>
#include <cpk_io.h>
#include <cpk_debug.h>
#include <cpk_valid.h>

/* INTERNAL PROTOS */

double cpk_compute_subset_value(uint32_vec_t *S, double_vec_t *W);
  /* Returns the total weight of the vertex set {S},
    that is, the sum of {W[S[i]]} for all {i} in {0..S.ne-1}. */

/* IMPLS */

void cpk_choose_auction_points 
  ( cpk_domain_t *C,   
    cpk_policy_t *P,     
    cpk_options_t *opt,  
    r2_vec_t *VJ,        
    double *WJ           
  )
  { 
    bool_t verbose = opt->verbose;
    if (verbose) { cpk_trace_entry(); }
    
    if (verbose) 
      { fprintf(stderr, "Policy parameters:\n");
        fprintf(stderr, "  min station-station distance = " XY_FFMT "\n", P->dMin);
        fprintf(stderr, "  min auction ctr dist to other municipalities = " XY_FFMT "\n", P->dMun);
        fprintf(stderr, "  min auction ctr dist to other nations = " XY_FFMT "\n", P->dNat);
        fprintf(stderr, "  radius of auction area = " XY_FFMT "\n", P->rAuc);
      }
    
    if (verbose)
      { /* Find Lon/Lat bounding box, for documentation onl: */
        interval_t LLB[2];
        r2_bbox(C->Urb.ne, C->Urb.e, LLB, TRUE);
        fprintf
          ( stderr, 
            "LL bounding box = [" LL_FFMT " _ " LL_FFMT "] × [" LL_FFMT " _ " LL_FFMT "]\n",
            LO(LLB[0]), HI(LLB[0]),   LO(LLB[1]), HI(LLB[1])
          );
      }        

    /* Compute mean longitude and mean Y for EUTM conversion: */
    double refLon, refY;
    cpk_pick_ref_coords(&(C->Urb), &refLon, &refY);
    if (verbose) { fprintf(stderr, "refLon = " LL_FMT "  refY = " XY_FMT "\n", refLon, refY); }

    /* Domain in EUTM XY coordinates: */
    cpk_domain_t *CXY;
    CXY = (cpk_domain_t *)notnull(malloc(sizeof(cpk_domain_t)), "no mem");

    /* Convert all geographic data from LL (degrees) to EUTM (XY): */
    CXY->Urb = cpk_LL_to_EUTM(&(C->Urb), refLon, refY, opt->magnify);
    CXY->Exs = cpk_LL_to_EUTM(&(C->Exs), refLon, refY, opt->magnify);
    CXY->Auc = cpk_LL_to_EUTM(&(C->Auc), refLon, refY, opt->magnify);
    CXY->Mun = cpk_LL_to_EUTM(&(C->Mun), refLon, refY, opt->magnify);
    CXY->Nat = cpk_LL_to_EUTM(&(C->Nat), refLon, refY, opt->magnify);
    CXY->Dem = cpk_LL_to_EUTM(&(C->Dem), refLon, refY, opt->magnify);

    /* Find bounding box and general extent of polygon: */
    interval_t B[2];
    r2_bbox(CXY->Urb.ne, CXY->Urb.e, B, TRUE);
    if (verbose)
      { fprintf
          ( stderr, 
            "XY bounding box = [" XY_FFMT " _ " XY_FFMT "] × [" XY_FFMT " _ " XY_FFMT "]\n",
            LO(B[0]), HI(B[0]),   LO(B[1]), HI(B[1])
          );
      }        
    /* double zx = HI(B[0]) - LO(B[0]), zy = HI(B[1]) - LO(B[1]); */
    /* double apmm = hypot(zx,zy)/350; */ /* Approx 1 mm in the plot */

    /* Get the candidates and their incompatibility edges: */
    r2_vec_t V;
    double_vec_t W;
    ui2_vec_t E;
    cpk_build_graph(CXY, P, &V, &W, &E, verbose);
    uint32_t nV = V.ne;

    /* Print the candidates: */
    if (verbose)
      { cpk_print_vertices(stderr, "candidates for auctions", &V, &W, 25);
        cpk_print_edges(stderr, "incompatible pairs", &E, &V, 25);
      }

    if (opt->validate)
      { assert(cpk_check_edges(&V, &E, P->dMin + 2*P->rAuc)); }
      
    if (nV == 0)
      { if (verbose) { fprintf(stderr, "Empty graph!\n"); }
        (*VJ) = r2_vec_new(0);
        (*WJ) = 0.0;
      }
    else
      { double start_time;
      
      /* Convert incompatibility graph to neighbor list format: */
        cpk_graph_t G = cpk_gather_neighbors(nV, E);    

        /* Demand suffix {demtag} for file name"*/
        char *demTag = (P->for_demand > 0 ? "wd" : "nd");

        /* Compute a quick greedy solution: */
        char *alg0 = "greed";
        start_time = cpk_cpu_time_1();
        uint32_vec_t J0 = cpk_greedy_find_indep_set(&G, W.e, verbose);
        double greedy_time = (cpk_cpu_time_1() - start_time)/1000000;
        if (verbose) { fprintf(stderr, "greedy time = %.2f sec\n", greedy_time); }
        double WJ0 = cpk_compute_subset_value(&J0, &W);

        if (opt->validate)
          { assert(cpk_check_solution(CXY, P, &V, &J0)); }

        if (opt->plot) 
          { cpk_plot_solution(opt->outDir, alg0, demTag, B, CXY, P, &V, &W, &J0); }

        /* Compute the GRASP solution: */
        char *alg1 = "grasp";
        start_time = cpk_cpu_time_1();
        double maxClock = start_time + opt->maxSecsGRASP*1000000;
        uint32_vec_t J1 = cpk_grasp_find_indep_set(&G, W.e, maxClock, opt->seed, verbose);
        double grasp_time = (cpk_cpu_time_1() - start_time)/1000000;
        if (verbose) { fprintf(stderr, "grasp time = %.2f sec\n", grasp_time); }
        double WJ1 = cpk_compute_subset_value(&J1, &W);

        if (opt->validate)
          { assert(cpk_check_solution(CXY, P, &V, &J1)); }

        if (opt->plot) 
          { cpk_plot_solution(opt->outDir, alg1, demTag, B, CXY, P, &V, &W, &J1); }

        /* Ensure that the best solution is in {J1}: */
        if (WJ0 > WJ1)
          { if (verbose) { fprintf(stderr, "Swapping solutions...\n"); }
            { uint32_vec_t JT = J0; J0 = J1; J1 = JT; }
            { double WJT = WJ0; WJ0 = WJ1; WJ1 = WJT; }
            { char *algT = alg0; alg0 = alg1; alg1 = algT; }
          }
        if (verbose) 
          { fprintf(stderr, "Algorithm \"%s\" produced the best solution\n", alg1); }
        
        /* Extract proposed auction centers: */
        r2_vec_t VJXY = r2_vec_new(J1.ne);
        for (int32_t i = 0; i < J1.ne; i++) { VJXY.e[i] = V.e[J1.e[i]]; }
        
        /* Map proposed auction centers back to Lon/Lat pairs: */
        (*VJ) = cpk_EUTM_to_LL(&VJXY, refLon, refY, opt->magnify);
        
        /* Return score: */
        (*WJ) = WJ1;
        
        /* Print the result: */
        if (verbose) { cpk_print_vertices(stderr, "proposed auctions", VJ, NULL, 25); }

        /* Housecleaning: */
        free(J0.e); free(J1.e);
        free(G.deg); free(G.fnb); free(G.nbr); 
        free(VJXY.e); 
      }
    
    /* Housecleaning: */
    free(V.e); free(W.e); free(E.e);
    cpk_domain_free(CXY); free(CXY);
    if (verbose) { cpk_trace_exit(); }
  }

void cpk_domain_free(cpk_domain_t *C)
  { free(C->Urb.e);
    free(C->Exs.e);
    free(C->Auc.e);
    free(C->Mun.e);
    free(C->Nat.e);
    free(C->Dem.e);
  }

void cpk_policy_free(cpk_policy_t *P)
  {  }

double cpk_compute_subset_value(uint32_vec_t *S, double_vec_t *W)
  { double s = 0;
    for (int32_t i = 0; i < S->ne; i++) { s += W->e[S->e[i]]; }
    return s;
  }
