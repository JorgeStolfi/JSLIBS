/*  Last edited on 2024-10-06 17:00:05 by stolfi */
/* Test of {tosl_mesh_slice} in {totposlic.h} */

#define _GNU_SOURCE
#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <tosl.h>
#include <tosl_slice.h>
#include <tosl_mesh.h>
#include <tosl_mesh_slice.h>
#include <tosl_mesh_make_keg.h>
#include <tosl_build_lists_hash.h>

int32_t main(int32_t argc, char **argv);

void do_test(int32_t NS, int32_t NR, int32_t NB, int32_t NP, int32_t NT);
  /* Tests the procedure {tosl_mesh_slice} {NT} times, all with the 
    mesh created by {tosl_mesh_make_keg(NS,NR,NB,...)} but different sets
    sets of {NP} slicing planes.  Prints the average and deviation of
    the user CPU time spent slicing, excluding the list building time,
    and the slice post-processing time. 
    
    The post-processing of each slice simply does {free} on it. */

tosl_coord_t *pick_planes(int32_t NP, tosl_coord_t Zmax, int8_t debug);
  /* Returns a vector of {NP} {Z}-coordinates of planes, all odd, 
    strictly increasing, equally spaced, strictly inside the range
    {-Zmax .. +Zmax}. */

double check_time
  ( int32_t NE,
    tosl_arc_t Arc[],
    int32_t NV,
    tosl_point_t Vpos[],
    int32_t NP,
    tosl_coord_t Zplane[],
    tosl_arc_id_t L[] 
  );
  /* Returns the time of {tosl_mesh_slice}, in seconds,
    for the mesh with arcs {Arc[0..2*NE-1]} and vertices {Vpos[0..NV-1]}
    and the slicing planes with Z-coordinates {Zplane[0..NP-1]}.
    Assumes that the 
    
    If either or both of the boolean arguments {use_bin} and {use_sec} are 1 (true),
    uses the binary and/or secant search method with those parameters.  If both
    {use_bin} and {use_sec} are 0 (false), uses the hash method. */

int32_t main(int32_t argc, char **argv)
  {
    srandom(4615);
    
    fprintf(stderr, "%6s %6s %6s  %10s %10s %12s %8s %8s\n", "NS", "NR", "NB", "NE", "NV", "NM", "avg", "dev");
    do_test(10, 20, 2, 1000, 10);
    
    return 0;
  }
  
void do_test(int32_t NS, int32_t NR, int32_t NB, int32_t NP, int32_t NT)
  {
    int8_t debug = 1;
    
    /* Decide the size of the barrel: */
    tosl_coord_t Zmax_plane = 100*NP; /* So that the plane spacing is about 100. */
    tosl_coord_t Zmax_edge = 100*NR*(NB+1); /* So that the vertex spacing is about 100. */
    tosl_coord_t Zmax = (Zmax_plane > Zmax_edge ? Zmax_plane : Zmax_edge);
    
    fprintf(stderr, "creating test mesh...\n");
    tosl_mesh_t *mesh = tosl_mesh_make_keg(NS, NR, NB, Zmax);
    
    fprintf(stderr, "checking test mesh...\n");
    tosl_mesh_check(mesh);
    
    /* Estimate number of plane-edge intersections: */
    int32_t NM_exp = NP*NS;

    /* Timing tests: */
    double time[NT]; /* CPU times of runs. */
    for (int32_t it = 0; it < NT; it++)
      { fprintf(stderr, "choosing slicing planes...\n");
        tosl_coord_t *Zplane = pick_planes(NP, Zmax, debug);
        
        fprintf(stderr, "creating bucket lists...\n");
        tosl_arc_id_t *L = tosl_build_lists_hash(NP, Zplane, mesh, 0);
    
        /* Run the slicer: */
        int32_t NM_cmp = 0;  /* Total number of intersections actualy found. */
        
        auto void post(tosl_plane_id_t ip, tosl_slice_t *S);
          /* Consumes a slice.  Mainly, runs {free} on it. */
        
        fprintf(stderr, "slicing mesh...\n");
        tosl_mesh_slice(mesh, NP, Zplane, L, &post, &(time[it]));
        
        assert(NM_cmp == NM_exp);
        free(Zplane);
        free(L);
        
        void post(tosl_plane_id_t ip, tosl_slice_t *S)
          { assert(S->Z == Zplane[ip]);
            NM_cmp += S->NV;
            tosl_slice_free(S);
          }
      }
    /* Compute average and deviation of times: */
    double avg_t, dev_t;
    tosl_compute_avg_dev(NT, time, &avg_t, &dev_t);
    fprintf(stderr, "%6d %6d %6d  %10d %10d %12d %8.3f %8.3f", NS, NR, NB, mesh->NE, mesh->NV, NM_exp, avg_t, dev_t);
  }
