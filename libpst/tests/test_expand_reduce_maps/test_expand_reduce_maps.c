#define PROG_NAME "test_expand_reduce_maps"
#define PROG_DESC "checks the {sample_conv.h} routines"
#define PROG_VERS "1.0"

/* Last edited on 2024-12-23 07:13:00 by stolfi */
/* Created on 2007-07-11 by J. Stolfi, UNICAMP */

#define test_expand_reduce_maps_C_COPYRIGHT \
  "Copyright © 2010  by the State University of Campinas (UNICAMP)"

#define PROG_HELP \
  "  " PROG_NAME " \\\n" \
  "    -function {ZFUNC} \\\n" \
  "    " argparser_help_info_HELP ""

#define PROG_INFO \
  "NAME\n" \
  "  " PROG_NAME " - " PROG_DESC "\n" \
  "\n" \
  "SYNOPSIS\n" \
  PROG_HELP "\n" \
  "\n" \
  "DESCRIPTION\n" \
  PROG_INFO_DESC "\n" \
  "OPTIONS\n" \
  PROG_INFO_OPTS "\n" \
  "DOCUMENTATION OPTIONS\n" \
  argparser_help_info_HELP_INFO "\n" \
  "\n" \
  "SEE ALSO\n" \
  "  test_encode_gamma(1).\n" \
  "\n" \
  "AUTHOR\n" \
  "  Created 2010-08-14 by J. Stolfi, IC-UNICAMP.\n" \
  "\n" \
  "WARRANTY\n" \
  argparser_help_info_NO_WARRANTY "\n" \
  "\n" \
  "RIGHTS\n" \
  "  " test_expand_reduce_maps_C_COPYRIGHT ".\n" \
  "\n" \
  argparser_help_info_STANDARD_RIGHTS

#define PROG_INFO_DESC \
  "  The program checks the functions {pst_height_map_expand}, {pst_height_map_reduce}," \
  " {pst_slope_map_expand}, {pst_slope_map_reduce}."

#define PROG_INFO_OPTS \
  "  -function {ZFUNC}\n" \
  "    Specifies the number of the function to use. See {pst_proc_map.h}"

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <float_image.h>
#include <pst_slope_map.h>
#include <pst_height_map.h>
#include <pst_proc_map.h>
#include <bool.h>
#include <jsmath.h>
#include <jsrandom.h>
#include <jsfile.h>
#include <jsprintf.h>
#include <affirm.h>
#include <vec.h>
#include <argparser.h>

typedef struct options_t
  { int32_t function;  /* Number of function to use. */
  } options_t;
  
typedef pst_proc_map_zfunc_t zfunc_t;

options_t *tb_parse_options(int32_t argc, char **argv);
  /* Parses the command line arguments and returns them as
    an {options_t} record. */

int32_t main(int32_t argc, char **argv);

void test_r_e_slope_map(int32_t nfunc);
  /* Tests {pst_slope_map_expand}, {pst_slope_map_reduce}.  The integer {nfunc} is 
    the iteration number. Output files will be called "out/{nfunc}-{NX}-{NY}-{tag}.fni"
    where {nfunc} is a procedural image ID number (see {pst_proc_map.h}),
    {NX,NY} are the dimensions of the test map, and {tag} is "IG", "IW", "JG", etc. */

void make_maps(int32_t NX, int32_t NY, zfunc_t func, double maxGDiff, float_image_t **ZP, float_image_t **GP, float_image_t **WP);
  /* Creates a height map {*ZP} of size {(NX+1)×(NY+1)} from a procedural image {func}. 
    Also creates its gradient map {*GP} and the corresponding weight map {*WP},
    both with size {NX×NY}. The {maxGDiff} is the threshold that defines
    the weight map (see {pst_proc_map_make_images}).  Allocates the images internally. If any of
    {ZP,GP,WP} is null, omits the corresponding map. */

void compare_and_write_slope_maps
  ( int32_t nfunc,
    float_image_t *IZ,
    float_image_t *IG, 
    float_image_t *IW, 
    float_image_t *JG, 
    float_image_t *JW, 
    float_image_t *KZ, 
    float_image_t *KG,
    float_image_t *KW
  );
  /* Takes a gradient map {IG} and its weight {IW}, a reduced version {JG,JW}
     of the same, and a reference image {KG,KW} for {JG,JW}. Writes all
     images and also the difference maps {DG=JG-KG} and {DW=JW-KW}.  The parameter {nfunc} is the 
     function ID used for file names. All images are written to files 
     "out/{nfunc}-{NX}-{NY}-{tag}.fni" where {NX,NY} are the dimensions of {IG},
     and {tag} identifies the map. */

void write_map(float_image_t *M, int32_t nfunc, int32_t NX,int32_t NY, char *tag);
  /* Writes the image {M} as a FNI file called "out/{nfunc}-{NX}-{NY}-{tag}.fni", where
    {nfunc} is formatted as 4 digits zero-filled.  The parameters {NX} and {NY} 
    need not be the dimensions of {M}. */

void test_r_e_height_map(int32_t nfunc);
  /* Tests {pst_height_map_expand}, {pst_height_map_reduce}.  The integer {nfunc} is 
    the iteration number. Output files will be called "out/{nfunc}-{tag}.fni"
    where {tag} is "IZ", "IG", etc. */

int32_t main(int32_t argc, char **argv)
  {
   /* Parse the command line options: */
    options_t *o = tb_parse_options(argc, argv);
    
    /* Create a test map: */
    float_image_t *IZ = NULL;
    float_image_t *IG = NULL;
    float_image_t *IW = NULL;
    test_r_e_slope_map(o->function);
    test_r_e_height_map(o->function);
    fprintf(stderr, "\n");
    float_image_free(IZ);
    float_image_free(IG);
    float_image_free(IW);

    /* Cleanup: */
    free(o); o = NULL;
    return 0;
  }

void test_r_e_slope_map(int32_t nfunc)
  {
    fprintf(stderr, "--- testing slopes ---------------------------------------------------\n");
    /* Choose size: */
    int32_t NX_I = int32_abrandom(10, 20);
    int32_t NY_I = int32_abrandom(10, 20);
    /* Choose the function: */
    zfunc_t *func = pst_proc_map_function_generic(nfunc);  /* The function that defines the map. */
    double maxGDiff = NAN; /* Threshold for weight map definition. */
    switch(nfunc)
      { 
        case 0:
        case 1:
        case 2:
        case 3:
        case 4:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
        case 16:
          maxGDiff = +INF;
          break;
        case 8:
          maxGDiff = 0.01;
          break;
        case 5:
        case 6:
        case 7:
        case 9:
        case 10:
        case 17:
        case 18:
        case 19:
        case 20:
          maxGDiff = 0.001;
          break;
        default:
          maxGDiff = +INF;
          break;
      }
    /* Create full-size height, slope and weight maps: */
    float_image_t *IZ = NULL;
    float_image_t *IG = NULL;
    float_image_t *IW = NULL;
    make_maps(NX_I, NY_I, func, maxGDiff, &IZ, &IG, &IW);
    /* Shrink it: */
    float_image_t *JG;
    float_image_t *JW;
    pst_slope_and_weight_map_shrink(IG, IW, &JG, &JW);
    /* Make a reduced version of the height, slope, and weight maps: */
    int32_t NX_J = (NX_I + 1)/2;
    int32_t NY_J = (NY_I + 1)/2;
    float_image_t *KZ = NULL;
    float_image_t *KG = NULL;
    float_image_t *KW = NULL;
    make_maps(NX_J, NY_J, func, maxGDiff, &KZ, &KG, &KW);
    /* Compare and write: */
    compare_and_write_slope_maps(nfunc, IZ, IG, IW, JG, JW, KZ, KG, KW);
    fprintf(stderr, "----------------------------------------------------------------------\n");
    float_image_free(JG);
    float_image_free(JW);
    float_image_free(KG);
    float_image_free(KW);
  }
  
void test_r_e_height_map(int32_t nfunc) { }

void compare_and_write_slope_maps
  ( int32_t nfunc,
    float_image_t *IZ, 
    float_image_t *IG, 
    float_image_t *IW, 
    float_image_t *JG, 
    float_image_t *JW, 
    float_image_t *KZ,
    float_image_t *KG,
    float_image_t *KW
  )  
  {
    /* Dimensions of full-size maps: */
    int32_t NX_I = (int32_t)IG->sz[1];
    int32_t NY_I = (int32_t)IG->sz[2];
    
    /* Dimensions of reduced reference map: */
    int32_t NX_K = (int32_t)KG->sz[1];
    int32_t NY_K = (int32_t)KG->sz[2];
    
    /* Write the given height, slope, and weight maps: */
    /* Use the full dimensions to compose the filename: */
    write_map(IZ, nfunc, NX_I, NY_I, "IZ");
    write_map(IG, nfunc, NX_I, NY_I, "IG");
    write_map(IW, nfunc, NX_I, NY_I, "IW");
    write_map(JG, nfunc, NX_I, NY_I, "JG");
    write_map(JW, nfunc, NX_I, NY_I, "JW");
    write_map(KZ, nfunc, NX_I, NY_I, "KZ");
    write_map(KG, nfunc, NX_I, NY_I, "KG");
    write_map(KW, nfunc, NX_I, NY_I, "KW");
    
    /* Compare the reduced gradient maps (shrunk and pristine): */
    float_image_t *DG = float_image_new(2, NX_K, NY_K);
    float_image_mix_channels(1.0,JG,0, -1.0,KG,0, DG,0);
    float_image_mix_channels(1.0,JG,1, -1.0,KG,1, DG,1);
    write_map(DG, nfunc, NX_I, NY_I, "DG");
    
    /* Compare the reduced weight maps (shrunk and pristine): */
    float_image_t *DW = float_image_new(1, NX_K, NY_K);
    float_image_mix_channels(1.0,JW,0, -1.0,KW,0, DW,0);
    write_map(DW, nfunc, NX_I, NY_I, "DW");
    
    float_image_free(DG);
    float_image_free(DW);
  }
  
void write_map(float_image_t *M, int32_t nfunc, int32_t NX,int32_t NY, char *tag)
  { 
    char *fname = jsprintf("out/%02d-%04d-%04d-%s.fni", nfunc, NX, NY, tag);
    FILE *wr = open_write(fname, TRUE);
    float_image_write(wr, M);
    fclose(wr);
    free(fname);
  }

void make_maps(int32_t NX, int32_t NY, zfunc_t func, double maxGDiff, float_image_t **ZP, float_image_t **GP, float_image_t **WP)
  {
    int32_t smoothGL = 2;
    int32_t smoothGN = 11;
    pst_proc_map_sampling_t smpG = pst_proc_map_make_sampling_tables(smoothGL, smoothGN);
    int32_t smoothZL = 1;
    int32_t smoothZN = 11;
    pst_proc_map_sampling_t smpZ = pst_proc_map_make_sampling_tables(smoothZL, smoothZN);
    bool_t numGrad = FALSE;
    double sigmaG = 0.0;
    double sigmaW = 0.0;
    float_image_t *Z = float_image_new(1, NX+1, NY+1); /* Must be non-null. */
    float_image_t *G = (GP == NULL ? NULL : float_image_new(2, NX, NY));
    float_image_t *W = (WP == NULL ? NULL : float_image_new(1, NX, NY));
    pst_proc_map_make_images(func, NX, NY, smpZ, smpG, numGrad, maxGDiff, sigmaG, sigmaW, Z, G, W, NULL);
    if (ZP != NULL) { (*ZP) = Z; }
    if (GP != NULL) { (*GP) = G; }
    if (WP != NULL) { (*WP) = W; }
    pst_proc_map_free_sampling_tables(&smpZ);
    pst_proc_map_free_sampling_tables(&smpG);
  }

options_t *tb_parse_options(int32_t argc, char **argv)
  {
    /* INITIALIZATION: */

    /* Start the command line analyzer {pp}: */
    argparser_t *pp = argparser_new(stderr, argc, argv);
    argparser_set_help(pp, PROG_NAME " version " PROG_VERS ", usage:\n" PROG_HELP);
    argparser_set_info(pp, PROG_INFO);

    /* Process "-help" and "-info" options: */
    argparser_process_help_info_options(pp);

    /* Allocate the return record: */
    options_t *o = (options_t *)malloc(sizeof(options_t));

    /* PARSE KEYWORD ARGUMENTS: */
    argparser_get_keyword(pp, "-function");
    o->function = (int32_t)argparser_get_next_int(pp, pst_proc_map_MIN_ZFUNC, pst_proc_map_MAX_ZFUNC);
    
    /* PARSE POSITIONAL ARGUMENTS: */

    /* Skip to first positional argument: */
    argparser_skip_parsed(pp);

    /* FINALIZATION: */

    /* Check for leftover arguments: */
    argparser_finish(pp);

    return o;
  }
