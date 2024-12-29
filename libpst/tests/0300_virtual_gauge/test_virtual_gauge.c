#define PROG_NAME "compare_gauge_virtual"
#define PROG_DESC "compares a simple lighting model to a spherical light gauge image"

/* Last edited on 2024-12-28 23:58:55 by stolfi */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <jsfile.h>
#include <jsprintf.h>
#include <math.h>
#include <r3.h>
#include <rn.h>
#include <rmxn.h>
#include <ellipse_crs.h>
#include <ellipse_crs_args.h>
#include <affirm.h>
#include <float_image.h>
#include <float_image_read_pnm.h>
#include <float_image_write_pnm.h>
#include <argparser.h>
#include <jsrandom.h>

#include <pst_virtual_gauge.h>

#define PROG_HELP \
  PROG_NAME " \\\n" \
  "      -gauge \\\n" \
  "          image {IMAGE_FNAME} \\\n" \
  "          " ellipse_crs_args_center_HELP " \\\n" \
  "          " ellipse_crs_args_radius_HELP " \\\n" \
  "          " ellipse_crs_args_stretch_HELP " \\\n" \
  "          [  ] \\\n" \
  "-inprefix {FILENAME} \\\n" \
  "           {pointlike,harmonic,radial,glossy,harmonicSG} \\\n" \
  "-outprefix {FILENAME} \\\n" \
  "    " argparser_help_info_HELP ""
  
#define PROG_INFO \
  "    The \"view\" option, if present, must immediately follow the ellipse parameters."

typedef struct options_t
  { pst_virtual_gauge_data_t* gd;
    char* inprefix;
    char* image_in;          /* Filename of image containing the gauge. */
    char* outprefix;
    char* image_out;          /* Filename of image containing the gauge. */
    int32_t magnify;          /* Magnification factor for output image. */
    int32_t margin;           /* Extra margin in pixels for output image. */
    double gamma;
  } options_t;

options_t* parse_options(int32_t argc, char** argv);

float_image_t* read_gauge_image(char* filename, double gamma);

int32_t pixel_position_in_gauge_image(int32_t x, int32_t y, ellipse_crs_t* E);
  /* Position of the pixel whose lower left corner is {(x,y)} relative to the gauge's projection.
    Returns {+1} if totally inside, {-1} if totally outside, and 0 if the pixel straddles
    the projection's outline (or too close to tell). */

// float_image_t* computeWeightGauge(float_image_t* img_in,  pst_virtual_gauge_data_t * gdata,void** l_data, pst_light_t* lm);
// float_image_t* computeWeightGauge(float_image_t* img_in,  pst_virtual_gauge_data_t * gdata,void** l_data, pst_light_t* lm){
//   int32_t NC   = img_in->sz[0];
//   int32_t NXot = img_in->sz[1];
//   int32_t NYot = img_in->sz[2];
//   float_image_t* img_sy = float_image_new(NC, NXot, NYot);
//   pst_vitual_gauge_paint(img_sy, &(gdata->E), &(gdata->view), l_data,lm->data_weight);
//   return img_sy;
// }

float_image_t* read_gauge_image(char* filename, double gamma)
  { char* ext = strrchr(filename, '.');
    demand(ext != NULL, "image name has no extension");
    float_image_t* img;
    if ((strcmp(ext, ".ppm") == 0) || (strcmp(ext, ".pgm") == 0))
      { bool_t isMask = FALSE;
        double bias = 0.0; /* For now */
        bool_t yup = TRUE;
        bool_t warn = FALSE;
        bool_t verbose = TRUE;
        img = float_image_read_pnm_named(filename, isMask, gamma, bias, yup, warn, verbose); }
    else if (strcmp(ext, ".fni") == 0)
      { FILE* rd = open_read(filename, TRUE);
        img = float_image_read(rd);
        if (rd != stdin) { fclose(rd); }
      }
    else
     { demand(FALSE, "unknown image name extension"); }
    return img;
  }

int32_t pixel_position_in_gauge_image(int32_t x, int32_t y, ellipse_crs_t* E)
  {
    double debug = FALSE;
    if (debug) { fprintf(stderr, "position of pixel [%d %d]", x, y); }
    
    /* Find the center {p} of pixel {x,y}: */
    r2_t p = (r2_t){{ x + 0.5, y + 0.5 }};
    /* If the four corners are inside, the pixel is inside: */
    bool_t inside = TRUE;
    int32_t kx, ky;
    for (kx = -1; (kx <= 1) && inside; kx += 2)
      { for (ky = -1; (ky <= 1) && inside; ky += 2)
          { /* Get corner, slightly outwards of pixel: */
            r2_t q = p;
            q.c[0] += 0.50001*kx;
            q.c[1] += 0.50001*ky;
            if (debug) { r2_gen_print(stderr, &q, "%12.7f", "  q = (", " ", ")\n"); }
            inside &= ellipse_crs_inside(E, &q);
          }
      }
    if (inside) { return +1; }
    
    /* Else, if {p} is outside the fattened ellipse, the pixel is outside: */
    ellipse_crs_t EF = (ellipse_crs_t){ .ctr = E->ctr, .rad = E->rad + 0.707107, .str = E->str };
    bool_t outside = ! ellipse_crs_inside(&EF, &p);
    if (outside) { return -1; }
      
    /* Else give up: */
    return 0;
  }

double gauge_coverage_of_pixel(int32_t x, int32_t y, ellipse_crs_t* E);
  /* Returns 1.0 if the pixel with lower left corer {(x,y)} is entirely inside the gauge's
    projection; 0.0 if it is entirely outside; and a number strictly between 0 and 1 
    if the pixel straddles the projection's boundary. */
    
double gauge_coverage_of_pixel(int32_t x, int32_t y, ellipse_crs_t* E)
  {
    double debug = FALSE;
    int32_t pos = pixel_position_in_gauge_image(x, y, E);
    if (debug) { fprintf(stderr, "pixel [%d %d]  position = %+d ", x, y, pos); }
    double cov;
    if (pos < 0)
      { /* Fully outside */
        if (debug) { fprintf(stderr, " (outside)\n"); }
        cov = 0.0;
      }
    else if (pos > 0)
      { /* Fully inside */
        if (debug) { fprintf(stderr, " (inside)\n"); }
        cov = 1.0;
      }
    else
      { /* Straddling the border; find {cov} by sampling: */
        if (debug) { fprintf(stderr, " (straddles)"); }
        int32_t NS = 15;  /* Pixel subsampling points along each axis for coverage. */
        int32_t kx, ky;
        int32_t nin = 0; /* Counts subsampling points inside the pixel. */
        for (kx = 0; kx < NS; kx++)
          { for (ky = 0; ky < NS; ky++)
              { /* Pick a subsampling point in pixel: */
                r2_t q = (r2_t){{ x + (kx + 0.5)/NS, y + (ky + 0.5)/NS }};
                /* Check whether {q} is inside the sphere's projection: */
                bool_t inside = ellipse_crs_inside(E, &q);
                if (inside) { nin++; }
              }
          }
        if (debug) { fprintf(stderr, " nin = %d", nin); }
        /* Make sure that {cov} is sufficiently away from 0 and 1 to avoid confusion: */
        double cov_min = 1.0/255.0;
        double cov_max = 254.0/255.0;
        cov = ((double)nin + 1)/(NS*NS + 2);
        if (cov < cov_min) { cov = cov_min; }
        if (cov > cov_max) { cov = cov_max; }
      }
    if (debug) { fprintf(stderr, " cov = %8.6f\n", cov); }
    return cov;
  }


float_image_t* make_gauge_mask ( float_image_t* img, ellipse_crs_t* E, r3_t* view);
float_image_t* make_gauge_mask ( float_image_t* img, ellipse_crs_t* E, r3_t* view)
  {
    /* Get the image dimensions: */
    int32_t NX = (int32_t)img->sz[1];
    int32_t NY = (int32_t)img->sz[2];
    
    float_image_t* xtr = float_image_new(1,NX,NY);
   
    
   
    int32_t x, y;
    for (x = 0; x < NX; x++) 
      { for (y = 0; y < NY; y++) 
          { /* Compute fractional coverage {cov}: */
            double cov = gauge_coverage_of_pixel(x, y, E);
            float_image_set_sample(xtr, 0,x,y, (float)cov);
	    
          }
      }
              
    return xtr;
  }

int32_t main(int32_t argc, char** argv){
  
  options_t* o = parse_options(argc,argv);
  
  /*First we init the light model to be tested*/
 
  
  fprintf(stderr,"Reading Gauge Image...");
  float_image_t* img_in = read_gauge_image(o->image_in, o->gamma);
  fprintf(stderr,"OK.\n");
  
  int32_t  NC, NX, NY;
  float_image_get_size(img_in, &NC, &NX, &NY);
  
  /* pst_light_t* lht_test = create_approx_lighting_model(o->in_lht_data.lightingModelType); */
  
  //char* arq_light_param_name = jsprintf("%s",o->inprefix);
  //FILE* arq_light_param = open_read(arq_light_param_name,TRUE);
  /* void** l_datatest = (void**)malloc(sizeof(void*)*NC); */
  //for(int32_t c = 0; c < NC; c++){
  //  l_datatest[c] = lht_test->read_parameters(arq_light_param);
  //  fprintf(stderr,"Lighting model chanel %d\n",c);
  //  fprintf(stderr,"************************************\n");
  //  lht_test->write_parameters(stderr,l_datatest[c]);
  //  fprintf(stderr,"************************************\n");
  //  
  //}
  //fclose(arq_light_param);
  
  fprintf(stderr,"making mask image...\n");
  float_image_t* img_mk = make_gauge_mask(img_in,&(o->gd->E),&(o->gd->view));
  
  auto double shading(r3_t* nrm, int32_t c);
    /* Shading function for the virtual gauge. */
  
  fprintf(stderr,"making virtual gauge image...\n");
  float_image_t* img_sy = pst_shading_make_image(NC, NX, NY, o->gd, &shading, &nbad, &nneg, &nbig);
  if (nneg[c] > 0) { fprintf(stderr, "** warning: %d NAN/INF samples in synthetic image.\n", nnbad; }
  if (nneg[c] > 0) { fprintf(stderr, "** warning: %d negative samples in synthetic image.\n", nneg); }
  if (nbig[c] > 0) { fprintf(stderr, "** warning: %d samples over 1.0 synthetic image.\n", nbig); }
  
  /*  fprintf(stderr,"making weight image...\n");
  float_image_t* img_wgt = computeWeightGauge( img_in, o->gd,l_datatest, lht_test);
  fprintf(stderr,"OK.\n");*/
  
  fprintf(stderr,"computing error image...");
  float_image_t* img_nr = NULL; /* For now */
  float_image_t* img_err = pst_shading_make_image_error(img_in, img_sy, img_mk, img_nr);
  
  char* wr_err_name = jsprintf("%s_err.fni", o->outprefix);
  float_image_write_named(wr_fni_err,img_err);
  free(wr_err_name);
  
//   char* _wr_wgt_name = NULL;
//   char* _wr_wgt_name = jsprintf("%s_wgt.fni",o->outprefix);
//   FILE* wr_fni_wgt = open_write(_wr_wgt_name,TRUE);
//   float_image_write(wr_fni_wgt,img_wgt);
//   fclose(wr_fni_wgt);
  
  char* wr_sy_fni_name = jsprintf("%s_sy.fni",o->outprefix);
  float_image_write_named(wr_sy_fni, img_sy);
  free(wr_sy_fni_name);
  
  char* wr_sy_pnm_name = jsprintf("%s_sy.ppm",o->outprefix);
  float_image_pnm_write_named(wr_sy_pnm_name, img_sy, FALSE, o->gamma, 0.0, TRUE,TRUE,TRUE);
  free(wr_sy_pnm_name);
 
  return 0;
}

void check_dup(argparser_t* pp, bool_t dup, char* key);
  /* If {dup} is true, prints an error "duplicate keyword {key}" and aborts. */
  

pst_virtual_gauge_data_t* parse_gauge_args(argparser_t* pp, bool_t input)
  {
    pst_virtual_gauge_data_t* gd = (pst_virtual_gauge_data_t*)malloc(sizeof(pst_virtual_gauge_data_t));
    ellipse_crs_t* E = &(gd->E);

    int32_t c; 
    for (c = 0; c < 3; c++) { gd->albedo[c] = NAN; }
    
    /* Parse fields of a gauge spec: */
    bool_t done;
    do 
      { 
        done = FALSE;
        if (argparser_keyword_present_next(pp, "image"))
          { check_dup(pp, gd->image_{in|out} != NULL, "image");
            gd->image_{in|out} = argparser_get_next(pp);
          }
        else 
        else if (input)
          { /* Fields allowed only in input gauge spec: */
            pst_virtual_gauge_args_parse(pp, &gd);
            
          }
        else
          { /* Fields allowed only in output gauge spec: */
            if (argparser_keyword_present_next(pp, "magnify"))
              {
                check_dup(pp, gd->magnify > 0, "magnify");
                gd->magnify = argparser_get_next_int32_t(pp, 1, 64);
              }
            else if (argparser_keyword_present_next(pp, "margin"))
              {
                check_dup(pp, gd->margin >= 0, "margin");
                gd->margin = argparser_get_next_int32_t(pp, 0, 1024);
              }
            else
              { done = TRUE; }
          }
      }
    while (! done);
      
    /* Check for required fields and provide defaults: */
    if (gd->image_{in|out} == NULL) { argparser_error(pp, "missing \"image\" for gauge."); }
    if (isnan(gd->albedo[0])) 
      { for (c = 0; c < 3; c++) { gd->albedo[c] = 1.0; } }
    
    if (input)
      { /* Fields allowed/required only in input gauge spec: */
        if (isnan(E->rad)) { argparser_error(pp, "missing \"radius\" for gauge."); }
        if (isnan(E->ctr.c[0])) { argparser_error(pp, "missing \"center\" for gauge."); }
        if (isnan(E->str.c[0])) { E->str = (r2_t) {{ 0.0, 0.0}}; }
        if (isnan(gd->view.c[0])) { gd->view = (r3_t) {{ 0.0, 0.0, 1.0}}; }
      }
    else
      { /* Fields allowed/required only in output gauge spec: */
        if (gd->magnify < 0) { gd->magnify = 1; }
        if (gd->margin < 0) { gd->margin = 0; }
      }
      
    return gd;
  }


options_t* parse_options(int32_t argc, char** argv){
  options_t* o = (options_t*)malloc(sizeof(options_t));
  argparser_t* pp = argparser_new(stderr, argc, argv);
  argparser_set_help(pp, PROG_HELP);
  argparser_set_info(pp, PROG_INFO);
  argparser_process_help_info_options(pp);
    
  o->image_{in|out} = NULL;
  o->magnify = -1;
  o->margin = -1;
  
  argparser_get_keyword(pp, "-gauge");
    o->gd = parse_gauge_args(pp, TRUE);
  
  argparser_get_keyword(pp, "-appModel");
  o->inprefix = argparser_get_next(pp);
  fprintf(stderr,"INPREFIX %s\n",o->inprefix);
  o->in_lht_data = lighting_parse(pp,FALSE);
  
  
  argparser_get_keyword(pp, "-outprefix");
  o->outprefix = argparser_get_next(pp);
  fprintf(stderr,"OUTPREFIX %s\n",o->outprefix);
  
  if (argparser_keyword_present(pp, "-gamma"))
      { o->gamma = argparser_get_next_double(pp, 0.0, 10.0); }
    else
      { o->gamma = 1.0; }
  fprintf(stderr,"GAMMA %lf\n",o->gamma);
   

  argparser_finish(pp);
  return o;

}
 
void check_dup(argparser_t* pp, bool_t dup, char* key)
  { if (dup)
      { char* msg = jsprintf("duplicate keyword \"%s\"", key);
        argparser_error(pp, msg);
        free(msg);
      }
  }
