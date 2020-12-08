/* Last edited on 2017-06-22 18:18:49 by stolfilocal */

#define PROG_NAME "test_povray_camera"
#define PROG_DESC "tests the conversion from Tsai camera matrix to POV-Ray camera spec"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <r2x2.h> 
#include <affirm.h> 
#include <r2.h>
#include <jsfile.h>
#include <pswr.h>
#include <float_image.h>
#include <float_image_io_pnm.h>

#include <tf_calib.h> 
#include <tf_camera.h> 
#include <tf_targets.h> 


void my_plot_marks
  ( PSStream *ps, 
    tf_camera_params_t *cpar, 
    tf_calib_data_t * tdat,
    bool_t cross, /*TRUE paints crosses, FALSE paints dots.*/
    bool_t red /*TRUE paints red crosses, FALSE paints black ones.*/
  );
/* 
  Plots all marks to the Postscript stream {ps}.
  If {cpar} is NULL, plots the given image coordinates {cdat->image},
  else computes the image coords from {cdat->world} using {cpar}
  and plots the result. */

void pswr_cross(PSStream *ps, double xc, double yc, double radius);
  /* Draws a cross at {(xc,yc)} with the given radius in mm.
   !!! Should be in {pswr.h}. !!! */



int main (int argc, char *argv[])
{
  /*Verify memory*/
  void *trash = malloc(1);
  struct mallinfo info;
  int MemDinInicial, MemDinFinal;
  free(trash);
  info = mallinfo();
  MemDinInicial = info.uordblks;
  /*End*/

  char *frame_prefix = argv[1];

  /* Get the fixed camera parameters: */
  tf_camera_specs_t *cspec;

  if (strcmp(argv[2],"optura") == 0) {
      cspec = tf_camera_specs_for_canon_optura ();
  }
  else if (strcmp(argv[2],"svga") == 0) {
      cspec = tf_camera_specs_for_povray_svga (); 
  }
  else if (strcmp(argv[2],"hvga") == 0) {
      cspec = tf_camera_specs_for_povray_hvga ();
  }
  else {
     fprintf(stderr, "error: define the camera specifications\n");
     exit(1);
  }  

  tf_camera_params_t *cpar = tf_camera_specs_get_new_mean_params(cspec);
  tf_camera_params_t *ctrue = tf_camera_specs_get_new_mean_params(cspec);

  FILE *f_cpars = fopen(argv[3], "r");

  int n_p_w;
  tf_calib_data_t *cdat = tf_calib_data_new();
  tf_calib_data_read_world_points (argv[4], &n_p_w, &(cdat->world));


  FILE *f_true;
  if (strcmp(argv[5], "none") != 0) {
      f_true = fopen(argv[5], "r");
  }    

  while (!feof(f_cpars)) {
      char *frame_name = NULL; 
      //int frame = tf_get_cpar_from_file (cpar,  f_cpars);
      //int frame = tf_camera_params_read_mutable (f_cpars, cpar);
      int frame = tf_camera_params_read (f_cpars, cpar);
      fprintf(stderr, "O frame lido foi %d\n", frame);
      asprintf(&frame_name, "%s/%05d/frame.pgm", frame_prefix, frame);
      float_image_t *img = float_image_read_pnm_named(frame_name);

      fprintf(stderr, "In camera parameters:\n");
      tf_camera_params_print (cpar, stderr);

      int i;

      mkdir("out", 0777);
      char *out_dir = NULL;
      asprintf(&out_dir, "points/%05d", frame);
      /* Creating the frame directory */
      assert(mkdir(out_dir, 0777)==0);

      if (strcmp(argv[5], "none") != 0) {
          tf_camera_params_read_mutable (f_true, ctrue);

          FILE *f_true_out;
          char *f_true_name;
          asprintf(&f_true_name, "%s/true.cpar", out_dir);
          f_true_out = fopen(f_true_name, "w");
          free(f_true_name); 
	  tf_camera_params_write_mutable (f_true_out, frame, ctrue);
	  fclose(f_true_out);
      }

      FILE *f_pi;
      char *out_f_pi_name;
      asprintf(&out_f_pi_name, "%s/p_i.txt", out_dir);
      f_pi = fopen(out_f_pi_name, "w");
      free(out_f_pi_name); 

      FILE *f_pw;
      char *out_f_pw_name;
      asprintf(&out_f_pw_name, "%s/p_w.txt", out_dir);
      f_pw = fopen(out_f_pw_name, "w");
      free(out_f_pw_name); 

      FILE *f_p_wgt;
      char *out_f_p_wgt_name;
      asprintf(&out_f_p_wgt_name, "%s/p_wgt.txt", out_dir);
      f_p_wgt = fopen(out_f_p_wgt_name, "w");
      free(out_f_p_wgt_name); 

      fprintf(f_pi, "%d\n", n_p_w);
      fprintf(f_pw, "%d\n", n_p_w);
      fprintf(f_p_wgt, "%d\n", n_p_w);

      fprintf(stdout, "######################### Frame: %d\n", frame); 
      fprintf(stdout, "tdat->ntargets: %d\n", n_p_w);
      /* Output Postscript file: */

      /* Camera image dimensions and aspect ratio: */
      double Nx = interval_mid(&cspec->Npx);
      double Ny = interval_mid(&cspec->Npy);
      double aspect = Nx/Ny;

      /* Postscript figure sizeand margin, in pt: */
      double hSize = 300.0*sqrt(aspect);
      double vSize = 300.0/sqrt(aspect);
      double hMrg = 4;
      double vMrg = 3;

      char *ps_prefix = NULL;
      asprintf(&ps_prefix, "out/fig");
      PSStream *ps = pswr_new_stream(TRUE, ps_prefix, NULL, NULL, FALSE, hSize + 2*hMrg, vSize + 2*vMrg);
      free(out_dir);
   
      pswr_new_canvas(ps, "all");
      /* Negate all Y coordinates to get the proper orientation for the Y axis: */
      pswr_set_window(ps, 0.0, Nx,  -Ny, 0.0, hMrg, hSize+hMrg, vMrg, vSize+vMrg);
      /* Draw a frame around the image: */
      pswr_set_pen(ps, 0.0,0.0,0.0, 0.25, 0,0); /* Black lines */
      pswr_frame(ps); 

      for (i = 0; i < n_p_w; i++) {
          fprintf(stdout, "before tf_world_coords_to_image_coords\n"); 
          r2_t p_i = tf_world_coords_to_image_coords (cpar, cdat->world[i]);
	  cdat->image[i] = p_i;
          fprintf(stdout, "after tf_world_coords_to_image_coords\n"); 

          int mark_radius = 5;
	  fprintf(stdout, "cross in %f %f image position\n", p_i.c[0], p_i.c[1]);

	  if (   (p_i.c[0] > 5) && (p_i.c[1] > 5)
	      && (p_i.c[0] < img->nx - mark_radius) && (p_i.c[1] < img->ny - mark_radius))
	  image_print_cross (img, mark_radius, p_i); 

	  fprintf(f_pi, "%f %f\n", p_i.c[0], p_i.c[1]);
	  fprintf(f_pw, "%f %f %f\n", cdat->world[i].c[0], cdat->world[i].c[1], cdat->world[i].c[2]);
	  fprintf(f_p_wgt, "%f\n", 1.00000);
 
      }

     my_plot_marks(ps, NULL, cdat, FALSE, FALSE); /* Given coordinates. */

      pswr_close_stream(ps);


      fclose(f_pi);
      fclose(f_pw);
      fclose(f_p_wgt);

      fprintf(stdout, "Loop outside\n"); 

      char *cross_fname = NULL; 
      asprintf(&cross_fname, "out/cross_%05d.pgm", frame);
      image_write_pgm (img, PGM_RAW_MAGIC, MAXVAL, cross_fname);
      free(cross_fname);

      free(frame_name);
      image_free(img);
  }

  if (strcmp(argv[5], "none") != 0) {
     fclose(f_true);
  }

  info = mallinfo();
  fclose(f_cpars);
  MemDinFinal = info.uordblks;
  if (MemDinInicial != MemDinFinal)
    { fprintf(stderr, "*** main: memory leak %d blocks\n", MemDinFinal - MemDinInicial); }
  return 0;
}


void my_plot_marks
  ( PSStream *ps, 
    tf_camera_params_t *cpar, 
    tf_calib_data_t * cdat,
    bool_t cross, /*TRUE paints crosses, FALSE paints dots.*/
    bool_t red /*TRUE paints red crosses, FALSE paints black ones.*/
  )
  {    
    double cross_radius = 1.5;
    double dot_radius = 1.5;
    if (red) {
      pswr_set_pen(ps, 1.0,0.0,0.0, 0.5, 0,0); /* Red lines */
    } else {
      pswr_set_pen(ps, 0.0,0.0,0.0, 0.5, 0,0); /* Black lines */
    }
    pswr_set_fill_color(ps, 1.0,1.0,0.0);    /* Yellow fill. */
    int mark;
    for (mark = 0; mark < cdat->np; mark++) {
      r3_t p_w = cdat->world[mark];
      r2_t p_i;
      if (cpar == NULL) {
        p_i = cdat->image[mark];
      } else {
         p_i = tf_world_coords_to_image_coords(cpar, p_w);
      }
      if (cross) {
        pswr_cross(ps, p_i.c[0], -p_i.c[1], cross_radius);
      } else {
        pswr_dot(ps, p_i.c[0], -p_i.c[1], dot_radius, TRUE, TRUE);
      }
    }
  }

void pswr_cross(PSStream *ps, double xc, double yc, double radius)
  {
    pswr_tic(ps, 0, xc, yc, 2*radius, 0.5);
    pswr_tic(ps, 1, xc, yc, 2*radius, 0.5);
  }

