#define PROG_NAME "test_povray_camera"
#define PROG_DESC "tests the conversion from Tsai camera matrix to POV-Ray camera spec"

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <image.h>
#include <tf_calib.h> 
#include <tf_camera.h> 
#include <tf_targets.h> 
#include <r2x2.h> 
#include <r2.h>
#include <jsfile.h>
#include <sys/stat.h>
#include <sys/types.h>

int32_t main (int32_t argc, char *argv[])
{
  /*Verify memory*/
  void *trash = malloc(1);
  struct mallinfo info;
  int32_t MemDinInicial, MemDinFinal;
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

  FILE *f_cpars = fopen(argv[3], "r");

  int32_t n_p_w;
  tf_calib_data_t *cdat = tf_calib_data_new();
  tf_calib_data_read_world_points (argv[4], &n_p_w, &(cdat->world));


  while (!feof(f_cpars)) {
      char *frame_name = NULL; 
      int32_t frame = tf_camera_params_read (f_cpars, cpar);
      fprintf(stderr, "O frame lido foi %d\n", frame);
      asprintf(&frame_name, "%s%05d.pgm", frame_prefix, frame);
      image_t img = image_read_pgm (frame_name);

      fprintf(stderr, "In camera parameters:\n");
      tf_camera_params_print (cpar, stderr);

      int32_t i;

      mkdir("out", 0777);
      char *out_dir = NULL;
      asprintf(&out_dir, "points/%05d", frame);
      /* Creating the frame directory */
      assert(mkdir(out_dir, 0777)==0);

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

      fprintf(f_pw, "%d\n", n_p_w);
      fprintf(f_p_wgt, "%d\n", n_p_w);

      fprintf(stdout, "######################### Frame: %d\n", frame); 
      fprintf(stdout, "tdat->ntargets: %d\n", n_p_w);
      /* Output Postscript file: */

      for (i = 0; i < n_p_w; i++) {
          fprintf(stdout, "before tf_world_coords_to_image_coords\n"); 
          r2_t p_i = tf_world_coords_to_image_coords (cpar, cdat->world[i]);
          fprintf(stdout, "after tf_world_coords_to_image_coords\n"); 

          int32_t mark_radius = 5;
	  fprintf(stdout, "cross in %f %f image position\n", p_i.c[0], p_i.c[1]);

	  if (   (p_i.c[0] > 5) && (p_i.c[1] > 5)
	      && (p_i.c[0] < img->nx - mark_radius) && (p_i.c[1] < img->ny - mark_radius))
	  image_print_cross (img, mark_radius, p_i); 

	  fprintf(f_pi, "%f %f\n", p_i.c[0], p_i.c[1]);
	  fprintf(f_pw, "%f %f %f\n", cdat->world[i].c[0], cdat->world[i].c[1], cdat->world[i].c[2]);
	  fprintf(f_p_wgt, "%f\n", 1.00000);
 
      }

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

