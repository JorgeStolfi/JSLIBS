/* See {tf_camera.h}. */
/* Last edited on 2022-10-20 05:53:22 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <jsfile.h>
#include <affirm.h>
#include <ctype.h>
#include <r3.h>
#include <r2.h>
#include <r4x4.h>
#include <tf_camera.h>

tf_camera_params_t *tf_camera_params_new (void)
{
  return (tf_camera_params_t *)notnull(malloc(sizeof(struct tf_camera_params_t)), "no mem");
}

tf_camera_params_t *tf_camera_params_copy (tf_camera_params_t *cpar)
{
  tf_camera_params_t *cpar_new = tf_camera_params_new();
  (*cpar_new) = (*cpar);
  return cpar_new;
}

void tf_camera_params_print (tf_camera_params_t *cpar, FILE *ferr)
{
    fprintf(ferr, "// ------------ camera parameters --------------------\n");
    fprintf(ferr, "// Npx = %f  Npy = %f\n", cpar->Npx, cpar->Npy);
    fprintf(ferr, "// dpx = %.5f  dpy = %.5f\n", cpar->dpx, cpar->dpy);
    fprintf(ferr, "// Cx = %.3f  Cy = %.3f\n", cpar->Cx, cpar->Cy);
    fprintf(ferr, "// sx = %.6f\n", cpar->sx);
    tf_camera_matrix_print(&(cpar->S), ferr);
    fprintf(ferr, "// f     = %15.15f\n", cpar->f);
    fprintf(ferr, "// kappa = %15.15f\n", cpar->kappa);
    fprintf(ferr, "// -----------------------------------------------------------\n");
}

void tf_camera_matrix_print(r4x4_t *S, FILE *ferr)
{
    fprintf(ferr, "// Transformatiom Matrix =\n");
    int32_t i, j;
    for (i = 0; i < 4; i++) {
      fprintf(ferr, "//   ");
      for (j = 0; j < 4; j++) {
            fprintf(ferr, "%+10.15f \t", S->c[i][j]);
        }
      fprintf(ferr, "\n");
    }
}

void tf_camera_matrix_print_rotation (r4x4_t *S, FILE *ferr)
{
    int32_t i, j;
    for (i = 1; i < 4; i++) {
      fprintf(ferr, "//   ");
      for (j = 1; j < 4; j++) {
          fprintf(ferr, "%+10.15f \t", S->c[i][j]);
      }
      fprintf(ferr, "\n");
   }  
}
void tf_camera_params_write (FILE *wr, int32_t index, tf_camera_params_t *cpar)
{
  fprintf(wr, "%d", index);
    
  /*writing rotating parameters*/
  fprintf(wr, "  %22.17f %22.17f %22.17f", cpar->S.c[1][1], cpar->S.c[1][2], cpar->S.c[1][3]);
  fprintf(wr, "  %22.17f %22.17f %22.17f", cpar->S.c[2][1], cpar->S.c[2][2], cpar->S.c[2][3]);
  fprintf(wr, "  %22.17f %22.17f %22.17f", cpar->S.c[3][1], cpar->S.c[3][2], cpar->S.c[3][3]);

  /*writing tranlating parameters*/
  fprintf(wr, "  %30.17f %30.17f %30.17f", cpar->S.c[1][0], cpar->S.c[2][0], cpar->S.c[3][0]);

  /*writing focal distance and kappa distortion*/ 
  fprintf(wr, "  %30.17f  %22.17f", cpar->f, cpar->kappa);

  /*writing the horizontal scale factor*/
  fprintf(wr, "  %22.17f", cpar->sx);

  /*writing the sensor aray parameters*/
  fprintf(wr, "  %15f", cpar->Npx);
  fprintf(wr, "  %15f", cpar->Npy);
  fprintf(wr, "  %22.17f", cpar->dpx);
  fprintf(wr, "  %22.17f", cpar->dpy);
  fprintf(wr, "  %30.17f", cpar->Cx);
  fprintf(wr, "  %30.17f", cpar->Cy);

  fprintf(wr, "\n");

  fflush(wr);

}

int32_t tf_camera_params_read (FILE *rd, tf_camera_params_t *cpar)
{ 
    int32_t frame;
    /*Basic cpar initializations*/ 
    cpar->S.c[0][0] = 1.0; cpar->S.c[0][1] = 0.0;  
    cpar->S.c[0][2] = 0.0; cpar->S.c[0][3] = 0.0; 
 
    /* Getting the camera rotation parameters */ 
    fscanf(rd, "%d", &frame); /*r1*/ 
    fscanf(rd, "%lf", &cpar->S.c[1][1]); /*r1*/ 
    fscanf(rd, "%lf", &cpar->S.c[1][2]); /*r2*/ 
    fscanf(rd, "%lf", &cpar->S.c[1][3]); /*r3*/  
    fscanf(rd, "%lf", &cpar->S.c[2][1]); /*r4*/  
    fscanf(rd, "%lf", &cpar->S.c[2][2]); /*r5*/  
    fscanf(rd, "%lf", &cpar->S.c[2][3]); /*r6*/  
    fscanf(rd, "%lf", &cpar->S.c[3][1]); /*r7*/  
    fscanf(rd, "%lf", &cpar->S.c[3][2]); /*r8*/  
    fscanf(rd, "%lf", &cpar->S.c[3][3]); /*r9*/  
 
    /* Getting the camera translation parameters */ 
    fscanf(rd, "%lf", &cpar->S.c[1][0]); /*Tx*/  
    fscanf(rd, "%lf", &cpar->S.c[2][0]); /*Ty*/  
    fscanf(rd, "%lf", &cpar->S.c[3][0]); /*Tz*/  
         
    /* Getting focal distance */ 
    fscanf(rd, "%lf", &cpar->f);  
 
    /* Getting radial distortion */ 
    fscanf(rd, "%lf", &cpar->kappa);  
 
    /* Getting sx */ 
    fscanf(rd, "%lf", &cpar->sx); 
 
    /* Getting the sensor array parameters */ 
    fscanf(rd, "%lf", &cpar->Npx); 
    fscanf(rd, "%lf", &cpar->Npy); 
    fscanf(rd, "%lf", &cpar->dpx); 
    fscanf(rd, "%lf", &cpar->dpy); 
    fscanf(rd, "%lf", &cpar->Cx); 
    fscanf(rd, "%lf", &cpar->Cy); 
    return frame;
}

void tf_camera_params_write_mutable (FILE *wr, int32_t index, tf_camera_params_t *cpar)
{
  fprintf(wr, "%d", index);
    
  /*writing rotating parameters*/
  fprintf(wr, "  %22.17f %22.17f %22.17f", cpar->S.c[1][1], cpar->S.c[1][2], cpar->S.c[1][3]);
  fprintf(wr, "  %22.17f %22.17f %22.17f", cpar->S.c[2][1], cpar->S.c[2][2], cpar->S.c[2][3]);
  fprintf(wr, "  %22.17f %22.17f %22.17f", cpar->S.c[3][1], cpar->S.c[3][2], cpar->S.c[3][3]);

  /*writing tranlating parameters*/
  fprintf(wr, "  %30.17f %30.17f %30.17f", cpar->S.c[1][0], cpar->S.c[2][0], cpar->S.c[3][0]);

  /*writing focal distance and kappa distortion*/ 
  fprintf(wr, "  %30.17f  %22.17f", cpar->f, cpar->kappa);

  /*writing the horizontal scale factor*/
  fprintf(wr, "  %22.17f", cpar->sx);

  fprintf(wr, "\n");

  fflush(wr);

}

int32_t tf_camera_params_read_mutable (FILE *rd, tf_camera_params_t *cpar)
{ 
    int32_t frame;
    /*Basic cpar initializations*/ 
    cpar->S.c[0][0] = 1.0; cpar->S.c[0][1] = 0.0;  
    cpar->S.c[0][2] = 0.0; cpar->S.c[0][3] = 0.0; 
 
    /* Getting the camera rotation parameters */ 
    fscanf(rd, "%d", &frame); /*r1*/ 
    fscanf(rd, "%lf", &cpar->S.c[1][1]); /*r1*/ 
    fscanf(rd, "%lf", &cpar->S.c[1][2]); /*r2*/ 
    fscanf(rd, "%lf", &cpar->S.c[1][3]); /*r3*/  
    fscanf(rd, "%lf", &cpar->S.c[2][1]); /*r4*/  
    fscanf(rd, "%lf", &cpar->S.c[2][2]); /*r5*/  
    fscanf(rd, "%lf", &cpar->S.c[2][3]); /*r6*/  
    fscanf(rd, "%lf", &cpar->S.c[3][1]); /*r7*/  
    fscanf(rd, "%lf", &cpar->S.c[3][2]); /*r8*/  
    fscanf(rd, "%lf", &cpar->S.c[3][3]); /*r9*/  
 
    /* Getting the camera translation parameters */ 
    fscanf(rd, "%lf", &cpar->S.c[1][0]); /*Tx*/  
    fscanf(rd, "%lf", &cpar->S.c[2][0]); /*Ty*/  
    fscanf(rd, "%lf", &cpar->S.c[3][0]); /*Tz*/  
         
    /* Getting focal distance */ 
    fscanf(rd, "%lf", &cpar->f);  
 
    /* Getting radial distortion */ 
    fscanf(rd, "%lf", &cpar->kappa);  
 
    /* Getting sx */ 
    fscanf(rd, "%lf", &cpar->sx); 
    return frame;
}

void tf_camera_params_print_changes (r3_t *Do, r3x3_t *DR, double Df, double Dkappa, FILE *ferr)
{
    int32_t i, j;
    fprintf(ferr, "%s\n", "Camera displacement:");
    for (i = 0; i < 3; i++) { fprintf(ferr, "%+10.15f \t", Do->c[i]); }
    fprintf(ferr, "\n");
    fprintf(ferr, "%s\n", "Camera orientation matrix:");
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            fprintf(ferr, "%+10.15f \t", DR->c[i][j]);
        }
        fprintf(ferr, "\n");
    }
    fprintf(ferr, "Df     = %15.15f\n", Df);
    fprintf(ferr, "Dkappa = %15.15f\n", Dkappa);
}

r3_t tf_world_coords_to_camera_coords (tf_camera_params_t *cpar, r3_t pw)
{
    r3_t pc;

    pc.c[0] = cpar->S.c[1][1] * pw.c[0] +
              cpar->S.c[1][2] * pw.c[1] +
              cpar->S.c[1][3] * pw.c[2] + cpar->S.c[1][0];

    pc.c[1] = cpar->S.c[2][1] * pw.c[0] +
              cpar->S.c[2][2] * pw.c[1] +
              cpar->S.c[2][3] * pw.c[2] + cpar->S.c[2][0];

    pc.c[2] = cpar->S.c[3][1] * pw.c[0] +
              cpar->S.c[3][2] * pw.c[1] +
              cpar->S.c[3][3] * pw.c[2] + cpar->S.c[3][0];
    return pc;
}

r2_t tf_camera_coords_to_und_sensor_coords (tf_camera_params_t *cpar, r3_t pc)
{
    r2_t pu; 
    pu.c[0] = cpar->f * pc.c[0] / pc.c[2];
    pu.c[1] = cpar->f * pc.c[1] / pc.c[2];
    return pu; 
}

r2_t tf_und_sensor_coords_to_dis_sensor_coords (tf_camera_params_t *cpar, r2_t pu)
{
  r2_t pd = tf_camera_apply_kappa(&pu, cpar->kappa);
  return pd;
}

r2_t tf_dis_sensor_coords_to_und_sensor_coords (tf_camera_params_t *cpar, r2_t pd)
{
  r2_t pu = tf_camera_apply_kappa(&pd, -cpar->kappa);
  return pu;
}

r2_t tf_camera_apply_kappa(r2_t *p, double kappa)
{
  if (kappa == 0.0) {
    /*Ver sem sentido !!!!*/
    //r2_t *q;
    //(*q) = (*p);
    return *p;
  } else {
    /* Get the coordinates of {p}: */
    double pX = p->c[0];
    double pY = p->c[1];

    /* Compute distance squared {s} from optical axis: */
    double s = pX*pX + pY*pY;

    /* Compute the scaling factor {f}: */
    double m = 2*kappa*s;
    if (m >= 0.995) {
      /* Point is outside the maximum radius {R}: */
      m = 0.995;
    }
    double h = 1.0/(1.0 - m);
    double f = sqrt(h);

    /* Apply the radial scaling by {f}: */
    return (r2_t){{ pX*f, pY*f }};
  }
}
  
double tf_camera_maximum_safe_kappa (tf_camera_params_t *cpar)
  { 
    /* Find the max distance squared {R2} from any corner to optical axis: */
    double R2 = 0.0;
    int32_t sx, sy;
    for (sx = 0; sx <= 1; sx++) {
      for (sy = 0; sy <= 1; sy ++) {
        double duX = sx*cpar->Npx*cpar->dpx - cpar->Cx;
        double duY = sy*cpar->Npy*cpar->dpy - cpar->Cy;
        double du2 = duX*duX + duY*duY;
        if (du2 > R2) { R2 = du2; }
      }
    }

    /* Throw in a safety factor of 2.0 in the radius: */
    R2 = 2*2*R2;

    /* Compute the max {kappa} such that {2*kappa*R2 <= 0.995}: */
    double max_kappa = 0.995/R2/2;

    return max_kappa;
  }


r2_t tf_sensor_coords_to_image_coords (tf_camera_params_t *cpar, r2_t pd)
{
    r2_t pi;
    pi.c[0] = pd.c[0] * cpar->sx/cpar->dpx + cpar->Cx;
    pi.c[1] = pd.c[1]/cpar->dpy + cpar->Cy;
    return pi;
}

r2_t tf_image_coords_to_sensor_coords (tf_camera_params_t *cpar, r2_t pi)
{
    r2_t pd;
    pd.c[0] = cpar->dpx * (pi.c[0] - cpar->Cx) / cpar->sx;
    pd.c[1] = cpar->dpy * (pi.c[1] - cpar->Cy);
    return pd;
}

r2_t tf_world_coords_to_image_coords (tf_camera_params_t *cpar, r3_t pw)
{
   bool_t debug = FALSE;
   if (debug) { tf_camera_params_print (cpar, stdout); }
   if (debug) { fprintf(stdout, "  world (mm):       ( %10.6f %10.6f %10.6f )\n", pw.c[0], pw.c[1], pw.c[2]); }
   r3_t pc = tf_world_coords_to_camera_coords (cpar, pw);
   if (debug) { fprintf(stdout, "  camera (mm):      ( %10.6f %10.6f %10.6f )\n", pc.c[0], pc.c[1], pc.c[2]); }
   r2_t pu = tf_camera_coords_to_und_sensor_coords (cpar, pc);
   if (debug) { fprintf(stdout, "  undistorted (mm): ( %10.6f %10.6f )\n", pu.c[0], pu.c[1]); }
   r2_t pd = tf_und_sensor_coords_to_dis_sensor_coords (cpar, pu);   
   if (debug) { fprintf(stdout, "  distorted (mm):   ( %10.6f %10.6f )\n", pd.c[0], pd.c[1]); }
   r2_t pi = tf_sensor_coords_to_image_coords (cpar, pd);
   if (debug) { fprintf(stdout, "  image (pix):      ( %10.6f %10.6f )\n", pi.c[0], pi.c[1]); }
   return pi;
}

void tf_camera_get_mark_position_and_shape
  ( tf_camera_params_t *cpar, 
    r3_t p_w, 
    r3_t u_w, 
    r3_t v_w, 
    r2_t *p_i, 
    r2_t *u_i, 
    r2_t *v_i )
{
    r3_t pu_w, pv_w;

    *p_i  = tf_world_coords_to_image_coords (cpar, p_w);

    r3_add (&p_w, &u_w, &pu_w);
    r2_t pu_i = tf_world_coords_to_image_coords (cpar, pu_w);
    r2_sub (&pu_i, p_i, u_i);

    r3_add (&p_w, &v_w, &pv_w);
    r2_t pv_i = tf_world_coords_to_image_coords (cpar, pv_w);
    r2_sub (&pv_i, p_i, v_i);

}

void tf_camera_matrix_split(r4x4_t *S, r3x3_t *R, r3_t *o)
{
  r4x4_t iS;
  r4x4_inv(S, &iS);

  /* Get camera position: */
  o->c[0] = iS.c[1][0];
  o->c[1] = iS.c[2][0];
  o->c[2] = iS.c[3][0];
  
  /* Get camera rotation: */
  R->c[0][0] = iS.c[1][1]; R->c[0][1] = iS.c[1][2]; R->c[0][2] = iS.c[1][3];
  R->c[1][0] = iS.c[2][1]; R->c[1][1] = iS.c[2][2]; R->c[1][2] = iS.c[2][3];
  R->c[2][0] = iS.c[3][1]; R->c[2][1] = iS.c[3][2]; R->c[2][2] = iS.c[3][3];
}

void tf_camera_matrix_assemble(r3x3_t *R, r3_t *o, r4x4_t *S)
{ r4x4_t iS;
  
  /* Set row 0: */
  iS.c[0][0] = 1.0; iS.c[0][1] = iS.c[0][2] = iS.c[0][3] = 0.0;

  /* Store camera position: */
  iS.c[1][0] = o->c[0];
  iS.c[2][0] = o->c[1];
  iS.c[3][0] = o->c[2];
  
  /* Store camera rotation: */
  iS.c[1][1] = R->c[0][0]; iS.c[1][2] = R->c[0][1]; iS.c[1][3] = R->c[0][2];
  iS.c[2][1] = R->c[1][0]; iS.c[2][2] = R->c[1][1]; iS.c[2][3] = R->c[1][2];
  iS.c[3][1] = R->c[2][0]; iS.c[3][2] = R->c[2][1]; iS.c[3][3] = R->c[2][2];
  
  r4x4_inv(&iS, S);
}

void tf_camera_params_free (tf_camera_params_t *cpar)
{
    free(cpar);
}

void tf_camera_params_write_povray (FILE *wr, tf_camera_params_t *cpar, bool_t flip)
{
  /* For documentation: */
  tf_camera_params_print (cpar, wr);

  /* Compute the camera-to-world transformation matrix {T}: */
  r4x4_t S = cpar->S;
  r4x4_t T;
  r4x4_inv(&S, &T);
  
  /* Compute the POV-Ray angle-of-view parameter: */
  double dtan = (cpar->Npx/2.0)*cpar->dpx/cpar->f; /* Tangent of horizontal half-angle. */
  double angle = 2.0 * atan(dtan) * 180.0/M_PI;

  fprintf(wr, "\n");
  fprintf(wr, "camera {\n");
  fprintf(wr, "  right  %s sqrt(image_width/image_height)*x\n", (flip?"+":"-"));
  fprintf(wr, "  up     sqrt(image_height/image_width)*y\n");
  fprintf(wr, "  angle  %7.3f\n", angle);
  fprintf(wr, "  matrix\n");
  fprintf(wr, "    < %+18.15f, %+18.15f, %+18.15f,\n", -T.c[1][1], -T.c[2][1], -T.c[3][1]);
  fprintf(wr, "      %+18.15f, %+18.15f, %+18.15f,\n", -T.c[1][2], -T.c[2][2], -T.c[3][2]);
  fprintf(wr, "      %+18.15f, %+18.15f, %+18.15f,\n", T.c[1][3], T.c[2][3], T.c[3][3]);
  fprintf(wr, "      %8.2f, %8.2f, %8.2f\n", T.c[1][0], T.c[2][0], T.c[3][0]);
  fprintf(wr, "    >\n");
  fprintf(wr, "  translate ctr\n");
  fprintf(wr, "}\n");

  // /* Compute the camera position in world coordinates: */
  // r4_t org = (r4_t){{ 1.0, 0.0, 0.0, 0.0 }};  /* Camera coords of camera. */
  // r4_t obs;
  // r4x4_map_col(&T, &org, &obs);
  // assert(obs.c[0] == 1.0);
  // /* Compute the look-at point (a point 1m ahead of the camera, on the optical axis): */
  // r4_t cmz = (r4_t){{ 1.0, 0.0, 0.0, 1000.0 }};  /* Camera coords of look-at point. */
  // r4_t see;
  // r4x4_map_col(&T, &cmz, &see); 
  // assert(see.c[0] == 1.0);
  // /* Compute the top-of-camera point (a point 1m above the camera): */
  // r4_t cmy = (r4_t){{ 1.0, 0.0, 1000.0, 0.0 }};  /* Camera coords of top-of-camera point. */
  // r4_t top;
  // r4x4_map_col(&T, &cmy, &top); 
  // assert(top.c[0] == 1.0);

  // fprintf(wr, "  location < %7.1f, %7.1f, %7.1f >\n", obs.c[1], obs.c[2], obs.c[3]);
  // fprintf(wr, "  look_at  < %7.1f, %7.1f, %7.1f >\n", see.c[1], see.c[2], see.c[3]);

  fflush(wr);
}

bool_t tf_camera_point_is_inside_image(r2_t p_i, tf_camera_params_t *cpar)
{
  return 
    (p_i.c[0] >= 0) && (p_i.c[0] <= cpar->Npx) &&
    (p_i.c[1] >= 0) && (p_i.c[1] <= cpar->Npy);
}

void tf_camera_matrix_inverse_from_v_w_and_R (r3_t v_w, r4x4_t *S)
{
  S->c[1][0] = -(S->c[1][1]*v_w.c[0] + S->c[1][2]*v_w.c[1] + S->c[1][3]*v_w.c[2]);
  S->c[2][0] = -(S->c[2][1]*v_w.c[0] + S->c[2][2]*v_w.c[1] + S->c[2][3]*v_w.c[2]);
  S->c[3][0] = -(S->c[3][1]*v_w.c[0] + S->c[3][2]*v_w.c[1] + S->c[3][3]*v_w.c[2]);
}

r3_t tf_camera_matrix_to_v_w (r4x4_t *S)
{
  r3_t  v_w;
  v_w.c[0] =  -(S->c[1][1]*S->c[1][0] + S->c[2][1]*S->c[2][0] + S->c[3][1]*S->c[3][0]);
  v_w.c[1] =  -(S->c[1][2]*S->c[1][0] + S->c[2][2]*S->c[2][0] + S->c[3][2]*S->c[3][0]);
  v_w.c[2] =  -(S->c[1][3]*S->c[1][0] + S->c[2][3]*S->c[2][0] + S->c[3][3]*S->c[3][0]);
  return v_w;
}

void tf_camera_extrapolate 
  ( r4x4_t *S_2, double f_2, double kappa_2, 
    r4x4_t *S_1, double f_1, double kappa_1, 
    r4x4_t *S_0, double *f_0, double *kappa_0
  )
{
  auto void print_params(r4x4_t *S, double f, double kappa, char *title);

  void print_params(r4x4_t *S, double f, double kappa, char *title) {
    fprintf(stderr, "// ---------------------------\n");
    fprintf(stderr, "// %s\n", title);
    fprintf(stderr, "// S =\n");
    int32_t i, j;
    for (i = 0; i < 4; i++) {
      fprintf(stderr, "//   ");
      for (j = 0; j < 4; j++) {
            fprintf(stderr, "%+10.15f \t", S->c[i][j]);
        }
      fprintf(stderr, "\n");
    }
    fprintf(stderr, "// f     = %15.15f\n", f);
    fprintf(stderr, "// kappa = %15.15f\n", kappa);
    fprintf(stderr, "// ---------------------------\n");
 
  }

  print_params(S_2, f_2, kappa_2, "Camera parameters for frame {iframe-2}");
  print_params(S_1, f_1, kappa_1, "Camera parameters for frame {iframe-1}");
 
  /* Split the camera matrices {S_1,S_2} into {R_1,o_1}, {R_2, o_2}: */
  r3x3_t R_1, R_2;
  r3_t o_1, o_2;
  tf_camera_matrix_split(S_1, &R_1, &o_1);
  tf_camera_matrix_split(S_2, &R_2, &o_2);
       
  /* Compute camera translation {Do} in world coords: */
  r3_t Do;
  r3_sub(&o_1, &o_2, &Do);

  fprintf(stderr, "Do: Tx = %f, Ty = %f, Tf = %f\n", Do.c[0], Do.c[1], Do.c[2]);

  /* Compute the incremental camera rotation matrix {DR}: */
  r3x3_t I;
  r3x3_inv (&R_2, &I);
  r3x3_t DR; /*Change in camera matrix {inv(R_2)*R_1}*/
  r3x3_mul(&I, &R_1, &DR);
       
  /* Compute the changes {Df,Dkappa} in focal distance and radial distortion: */
  double Df = f_1/f_2; /* Change in focal distance. */
  double Dkappa = kappa_1 - kappa_2; /* Change in radial distortion. */
       
  fprintf(stderr, "---------------------------\n");
  fprintf(stderr, "changes between the two frames\n");
  tf_camera_params_print_changes (&Do, &DR, Df, Dkappa, stderr);
  fprintf(stderr, "---------------------------\n");
 
  /*avoid excessive changes:*/
  double max_DR = 0.2; /*max camera rotation (radians) */
  double max_Do = 100; /*max camera translation (mm) */
  double max_Df = 1.05; /*max relative change*/
  double max_Dkappa = 0.00001; /*max absolute change*/
  if (sqrt(r3x3_mod_norm_sqr(&DR)) > max_DR) { 
    r3x3_ident (&DR); /* !!! not good */
    fprintf(stderr, "[DR set to ident!]\n");
  }          
  if (r3_norm(&Do) > max_Do) { 
    r3_scale (max_Do/(r3_norm(&Do) + 1.0e-100), &Do, &Do); 
    fprintf(stderr, "[Do clipped to max_Do!]\n");
  }          
  if (Df > max_Df) Df = max_Df;
  if (Df < 1/max_Df) Df = 1/max_Df;
  if (Dkappa < -max_Dkappa) Dkappa = -max_Dkappa;
  if (Dkappa > +max_Dkappa) Dkappa = +max_Dkappa;
       
  fprintf(stderr, "---------------------------\n");
  fprintf(stderr, "max-clipped changes\n");
  tf_camera_params_print_changes (&Do, &DR, Df, Dkappa, stderr);
  fprintf(stderr, "---------------------------\n");

  /* Apply changes to last camera parameters to get the extrapolated guess: */       
  r3x3_t R_0;   /* Extrapolated camera rotation. */
  r3x3_mul(&R_1, &DR, &R_0);
  r3_t o_0;  /* Extrapolated camera position. */
  r3_add (&o_1, &Do, &o_0);
  fprintf(stderr, "o_0: Tx = %f, Ty = %f, Tf = %f\n", o_0.c[0], o_0.c[1], o_0.c[2]);
  tf_camera_matrix_assemble(&R_0, &o_0, S_0);
  (*f_0) = f_1 * Df;
  (*kappa_0) = kappa_1 + Dkappa;

  /* Avoid excessive {kappa} values: */
  double min_kappa = 0.0000;
  double max_kappa = 0.0002;
  if ((*kappa_0) < min_kappa) { (*kappa_0) = min_kappa; }
  if ((*kappa_0) > max_kappa) { (*kappa_0) = max_kappa; }

  print_params(S_0, *f_0, *kappa_0, "Camera parameters extrapolated for frame {iframe}");

}



void tf_camera_matrix_to_euler_angles (r4x4_t *S, r3_t *R)
{
    R->c[2] = atan2 (S->c[2][1], S->c[1][1]);
    double sg = sin(R->c[2]);
    double cg = cos(R->c[2]);
    R->c[1] = atan2 (-S->c[3][1], S->c[1][1] * cg + S->c[2][1] * sg);
    R->c[0] = atan2 (S->c[1][3] * sg - S->c[2][3] * cg, S->c[2][2] * cg - S->c[1][2] * sg);
}


void tf_camera_matrix_from_euler_angles ( r3_t *R, r4x4_t *S )
  {
    double sa = sin(R->c[0]);
    double ca = cos(R->c[0]);
    
    double sb = sin(R->c[1]);
    double cb = cos(R->c[1]);
    
    double sg = sin(R->c[2]);
    double cg = cos(R->c[2]);
    
    S->c[1][1] = cb * cg;                  /* r1 */ 
    S->c[1][2] = cg * sa * sb - ca * sg;   /* r2 */ 
    S->c[1][3] = sa * sg + ca * cg * sb;   /* r3 */ 
    S->c[2][1] = cb * sg;		   /* r4 */ 
    S->c[2][2] = sa * sb * sg + ca * cg;   /* r5 */ 
    S->c[2][3] = ca * sb * sg - cg * sa;   /* r6 */ 
    S->c[3][1] = -sb;			   /* r7 */ 
    S->c[3][2] = cb * sa;		   /* r8 */ 
    S->c[3][3] = ca * cb;		   /* r9 */ 
  }

double tf_camera_adjust_angle (double ang, double aref)
{
  if (isnan(aref)) { return ang; }
  int32_t i = (int32_t)floor((ang-aref)/(2*M_PI) + 0.5);
  ang = ang - i*2*M_PI;
  if (ang < (aref - M_PI)) { ang += 2*M_PI; }
  if (ang > (aref + M_PI)) { ang -= 2*M_PI; }
  return ang;
}

interval_t tf_camera_adjust_angle_range (interval_t ang, double aref)
{
  if (isnan(aref)) { return ang; }
  if ((! isfinite(LO(ang))) || (! isfinite(HI(ang)))) { return ang; }
  /* Convert {ang} to {amid +/- arad} format: */
  double amid = interval_mid(&ang);
  double arad = interval_rad(&ang);
  int32_t i = (int32_t)floor((amid-aref)/(2*M_PI) + 0.5);
  /* Adjust {amid} so that it lies as close as possible to {aref}: */
  amid = amid - i*2*M_PI;
  if (amid < (aref - M_PI)) { amid += 2*M_PI; }
  if (amid > (aref + M_PI)) { amid -= 2*M_PI; }
  /* Return the adjusted interval: */
  return (interval_t){{ amid-arad, amid+arad }};
}

#define MAX_X (1.0e+100)
#define MIN_X (1.0e-100)
#define MAX_LOG_X (+230.25850929940456840179)
#define MIN_LOG_X (-230.25850929940456840179)

double tf_camera_safe_log (double x)
{
  if (x < MIN_X) { return MIN_LOG_X;}
  if (x > MAX_X) { return MAX_LOG_X;}
  return log(x);
}

double tf_camera_safe_exp (double e)
{
  if (e < MIN_LOG_X) {e = MIN_LOG_X;}
  if (e > MAX_LOG_X) {e = MAX_LOG_X;}
  return exp(e);
}

interval_t tf_camera_interval_safe_log (interval_t *xr)
{
  return (interval_t){{ tf_camera_safe_log(LO(*xr)), tf_camera_safe_log(HI(*xr)) }};
}

interval_t tf_camera_interval_safe_exp (interval_t *er)
{
  return (interval_t){{ tf_camera_safe_exp(LO(*er)), tf_camera_safe_exp(HI(*er)) }};
}

double tf_camera_stretch_param (double x, interval_t *xr)
{
  if (x <= LO(*xr)) {return MIN_LOG_X;}
  if (x >= HI(*xr)) {return MAX_LOG_X;}
  return log( (x - LO(*xr) )/( HI(*xr) - x)  );
}

double tf_camera_squeeze_param (double e, interval_t *er)
{
  if (e < MIN_LOG_X) {return LO(*er);}
  if (e > MAX_LOG_X) {return HI(*er);}
  double z = exp(e);
  return (z * HI(*er) + LO(*er))/(1 + z);
}

double tf_camera_params_get_value_from_index (tf_camera_params_t *cpar, int32_t iparam, double vref)
{
  r3_t R;
  if ((iparam == 0) || (iparam == 1) || (iparam == 2))
    { tf_camera_matrix_to_euler_angles (&(cpar->S), &R); }
  r3_t v_w;
  if ((iparam == 3) || (iparam == 4) || (iparam == 5))
    {  v_w = tf_camera_matrix_to_v_w (&(cpar->S)); }
  switch (iparam) {
  case 0:
    return tf_camera_adjust_angle(R.c[0], vref);
  case 1:
    return tf_camera_adjust_angle(R.c[1], vref);
  case 2:
    return tf_camera_adjust_angle(R.c[2], vref);
  case 3:
    return v_w.c[0];
  case 4:
    return v_w.c[1];
  case 5:
    return v_w.c[2];
  case 6:
    return tf_camera_safe_log(cpar->f);
  case 7:
    return cpar->kappa;
  case 8:
    return cpar->sx;
  case 9:
    return cpar->Cx;
  case 10:
    return cpar->Cy;
  case 11:
    return cpar->dpx;
  case 12:
    return cpar->dpy;
  case 13:
    return cpar->Npx;
  case 14:
    return cpar->Npy;
  default:
    fprintf(stderr, "error: the camera parameter doesn't exist\n");
    assert(FALSE);
  }
}

char * tf_camera_params_get_name_from_index (int32_t iparam)
{
  switch (iparam) {
  case 0:
    return "Rx";
  case 1:
    return "Ry";
  case 2:
    return "Rz";
  case 3:
    return "Vx";
  case 4:
    return "Vy";
  case 5:
    return "Vz";
  case 6:
    return "logf";
  case 7:
    return "kappa";
  case 8:
    return "sx";
  case 9:
    return "Cx";
  case 10:
    return "Cy";
  case 11:
    return "dpx";
  case 12: 
    return "dpy";
  case 13: 
    return "Npx";
  case 14: 
    return "Npy";
  default:
    fprintf(stderr, "error: the camera parameter doesn't exist\n");
    assert(FALSE);
  }
}

void tf_camera_params_set_from_vector (tf_camera_params_t *cpar, double params[])
{
   r3_t R;
   R.c[0] = params[0];
   R.c[1] = params[1];
   R.c[2] = params[2];
   tf_camera_matrix_from_euler_angles(&R, &(cpar->S));
   r3_t v_w;
   v_w.c[0] = params[3];
   v_w.c[1] = params[4];
   v_w.c[2] = params[5];
   tf_camera_matrix_inverse_from_v_w_and_R (v_w, &(cpar->S));
   cpar->f = tf_camera_safe_exp(params[6]);
   cpar->kappa = params[7];
   cpar->sx = params[8];
   cpar->Cx = params[9];
   cpar->Cy = params[10];
   cpar->dpx = params[11];
   cpar->dpy = params[12];
   cpar->Npx = params[13];
   cpar->Npx = params[14];
}

r2_t tf_camera_compute_image_error(tf_camera_params_t *cpar, r3_t p_w, r2_t p_i)
{
  r2_t p_i_pre = tf_world_coords_to_image_coords (cpar, p_w);
  double error_X = p_i.c[0] - p_i_pre.c[0];
  double error_Y = p_i.c[1] - p_i_pre.c[1];
  return (r2_t){{ error_X, error_Y }};
}

void tf_camera_compute_all_image_errors
  ( tf_camera_params_t *cpar,
    int32_t nmarks,
    r3_t p_w[],
    r2_t p_i[],
    r2_t e_i[] )
{
  int32_t i;
  for (i = 0; i < nmarks; i++) {
    e_i[i] = tf_camera_compute_image_error(cpar, p_w[i], p_i[i]);
  }
}

double tf_camera_compute_world_error_sqr(tf_camera_params_t *cpar, r3_t p_w, r2_t p_i)
{
  /* Determine the position of the 3D object space point in camera coordinates */
  r3_t p_c = tf_world_coords_to_camera_coords (cpar, p_w);

  /* Convert the measured 2D image coordinates into distorted sensor coordinates */
  r2_t p_d = tf_image_coords_to_sensor_coords (cpar, p_i); 

  /* Convert from distorted sensor coordinates into undistorted sensor plane coordinates */
  r2_t p_u = tf_dis_sensor_coords_to_und_sensor_coords (cpar, p_d);

  /* Find the point of closest approach */
  /* between the undistorted line of sight and the point in 3 space */
  double num = p_c.c[0]*p_u.c[0] + p_c.c[1]*p_u.c[1] + p_c.c[2]*cpar->f;
  double den = p_u.c[0]*p_u.c[0] + p_u.c[1]*p_u.c[1] + cpar->f*cpar->f;
  double t = num/den ;
  r3_t q_c = (r3_t){{ t*p_u.c[0], t*p_u.c[1], t*cpar->f }};

  /* Return the difference between {p_w} and that point: */
  double squared_error = r3_dist_sqr(&p_c, &q_c);
  return squared_error;
}

