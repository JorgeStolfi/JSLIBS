#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include <f2c.h>
#include <lmdif.h>
#include <r3.h>
#include <r2.h>
#include <r4x4.h>
#include <rmxn.h>
#include <string.h>
#include <affirm.h>
#include <tf_camera.h>
#include <tf_math.h>

camera_parameters_t tf_new_camera_parameters (void)
{
  return (camera_parameters_t)notnull(malloc(sizeof(struct _camera_parameters_t)), "no mem");
}

camera_parameters_t tf_copy_camera_parameters (camera_parameters_t cpar)
{
  camera_parameters_t cpar_new = tf_new_camera_parameters();
  (*cpar_new) = (*cpar);
  return cpar_new;
}

void tf_show_camera_parameters (camera_parameters_t cpar, FILE *ferr)
{
    fprintf(ferr, "// ------------ camera parameters --------------------\n");
    fprintf(ferr, "// Npx = %f  Npy = %f\n", cpar->Npx, cpar->Npy);
    fprintf(ferr, "// dpx = %.5f  dpy = %.5f\n", cpar->dpx, cpar->dpy);
    fprintf(ferr, "// Cx = %.3f  Cy = %.3f\n", cpar->Cx, cpar->Cy);
    fprintf(ferr, "// sx = %.6f\n", cpar->sx);
    fprintf(ferr, "// Transformatiom Matrix =\n");
    int32_t i, j;
    for (i = 0; i < 4; i++) {
      fprintf(ferr, "//   ");
      for (j = 0; j < 4; j++) {
            fprintf(ferr, "%+10.15f \t", cpar->S.c[i][j]);
        }
      fprintf(ferr, "\n");
    }
    fprintf(ferr, "// f     = %15.15f\n", cpar->f);
    fprintf(ferr, "// kappa = %15.15f\n", cpar->kappa);
    fprintf(ferr, "// -----------------------------------------------------------\n");
}

void tf_show_rotation_matrix (r4x4_t *S, FILE *ferr)
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
void tf_write_camera_parameters (FILE *wr, int32_t index, camera_parameters_t cpar)
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

int32_t tf_read_camera_parameters (FILE *rd, camera_parameters_t cpar)
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

void tf_write_mutable_camera_parameters (FILE *wr, int32_t index, camera_parameters_t cpar)
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

int32_t tf_read_mutable_camera_parameters (FILE *rd, camera_parameters_t cpar)
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

void tf_print_changes_in_camera_parameters (r3_t *Do, r3x3_t *DR, double Df, double Dkappa, FILE *ferr)
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

r3_t tf_world_coords_to_camera_coords (camera_parameters_t cpar, r3_t pw)
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

r2_t tf_camera_coords_to_undistorted_sensor_coords (camera_parameters_t cpar, r3_t pc)
{
    r2_t pu; 
    pu.c[0] = cpar->f * pc.c[0] / pc.c[2];
    pu.c[1] = cpar->f * pc.c[1] / pc.c[2];
    return pu; 
}

r2_t tf_undistorted_sensor_coords_to_distorted_sensor_coords (camera_parameters_t cpar, r2_t pu)
{
    double Xu = pu.c[0], Yu = pu.c[1];
    double Ru, Rd, lambda, c, d;
    double Q, R, D, S, T, sinT, cosT;

    if (((Xu == 0) && (Yu == 0)) || (cpar->kappa == 0)) {
        return pu;
    }
    
    Ru = hypot (Xu, Yu);        /* sqrt(Xu*Xu+Yu*Yu) */

    c = 1 / cpar->kappa;
    d = -c * Ru;
    
    Q = c / 3;
    R = -d / 2;
    D = CUB (Q) + SQR (R);

    if (D >= 0) {               /* one real root */
        D = sqrt (D);
        S = cbrt (R + D);
        T = cbrt (R - D);
        Rd = S + T;

        if (Rd < 0) {
            Rd = sqrt (-1 / (3 * cpar->kappa));
            fprintf (stderr, "\nWarning: undistorted image point to distorted image point mapping limited by\n");
            fprintf (stderr, "         maximum barrel distortion radius of %lf\n", Rd);
            fprintf (stderr, "         (Xu = %lf, Yu = %lf) -> (Xd = %lf, Yd = %lf)\n\n",
                     Xu, Yu, Xu * Rd / Ru, Yu * Rd / Ru);
        }
        } else {                    /* three real roots */
        D = sqrt (-D);
        S = cbrt (hypot (R, D));
        T = atan2 (D, R) / 3;
        SINCOS (T, sinT, cosT);

        /* the larger positive root is    2*S*cos(T)                   */
        /* the smaller positive root is   -S*cos(T) + sqrt(3)*S*sin(T) */
        /* the negative root is           -S*cos(T) - sqrt(3)*S*sin(T) */

        Rd = -S * cosT + SQRT3 * S * sinT;      /* use the smaller positive root */
    }

    lambda = Rd / Ru;

    r2_t pd;

    pd.c[0] = Xu * lambda;
    pd.c[1] = Yu * lambda;
    return pd;
}


r2_t tf_sensor_coords_to_image_coords (camera_parameters_t cpar, r2_t pd)
{
    r2_t pi;
    pi.c[0] = pd.c[0] * cpar->sx/cpar->dpx + cpar->Cx;
    pi.c[1] = pd.c[1]/cpar->dpy + cpar->Cy;
    return pi;
}

r2_t tf_image_coords_to_sensor_coords (camera_parameters_t cpar, r2_t pi)
{
    r2_t pd;
    pd.c[0] = cpar->dpx * (pi.c[0] - cpar->Cx) / cpar->sx;
    pd.c[1] = cpar->dpy * (pi.c[1] - cpar->Cy);
    return pd;
}


r2_t tf_distorted_sensor_coords_to_undistorted_sensor_coords (camera_parameters_t cpar, r2_t pd)
{
   r2_t pu;

   double distortion_factor = 1 + cpar->kappa * r2_norm_sqr (&pd);
   pu.c[0] = pd.c[0] * distortion_factor;
   pu.c[1] = pd.c[1] * distortion_factor;
   return pu;
}

r2_t tf_world_coords_to_image_coords (camera_parameters_t cpar, r3_t pw)
{
   bool_t debug = FALSE;
   if (debug) { tf_show_camera_parameters (cpar, stdout); }
   if (debug) { fprintf(stdout, "  world (mm):       ( %10.6f %10.6f %10.6f )\n", pw.c[0], pw.c[1], pw.c[2]); }
   r3_t pc = tf_world_coords_to_camera_coords (cpar, pw);
   if (debug) { fprintf(stdout, "  camera (mm):      ( %10.6f %10.6f %10.6f )\n", pc.c[0], pc.c[1], pc.c[2]); }
   r2_t pu = tf_camera_coords_to_undistorted_sensor_coords (cpar, pc);
   if (debug) { fprintf(stdout, "  undistorted (mm): ( %10.6f %10.6f )\n", pu.c[0], pu.c[1]); }
   r2_t pd = tf_undistorted_sensor_coords_to_distorted_sensor_coords (cpar, pu);   
   if (debug) { fprintf(stdout, "  distorted (mm):   ( %10.6f %10.6f )\n", pd.c[0], pd.c[1]); }
   r2_t pi = tf_sensor_coords_to_image_coords (cpar, pd);
   if (debug) { fprintf(stdout, "  image (pix):      ( %10.6f %10.6f )\n", pi.c[0], pi.c[1]); }
   return pi;
}

void tf_mark_estimate_position_and_shape
  ( camera_parameters_t cpar, 
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

void tf_split_camera_matrix(r4x4_t *S, r3x3_t *R, r3_t *o)
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

void tf_assemble_camera_matrix(r3x3_t *R, r3_t *o, r4x4_t *S)
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

void tf_free_camera_parameters_structure (camera_parameters_t cpar)
{
    free(cpar);
}

void tf_write_pov_ray_camera (FILE *wr, camera_parameters_t cpar, bool_t flip)
{
  /* For documentation: */
  tf_show_camera_parameters (cpar, wr);

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

bool_t tf_camera_point_is_inside_image(r2_t p_i, camera_parameters_t cpar)
{
  return 
    (p_i.c[0] >= 0) && (p_i.c[0] <= cpar->Npx) &&
    (p_i.c[1] >= 0) && (p_i.c[1] <= cpar->Npy);
}


void tf_show_camera_specs (camera_specs_t cspec, FILE *ferr)
{
    fprintf(ferr, "// ------------ camera specs ---------------------\n");
    fprintf(ferr, "// Npx =    [ %f _ %f ]\n", LO(cspec->Npx), HI(cspec->Npx));
    fprintf(ferr, "// Npy =    [ %f _ %f ]\n", LO(cspec->Npy), HI(cspec->Npy));
    fprintf(ferr, "// dpx =    [ %.5f _ %.5f ]\n", LO(cspec->dpx), HI(cspec->dpx));
    fprintf(ferr, "// dpy =    [ %.5f _ %.5f ]\n", LO(cspec->dpy), HI(cspec->dpy));
    fprintf(ferr, "// Cx =     [ %.3f _ %.3f ]\n", LO(cspec->Cx), HI(cspec->Cx));
    fprintf(ferr, "// Cy =     [ %.3f _ %.3f ]\n", LO(cspec->Cy), HI(cspec->Cy));
    fprintf(ferr, "// sx =     [ %.6f _ %.6f ]\n", LO(cspec->sx), HI(cspec->sx));
    fprintf(ferr, "// f =      [ %17.15f _ %17.15f ]\n", LO(cspec->f), HI(cspec->f));
    fprintf(ferr, "// kappa =  [ %+17.15f _ %+17.15f ]\n", LO(cspec->kappa), HI(cspec->kappa));
    int32_t i;
    for (i = 0; i < 3; i++) {
      fprintf(ferr, "// v_w[%d] = [ %+10.3f _ %+10.3f ]\n", i, LO(cspec->v_w[i]), HI(cspec->v_w[i]));
    }
    for (i = 0; i < 3; i++) {
      fprintf(ferr, "// R[%d] =   [ %+10.7f _ %+10.7f ]\n", i, LO(cspec->R[i]), HI(cspec->R[i]));
    }
    
  
    fprintf(ferr, "// -----------------------------------------------------------\n");
}

void tf_write_camera_specs (camera_specs_t cspec, FILE *f)
{
    fprintf(f, "Npx         %f \t\t %f\n", LO(cspec->Npx), HI(cspec->Npx));
    fprintf(f, "Npy         %f \t\t %f\n", LO(cspec->Npy), HI(cspec->Npy));
    fprintf(f, "dpx         %.5f \t\t %.5f\n", LO(cspec->dpx), HI(cspec->dpx));
    fprintf(f, "dpy         %.5f \t\t %.5f\n", LO(cspec->dpy), HI(cspec->dpy));
    fprintf(f, "Cx          %.3f \t\t %.3f\n", LO(cspec->Cx), HI(cspec->Cx));
    fprintf(f, "Cy          %.3f \t\t %.3f\n", LO(cspec->Cy), HI(cspec->Cy));
    fprintf(f, "sx          %.6f \t\t %.6f\n", LO(cspec->sx), HI(cspec->sx));
    fprintf(f, "f           %17.15f \t\t %17.15f\n", LO(cspec->f), HI(cspec->f));
    fprintf(f, "kappa       %+17.15f \t\t %+17.15f\n", LO(cspec->kappa), HI(cspec->kappa));
    fprintf(f, "v_w_x       %+10.3f \t\t %+10.3f\n", LO(cspec->v_w[0]), HI(cspec->v_w[0]));
    fprintf(f, "v_w_y       %+10.3f \t\t %+10.3f\n", LO(cspec->v_w[1]), HI(cspec->v_w[1]));
    fprintf(f, "v_w_z       %+10.3f \t\t %+10.3f\n", LO(cspec->v_w[2]), HI(cspec->v_w[2]));
    fprintf(f, "R_x         %+10.7f \t\t %+10.7f\n", LO(cspec->R[0]), HI(cspec->R[0]));
    fprintf(f, "R_y         %+10.7f \t\t %+10.7f\n", LO(cspec->R[1]), HI(cspec->R[1]));
    fprintf(f, "R_z         %+10.7f \t\t %+10.7f\n", LO(cspec->R[2]), HI(cspec->R[2]));
}

void tf_read_camera_specs (camera_specs_t cspec, FILE *f)
{ 
    tf_read_camera_specs_param (f, "Npx",   &(cspec->Npx.end[0]), &(cspec->Npx.end[1]));
    tf_read_camera_specs_param (f, "Npy",   &(cspec->Npy.end[0]), &(cspec->Npy.end[1]));
    tf_read_camera_specs_param (f, "dpx",   &(cspec->dpx.end[0]), &(cspec->dpx.end[1]));
    tf_read_camera_specs_param (f, "dpy",   &(cspec->dpy.end[0]), &(cspec->dpy.end[1]));
    tf_read_camera_specs_param (f, "Cx",    &(cspec->Cx.end[0]), &(cspec->Cx.end[1]));
    tf_read_camera_specs_param (f, "Cy",    &(cspec->Cy.end[0]), &(cspec->Cy.end[1]));
    tf_read_camera_specs_param (f, "sx",    &(cspec->sx.end[0]), &(cspec->sx.end[1]));
    tf_read_camera_specs_param (f, "f",     &(cspec->f.end[0]), &(cspec->f.end[1]));
    tf_read_camera_specs_param (f, "kappa", &(cspec->kappa.end[0]), &(cspec->kappa.end[1]));
    tf_read_camera_specs_param (f, "v_w_x", &(cspec->v_w[0].end[0]), &(cspec->v_w[0].end[1]));
    tf_read_camera_specs_param (f, "v_w_y", &(cspec->v_w[1].end[0]), &(cspec->v_w[1].end[1]));
    tf_read_camera_specs_param (f, "v_w_z", &(cspec->v_w[2].end[0]), &(cspec->v_w[2].end[1]));
    tf_read_camera_specs_param (f, "R_x",   &(cspec->R[0].end[0]), &(cspec->R[0].end[1]));
    tf_read_camera_specs_param (f, "R_y",   &(cspec->R[1].end[0]), &(cspec->R[1].end[1]));
    tf_read_camera_specs_param (f, "R_z",   &(cspec->R[2].end[0]), &(cspec->R[2].end[1]));
}

void tf_read_camera_specs_param (FILE *f, char *param, double *end0, double *end1)
{
    char *str = (char *)malloc(512*sizeof(char *)); 
    fscanf(f, "%s", str);
    if (strcmp(str, param) == 0) {
      fscanf(f, "%lf", end0);
      fscanf(f, "%lf", end1);
    }	
    else {
        fprintf(stderr, "error: missing the %s camera specs parameter\n", param);
        assert(FALSE);
    }
    free(str);
}


camera_specs_t tf_camera_copy_specs (camera_specs_t cspec)
{
    camera_specs_t copy_cspec = (camera_specs_t)malloc(sizeof(struct _camera_specs_t));
    copy_cspec->Npx =    cspec->Npx;
    copy_cspec->Npy =    cspec->Npy;
    copy_cspec->dpx =    cspec->dpx;
    copy_cspec->dpy =    cspec->dpy;
    copy_cspec->Cx =     cspec->Cx;
    copy_cspec->Cy =     cspec->Cy;
    copy_cspec->sx =     cspec->sx;
    copy_cspec->f =      cspec->f;
    copy_cspec->kappa =  cspec->kappa;
    copy_cspec->v_w[0] = cspec->v_w[0];  
    copy_cspec->v_w[1] = cspec->v_w[1];  
    copy_cspec->v_w[2] = cspec->v_w[2]; 
    copy_cspec->R[0] =   cspec->R[0];   
    copy_cspec->R[1] =   cspec->R[1];   
    copy_cspec->R[2] =   cspec->R[2];   
    return copy_cspec;
}

camera_specs_t tf_get_canon_optura_specs (void)
{
    camera_specs_t cspec = (camera_specs_t)malloc(sizeof(struct _camera_specs_t));
    cspec->Npx =   (interval_t) {{ 320, 320}};
    cspec->Npy =   (interval_t) {{ 240, 240}};
    cspec->dpx =   (interval_t) {{ 0.04, 0.04}}; // Guessed
    cspec->dpy =   (interval_t) {{ 0.04, 0.04}}; // Guessed
    cspec->Cx =    (interval_t) {{ LO(cspec->Npx) / 2, HI(cspec->Npx) / 2 }};
    cspec->Cy =    (interval_t) {{ LO(cspec->Npy) / 2, HI(cspec->Npy) / 2 }};
    cspec->sx =    (interval_t) {{ 1.00, 1.00 }};
    cspec->f =     (interval_t) {{ 1.0, 1000.0}}; // To check
    cspec->kappa = (interval_t) {{ -0.1, 0.1 }}; // To check

    cspec->v_w[0] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 100 meters */
    cspec->v_w[1] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 100 meters */
    cspec->v_w[2] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 100 meters */
    //cspec->v_w[0] =  (interval_t) {{ -INF, +INF }}; /* 10 meters */
    //cspec->v_w[1] =  (interval_t) {{ -INF, +INF }}; /* 10 meters */
    //cspec->v_w[2] =  (interval_t) {{ -INF, +INF }}; /* 10 meters */
    cspec->R[0]   =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */
    cspec->R[1]   =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */
    cspec->R[2]   =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */

    return cspec;
}

camera_specs_t tf_get_sony_dv40_specs (void)
{
    camera_specs_t cspec = (camera_specs_t)malloc(sizeof(struct _camera_specs_t));
    cspec->Npx =   (interval_t) {{ 720, 720}};
    cspec->Npy =   (interval_t) {{ 540, 540}};
    cspec->dpx =   (interval_t) {{ 0.04, 0.04}}; // Guessed
    cspec->dpy =   (interval_t) {{ 0.04, 0.04}}; // Guessed
    cspec->Cx =    (interval_t) {{ LO(cspec->Npx) / 2, HI(cspec->Npx) / 2 }};
    cspec->Cy =    (interval_t) {{ LO(cspec->Npy) / 2, HI(cspec->Npy) / 2 }};
    cspec->sx =    (interval_t) {{ 1.00, 1.00 }};
    cspec->f =     (interval_t) {{ 1.0, 1000.0}}; // To check
    cspec->kappa = (interval_t) {{ -0.1, 0.1 }}; // To check
    cspec->v_w[0] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 100 meters */
    cspec->v_w[1] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 100 meters */
    cspec->v_w[2] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 100 meters */
    cspec->R[0]   =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */
    cspec->R[1]   =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */
    cspec->R[2]   =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */

    return cspec;
}

camera_specs_t tf_get_svga_povray_specs (void)
{
    camera_specs_t cspec = (camera_specs_t)malloc(sizeof(struct _camera_specs_t));
    cspec->Npx =   (interval_t) {{ 640, 640}};
    cspec->Npy =   (interval_t) {{ 480, 480}};
    cspec->dpx =   (interval_t) {{ 0.04, 0.04}}; // Guessed
    cspec->dpy =   (interval_t) {{ 0.04, 0.04}}; // Guessed
    cspec->Cx =    (interval_t) {{ LO(cspec->Npx) / 2, HI(cspec->Npx) / 2 }};
    cspec->Cy =    (interval_t) {{ LO(cspec->Npy) / 2, HI(cspec->Npy) / 2 }};
    cspec->sx =    (interval_t) {{ 1.00, 1.00 }};
    cspec->f =     (interval_t) {{ 0.10, 1000000.0}};
    cspec->kappa = (interval_t) {{ 0.00, 0.00 }}; // No radial distortion

    cspec->v_w[0] =  (interval_t) {{ -1000000.0, +1000000.0 }}; /* 100 meters */
    cspec->v_w[1] =  (interval_t) {{ -1000000.0, +1000000.0 }}; /* 100 meters */
    cspec->v_w[2] =  (interval_t) {{ -1000000.0, +1000000.0 }}; /* 100 meters */
    cspec->R[0] =  (interval_t) {{ -INF, +INF}}; /* Unconstrained */
    cspec->R[1] =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */
    cspec->R[2] =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */

    return cspec;
}

camera_specs_t tf_get_svga_povray_webcam_specs (void)
{
    camera_specs_t cspec = (camera_specs_t)malloc(sizeof(struct _camera_specs_t));
    cspec->Npx =   (interval_t) {{ 640, 640}};
    cspec->Npy =   (interval_t) {{ 480, 480}};
    cspec->dpx =   (interval_t) {{ 0.04, 0.04}}; // Guessed
    cspec->dpy =   (interval_t) {{ 0.04, 0.04}}; // Guessed
    cspec->Cx =    (interval_t) {{ LO(cspec->Npx) / 2, HI(cspec->Npx) / 2 }};
    cspec->Cy =    (interval_t) {{ LO(cspec->Npy) / 2, HI(cspec->Npy) / 2 }};
    cspec->sx =    (interval_t) {{ 1.0, 1.0 }};
    //cspec->f =     (interval_t) {{ 46.50, 46.50}};
    cspec->f =     (interval_t) {{ 100.0, 100.0}};
    cspec->kappa = (interval_t) {{ 0.00, 0.00 }}; 

    cspec->v_w[0] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 10 meters */
    cspec->v_w[1] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 10 meters */
    cspec->v_w[2] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 10 meters */
    cspec->R[0] =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */
    cspec->R[1] =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */
    cspec->R[2] =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */

    return cspec;
}

camera_specs_t tf_get_hvga_povray_specs (void)
{
    camera_specs_t cspec = (camera_specs_t)malloc(sizeof(struct _camera_specs_t));
    cspec->Npx =   (interval_t) {{ 320, 320}};
    cspec->Npy =   (interval_t) {{ 240, 240}};
    cspec->dpx =   (interval_t) {{ 0.04, 0.04}}; // Guessed
    cspec->dpy =   (interval_t) {{ 0.04, 0.04}}; // Guessed
    cspec->Cx =    (interval_t) {{ LO(cspec->Npx) / 2, HI(cspec->Npx) / 2 }};
    cspec->Cy =    (interval_t) {{ LO(cspec->Npy) / 2, HI(cspec->Npy) / 2 }};
    cspec->sx =    (interval_t) {{ 1.00, 1.00 }};
    cspec->f =     (interval_t) {{ 10.0, 10000.0}};
    /*Cuidado mudar*/
    //cspec->kappa = (interval_t) {{ 0.00, 0.00 }}; // No radial distortion
    cspec->kappa = (interval_t) {{ 0.00, 0.00 }}; // No radial distortion

    cspec->v_w[0] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 10 meters */
    cspec->v_w[1] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 10 meters */
    cspec->v_w[2] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 10 meters */
    cspec->R[0] =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */
    cspec->R[1] =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */
    cspec->R[2] =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */

    return cspec;
}

camera_parameters_t tf_get_mean_camera_parameters(camera_specs_t cspec)
{ 
  camera_parameters_t cpar = (camera_parameters_t)malloc(sizeof(struct _camera_parameters_t));
  tf_store_mean_camera_parameters(cspec, cpar);
    return cpar;
}

void tf_store_mean_camera_parameters(camera_specs_t cspec, camera_parameters_t cpar)
{
    cpar->Npx =   interval_mid(&(cspec->Npx));   
    cpar->Npy =   interval_mid(&(cspec->Npy));   
    cpar->dpx =   interval_mid(&(cspec->dpx));   
    cpar->dpy =   interval_mid(&(cspec->dpy));
    cpar->Cx =    interval_mid(&(cspec->Cx));    
    cpar->Cy =    interval_mid(&(cspec->Cy));    
    cpar->sx =    interval_mid(&(cspec->sx));    
    cpar->f =     sqrt(LO(cspec->f)*HI(cspec->f)); /* geometric mean */
    cpar->kappa = interval_mid(&(cspec->kappa)); 

    /* Initialize the first row of {S}: */
    cpar->S.c[0][0] = 1.0;
    cpar->S.c[0][1] = 0.0;
    cpar->S.c[0][2] = 0.0;
    cpar->S.c[0][3] = 0.0;
    if ( interval_is_full (&(cspec->R[0])) && 
         interval_is_full (&(cspec->R[1])) && 
         interval_is_full (&(cspec->R[2])) ) {
      /* Default: camera on world's X axis, looking at origin. */
      cpar->S = (r4x4_t){{
        { 1.0,  00.00,  00.00, 00.00 },
        { 0.0,  00.00,  +1.00, 00.00 },
        { 0.0,  00.00,  00.00, -1.00 },
        { 0.0,  -1.00,  00.00, 00.00 }
      }};
    } else {
      /* Get the rotation matrix from the mean Euler angles: */
      r3_t R;
      R.c[0] = interval_mid(&(cspec->R[0]));
      R.c[1] = interval_mid(&(cspec->R[1]));
      R.c[2] = interval_mid(&(cspec->R[2]));
      tf_set_S_matrix_from_euler_angles(&R, &(cpar->S));
    }
  
    if ( interval_is_full (&(cspec->v_w[0])) && 
         interval_is_full (&(cspec->v_w[1])) &&  
         interval_is_full (&(cspec->v_w[2]))) {
      /* Place the camera 1 m away from the origin, looking into it: */
      cpar->S.c[1][0] = 0.0;
      cpar->S.c[2][0] = 0.0;
      cpar->S.c[3][0] = 1000.0;
    } else {
      /* Compute {Tx,Ty,Tz} from midpoint of {Vx_w,Vy_w,Vz_w}: */
      r3_t v_w;
      v_w.c[0] = interval_mid(&(cspec->v_w[0]));
      v_w.c[1] = interval_mid(&(cspec->v_w[1]));
      v_w.c[2] = interval_mid(&(cspec->v_w[2]));
    
      tf_compute_T_from_v_w_and_R (v_w, &(cpar->S));
    }
}


void tf_compute_T_from_v_w_and_R (r3_t v_w, r4x4_t *S)
{
  S->c[1][0] = -(S->c[1][1]*v_w.c[0] + S->c[1][2]*v_w.c[1] + S->c[1][3]*v_w.c[2]);
  S->c[2][0] = -(S->c[2][1]*v_w.c[0] + S->c[2][2]*v_w.c[1] + S->c[2][3]*v_w.c[2]);
  S->c[3][0] = -(S->c[3][1]*v_w.c[0] + S->c[3][2]*v_w.c[1] + S->c[3][3]*v_w.c[2]);
}

r3_t tf_compute_v_w_from_S (r4x4_t *S)
{
  r3_t  v_w;
  v_w.c[0] =  -(S->c[1][1]*S->c[1][0] + S->c[2][1]*S->c[2][0] + S->c[3][1]*S->c[3][0]);
  v_w.c[1] =  -(S->c[1][2]*S->c[1][0] + S->c[2][2]*S->c[2][0] + S->c[3][2]*S->c[3][0]);
  v_w.c[2] =  -(S->c[1][3]*S->c[1][0] + S->c[2][3]*S->c[2][0] + S->c[3][3]*S->c[3][0]);
  return v_w;
}

void tf_extrapolate_camera_parameters 
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
  tf_split_camera_matrix(S_1, &R_1, &o_1);
  tf_split_camera_matrix(S_2, &R_2, &o_2);
       
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
  tf_print_changes_in_camera_parameters (&Do, &DR, Df, Dkappa, stderr);
  fprintf(stderr, "---------------------------\n");
 
  /*avoid excessive changes:*/
  double max_DR = 0.2; /*max camera rotation (radians) */
  double max_Do = 100; /*max camera translation (mm) */
  double max_Df = 1.05; /*max relative change*/
  double max_Dkappa = 0.00001; /*max absolute change*/
  if (sqrt(r3x3_mag_sqr(&DR)) > max_DR) { 
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
  tf_print_changes_in_camera_parameters (&Do, &DR, Df, Dkappa, stderr);
  fprintf(stderr, "---------------------------\n");

  /* Apply changes to last camera parameters to get the extrapolated guess: */       
  r3x3_t R_0;   /* Extrapolated camera rotation. */
  r3x3_mul(&R_1, &DR, &R_0);
  r3_t o_0;  /* Extrapolated camera position. */
  r3_add (&o_1, &Do, &o_0);
  fprintf(stderr, "o_0: Tx = %f, Ty = %f, Tf = %f\n", o_0.c[0], o_0.c[1], o_0.c[2]);
  tf_assemble_camera_matrix(&R_0, &o_0, S_0);
  (*f_0) = f_1 * Df;
  (*kappa_0) = kappa_1 + Dkappa;

  /* Avoid excessive {kappa} values: */
  double min_kappa = 0.0000;
  double max_kappa = 0.0002;
  if ((*kappa_0) < min_kappa) { (*kappa_0) = min_kappa; }
  if ((*kappa_0) > max_kappa) { (*kappa_0) = max_kappa; }

  print_params(S_0, *f_0, *kappa_0, "Camera parameters extrapolated for frame {iframe}");

}



void tf_get_euler_angles_from_S_matrix (r4x4_t *S, r3_t *R)
{
    double    sg, cg;
    R->c[2] = atan2 (S->c[2][1], S->c[1][1]);
    SINCOS (R->c[2], sg, cg);
    R->c[1] = atan2 (-S->c[3][1], S->c[1][1] * cg + S->c[2][1] * sg);
    R->c[0] = atan2 (S->c[1][3] * sg - S->c[2][3] * cg, S->c[2][2] * cg - S->c[1][2] * sg);
}


void tf_set_S_matrix_from_euler_angles ( r3_t *R, r4x4_t *S )
{
    double sa, ca, sb, cb, sg, cg;
    SINCOS (R->c[0], sa, ca);
    SINCOS (R->c[1], sb, cb);
    SINCOS (R->c[2], sg, cg);
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

double tf_safe_log (double x)
{
  if (x < MIN_X) { return MIN_LOG_X;}
  if (x > MAX_X) { return MAX_LOG_X;}
  return log(x);
}

double tf_safe_exp (double e)
{
  if (e < MIN_LOG_X) {e = MIN_LOG_X;}
  if (e > MAX_LOG_X) {e = MAX_LOG_X;}
  return exp(e);
}

interval_t tf_interval_safe_log (interval_t *xr)
{
  return (interval_t){{ tf_safe_log(LO(*xr)), tf_safe_log(HI(*xr)) }};
}

interval_t tf_interval_safe_exp (interval_t *er)
{
  return (interval_t){{ tf_safe_exp(LO(*er)), tf_safe_exp(HI(*er)) }};
}

double tf_stretch (double x, interval_t *xr)
{
  if (x <= LO(*xr)) {return MIN_LOG_X;}
  if (x >= HI(*xr)) {return MAX_LOG_X;}
  return log( (x - LO(*xr) )/( HI(*xr) - x)  );
}

double tf_squeeze (double e, interval_t *er)
{
  if (e < MIN_LOG_X) {return LO(*er);}
  if (e > MAX_LOG_X) {return HI(*er);}
  double z = exp(e);
  return (z * HI(*er) + LO(*er))/(1 + z);
}

double tf_camera_get_param (camera_parameters_t cpar, int32_t iparam, double vref)
{
  r3_t R;
  if ((iparam == 0) || (iparam == 1) || (iparam == 2))
    { tf_get_euler_angles_from_S_matrix (&(cpar->S), &R); }
  r3_t v_w;
  if ((iparam == 3) || (iparam == 4) || (iparam == 5))
    {  v_w = tf_compute_v_w_from_S (&(cpar->S)); }
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
    return tf_safe_log(cpar->f);
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

interval_t tf_camera_get_param_range (camera_specs_t cspec, int32_t iparam, double vref)
{
  switch (iparam) {
  case 0:
    return tf_camera_adjust_angle_range(cspec->R[0], vref);
  case 1:
    return tf_camera_adjust_angle_range(cspec->R[1], vref);
  case 2:
    return tf_camera_adjust_angle_range(cspec->R[2], vref);
  case 3:
    return cspec->v_w[0];
  case 4:
    return cspec->v_w[1];
  case 5:
    return cspec->v_w[2];
  case 6:
    return tf_interval_safe_log(&(cspec->f));
  case 7:
    return cspec->kappa;
  case 8:
    return cspec->sx;
  case 9:
    return cspec->Cx;
  case 10:
    return cspec->Cy;
  case 11:
    return cspec->dpx;
  case 12:
    return cspec->dpy;
  case 13:
    return cspec->Npx;
  case 14:
    return cspec->Npy;
  default:
    fprintf(stderr, "error: the camera parameter doesn't exist\n");
    assert(FALSE);
  }
}

void tf_camera_set_param_range (camera_specs_t cspec, int32_t iparam, interval_t *range)
{
  switch (iparam) {
  case 0:
    cspec->R[0] = *range; break;
  case 1:
    cspec->R[1] = *range; break;
  case 2:
    cspec->R[2] = *range; break;
  case 3:
    cspec->v_w[0] = *range; break;
  case 4:
    cspec->v_w[1] = *range; break;
  case 5:
    cspec->v_w[2] = *range; break;
  case 6:
    cspec->f = tf_interval_safe_exp(range); break;
  case 7:
    cspec->kappa = *range; break;
  case 8:
    cspec->sx = *range; break;
  case 9:
    cspec->Cx = *range; break;
  case 10:
    cspec->Cy = *range; break;
  case 11:
    cspec->dpx = *range; break;
  case 12:
    cspec->dpy = *range; break;
  case 13:
    cspec->Npx = *range; break;
  case 14:
    cspec->Npy = *range; break;
  default:
    fprintf(stderr, "error: the camera parameter doesn't exist\n");
    assert(FALSE);
  }
}

char * tf_camera_param_name (int32_t iparam)
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

void tf_camera_set_params_from_vector (camera_parameters_t cpar, double params[])
{
   r3_t R;
   R.c[0] = params[0];
   R.c[1] = params[1];
   R.c[2] = params[2];
   tf_set_S_matrix_from_euler_angles(&R, &(cpar->S));
   r3_t v_w;
   v_w.c[0] = params[3];
   v_w.c[1] = params[4];
   v_w.c[2] = params[5];
   tf_compute_T_from_v_w_and_R (v_w, &(cpar->S));
   cpar->f = tf_safe_exp(params[6]);
   cpar->kappa = params[7];
   cpar->sx = params[8];
   cpar->Cx = params[9];
   cpar->Cy = params[10];
   cpar->dpx = params[11];
   cpar->dpy = params[12];
   cpar->Npx = params[13];
   cpar->Npx = params[14];
}

r2_t tf_camera_compute_image_error(camera_parameters_t cpar, r3_t p_w, r2_t p_i)
{
  r2_t p_i_pre = tf_world_coords_to_image_coords (cpar, p_w);
  double error_X = p_i.c[0] - p_i_pre.c[0];
  double error_Y = p_i.c[1] - p_i_pre.c[1];
  return (r2_t){{ error_X, error_Y }};
}

void tf_camera_compute_all_image_errors
  ( camera_parameters_t cpar,
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

double tf_camera_compute_world_error_sqr(camera_parameters_t cpar, r3_t p_w, r2_t p_i)
{
  /* Determine the position of the 3D object space point in camera coordinates */
  r3_t p_c = tf_world_coords_to_camera_coords (cpar, p_w);

  /* Convert the measured 2D image coordinates into distorted sensor coordinates */
  r2_t p_d = tf_image_coords_to_sensor_coords (cpar, p_i); 

  /* Convert from distorted sensor coordinates into undistorted sensor plane coordinates */
  r2_t p_u = tf_distorted_sensor_coords_to_undistorted_sensor_coords (cpar, p_d);

  /* Find the point of closest approach */
  /* between the undistorted line of sight and the point in 3 space */
  double num = p_c.c[0]*p_u.c[0] + p_c.c[1]*p_u.c[1] + p_c.c[2]*cpar->f;
  double den = SQR(p_u.c[0]) + SQR(p_u.c[1]) + SQR(cpar->f);
  double t = num/den ;
  r3_t q_c = (r3_t){{ t*p_u.c[0], t*p_u.c[1], t*cpar->f }};

  /* Return the difference between {p_w} and that point: */
  double squared_error = SQR(p_c.c[0] - q_c.c[0]) + SQR (p_c.c[1] - q_c.c[1]) + SQR (p_c.c[2] - q_c.c[2]);
  return squared_error;
}

