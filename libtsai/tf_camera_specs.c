/* See {tf_camera_specs.h}. */
/* Last edited on 2022-10-20 05:52:13 by stolfi */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>

#include <jsfile.h>
#include <interval.h>
#include <affirm.h>
#include <r3.h>
#include <fget.h>
#include <r2.h>
#include <r4x4.h>

#include <tf_camera.h>
#include <tf_camera_specs.h>

tf_camera_specs_t *tf_camera_specs_new (void)
  {
    tf_camera_specs_t *cspec = notnull(malloc(sizeof(tf_camera_specs_t)), "no mem");
    return cspec;
  }

tf_camera_specs_t *tf_camera_specs_copy (tf_camera_specs_t *cspec)
  {
    tf_camera_specs_t *copy = tf_camera_specs_new ();
    copy->Npx =    cspec->Npx;
    copy->Npy =    cspec->Npy;
    copy->dpx =    cspec->dpx;
    copy->dpy =    cspec->dpy;
    copy->Cx =     cspec->Cx;
    copy->Cy =     cspec->Cy;
    copy->sx =     cspec->sx;
    copy->f =      cspec->f;
    copy->kappa =  cspec->kappa;
    copy->v_w[0] = cspec->v_w[0];  
    copy->v_w[1] = cspec->v_w[1];  
    copy->v_w[2] = cspec->v_w[2]; 
    copy->R[0] =   cspec->R[0];   
    copy->R[1] =   cspec->R[1];   
    copy->R[2] =   cspec->R[2];   
    return copy;
  }

void tf_camera_specs_print (FILE *wr, tf_camera_specs_t *cspec)
  {
    fprintf(wr, "// ------------ camera specs ---------------------\n");
    fprintf(wr, "// Npx =    [ %f _ %f ]\n", LO(cspec->Npx), HI(cspec->Npx));
    fprintf(wr, "// Npy =    [ %f _ %f ]\n", LO(cspec->Npy), HI(cspec->Npy));
    fprintf(wr, "// dpx =    [ %.5f _ %.5f ]\n", LO(cspec->dpx), HI(cspec->dpx));
    fprintf(wr, "// dpy =    [ %.5f _ %.5f ]\n", LO(cspec->dpy), HI(cspec->dpy));
    fprintf(wr, "// Cx =     [ %.3f _ %.3f ]\n", LO(cspec->Cx), HI(cspec->Cx));
    fprintf(wr, "// Cy =     [ %.3f _ %.3f ]\n", LO(cspec->Cy), HI(cspec->Cy));
    fprintf(wr, "// sx =     [ %.6f _ %.6f ]\n", LO(cspec->sx), HI(cspec->sx));
    fprintf(wr, "// f =      [ %17.15f _ %17.15f ]\n", LO(cspec->f), HI(cspec->f));
    fprintf(wr, "// kappa =  [ %+17.15f _ %+17.15f ]\n", LO(cspec->kappa), HI(cspec->kappa));
    int32_t i;
    for (i = 0; i < 3; i++) {
      fprintf(wr, "// v_w[%d] = [ %+10.3f _ %+10.3f ]\n", i, LO(cspec->v_w[i]), HI(cspec->v_w[i]));
    }
    for (i = 0; i < 3; i++) {
      fprintf(wr, "// R[%d] =   [ %+10.7f _ %+10.7f ]\n", i, LO(cspec->R[i]), HI(cspec->R[i]));
    }
  
    fprintf(wr, "// -----------------------------------------------------------\n");
  }

void tf_camera_specs_write (FILE *wr, tf_camera_specs_t *cspec)
  {
    fprintf(wr, "Npx         %f \t\t %f\n", LO(cspec->Npx), HI(cspec->Npx));
    fprintf(wr, "Npy         %f \t\t %f\n", LO(cspec->Npy), HI(cspec->Npy));
    fprintf(wr, "dpx         %.5f \t\t %.5f\n", LO(cspec->dpx), HI(cspec->dpx));
    fprintf(wr, "dpy         %.5f \t\t %.5f\n", LO(cspec->dpy), HI(cspec->dpy));
    fprintf(wr, "Cx          %.3f \t\t %.3f\n", LO(cspec->Cx), HI(cspec->Cx));
    fprintf(wr, "Cy          %.3f \t\t %.3f\n", LO(cspec->Cy), HI(cspec->Cy));
    fprintf(wr, "sx          %.6f \t\t %.6f\n", LO(cspec->sx), HI(cspec->sx));
    fprintf(wr, "f           %17.15f \t\t %17.15f\n", LO(cspec->f), HI(cspec->f));
    fprintf(wr, "kappa       %+17.15f \t\t %+17.15f\n", LO(cspec->kappa), HI(cspec->kappa));
    fprintf(wr, "v_w_x       %+10.3f \t\t %+10.3f\n", LO(cspec->v_w[0]), HI(cspec->v_w[0]));
    fprintf(wr, "v_w_y       %+10.3f \t\t %+10.3f\n", LO(cspec->v_w[1]), HI(cspec->v_w[1]));
    fprintf(wr, "v_w_z       %+10.3f \t\t %+10.3f\n", LO(cspec->v_w[2]), HI(cspec->v_w[2]));
    fprintf(wr, "R_x         %+10.7f \t\t %+10.7f\n", LO(cspec->R[0]), HI(cspec->R[0]));
    fprintf(wr, "R_y         %+10.7f \t\t %+10.7f\n", LO(cspec->R[1]), HI(cspec->R[1]));
    fprintf(wr, "R_z         %+10.7f \t\t %+10.7f\n", LO(cspec->R[2]), HI(cspec->R[2]));
  }

tf_camera_specs_t *tf_camera_specs_read (FILE *rd)
  { 
    tf_camera_specs_t *cspec = tf_camera_specs_new();
    cspec->Npx    = tf_camera_specs_read_range (rd, "Npx");
    cspec->Npy    = tf_camera_specs_read_range (rd, "Npy");
    cspec->dpx    = tf_camera_specs_read_range (rd, "dpx");
    cspec->dpy    = tf_camera_specs_read_range (rd, "dpy");
    cspec->Cx     = tf_camera_specs_read_range (rd, "Cx");
    cspec->Cy     = tf_camera_specs_read_range (rd, "Cy");
    cspec->sx     = tf_camera_specs_read_range (rd, "sx");
    cspec->f      = tf_camera_specs_read_range (rd, "f");
    cspec->kappa  = tf_camera_specs_read_range (rd, "kappa");
    cspec->v_w[0] = tf_camera_specs_read_range (rd, "v_w_x");
    cspec->v_w[1] = tf_camera_specs_read_range (rd, "v_w_y");
    cspec->v_w[2] = tf_camera_specs_read_range (rd, "v_w_z");
    cspec->R[0]   = tf_camera_specs_read_range (rd, "R_x");
    cspec->R[1]   = tf_camera_specs_read_range (rd, "R_y");
    cspec->R[2]   = tf_camera_specs_read_range (rd, "R_z");
    return cspec;
  }

interval_t tf_camera_specs_read_range (FILE *rd, char *name)
  {
    char *str = fget_string(rd);
    demand(strcmp(str, name) == 0, "unexpected parameter name");
    interval_t range = (interval_t){{ fget_double(rd), fget_double(rd) }};
    demand(LO(range) <= HI(range), "empty parameter range");
    free(str);
    return range;
  }

tf_camera_specs_t *tf_camera_specs_for_canon_optura (void)
  {
    tf_camera_specs_t *cspec = tf_camera_specs_new ();
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

tf_camera_specs_t *tf_camera_specs_for_sony_dv40 (void)
  {
    tf_camera_specs_t *cspec = tf_camera_specs_new ();
    cspec->Npx =   (interval_t) {{ 720, 720}};
    cspec->Npy =   (interval_t) {{ 480, 480}};
    cspec->dpx =   (interval_t) {{ 0.04, 0.04}}; // Guessed
    cspec->dpy =   (interval_t) {{ 0.04, 0.04}}; // Guessed
    cspec->Cx =    (interval_t) {{ LO(cspec->Npx) / 2, HI(cspec->Npx) / 2 }};
    cspec->Cy =    (interval_t) {{ LO(cspec->Npy) / 2, HI(cspec->Npy) / 2 }};
    cspec->sx =    (interval_t) {{ 1.00, 1.00 }};
    cspec->f =     (interval_t) {{ 1.0, 10000.0}}; // To check
    cspec->kappa = (interval_t) {{ -0.1, 0.1 }}; // To check
    cspec->v_w[0] =  (interval_t) {{ -1000000.0, +1000000.0 }}; /* 100 meters */
    cspec->v_w[1] =  (interval_t) {{ -1000000.0, +1000000.0 }}; /* 100 meters */
    cspec->v_w[2] =  (interval_t) {{ -1000000.0, +1000000.0 }}; /* 100 meters */
    cspec->R[0]   =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */
    cspec->R[1]   =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */
    cspec->R[2]   =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */

    return cspec;
  }

tf_camera_specs_t *tf_camera_specs_for_povray_svga (void)
  {
    tf_camera_specs_t *cspec = tf_camera_specs_new ();
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

tf_camera_specs_t *tf_camera_specs_for_povray_svga_webcam (void)
  {
    tf_camera_specs_t *cspec = tf_camera_specs_new ();
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

tf_camera_specs_t *tf_camera_specs_for_povray_hvga (void)
  {
    tf_camera_specs_t *cspec = tf_camera_specs_new ();
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

tf_camera_specs_t *tf_camera_specs_for_povray_hvga_distorted (void)
  {
    tf_camera_specs_t *cspec = tf_camera_specs_new ();
    cspec->Npx =   (interval_t) {{ 320, 320}};
    cspec->Npy =   (interval_t) {{ 240, 240}};
    cspec->dpx =   (interval_t) {{ 0.04, 0.04}}; // Guessed
    cspec->dpy =   (interval_t) {{ 0.04, 0.04}}; // Guessed
    cspec->Cx =    (interval_t) {{ LO(cspec->Npx) / 2, HI(cspec->Npx) / 2 }};
    cspec->Cy =    (interval_t) {{ LO(cspec->Npy) / 2, HI(cspec->Npy) / 2 }};
    cspec->sx =    (interval_t) {{ 1.00, 1.00 }};
    cspec->f =     (interval_t) {{ 10.0, 10000.0}};
    cspec->kappa = (interval_t) {{ -0.01, +0.01 }}; 
    cspec->v_w[0] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 10 meters */
    cspec->v_w[1] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 10 meters */
    cspec->v_w[2] =  (interval_t) {{ -100000.0, +100000.0 }}; /* 10 meters */
    cspec->R[0] =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */
    cspec->R[1] =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */
    cspec->R[2] =  (interval_t) {{ -INF, +INF }}; /* Unconstrained */

    return cspec;
  }

tf_camera_params_t *tf_camera_specs_get_new_mean_params(tf_camera_specs_t *cspec)
  { 
    tf_camera_params_t *cpar = tf_camera_params_new ();
    tf_camera_specs_get_mean_params(cspec, cpar);
    return cpar;
  }

void tf_camera_specs_get_mean_params(tf_camera_specs_t *cspec, tf_camera_params_t *cpar)
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
      tf_camera_matrix_from_euler_angles(&R, &(cpar->S));
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
    
      tf_camera_matrix_inverse_from_v_w_and_R (v_w, &(cpar->S));
    }
  }

interval_t tf_camera_specs_get_param_range (tf_camera_specs_t *cspec, int32_t iparam, double vref)
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
      return tf_camera_interval_safe_log(&(cspec->f));
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

void tf_camera_specs_set_param_range (tf_camera_specs_t *cspec, int32_t iparam, interval_t *range)
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
      cspec->f = tf_camera_interval_safe_exp(range); break;
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

