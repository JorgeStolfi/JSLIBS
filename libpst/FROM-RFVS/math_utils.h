#ifndef __MATHUTILS_H__
#define __MATHUTILS_H__

#include <r3.h>
#include <r3x3.h>

r3x3_t compute_normal_correction_matrix(r3_t view_dir);
/*Computes transformation matrix from view dir*/
r3_t compute_bissetriz(r3_t u, r3_t v);
/*Computes the bisectrix of {u} and {v} */
void convert_r3_to_stereographic(r3_t v, double* X, double* Y,r3_t view_dir);
/*Converts from a r3 vector to a 2d coordinates of stereographic projection with center in {view_dir}*/
r3_t find_ortho_vector(r3_t view_dir,r3_t u);
/*Find a vector orthogonal {t} to view_dir and coplanar with u*/
r3_t bend_towards(r3_t v, r3_t u, double ang);
/*Bends a vector {u} towards unit vector {v} by {ang} radians */
void r3_and_error_to_stereographic_box(r3_t u, double err_u, r3_t view_dir, double* X_min, double* X_max, double* Y_min, double* Y_max);
/*Given a r3 {u} computes the box of a stereographic projection with with center in {view_dir}
  where the normals with error {err_u} (in radians) are contained*/
r3_t convert_stereographic_to_r3(double X, double Y, r3_t view_dir);
/*Converts from 2d coordinates of stereographic projection with center in {view_dir} to a r3 vector*/
double gaussian_distribution(double x, double mean, double sigma);
/*Computes value of f(x) of a gaussian distribution with mean {mean} and standard deviation {sigma} */
void clip_direction_to_circle(r3_t* u, r3_t u_old, double err_old);
/*If the angle {u} and ${u_old} is greater than {err_old}, replaces {u} by the nearest vector
  that makes an angle equal to {err_old}.
*/
void clip_scalar_to_interval(double* v, double v_old, double err_old);

#endif