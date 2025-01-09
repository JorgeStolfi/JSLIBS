#ifndef __LIBRAAB_NORMINTERP_H__
#define __LIBRAAB_NORMINTERP__

#include <r3.h>
#include <r2.h>
#include <float_image.h> 

typedef enum{ PROB_BEST,PROB_AVERAGE,POS_POLY, HIGHLIGHT_12PS } normal_interp_t;

float_image_t* normal_interpolate_prob(float_image_t** fim_normals, float_image_t** fim_logprobs, int n, normal_interp_t interp_opt,float_image_t* select_map);
float_image_t* normal_interpolate_pos(float_image_t** fim_normals, r2_t* gauge_positions, int n, int degree);
float_image_t* normal_interpolate_hightlight12(float_image_t** fim_normals,float_image_t** fim_weights, int n,r3_t* cluster_dir,double k);
/*Interpolates averaged normals using a hightlight correcction factor where {k} is the  specular factor. recommended k = 5 */
#endif 