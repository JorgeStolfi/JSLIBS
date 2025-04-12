/* Last edited on 2025-03-19 12:47:27 by stolfi */

#define _GNU_SOURCE
#include <stdlib.h>
#include <math.h>
#include <analytic_alpha.h>
#include <r3x3.h>
#include <tabela.h>
#include <assert.h>
#include <sve_minn.h>
#include <sve_minn_iterate.h>
#include <jsfile.h>
#include <math_utils.h>
#include <estlogprob.h>

#include <nonlinear_optimization.h>


analytic_alpha_log_prob_function_t* escolhe_analytic_alpha_log_prob_function_t(int nLogProb){
    switch(nLogProb) {
    case 0: return NULL;
    case 7: return analytic_alpha_log_prob_noalbedo07;
    case 8: return analytic_alpha_log_prob_noalbedo08;
    case 9: return analytic_alpha_log_prob_noalbedo09;
   default: demand(FALSE, "numero invalido da funcao de probabilidade"); return NULL;
  }
}

double analytic_alpha_log_prob_noalbedo07(double a,double S[],double G[],int m,double sigma, double omg0, double omg1 ){
//    assert(isnan(a));
//    double so[m];
//    double go[m];
//    double Smag = rn_dir(m,S,so);
//    double Gmag = rn_dir(m,G,go);
//    double val = EstLogPrSG_07(so,Smag,go, Gmag,m, sigma, omg0,omg1);
   int i;
   double val = 0;
   for(i = 0; i < m; i++){
     val+=LogProbSiGialb(S[i], G[i],a,sigma,omg0,omg1);
   }
//    double val = EstLogPrSG_07(so,Smag,go, Gmag,m, sigma, omg0,omg1);
//    assert(!isnan(val));
   if(isnan(val)) return 0;
     
   return val;


}


double analytic_alpha_log_prob_noalbedo08(double a,double S[],double G[],int m,double sigma, double omg0, double omg1 ){
   assert(isnan(a));
   double so[m];
   double go[m];
   double Smag = rn_dir(m,S,so);
   double Gmag = rn_dir(m,G,go);
   return EstLogPrSG_08(so,Smag,go, Gmag,m, sigma, omg0,omg1);


}

double analytic_alpha_log_prob_noalbedo09(double a,double S[],double G[],int m,double sigma, double omg0, double omg1 ){
   assert(isnan(a));
   double so[m];
   double go[m];
   double Smag = rn_dir(m,S,so);
   double Gmag = rn_dir(m,G,go);
   return EstLogPrSG_09(so,Smag,go, Gmag,m, sigma, omg0,omg1);


}

double analytic_alpha_log_prob_withalbedo(double a, double S[],double G[],int m,double sigma, double omg0, double omg1 ){
  double val = 0;
   int i;
   assert(!isnan(a));
   for(i = 0; i < m; i++){
     double p_outlier = (S[i] < a*G[i] ? omg0/(a*G[i] + 0.001) : omg1/(1.001 - a*G[i]) ); /*Prob. density due to outliers*/
     double p_inlier  = (1 - omg0 - omg1)*gaussian_distribution(S[i] - a*G[i],0.0,sigma); /*Prob. density due to inliers*/
     val+= log(p_outlier + p_inlier);
   }

  return val;

}

void analytic_alpha_get_nl_centers_and_deltas(
  r3_t s, double err_s,
  double a,double err_a,
  double* center,double* delta,
  int num_parms)
{
    double X_s_min,X_s_max,Y_s_min,Y_s_max;
    r3_t view_dir = (r3_t){{ 0,0,1}};
    r3_and_error_to_stereographic_box(s,err_s,view_dir,&X_s_min,&X_s_max,&Y_s_min,&Y_s_max);
    center[0] = (X_s_min + X_s_max)/2.0;
    delta[0] = (X_s_max - X_s_min)/2.0;
    center[1] = (Y_s_min + Y_s_max)/2.0;
    delta[1] = (Y_s_max - Y_s_min)/2.0;
    if(num_parms > 2){
      center[2] = a;
      delta[2] = err_a;
    }

}


void analytic_alpha_generate_simplex( double center[],double delta[], int n, double* S){
  rmxn_regular_simplex(n, S);
  rmxn_spin_rows(n+1, n,S,S);
  double simplex_radius = rmxn_regular_simplex_radius(n);
  /*now modify S so that it lies inside the uncertainty intervals*/
  int i;
  for(i = 0; i <=n ; i++){
    double* Si = &(S[i*n]);
    int j;
    for(j = 0; j < n; j++){
      double rj = (0.125)*delta[j];
      Si[j] = (Si[j]/simplex_radius)*rj + center[j];
      assert(fabs(Si[j] - center[j]) <= 1.00001*rj);
    }
  }
  
}

double  analytic_alpha_compute_new_scalar_error(double t_old, double eps_old, double t, double alpha, double beta, double gamma);
double  analytic_alpha_compute_new_scalar_error(double t_old, double eps_old, double t, double alpha, double beta, double gamma){
  double delta = fabs(t - t_old); /*parameter change*/
  double eps = fmax(fmin(alpha*delta,beta*eps_old), gamma*eps_old); /*New uncertainty estimate*/ 
  return eps;
}

double analytic_alpha_compute_new_direction_error(r3_t u_old, double eps_old, r3_t u, double alpha, double beta, double gamma);
double analytic_alpha_compute_new_direction_error(r3_t u_old, double eps_old, r3_t u, double alpha, double beta, double gamma){
  double delta = acos(r3_dot(&u_old,&u));
  return fmin(fmax(fmin(alpha*delta,beta*eps_old), gamma*eps_old),M_PI/2.0);
}



double analytic_alpha_goal_function(analytic_alpha_log_prob_function_t* est_logprob ,r3_t normal,bool_t known_albedo,double albedo, double S[] ,int m, approx_model_t* lm, void** l_data, double sigma, double omg0, double omg1 ){
  double G_vec[m];
  int i;
  double val = 0;
  for(i = 0; i < m; i++){
      G_vec[i] = lm->evaluate(normal.c,l_data[i]);
      if(G_vec[i] < 0) { G_vec[i] = 0; }
      val+= G_vec[i];
  }

  
  if(known_albedo){
     return est_logprob(albedo,S,G_vec,m,sigma,omg0,omg1);
//       return analytic_alpha_log_prob_noalbedo02(albedo,S,G_vec,m,sigma,omg0,omg1);
  }else{
    return est_logprob(NAN,S,G_vec,m,sigma,omg0,omg1);
  }
}




void analytic_alpha_compute_normal(
  double S[],int m,analytic_inversion_t* ai, bool_t optimize_albedo,/*Input*/
  r3_t* normal, double* albedo,double* dist,double G[], /*Output*/
  int nLogProb,double sigma,double omg0, double omg1, int update_steps,double tol_albedo,double tol_normal, r3_t view_dir,char* debug_prefix
){

//   fprintf(stderr,"Start analytic_alpha_compute_normal\n");
//    int i;
//   fprintf(stderr,"S - ( ");
//   for(i = 0; i < m; i++){
//     fprintf(stderr," %9.6lf ",S[i]);
//   }
//   fprintf(stderr,")\n");
  
  
  /*Initial guess by analytic inversion*/
  double albedo_guess;
  r3_t normal_guess;
  double f_guess;
  double G_test[ai->m]; 
  
  analytic_inversion_compute_normal(S,ai,&normal_guess, &albedo_guess,G_test,&f_guess);
  if(debug_prefix){
    fprintf(stderr,"Initial Guesses Normal - ( %9.6lf %9.6lf %9.6lf ) Albedo ( %9.6lf ) \n",normal_guess.c[0],normal_guess.c[1],normal_guess.c[2],albedo_guess);
  }
  
  double rad_albedo = (optimize_albedo ? 0.3: 0.0) ;  /*Estimated distance from guessed albedo to optimal albedo*/
  double rad_normal = M_PI/9.0; /*Estimated distance from guessed normal to optimal normal*/
  
  
  double max_scale = 1.0;
  int nscales = 1;
  double alpha = sqrt(2.0);
  while (sigma*max_scale  < 0.2){
    max_scale*=alpha;
    nscales++;
  }
//   
  int iscale = nscales;
  
  r3_t normal_guess_scale = normal_guess;
  double albedo_guess_scale = albedo_guess;
  double scale = max_scale;
  
  while(iscale >= 1){
      char* debug_prefix_scale = NULL;
      if(debug_prefix){
	char *debug_prefix_scale = jsprintf("%s-S%02d",debug_prefix,iscale);
      }
      double sigma_scale = sigma*scale;
      double rad_normal_scale = fmin(rad_normal*scale,M_PI/4.0);
      double rad_albedo_scale = fmin(rad_albedo*scale,0.5);
      int update_steps_iter = (int)fmax(6,ceil(update_steps/scale));
      analytic_alpha_compute_normal_aux(S,ai, m,normal_guess_scale,rad_normal_scale,optimize_albedo, albedo_guess_scale,rad_albedo_scale,normal,albedo,dist,G,
					 nLogProb,sigma_scale,omg0,omg1,update_steps_iter,tol_albedo,tol_normal, view_dir,debug_prefix_scale);
      
      if(debug_prefix){
	fprintf(stderr,"Optimized Scale %02d Normal - ( %9.6lf %9.6lf %9.6lf ) Albedo ( %9.6lf ) \n",iscale,normal->c[0],normal->c[1],normal->c[2],*albedo);
      }	
      normal_guess_scale= *normal;
      albedo_guess_scale= *albedo;
      
      scale/=alpha;
      iscale--;
      free(debug_prefix_scale);
  }
//   fprintf(stderr,"end analytic_alpha_compute_normal\n");

  
}

void analytic_alpha_compute_normal_aux(
  double S[],analytic_inversion_t* ai,int m,r3_t normal_guess,double rad_normal,  bool_t optimize_albedo,double albedo_guess,double rad_albedo,/*Input*/
  r3_t* normal, double* albedo,double* dist,double G[], /*Output*/
  int nLogProb,double sigma,double omg0, double omg1, int update_steps,double tol_albedo,double tol_normal, r3_t view_dir,char* debug_prefix
){

//   fprintf(stderr,"Normal - ( %9.6lf %9.6lf %9.6lf ) Albedo %9.6lf Sigma %9.6lf\n",normal_guess.c[0],normal_guess.c[1],normal_guess.c[2],albedo_guess,sigma);
  
  analytic_alpha_log_prob_function_t* est_logprob = escolhe_analytic_alpha_log_prob_function_t(nLogProb);
    
  /*Validating the normal guess*/
  double thetaMAX = 0.50*M_PI;
  double ctMIN = cos(thetaMAX);
  if(r3_dot(&view_dir,&normal_guess) < ctMIN ){
     fprintf(stderr,"Impossible normal guess (%9.6lf %9.6lf %9.6lf) has angle too large with  with view dir (%9.6lf,%9.6lf,%9.6lf)\n",normal_guess.c[0],normal_guess.c[1],normal_guess.c[2],view_dir.c[0],view_dir.c[1],view_dir.c[2]);
     normal_guess = view_dir;
  }
  /*validating the albedo guess*/
  if((albedo_guess < 0) || (albedo_guess > 1.0) ){
    fprintf(stderr,"Impossible albedo guess - %9.6lf\n",albedo_guess);
    albedo_guess = fmin(1.0,fmax(0.0,albedo_guess));
  }
  
  int num_parms = 2 + (optimize_albedo);

  /*Albedo encoding/decoding*/
  /*The range [min_albedo,max_albedo] is mapped o [-INF,+INF] */
  
  double min_albedo = (optimize_albedo ? 0.0 : albedo_guess );
  double max_albedo = (optimize_albedo ? 1.0 : albedo_guess );
  double ALB_SCALE = 0.499*M_PI;
  
  auto double compute_rawz2_from_albedo(double alb);
  double compute_rawz2_from_albedo(double alb){
    assert(optimize_albedo);
    if( alb <= min_albedo) return -INF;
    if( alb >= max_albedo) return +INF;
    double rel = (alb - min_albedo)/(max_albedo - min_albedo);
    return tan(ALB_SCALE*(2*rel - 1));
  }
  
  auto double compute_albedo_from_rawz2(double z2);
  double compute_albedo_from_rawz2(double z2){
    assert(optimize_albedo);
    if( z2 <= -INF) return min_albedo;
    if( z2 >= INF ) return max_albedo;
    double rel = 0.5*(atan(z2)/ALB_SCALE +1);
    return  min_albedo + (rel*(max_albedo - min_albedo));
  }
    
  double rad_z2 = NAN; /*scale factor for  z2*/
  if(optimize_albedo){
    /*Compute the z2 scale factor*/
   double eps_albedo = 0.01*(max_albedo - min_albedo);
   double lo_alb = fmax(albedo_guess - rad_albedo,min_albedo + eps_albedo); 
   double hi_alb = fmin(albedo_guess + rad_albedo,max_albedo - eps_albedo);
   double lo_z2 = compute_rawz2_from_albedo(lo_alb);
   double hi_z2 = compute_rawz2_from_albedo(hi_alb);
   rad_z2 = (hi_z2 - lo_z2)/2.0; 
   if (debug_prefix){
     fprintf(stderr,"rad_z2 = %24.15e\n",rad_z2);
   }
  }
   
  /*Normal encoding-decoding*/
  double eps_z = 1.0e-6;
  r3x3_t M_guess_correction = compute_normal_correction_matrix(view_dir); 
  
  auto r2_t compute_rawz0z1_from_normal(r3_t normal);
  r2_t compute_rawz0z1_from_normal(r3_t normal){
    
    r3_t p;
    
    r3x3_map_col (&M_guess_correction,&normal, &p);
    r2_t s = (r2_t){{p.c[0],p.c[1]}};
    
    double ct = p.c[2] + eps_z;
    
    if(ct > eps_z ){
      r2_scale(1.0/ct,&s,&s);
    }else{
      r2_dir(&s,&s);
      r2_scale(1/eps_z,&s,&s);
    }
    return s;
  }
  
  auto r3_t compute_normal_from_rawz0z1(double z0,double z1);
  r3_t compute_normal_from_rawz0z1(double z0,double z1){
    r3_t p = (r3_t){{ z0, z1,1.0}};
    r3_dir(&p,&p);
    if(p.c[2] < eps_z ){
      p.c[2] = eps_z;
      r3_dir(&p,&p);
    }
    r3_t normal;
    r3x3_map_row(&p,&M_guess_correction,&normal);
    r3_dir(&normal,&normal);
    
    return normal;
  }
  
  double rad_z0_z1 = NAN; /*scale factor for z0 and z1 */
  {
    r3_t u1  = bend_towards(view_dir,normal_guess,-rad_normal);
    r3_t u2  = bend_towards(view_dir,normal_guess,+rad_normal);
    
    r2_t zu1 = compute_rawz0z1_from_normal(u1);
    r2_t zu2 = compute_rawz0z1_from_normal(u2);
    
    double dist_zu12 = r2_dist(&zu1,&zu2);
    
    assert(dist_zu12 > eps_z);
    rad_z0_z1 = dist_zu12/2.0;
    
  }

  
  
  auto void convert_normal_and_albedo_to_parms(r3_t normal, double alb,double z[], int num_parameters);
  void convert_normal_and_albedo_to_parms(r3_t normal, double alb,double z[], int num_parameters){
    assert(r3_dot(&normal,&view_dir) >= 0 );
    r2_t s = compute_rawz0z1_from_normal(normal);
    z[0] = s.c[0]/rad_z0_z1;
    z[1] = s.c[1]/rad_z0_z1;
    if(optimize_albedo){
       assert((alb <= max_albedo) && (alb >= min_albedo));
      z[2] = compute_rawz2_from_albedo(alb)/rad_z2;
    }
  }
  
  auto void convert_parm_to_normal_and_albedo(double z[], int num_parameters, r3_t* normal, double* alb);
  void convert_parm_to_normal_and_albedo(double z[], int num_parameters, r3_t* normal, double* alb){
    double z0 = z[0]*rad_z0_z1;     
    double z1 = z[1]*rad_z0_z1;
    *normal = compute_normal_from_rawz0z1(z0,z1);
    assert(r3_dot(normal,&view_dir) >=0 );
    if(optimize_albedo){
       *alb = compute_albedo_from_rawz2(z[2]*rad_z2);
       if( isnan(*alb) || (*alb > max_albedo) || (*alb < min_albedo)){
	 fprintf(stderr,"Impossible optimized albedo - %9.6lf \n",*alb);
       }
       
       assert((*alb <= max_albedo) && (*alb >= min_albedo));
    }
  }
  
  
  auto double  goal_function_alpha(double parameters[],int n);
  double  goal_function_alpha(double parameters[],int n){
    r3_t normal_estimated;
    double alb_estimated = NAN;
    convert_parm_to_normal_and_albedo(parameters, n, &normal_estimated,&alb_estimated);
    double val = 0;
    val =  analytic_alpha_goal_function(est_logprob,normal_estimated,optimize_albedo,alb_estimated, S,m,ai->gauge_model,ai->gauge_data,sigma,omg0,omg1 );

    return val;
  }
  
  
  double z_guess[num_parms];
  convert_normal_and_albedo_to_parms(normal_guess, albedo_guess,z_guess, num_parms); /*It is supposed to fill parms with zeros*/
  
  double r = 1.0;
  double theta = 0.5; 
  double alpha = 1.25;
  double epsilon = 1.0e-5;
  
  
//   fprintf(stderr,"--------------------------------------\n");
  double optmized_parameters[num_parms];
  nonlinear_optimize_z(z_guess,num_parms,r,theta,alpha,update_steps,epsilon,+1,goal_function_alpha,optmized_parameters,debug_prefix);
  //fprintf(stderr,"--------------------------------------\n");
  r3_t normal_approx ;
  double albedo_approx = NAN;
  convert_parm_to_normal_and_albedo(optmized_parameters, num_parms,&normal_approx, &albedo_approx);
  
  double go[m];
  double so[m];
  int i;
  for(i = 0; i < m; i++){
    G[i] = ai->gauge_model->evaluate(normal_approx.c,ai->gauge_data[i]);
  }
  
  double Smag  = rn_dir(m,S,so);
  double Gmag =  rn_dir(m,G,go);
  albedo_approx = (optimize_albedo ? albedo_approx : EstAlbedo_00(so, Smag,go,Gmag,m,sigma,omg0,omg1));
  if(r3_dot(&view_dir,&normal_approx) < ctMIN ){
     fprintf(stderr,"Impossible optimized (%9.6lf %9.6lf %9.6lf) has angle too large with view dir (%9.6lf,%9.6lf,%9.6lf)\n",normal_approx.c[0],normal_approx.c[1],normal_approx.c[2],view_dir.c[0],view_dir.c[1],view_dir.c[2]);
    albedo_approx = 0.63;
    r3_t para,perp;
    r3_decomp (&normal_approx, &view_dir, &para, &perp);
    normal_approx = perp;
    r3_dir(&normal_approx,&normal_approx);
     for(i = 0; i < m; i++){
      G[i] = ai->gauge_model->evaluate(normal_approx.c,ai->gauge_data[i]);
     }
//     assert(FALSE);
  }
  
  if(debug_prefix){
    plot_nonlinear_goalfunction(debug_prefix,optmized_parameters,NULL, num_parms, goal_function_alpha);
    float_image_t* im = analytic_alpha_plot_goal_2D(est_logprob, normal_approx, optimize_albedo,albedo_approx,S, m ,ai->gauge_model, ai->gauge_data,sigma,omg0,omg1,view_dir,100);
    char* filename = NULL;
    char *filename = jsprintf("%s-Goal.fni",debug_prefix);
    FILE* arq = open_write(filename,TRUE);
    float_image_write(arq,im);
    fclose(arq);
    float_image_free(im);
    free(filename);
  }
  *normal = normal_approx;
  *albedo = albedo_approx;
  *dist  = analytic_alpha_goal_function(est_logprob,normal_approx,optimize_albedo,albedo_approx, S,m,ai->gauge_model,ai->gauge_data, sigma,omg0,omg1 );
  
}

float_image_t* analytic_alpha_plot_goal_2D(analytic_alpha_log_prob_function_t* est_logprob, r3_t correct_normal, bool_t known_albedo, double albedo,double S[], int m , approx_model_t* gauge_model, void** gauge_data,double sigma,double omg0, double omg1,r3_t view_dir, int imSize){
  
  float_image_t* im = float_image_new(3,imSize,imSize);
  double Smax = 1.0;
  double radius = imSize/2.0;
  int x,y;
  double so[m];
  double Smag = rn_dir(m,S,so);
  for(x = 0; x < imSize; x++){
    for(y = 0; y < imSize; y++){


      double dX = Smax*(x - radius)/radius;
      double dY = Smax*(y - radius)/radius;
       double nx = 2*dX/((dX*dX) + (dY*dY) + 1);
       double ny = 2*dY/((dX*dX) + (dY*dY) + 1);
       double nz = (1 - (dX*dX) - (dY*dY))/((dX*dX) + (dY*dY) + 1);
       r3_t normal = (r3_t){{nx,ny,nz}};
       bool_t valid_hemisphere = (r3_dot(&normal,&view_dir) > 0 ? TRUE : FALSE);
      
      float_image_set_sample(im,0,x,y,0);
      if( valid_hemisphere){
	double val = analytic_alpha_goal_function(est_logprob,normal,known_albedo, albedo,S,m,gauge_model, gauge_data,sigma, omg0,omg1);
	float_image_set_sample(im,0,x,y,(float)val);
	double alb = albedo;
	if(!known_albedo){
	  double G[m];
	  double go[m];
	  double Gmag;
	  int i;
	  for(i = 0; i < m; i++){
	    G[i] = gauge_model->evaluate(normal.c,gauge_data[i]);
	    if(G[i] < 0) { G[i] = 0; }
	  }
	  Gmag = rn_dir(m,G,go);
	  alb = EstAlbedo_00(so, Smag,go,Gmag,m,sigma,omg0,omg1);
	}
	float_image_set_sample(im,1,x,y,(float)alb);
      }
      float_image_set_sample(im,2,x,y,0);
    }
  }
  double nx,ny;
  convert_r3_to_stereographic(correct_normal,&nx, &ny,view_dir);
  nx = (nx*radius) + radius;
  ny = (ny*radius) + radius;
  if( (nx >=0) && (nx < imSize) && (ny >= 0) && (ny < imSize) ){
    float_image_set_sample(im,2,(int)nx,(int)ny,1);
  }
  
  return im;
}

// float_image_t* analytic_alpha_plot_goal(r3_t correct_normal, bool_t known_albedo, double albedo,double S[], int m , approx_model_t* gauge_model, void** gauge_data,double sigma,double omg0, double omg1,r3_t view_dir, int imSize){
//   
//   float_image_t* im = float_image_new(2,imSize,imSize);
//   double Smax = 1.0;
//   double radius = imSize/2.0;
//   int x,y;
//   for(x = 0; x < imSize; x++){
//     for(y = 0; y < imSize; y++){
// 
// 
//       double dX = Smax*(x - radius)/radius;
//       double dY = Smax*(y - radius)/radius;
//        double nx = 2*dX/((dX*dX) + (dY*dY) + 1);
//        double ny = 2*dY/((dX*dX) + (dY*dY) + 1);
//        double nz = (1 - (dX*dX) - (dY*dY))/((dX*dX) + (dY*dY) + 1);
//        r3_t normal = (r3_t){{nx,ny,nz}};
//        bool_t valid_hemisphere = (r3_dot(&normal,&view_dir) > 0 ? TRUE : FALSE);
//       
//       float_image_set_sample(im,0,x,y,0);
//       if( valid_hemisphere){
// 	double val = analytic_alpha_goal_function(normal,known_albedo, albedo,S,m,gauge_model, gauge_data,sigma, omg0,omg1);
// 	float_image_set_sample(im,0,x,y,val);
//       }
//       float_image_set_sample(im,1,x,y,0);
//     }
//   }
//   double nx,ny;
//   convert_r3_to_stereographic(correct_normal,&nx, &ny,view_dir);
//   nx = (nx*radius) + radius;
//   ny = (ny*radius) + radius;
//   if( (nx >=0) && (nx < imSize) && (ny >= 0) && (ny < imSize) ){
//     float_image_set_sample(im,1,nx,ny,1);
//   }
//   
//   return im;
// }
