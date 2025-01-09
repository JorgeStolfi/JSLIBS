#define _GNU_SOURCE
#include <stdio.h>
#include <rn.h>
#include <stdlib.h>
#include <math.h>
#include <rmxn_extra.h>
#include <assert.h>
#include <sve_minn.h>
#include <nonlinear_optimization.h>
#include <jsfile.h>

void nonlinear_optimize(
  double center[], double delta[], int num_parms, /*Input*/
  double r, double rmin, double alpha, int max_iters,double epsilon,int which, /*Adjust parameters*/
  goal_function_t* goal_function, clip_function_t* clip_function, /*Callback functions*/
  double optmized_parameters[], /*Output*/
  char* debug_prefix
){
  if(debug_prefix){ fprintf(stderr,"\n MAX ITERS IS %d - %s\n",max_iters,debug_prefix);}
  auto void debug_parameters(FILE* arq,char* tag,int i, int j, double p[],double f);
  void debug_parameters(FILE* arq,char* tag,int i, int j, double p[],double f){
    assert(arq != NULL);
    fprintf(arq,"%s ",tag);
    fprintf(arq,"%5d %5d ",i,j);
    int rr;
    for(rr = 0; rr < num_parms; rr++){
      fprintf( arq,"%24.17lf ",p[rr]);
    }
    if(!isnan(f)){ fprintf( arq,"  %24.17lf ",f);}
    fprintf(arq,"\n");
  }
  
  assert(r < 1);
  assert(alpha > 1);
  int num_vertices = num_parms +1;
   /*Arrays for non linear optimization*/
  int num_samples = (num_vertices*(num_vertices+1))/2;
  double parameters[num_parms];
  double z[num_parms];
  double uline;
  double* Fv = rn_alloc(num_samples);
  double* S = rmxn_alloc(num_vertices,num_parms);
  
  /*init deltas*/
  double old_center[num_parms];
  double old_delta[num_parms];
  int i;
  
  int num_iters = 0;
  double diff_centers = INFINITY;
  while( (num_iters < max_iters) && (diff_centers > epsilon) ){
    
    FILE* arq_debug = NULL;
    if(debug_prefix){ 
      char* fname_debug = NULL;
      char *fname_debug = jsprintf("%s-NL-I%05d.txt",debug_prefix,num_iters);
      fprintf(stderr,"we have to open %s here\n",fname_debug);
      arq_debug = open_write(fname_debug,TRUE);
      free(fname_debug);
    }
    
    if(arq_debug != NULL){  debug_parameters(arq_debug,"GC",-1,-1,center,NAN); }
    if(arq_debug != NULL){  debug_parameters(arq_debug,"GD",-1,-1,delta,NAN); }
    /*Copy previous iteration of centers*/
    rn_copy(num_parms,center,old_center);
    rn_copy(num_parms,delta,old_delta);
    /*Creates a regular simplex*/
    rmxn_regular_simplex(num_parms, S);
    rmxn_spin_rows(num_vertices, num_parms,S,S);
    double radius_simplex = rmxn_regular_simplex_radius(num_parms);
        
    /*Compresses the simplex*/
    for(i = 0; i < num_vertices ; i++){
      double* Si = &(S[i*num_parms]);
      int j;
      for(j = 0; j < num_parms; j++){
	Si[j] = (r/radius_simplex)*(Si[j]*delta[j]) + center[j];
      }
    }
        
    /*Optmize*/
    int k = 0;
    
    int ibest,jbest;
    ibest = jbest = -1;
    double fbest = 0;
    
    for(i = 0; i < num_vertices; i++){
      double* Si = &(S[i*num_parms]);
      int j;
      for(j = 0; j <= i; j++){
	double* Sj = &(S[j*num_parms]);
	int rr;
	for(rr = 0; rr < num_parms; rr++){
	  parameters[rr] = (i == j ? Si[rr] : (Si[rr] + Sj[rr])/2.0);
	}
 	
	Fv[k] = goal_function(parameters,num_parms);
	if(arq_debug != NULL){ debug_parameters(arq_debug,"SP",i,j,parameters,Fv[k]); }
	if( which != 0){
	  if( ibest == -1){
	    ibest = i;
	    jbest = j;
	    fbest = Fv[k];
	  }else {
	    bool_t update = (which == -1 ? fbest > Fv[k] : fbest < Fv[k]);
	    fbest = (update ? Fv[k] : fbest);
	    ibest = (update ? i : ibest);
	    jbest = (update ? j : jbest);
	  }
	}
	
	k++;
      }
    }
    double baricentric[num_vertices];
    sve_minn_step(num_parms,Fv, baricentric);
    rmxn_map_row (num_vertices, num_parms, baricentric,S, parameters);
    
    if(arq_debug != NULL){  debug_parameters(arq_debug,"OP",-1,-1,parameters,NAN); }
    
    clip_function(parameters,delta,num_parms);
    if(arq_debug != NULL){  debug_parameters(arq_debug,"CP",-1,-1,parameters,NAN); }
    /*Decompresses the simplex parameters*/
    for(i = 0; i < num_parms; i++){
      z[i] = (parameters[i] - center[i])/delta[i];
    }
    
    /*trunkate z if needed*/
    double modZ = rn_norm(num_parms,z);
    if(modZ > 1){
      rn_dir(num_parms,z,z); /*It truncates z, but i'm not sure if it is correct ...*/
      modZ = 1.0;
    }
    if(debug_prefix) fprintf(stderr,"modz = %lf and thing is %lf\n",modZ,alpha*fmin(modZ,1 - modZ));
    /*Compute U line*/
//     uline = fmax(alpha*fmin(modZ,1 - modZ),rmin);
      uline = fmax(alpha*modZ,rmin);
    /*Update center and deltas*/
    for(i = 0; i < num_parms; i++){
      center[i] = (z[i]*old_delta[i]) + old_center[i];
      delta[i] = uline*old_delta[i];
    }
    clip_function(center, delta,num_parms);
    double fcurrent = goal_function(center,num_parms);
    
    if(arq_debug != NULL){  debug_parameters(arq_debug,"OC",-1,-1,center,fcurrent); }
    if(arq_debug != NULL){  debug_parameters(arq_debug,"OD",-1,-1,delta,uline); }
    /*Choose the minimum or max instead*/
    if(which != 0){
      double* Si = &(S[ibest*num_parms]);
      double* Sj = &(S[jbest*num_parms]);
      
      bool_t update = (which == -1 ? fcurrent > fbest : fcurrent < fbest);
      if(update){
	int rr;
	for(rr = 0; rr < num_parms; rr++){
	  center[rr] = (ibest == jbest ? Si[rr] : (Si[rr] + Sj[rr])/2.0);
	}
	if(arq_debug != NULL){  debug_parameters(arq_debug,"UC",ibest,jbest,center,fbest);} 
	clip_function(center, delta,num_parms);
	fcurrent = goal_function(center,num_parms);
	if(arq_debug != NULL){  debug_parameters(arq_debug,"CC",ibest,jbest,center,fcurrent); }
      }
    }

    diff_centers = rn_dist(num_parms,center, old_center);
    if(debug_prefix) { fclose(arq_debug);}
    num_iters++;
    
  }
  /*Clean the arrays*/
 free(S);
 free(Fv);
 
 rn_copy(num_parms,center,optmized_parameters);

}


void nonlinear_optimize_z(
  double initial_z[], int num_parms, /*Input*/
  double r,double theta, double alpha, int max_iters,double epsilon,int which, /*Adjust parameters*/
  goal_function_t* goal_function, /*Callback functions*/
  double opt_z[], /*Output*/
  char* debug_prefix
){
  
  if(debug_prefix){ fprintf(stderr,"\n MAX ITERS IS %d - %s\n",max_iters,debug_prefix);}
  auto void debug_parameters(FILE* arq,char* tag,int i, int j, double p[],double f);
  void debug_parameters(FILE* arq,char* tag,int i, int j, double p[],double f){
    assert(arq != NULL);
    fprintf(arq,"%s ",tag);
    fprintf(arq,"%5d %5d ",i,j);
    int rr;
    for(rr = 0; rr < num_parms; rr++){
      fprintf( arq,"%24.17lf ",p[rr]);
    }
    if(!isnan(f)){ fprintf( arq,"  %24.17lf ",f);}
    fprintf(arq,"\n");
  }
  
//   assert(r < 1);
  assert(theta < 1);
  int num_vertices = num_parms +1;
   /*Arrays for non linear optimization*/
  int num_samples = (num_vertices*(num_vertices+1))/2;
  double* Fv = rn_alloc(num_samples);
  double* S = rmxn_alloc(num_vertices,num_parms);
  
  /*init deltas*/
  double z[num_parms];
  rn_copy(num_parms,initial_z,z);
  double z_star[num_parms];
  int i;
  
  int num_iters = 0;
  double diff_centers = INFINITY;
  while( (num_iters < max_iters) && (diff_centers > epsilon) ){
    
    FILE* arq_debug = NULL;
    if(debug_prefix){ 
      char* fname_debug = NULL;
      char *fname_debug = jsprintf("%s-NL-I%05d.txt",debug_prefix,num_iters);
      arq_debug = open_write(fname_debug,TRUE);
      free(fname_debug);
    }
    
    if(arq_debug != NULL){  debug_parameters(arq_debug,"GZ",-1,-1,z,NAN); }
    
    /*Creates a regular simplex*/
    rmxn_regular_simplex(num_parms, S);
    rmxn_spin_rows(num_vertices, num_parms,S,S);
    double radius_simplex = rmxn_regular_simplex_radius(num_parms);
        
    /*Compresses the simplex*/
    for(i = 0; i < num_vertices ; i++){
      double* Si = &(S[i*num_parms]);
      int j;
      for(j = 0; j < num_parms; j++){
	Si[j] = (theta*r/radius_simplex)*(Si[j]) + z[j];
      }
    }
        
    /*Optmize*/
    int k = 0;
    
    int ibest,jbest;
    ibest = jbest = -1;
    double fbest = 0;
    
    for(i = 0; i < num_vertices; i++){
      double* Si = &(S[i*num_parms]);
      int j;
      for(j = 0; j <= i; j++){
	double* Sj = &(S[j*num_parms]);
	int rr;
	for(rr = 0; rr < num_parms; rr++){
	  z_star[rr] = (i == j ? Si[rr] : (Si[rr] + Sj[rr])/2.0);
	}
 	
	Fv[k] = goal_function(z_star,num_parms);
	if(arq_debug != NULL){ debug_parameters(arq_debug,"SP",i,j,z_star,Fv[k]); }
	if( which != 0){
	  if( ibest == -1){
	    ibest = i;
	    jbest = j;
	    fbest = Fv[k];
	  }else {
	    bool_t update = (which == -1 ? fbest > Fv[k] : fbest < Fv[k]);
	    fbest = (update ? Fv[k] : fbest);
	    ibest = (update ? i : ibest);
	    jbest = (update ? j : jbest);
	  }
	}
	
	k++;
      }
    }
    double baricentric[num_vertices];
    sve_minn_step(num_parms,Fv, baricentric);
    rmxn_map_row (num_vertices, num_parms, baricentric,S, z_star);
    
    if(arq_debug != NULL){  debug_parameters(arq_debug,"OP",-1,-1,z_star,NAN); }
    /*Decompresses the simplex parameters*/
    /*trunkate z if needed*/
    double distZ = rn_dist(num_parms,z_star,z);
    if(distZ > r){
      double diff_z[num_parms];
      rn_sub(num_parms,z_star,z,diff_z);
      rn_mix(num_parms, 1.0, z,r/distZ,diff_z,z_star);
      distZ = r;
    }
    
      
    /*Update z and r*/
    rn_copy(num_parms,z_star,z);
    r = alpha*distZ;
    double fcurrent = goal_function(z,num_parms);
    if(debug_prefix) fprintf(stderr,"DIFFZ is %12.7lf\nR is %12.7lf\n",distZ,r);
    
    if(arq_debug != NULL){  debug_parameters(arq_debug,"OZ",-1,-1,z,fcurrent); }
    /*Choose the minimum or max instead*/
    if(which != 0){
      double* Si = &(S[ibest*num_parms]);
      double* Sj = &(S[jbest*num_parms]);
      
      bool_t update = (which == -1 ? fcurrent > fbest : fcurrent < fbest);
      if(update){
	int rr;
	for(rr = 0; rr < num_parms; rr++){
	  z[rr] = (ibest == jbest ? Si[rr] : (Si[rr] + Sj[rr])/2.0);
	}
	if(arq_debug != NULL){  debug_parameters(arq_debug,"UZ",ibest,jbest,z,fbest);} 
	fcurrent = goal_function(z_star,num_parms);
	if(arq_debug != NULL){  debug_parameters(arq_debug,"CC",ibest,jbest,z,fcurrent); }
      }
    }
    
    
    diff_centers = distZ;
    
    if(debug_prefix) { fclose(arq_debug);}
    num_iters++;
    
  }
  /*Clean the arrays*/
 free(S);
 free(Fv);
 
 rn_copy(num_parms,z,opt_z);

}


void plot_nonlinear_goalfunction(char* prefix,double center[],double delta[], int num_parms, goal_function_t goal_function)
{
  
  double parameter[num_parms];
 
  int nsteps = 20;
  
  
   
  int i_par;
  for(i_par = 0; i_par < num_parms; i_par++){
    parameter[i_par] = center[i_par];
  }
  for(i_par = 0; i_par < num_parms; i_par++){
    /*plots variation with respect to parameter i_par*/
    char* filename = NULL;
    char *filename = jsprintf("%s-P%02d.txt",prefix,i_par);
    FILE* arq = open_write(filename,TRUE);
    double v_old = parameter[i_par];
    int j;
    for( j = -nsteps; j<= +nsteps; j++){
      if(delta != NULL){
	parameter[i_par] = center[i_par]+delta[i_par]*(j/(double)nsteps);
      }else{
	parameter[i_par] = center[i_par]+ (j/(double)nsteps);
      }
      double f = goal_function(parameter,num_parms);
      fprintf(arq,"%15.9lf %15.9lf\n",parameter[i_par],f);
    }
    parameter[i_par] = v_old;
    fclose(arq);
    free(filename);
  }
  
}