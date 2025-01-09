#include <normal_interpolation.h>
#include <polynomial_functions.h>
#include <r2.h>
#include <affirm.h>

r3_t float_image_get_normal(float_image_t* fim, int x,int y);
void float_image_set_normal(float_image_t* fim, int x,int y,r3_t normal);


r3_t float_image_get_normal(float_image_t* fim, int x,int y){
  r3_t normal;
  normal.c[0] = float_image_get_sample(fim,0,x,y);
  normal.c[1] = float_image_get_sample(fim,1,x,y);
  normal.c[2] = float_image_get_sample(fim,2,x,y);
  return normal;
}

void float_image_set_normal(float_image_t* fim, int x,int y,r3_t normal){
  float_image_set_sample(fim,0,x,y,(float)(normal.c[0]));
  float_image_set_sample(fim,1,x,y,(float)(normal.c[1]));
  float_image_set_sample(fim,2,x,y,(float)(normal.c[2]));
}

float_image_t* normal_interpolate_prob(float_image_t** fim_normals, float_image_t** fim_logprobs, int n, normal_interp_t interp_opt,float_image_t* select_map){
  int NX = (int)(fim_normals[0]->sz[1]);
  int NY = (int)(fim_normals[0]->sz[2]);
  float_image_t* fim_normal_it = float_image_new(3,NX,NY);
  int x,y;
  for(x = 0 ;  x < NX; x++){
    for(y = 0; y < NY; y++){

      r3_t normal_it;
      r3_t normals[n];
      double logprobs[n];
      int i;
      for(i = 0; i < n;i++){
	normals[i] = float_image_get_normal(fim_normals[i],x,y);
	logprobs[i] = float_image_get_sample(fim_logprobs[i],0,x,y);
      }

      /*Option 1 - choose the minor logprob*/
      if(interp_opt == PROB_BEST){
	int best_prob = 0;
	for(i = 1; i < n; i++){
	  if(logprobs[i] > logprobs[best_prob]){
	    best_prob = i;
	  }
	}
	if(select_map != NULL){
          float val = (float)((double)best_prob/(double)n);
	  float_image_set_sample(select_map,0,x,y,val);
	}
	normal_it = normals[best_prob];
      }else if(interp_opt == PROB_AVERAGE){
	normal_it = (r3_t){{0,0,0}};
	for(i = 0; i <n;i++){
	  double w = exp(logprobs[i]);
// 	  double w = exp(-logprob[i]);
	  normal_it.c[0]+=w*normals[i].c[0];
	  normal_it.c[1]+=w*normals[i].c[1];
	  normal_it.c[2]+=w*normals[i].c[2];
	}
	double len = r3_dir(&normal_it,&normal_it);
	if(len < 0.01) r3_zero(&normal_it);
	
      }

      float_image_set_normal(fim_normal_it,x,y,normal_it);
      
    }
  }
  return fim_normal_it;
}


float_image_t* normal_interpolate_pos(float_image_t** fim_normals, r2_t* gauge_positions, int n, int degree){
  int NX = (int)fim_normals[0]->sz[1];
  int NY = (int)fim_normals[0]->sz[2];
  float_image_t* fim_normal_it = float_image_new(3,NX,NY);
  int x,y;
  
  poly_function_t* Vx;
  poly_function_t* Vy;
  poly_function_t* Vz;
  
  Vx = poly_init_components(2,degree,FALSE);
  Vy = poly_init_components(2,degree,FALSE);
  Vz = poly_init_components(2,degree,FALSE);
  approx_model_t* lm = create_ls_polynomial_model();
  
  demand( Vx->num_coefs <= n,"Invalid number of gauges for degree !");
  
  double** X = (double**)malloc(sizeof(double*)*n);
  double Fx[n];
  double Fy[n];
  double Fz[n];
  int i;
  for(i = 0; i < n; i++){
    X[i] = (double*)malloc(sizeof(double)*2);
    X[i][0] = gauge_positions[i].c[0];
    X[i][1] = gauge_positions[i].c[1];
  }
  fprintf(stderr,"\n\n");
  int count = 0;
  for(x = 0 ;  x < NX; x++){
    for(y = 0; y < NY; y++){
      /*First arrange the system*/
      /*We dont want realloc every time that we solve a system, so we just reset the weights*/
      fprintf(stderr,"\033[1A");
      fprintf(stderr,"processed %d pixels of %d\n",count,NX*NY);
      count++;
      for(i = 0; i < n;i++){
	r3_t normal = float_image_get_normal(fim_normals[i],x,y);
	Fx[i] = normal.c[0];
	Fy[i] = normal.c[1];
	Fz[i] = normal.c[2];
      }
      /*solve the system*/
     approx_model_fit_linear_model(X,Fx,NULL,NULL,n,lm,Vx);
     approx_model_fit_linear_model(X,Fy,NULL,NULL,n,lm,Vy);
     approx_model_fit_linear_model(X,Fz,NULL,NULL,n,lm,Vz);
      /*
      fitModelToFunction(X,Fx,n,lm,Vx,3);
      fitModelToFunction(X,Fy,n,lm,Vy,3);
      fitModelToFunction(X,Fz,n,lm,Vz,3);*/
      
      r3_t normal_it;
      double p[2];
      p[0] = x;
      p[1] = y;
      normal_it.c[0] = lm->evaluate(p,Vx);
      normal_it.c[1] = lm->evaluate(p,Vy);
      normal_it.c[2] = lm->evaluate(p,Vz);
      
      r3_dir(&normal_it,&normal_it);
      
      float_image_set_normal(fim_normal_it,x,y,normal_it);
      
      for(i = 0; i < Vx->num_coefs;i++){
	Vx->weights[i] = 1.0;
	Vy->weights[i] = 1.0;
	Vz->weights[i] = 1.0;
      }
    }
  }
  return fim_normal_it; 
}

float_image_t* normal_interpolate_hightlight12(float_image_t** fim_normals,float_image_t** fim_weights, int n,r3_t* cluster_dir,double k){
  int NX,NY;
  int NX = (int)(fim_normals[0]->sz[1]);
  int NY = (int)(fim_normals[0]->sz[2]);
  float_image_t* fim_normal_it = float_image_new(3,NX,NY);
  int x,y,i;
  
  double r = 5.0; //TEMPORARY 
  
  //calcula NHIGH = r3_dir(DIRCLUSTER + (0,0,1))
  r3_t Nhight[n];
  
  for(i = 0; i < n; i++){
    Nhight[i] = cluster_dir[i];
    Nhight[i].c[2]+= 1;
    r3_dir(&(Nhight[i]),&(Nhight[i]));
  }
  
  for(x = 0 ;  x < NX; x++){
    for(y = 0; y < NY; y++){
		
	r3_t norms[n];
	double w[n];
	double wh[n];
	//int i;
	r3_t norm_it;
	r3_zero(&norm_it);
	for(i = 0; i < n; i++){
	    norms[i] = float_image_get_normal(fim_normals[i],x,y);
	    wh[i] = w[i] = float_image_get_sample(fim_weights[i],0,x,y);
	    r3_mix_in (w[i],&(norms[i]),&norm_it);
	    /* Sets {r := r + s * a}. */
	}
	
	double len  = r3_dir(&norm_it,&norm_it);
	if(len < 0.01) {
	  r3_zero(&norm_it);
	  float_image_set_normal(fim_normal_it,x,y,norm_it);
	  continue;
	}
	
	r3_t norm_previous;
	double w_previous[n];
	double diff = 1.0;
	double epsilon = 10e-10;
		
	int num_iter = 0;
	while((diff > epsilon) && (num_iter < 200)){
	  //save previous data
	  norm_previous = norm_it;
	  rn_copy(n,wh,w_previous);
	  //update weights
	  for(i = 0; i < n; i++){
	    double T = r3_dot(&(Nhight[i]),&norm_it);
	    T = cos(acos(T) + (r*M_PI/180.0));
// 	    double wh = (T <= 0 ? 1 : pow(1-T,k));
	    double whigh = (T <= 0 ? 1 : 1 - pow(T,k));
	    wh[i] = whigh*w[i];
	  }
	  //update  computed normal
	  r3_zero(&norm_it);
	  for(i = 0; i < n; i++){
	    r3_mix_in (wh[i],&(norms[i]),&norm_it);
	  }
	  r3_dir(&norm_it,&norm_it);

	  //compute the difference
	  diff = r3_dist(&norm_it,&norm_previous);
	 // diff = rn_dist(n,wh,w_previous);
// 	 fprintf(stderr,"[%d,%d] - ");
// 	 rn_gen_print(stderr,n,wh, "%8.6f","("," , ",") - ");
// 	 r3_gen_print(stderr,&norm_it,
	num_iter++;
      }
	
	float_image_set_normal(fim_normal_it,x,y,norm_it);
	
   }
  }
  return fim_normal_it;
}
