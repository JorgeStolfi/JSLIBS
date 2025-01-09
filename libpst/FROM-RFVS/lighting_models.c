#include <lighting_models.h>
#include <assert.h>
#include <rn.h>
#include <pst_lamp.h>
#include <lighting_harmonic.h>
#include <lighting_radial.h>
#include <lighting_stereopoly.h>
#include <lighting_compact.h>
#include <lighting_nearest.h>
#include <jsfile.h>



/*Argpareser reader*/

double lighting_error_parse(argparser_t* pp,double emax){
  if(argparser_keyword_present_next(pp,"error")){
    return argparser_get_next_double(pp,0,emax);
  }else{
    return 0;
  }
}



/*Generic model type*/

approx_model_t* create_approx_lighting_model(model_type_t type){
  
  
  if( type == COMPACT_MODEL ) return lighting_compact_create_approx_lighting_model();
  else if( type == HARMONIC_MODEL ) return lighting_harmonic_create_approx_lighting_model();
  else if( type == RADIALBASIS_MODEL ) return lighting_radial_create_approx_lighting_model();
  else if( type == STEREOPOLY_MODEL ) return lighting_stereopoly_create_approx_lighting_model();
  else if( type == NEAREST_MODEL ) return lighting_nearest_create_approx_lighting_model();
  demand(FALSE, "unknown lighting model");
    
  return NULL;
}


void generate_compare_plot(FILE* arq, approx_model_t* lm_test,void* l_datatest, approx_model_t* lm_approx, void* l_datapprox, long int num_samples, r3_t view_dir, r3_t ref_dir){
  long int i;
  r3_t para;
  r3_decomp(&ref_dir,&view_dir,&para,&ref_dir);
  double temp= r3_dir(&ref_dir,&ref_dir); 
  demand(temp > 0.0001,"Invalid ref_dir");
  for(i = 0; i < num_samples;i++){
    
    double theta = ((2*i)/(double)num_samples -1)*M_PI;
    double st = sin(theta);
    double ct = cos(theta);
    r3_t normal;
    r3_mix(ct,&view_dir,st,&ref_dir,&normal);
    double temp2 = r3_dir(&normal,&normal);
    assert(fabs(temp2 -1) < 0.00001);
    fprintf(arq,"%+9.6lf %+1.6lf %+1.6lf %+1.6lf ",theta,normal.c[0],normal.c[1],normal.c[2]);
    double f = lm_test->evaluate(normal.c,l_datatest);
    double fi = lm_approx->evaluate(normal.c,l_datapprox);
     if(isnan(f)){
      fprintf(stderr,"NAN value for input ! \n");
      f = 0;
    }
    if(isnan(fi)){
      fprintf(stderr,"NAN value for approximation ! \n");
      fi = 0;
    }
    fprintf(arq,"%+9.6lf %+9.6lf  ",f,fi);
    
     /*First Phi values of test model*/
    long int j;
    long int components_fi = lm_approx->get_num_components(l_datapprox);
    double s = 0;
    double alphas[components_fi];
    lm_approx->get_alphas(l_datapprox,alphas,components_fi);
    for(j = 0; j < components_fi; j++){
      double alpha = alphas[j];
      double val_phi = lm_approx->phi(j,normal.c,l_datapprox);
      fprintf(arq,"%+9.6lf ",alpha*val_phi);
      s+= alpha*val_phi;
    }
    if( (fabs(s - fi)/fabs(s + fi + 1.0e-200))  > 1.0e-12){
      fprintf(stderr,"Approx for theta = %+9.6lf inconsistent. %+21.15lf %+21.15lf\n",theta,fi,s);
    }
    fprintf(arq,"\n");
  }
  
}


void   generate_SG_plot(char* prefix, approx_model_t* lm_test,void* l_datatest, approx_model_t* lm_approx, void* l_datapprox, long int num_samples_axis,r3_t view_dir){
  long int X,Y;
  char* filename = NULL;
  r3_t normal;
  auto double func_plot(int id, int component);
  auto void init_plot(int id, int component);
   int FUNC_F = 0;
   int FUNC_FI = 1;
   int FUNC_E = 2;
   int FUNC_PHI = 3;
   int FUNC_NONE = 4;
  
  
   void init_plot(int id, int component){
      char* prefix_name;
      char* component_name  = "\0";
      if(id == FUNC_F){
	prefix_name = "F";
      }else if(id == FUNC_FI){
	prefix_name = "FI";
      }else if (id == FUNC_E){
	prefix_name = "E";
      }else if(id == FUNC_PHI){
	prefix_name = "A";
	component_name = NULL;
	char *component_name = jsprintf("_PHI%03d",component);
      }else if (id == FUNC_NONE){
	prefix_name = "NONE";
      }else{
	fprintf(stderr,"ERROR -NOT A DEFINED FUNCTION FOR PLOT - %d\n",id);
	assert(FALSE);
      }
      char *filename = jsprintf("%s-%s%s.fni",prefix,prefix_name,component_name);
    }
   
  double func_plot(int id, int component){  
    double val;
    if(id == FUNC_F){
      val = lm_test->evaluate(normal.c,l_datatest);
    }else if(id == FUNC_FI){
      val = lm_approx->evaluate(normal.c,l_datapprox);
    }else if (id == FUNC_E){
      double v_fi = lm_approx->evaluate(normal.c,l_datapprox);
      double v_f = lm_test->evaluate(normal.c,l_datatest);
      val = v_fi - v_f;
    }else if(id == FUNC_PHI){
      int n_compe = lm_approx->get_num_components(l_datapprox);
      double alphas[n_compe];
      lm_approx->get_alphas(l_datapprox,alphas,n_compe);
      double alpha = alphas[component];
      double val_phi = lm_approx->phi(component,normal.c,l_datapprox);
      val = val_phi;
    }else if(id == FUNC_NONE){
      val = 1;
    }else{
      fprintf(stderr,"ERROR -NOT A DEFINED FUNCTION FOR PLOT - %d\n",id);
      assert(FALSE);
    }
    return val;
  }
  
  
  
  int funcao;
  double Smax = 1.0;
  double thetaMAX = M_PI/2.0;
  for(funcao = FUNC_F; funcao < FUNC_NONE; funcao++){
    
    long int num_components = 1;
    if (funcao == FUNC_PHI) num_components = lm_approx->get_num_components(l_datapprox);
    int j;
    for( j = 0; j < num_components; j++){
      float_image_t* plot_im = float_image_new(1,num_samples_axis,num_samples_axis);
      init_plot(funcao,j);
      for(X = 0; X < num_samples_axis;X++){
	for(Y = 0; Y < num_samples_axis; Y++){
	  double radius = num_samples_axis/2.0;
	  double dX = Smax*(X - radius)/radius;
	  double dY = Smax*(Y - radius)/radius;
	  
//           bool_t north_hemisphere = sqrt((dX*dX) + (dY*dY)) <  1 ;
	  
	  
	  double func_val = 0;
	  double x = 2*dX/((dX*dX) + (dY*dY) + 1);
	  double y = 2*dY/((dX*dX) + (dY*dY) + 1);
	  double z = (1 - (dX*dX) - (dY*dY))/((dX*dX) + (dY*dY) + 1);
	  normal = (r3_t){{x,y,z}};
	  bool_t valid_hemisphere = (r3_dot(&normal,&view_dir) > cos(thetaMAX) ? TRUE : FALSE);
// 	  valid_hemisphere = TRUE;
	  func_val = func_plot(funcao,j);
	  if(!valid_hemisphere) func_val = 0;
	  float_image_set_sample(plot_im,0,X,Y,func_val);
	}
      }
      FILE* arq = open_write(filename,TRUE);
      float_image_write(arq,plot_im);
      fclose(arq);
      float_image_free(plot_im);
      free(filename);
      filename=NULL;
   }
  }
  
 
  
  
}





void lighting_model_generate_plot(FILE* arq, approx_model_t* lm_test,void* l_datatest,long int num_samples, r3_t view_dir, r3_t ref_dir){
  long int i;
  r3_t para;
  r3_decomp(&ref_dir,&view_dir,&para,&ref_dir);
  double temp= r3_dir(&ref_dir,&ref_dir); 
  demand(temp > 0.0001,"Invalid ref_dir");
  for(i = 0; i < num_samples;i++){
    
    double theta = ((2*i)/(double)num_samples -1)*M_PI;
    double st = sin(theta);
    double ct = cos(theta);
    r3_t normal;
    r3_mix(ct,&view_dir,st,&ref_dir,&normal);
    double temp2 = r3_dir(&normal,&normal);
    assert(fabs(temp2 -1) < 0.00001);
    fprintf(arq,"%+9.6lf %+1.6lf %+1.6lf %+1.6lf ",theta,normal.c[0],normal.c[1],normal.c[2]);
    double f = lm_test->evaluate(normal.c,l_datatest);
    if(isnan(f)){
      fprintf(stderr,"NAN value for input ! \n");
      f = 0;
    }
    fprintf(arq,"%+9.6lf ",f);
     /*First Phi values of test model*/
    long int j;
    long int components_f = lm_test->get_num_components(l_datatest);
    double s = 0;
    double alphas[components_f];
    lm_test->get_alphas(l_datatest,alphas,components_f);
    for(j = 0; j < components_f; j++){
      double alpha = alphas[j];
      double val_phi = lm_test->phi(j,normal.c,l_datatest);
      fprintf(arq,"%+9.6lf ",alpha*val_phi);
      s+= alpha*val_phi;
    }
    if( (fabs(s - f)/fabs(s + f + 1.0e-200))  > 1.0e-12){
      fprintf(stderr,"Theta = %+9.6lf inconsistent. %+21.15lf %+21.15lf\n",theta,f,s);
    }
    fprintf(arq,"\n");
  }
  
}


void   lighting_model_generate_SG_plot(char* prefix, approx_model_t* lm_test,void* l_datatest,long int num_samples_axis,r3_t view_dir){
  long int X,Y;
  char* filename = NULL;
  r3_t normal;
  auto double func_plot(int id, int component);
  auto void init_plot(int id, int component);
   int FUNC_F = 0;
   int FUNC_PHI = 1;
   int FUNC_NONE = 2;
  
  
   void init_plot(int id, int component){
      char* prefix_name;
      char* component_name  = "\0";
      if(id == FUNC_F){
	prefix_name = "F";
      }else if(id == FUNC_PHI){
	prefix_name = "A";
	component_name = NULL;
	char *component_name = jsprintf("_PHI%03d",component);
      }else if (id == FUNC_NONE){
	prefix_name = "NONE";
      }else{
	fprintf(stderr,"ERROR -NOT A DEFINED FUNCTION FOR PLOT - %d\n",id);
	assert(FALSE);
      }
      char *filename = jsprintf("%s-%s%s.fni",prefix,prefix_name,component_name);
    }
   
  double func_plot(int id, int component){  
    double val;
    if(id == FUNC_F){
      val = lm_test->evaluate(normal.c,l_datatest);
    }else if(id == FUNC_PHI){
      int n_compe = lm_test->get_num_components(l_datatest);
      double alphas[n_compe];
      lm_test->get_alphas(l_datatest,alphas,n_compe);
      double alpha = alphas[component];
      double val_phi = lm_test->phi(component,normal.c,l_datatest);
      val = alpha*val_phi;
    }else if(id == FUNC_NONE){
      val = 1;
    }else{
      fprintf(stderr,"ERROR -NOT A DEFINED FUNCTION FOR PLOT - %d\n",id);
      assert(FALSE);
    }
    return val;
  }
  
  
  
  int funcao;
  double Smax = 1.5;
  double thetaMAX = M_PI/2.0;
  for(funcao = FUNC_F; funcao < FUNC_NONE; funcao++){
    
    long int num_components = 1;
    if (funcao == FUNC_PHI) num_components = lm_test->get_num_components(l_datatest);
    int j;
    for( j = 0; j < num_components; j++){
      float_image_t* plot_im = float_image_new(1,num_samples_axis,num_samples_axis);
      init_plot(funcao,j);
      for(X = 0; X < num_samples_axis;X++){
	for(Y = 0; Y < num_samples_axis; Y++){
	  double radius = num_samples_axis/2.0;
	  double dX = Smax*(X - radius)/radius;
	  double dY = Smax*(Y - radius)/radius;
	  
//           bool_t north_hemisphere = sqrt((dX*dX) + (dY*dY)) <  1 ;
	  
	  
	  double func_val = 0;
	  double x = 2*dX/((dX*dX) + (dY*dY) + 1);
	  double y = 2*dY/((dX*dX) + (dY*dY) + 1);
	  double z = (1 - (dX*dX) - (dY*dY))/((dX*dX) + (dY*dY) + 1);
	  normal = (r3_t){{x,y,z}};
	  bool_t valid_hemisphere = (r3_dot(&normal,&view_dir) > cos(thetaMAX) ? TRUE : FALSE);
	  
	  func_val = func_plot(funcao,j);
	  if(!valid_hemisphere) func_val = 0;
	  float_image_set_sample(plot_im,0,X,Y,func_val);
	}
      }
      FILE* arq = open_write(filename,TRUE);
      float_image_write(arq,plot_im);
      fclose(arq);
      float_image_free(plot_im);
      free(filename);
      filename=NULL;
   }
  }
  
 
  
  
}

