#define _GNU_SOURCE
#include <jsfile.h>
#include <stdio.h>
#include <lighting_models.h>
#include <polynomial_functions.h>
#include <analytic_inversion.h>
#include <tabela.h>
#include <hash.h>
#include <rn.h>
#include <assert.h>
#include <quad.h>
#include <delaunay.h>
#include <delaunay_plot_POV.h>


void analytic_inversion_plot_delaunay_POV (quad_arc_t e, delaunay_site_t *st, int nsites, double height[] ,FILE* wr,double scale,double offset);
void analytic_inversion_plot_delaunay_POV (quad_arc_t e, delaunay_site_t *st, int nsites, double height[] ,FILE* wr,double scale,double offset)
  { 
    float minX,minY,maxX,maxY;
    minX = minY = -1;
    maxX = maxY = +1;
    float maxZ = offset;
    float minZ = offset;
    int i;
    for(i = 0; i < nsites; i++){
	maxZ = fmax(maxZ,height[i]);
	minZ = fmin(minZ,height[i]);
    }
    fprintf(wr,"#declare minX = %f;\n",minX);
    fprintf(wr,"#declare maxX = %f;\n",maxX);
    fprintf(wr,"#declare minY = %f;\n",minY);
    fprintf(wr,"#declare maxY = %f;\n",maxY);
    fprintf(wr,"#declare maxZ = %f;\n",maxZ);
    fprintf(wr,"#declare minZ = %f;\n",minZ);
    fprintf(wr,"#declare zeroLevel = %f;\n",offset);
    fprintf(wr,"#declare graph = \n");
    fprintf(wr,"  mesh{\n");
    delaunay_plot_POV_triangles(wr, e,height,"graph_tx_out");
    fprintf(wr,"  }\n");
    fprintf(wr,"#declare skirt = \n");
    fprintf(wr,"  mesh{\n");
    delaunay_plot_POV_skirt(wr,e,height,"skirt_tx_out");
    fprintf(wr,"  }\n");
  }
  
  void analytic_inversion_generatePOVSceneTri(FILE* arq,double** X ,double* F,int n,double scale,double offset);
  void analytic_inversion_generatePOVSceneTri(FILE* arq,double** X ,double* F,int n,double scale,double offset){
	//the main idea is use an workspace with 2 decimal digits of precision, so we multiply all value by escala
	
	int nsites = n;
	delaunay_site_t *st = notnull(malloc(nsites*sizeof(delaunay_site_t)), "no mem");
	int i;
	for (i = 0; i < nsites; i++) {
	    st[i].index = i;
      	    st[i].p.c[0] = X[i][0];
            st[i].p.c[1] = X[i][1];
        }
	quad_arc_t e = delaunay_build (st, nsites);
	analytic_inversion_plot_delaunay_POV(e, st, nsites, F, arq,scale,offset);
	free(st);
}


analytic_inversion_t* analytic_inversion_create_model
( 
 int m,
 int type_gauge,
 approx_model_t* gauge_model,
 void** gauge_data,
 r3_t view_dir,
 double thetaMax,
 int resolution,
 int type_inv,
 approx_model_t* inv_model,
 void** inv_data,
 void* mag_data,
 char* plot_pov_prefix,
 char* plot_txt_prefix
)
{
  analytic_inversion_t* ai = (analytic_inversion_t*) malloc(sizeof(analytic_inversion_t));
  ai->type_gauge = type_gauge;
  ai->m = m;
  ai->gauge_model = gauge_model;
  ai->gauge_data = gauge_data;
  
  Tabela* tab = cria_tabela_from_model(m,resolution,view_dir,thetaMax,gauge_model,gauge_data);
  ai->u = (double*)malloc(sizeof(double)*m);
  ai->v = (double*)malloc(sizeof(double)*m);
  ai->b = (double*)malloc(sizeof(double)*m);
  calcula_sistema_de_coordenadas(tab, ai->u, ai->v, ai->b, &(ai->bu), &(ai->bv));
  int n_linhas = get_num_linhas(tab);
  
  
  
  double** X = (double**)malloc(sizeof(double*)*n_linhas);
  double* F =  (double*)malloc(sizeof(double)*n_linhas);
  int i;
  double radius = 0;
  
  for(i = 0; i < n_linhas; i++){
    X[i] = (double*)malloc(sizeof(double)*2);
    double* g = (double*) get_intdir(tab,i);
    X[i][0] = rn_dot(m,g,ai->u) - ai->bu;
    X[i][1] = rn_dot(m,g,ai->v) - ai->bv;
    double dist_center = rn_norm(2,X[i]);
    if(dist_center > radius) radius = dist_center;
  }
  assert(radius > 0.0000001);
  for(i = 0; i < n_linhas; i++){
    X[i][0]/= radius;
    X[i][1]/= radius;
  }
  ai->radius = radius;
  ai->type_inv = type_inv;
  ai->inv_model = inv_model;
  
  
  int k;
  for(k = 0; k < 3 ; k++){
    
    
    for(i = 0; i < n_linhas; i++){
      r3_t normal = get_normal(tab,i);
      F[i] = normal.c[k];
    }
    if(plot_txt_prefix != NULL){
      char*ref_name = NULL;
      char *ref_name = jsprintf("%s_InvData_%d.txt",plot_txt_prefix,k);
      FILE* arq_debug = open_write(ref_name,TRUE);
      for(i = 0; i < n_linhas; i++){
	fprintf(arq_debug,"%15.12lf %15.12lf %15.12lf\n",X[i][0],X[i][1],F[i]);
      }
      free(ref_name);
      fclose(arq_debug);
    } 
    
    if(plot_pov_prefix != NULL){
      char*ref_name = NULL;
      char *ref_name = jsprintf("%s_InvData_%d.inc",plot_pov_prefix,k);
      FILE* arq_debug = open_write(ref_name,TRUE);
      analytic_inversion_generatePOVSceneTri(arq_debug,X,F,n_linhas,1.0,0.0);
      free(ref_name);
      fclose(arq_debug);
    } 
    
    
    assert(inv_data != NULL);
    approx_model_fit_linear_model(X, F,NULL,NULL, n_linhas, inv_model, inv_data[k]);
    ai->inv_data[k] = inv_data[k];
  }
  
  /*Now the albedo*/
  for(i = 0; i < n_linhas; i++){
      r3_t normal = get_normal(tab,i);
      double G[m];
      int l;
      for(l = 0; l < m; l++){ G[l] = gauge_model->evaluate(normal.c,gauge_data[l]);}
      double Gmag = rn_norm(m,G);
      F[i] = Gmag;
  }
  if(plot_txt_prefix != NULL){
      char*ref_name = NULL;
      char *ref_name = jsprintf("%s_MagData_%d.txt",plot_txt_prefix,0);
      FILE* arq_debug = open_write(ref_name,TRUE);
      for(i = 0; i < n_linhas; i++){
	fprintf(arq_debug,"%15.12lf %15.12lf %15.12lf\n",X[i][0],X[i][1],F[i]);
      }
      free(ref_name);
      fclose(arq_debug);
  } 
  assert(mag_data != NULL);
  approx_model_fit_linear_model(X, F,NULL,NULL, n_linhas, inv_model, mag_data);
  ai->mag_data = mag_data;
  
  return ai;
  
}
/* Returns a inverse model of a set of ${m} gauges adjusted to the hemisphere with center {view_dir} and angular radius ${thetaMax},
using {nSamples} sampling normals */

void analytic_inversion_compute_gauge_OV(analytic_inversion_t* ai, r3_t normal, double G[])
{
  int m = ai->m;
  int i;
  for(i = 0; i < m; i++){
    G[i] = ai->gauge_model->evaluate(normal.c,ai->gauge_data[i]);
  }
}

void analytic_inversion_compute_normal(double* S, analytic_inversion_t* ai,r3_t* normal, double* albedo, double G[], double* dist)
{
  int m = ai->m;
  double s[m];
  double smag = rn_dir(m,S,s);
  
   
  double z[2];
  z[0] = (rn_dot(m,s,ai->u) - ai->bu)/ai->radius;
  z[1] = (rn_dot(m,s,ai->v) - ai->bv)/ai->radius;

  r3_t normal_v;
  int i;
  for(i = 0; i < 3 ; i++){
    normal_v.c[i] = ai->inv_model->evaluate(z,ai->inv_data[i]);
  }
  double mag  = r3_dir(&normal_v,&normal_v);
  
  if(isnan(mag)){
    for(i = 0; i < 3 ; i++){
      normal_v.c[i] = ai->inv_model->evaluate(z,ai->inv_data[i]);
    } 
  }
  assert(mag > 0.01);
//   analytic_inversion_compute_gauge_OV(ai,normal_v,G);
//   double g[m];
//   double gmag = rn_dir(m,G,g);
  double gmag = ai->inv_model->evaluate(z,ai->mag_data);
  if(gmag > 0.0){
  *normal = normal_v;
  *albedo = smag/gmag;
  }else{
    *normal = normal_v;
    *albedo = 0;
  }
//   *dist = rn_dist(m,s,g);  
  *dist = 0.5;
  
}


void analytic_inversion_write(FILE* arq,analytic_inversion_t* ai){
  fprintf(arq,"%d %d\n",ai->m,ai->type_gauge);
  int m = ai->m;
  int i;
  for(i = 0; i < m; i++){
    ai->gauge_model->write_parameters(arq,ai->gauge_data[i]);
  }
  for(i = 0; i < m; i++){
    fprintf(arq,"%9.6lf ",ai->u[i]);
  }
  fprintf(arq,"\n");
  for(i = 0; i < m; i++){
    fprintf(arq,"%9.6lf ",ai->v[i]);
  }
  fprintf(arq,"\n");
  for(i = 0; i < m; i++){
    fprintf(arq,"%9.6lf ",ai->b[i]);
  }
  fprintf(arq,"\n");
  fprintf(arq,"%9.6lf %9.6lf %9.6lf\n",ai->bu,ai->bv,ai->radius);
  
  fprintf(arq,"%d \n",ai->type_inv);
  for(i = 0; i < 3; i++){
    ai->inv_model->write_parameters(arq,ai->inv_data[i]);
  }
  ai->inv_model->write_parameters(arq,ai->mag_data);

  
}

analytic_inversion_t* analytic_inversion_read(FILE* arq){
  analytic_inversion_t* ai = (analytic_inversion_t*)malloc(sizeof(analytic_inversion_t));
  fscanf(arq,"%d %d",&(ai->m),&(ai->type_gauge));
  int m = ai->m;
  int i;
  ai->gauge_model = create_approx_lighting_model(ai->type_gauge);
  ai->gauge_data = (void*)malloc(sizeof(void*)*m);
  for(i = 0; i < m; i++){
    fprintf(stderr,"Gauge %02d\n",i);
    ai->gauge_data[i] = ai->gauge_model->read_parameters(arq);
//     ai->gauge_model->write_parameters(stderr,ai->gauge_data[i]);
  }
  
  ai->u = (double*)malloc(sizeof(double)*m);
  for(i = 0; i < m; i++){
    fscanf(arq,"%lf",&(ai->u[i]));
  }
  
  ai->v = (double*)malloc(sizeof(double)*m);
  for(i = 0; i < m; i++){
    fscanf(arq,"%lf",&(ai->v[i]));
  }
  
  ai->b = (double*)malloc(sizeof(double)*m);
  for(i = 0; i < m; i++){
    fscanf(arq,"%lf",&(ai->b[i]));
  }
  
  fscanf(arq,"%lf %lf %lf",&(ai->bu),&(ai->bv),&(ai->radius));
  fscanf(arq,"%d",&(ai->type_inv));
  if(ai->type_inv == GENERICPOLY_MODEL){
    ai->inv_model = create_ls_polynomial_model();
  }
  else{
    fprintf(stderr,"analytic_inversion_read - ERROR : invalid model for 2D inversion\n");
    assert(FALSE);
  }
  
  for(i = 0; i < 3; i++){
      ai->inv_data[i] = ai->inv_model->read_parameters(arq);
//       ai->inv_model->write_parameters(stderr,ai->inv_data[i]);
  }
  ai->mag_data = ai->inv_model->read_parameters(arq);
//   ai->inv_model->write_parameters(stderr,ai->mag_data);
  
  
  
  return ai;
}

float_image_t* analytic_inversion_plot_map(analytic_inversion_t* ai,int imSize){
  assert(imSize > 0);
  float_image_t* im = float_image_new(4,imSize,imSize);
  int X,Y;
  double imRad = imSize/2.0;
  for(X = 0; X < imSize; X++){
    for(Y = 0; Y < imSize; Y++){
     
      double coords[2];
      coords[0] = (X-imRad)/imRad;
      coords[1] = (Y-imRad)/imRad;
      bool_t test_inside = hypot(coords[0],coords[1]) <=1;
      double func_val[4];
      int c;
      for(c = 0; c < 4; c++){
	void* inv_data_used = (c < 3 ? ai->inv_data[c] : ai->mag_data);
	func_val[c] = (test_inside ? ai->inv_model->evaluate(coords,inv_data_used) : 0);
	float_image_set_sample(im,c,X,Y,func_val[c]);
      }
      
    }
  }
  
  return im;
}


void analytic_inversion_generate_UV_plot(char* prefix, analytic_inversion_t* ai, long int num_samples_axis,r3_t view_dir){
  double theta,phi;
  r3_t normal;
  char* filename = NULL;
  char *filename = jsprintf("%s-UVplot.fni",prefix);
  FILE* arq = open_write(filename, TRUE);
  double theta_interval = ((M_PI)/2.0)/100.0;
  double phi_interval = ((M_PI)*2.0)/100.0;
  for(theta = 0; theta < (M_PI)/2.0 ; theta+= theta_interval){
    for(phi = 0; phi < (M_PI*2.0); phi+=phi_interval){
   
//       fprintf(stderr,"%9.6lf %9.6lf\n",theta,phi);
      double x = sin(theta)*cos(phi);
      double y = sin(theta)*sin(phi);
      double z = cos(theta);
      normal = (r3_t){{x,y,z}};
      double G[ai->m];
      int i;
      for(i = 0; i < ai->m; i++){
	G[i] = ai->gauge_model->evaluate(normal.c,ai->gauge_data[i]);
      }
      double g[ai->m];
      rn_dir(ai->m,G,g); 
      
      double Z[2];
      Z[0] = (rn_dot(ai->m,g,ai->u) - ai->bu)/ai->radius;
      Z[1] = (rn_dot(ai->m,g,ai->v) - ai->bv)/ai->radius;
      
      double G_aux[ai->m];
      double dist,albedo;
      r3_t normal_approx;
      analytic_inversion_compute_normal(G, ai,&normal_approx,&albedo, G_aux,&dist);
      
      double func_val = acos(r3_dot(&normal_approx,&normal));
      
      fprintf(arq,"%9.6lf %9.6lf %9.6lf %9.6lf %9.6lf %9.6lf \n",Z[0],Z[1],x,y,z,func_val);
     
    }
  }
  
  fclose(arq);
  free(filename);
  filename=NULL;
  
}

void   analytic_inversion_generate_SG_plot(char* prefix, analytic_inversion_t* ai, long int num_samples_axis,r3_t view_dir){
  long int X,Y;
  r3_t normal;
  double Smax = 1.0;
  double thetaMAX = M_PI/2.0;
  float_image_t* plot_im = float_image_new(1,num_samples_axis,num_samples_axis);
  for(X = 0; X < num_samples_axis;X++){
    for(Y = 0; Y < num_samples_axis; Y++){
      double radius = num_samples_axis/2.0;
      double dX = Smax*(X - radius)/radius;
      double dY = Smax*(Y - radius)/radius;
      double func_val = 0;
      double x = 2*dX/((dX*dX) + (dY*dY) + 1);
      double y = 2*dY/((dX*dX) + (dY*dY) + 1);
      double z = (1 - (dX*dX) - (dY*dY))/((dX*dX) + (dY*dY) + 1);
      normal = (r3_t){{x,y,z}};
      bool_t valid_hemisphere = (r3_dot(&normal,&view_dir) > cos(thetaMAX) ? TRUE : FALSE);
      if(!valid_hemisphere) {
	func_val = 0;
      }else{
	int i;
	double G[ai->m];
	for(i = 0; i < ai->m; i++){
	  G[i] = ai->gauge_model->evaluate(normal.c,ai->gauge_data[i]);
	}

	
	
	r3_t normal_approx;
	double G_aux[ai->m];
	double dist,albedo;
	analytic_inversion_compute_normal(G, ai,&normal_approx,&albedo, G_aux,&dist);
// 	func_val = r3_dist(&normal_approx,&normal);
	func_val = acos(r3_dot(&normal_approx,&normal));
// 	fprintf(stderr,"Acos");
      }
      float_image_set_sample(plot_im,0,X,Y,func_val);
    }
  }
  char* filename = NULL;
  char *filename = jsprintf("%s-ApproxError.fni",prefix);
  FILE* arq = open_write(filename,TRUE);
  float_image_write(arq,plot_im);
  fclose(arq);
  float_image_free(plot_im);
  free(filename);
  filename=NULL;
  
}
