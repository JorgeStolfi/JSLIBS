#include <fast_hash.h>
#include <hash.h>
#include <tabela.h>
#include <polynomial_functions.h>
#include <assert.h>

r2_t fasthash_mapeiaHash
  ( const double ass[], 
    double u[], 
    double v[], 
    const double baricentro[],
    int num_luzes,
    int tam_grid,
    double R
    );

r2_t fasthash_mapeiaHash
  ( const double ass[], 
    double u[], 
    double v[], 
    const double baricentro[],
    int num_luzes,
    int tam_grid,
    double R
    )
{
  int N = tam_grid;
  r2_t indice;
  double s[num_luzes];
  int i;
  for (i = 0; i < num_luzes;i++) s[i] = ass[i] - baricentro[i];
  double su = rn_dot(num_luzes, s,u);
  double sv = rn_dot(num_luzes, s,v);
  assert(fabs(su) <= 1.0001); 
  assert(fabs(sv) <= 1.0001); 
  indice.c[0] = (su/R + 1)*(N/2.0);
 // if(indice.c[0] < 0) { indice.c[0] = 0; }
 // if(indice.c[0] > (N-1)) { indice.c[0]= N-1; }
  indice.c[1] = (sv/R + 1)*(N/2.0);
 // if(indice.c[1] < 0) i{ ndice.c[1] = 0; }
 // if(indice.c[1] > (N-1)) { indice.c[1]= N-1; }
  // fprintf(stderr,  "XR: %f %d YR: %f %d \n",su,indice.c[0],sv,indice.c[1]);
  return indice;
}




fast_hash_t* create_fasthash(int N,Tabela* tab,double sigma,int normal_degree, int albedo_degree){
  int nLights = get_num_luzes(tab);
  fast_hash_t* fh = (fast_hash_t*)malloc(sizeof(fast_hash_t));
  fh->N = N;
  fh->nLights = nLights;
  fh->baricenter = (double*)malloc(sizeof(double)*nLights);
  fh->u = (double*)malloc(sizeof(double)*nLights);
  fh->v = (double*)malloc(sizeof(double)*nLights);
  
  
  double* u = fh->u;
  double* v = fh->v;
  double* bar = fh->baricenter;
  double R ;
  calcula_sistema_de_coordenadas(tab, u, v, bar, &fh->bu, &fh->bv);
  int num_lines = get_num_linhas(tab);
  
//   bucketGrid* bg = CriaBucketGrid(tab,N);
  
  /*Init mass of data*/
 
  
  int i;
  
  R = -1;
  for(i = 0; i < num_lines; i++){
    const double *go = get_intdir(tab,i);
    double dgo[nLights];
    int t;
    for (t = 0; t < nLights; t++) dgo[t] = go[t] -  bar[t];
    double su = fabs(rn_dot(nLights,dgo,u));
    double sv = fabs(rn_dot(nLights,dgo,v));
    // fprintf(stderr,  "linha %5d SU = %f SV = %f\n",i,su,sv);
    if(R < su) R = su;
    if(R < sv) R = sv;
  }
  
  fh->R = R;
  
  double** OX = (double**)malloc(sizeof(double*)*num_lines);
  double** X = (double**)malloc(sizeof(double*)*num_lines);
  
  double* OFx = (double*)malloc(sizeof(double)*num_lines);
  double* OFy = (double*)malloc(sizeof(double)*num_lines);
  double* OFz = (double*)malloc(sizeof(double)*num_lines);
  double* OFa = (double*)malloc(sizeof(double)*num_lines);
  
  double* Fx = (double*)malloc(sizeof(double)*num_lines);
  double* Fy = (double*)malloc(sizeof(double)*num_lines);
  double* Fz = (double*)malloc(sizeof(double)*num_lines);
  double* Fa = (double*)malloc(sizeof(double)*num_lines);
  
  
  double* w = (double*)malloc(sizeof(double)*num_lines);
  double* wpos = (double*)malloc(sizeof(double)*num_lines);
  double* wcopy = (double*)malloc(sizeof(double)*num_lines);
    
  for(i = 0 ; i < num_lines;i++){
    
    r3_t normal = get_normal(tab,i);
    const double* go = get_intdir(tab,i);
    r2_t mapped_cell= fasthash_mapeiaHash(go,u,v,bar,nLights,N,R);
    OX[i] = (double*)malloc(sizeof(double)*2);
    X[i] = (double*)malloc(sizeof(double)*2);
    
    OX[i][0] = mapped_cell.c[0];
    OX[i][1] = mapped_cell.c[1];
    OFx[i] = normal.c[0];
    
//     OFx[i] = normal.c[0]/(1 + normal.c[2]);
    
    OFy[i] = normal.c[1];
//     OFy[i] = normal.c[1]/(1 + normal.c[2]);
    
    OFz[i] = normal.c[2];
    OFa[i] = get_intmag(tab,i);
    wpos[i] = 1.0;
  }
  
  approx_model_t* lm = create_ls_polynomial_model();
  
  int count = 0;
  
  fprintf(stderr,"\n\n");
  int iu,iv;
  
  fh->hash_table = (r3_t**)malloc(sizeof(r3_t*)*N);
  fh->albedo_table = (double**)malloc(sizeof(double*)*N);
  fh->weight_table = (double**)malloc(sizeof(double*)*N);
  for(iu = 0; iu < N; iu++){
    fh->hash_table[iu] = (r3_t*)malloc(sizeof(r3_t)*N);
    fh->albedo_table[iu] = (double*)malloc(sizeof(double)*N);
    fh->weight_table[iu] = (double*)malloc(sizeof(double)*N);
    
    for(iv = 0; iv < N; iv++){
      
      count++;
      
      double PESO_LIMITE = 10e-8;
      int l;
 //     double sigma = 0.5;
      int nComp = 0;
      double sumW = 0;
      
      for(l = 0; l < num_lines; l++){
	//double dist2u = (iu - X[l][0])*(iu - X[l][0]);
	double dist2u = (iu  - OX[l][0])/(double)N;
	dist2u = dist2u*dist2u;
	
	double dist2v = (iv  - OX[l][1])/(double)N;
	dist2v = dist2v*dist2v;
	
// 	double dist2v = (iv - X[l][1])*(iv - X[l][1]);
	double peso = exp(-(dist2u+dist2v)/(2.0*(sigma*sigma)));
	
// 	if(sqrt(dist2u + dist2v) < 0.02){
// 	  if(bg->buckets[iu][iv].itens != NULL){
// 	    double* centro = bg->buckets[iu][iv].centro;
// 	    
// 	    double dist = rn_dist(nLights,centro,
// 	  }
// 	}
	
	if(peso > PESO_LIMITE){
	  //w[l]= 
	  w[nComp] = peso;
	  X[nComp][0] = OX[l][0];
	  X[nComp][1] = OX[l][1];
	  Fx[nComp] = OFx[l];
	  Fy[nComp] = OFy[l];
	  Fz[nComp] = OFz[l];
	  Fa[nComp] = OFa[l];
	  sumW+=peso;
	  nComp++;
	}
	
	
	
      }
      
      if(nComp > 0){
	sumW = sumW/(double)nComp;
      }
      fh->weight_table[iu][iv] = sumW;
      
      fprintf(stderr,"\033[1A");
      fprintf(stderr,"[%d][%d] - %d components-  %6.6f%% \n",iu,iv,nComp,count*100/(float)(N*N));
      
     
      
      
      
      
      poly_function_t* l_data_x = poly_init_components(2, normal_degree,TRUE);
      poly_function_t* l_data_y = poly_init_components(2, normal_degree,TRUE);
      poly_function_t* l_data_z = poly_init_components(2, normal_degree,TRUE);
      poly_function_t* l_data_a = poly_init_components(2, albedo_degree,TRUE);
     
            
      rn_copy(nComp,w,wcopy);
      approx_model_fit_linear_model(X,Fx,wcopy,wpos,nComp,lm,l_data_x);
     
      
//       fitModelToFunction(X,Fx, nComp,lm,l_data_x,1);
      
      rn_copy(num_lines,w,wcopy);
//       fitModelToFunction(X,Fy, nComp,lm,l_data_y,1);
      approx_model_fit_linear_model(X,Fy,wcopy,wpos,nComp,lm,l_data_y);
      
      rn_copy(num_lines,w,wcopy);
//       lm->weights = wcopy;
   //   lm->wpos = NULL;
//      fitModelToFunction(X,Fz, nComp,lm,l_data_z,1);
       approx_model_fit_linear_model(X,Fz,wcopy,wpos,nComp,lm,l_data_z);
      
      rn_copy(num_lines,w,wcopy);
//       lm->weights = wcopy;
   //   lm->wpos = NULL;
//       fitModelToFunction(X,Fa, nComp,lm,l_data_a,1);
    approx_model_fit_linear_model(X,Fa,wcopy,wpos,nComp,lm,l_data_a);

      r3_t normal_res;
      double pt_coord[] = {iu,iv};
      normal_res.c[0] = lm->evaluate(pt_coord,l_data_x);
      
      normal_res.c[1] = lm->evaluate(pt_coord,l_data_y);
      
      
      normal_res.c[2] = lm->evaluate(pt_coord,l_data_z);
      double albedo = lm->evaluate(pt_coord,l_data_a);
      lm->release_data(l_data_x);
      lm->release_data(l_data_y);
      lm->release_data(l_data_z);
      lm->release_data(l_data_a);
	  
      if(isnan(albedo)){
	albedo = 0;
	 fh->weight_table[iu][iv] = 0;
      }
      
      
      
      r3_dir(&normal_res,&normal_res);
      
      if(isnan(normal_res.c[0]) && isnan(normal_res.c[1]) && isnan(normal_res.c[2])){
	r3_zero(&normal_res);
      }
      fh->hash_table[iu][iv] = normal_res;
      fh->albedo_table[iu][iv] = albedo;
   }
  }
  
  free(Fx);
  free(Fy);
  free(Fz);
  free(Fa);
  
  free(OFx);
  free(OFy);
  free(OFz);
  free(OFa);
  
  free(w);
  free(wcopy);
  free(wpos);
  
  return fh;
}


void SaveFastHash(FILE* arq, fast_hash_t* fh){
  fprintf(arq,"%d\n", fh->N); 
  fprintf(arq,"%d\n",fh->nLights);
 
  int i;
  
  for(i = 0; i < fh->nLights; i++){
    fprintf(arq,"%8.6lf ",fh->u[i]);
  }
  fprintf(arq,"\n");
  
  
  for(i = 0; i < fh->nLights; i++){
    fprintf(arq,"%8.6lf ",fh->v[i]);
  }
  fprintf(arq,"\n");
  
  fprintf(arq,"%8.6lf %8.6lf %8.6lf\n",fh->bu,fh->bv,fh->R);

  for(i = 0; i < fh->nLights; i++){
    fprintf(arq,"%8.6lf ",fh->baricenter[i]);
  }
  fprintf(arq,"\n");
  
  int j;
  for(i = 0; i < fh->N; i++){
    for(j = 0; j < fh->N; j++){
      fprintf(arq,"%8.6lf %8.6lf %8.6lf\n",fh->hash_table[i][j].c[0],fh->hash_table[i][j].c[1],fh->hash_table[i][j].c[2]);
    }
  }
  
  for(i = 0; i < fh->N; i++){
    for(j = 0; j < fh->N; j++){
      fprintf(arq,"%8.6lf\n",fh->albedo_table[i][j]);
    }
  }
  
  for(i = 0; i < fh->N; i++){
    for(j = 0; j < fh->N; j++){
      fprintf(arq,"%8.6lf\n",fh->weight_table[i][j]);
    }
  }
  
 
}



void PrintFastHash(FILE* arq, fast_hash_t* fh){
  fprintf(arq,"********************************\n");
  fprintf(arq,"Hash Size: %d\n", fh->N); 
  fprintf(arq,"Number of lights: %d\n",fh->nLights);
 
  rn_gen_print(arq,fh->nLights,fh->u,"%8.6lf","U (",",",")\n");
  rn_gen_print(arq,fh->nLights,fh->v,"%8.6lf","V (",",",")\n");
  fprintf(arq,"BU: %8.6lf\nBV: %8.6lf\nR:%8.6lf\n",fh->bu,fh->bv,fh->R);
  rn_gen_print(arq,fh->nLights,fh->baricenter,"%8.6lf","BAR (",",",")\n");
  fprintf(arq,"********************************\n");

 
}



float_image_t* FastHashNormalsToFNI(fast_hash_t* fh){
  int N = fh->N;
  float_image_t* im = float_image_new(3,N,N);
  int x,y;
  for(x = 0; x < N; x++){
    for(y = 0; y < N; y++){
	r3_t normal = fh->hash_table[x][y];
	int c;
	for(c = 0; c < 3; c++){
	  float_image_set_sample(im,c,x,y,normal.c[c]);
	}
    }
  }
  
  return im;
  
}

float_image_t* FastHashAlbedoToFNI(fast_hash_t* fh){
  int N = fh->N;
  float_image_t* im = float_image_new(1,N,N);
  int x,y;
  for(x = 0; x < N; x++){
    for(y = 0; y < N; y++){
	float_image_set_sample(im,0,x,y,fh->albedo_table[x][y]);
    }
  }
  
  return im;
  
}


float_image_t* FastHashWeightsToFNI(fast_hash_t* fh){
  int N = fh->N;
  float_image_t* im = float_image_new(1,N,N);
  int x,y;
  for(x = 0; x < N; x++){
    for(y = 0; y < N; y++){
	float_image_set_sample(im,0,x,y,fh->weight_table[x][y]);
    }
  }
  
  return im;
  
}

fast_hash_t* LoadFastHash(FILE* arq){
  
  fast_hash_t* fh = (fast_hash_t*)malloc(sizeof(fast_hash_t));
  fscanf(arq,"%d", &(fh->N)); 
  fscanf(arq,"%d",&(fh->nLights));
 
  int i;
  fh->u = (double*)malloc(sizeof(double)*(fh->nLights));  
  for(i = 0; i < fh->nLights; i++){
    fscanf(arq,"%lf",&(fh->u[i]));
  }
  
  fh->v = (double*)malloc(sizeof(double)*(fh->nLights));
  for(i = 0; i < fh->nLights; i++){
    fscanf(arq,"%lf",&(fh->v[i]));
  }
    
  fscanf(arq,"%lf %lf %lf",&(fh->bu),&(fh->bv),&(fh->R));

  fh->baricenter = (double*)malloc(sizeof(double)*(fh->nLights));
  for(i = 0; i < fh->nLights; i++){
    fscanf(arq,"%lf",&(fh->baricenter[i]));
  }
  
  int j;
  fh->hash_table = (r3_t**)malloc(sizeof(r3_t*)*(fh->N));
  for(i = 0; i < fh->N; i++){
    fh->hash_table[i] = (r3_t*)malloc(sizeof(r3_t)*(fh->N));
    for(j = 0; j < fh->N; j++){
      fscanf(arq,"%lf %lf %lf",&(fh->hash_table[i][j].c[0]),&(fh->hash_table[i][j].c[1]),&(fh->hash_table[i][j].c[2]));
    }
  }
  
  fh->albedo_table = (double**)malloc(sizeof(double*)*(fh->N));
  for(i = 0; i < fh->N; i++){
    fh->albedo_table[i] = (double*)malloc(sizeof(double)*(fh->N));
    for(j = 0; j < fh->N; j++){
      fscanf(arq,"%lf",&(fh->albedo_table[i][j]));
    }
  } 
  
  fh->weight_table = (double**)malloc(sizeof(double*)*(fh->N));
  for(i = 0; i < fh->N; i++){
    fh->weight_table[i] = (double*)malloc(sizeof(double)*(fh->N));
    for(j = 0; j < fh->N; j++){
      fscanf(arq,"%lf",&(fh->weight_table[i][j]));
    }
  } 
  
  return fh;
}


r3_t fast_hash_compute_normal( fast_hash_t* fh, const double SO[], double sigma,double omg0,double omg1, double *albedo,double* logProb ){
  
  int nLights = fh->nLights;
  double so[nLights];
  double Smag = rn_dir(nLights,(double*)SO,so);
  r2_t pt_coords = fasthash_mapeiaHash(so,fh->u,fh->v,fh->baricenter,nLights,fh->N,fh->R);
    
  double u = pt_coords.c[0];
  double v = pt_coords.c[1];
  
  if( (u < 1) || (u >= fh->N ) ){
    *albedo = 0;
    *logProb = 0;
    return (r3_t){{0,0,0}};
  }
  
  r3_t normal_res;
 
//   int u0 = floor(pt_coords.c[0]);
//   int u1 = ceil(pt_coords.c[0]);
//   
//   int v0 = floor(pt_coords.c[1]);
//   int v1 = ceil(pt_coords.c[1]);
//   
//   r3_t normal_00 = fh->hash_table[u0][v0];
//   r3_t normal_01 = fh->hash_table[u0][v1];
//   r3_t normal_10 = fh->hash_table[u1][v0];
//   r3_t normal_11 = fh->hash_table[u1][v1];
//   
//   r3_t R1,R2;
//   double qu0 = (u1 - u)/(float)(u1 - u0);
//   double qu1 = (u - u0)/(float)(u1 - u0);
//   double qv0 = (v1 - v)/(float)(v1 - v0);
//   double qv1 = (v - v0)/(float)(v1 - v0);
//   
//   r3_mix(qu0, &normal_00,qu1,&normal_10, &R1);
//   r3_mix(qu0, &normal_01,qu1,&normal_11, &R2);
//   
//   r3_mix(qv0, &R1,qv1,&R2, &normal_res);
//   
//   r3_dir(&normal_res,&normal_res);
//   
//   /* Sets {r := s * a + t * b}. */ 
//   double a00 = fh->albedo_table[u0][v0];
//   double a01 = fh->albedo_table[u0][v1];
//   double a10 = fh->albedo_table[u1][v0];
//   double a11 = fh->albedo_table[u1][v1];
//   
//   double ar1 = (qu0*a00) + (qu1*a10) ;
//   double ar2 = (qu0*a01) + (qu1*a11) ;
//   
//   double Gmag =   (qv0*ar1) + (qv1*ar2);
//   *albedo =  Smag/(Gmag + 10e-16);
//   
//   double w00 = fh->albedo_table[u0][v0];
//   double w01 = fh->albedo_table[u0][v1];
//   double w10 = fh->albedo_table[u1][v0];
//   double w11 = fh->albedo_table[u1][v1];
//   
//   double wr1 = (qu0*w00) + (qu1*w10) ;
//   double wr2 = (qu0*w01) + (qu1*w11) ;
//   double weight = (qv0*wr1) + (qv1*wr2);
  
  int uu,vv;
  uu = u;
  vv = v;
  normal_res = fh->hash_table[uu][vv];
  double Gmag =   fh->albedo_table[uu][vv];
  double weight = fh->weight_table[uu][vv];
  
  *albedo = Smag/(Gmag + 10e-16);
  
  *logProb = 2.0*log(Smag)*weight;
  
  return normal_res;
}
