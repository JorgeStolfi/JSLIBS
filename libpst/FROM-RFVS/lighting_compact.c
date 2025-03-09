
#include <assert.h>

#include <math_utils.h>
#include <jsfile.h>
#include <lighting_models.h>
#include <r3x3.h>
#include <rn.h>
#include <pst_lamp.h>
#include <lighting_models.h>
#include <lighting_compact.h>

lighting_compact_options_t* lighting_compact_parse(argparser_t* pp){
  lighting_compact_options_t* o = (lighting_compact_options_t*)malloc(sizeof(lighting_compact_options_t));
  fprintf(stderr,"  COMPACT_MODEL\n");
  argparser_get_keyword_next(pp, "dirlight");
  o->dirlight.c[0] = argparser_get_next_double(pp, -10000, 10000);
  o->dirlight.c[1] = argparser_get_next_double(pp, -10000, 10000);
  o->dirlight.c[2] = argparser_get_next_double(pp, -10000, 10000);
  r3_t normU;
  double lenghtU = r3_dir(&(o->dirlight),&normU);
  if(lenghtU != 1.0){
	fprintf(stderr,"Non normalized U (%9.6lf %9.6lf %9.6lf), normalized to (%9.6lf %9.6lf %9.6lf)\n",o->dirlight.c[0],o->dirlight.c[1],o->dirlight.c[2],normU.c[0],normU.c[1],normU.c[2]);
  }
  o->dirlight = normU;
  o->err_dirlight = lighting_error_parse(pp,M_PI/2.0);
  argparser_get_keyword_next(pp, "radlight");
  o->radlight = argparser_get_next_double(pp, 0.0, M_PI/2.0);
  o->err_radlight = lighting_error_parse(pp,o->radlight);
    

  if(argparser_keyword_present_next(pp, "shine")){
    /*INFINITY means a mirror, maybe someday...*/
    o->shine = argparser_get_next_double(pp, 1.0, 1e5);
    o->err_shine = lighting_error_parse(pp,o->shine - 1);
  }else{
    o->shine = INFINITY;
    o->err_shine = 0;
  }

  o->has_isotropic = argparser_keyword_present_next(pp,"isotropic");
  o->has_backlight = argparser_keyword_present_next(pp,"backlight");
  if(o->has_backlight){
    o->backnormal.c[0] = argparser_get_next_double(pp, -10000, 10000);
    o->backnormal.c[1] = argparser_get_next_double(pp, -10000, 10000);
    o->backnormal.c[2] = argparser_get_next_double(pp, -10000, 10000);
    r3_t normXi;
    double lenghtXi = r3_dir(&(o->backnormal),&normXi);
    if(lenghtXi != 1.0){
	  fprintf(stderr,"Non normalized Xi (%9.6lf %9.6lf %9.6lf), normalized to (%9.6lf %9.6lf %9.6lf)\n",o->backnormal.c[0],o->backnormal.c[1],o->backnormal.c[2],normXi.c[0],normXi.c[1],normXi.c[2]);
    }
    o->backnormal = normXi;   
    o->err_backnormal = lighting_error_parse(pp,M_PI/2.0);    
  }else{
    o->backnormal = (r3_t){{0,0,0}};
    o->err_backnormal = 0;
  }
  return o;
}

approx_model_t* lighting_compact_create_approx_lighting_model(void){
  approx_model_t* am = (approx_model_t*)malloc(sizeof(approx_model_t));
  am->type = COMPACT_MODEL;

  am->phi = lighting_compact_phi;
  am->get_num_components = lighting_compact_get_number_components;
  am->set_alphas = lighting_compact_set_alphas;
  am->get_alphas = lighting_compact_get_alphas;
  am->copy_data = lighting_compact_copy_lighting_data;
  am->release_data =  lighting_compact_release_lighting_data;
  
  /*Optional members - they wont be called by the fitting functions and can be NULL*/
  am->evaluate = lighting_compact_shading; 
  am->write_parameters = lighting_compact_write_parameters;
  am->read_parameters = lighting_compact_read_parameters; 
  /*for non linear models - NOT USED*/
  am->compare = NULL; 
  am->get_num_nl_parameters = lighting_compact_num_nl_parameters;
  am->pack_nl_parameters =  lighting_compact_pack_nl_parameters; 
  am->get_nl_centers_and_deltas = lighting_compact_get_nl_centers_and_deltas;
  am->unpack_nl_parameters = lighting_compact_unpack_nl_parameters;
  am->update_nl_parameters_and_errors = lighting_compact_update_params_and_errors;
  am->generate_nl_parameter_simplex = lighting_compact_generate_simplex; 
  return am;
}


/*Compact Light model*/

lighting_compact_data_t* lighting_compact_init_components(lighting_compact_options_t* o, r3_t view_dir){
  
  lighting_compact_data_t* cl = (lighting_compact_data_t*)malloc(sizeof(lighting_compact_data_t));
  cl->view_dir = view_dir;
  cl->dirlight = o->dirlight;
  cl->err_dirlight =  o->err_dirlight;
  
  if( r3_norm(&(o->dirlight)) == 0){
    cl->Hc = FALSE;
    cl->dirlight.c[0] = 0;
    cl->dirlight.c[1] = 0;
    cl->dirlight.c[2] = 1;
  }else{
      cl->Hc = TRUE;
  }
  
  cl->radlight = o->radlight;
  cl->err_radlight = o->err_radlight;
  
  
  cl->shine = o->shine;
  cl->err_shine = o->err_shine;
  if(o->shine == INFINITY){
    cl->Hl = FALSE;
  }else{
    cl->Hl = TRUE;
  }
  
  cl->backnormal = o->backnormal;
  cl->err_backnormal = o->err_backnormal;
  if( !o->has_backlight){
    cl->Hf = FALSE;
    cl->backnormal.c[0] = 0;
    cl->backnormal.c[1] = 0;
    cl->backnormal.c[2] = 1;
  }else{
    cl->Hf = TRUE;
  }
  
  cl->Hi = o->has_isotropic;
  cl->Ec = cl->El =  cl->Ei = cl->Ef = 0;
  return cl;
}



void lighting_compact_determine_basis_indices(void* l_data,int* ne, int* num_compac, int* num_phong, int* num_isotropic, int* num_backplane){
  /*
  Returns the effective number elements {ne} in the basis and the indices of each element if they are included in the basis, if
  it is not included sets its index to -1.
  */
  lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
  int n = 0;
  
  int num_c,num_p,num_i,num_b;
  num_c = num_p = num_i = num_b = -1;
  
  if(cl->Hc){num_c = n; n++;}
  if(cl->Hl){num_p = n; n++;}
  if(cl->Hi){num_i = n; n++;}
  if(cl->Hf){num_b = n; n++;}
  
  if(num_compac != NULL)    {*num_compac    = num_c;   }
  if(num_phong != NULL)     {*num_phong     = num_p;   } 
  if(num_isotropic != NULL) {*num_isotropic = num_i;   }
  if(num_backplane != NULL) {*num_backplane = num_b;   }
  
  *ne = n;
}



double lighting_compact_phi(int r, double* x,void* l_data){
  lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
  r3_t normal = (r3_t){{x[0],x[1],x[2]}};
  
  auto double Phong(void);
  double Phong(void){
    r3_t w = compute_bissetriz(cl->dirlight,cl->view_dir);
    double phong_base = fmax(0,r3_dot(&(cl->dirlight),&normal))*fmax(0, r3_dot(&w,&normal));
    return ( phong_base <= 0 ? 0 : pow(phong_base,cl->shine) ) ;
  }
  
  int n_comp;
  int num_compac,num_phong,num_isotropic,num_backplane;
  lighting_compact_determine_basis_indices(l_data,&n_comp,&num_compac,&num_phong,&num_isotropic,&num_backplane);
   
  if(num_compac == r){
    double val = pst_lamp_geom_factor(&normal, &(cl->dirlight),cos(cl->radlight));
//     val = r3_dot(&normal,&(cl->dirlight));
    return val;
  }else if (num_phong == r){
    return Phong();
  }else if (num_isotropic == r){
    return 1;
  }else if( num_backplane == r){
    return 1 - r3_dot(&normal,&(cl->backnormal));
  }
  
  fprintf(stderr,"libraab - lighting_compact_phi: ERROR - invalid value for r = %d within range [0-%d]\n",r,n_comp -1 );
  assert(FALSE);
 
}


int lighting_compact_get_number_components(void* l_data){
  int n_comp;
  lighting_compact_determine_basis_indices(l_data,&n_comp,NULL,NULL,NULL,NULL);
  return n_comp;
}


void lighting_compact_set_alphas(double* C,int n, void* l_data){
  lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
  int n_comp;
  int num_compac,num_phong,num_isotropic,num_backplane;
  lighting_compact_determine_basis_indices(l_data,&n_comp,&num_compac,&num_phong,&num_isotropic,&num_backplane);
  assert(n_comp == n);
  cl->Ec = (num_compac    >= 0 ?   C[num_compac]    : 0);
  cl->El = (num_phong     >= 0 ?   C[num_phong]     : 0);
  cl->Ei = (num_isotropic >= 0 ?   C[num_isotropic] : 0);
  cl->Ef = (num_backplane >= 0 ?   C[num_backplane] : 0);

}

void lighting_compact_get_alphas(void* l_data,double* C,int n){
  lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
  int n_comp;
  int num_compac,num_phong,num_isotropic,num_backplane;
  lighting_compact_determine_basis_indices(l_data,&n_comp,&num_compac,&num_phong,&num_isotropic,&num_backplane);
  assert(n_comp == n);
  if( num_compac >= 0 ) C[num_compac] = cl->Ec;
  if( num_phong  >= 0 ) C[num_phong]  = cl->El;
  if( num_isotropic >= 0 ) C[num_isotropic] = cl->Ei;
  if( num_backplane >= 0 ) C[num_backplane] = cl->Ef;
}

void* lighting_compact_copy_lighting_data(void* l_data){
  lighting_compact_data_t* cl1 = (lighting_compact_data_t*)l_data;
  lighting_compact_data_t* cl2 =  (lighting_compact_data_t*) malloc(sizeof(lighting_compact_data_t));
  *cl2 = *cl1;
  
  return cl2;
}

void lighting_compact_release_lighting_data(void* l_data){
  free(l_data);
}

double lighting_compact_shading(double* x,void* l_data){
  int n_comp = lighting_compact_get_number_components(l_data);
  int i;
  double val = 0;
  double alphas[n_comp];
  lighting_compact_get_alphas(l_data,alphas,n_comp);
  for(i = 0; i < n_comp; i++){
    val+= lighting_compact_phi(i,x,l_data)*alphas[i];
  }
  if(val < 0 ) val = 0;
  return val ;
}

void lighting_compact_write_parameters(FILE* arq,void* l_data){
  lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
  fprintf(arq,"%9.6lf %9.6lf %9.6lf %9.6lf %9.6lf %9.6lf %9.6lf        ",cl->Ec,cl->dirlight.c[0],cl->dirlight.c[1],cl->dirlight.c[2],cl->err_dirlight,cl->radlight, cl->err_radlight );
  fprintf(arq,"%9.6lf %9.6lf %9.6lf         ",cl->El,cl->shine,cl->err_shine);
  fprintf(arq,"%9.6lf         ",cl->Ei);
  fprintf(arq,"%9.6lf %9.6lf %9.6lf %9.6lf %9.6lf       ",cl->Ef,cl->backnormal.c[0],cl->backnormal.c[1],cl->backnormal.c[2],cl->err_backnormal);
  fprintf(arq,"%d %d %d %d        ", cl->Hc, cl->Hl, cl->Hi, cl->Hf);
  fprintf(arq,"%9.6lf %9.6lf %9.6lf",cl->view_dir.c[0],cl->view_dir.c[1],cl->view_dir.c[2]);
  fprintf(arq,"\n");
}

void* lighting_compact_read_parameters(FILE* arq){
  lighting_compact_data_t* cl = (lighting_compact_data_t*)malloc(sizeof(lighting_compact_data_t));
  fscanf(arq,"%lf %lf %lf %lf %lf %lf %lf",&(cl->Ec),&(cl->dirlight.c[0]),&(cl->dirlight.c[1]),&(cl->dirlight.c[2]),&(cl->err_dirlight),&(cl->radlight),&(cl->err_radlight) );
  fscanf(arq,"%lf %lf %lf",&(cl->El),&(cl->shine),&(cl->err_shine));
  fscanf(arq,"%lf",&(cl->Ei));
  fscanf(arq,"%lf %lf %lf %lf %lf",&(cl->Ef),&(cl->backnormal.c[0]),&(cl->backnormal.c[1]),&(cl->backnormal.c[2]),&(cl->err_backnormal));
  int Hc,Hl,Hi,Hf;
  fscanf(arq,"%d %d %d %d",&Hc,&Hl,&Hi,&Hf);
  fscanf(arq,"%lf %lf %lf",&(cl->view_dir.c[0]),&(cl->view_dir.c[1]),&(cl->view_dir.c[2]));
  cl->Hc = (Hc == 1 ? TRUE : FALSE);
  cl->Hl = (Hl == 1 ? TRUE : FALSE);
  cl->Hi = (Hi == 1 ? TRUE : FALSE);
  cl->Hf = (Hf == 1 ? TRUE : FALSE);
  
  return cl;
}

// double lighting_compact_compute_weight(lighting_compact_data_t* cl ,r3_t normal){
//   int i;
//   r3_t u = cl->dirlight;
//   double K = cl->shine;
//   double rho = cl->radlight;
//   double eps = cl->err_dirlight;
//   
//   /*Verifies if it is a region of self-shadow or terminator*/
//   double wt = (r3_dot(&u,&normal) < cos(rho + eps) ? 0 : 1);
//   wt = 1;
//   /*Check if it is in region of phong highlight*/
//   r3_t w_bissetriz = compute_bissetriz(u,normal);
//   double tol = 10e-5;
//   double lambda = acos(pow(tol,1/K));
//   double wl = ( r3_dot(&normal,&w_bissetriz) > cos( (eps/2.0) + lambda) ? 0 : 1);
//   fprintf(stderr,"ErrR %9.6lf %9.6lf  ErrK %9.6lf - %9.6lf \n",rho,wt,eps,wl);
//   return  wl*wt;
// }

// void lighting_compact_update_weights(void* l_data,double** X, double* F,double* weights, int n,double* sigma){
//   lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
//   int i;
//   r3_t u = cl->dirlight;
//   double K = cl->shine;
//   double rho = cl->radlight;
//   double eps = cl->err_dirlight;
//   for(i = 0; i < n; i++){
//     r3_t normal = (r3_t){{X[i][0],X[i][1],X[i][2]}};
//     /*Verifies if it is a region of self-shadow or terminator*/
// //     double wt = (r3_dot(&u,&normal) < cos(rho + eps) ? 0 : 1);
// //     /*Check if it is in region of phong highlight*/
// //     r3_t w_bissetriz = compute_bissetriz(u,normal);
// //     double tol = 10e-5;
// //     double lambda = acos(pow(tol,1/K));
// //     double wl = ( r3_dot(&normal,&w_bissetriz) > cos( (eps/2.0) + lambda) ? 0 : 1);
// //     weights[i] = wl*wt;
//     weights[i] = lighting_compact_compute_weight(cl,normal) ;
//   }
// }


double lighting_compact_compare_lightings(void* l_data1,void* l_data2){
  lighting_compact_data_t* cl1 = (lighting_compact_data_t*)l_data1;
  lighting_compact_data_t* cl2 = (lighting_compact_data_t*)l_data2;
  
  if(cl1->Hc != cl2->Hc){ return INFINITY;}
  if(cl1->Hl != cl2->Hl){ return INFINITY;}
  if(cl1->Hi != cl2->Hi){ return INFINITY;}
  if(cl1->Hf != cl2->Hf){ return INFINITY;}
    
  int n_comp1 = lighting_compact_get_number_components(l_data1);
  double v1[n_comp1];
  lighting_compact_get_alphas(l_data1,v1,n_comp1);
  
  int n_comp2 = lighting_compact_get_number_components(l_data2);
  double v2[n_comp2];
  lighting_compact_get_alphas(l_data2,v2,n_comp2);
  
  assert(n_comp1 == n_comp2);
  
  return rn_dist(n_comp1,v1,v2);
  
}


int lighting_compact_num_nl_parameters(void* l_data){
  lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
  /*Only valid params*/
  int num_parms = 0;
  if(cl->Hc && (cl->err_dirlight   > 0)) num_parms+= 2;
  if(cl->Hc && (cl->err_radlight   > 0)) num_parms+= 1;
  if(cl->Hl && (cl->err_shine      > 0)) num_parms+= 1;
  if(cl->Hf && (cl->err_backnormal > 0)) num_parms+= 2;
  return num_parms;  
}


// void convert_r3_to_stereographic(r3_t v, double* X, double* Y,r3_t view_dir);
// void convert_r3_to_stereographic(r3_t v, double* X, double* Y,r3_t view_dir){
//      r3_t v_novo;
//      r3x3_t VM = compute_normal_correction_matrix(view_dir);
//      r3x3_map_row(&v, &VM, &v_novo);
//      
//      *X = v_novo.c[0]/(1 + v_novo.c[2]);
//      *Y = v_novo.c[1]/(1 + v_novo.c[2]);
//   
// }



// r3_t convert_stereographic_to_r3(double X, double Y, r3_t view_dir);
// r3_t convert_stereographic_to_r3(double X, double Y, r3_t view_dir){
//     r3_t v;
//    
//     r3_t v_novo;
//     double den = (1+ X*X + Y*Y);
//     v_novo.c[0] = 2*X/den;
//     v_novo.c[1] = 2*Y/den;
//     v_novo.c[2] = (1 - X*X - Y*Y)/den;
//     
//     r3x3_t VM = compute_normal_correction_matrix(view_dir);
//     r3x3_inv(&VM,&VM);
//     r3x3_map_row(&v_novo, &VM, &v);
//    
//     return v;
//   }


void lighting_compact_pack_nl_parameters(void* l_data,double* parameters,int n){
  
  int num_parms = lighting_compact_num_nl_parameters(l_data);
  lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
  r3_t view_dir = cl->view_dir;
  
  
  
  int ind_parm = 0;
  /*primeiro o vetor de direção de luz*/
  if(cl->Hc){
    if(cl->err_dirlight > 0){
      r3_t u = cl->dirlight;
      double X_u,Y_u;
      convert_r3_to_stereographic(u,&X_u,&Y_u,view_dir);
      parameters[ind_parm] = X_u;
      parameters[ind_parm+1] = Y_u;
      ind_parm+=2;
    }
    if(cl->err_radlight > 0){
      parameters[ind_parm] = cl->radlight;
      ind_parm+=1;
    }
    
  }
  /*Now K and rho*/
  if(cl->Hl && (cl->err_shine > 0)){
    parameters[ind_parm] = cl->shine;
    ind_parm+=1;
  }
  
  /*Now Xi*/
  if( cl->Hf && (cl->err_backnormal > 0)){
    r3_t xi = cl->backnormal;
    double X_xi,Y_xi;
    convert_r3_to_stereographic(xi,&X_xi,&Y_xi,view_dir);
    parameters[ind_parm] = X_xi;
    parameters[ind_parm+1] = Y_xi;
    ind_parm+=2;
  }
  int i;
  for(i = 0; i < num_parms ; i++){
    assert(!isnan(parameters[i]));
  }

}

// r3_t find_ortho_vector(r3_t view_dir,r3_t u);
// /*Find a vector orthogonal {t} to view_dir and coplanar with u*/
// r3_t find_ortho_vector(r3_t view_dir,r3_t u){
//   r3_t t,w;
//       r3_cross(&view_dir,&u,&w);
//       double nw = r3_dir(&w,&w);
//       if(nw == 0){
// 	t = (r3_t){{1,0,0}};
//       }else{
// 	r3_cross(&view_dir,&w,&t);
// 	r3_dir(&t,&t);
//       }
//       return t;
// }
// 
// r3_t bend_towards(r3_t v, r3_t u, double ang);
// /*Bends a vector {u} towards unit vector {v} by {ang} radians */
// r3_t bend_towards(r3_t v, r3_t u, double ang){
//   r3_t t = find_ortho_vector(v,u);
//   double uv = r3_dot(&v,&u);
//   double ut = r3_dot(&t,&u);
//   double u1v = +uv*cos(ang) + ut*sin(ang);
//   double u1t = -uv*sin(ang) + ut*cos(ang);
//   r3_t r;
//   r3_mix(u1v,&v,u1t,&t,&r);
//   return r;
// }
/*
void r3_and_error_to_stereographic_box(r3_t u, double err_u, r3_t view_dir, double* X_min, double* X_max, double* Y_min, double* Y_max);
void r3_and_error_to_stereographic_box(r3_t u, double err_u, r3_t view_dir, double* X_min, double* X_max, double* Y_min, double* Y_max){
 if( err_u == 0){
   double X,Y;
   convert_r3_to_stereographic(u,&X,&Y,view_dir);
   *X_min = X;
   *X_max = X;
   *Y_min = Y;
   *Y_max = Y;
   return;
 }
 
 r3_t para,perp;
 r3_decomp (&u, &view_dir, &para, &perp);
 double X1,Y1,X2,Y2,rad;
 double Xc,Yc;
 
 if( r3_norm(&perp) > 1e-6 ){
   //Bend towards view_dir, i guess that i can consider viewdir = 0,0,1 
   r3_t u1  = bend_towards(view_dir,u,-err_u);
   convert_r3_to_stereographic(u1,&X1,&Y1,view_dir);
   r3_t u2  = bend_towards(view_dir,u,+err_u);
   convert_r3_to_stereographic(u2,&X2,&Y2,view_dir);
   Xc = (X1+X2)/2.0;
   Yc = (Y1+Y2)/2.0;
   rad = sqrt( (Xc-X1)*(Xc-X1) + (Yc-Y1)*(Yc-Y1));
 }else{
   X1 = Y1 = X2 = Y2 = Xc = Yc = 0;
   rad = 2*tan(err_u/2.0);
 }
 
 
 
 *X_min = Xc - rad;
 *X_max = Xc + rad;
 *Y_min = Yc - rad;
 *Y_max = Yc + rad;
 
 return;
 
}*/

void lighting_compact_get_nl_centers_and_deltas(void* l_data,double* center,double* delta,int n){
  lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
  int num_parms = lighting_compact_num_nl_parameters(l_data);
  int ind_parms = 0;
  
  if(cl->Hc){
    if(cl->err_dirlight > 0){
    double X_u_min,X_u_max,Y_u_min,Y_u_max;
    r3_and_error_to_stereographic_box(cl->dirlight,cl->err_dirlight,cl->view_dir,&X_u_min,&X_u_max,&Y_u_min,&Y_u_max);
    center[ind_parms] = (X_u_min + X_u_max)/2.0;
    delta[ind_parms] = (X_u_max - X_u_min)/2.0;
    center[ind_parms+1] = (Y_u_min + Y_u_max)/2.0;
    delta[ind_parms+1] = (Y_u_max - Y_u_min)/2.0;
    ind_parms+=2;
   }
   if(cl->err_radlight > 0){
    center[ind_parms] = cl->radlight;
    delta[ind_parms] = cl->err_radlight;
    ind_parms+=1;
   }
  }
  
  if(cl->Hl && (cl->err_shine > 0)){
    center[ind_parms] = cl->shine;
    delta[ind_parms] = cl->err_shine;
    ind_parms+=1;
  }
  
  if(cl->Hf && (cl->err_backnormal > 0)){
    double X_xi_min,X_xi_max,Y_xi_min,Y_xi_max;
    r3_and_error_to_stereographic_box(cl->backnormal,cl->err_backnormal,cl->view_dir,&X_xi_min,&X_xi_max,&Y_xi_min,&Y_xi_max);
    center[ind_parms] = (X_xi_min + X_xi_max)/2.0;
    delta[ind_parms] = (X_xi_max - X_xi_min)/2.0;
    center[ind_parms+1] = (Y_xi_min + Y_xi_max)/2.0;
    delta[ind_parms+1] = (Y_xi_max - Y_xi_min)/2.0;
    ind_parms+=2;
  }
  
  assert(ind_parms == num_parms);
}

void lighting_compact_unpack_nl_parameters(double* parameters,int n,void* l_data){
  lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
  int num_parms = lighting_compact_num_nl_parameters(l_data);
  assert(n == num_parms);
  int i;
  r3_t view_dir = cl->view_dir;
  for(i = 0; i < num_parms; i++){
    assert(!isnan(parameters[i]));
  }
  /*Copy unpacked data*/

  int ind_parm = 0;
  if(cl->Hc){  
    if(cl->err_dirlight > 0 ){
      double X_u = parameters[ind_parm];
      double Y_u = parameters[ind_parm+1];
      cl->dirlight = convert_stereographic_to_r3(X_u,Y_u,view_dir);
      ind_parm+=2;
    }
    if(cl->err_radlight > 0){
      cl->radlight = parameters[ind_parm];
      ind_parm+=1;
    }
  }
  
  if( cl->Hl && (cl->err_shine > 0)){
    cl->shine = parameters[ind_parm];
    ind_parm+=1;
  }
  
  if( cl->Hf && (cl->err_backnormal > 0)){
    double X_xi = parameters[ind_parm];
    double Y_xi = parameters[ind_parm+1];
    cl->backnormal = convert_stereographic_to_r3(X_xi,Y_xi,view_dir);
    ind_parm+=2;
  }
    
  return;
}




double lighting_compact_compute_new_scalar_error(double t_old, double eps_old, double t, double alpha, double beta, double gamma);
double lighting_compact_compute_new_scalar_error(double t_old, double eps_old, double t, double alpha, double beta, double gamma){
  double delta = fabs(t - t_old); /*parameter change*/
  double eps = fmax(fmin(alpha*delta,beta*eps_old), gamma*eps_old); /*New uncertainty estimate*/ 
  return eps;
}

double lighting_compact_compute_new_direction_error(r3_t u_old, double eps_old, r3_t u, double alpha, double beta, double gamma);
double lighting_compact_compute_new_direction_error(r3_t u_old, double eps_old, r3_t u, double alpha, double beta, double gamma){
  double delta = acos(r3_dot(&u_old,&u));
  return fmin(fmax(fmin(alpha*delta,beta*eps_old), gamma*eps_old),M_PI/2.0);
}




void lighting_compact_update_params_and_errors(double* parameters,int n,void* l_data, double alpha, double beta, double gamma){
  lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
  lighting_compact_data_t cl_old = *cl;
  lighting_compact_unpack_nl_parameters(parameters,n,l_data);
  
  if(cl->Hc){
    if(cl->err_dirlight > 0){
      clip_direction_to_circle(&(cl->dirlight), cl_old.dirlight, cl_old.err_dirlight);
      cl->err_dirlight = lighting_compact_compute_new_direction_error(cl_old.dirlight, cl_old.err_dirlight,cl->dirlight,alpha,beta,gamma);
      if(cl->err_dirlight < 0.005){ fprintf(stderr,"!! Err_dirlight too small\n"); cl->err_dirlight = 0;}
    }
    if(cl->err_radlight > 0){
      clip_scalar_to_interval(&(cl->radlight),cl_old.radlight,cl_old.err_radlight);
      cl->err_radlight = fmin(cl->radlight, lighting_compact_compute_new_scalar_error(cl_old.radlight,cl_old.err_radlight,cl->radlight,alpha,beta,gamma));
      if(cl->err_radlight <= 0.005){ fprintf(stderr,"!! Err_radlight too small\n"); cl->err_radlight = 0;}
    }
  }
  if(cl->Hl && (cl->err_shine > 0)){
    clip_scalar_to_interval(&(cl->shine),cl_old.shine,cl_old.err_shine);
    cl->err_shine = fmin(cl->shine, lighting_compact_compute_new_scalar_error(cl_old.shine,cl_old.err_shine,cl->shine,alpha,beta,gamma));
    if(cl->err_shine <= 0.05*cl->shine){ fprintf(stderr,"!! Err_shine too small\n"); cl->err_shine = 0;}
  }
  
  if(cl->Hf && (cl->err_backnormal > 0)){
    clip_direction_to_circle(&(cl->backnormal), cl_old.backnormal, cl_old.err_backnormal);
    cl->err_backnormal = lighting_compact_compute_new_direction_error(cl_old.backnormal, cl_old.err_backnormal,cl->backnormal,alpha,beta,gamma);
    if(cl->err_backnormal <= 0.005){ fprintf(stderr,"!! Err_backnormal too small\n"); cl->err_backnormal = 0;}
  }
  
}




void lighting_compact_generate_simplex(void* l_data, int n, double* S){
  int num_parms = lighting_compact_num_nl_parameters(l_data);
  
  assert(n == num_parms);
  double delta[num_parms];
  double center[num_parms];
  
  lighting_compact_get_nl_centers_and_deltas(l_data,center,delta,n);
  
  rmxn_regular_simplex(n, S);
  rmxn_spin_rows(n+1, n,S,S);
  double simplex_radius = rmxn_regular_simplex_radius(n);
  /*now modify S so that it lies inside the uncertainty intervals*/
  int i;
  for(i = 0; i <=n ; i++){
    double* Si = &(S[i*n]);
    int j;
    for(j = 0; j < n; j++){
      Si[j] = (Si[j]/simplex_radius)*delta[j] + center[j];
      assert(fabs(Si[j] - center[j]) <= 1.00001*delta[j]);
    }
  }
  
}

/*Functions for light source direction estimation !*/

 double lighting_compact_simple_phi(int r, double* x,void* l_data){
   lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
   r3_t normal = (r3_t){{x[0],x[1],x[2]}};
   if((r >=0) && (r < 3)){
      assert(cl->Hc);
      double s = r3_dot(&(cl->dirlight), &normal);
      return (s <= 0 ? 0 : normal.c[r]);
   }else if( r == 3){
      return (cl-> Hf ? 0.5*(1 - normal.c[2]) : 1);
   }
   assert(FALSE);
 }
 
 int lighting_compact_simple_get_number_components(void* l_data){
   lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
   if(cl->Hi || cl->Hf) return 4;
   return 3;
 }
 
 void lighting_compact_simple_set_alphas(double* C,int n, void* l_data){
   int num_comp = lighting_compact_simple_get_number_components(l_data);
   assert(n == num_comp);
   lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
   r3_t u = (r3_t){{ C[0],C[1],C[2]}};
   double  Ec = r3_dir(&u,&u);
   cl->dirlight = u;
   cl->Ec = Ec;
   if(cl->Hi || cl->Hf){
     double Ef,Ei;
     Ei = Ef = 0;
     if(cl->Hf){
       Ef = C[3];
     }else{
       Ei = C[3];
     }
     cl->Ef = Ef;
     cl->Ei = Ei;
   }
 }

void lighting_compact_simple_get_alphas(void* l_data,double* C,int n){
  int num_comp = lighting_compact_simple_get_number_components(l_data);
   assert(n == num_comp);
   lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
   C[0] = cl->Ec*cl->dirlight.c[0];
   C[1] = cl->Ec*cl->dirlight.c[1];
   C[2] = cl->Ec*cl->dirlight.c[2];
   if(cl->Hi || cl->Hf){
     C[3] = (cl->Hf ? cl->Ef : cl->Ei);
   }
}


approx_model_t* lighting_compact_simple_create_approx_lighting_model(void){
  approx_model_t* am = (approx_model_t*)malloc(sizeof(approx_model_t));
  am->type = COMPACT_MODEL;

  am->phi = lighting_compact_simple_phi;
  am->get_num_components = lighting_compact_simple_get_number_components;
  am->set_alphas = lighting_compact_simple_set_alphas;
  am->get_alphas = lighting_compact_simple_get_alphas;
  am->copy_data = lighting_compact_copy_lighting_data;
  am->release_data =  lighting_compact_release_lighting_data;
  
  /*Optional members - they wont be called by the fitting functions and can be NULL*/
  am->evaluate = lighting_compact_shading; 
  am->write_parameters = lighting_compact_write_parameters;
  am->read_parameters = lighting_compact_read_parameters; 
  /*for non linear models - NOT USED*/
  am->compare = NULL; 
  am->get_num_nl_parameters = NULL;
  am->pack_nl_parameters =  NULL; 
  am->get_nl_centers_and_deltas = NULL;
  am->unpack_nl_parameters = NULL;
  am->update_nl_parameters_and_errors = NULL;
  am->generate_nl_parameter_simplex = NULL; 
  return am;
}


double  light_compact_simple_compute_weight(lighting_compact_data_t* cl,r3_t normal){
//   assert(cl->Hc);
  double eps = cl->err_dirlight;
  double rhomax = cl->radlight + cl->err_radlight;
  
  r3_t w_r; /*Bisectrix of {dirlight} with {view_dir}*/
  r3_add(&(cl->dirlight),&(cl->view_dir),&w_r);
  r3_dir(&w_r,&w_r);
  double tol  = 0.1;
  double Kmin = cl->shine - cl->err_shine;
  double lambda = acos(pow(tol,1.0/Kmin));
  double s = r3_dot(&(cl->dirlight), &normal);
  /*Eliminate terimnator points*/
  double weight_terminator = 1.0;
  if(cl->Hc){
    weight_terminator = (s < sin(rhomax + eps) ? 0 : 1.0);
  }
  double weight_highlight  = 1.0;
  if(cl->Hl){
    weight_highlight = ( r3_dot(&normal,&w_r) > cos((eps/2.0) + lambda) ? 0 : 1.0);
  }
  return  weight_terminator*weight_highlight;
}

void light_compact_simple_update_weights(void* l_data,double** X,double* weights, int n){
  lighting_compact_data_t* cl = (lighting_compact_data_t*)l_data;
  int i;
  assert(cl->Hc);
  double eps = cl->err_dirlight;
  double rhomax = cl->radlight + cl->err_radlight;
  
  r3_t w_r; /*Bisectrix of {dirlight} with {view_dir}*/
  r3_add(&(cl->dirlight),&(cl->view_dir),&w_r);
  r3_dir(&w_r,&w_r);
  double tol  = 0.2;
  double Kmin = cl->shine - cl->err_shine;
  double lambda = acos(pow(tol,1.0/Kmin));
  
  
  
  fprintf(stderr,"RHOMAX %9.6lf (%9.6lf) EPS %9.6lf (%9.6lf) KMIN %9.6lf LAMBDA %9.6lf (%9.6lf)\n",rhomax, rhomax*180/M_PI,eps, eps*180/M_PI,Kmin,lambda,lambda*180/M_PI);
  for(i = 0; i < n; i++){
    double*x = X[i];
    r3_t normal = (r3_t){{x[0],x[1],x[2]}};
//     double s = r3_dot(&(cl->dirlight), &normal);
//     /*Eliminate terimnator points*/
//     double weight_terminator = (s < sin(rhomax + eps) ? 0 : 1.0);
//     
//     double weight_highlight  = 1.0;
//     if(cl->Hl){
//       weight_highlight = ( r3_dot(&normal,&w_r) > cos((eps/2.0) + lambda) ? 0 : 1.0);
//     }
//     
//     weights[i] = weight_terminator*weight_highlight;
      weights[i] = light_compact_simple_compute_weight(cl,normal);
  }
}
/*
void dump_weights(FILE* arq,double** X,double* w, int n){
  int i;
  for(i = 0; i < n; i++){
    fprintf(arq,"%9.6lf %9.6lf %9.6lf %9.6lf\n",X[i][0],X[i][1],X[i][2],w[i]);
  }
}*/

void light_compact_estimate_lightdir(double** X, double* F,double* wpos, int n,void* l_data,int update_steps){
  approx_model_t* am = lighting_compact_simple_create_approx_lighting_model();
  double* weights = (double*)malloc(sizeof(double)*n);
  rn_all(n,1,weights);
  int num_steps = 0;
  do{
//    char* filename= NULL;
//     char *filename = jsprintf("dump_%02d.txt",num_steps);
//     FILE* dump_f = open_write(filename,TRUE);
    approx_model_fit_linear_model(X,F,weights,wpos,n, am,l_data);
    light_compact_simple_update_weights(l_data,X,weights,n);
//     dump_weights(dump_f,X,weights,n);
//     fclose(dump_f);
    num_steps++;
  }while(num_steps < update_steps);
  
  free(weights);
  free(am);
}
/*

void light_compact_estimate_lightdir(double** X, double* F, int n,
		    ls_model_t* lm,
		     void* l_data,
		      int update_steps
		     ){
	double* A;
	double* b;
	double* c;
	double* w; //dynamically adjusted weight, probability of goodness of a point of table
	double* wpos; //a priori weights
	
	int basis_size = lm->get_num_components(l_data);
	
	A = (double*)malloc(sizeof(double)*(basis_size*basis_size));
	b = (double*)malloc(sizeof(double)*basis_size);
	c = (double*)malloc(sizeof(double)*basis_size);
	//w = (double*)malloc(sizeof(double)*n);
	int i;

	if(lm->wpos != NULL){
	   wpos = lm->wpos;
	}else{
	  wpos = (double*)malloc(sizeof(double)*n);
	}
	
	w = lm->weights;
	assert(w != NULL);
	 
	 

        // Compute matrix that maps true normal to view-relative normal:
	double sigma = 0.2;
	//fprintf(stderr,"Basis size is %d\n",lm->get_num_components(l_data));
	//fprintf(stderr,"Starting weight-adjusting iteration SIGMA = %9.6f\n",sigma);
	int passo = 0;
	//int npassos = 3;
	
//	int luz;
	double threshold = 0.01;
// 	if(lm->type == 0){
// 	    fprintf(stderr," Start : ");
// 	    lm->write_param(stderr,l_data);
// 	}
	for(passo = 0; passo < update_steps; passo++){
//		double sigma2 = sigma*sigma;
		//fprintf(stderr,"Weight adjustment step %d\n",passo);
		
		int num_iters = 0;
		int MAX_ITER = 100;
		double epsilon = 10e-5;
		double diff = 1.0; 
		while( (num_iters < MAX_ITER) && (diff > epsilon)){
		//  fprintf(stderr,"Generating system matrix...");
		
	
		  computeLeastSquaresMatrix(A,lm->phi,basis_size,X,n,w,wpos,l_data);
		//  fprintf(stderr,"OK\n");
		  void* l_data_old = lm->copy_data(l_data);
		  computeLSTerms(A,
			      b,
			      c,
			      lm->phi,
			      basis_size,
			      lm->validate_results,
			      X,F,n,
			      w,
			      wpos,
			      l_data
		  );
	
		  lm->retrieve_components(c,basis_size,l_data);
		  diff = lm->compare(l_data,l_data_old);
// 
		//  fprintf(stderr,"Finished Iteration %d - (diff %lf vs %lf) \n",num_iters,diff,epsilon);
		if(lm->type == 3){
		  lm->write_param(stderr,l_data);
		}
		  lm->release_data(l_data_old);
		  num_iters++;
		}
		
		//fprintf(stderr,"Computing errors and adjusting weight...");
		lm->update_weights(l_data,X,F,w,n,&sigma);
		
		int bad_lines = 0;
		double sumw = 1.0e-200;
		for(i = 0; i<  n; i++){
		  sumw+=w[i];
		}
		for(i = 0; i<  n; i++){
			w[i] = (w[i]/sumw)*n;
			if(w[i] < threshold) bad_lines++;
		}
// 			if(lm->type == 0){
// 			    fprintf(stderr," Step : ");
// 			     lm->write_param(stderr,l_data);
// 			 }
		//fprintf(stderr,"OK - %d  bad lines in table (%9.6f %%)\n",bad_lines,100.0*bad_lines/(double)n);
		
	}

	lm->num_weights =  n;
	lm->weights = w;
	lm->wpos = wpos;
	
        free(A);
	free(b);

	free(c);
	
} */
