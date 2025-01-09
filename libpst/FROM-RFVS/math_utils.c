#include <math_utils.h>
#include <r3x3.h>
#include <tabela.h>



r3x3_t compute_normal_correction_matrix(r3_t view_dir){
	r3x3_t roda_normal;
	double sR = hypot(view_dir.c[0],view_dir.c[1]);
	double cR = view_dir.c[2];
	if(sR != 0) {
		
		r3x3_t A,R;
		double ux = view_dir.c[0]/sR;
		double uy = view_dir.c[1]/sR;
	
	
		A = (r3x3_t){
			{
			{ux,0.0,-uy},
			{uy, 0, ux},
			{0.0,1.0,0.0}
			}
		};
// 		fprintf(stderr,"A\n");
// 		r3x3_print(stderr,&A);
// 		fprintf(stderr,"Determinant = %lf \n",r3x3_det(&A));
		r3x3_t temp;
		r3x3_mul_tr(&A,&A,&temp);
/*		fprintf(stderr,"A*Atr\n");
		r3x3_print(stderr,&temp);*/
		R = (r3x3_t){
			{
				{cR,-sR, 0.0},
				{sR,cR,0.0},
				{0.0,0.0,1.0}
			}
		};
		r3x3_mul(&A,&R,&roda_normal);
		r3x3_mul_tr(&roda_normal,&A,&roda_normal);
// 		fprintf(stderr,"Roda Normal\n");
// 		r3x3_print(stderr,&roda_normal);
// 		fprintf(stderr,"Determinant = %lf \n",r3x3_det(&roda_normal));
// 		r3x3_mul_tr(&roda_normal,&roda_normal,&temp);
// 		fprintf(stderr,"M*Mtr\n");
// 		r3x3_print(stderr,&temp);
	}else {
		r3x3_ident(&roda_normal);
	}
	return roda_normal;
}

r3_t compute_bissetriz(r3_t u, r3_t v);
r3_t compute_bissetriz(r3_t u, r3_t v){
  r3_t w;
  r3_add(&u,&v,&w);
  (void) r3_dir(&w,&w);
  return w;
}

void convert_r3_to_stereographic(r3_t v, double* X, double* Y,r3_t view_dir){
     r3_t v_novo;
     r3x3_t VM = compute_normal_correction_matrix(view_dir);
     r3x3_map_row(&v, &VM, &v_novo);
     
     *X = v_novo.c[0]/(1 + v_novo.c[2]);
     *Y = v_novo.c[1]/(1 + v_novo.c[2]);
  
}


r3_t find_ortho_vector(r3_t view_dir,r3_t u){
  r3_t t,w;
      r3_cross(&view_dir,&u,&w);
      double nw = r3_dir(&w,&w);
      if(nw == 0){
	t = (r3_t){{1,0,0}};
      }else{
	r3_cross(&view_dir,&w,&t);
	r3_dir(&t,&t);
      }
      return t;
}


r3_t bend_towards(r3_t v, r3_t u, double ang){
  r3_t t = find_ortho_vector(v,u);
  double uv = r3_dot(&v,&u);
  double ut = r3_dot(&t,&u);
  double u1v = +uv*cos(ang) + ut*sin(ang);
  double u1t = -uv*sin(ang) + ut*cos(ang);
  r3_t r;
  r3_mix(u1v,&v,u1t,&t,&r);
  return r;
}

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
   /*Bend towards view_dir, i guess that i can consider viewdir = 0,0,1*/
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
 
}

r3_t convert_stereographic_to_r3(double X, double Y, r3_t view_dir){
    r3_t v;
   
    r3_t v_novo;
    double den = (1+ X*X + Y*Y);
    v_novo.c[0] = 2*X/den;
    v_novo.c[1] = 2*Y/den;
    v_novo.c[2] = (1 - X*X - Y*Y)/den;
    
    r3x3_t VM = compute_normal_correction_matrix(view_dir);
    r3x3_inv(&VM,&VM);
    r3x3_map_row(&v_novo, &VM, &v);
   
    return v;
  }



void clip_direction_to_circle(r3_t* u, r3_t u_old, double err_old){
   r3_t para,perp;
   r3_dir(&u_old,&u_old);
   r3_decomp (u, &u_old, &para, &perp);
   double ct = r3_norm(&para)/r3_norm(u); /*sine of angle between u and u_old*/
   double ce = cos(err_old);
   if( ct < ce ){
     if(err_old < 1.0e-8){ *u = u_old; return; };
     double st = r3_norm(&perp)/r3_norm(u);
     if( st < 1.0e-8){ *u = u_old; return; }
     double se = sin(err_old);
     r3_scale(se/st,&perp,&perp);
     r3_scale(ce/ct,&para,&para);
     r3_add(&perp,&para,u);
     r3_dir(u,u);
   }
  
}


void clip_scalar_to_interval(double* v, double v_old, double err_old){
  *v = fmin(v_old + err_old, fmax(v_old - err_old,*v));
}

double gaussian_distribution(double x, double mean, double sigma){
  auto double normal_dist(double x);
  double normal_dist(double x){
    return exp(-x*x/2.0)/sqrt(2.0*M_PI);
  }
  return normal_dist((x - mean)/sigma)/sigma;
}