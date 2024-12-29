#include <stdio.h>
#include <jsfile.h>
#include <float_image.h>
#include <lighting_compact.h>
#include <estlogprob.h>
#include <r3.h>

 
r3_t normal_from_azi_inc(double azi_rad, double inc_rad);
r3_t normal_from_azi_inc(double azi_rad, double inc_rad){
  return (r3_t){{sin(inc_rad)*cos(azi_rad),sin(inc_rad)*sin(azi_rad),cos(inc_rad)}};
}

lighting_compact_data_t create_ldata(r3_t light_dir);
lighting_compact_data_t create_ldata(r3_t light_dir){
  
  lighting_compact_data_t  ldata;
  ldata.view_dir = (r3_t){{0,0,1}};
  ldata.dirlight = light_dir;
  ldata.radlight = 0;
  ldata.shine = 0;
  ldata.backnormal = (r3_t){{0,0,0}}; /*Background*/
  
  /*Uncertainties of the non linear coefficients*/
  ldata.err_dirlight = ldata.err_radlight = ldata.err_shine = ldata.err_backnormal = 0;
  
  ldata.Ec = 0.5;
  ldata.El = ldata.Ei = ldata.Ef = 0; /*Backplane light source intensity*/
  
  /*Which basis elements to include*/ 
  ldata.Hc = TRUE;
  ldata.Hl = ldata.Hi = ldata.Hf = FALSE;
  
  return ldata;
}

r3_t create_random_normal(void);
r3_t create_random_normal(void){
  double rand_azi = (rand()/(double)RAND_MAX)*M_PI;
  double rand_inc = (rand()/(double)RAND_MAX);
  rand_inc*=rand_inc;
  rand_inc*=(M_PI/2);
  fprintf(stderr,"Random AZ is %9.6lf and IC is %9.6lf\n",rand_azi,rand_inc);
  return normal_from_azi_inc(rand_azi, rand_inc);
}


double ProbSiGialb_CBP_shadow(double Si, double Gi, double sigma, double omg0,double omg1);
double ProbSiGialb_CBP_shadow(double Si, double Gi, double sigma, double omg0,double omg1){
  double R = K_09*sigma;
  double alo = (Si < 2*R ? 0.0 : (Si > Gi + R ? 1.0 : (Si - R)/Gi ));
  return ProbSiGialb(Si,Gi,0.5*(0 + alo),sigma,omg0,omg1);
}

double ProbSiGialb_CBP_highlight(double Si, double Gi, double sigma, double omg0,double omg1);
double ProbSiGialb_CBP_highlight(double Si, double Gi, double sigma, double omg0,double omg1){
  double R = K_09*sigma;
  double ahi = ((Si > Gi - R) || (Si > 1 - 2*R) ? 1.0 : (Si + R)/Gi);
  return  ProbSiGialb(Si,Gi,0.5*(1 + ahi),sigma,omg0,omg1);
}


double ProbSiGialb_CBP_canonical(double Si, double Gi, double sigma, double omg0,double omg1);
double ProbSiGialb_CBP_canonical(double Si, double Gi, double sigma, double omg0,double omg1){
  double R = K_09*sigma;
  double alo = (Si < 2*R ? 0.0 : (Si > Gi + R ? 1.0 : (Si - R)/Gi ));
  double ahi = ((Si > Gi - R) || (Si > 1 - 2*R) ? 1.0 : (Si + R)/Gi);
  double amd = 0.5*(alo + ahi);
  double lam = 0.57735026918962576450; /* sqrt(1/3) para quadratura de Gauss. */
  double grd = lam * (ahi - amd);
  
  double P0 = ProbSiGialb(Si,Gi,amd - grd,sigma,omg0,omg1);
  double P1 = ProbSiGialb(Si,Gi,amd + grd,sigma,omg0,omg1);
  return  0.5*(P0 + P1);
  
}


double ProbSiGialb_CBP(double Si, double Gi,double alb, double sigma, double omg0,double omg1);
double ProbSiGialb_CBP(double Si, double Gi,double alb, double sigma, double omg0,double omg1){
  double R = K_09*sigma;
  double alo = (Si < 2*R ? 0.0 : (Si > Gi + R ? 1.0 : (Si - R)/Gi ));
  double ahi = ((Si > Gi - R) || (Si > 1 - 2*R) ? 1.0 : (Si + R)/Gi);
  if(alb < alo) return ProbSiGialb_CBP_shadow(Si,Gi,sigma,omg0,omg1);
  if(alb > ahi) return ProbSiGialb_CBP_highlight(Si,Gi,sigma,omg0,omg1);
    
  return  ProbSiGialb_CBP_canonical(Si,Gi,sigma,omg0,omg1);
  
}

int main(int argc, char** argv){
  
//   srand ( time(NULL) );
  
  
  
  
  


  double omg0 = 0.1;
  double omg1 = 0.15;
  double sigma = 0.03;
  double delta=0.001;
  
  int i;
  
  double Si = 0;
  double Gi= 0.5;
  double alb = 0.7;
  
  
  
  FILE* arq = open_write("graph_SG.txt",TRUE);
  while( Si <= 1.0){
    
    double PR,PRSG,PRW0,PRW1;
    PR = PRSG = PRW0 = PRW1 = 0.0;
    PR=ProbSiGialb(Si, Gi,alb,sigma,omg0, omg1);
    PRSG=ProbSiGialb_canonical(Si, Gi,alb,sigma);
    PRW0=ProbSiGialb_shadow(Si, Gi,alb,sigma);
    PRW1=ProbSiGialb_highlight(Si, Gi,alb,sigma);
    
    
    fprintf(arq,"%22.15e %22.15e %22.15e %22.15e %22.15e \n",Si, PR,PRSG,PRW0,PRW1);
    Si+=delta;
  }
  fclose(arq);
  
  
  arq = open_write("graph_SP.txt",TRUE);
  Si = 0;
  while( Si <= 1.0){
    
    double PR,PRSG,PRW0,PRW1;
    PR = PRSG = PRW0 = PRW1 = 0.0;
    PR=ProbSiGialb_CBP(Si, Gi,alb,sigma,omg0, omg1);
    PRSG=ProbSiGialb_CBP_canonical(Si, Gi,sigma,omg0, omg1);
    PRW0=ProbSiGialb_CBP_shadow(Si, Gi,sigma,omg0, omg1);
    PRW1=ProbSiGialb_CBP_highlight(Si, Gi,sigma,omg0, omg1);
    
    
    fprintf(arq,"%22.15e %22.15e %22.15e %22.15e %22.15e \n",Si, PR,PRSG,PRW0,PRW1);
    Si+=delta;
  }
  fclose(arq);
  
  int m = 10;
  double S[m];
  double G[m];
  r3_t normalS = (r3_t){{0,0,1}};
  r3_t normalG = normalS;
  double albG = 0.9;
  double albS = 0.6;
  double albRel = albS/albG;
  
  fprintf(stderr,"Albedo S %9.6lf Albedo G %9.6lf Albedo Rel %9.6lf\n",albS,albG,albRel);
  
  lighting_compact_data_t ldata[m];
  
  for(i = 0; i < m; i++){
    r3_t ldir = create_random_normal();
    ldata[i] = create_ldata(ldir);
    S[i] = albS*lighting_compact_shading(normalS.c,&(ldata[i]));
    G[i] = albG*lighting_compact_shading(normalG.c,&(ldata[i]));;
    fprintf(stderr,"Light Dir[%02d]: (%9.6lf %9.6lf %9.6lf)\n",i,ldir.c[0],ldir.c[1],ldir.c[2] );
  }
  
  
  
  double defect_qt = 5;
  for(i = 0; i < defect_qt; i++){
      int ind = rand()%m;
      double new = rand()/(double)RAND_MAX;
      fprintf(stderr,"Point S[%d] is has anomaly: %9.6lf -> %9.6lf\n",ind,S[ind],new);
      S[ind] = new;
  }
  
  fprintf(stderr,"G - ( ");
  for(i = 0; i < m; i++){
    fprintf(stderr," %9.6lf ",G[i]);
  }
  fprintf(stderr,")\n");
  
  fprintf(stderr,"S - ( ");
  for(i = 0; i < m; i++){
    fprintf(stderr," %9.6lf ",S[i]);
  }
  fprintf(stderr,")\n");
  
  
  FILE* arq_sag = open_write("graph_SAG.txt",TRUE);
    
  double delta_alb = 0.001;
  for( alb = 0; alb <=1; alb+=delta_alb){
    double PRSAG = 1.0;
    fprintf(arq_sag,"%22.15e ",alb);
    for(i =0; i < m; i++){
      double PrSiGi = ProbSiGialb(S[i], G[i],alb,sigma,omg0, omg1);
      fprintf(arq_sag,"%22.15e ", PrSiGi);
      PRSAG*= PrSiGi;
    }
    fprintf(arq_sag,"%22.15e \n", PRSAG);
  }
  
  fclose(arq_sag);  
    
  
  return 0;
}


