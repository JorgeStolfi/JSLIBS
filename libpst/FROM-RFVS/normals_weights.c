
#define _GNU_SOURCE
#include <math.h>
#include <normals_weights.h>



double computeWalb(double albedo, double Amax){
  double t = albedo/Amax;
  double s = fmax(0,4*t*(1-t));
  double walb = 1 - ((1-s)*(1-s));
  return walb;
}

double computeWnrm(r3_t norm, r3_t avg_dir, double r){
  double s = r3_dot(&avg_dir,&norm);
  double wnrm = (s < cos((90-2*r)*M_PI/180.0) ? 0 : 1-s);
  return wnrm;
}

double computeWsmp(double S[], int nLights){
  int i;
  double prod = 1;
  for(i = 0; i < nLights; i++){
	  prod*=(1-S[i]);
  }
  double Smsq = 1 - pow(prod,1/(double)nLights);
  double A = 4*(1 - Smsq)*Smsq;
  double wsmp = 1 - ((1-A)*(1-A));
  return wsmp;
}

double computeWhigh(r3_t norm,r3_t cluster_dir,double k,double r){
  r3_t Nhigh = cluster_dir;
  Nhigh.c[2] += 1;
  r3_dir(&Nhigh,&Nhigh);
  double T = r3_dot(&Nhigh,&norm);
  T =  cos(acos(T) + (r*M_PI/180.0));
  double whigh = (T <= 0 ? 1 : 1 - pow(T,k));
  return whigh;
}
