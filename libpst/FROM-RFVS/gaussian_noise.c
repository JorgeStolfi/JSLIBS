/* Last edited on 2017-01-04 19:25:46 by stolfilocal */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <gaussian_noise.h>

double uniform(int* seed){
  if(*seed == 1){
    *seed = (int)time(NULL);
    srand(*seed);
  }
  return (rand()/(double)RAND_MAX);
}

double BoxMuller(double mean, double sigma,int* seed){
  /*
    The Box-Muller method uses the technique of inverse transformation to turn two uniformly distributed randoms,
    U1 and U2, into two unit normal randoms, X and Y.
  */
  double X,Y,U1,U2,V1,V2,S;
  do{
    U1=uniform(seed);            /* U1=[0,1] */
    U2=uniform(seed);           /* U2=[0,1] */
    V1=2.0 * U1 -1.0;            /* V1=[-1,1] */
    V2=2.0 * U2 -1.0;           /* V2=[double uniform(void){
  if(seed == 1){
    seed = time(NULL);
    srand(seed);
  }
  return (rand()/(double)RAND_MAX);
}
-1,1] */
    S=V1 * V1 + V2 * V2;
  }while( S >=1.0);
	
  X=sqrt(-2.0 * log(S) / S) * V1;
  Y=sqrt(-2.0 * log(S) / S) * V2;
  /*
    The above is called the Polar Method and is fully described in the Ross book [Ros88].
    X and Y will be unit normal random variables (mean = 0 and (sigma*sigma) = 1), and can be easilly modified for
    different mean and (sigma*sigma).
  */
  X = mean + sigma * X;
  Y = mean + sigma * Y;
  return X;
 }


double CLTM(int N, double mean, double sigma,int* seed){
  /*
    The Central Limit Theorm Method (CLTM) states that the sum
    of N randoms will approach normal distribution as N
    approaches infinity. We can outline an algorithm that uses this approach as:
  */
  double X,U;
  int i;
  X=0;
  for (i = 1; i <=N;i++){
    //U = uniform()
		
    U = uniform(seed);
    X = X + U;
  }
	
  /* for uniform randoms in [0,1], mu = 0.5 and var = 1/12 */
  /* adjust X so mu = 0 and var = 1 */
  X = X - ((double)N/2);            /* set mean to 0 */
  X = X * sqrt(12/(double)N);       /* adjust variance to 1 */	
  /*
    When the algorithm finishes, X will be our unit normal
     random. X can be further modified to have the given mean value 
    {mean} and  deviation {sigma}, e.g.:
  */
  X = mean + sigma * X;
  return X;
  /*
    The drawback to this method is that X will be in the range
    [-N/2, +N/2]*sqrt(12/(double)N), instead of (-Infinity, Infinity) and if
    the calls to uniform are not truly independent, then the 
    noise will no longer be white. Jeruchim, et. al.,
    recommend N >=20 for good results.
  */
}
