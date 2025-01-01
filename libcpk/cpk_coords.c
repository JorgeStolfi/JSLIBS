/* See cpk_coords.h */
/* Last edited on 2024-12-31 10:46:20 by stolfi */ 

#include <stdint.h>
#include <math.h>

#include <cpk_coords.h>

#define deg2rad (M_PI / 180)
#define rad2deg (180.0 / M_PI)

/* Parameters for the "South American 1969" reference ellipsoid: */
#define EquatorialRadius (6378160)
#define EccentricitySquared (0.006694542)

void LL_to_EUTM(double Lat, double Lon, double *dx, double *dy, double refLon)
  {
    /* Equaçoes provenientes do USGS Bulletin 1532. */
    /* Renamed {eccSquared -> e2, eccPrimeSquared -> e2P}. */
    /* Renamed {a -> ER}. */

    double LatRad = Lat * deg2rad;
    double LonRad = (Lon - refLon) * deg2rad;
    
    double ER = EquatorialRadius;
    double e2 = EccentricitySquared, e4 = e2*e2, e6 = e4*e2;
    double e2P = e2/(1 - e2);

    double MBas = LatRad * (1 - e2/4 - 3*e4/64 - 5*e6/256);
    double MCorr = 
      - (3*e2/8 + 3*e4/32  +  45*e6/1024) * sin(2*LatRad)
      + (15*e4/256 + 45*e6/1024) * sin(4*LatRad)
      - (35*e6/3072) * sin(6*LatRad);
    double M = ER * (MBas + MCorr);

    double sinLat = sin(LatRad);
    double cosLat = cos(LatRad);
    double tanLat = tan(LatRad);
    double N = ER / sqrt(1 - e2*sinLat*sinLat);
    double T = tanLat*tanLat, T2 = T*T;
    double C = e2P*cosLat*cosLat, C2 = C*C;
    double A = cosLat * LonRad;
    double A2 = A*A, A3 = A2*A, A4 = A2*A2, A5 = A3*A2, A6 = A3*A3;

    double k0 = 0.9996;
    double ZX = 
      + A 
      + (1 - T + C) * A3/6 
      + (5 - 18*T + T2 + 72*C - 58*e2P) * A5/120;

    double ZY = 
      + A2/2 
      + (5 - T + 9*C + 4*C2)*A4/24 
      + (61 - 58*T + T2 + 600*C - 330*e2P) * A6/720;

    *dx = k0 * N * ZX;
    *dy = k0 * (M + N * tanLat*ZY);
  }

void EUTM_to_LL(double dx, double dy, double *Lat, double *Lon, double refLon)
  {
    /* Equaçoes provenientes do USGS Bulletin 1532. */
    /* Renamed {eccSquared -> e2, eccPrimeSquared -> e2P}. */
    /* Renamed {a -> ER}, {e1 -> f}, {phi1Rad -> p}. */
    /* Renamed {T1 -> T}, {R1 -> R}, {C1 -> C}, {N1 -> N}. */

    double ER = EquatorialRadius;
    double e2 = EccentricitySquared, e4 = e2*e2, e6 = e4*e2;
    double e2P = e2/(1 - e2);

    double f = (1 - sqrt(1 - e2))/(1 + sqrt(1 - e2));
    double f2 = f*f, f3 = f2*f, f4 = f2*f2;

    double k0 = 0.9996;
    double mu = dy / k0 / (ER * (1 - e2/4 - 3*e4/64 - 5*e6/256));

    double p = 
      + mu 
      + (3*f/2 - 27*f3/32) * sin(2*mu)
      + (21*f2/16 - 55*f4/32) * sin(4*mu)
      + (151*f3/96) * sin(6*mu);
    double sinp = sin(p);
    double cosp = cos(p);
    double tanp = tan(p);
    
    double N = ER / sqrt(1 - e2*sinp*sinp);
    double T = tanp*tanp, T2 = T*T;
    double C = e2P * cosp*cosp, C2 = C*C;
    double h = 1 - e2*sinp*sinp; 
    double R = ER*(1 - e2)/h/sqrt(h);
    double D = dx / k0 / N;
    double D2 = D*D, D3 = D2*D, D4 = D2*D2, D5 = D3*D2, D6 = D3*D3;

    double ZLat = 
      + D2/2 
      - (5 + 3*T + 10*C - 4*C2 - 9*e2P) * D4/24
      + (61 + 90*T + 298*C + 45*T2 - 252*e2P - 3*C2) * D6/720;
    double LatRad = p - (N * tanp / R) * ZLat;

    double ZLon = 
      + D 
      - (1 + 2*T + C) * D3/6 
      + (5 - 2*C + 28*T - 3*C2 + 8*e2P + 24*T2) * D5/120;
    double LonRad = ZLon / cosp;
      
    (*Lat) = LatRad * rad2deg;
    (*Lon) = LonRad * rad2deg + refLon;
  }

//Northings sao positivos no norte e negativas no sul
//Easting sao positivas no leste e negativas no oeste

//Se a diferenca entre a longitude maxima e longitude minima for superior a 6 grau, o programa retorna 0
//indicando que nao será possível fazer a conversao
//Caso contrario, ele retornará o codigo 1

//converte coordenadas EUTM para coordenadas lat/long. Equaçoes provenientes do USGS Bulletin 1532
//Longitudes Leste são positivas, Longitudes Oeste são negativas.
//Latitudes Norte são positivas. Latitudes Sul são negativas
//Lat e Lon são dadas em graus decimais
