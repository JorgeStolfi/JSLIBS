// Test problems for univariate minimizers
// Last edited on 2023-02-19 13:02:12 by stolfi

#ifndef minu_test_problems_H
#define minu_test_problems_H

#define _GNU_SOURCE
#include <stdio.h>

#include <minu_test.h>

#define QuadXStart (0.0)
#define QuadDist (1.0)

#define QuadXMin (-10.0)
#define QuadXMax (+10.0)
#define QuadYMin ( -1.0)
#define QuadYMax (+151.0)

void QuadEval(void *prb, double x, double *fx, double *dfx);
double QuadError (void *prb, double x, double fx);

#define BentXStart (0.0)
#define BentDist (1.0)

#define BentXMin (-10.0)
#define BentXMax (+10.0)
#define BentYMin ( -1.0)
#define BentYMax (+151.0)

void BentEval(void *prb, double x, double *fx, double *dfx);
double BentError (void *prb, double x, double fx);
     
#define CornXStart (0.0)
#define CornDist (1.0)

#define CornXMin (-10.0)
#define CornXMax (+10.0)
#define CornYMin ( -1.0)
#define CornYMax (+15.1)

void CornEval(void *prb, double x, double *fx, double *dfx);
double CornError (void *prb, double x, double fx);

#define WavyXStart (0.0)
#define WavyDist (1.0)

#define WavyXMin (-10.0)
#define WavyXMax (+10.0)
#define WavyYMin ( -1.0)
#define WavyYMax (+151.0)
     
void WavyEval(void *prb, double x, double *fx, double *dfx);
double WavyError (void *prb, double x, double fx);

#define BiquXStart (0.0)
#define BiquDist (1.0)

#define BiquXMin (-10.0)
#define BiquXMax (+10.0)
#define BiquYMin ( -1000.0)
#define BiquYMax (+15001.0)
            
void BiquEval (void *prb, double x, double *fx, double *dfx);
double BiquError (void *prb, double x, double fx);

#define HoleXStart (0.0)
#define HoleDist (1.0)

#define HoleXMin (-10.0)
#define HoleXMax (+10.0)
#define HoleYMin ( -1511.0)
#define HoleYMax (+100.0)

void HoleEval (void *prb, double x, double *fx, double *dfx);
double HoleError (void *prb, double x, double fx);

#endif
