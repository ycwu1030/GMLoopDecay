#ifndef GMHVVHHVChannel_H
#define GMHVVHHVChannel_H

#include <iostream>
#include "GMModelParameters.h"
#include "GMFormFactorFunction.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_monte_miser.h>

typedef struct
{
    double Mmother;
    double MV1;
    double GAV1;
    double MV2;
    double GAV2;
    SFunc S;
    SFunc Stilde;
    double eta;

    GMModel Mod;
} GMModel_GSL;

typedef struct
{
    double MH1;
    double MH2;
    double GAH2;
    double MV;
    double GAV;
    double CHHV;
} GMModel_HHV_GSL;

double lambda(double a, double b, double c);

double GammaHVV(double MH, double MV1, double MV2, GMModel &Mod, SFunc S, SFunc Stilde, double eta);

double GammaHVV_Integrand_GSL(double *Q2, size_t dim, void * modparams);

double GammaHVV_Flatten_Integrand_GSL(double *rho, size_t dim, void * modparams);

double GammaHHV(double MH1, double MH2, double MV, double Coupling);

double GammaHHVOFF_DIRECT(double MH1, double MH2, double MV, PID PIDV, double Coupling);

double GammaHHV_Integrand_GSL(double *Q2, size_t dim, void * modparams);

double GammaHHV_Flatten_Integrand_GSL(double *rho, size_t dim, void * modparams);

// Following are the functions that will be used in Numerical integral 

double Rho(double Q, double M, double Gamma);

double Q(double rho, double M, double Gamma);

double BW(double Q2, double M, double Gamma);

double deltaV(PID P1);

double Gij(double MH1, double MH2, double MV);

// Following are the functions for the evaluation of off-shell decays

// This one uses naive grid/box approximation to perform the integral over the two off-shell resonance
double GammaHVVOFF(double MH, double MV1, double GA1, double MV2, double GA2, GMModel &Mod, SFunc S, SFunc Stilde, double eta);

// This is the same as above but for one of the final state is gamma
double GammaHVAOFF(double MH, double MV1, double GA1, GMModel &Mod, SFunc S, SFunc Stilde, double eta);

// This one uses VEGAS algorithm to perform the integral
double GammaHVVOFF_VEGAS(double MH, double MV1, double GA1, double MV2, double GA2, GMModel &Mod, SFunc FS, SFunc FStilde, double eta);
double GammaHVVOFF_MISER(double MH, double MV1, double GA1, double MV2, double GA2, GMModel &Mod, SFunc FS, SFunc FStilde, double eta);
double GammaHHVOFF_VEGAS(double MH1, double MH2, double GAH2, double MV, double GAV, double CHHV);
double GammaHHVOFF_MISER(double MH1, double MH2, double GAH2, double MV, double GAV, double CHHV);

double GammaHVVOFF_Flatten_VEGAS(double MH, double MV1, double GA1, double MV2, double GA2, GMModel &Mod, SFunc FS, SFunc FStilde, double eta);
double GammaHVVOFF_Flatten_MISER(double MH, double MV1, double GA1, double MV2, double GA2, GMModel &Mod, SFunc FS, SFunc FStilde, double eta);
double GammaHHVOFF_Flatten_VEGAS(double MH1, double MH2, double GAH2, double MV, double GAV, double CHHV);
double GammaHHVOFF_Flatten_MISER(double MH1, double MH2, double GAH2, double MV, double GAV, double CHHV);
#endif