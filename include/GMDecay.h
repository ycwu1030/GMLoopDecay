#ifndef GMDECAY_H
#define GMDECAY_H

#include "GMModelParameters.h"
#include "GMFormFactorFunction.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_vegas.h>

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
} GMModel_VEGAS;

double lambda(double a, double b, double c);

double GammaHVV(double MH, double MV1, double MV2, GMModel &Mod, SFunc S, SFunc Stilde, double eta);

double GammaHVV_VEGAS(double *Q2, size_t dim, void * modparams);

// Following are the functions that will be used in Numerical integral 

double Rho(double Q, double M, double Gamma);

double Q(double rho, double M, double Gamma);

double BW(double Q2, double M, double Gamma);

// Following are the functions for the evaluation of off-shell decays

// This one uses naive grid/box approximation to perform the integral over the two off-shell resonance
double GammaHVVOFF(double MH, double MV1, double GA1, double MV2, double GA2, GMModel &Mod, SFunc S, SFunc Stilde, double eta);

// This is the same as above but for one of the final state is gamma
double GammaHVAOFF(double MH, double MV1, double GA1, GMModel &Mod, SFunc S, SFunc Stilde, double eta);

// This one uses VEGAS algorithm to perform the integral
double GammaHVVOFF_VEGAS(double MH, double MV1, double GA1, double MV2, double GA2, GMModel &Mod, SFunc FS, SFunc FStilde, double eta);

#define NChannel 18
class GMDecay
{
public:
    GMDecay();
    GMDecay(GMModel Mod);
    GMDecay(double MH, double MH3, double MH5, double M2);
    ~GMDecay(){};

    void SetModel(GMModel Mod);
    void SetModel(double MH, double MH3, double MH5, double M2);
    double GetGammaHVV(PID Mother, PID P1, PID P2);
    double GetGammaHtot(PID Mother);

    double GetGammaHVVEFT(PID Mother, PID P1, PID P2);
    double GetGammaHtotEFT(PID Mother);
private: 
    GMModel _Mod;
    void ResettingGammas();
    double GammaHVV[NChannel];
    double GammaHVVEFT[NChannel];
    int ToChannelID(PID Mother, PID P1, PID P2);
    void CalcGammaHVV(PID Mother, PID P1, PID P2);
    void CalcGammaHVVEFT(PID Mother, PID P1, PID P2);
    void GetMotherMass(PID Mother, double &mass);
    void GetVectorMassWidth(PID P1, double &mass, double &Gamma);
    enum ChannelID
    {
        HAA = 0,
        HWW = 1,
        HZA = 2,
        HZZ = 3,
        H30AA = 4,
        H30WW = 5,
        H30ZA = 6,
        H30ZZ = 7,
        H50AA = 8,
        H50WW = 9,
        H50ZA = 10,
        H50ZZ = 11,
        H3pWA = 12,
        H3pWZ = 13,
        H5pWA = 14,
        H5pWZ = 15,
        H5ppWW = 16,
        UNKOWN = 17
    };
/*
    // Gamma's for H 
    void CalcGammaHAA();
    double GammaHAA;

    void CalcGammaHWW();
    double GammaHWW;

    void CalcGammaHZA();
    double GammaHZA;

    void CalcGammaHZZ();
    double GammaHZZ;

    void CalcGammaHtot();
    double GammaHtot;

    // CalcGamma's for H30
    void CalcGammaH30AA();
    double GammaH30AA;

    void CalcGammaH30WW();
    double GammaH30WW;

    void CalcGammaH30ZA();
    double GammaH30ZA;

    void CalcGammaH30ZZ();
    double GammaH30ZZ;

    void CalcGammaH30tot();
    double GammaH30tot;

    // CalcGamma's for H50
    void CalcGammaH50AA();
    double GammaH50AA;

    void CalcGammaH50WW();
    double GammaH50WW;

    void CalcGammaH50ZA();
    double GammaH50ZA;

    void CalcGammaH50ZZ();
    double GammaH50ZZ;

    void CalcGammaH50tot();
    double GammaH50tot;

    // CalcGamma's for H3pm
    void CalcGammaH3pWA();
    double GammaH3pWA;

    void CalcGammaH3pWZ();
    double GammaH3pWZ;

    void CalcGammaH3ptot();
    double GammaH3ptot;

    // CalcGamma's for H5pm
    void CalcGammaH5pWA();
    double GammaH5pWA;

    void CalcGammaH5pWZ();
    double GammaH5pWZ;

    void CalcGammaH5ptot();
    double GammaH5ptot;

    // CalcGamma's for H5pp
    void CalcGammaH5ppWW();
    double GammaH5ppWW;

    void CalcGammaH5pptot();
    double GammaH5pptot;
*/
};
/*
// Gamma's for H 
double GammaHAA(GMModel &Mod);

double GammaHWW(GMModel &Mod);

double GammaHZA(GMModel &Mod);

double GammaHZZ(GMModel &Mod);

double GammaHtot(GMModel &Mod);

// Gamma's for H30
double GammaH30AA(GMModel &Mod);

double GammaH30WW(GMModel &Mod);

double GammaH30ZA(GMModel &Mod);

double GammaH30ZZ(GMModel &Mod);

double GammaH30tot(GMModel &Mod);

// Gamma's for H50
double GammaH50AA(GMModel &Mod);

double GammaH50WW(GMModel &Mod);

double GammaH50ZA(GMModel &Mod);

double GammaH50ZZ(GMModel &Mod);

double GammaH50tot(GMModel &Mod);

// Gamma's for H3pm
double GammaH3pWA(GMModel &Mod);

double GammaH3pWZ(GMModel &Mod);

// Gamma's for H5pm
double GammaH5pWA(GMModel &Mod);

double GammaH5pWZ(GMModel &Mod);

// Gamma's for H5pp
double GammaH5ppWW(GMModel &Mod);
*/
#endif

