#ifndef GMDECAY_H
#define GMDECAY_H

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
} GMModel_VEGAS;

typedef struct
{
    double MH1;
    double MH2;
    double GAH2;
    double MV;
    double GAV;
    double CHHV;
} GMModel_HHV_VEGAS;

double lambda(double a, double b, double c);

double GammaHVV(double MH, double MV1, double MV2, GMModel &Mod, SFunc S, SFunc Stilde, double eta);

double GammaHVV_VEGAS(double *Q2, size_t dim, void * modparams);

double GammaHHV(double MH1, double MH2, double MV, double Coupling);

double GammaHHV_VEGAS(double *Q2, size_t dim, void * modparams);

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

#define NChannel 32
class GMDecay
{
public:
    GMDecay();
    GMDecay(GMModel Mod);
    GMDecay(double MH, double MH3, double MH5, double M2);
    ~GMDecay(){};

    void SetModel(GMModel Mod);
    void SetModel(double MH, double MH3, double MH5, double M2);
    double GetGammaHPartial(PID Mother, PID P1, PID P2);
    double GetGammaHtot(PID Mother);

    double GetGammaHPartialEFT(PID Mother, PID P1, PID P2);
    double GetGammaHtotEFT(PID Mother);
private: 
    GMModel _Mod;
    void ResettingGammas();
    double GammaHPartial[NChannel];
    double GammaHPartialEFT[NChannel];
    int ToChannelID(PID Mother, PID P1, PID P2);
    void CalcGammaHVV(PID Mother, PID P1, PID P2);
    void CalcGammaHHV(PID Mother, PID P1, PID P2);
    void CalcGammaHVVEFT(PID Mother, PID P1, PID P2);
    void CalcGammaHHVEFT(PID Mother, PID P1, PID P2);// Just a copy of CalcGammaHHV, there is not EFT version of this calculation
    void GetScalarMass(PID P1, double &mass);
    void GetVectorMassWidth(PID P1, double &mass, double &Gamma);
    double GetMulti(PID Mother, PID P1, PID P2); // Check if we need to multiply by 2 for Hi->HjV channel;
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
        HH30Z = 17,
        HH3pW = 18,
        H3pHW = 19,
        H3pH5pZ = 20,
        H3pH5ppW = 21,
        H3pH50W = 22,
        H30HZ = 23,
        H30H5pW = 24,
        H30H50Z = 25,
        H5pH3pZ = 26,
        H5pH30W = 27,
        H5ppH3pW = 28,
        H50H3pW = 29,
        H50H30Z = 30,
        UNKOWN = 31
    };
};

#endif

