#ifndef GMDECAY_H
#define GMDECAY_H

#include "GMModelParameters.h"
#include "GMFormFactorFunction.h"
#include "GMHVVHHVChannel.h"

enum INTEGRALFLAGS
{
    Q2VEGAS = 0,
    Q2MISER = 1,
    RHOVEGAS = 2,
    RHOMISER = 3,
    NONE = 4
};

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
    void SetIntegralMethod(INTEGRALFLAGS flag);
    double GetGammaHPartial(PID Mother, PID P1, PID P2);
    double GetGammaHtot(PID Mother);

    double GetGammaHPartialEFT(PID Mother, PID P1, PID P2);
    double GetGammaHtotEFT(PID Mother);
private: 
    GMModel _Mod;
    INTEGRALFLAGS _flag;
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

