#include <iostream>
#include <fstream>
#include <time.h>
#include "GMModelParameters.h"
#include "clooptools.h"
#include "GMDecay.h"
#include "GMFormFactorFunction.h"

using namespace GM;
using namespace std;

int main() {
    std::ofstream outfile("GMtest.dat");
    ltini();
    RealType mphi;
    double GammaOSWZ,GammaOSWA;
    double GammaAA,GammaZA,GammaZZ,GammaWW;
    double Gammatot;
    GMModel Mod;
    clock_t T1, T2, T3;
    GMDecay decay;
    for (double i = log10(100); i <= log10(1000); i=i+0.05)
    {
        mphi = pow(10,i);
        Mod.MH = mphi;
        Mod.MH3 = mphi;
        Mod.MH5 = mphi;
        Mod.M2 = 80.0;
        decay.SetModel(Mod);
        GammaAA = decay.GetGammaHVV(H50,A,A);
        GammaZA = decay.GetGammaHVV(H50,Z,A);
        GammaZZ = decay.GetGammaHVV(H50,Z,Z);
        GammaWW = decay.GetGammaHVV(H50,W,W);
        Gammatot = decay.GetGammaHtot(H50);
        outfile<<mphi<<"  "<<GammaAA<<"  "<<GammaZA<<"  "<<GammaZZ<<"  "<<GammaWW<<"  "<<Gammatot<<endl;
        // T1 = clock();
        // GammaWZ = GammaHVVOFF(Mod.MH5,MW,GAW,MZ,GAZ,Mod,SH5pWZ,NOIMPLEMENTED,1);
        // T2 = clock();
        // GammaWZVEGAS = GammaHVVOFF_VEGAS(Mod.MH5,MW,GAW,MZ,GAZ,Mod,SH5pWZ,NOIMPLEMENTED,1);
        // T3 = clock();
        // std::cout<<(double)(T2-T1)/CLOCKS_PER_SEC<<"  "<<(double)(T3-T2)/CLOCKS_PER_SEC<<std::endl;
        // Gammatot = GammaWA+GammaWZ;
        // BRWA = GammaWA/Gammatot;
        // BRWZ = GammaWZ/Gammatot;
        // std::cout<<GammaOSWZ<<"  "<<GammaWZ<<"  "<<GammaWZVEGAS<<std::endl;
        // std::cout<<SWA<<"  "<<GammaOSWA<<"  "<<GammaWA<<std::endl;
    }
    ltexi();

    return 0;
}
