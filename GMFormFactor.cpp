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
    std::ofstream BRH("GMBRH.dat");
    std::ofstream BRH30("GMBRH30.dat");
    std::ofstream BRH50("GMBRH50.dat");
    std::ofstream BRH3p("GMBRH3p.dat");
    std::ofstream BRH5p("GMBRH5p.dat");
    std::ofstream BRH5pp("GMBRH5pp.dat");

    ltini();
    RealType mphi;
    RealType dM;
    double GHAA,GHZA,GHZZ,GHWW, GHtot;
    double GH30AA,GH30ZA,GH30ZZ,GH30WW,GH30tot;
    double GH50AA,GH50ZA,GH50ZZ,GH50WW,GH50tot;
    double GH3pWA,GH3pWZ,GH3ptot;
    double GH5pWA,GH5pWZ,GH5ptot;
    double GH5ppWW;
    GMModel Mod;
    clock_t T1, T2, T3;
    GMDecay decay;
    for (mphi = 80.0; mphi<=500 ; mphi+=5)
    {
        for (dM = 0; dM < 200; dM+=10)
        {
            Mod.MH = sqrt(mphi*mphi+3.0/2.0*dM*dM);
            Mod.MH3 = sqrt(mphi*mphi+dM*dM);
            Mod.MH5 = mphi;
            Mod.M2 = 80.0;
            decay.SetModel(Mod);
            GHAA = decay.GetGammaHVV(H,A,A);
            GHZA = decay.GetGammaHVV(H,Z,A);
            GHZZ = decay.GetGammaHVV(H,Z,Z);
            GHWW = decay.GetGammaHVV(H,W,W);
            GHtot = decay.GetGammaHtot(H);
            BRH<<mphi<<"  "<<dM<<"  "<<GHAA/GHtot<<"  "<<GHZA/GHtot<<"  "<<GHZZ/GHtot<<"  "<<GHWW/GHtot<<"  "<<GHtot<<endl;

            GH30AA = decay.GetGammaHVV(H30,A,A);
            GH30ZA = decay.GetGammaHVV(H30,Z,A);
            GH30ZZ = decay.GetGammaHVV(H30,Z,Z);
            GH30WW = decay.GetGammaHVV(H30,W,W);
            GH30tot = decay.GetGammaHtot(H30);
            BRH30<<mphi<<"  "<<dM<<"  "<<GH30AA/GH30tot<<"  "<<GH30ZA/GH30tot<<"  "<<GH30ZZ/GH30tot<<"  "<<GH30WW/GH30tot<<"  "<<GH30tot<<endl;

            GH50AA = decay.GetGammaHVV(H50,A,A);
            GH50ZA = decay.GetGammaHVV(H50,Z,A);
            GH50ZZ = decay.GetGammaHVV(H50,Z,Z);
            GH50WW = decay.GetGammaHVV(H50,W,W);
            GH50tot = decay.GetGammaHtot(H50);
            BRH50<<mphi<<"  "<<dM<<"  "<<GH50AA/GH50tot<<"  "<<GH50ZA/GH50tot<<"  "<<GH50ZZ/GH50tot<<"  "<<GH50WW/GH50tot<<"  "<<GH50tot<<endl;

            GH3pWA = decay.GetGammaHVV(H3p,W,A);
            GH3pWZ = decay.GetGammaHVV(H3p,W,Z);
            GH3ptot = decay.GetGammaHtot(H3p);
            BRH3p<<mphi<<"  "<<dM<<"  "<<GH3pWA/GH3ptot<<"  "<<GH3pWZ/GH3ptot<<"  "<<GH3ptot<<endl;

            GH5pWA = decay.GetGammaHVV(H5p,W,A);
            GH5pWZ = decay.GetGammaHVV(H5p,W,Z);
            GH5ptot = decay.GetGammaHtot(H5p);
            BRH5p<<mphi<<"  "<<dM<<"  "<<GH5pWA/GH5ptot<<"  "<<GH5pWZ/GH5ptot<<"  "<<GH5ptot<<endl;

            GH5ppWW = decay.GetGammaHVV(H5pp,W,W);
            BRH5pp<<mphi<<"  "<<dM<<"  "<<1.0<<"  "<<GH5ppWW<<endl;

        }
    }
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
    ltexi();

    return 0;
}
