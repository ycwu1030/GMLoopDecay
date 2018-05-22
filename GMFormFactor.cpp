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
    std::ofstream BRS("GMBRScalars.dat");
    // std::ofstream BRH30("GMBRH30.dat");
    // std::ofstream BRH50("GMBRH50.dat");
    // std::ofstream BRH3p("GMBRH3p.dat");
    // std::ofstream BRH5p("GMBRH5p.dat");
    // std::ofstream BRH5pp("GMBRH5pp.dat");

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
    BRS<<"# mphi dM BRHAA BRHZA BRHZZ BRHWW GammaH BRHAAEFT BRHZAEFT BRHZZEFT BRHWWEFT GammaHEFT";
    BRS<<" BRH30AA BRH30ZA BRH30ZZ BRH30WW GammaH30 BRH30AAEFT BRH30ZAEFT BRH30ZZEFT BRH30WWEFT GammaH30EFT";
    BRS<<" BRH50AA BRH50ZA BRH50ZZ BRH50WW GammaH50 BRH50AAEFT BRH50ZAEFT BRH50ZZEFT BRH50WWEFT GammaH50EFT";
    BRS<<" BRH3pWA BRH3pWZ GammaH3p BRH3pWAEFT BRH3pWZEFT GammaH3pEFT";
    BRS<<" BRH5pWA BRH5pWZ GammaH5p BRH5pWAEFT BRH5pWZEFT GammaH5pEFT";
    BRS<<" BRH5ppWW GammaH5pp BRH5ppWWEFT GammaH5ppEFT"<<endl;
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
            BRS<<mphi<<"  "<<dM<<"  "<<GHAA/GHtot<<"  "<<GHZA/GHtot<<"  "<<GHZZ/GHtot<<"  "<<GHWW/GHtot<<"  "<<GHtot;
            GHAA = decay.GetGammaHVVEFT(H,A,A);
            GHZA = decay.GetGammaHVVEFT(H,Z,A);
            GHZZ = decay.GetGammaHVVEFT(H,Z,Z);
            GHWW = decay.GetGammaHVVEFT(H,W,W);
            GHtot = decay.GetGammaHtotEFT(H);
            BRS<<"  "<<GHAA/GHtot<<"  "<<GHZA/GHtot<<"  "<<GHZZ/GHtot<<"  "<<GHWW/GHtot<<"  "<<GHtot;

            GH30AA = decay.GetGammaHVV(H30,A,A);
            GH30ZA = decay.GetGammaHVV(H30,Z,A);
            GH30ZZ = decay.GetGammaHVV(H30,Z,Z);
            GH30WW = decay.GetGammaHVV(H30,W,W);
            GH30tot = decay.GetGammaHtot(H30);
            BRS<<"  "<<GH30AA/GH30tot<<"  "<<GH30ZA/GH30tot<<"  "<<GH30ZZ/GH30tot<<"  "<<GH30WW/GH30tot<<"  "<<GH30tot;
            GH30AA = decay.GetGammaHVVEFT(H30,A,A);
            GH30ZA = decay.GetGammaHVVEFT(H30,Z,A);
            GH30ZZ = decay.GetGammaHVVEFT(H30,Z,Z);
            GH30WW = decay.GetGammaHVVEFT(H30,W,W);
            GH30tot = decay.GetGammaHtotEFT(H30);
            BRS<<"  "<<GH30AA/GH30tot<<"  "<<GH30ZA/GH30tot<<"  "<<GH30ZZ/GH30tot<<"  "<<GH30WW/GH30tot<<"  "<<GH30tot;

            GH50AA = decay.GetGammaHVV(H50,A,A);
            GH50ZA = decay.GetGammaHVV(H50,Z,A);
            GH50ZZ = decay.GetGammaHVV(H50,Z,Z);
            GH50WW = decay.GetGammaHVV(H50,W,W);
            GH50tot = decay.GetGammaHtot(H50);
            BRS<<"  "<<GH50AA/GH50tot<<"  "<<GH50ZA/GH50tot<<"  "<<GH50ZZ/GH50tot<<"  "<<GH50WW/GH50tot<<"  "<<GH50tot;
            GH50AA = decay.GetGammaHVVEFT(H50,A,A);
            GH50ZA = decay.GetGammaHVVEFT(H50,Z,A);
            GH50ZZ = decay.GetGammaHVVEFT(H50,Z,Z);
            GH50WW = decay.GetGammaHVVEFT(H50,W,W);
            GH50tot = decay.GetGammaHtotEFT(H50);
            BRS<<"  "<<GH50AA/GH50tot<<"  "<<GH50ZA/GH50tot<<"  "<<GH50ZZ/GH50tot<<"  "<<GH50WW/GH50tot<<"  "<<GH50tot;

            GH3pWA = decay.GetGammaHVV(H3p,W,A);
            GH3pWZ = decay.GetGammaHVV(H3p,W,Z);
            GH3ptot = decay.GetGammaHtot(H3p);
            BRS<<"  "<<GH3pWA/GH3ptot<<"  "<<GH3pWZ/GH3ptot<<"  "<<GH3ptot;
            GH3pWA = decay.GetGammaHVVEFT(H3p,W,A);
            GH3pWZ = decay.GetGammaHVVEFT(H3p,W,Z);
            GH3ptot = decay.GetGammaHtotEFT(H3p);
            BRS<<"  "<<GH3pWA/GH3ptot<<"  "<<GH3pWZ/GH3ptot<<"  "<<GH3ptot;

            GH5pWA = decay.GetGammaHVV(H5p,W,A);
            GH5pWZ = decay.GetGammaHVV(H5p,W,Z);
            GH5ptot = decay.GetGammaHtot(H5p);
            BRS<<"  "<<GH5pWA/GH5ptot<<"  "<<GH5pWZ/GH5ptot<<"  "<<GH5ptot;
            GH5pWA = decay.GetGammaHVVEFT(H5p,W,A);
            GH5pWZ = decay.GetGammaHVVEFT(H5p,W,Z);
            GH5ptot = decay.GetGammaHtotEFT(H5p);
            BRS<<"  "<<GH5pWA/GH5ptot<<"  "<<GH5pWZ/GH5ptot<<"  "<<GH5ptot;

            GH5ppWW = decay.GetGammaHVV(H5pp,W,W);
            BRS<<"  "<<1.0<<"  "<<GH5ppWW;
            GH5ppWW = decay.GetGammaHVVEFT(H5pp,W,W);
            BRS<<"  "<<1.0<<"  "<<GH5ppWW<<endl;

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
