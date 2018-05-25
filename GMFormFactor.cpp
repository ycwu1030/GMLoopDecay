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
    // std::ofstream BRS("GMBRScalars.dat");
    std::ofstream BRH("GMBRH.dat");
    std::ofstream BRH3("GMBRH3.dat");
    std::ofstream BRH5("GMBRH5.dat");
    // std::ofstream BRH3p("GMBRH3p.dat");
    // std::ofstream BRH5p("GMBRH5p.dat");
    // std::ofstream BRH5pp("GMBRH5pp.dat");

    ltini();
    RealType mphi;
    RealType dM;
    double GHAA,GHZA,GHZZ,GHWW,GHH30Z,GHH3pW,GHtot;
    double GH30AA,GH30ZA,GH30ZZ,GH30WW,GH30HZ,GH30H5pW,GH30H50Z,GH30tot;
    double GH50AA,GH50ZA,GH50ZZ,GH50WW,GH50H3pW,GH50H30Z,GH50tot;
    double GH3pWA,GH3pWZ,GH3pHW,GH3pH5pZ,GH3pH5ppW,GH3pH50W,GH3ptot;
    double GH5pWA,GH5pWZ,GH5pH3pZ,GH5pH30W,GH5ptot;
    double GH5ppWW,GH5ppH3pW,GH5pptot;
    GMModel Mod;
    clock_t T1, T2, T3;
    GMDecay decay;
    // BRS<<"# mphi dM BRHAA BRHZA BRHZZ BRHWW GammaH BRHAAEFT BRHZAEFT BRHZZEFT BRHWWEFT GammaHEFT";
    // BRS<<" BRH30AA BRH30ZA BRH30ZZ BRH30WW GammaH30 BRH30AAEFT BRH30ZAEFT BRH30ZZEFT BRH30WWEFT GammaH30EFT";
    // BRS<<" BRH50AA BRH50ZA BRH50ZZ BRH50WW GammaH50 BRH50AAEFT BRH50ZAEFT BRH50ZZEFT BRH50WWEFT GammaH50EFT";
    // BRS<<" BRH3pWA BRH3pWZ GammaH3p BRH3pWAEFT BRH3pWZEFT GammaH3pEFT";
    // BRS<<" BRH5pWA BRH5pWZ GammaH5p BRH5pWAEFT BRH5pWZEFT GammaH5pEFT";
    // BRS<<" BRH5ppWW GammaH5pp BRH5ppWWEFT GammaH5ppEFT"<<endl;
// For H5:
    BRH5<<"# mphi dM BRH50AA BRH50ZA BRH50ZZ BRH50WW BRH50H3pW BRH50H30Z GammaH50";
    BRH5<<" BRH50AAEFT BRH50ZAEFT BRH50ZZEFT BRH50WWEFT BRH50H3pWEFT BRH50H30ZEFT GammaH50EFT";
    BRH5<<" BRH5pWA BRH5pWZ BRH5pH3pZ BRH5pH30W GammaH5p";
    BRH5<<" BRH5pWAEFT BRH5pWZEFT BRH5pH3pZ BRH5pH30W GammaH5pEFT";
    BRH5<<" BRH5ppWW BRH5ppH3pW GammaH5pp BRH5ppWWEFT BRH5ppH3pW GammaH5ppEFT"<<endl;
    for (mphi = 80.0; mphi<=500 ; mphi+=5)
    {
        for (dM = -500; dM <= 500; dM+=10)
        {
            if ((mphi*mphi + 3.0/2.0*abs(dM)*dM - 0.1 < 0))
            {
                continue;
            }
            Mod.MH = sqrt(mphi*mphi+3.0/2.0*abs(dM)*dM);
            Mod.MH3 = sqrt(mphi*mphi+abs(dM)*dM);
            Mod.MH5 = mphi;
            Mod.M2 = 80.0;
            decay.SetModel(Mod);

            GH50AA = decay.GetGammaHPartial(H50,A,A);
            GH50ZA = decay.GetGammaHPartial(H50,Z,A);
            GH50ZZ = decay.GetGammaHPartial(H50,Z,Z);
            GH50WW = decay.GetGammaHPartial(H50,W,W);
            GH50H3pW = decay.GetGammaHPartial(H50,H3p,W);
            GH50H30Z = decay.GetGammaHPartial(H50,H30,Z);
            GH50tot = decay.GetGammaHtot(H50);
            BRH5<<mphi<<"  "<<dM<<"  "<<GH50AA/GH50tot<<"  "<<GH50ZA/GH50tot<<"  "<<GH50ZZ/GH50tot<<"  "<<GH50WW/GH50tot;
            BRH5<<"  "<<GH50H3pW/GH50tot<<"  "<<GH50H30Z/GH50tot<<"  "<<GH50tot;
            GH50AA = decay.GetGammaHPartialEFT(H50,A,A);
            GH50ZA = decay.GetGammaHPartialEFT(H50,Z,A);
            GH50ZZ = decay.GetGammaHPartialEFT(H50,Z,Z);
            GH50WW = decay.GetGammaHPartialEFT(H50,W,W);
            GH50H3pW = decay.GetGammaHPartialEFT(H50,H3p,W);
            GH50H30Z = decay.GetGammaHPartialEFT(H50,H30,Z);
            GH50tot = decay.GetGammaHtotEFT(H50);
            BRH5<<"  "<<GH50AA/GH50tot<<"  "<<GH50ZA/GH50tot<<"  "<<GH50ZZ/GH50tot<<"  "<<GH50WW/GH50tot;
            BRH5<<"  "<<GH50H3pW/GH50tot<<"  "<<GH50H30Z/GH50tot<<"  "<<GH50tot;

            GH5pWA = decay.GetGammaHPartial(H5p,W,A);
            GH5pWZ = decay.GetGammaHPartial(H5p,W,Z);
            GH5pH3pZ = decay.GetGammaHPartial(H5p,H3p,Z);
            GH5pH30W = decay.GetGammaHPartial(H5p,H30,W);
            GH5ptot = decay.GetGammaHtot(H5p);
            BRH5<<"  "<<GH5pWA/GH5ptot<<"  "<<GH5pWZ/GH5ptot;
            BRH5<<"  "<<GH5pH3pZ/GH5ptot<<"  "<<GH5pH30W/GH5ptot<<"  "<<GH5ptot;
            GH5pWA = decay.GetGammaHPartialEFT(H5p,W,A);
            GH5pWZ = decay.GetGammaHPartialEFT(H5p,W,Z);
            GH5pH3pZ = decay.GetGammaHPartial(H5p,H3p,Z);
            GH5pH30W = decay.GetGammaHPartial(H5p,H30,W);
            GH5ptot = decay.GetGammaHtotEFT(H5p);
            BRH5<<"  "<<GH5pWA/GH5ptot<<"  "<<GH5pWZ/GH5ptot;
            BRH5<<"  "<<GH5pH3pZ/GH5ptot<<"  "<<GH5pH30W/GH5ptot<<"  "<<GH5ptot;

            GH5ppWW = decay.GetGammaHPartial(H5pp,W,W);
            GH5ppH3pW = decay.GetGammaHPartial(H5pp,H3p,W);
            GH5pptot = decay.GetGammaHtot(H5pp);
            BRH5<<"  "<<GH5ppWW/GH5pptot<<"  "<<GH5ppH3pW/GH5pptot<<"  "<<GH5pptot;
            GH5ppWW = decay.GetGammaHPartialEFT(H5pp,W,W);
            GH5ppH3pW = decay.GetGammaHPartialEFT(H5pp,H3p,W);
            GH5pptot = decay.GetGammaHtotEFT(H5pp);
            BRH5<<"  "<<GH5ppWW/GH5pptot<<"  "<<GH5ppH3pW/GH5pptot<<"  "<<GH5pptot<<endl;

        }
    }
// For H:
    BRH<<"# mphi dM BRHAA BRHZA BRHZZ BRHWW BRHH30Z BRHH3pW GammaH";
    BRH<<" BRHAAEFT BRHZAEFT BRHZZEFT BRHWWEFT BRHH30ZEFT BRHH3pWEFT GammaHEFT"<<endl;
    for (mphi = 80.0; mphi<=500 ; mphi+=5)
    {
        for (dM = -500; dM <= 500; dM+=10)
        {
            Mod.MH = mphi;
            if (mphi*mphi < 3.0/2.0*abs(dM)*dM + 0.1)
            {
                continue;
            }
            Mod.MH3 = sqrt(mphi*mphi-1.0/2.0*abs(dM)*dM);
            Mod.MH5 = sqrt(mphi*mphi-3.0/2.0*abs(dM)*dM);
            Mod.M2 = 80.0;
            decay.SetModel(Mod);
            GHAA = decay.GetGammaHPartial(H,A,A);
            GHZA = decay.GetGammaHPartial(H,Z,A);
            GHZZ = decay.GetGammaHPartial(H,Z,Z);
            GHWW = decay.GetGammaHPartial(H,W,W);
            GHH30Z = decay.GetGammaHPartial(H,H30,Z);
            GHH3pW = decay.GetGammaHPartial(H,H3p,W);
            GHtot = decay.GetGammaHtot(H);
            BRH<<mphi<<"  "<<dM<<"  "<<GHAA/GHtot<<"  "<<GHZA/GHtot<<"  "<<GHZZ/GHtot<<"  "<<GHWW/GHtot;
            BRH<<"  "<<GHH30Z/GHtot<<"  "<<GHH3pW/GHtot<<"  "<<GHtot;
            GHAA = decay.GetGammaHPartialEFT(H,A,A);
            GHZA = decay.GetGammaHPartialEFT(H,Z,A);
            GHZZ = decay.GetGammaHPartialEFT(H,Z,Z);
            GHWW = decay.GetGammaHPartialEFT(H,W,W);
            GHH30Z = decay.GetGammaHPartialEFT(H,H30,Z);
            GHH3pW = decay.GetGammaHPartialEFT(H,H3p,W);
            GHtot = decay.GetGammaHtotEFT(H);
            BRH<<"  "<<GHAA/GHtot<<"  "<<GHZA/GHtot<<"  "<<GHZZ/GHtot<<"  "<<GHWW/GHtot;
            BRH<<"  "<<GHH30Z/GHtot<<"  "<<GHH3pW/GHtot<<"  "<<GHtot<<endl;

        }
    }
// For H3:
    BRH3<<"# mphi dM BRH30AA BRH30ZA BRH30ZZ BRH30WW BRH30HZ BRH30H5pW BRH30H50Z GammaH30";
    BRH3<<" BRH30AAEFT BRH30ZAEFT BRH30ZZEFT BRH30WWEFT BRH30HZEFT BRH30H5pWEFT BRH30H50ZEFT GammaH30EFT";
    BRH3<<" BRH3pWA BRH3pWZ BRH3pHW BRH3pH5pZ BRH3pH5ppW BRH3pH50W GammaH3p";
    BRH3<<" BRH3pWAEFT BRH3pWZEFT BRH3pHWEFT BRH3pH5pZEFT BRH3pH5ppWEFT BRH3pH50WEFT GammaH3pEFT"<<endl;
    for (mphi = 80.0; mphi<=500 ; mphi+=5)
    {
        for (dM = -500; dM <= 500; dM+=10)
        {
            if ((mphi*mphi+1.0/2.0*abs(dM)*dM - 0.1 < 0) || (mphi*mphi - abs(dM)*dM - 0.1< 0))
            {
                continue;
            }
            Mod.MH3 = mphi;
            Mod.MH = sqrt(mphi*mphi+1.0/2.0*abs(dM)*dM);
            Mod.MH5 = sqrt(mphi*mphi-abs(dM)*dM);
            Mod.M2 = 80.0;
            decay.SetModel(Mod);
            GH30AA = decay.GetGammaHPartial(H30,A,A);
            GH30ZA = decay.GetGammaHPartial(H30,Z,A);
            GH30ZZ = decay.GetGammaHPartial(H30,Z,Z);
            GH30WW = decay.GetGammaHPartial(H30,W,W);
            GH30HZ = decay.GetGammaHPartial(H30,H,Z);
            GH30H5pW = decay.GetGammaHPartial(H30,H5p,W);
            GH30H50Z = decay.GetGammaHPartial(H30,H50,Z);
            GH30tot = decay.GetGammaHtot(H30);
            BRH3<<mphi<<"  "<<dM<<"  "<<GH30AA/GH30tot<<"  "<<GH30ZA/GH30tot<<"  "<<GH30ZZ/GH30tot<<"  "<<GH30WW/GH30tot;
            BRH3<<"  "<<GH30HZ/GH30tot<<"  "<<GH30H5pW/GH30tot<<"  "<<GH30H50Z/GH30tot<<"  "<<GH30tot;
            GH30AA = decay.GetGammaHPartialEFT(H30,A,A);
            GH30ZA = decay.GetGammaHPartialEFT(H30,Z,A);
            GH30ZZ = decay.GetGammaHPartialEFT(H30,Z,Z);
            GH30WW = decay.GetGammaHPartialEFT(H30,W,W);
            GH30HZ = decay.GetGammaHPartialEFT(H30,H,Z);
            GH30H5pW = decay.GetGammaHPartialEFT(H30,H5p,W);
            GH30H50Z = decay.GetGammaHPartialEFT(H30,H50,Z);
            GH30tot = decay.GetGammaHtotEFT(H30);
            BRH3<<"  "<<GH30AA/GH30tot<<"  "<<GH30ZA/GH30tot<<"  "<<GH30ZZ/GH30tot<<"  "<<GH30WW/GH30tot;
            BRH3<<"  "<<GH30HZ/GH30tot<<"  "<<GH30H5pW/GH30tot<<"  "<<GH30H50Z/GH30tot<<"  "<<GH30tot;

            GH3pWA = decay.GetGammaHPartial(H3p,W,A);
            GH3pWZ = decay.GetGammaHPartial(H3p,W,Z);
            GH3pHW = decay.GetGammaHPartial(H3p,H,W);
            GH3pH5pZ = decay.GetGammaHPartial(H3p,H5p,Z);
            GH3pH5ppW = decay.GetGammaHPartial(H3p,H5pp,W);
            GH3pH50W = decay.GetGammaHPartial(H3p,H50,W);
            GH3ptot = decay.GetGammaHtot(H3p);
            BRH3<<"  "<<GH3pWA/GH3ptot<<"  "<<GH3pWZ/GH3ptot;
            BRH3<<"  "<<GH3pHW/GH3ptot<<"  "<<GH3pH5pZ/GH3ptot<<"  "<<GH3pH5ppW/GH3ptot<<"  "<<GH3pH50W/GH3ptot<<"  "<<GH3ptot;
            GH3pWA = decay.GetGammaHPartialEFT(H3p,W,A);
            GH3pWZ = decay.GetGammaHPartialEFT(H3p,W,Z);
            GH3pHW = decay.GetGammaHPartialEFT(H3p,H,W);
            GH3pH5pZ = decay.GetGammaHPartialEFT(H3p,H5p,Z);
            GH3pH5ppW = decay.GetGammaHPartialEFT(H3p,H5pp,W);
            GH3pH50W = decay.GetGammaHPartialEFT(H3p,H50,W);
            GH3ptot = decay.GetGammaHtotEFT(H3p);
            BRH3<<"  "<<GH3pWA/GH3ptot<<"  "<<GH3pWZ/GH3ptot;
            BRH3<<"  "<<GH3pHW/GH3ptot<<"  "<<GH3pH5pZ/GH3ptot<<"  "<<GH3pH5ppW/GH3ptot<<"  "<<GH3pH50W/GH3ptot<<"  "<<GH3ptot<<endl;
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
