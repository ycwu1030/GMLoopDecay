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
    ltini();
    RealType mphi = 250;
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
    for (dM = 0; dM < 250; dM+=5)
    {
        Mod.MH = mphi;
        if (mphi*mphi< 3.0/2.0*dM*dM) continue;
        Mod.MH3 = sqrt(mphi*mphi-1.0/2.0*dM*dM);
        Mod.MH5 = sqrt(mphi*mphi-3.0/2.0*dM*dM);
        Mod.M2 = 80.0;
        decay.SetModel(Mod);
        // GHAA = decay.GetGammaHPartial(H,A,A);
        // GHZA = decay.GetGammaHPartial(H,Z,A);
        // GHZZ = decay.GetGammaHPartial(H,Z,Z);
        // GHWW = decay.GetGammaHPartial(H,W,W);
        GHH30Z = decay.GetGammaHPartial(H,H30,Z);
        GHH3pW = decay.GetGammaHPartial(H,H3p,W);
        GHtot = decay.GetGammaHtot(H);
        GH30tot = decay.GetGammaHtot(H30);
        GH3ptot = decay.GetGammaHtot(H3p);
        cout<<dM<<"  "<<Mod.MH<<"  "<<Mod.MH3<<"  "<<Mod.MH5<<"  "<<GHH30Z<<"  "<<GHH3pW<<"  "<<GHtot<<"  "<<GH30tot<<"  "<<GH3ptot<<endl;   
    }
    ltexi();

    return 0;
}
