#include "GMDecay.h"
#include <iostream>


double lambda(double a, double b, double c)
{
    return a*a+b*b+c*c-2*a*b-2*b*c-2*c*a;
}

double GammaHVV(double Mm, double MV1, double MV2, GMModel &Mod, SFunc FS, SFunc FStilde, double eta)
{
    double Mm2 = Mm*Mm;
    double MV12 = MV1*MV1;
    double MV22 = MV2*MV2;
    double MH2 = Mod.MH*Mod.MH;
    double MH32 = Mod.MH3*Mod.MH3;
    double MH52 = Mod.MH5*Mod.MH5;
    double lam = lambda(Mm2,MV12,MV22);
    // std::cout<<lam<<std::endl;
    double S = std::abs(FS(MH2,MH32,MH52,Mod.M2,Mm2,MV12,MV22));
    double Stilde = std::abs(FStilde(MH2,MH32,MH52,Mod.M2,Mm2,MV12,MV22));
    clearcache(); // This is used to clear the cache produced by LoopTools, otherwise the Lookup table will be too large if you loop over some parameters
    return (pow(lam,1.5)*(S*S+Stilde*Stilde)+6*MV12*MV22*pow(lam,0.5)*S*S)/(32*PI*pow(Mm,3)*eta);
}

// For VEGAS Integral
double GammaHVV_Integrand_VEGAS(double *Q2, size_t dim, void * modparams)
{
    GMModel_VEGAS *fp = (GMModel_VEGAS *) modparams;
    // Integral under Q2, not the rho
    double result;
    if (Q2[1] >= pow(fp->Mmother-sqrt(Q2[0]),2))
    {
        result = 0.0;
    }
    else
    {
        result = BW(Q2[0],fp->MV1,fp->GAV1)*BW(Q2[1],fp->MV2,fp->GAV2)*GammaHVV(fp->Mmother,sqrt(Q2[0]),sqrt(Q2[1]),fp->Mod,fp->S,fp->Stilde,fp->eta);
    }
    return result;
}

// Following are the functions that will be used in Numerical integral 

double Rho(double Q, double M, double Gamma)
{
    return 1.0/PI*atan((Q*Q-M*M)/(M*Gamma));
}

double Q(double rho, double M, double Gamma)
{
    return sqrt(M*Gamma*tan(PI*rho)+M*M);
}

double BW(double Q2, double M, double Gamma)
{
    return Q2*Gamma/M/(pow(Q2-M*M,2)+pow(M*Gamma,2));
}

double GammaHVVOFF(double MH, double MV1, double GA1, double MV2, double GA2, GMModel &Mod, SFunc FS, SFunc FStilde, double eta)
{
    int steps;
    double R1I, R1F, R1F1, R1F2, R2I, R2F, R2F1, R2F2, R1, R2;
    double interval1, interval2;
    double GammaFinal=0.0;
    double SLICE;
    double Q1,Q2;
    double PIECE1,PIECE2;
    if (MH > MV1 + MV2)
    {
        steps = 100;

        R1I = Rho(0,MV1,GA1);
        R1F = Rho(MH,MV1,GA1);
        R2I = Rho(0,MV2,GA2);
        interval1 = (R1F-R1I)/steps;
        for (int i = 0; i < steps; ++i)
        {
            R1 = R1I + (i+0.5)*interval1;
            Q1 = Q(R1,MV1,GA1);
            R2F = Rho(MH-Q1,MV2,GA2);
            SLICE = 0.0;
            interval2 = (R2F-R2I)/steps;
            for (int j = 0; j < steps; ++j)
            {
                R2 = R2I + (j+0.5)*interval2;
                Q2 = Q(R2,MV2,GA2);
                SLICE += interval2*pow(Q1*Q2/(MV1*MV2),2)*GammaHVV(MH,Q1,Q2,Mod,FS,FStilde,eta);
            }
            GammaFinal += SLICE*interval1;
        }
    }
    else
    {
        steps = OFFSHELLSTEPS;
        // The Integral region is divided into two parts for better numerical behavior
        // First region:
        R1I = Rho(0,MV1,GA1);
        R1F = Rho(MH,MV1,GA1);
        R2I = Rho(0,MV2,GA2);
        PIECE1 = 0.0;
        interval1 = (R1F-R1I)/steps;
#ifdef DEBUG
        std::cout<<"Entering Loop of the first piece in off-shell calculation"<<std::endl;
#endif
        for (int i = 0; i < steps; ++i)
        {
            R1 = R1I + (i+0.5)*interval1;
            Q1 = Q(R1,MV1,GA1);
            R2F1 = Rho(MH-Q1,MV2,GA2);
            R2F2 = Rho(MV2*Q1/MV1,MV2,GA2);
            R2F = R2F1<R2F2?R2F1:R2F2; // This condition separate the whole region into two.
            SLICE = 0.0;
            interval2 = (R2F-R2I)/steps;
            for (int j = 0; j < steps; ++j)
            {
                R2 = R2I + (j+0.5)*interval2;
                Q2 = Q(R2,MV2,GA2);
                SLICE+=interval2*pow(Q1*Q2/(MV1*MV2),2)*GammaHVV(MH,Q1,Q2,Mod,FS,FStilde,eta);
#ifdef DEBUG
                // if ((i+1)%100==0&&(j+1)%100==0)
                {
                    std::cout<<"Looping i: "<<i<<" j: "<<j<<"\r";
                }
#endif
            }
            PIECE1 += SLICE*interval1;
        }
#ifdef DEBUG
        std::cout<<std::endl;
        std::cout<<"Exit Loop of the first piece in off-shell calculation;"<<std::endl;
#endif

        // Second region:
        R2I = Rho(0,MV2,GA2);
        R2F = Rho(MH,MV2,GA2);
        R1I = Rho(0,MV1,GA1);
        PIECE2 = 0.0;
        interval2 = (R2F-R2I)/steps;
        for (int i = 0; i < steps; ++i)
        {
            R2 = R2I + (i+0.5)*interval2;
            Q2 = Q(R2,MV2,GA2);
            R1F1 = Rho(MH-Q2,MV1,GA1);
            R1F2 = Rho(MV1*Q2/MV2,MV1,GA1);
            R1F = R1F1<R1F2?R1F1:R1F2;
            SLICE = 0.0;
            interval1 = (R1F-R1I)/steps;
            for (int j = 0; j < steps; ++j)
            {
                R1 = R1I + (j+0.5)*interval1;
                Q1 = Q(R1,MV1,GA1);
                SLICE += interval1*pow(Q1*Q2/(MV1*MV2),2)*GammaHVV(MH,Q1,Q2,Mod,FS,FStilde,eta);
            }
            PIECE2 += interval2*SLICE;
        }
        GammaFinal = PIECE1 + PIECE2;
    }
    return GammaFinal;
}

// For off-shell decay involving one gamma final states (The gamma can be treated as always on-shell)
double GammaHVAOFF(double MH, double MV1, double GA1, GMModel &Mod, SFunc FS, SFunc FStilde, double eta)
{
    if (MH > MV1)
    {
        return GammaHVV(MH, MV1, 0, Mod, FS, FStilde, eta);
    }
    int steps;
    double R1I, R1F, R1;
    double interval1;
    double GammaFinal=0.0;
    double SLICE;
    double Q1;
    steps = OFFSHELLSTEPS;
    R1I = Rho(0,MV1,GA1);
    R1F = Rho(MH,MV1,GA1);
    interval1 = (R1F-R1I)/steps;
    SLICE = 0.0;
    for (int i = 0; i < steps; ++i)
    {
        R1 = R1I + (i+0.5)*interval1;
        Q1 = Q(R1,MV1,GA1);
        SLICE += interval1*pow(Q1/MV1,2)*GammaHVV(MH,Q1,0,Mod,FS,FStilde,eta);
    }
    GammaFinal = SLICE;
    return GammaFinal;
}

double GammaHVVOFF_VEGAS(double MH, double MV1, double GA1, double MV2, double GA2, GMModel &Mod, SFunc FS, SFunc FStilde, double eta)
{
    if (MH > MV1 + MV2)
    {
        return GammaHVV(MH, MV1, MV2, Mod, FS, FStilde, eta);
    }

    double res, err;

    double Q2L[2] = {0,0};
    double Q2U[2] = {MH*MH,MH*MH};

    const gsl_rng_type *T;
    gsl_rng *r;

    GMModel_VEGAS fp = {MH,MV1,GA1,MV2,GA2,FS,FStilde,eta,Mod};
    gsl_monte_function G = {&GammaHVV_Integrand_VEGAS, 2, &fp};

    size_t calls = 1000;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);

    gsl_monte_vegas_integrate(&G,Q2L,Q2U,2,calls,r,s,&res,&err);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    return res/PI2;

}

GMDecay::GMDecay()
{
    _Mod.MH = 200.0;
    _Mod.MH3 = 200.0;
    _Mod.MH5 = 200.0;
    _Mod.M2 = 80.0;
    ResettingGammas();
}

GMDecay::GMDecay(GMModel Mod)
{
    _Mod = Mod;
    ResettingGammas();
}

GMDecay::GMDecay(double MH, double MH3, double MH5, double M2)
{
    _Mod.MH = MH;
    _Mod.MH3 = MH3;
    _Mod.MH5 = MH5;
    _Mod.M2 = M2;
    ResettingGammas();
}

void GMDecay::SetModel(GMModel Mod)
{
    _Mod = Mod;
    ResettingGammas();
}

void GMDecay::SetModel(double MH, double MH3, double MH5, double M2)
{
    _Mod.MH = MH;
    _Mod.MH3 = MH3;
    _Mod.MH5 = MH5;
    _Mod.M2 = M2;
    ResettingGammas();
}

void GMDecay::ResettingGammas()
{
/*
    GammaHAA = -1;
    GammaHWW = -1;
    GammaHZA = -1;
    GammaHZZ = -1;
    GammaHtot = -1;
    GammaH30AA = -1;
    GammaH30WW = -1;
    GammaH30ZA = -1;
    GammaH30ZZ = -1;
    GammaH30tot = -1;
    GammaH50AA = -1;
    GammaH50WW = -1;
    GammaH50ZA = -1;
    GammaH50ZZ = -1;
    GammaH50tot = -1;
    GammaH3pWA = -1;
    GammaH3pWZ = -1;
    GammaH3ptot = -1;
    GammaH5pWA = -1;
    GammaH5pWZ = -1;
    GammaH5ptot = -1;
    GammaH5ppWW = -1;
    GammaH5pptot = -1;
*/
    for (int i = 0; i < NChannel; ++i)
    {
        GammaHVV[i] = -1;
    }
    GammaHVV[NChannel-1] = 0;
}

int GMDecay::ToChannelID(PID Mother, PID P1, PID P2)
{
    if (Mother == H)
    {
        if (P1 == A && P2 == A)
        {
          return HAA;
        }
        else if (P1 == W && P2 == W)
        {
          return HWW;
        }
        else if ( (P1 == Z && P2 == A) || (P1 == A && P2 == Z) )
        {
          return HZA;
        }
        else if ( P1 == Z && P2 == Z)
        {
          return HZZ;
        }
        else
        {
          std::cout<<"Warning: Didn't Find the desired H channel, using UNKOWN instead"<<std::endl;
          return UNKOWN;
        }
    }
    else if (Mother == H30)
    {
        if (P1 == A && P2 == A)
        {
          return H30AA;
        }
        else if (P1 == W && P2 == W)
        {
          return H30WW;
        }
        else if ( (P1 == Z && P2 == A) || (P1 == A && P2 == Z) )
        {
          return H30ZA;
        }
        else if ( P1 == Z && P2 == Z)
        {
          return H30ZZ;
        }
        else
        {
          std::cout<<"Warning: Didn't Find the desired H30 channel, using UNKOWN instead"<<std::endl;
          return UNKOWN;
        }
    }
    else if (Mother == H50)
    {
        if (P1 == A && P2 == A)
        {
          return H50AA;
        }
        else if (P1 == W && P2 == W)
        {
          return H50WW;
        }
        else if ( (P1 == Z && P2 == A) || (P1 == A && P2 == Z) )
        {
          return H50ZA;
        }
        else if ( P1 == Z && P2 == Z)
        {
          return H50ZZ;
        }
        else
        {
          std::cout<<"Warning: Didn't Find the desired H50 channel, using UNKOWN instead"<<std::endl;
          return UNKOWN;
        }
    }
    else if (Mother == H3p)
    {
        if ((P1 == W && P2 == A) || (P1 == A && P2 == W))
        {
          return H3pWA;
        }
        else if ((P1 == W && P2 == Z) || (P1 == Z && P2 == W))
        {
          return H3pWZ;
        }
        else
        {
          std::cout<<"Warning: Didn't Find the desired H3p channel, using UNKOWN instead"<<std::endl;
          return UNKOWN;
        }
    }
    else if (Mother == H5p)
    {
        if ((P1 == W && P2 == A) || (P1 == A && P2 == W))
        {
          return H5pWA;
        }
        else if ((P1 == W && P2 == Z) || (P1 == Z && P2 == W))
        {
          return H5pWZ;
        }
        else
        {
          std::cout<<"Warning: Didn't Find the desired H5p channel, using UNKOWN instead"<<std::endl;
          return UNKOWN;
        }
    }
    else if (Mother == H5pp)
    {
        if (P1 == W && P2 == W)
        {
          return H5ppWW;
        }
        else
        {
          std::cout<<"Warning: Didn't Find the desired H5pp channel, using UNKOWN instead"<<std::endl;
          return UNKOWN;
        }
    }
    else
    {
        std::cout<<"Warning: Didn't Find the Mother Particle, using UNKOWN instead"<<std::endl;
        return UNKOWN;
    }
}

void GMDecay::GetMotherMass(PID Mother, double &mass)
{
    if (Mother == H)
    {
        mass = _Mod.MH;
    }
    else if (Mother == H3p || Mother == H30)
    {
        mass = _Mod.MH3;
    }
    else if (Mother == H5pp || Mother == H5p || Mother == H50)
    {
        mass = _Mod.MH5;
    }
    else
    {
        std::cout<<"Warning: No corresponding Particle, just giving 200 GeV"<<std::endl;
        mass = 200.0;
    }
}

void GMDecay::GetVectorMassWidth(PID P1, double &mass, double &Gamma)
{
    if (P1 == A)
    {
        mass = 0;
        Gamma = -1;
    }
    else if (P1 == Z)
    {
        mass = MZ;
        Gamma = GAZ;
    }
    else if (P1 == W)
    {
        mass = MW;
        Gamma = GAW;
    }
    else
    {
        std::cout<<"Warning: No corresponding Vector Boson, just give the Z boson"<<std::endl;
        mass = MZ;
        Gamma = GAZ;
    }

}

void GMDecay::CalcGammaHVV(PID Mother, PID P1, PID P2)
{
    int channel = ToChannelID(Mother, P1, P2);
    if (GammaHVV[channel] > 0 || channel == UNKOWN)
    {
        return;
    }
    SFunc SHVV = GetSFunc(Mother, P1, P2);
    SFunc STilde = GM::NOIMPLEMENTED;

    double eta = P1==P2?2:1;
    double m1, m2, m3, ga2, ga3;
    GetMotherMass(Mother,m1);
    GetVectorMassWidth(P1,m2,ga2);
    GetVectorMassWidth(P2,m3,ga3);
    if ((ga3 < 0 && ga2 > 0)||(ga3 < 0 && ga2 < 0))
    {
        GammaHVV[channel] = GammaHVAOFF(m1,m2,ga2,_Mod,SHVV,STilde,eta);
    }
    else if (ga2 < 0 && ga3 > 0)
    {
        GammaHVV[channel] = GammaHVAOFF(m1,m3,ga3,_Mod,SHVV,STilde,eta);
    }
    else
    {
        GammaHVV[channel] = GammaHVVOFF_VEGAS(m1,m2,ga2,m3,ga3,_Mod,SHVV,STilde,eta);
    }
}

double GMDecay::GetGammaHVV(PID Mother, PID P1, PID P2)
{
    int channel = ToChannelID(Mother, P1, P2);
    if (GammaHVV[channel] > 0 || channel == UNKOWN)
    {
        return GammaHVV[channel];
    }
    else
    {
        CalcGammaHVV(Mother, P1, P2);
    }
    return GammaHVV[channel];
}

double GMDecay::GetGammaHtot(PID Mother)
{
    double Gammatot=0.0;
    if (Mother == H)
    {
        Gammatot += GetGammaHVV(H,A,A);
        Gammatot += GetGammaHVV(H,W,W);
        Gammatot += GetGammaHVV(H,Z,A);
        Gammatot += GetGammaHVV(H,Z,Z);
    }
    else if (Mother == H30)
    {
        Gammatot += GetGammaHVV(H30,A,A);
        Gammatot += GetGammaHVV(H30,W,W);
        Gammatot += GetGammaHVV(H30,Z,A);
        Gammatot += GetGammaHVV(H30,Z,Z);
    }
    else if (Mother == H50)
    {
        Gammatot += GetGammaHVV(H50,A,A);
        Gammatot += GetGammaHVV(H50,W,W);
        Gammatot += GetGammaHVV(H50,Z,A);
        Gammatot += GetGammaHVV(H50,Z,Z);
    }
    else if (Mother == H3p)
    {
        Gammatot += GetGammaHVV(H3p,W,A);
        Gammatot += GetGammaHVV(H3p,W,Z);
    }
    else if (Mother == H5p)
    {
        Gammatot += GetGammaHVV(H5p,W,A);
        Gammatot += GetGammaHVV(H5p,W,Z);
    }
    else if (Mother == H5pp)
    {
        Gammatot += GetGammaHVV(H5pp,W,W);
    }
    return Gammatot;
}

























