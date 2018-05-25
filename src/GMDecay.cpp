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

double GammaHHV(double MH1, double MH2, double MV, double Coupling)
{
    double MH12 = MH1*MH1;
    double MH22 = MH2*MH2;
    double MV2 = MV*MV;
    double lam = lambda(MH12,MH22,MV2);
    return Coupling*Coupling*pow(lam,1.5)/(16*PI*pow(MH1,3)*MV2);
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

double GammaHHV_Integrand_VEGAS(double *Q2, size_t dim, void * modparams)
{
    GMModel_HHV_VEGAS *fp = (GMModel_HHV_VEGAS *) modparams;
    // Integral under Q2, not the rho
    double result;
    if (Q2[1] >= pow(fp->MH1-sqrt(Q2[0]),2))
    {
        result = 0.0;
    }
    else
    {
        result = BW(Q2[0],fp->MH2,fp->GAH2)*BW(Q2[1],fp->MV,fp->GAV)*GammaHHV(fp->MH1,sqrt(Q2[0]),sqrt(Q2[1]),fp->CHHV);
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

    size_t calls = GSLCALLS;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);

    gsl_monte_vegas_integrate(&G,Q2L,Q2U,2,calls,r,s,&res,&err);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    return res/PI2;

}

double GammaHVVOFF_MISER(double MH, double MV1, double GA1, double MV2, double GA2, GMModel &Mod, SFunc FS, SFunc FStilde, double eta)
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

    size_t calls = GSLCALLS;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_miser_state *s = gsl_monte_miser_alloc(2);

    gsl_monte_miser_integrate(&G,Q2L,Q2U,2,calls,r,s,&res,&err);

    gsl_monte_miser_free(s);
    gsl_rng_free(r);

    return res/PI2;

}

double GammaHHVOFF_VEGAS(double MH1, double MH2, double GAH2, double MV, double GAV, double CHHV)
{
    if (MH1 > MH2 + MV)
    {
        return GammaHHV(MH1, MH2, MV, CHHV);
    }
    if (MH1 <= MH2 + 0.1)
    {
        return 0;
    }

    double res, err;

    double Q2L[2] = {0,0};
    double Q2U[2] = {MH1*MH1,MH1*MH1};

    const gsl_rng_type *T;
    gsl_rng *r;

    GMModel_HHV_VEGAS fp = {MH1,MH2,GAH2,MV,GAV,CHHV};
    gsl_monte_function G = {&GammaHHV_Integrand_VEGAS, 2, &fp};

    size_t calls = GSLCALLS;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);

    gsl_monte_vegas_integrate(&G,Q2L,Q2U,2,calls,r,s,&res,&err);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    return res/PI2;

}

double GammaHHVOFF_MISER(double MH1, double MH2, double GAH2, double MV, double GAV, double CHHV)
{
    if (MH1 > MH2 + MV)
    {
        return GammaHHV(MH1, MH2, MV, CHHV);
    }
    if (MH1 <= MH2 + 0.1)
    {
        return 0;
    }

    double res, err;

    double Q2L[2] = {0,0};
    double Q2U[2] = {MH1*MH1,MH1*MH1};

    const gsl_rng_type *T;
    gsl_rng *r;

    GMModel_HHV_VEGAS fp = {MH1,MH2,GAH2,MV,GAV,CHHV};
    gsl_monte_function G = {&GammaHHV_Integrand_VEGAS, 2, &fp};

    size_t calls = GSLCALLS;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_miser_state *s = gsl_monte_miser_alloc(2);

    gsl_monte_miser_integrate(&G,Q2L,Q2U,2,calls,r,s,&res,&err);

    gsl_monte_miser_free(s);
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
    for (int i = 0; i < NChannel; ++i)
    {
        GammaHPartial[i] = -1;
        GammaHPartialEFT[i] = -1;
    }
    GammaHPartial[NChannel-1] = 0;
    GammaHPartialEFT[NChannel-1] = 0;
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
        else if ((P1 == H30 && P2 == Z) || (P1 == Z && P2 == H30))
        {
            return HH30Z;
        }
        else if ((P1 == H3p && P2 == W) || (P1 == W && P2 == H3p))
        {
            return HH3pW;
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
        else if ((P1 == H && P2 == Z) || (P1 == Z && P2 == H))
        {
            return H30HZ;
        }
        else if ((P1 == H5p && P2 == W) || (P1 == W && P2 == H5p))
        {
            return H30H5pW;
        }
        else if ((P1 == H50 && P2 == Z) || (P1 == Z && P2 == H50))
        {
            return H30H50Z;
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
        else if ((P1 == H3p && P2 == W) || (P1 == W && P2 == H3p))
        {
            return H50H3pW;
        }
        else if ((P1 == H30 && P2 == Z) || (P1 == Z && P2 == H30))
        {
            return H50H30Z;
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
        else if ((P1 == H && P2 == W) || (P1 == W && P2 == H))
        {
            return H3pHW;
        }
        else if ((P1 == H5p && P2 == Z) || (P1 == Z && P2 == H5p))
        {
            return H3pH5pZ;
        }
        else if ((P1 == H5pp && P2 == W) || (P1 == W && P2 == H5pp))
        {
            return H3pH5ppW;
        }
        else if ((P1 == H50 && P2 == W) || (P1 == W && P2 == H50))
        {
            return H3pH50W;
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
        else if ((P1 == H3p && P2 == Z) || (P1 == Z && P2 == H3p))
        {
            return H5pH3pZ;
        }
        else if ((P1 == H30 && P2 == W) || (P1 == W && P2 == H30))
        {
            return H5pH30W;
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
        else if ((P1 == H3p && P2 == W) || (P1 == W && P2 == H3p))
        {
            return H5ppH3pW;
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

void GMDecay::GetScalarMass(PID Mother, double &mass)
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
    if (GammaHPartial[channel] > 0 || channel == UNKOWN)
    {
        return;
    }
    SFunc SHVV = GetSFunc(Mother, P1, P2);
    SFunc STilde = GM::NOIMPLEMENTED;

    double eta = P1==P2?2:1;
    double m1, m2, m3, ga2, ga3;
    GetScalarMass(Mother,m1);
    GetVectorMassWidth(P1,m2,ga2);
    GetVectorMassWidth(P2,m3,ga3);
    if ((ga3 < 0 && ga2 > 0)||(ga3 < 0 && ga2 < 0))
    {
        GammaHPartial[channel] = GammaHVAOFF(m1,m2,ga2,_Mod,SHVV,STilde,eta);
    }
    else if (ga2 < 0 && ga3 > 0)
    {
        GammaHPartial[channel] = GammaHVAOFF(m1,m3,ga3,_Mod,SHVV,STilde,eta);
    }
    else
    {
        GammaHPartial[channel] = GammaHVVOFF_VEGAS(m1,m2,ga2,m3,ga3,_Mod,SHVV,STilde,eta);
    }
}

void GMDecay::CalcGammaHHV(PID Mother, PID P1, PID P2)
{
    PID temp;
    if (P2>P1)
    {
        temp = P1;
        P1 = P2;
        P2 = temp;
    }
    int channel = ToChannelID(Mother, P1, P2);
    if (GammaHPartial[channel] > 0 || channel == UNKOWN)
    {
        return;
    }
    double CHHV = GetHHVCoupling(Mother,P1,P2);
    double eta = GetMulti(Mother, P1, P2);
    double m1, m2, m3, ga2, ga3;
    GetScalarMass(Mother,m1);
    GetScalarMass(P1,m2);
    if (m1 < m2 + 0.1)
    {
        GammaHPartial[channel] = 0.0;
        return;
    }
    ga2 = GetGammaHtot(P1);
    GetVectorMassWidth(P2,m3,ga3);
    GammaHPartial[channel] = eta*GammaHHVOFF_MISER(m1,m2,ga2,m3,ga3,CHHV);
}

void GMDecay::CalcGammaHVVEFT(PID Mother, PID P1, PID P2)
{
    int channel = ToChannelID(Mother, P1, P2);
    if (GammaHPartialEFT[channel] > 0 || channel == UNKOWN)
    {
        return;
    }
    SFunc SHVV = GM::SH50AA; // Always using SH50AA
    double ratio = GetSFuncEFTRatio(Mother,P1,P2);
    SFunc STilde = GM::NOIMPLEMENTED;

    double eta = P1==P2?2:1;
    double m1, m2, m3, ga2, ga3;
    GetScalarMass(Mother,m1);
    GetVectorMassWidth(P1,m2,ga2);
    GetVectorMassWidth(P2,m3,ga3);
    if ((ga3 < 0 && ga2 > 0)||(ga3 < 0 && ga2 < 0))
    {
        GammaHPartialEFT[channel] = ratio*ratio*GammaHVAOFF(m1,m2,ga2,_Mod,SHVV,STilde,eta);
    }
    else if (ga2 < 0 && ga3 > 0)
    {
        GammaHPartialEFT[channel] = ratio*ratio*GammaHVAOFF(m1,m3,ga3,_Mod,SHVV,STilde,eta);
    }
    else
    {
        GammaHPartialEFT[channel] = ratio*ratio*GammaHVVOFF_VEGAS(m1,m2,ga2,m3,ga3,_Mod,SHVV,STilde,eta);
    }
}

void GMDecay::CalcGammaHHVEFT(PID Mother, PID P1, PID P2)
{
    PID temp;
    if (P2>P1)
    {
        temp = P1;
        P1 = P2;
        P2 = temp;
    }
    int channel = ToChannelID(Mother, P1, P2);
    if (GammaHPartialEFT[channel] > 0 || channel == UNKOWN)
    {
        return;
    }
    else if (GammaHPartial[channel] > 0)
    {
        GammaHPartialEFT[channel] = GammaHPartial[channel];
        return;
    }
    double CHHV = GetHHVCoupling(Mother,P1,P2);
    double eta = GetMulti(Mother, P1, P2);
    double m1, m2, m3, ga2, ga3;
    GetScalarMass(Mother,m1);
    GetScalarMass(P1,m2);
    if (m1 < m2 + 0.1)
    {
        GammaHPartialEFT[channel] = 0.0;
        return;
    }
    ga2 = GetGammaHtot(P1);
    GetVectorMassWidth(P2,m3,ga3);
    GammaHPartialEFT[channel] = eta*GammaHHVOFF_MISER(m1,m2,ga2,m3,ga3,CHHV);
}

double GMDecay::GetGammaHPartial(PID Mother, PID P1, PID P2)
{
    int channel = ToChannelID(Mother, P1, P2);
    if (GammaHPartial[channel] > 0 || channel == UNKOWN)
    {
        return GammaHPartial[channel];
    }
    if (P1 < 30 && P2 < 30)
    {
        CalcGammaHVV(Mother, P1, P2);
    }
    else
    {
        // std::cout<<"Calculating: HHV"<<Mother<<"  "<<P1<<" "<<P2<<std::endl;
        CalcGammaHHV(Mother, P1, P2);
    }
    return GammaHPartial[channel];
}

double GMDecay::GetGammaHPartialEFT(PID Mother, PID P1, PID P2)
{
    int channel = ToChannelID(Mother, P1, P2);
    if (GammaHPartialEFT[channel] > 0 || channel == UNKOWN)
    {
        return GammaHPartialEFT[channel];
    }
    if (P1 < 30 && P2 < 30)
    {
        CalcGammaHVVEFT(Mother, P1, P2);
    }
    else
    {
        CalcGammaHHVEFT(Mother, P1, P2);
    }
    return GammaHPartialEFT[channel];
}

double GMDecay::GetGammaHtot(PID Mother)
{
    double Gammatot=0.0;
    if (Mother == H)
    {
        Gammatot += GetGammaHPartial(H,A,A);
        Gammatot += GetGammaHPartial(H,W,W);
        Gammatot += GetGammaHPartial(H,Z,A);
        Gammatot += GetGammaHPartial(H,Z,Z);
        Gammatot += GetGammaHPartial(H,H30,Z);
        Gammatot += GetGammaHPartial(H,H3p,W);
    }
    else if (Mother == H30)
    {
        Gammatot += GetGammaHPartial(H30,A,A);
        Gammatot += GetGammaHPartial(H30,W,W);
        Gammatot += GetGammaHPartial(H30,Z,A);
        Gammatot += GetGammaHPartial(H30,Z,Z);
        Gammatot += GetGammaHPartial(H30,H,Z);
        Gammatot += GetGammaHPartial(H30,H5p,W);
        Gammatot += GetGammaHPartial(H30,H50,Z);
    }
    else if (Mother == H50)
    {
        Gammatot += GetGammaHPartial(H50,A,A);
        Gammatot += GetGammaHPartial(H50,W,W);
        Gammatot += GetGammaHPartial(H50,Z,A);
        Gammatot += GetGammaHPartial(H50,Z,Z);
        Gammatot += GetGammaHPartial(H50,H3p,W);
        Gammatot += GetGammaHPartial(H50,H30,Z);
    }
    else if (Mother == H3p)
    {
        Gammatot += GetGammaHPartial(H3p,W,A);
        Gammatot += GetGammaHPartial(H3p,W,Z);
        Gammatot += GetGammaHPartial(H3p,H,W);
        Gammatot += GetGammaHPartial(H3p,H5p,Z);
        Gammatot += GetGammaHPartial(H3p,H5pp,W);
        Gammatot += GetGammaHPartial(H3p,H50,W);
    }
    else if (Mother == H5p)
    {
        Gammatot += GetGammaHPartial(H5p,W,A);
        Gammatot += GetGammaHPartial(H5p,W,Z);
        Gammatot += GetGammaHPartial(H5p,H3p,Z);
        Gammatot += GetGammaHPartial(H5p,H30,W);
    }
    else if (Mother == H5pp)
    {
        Gammatot += GetGammaHPartial(H5pp,W,W);
        Gammatot += GetGammaHPartial(H5pp,H3p,W);
    }
    return Gammatot;
}

double GMDecay::GetGammaHtotEFT(PID Mother)
{
    double Gammatot=0.0;
    if (Mother == H)
    {
        Gammatot += GetGammaHPartialEFT(H,A,A);
        Gammatot += GetGammaHPartialEFT(H,W,W);
        Gammatot += GetGammaHPartialEFT(H,Z,A);
        Gammatot += GetGammaHPartialEFT(H,Z,Z);
        Gammatot += GetGammaHPartialEFT(H,H30,Z);
        Gammatot += GetGammaHPartialEFT(H,H3p,W);
    }
    else if (Mother == H30)
    {
        Gammatot += GetGammaHPartialEFT(H30,A,A);
        Gammatot += GetGammaHPartialEFT(H30,W,W);
        Gammatot += GetGammaHPartialEFT(H30,Z,A);
        Gammatot += GetGammaHPartialEFT(H30,Z,Z);
        Gammatot += GetGammaHPartialEFT(H30,H,Z);
        Gammatot += GetGammaHPartialEFT(H30,H5p,W);
        Gammatot += GetGammaHPartialEFT(H30,H50,Z);
    }
    else if (Mother == H50)
    {
        Gammatot += GetGammaHPartialEFT(H50,A,A);
        Gammatot += GetGammaHPartialEFT(H50,W,W);
        Gammatot += GetGammaHPartialEFT(H50,Z,A);
        Gammatot += GetGammaHPartialEFT(H50,Z,Z);
        Gammatot += GetGammaHPartialEFT(H50,H3p,W);
        Gammatot += GetGammaHPartialEFT(H50,H30,Z);
    }
    else if (Mother == H3p)
    {
        Gammatot += GetGammaHPartialEFT(H3p,W,A);
        Gammatot += GetGammaHPartialEFT(H3p,W,Z);
        Gammatot += GetGammaHPartialEFT(H3p,H,W);
        Gammatot += GetGammaHPartialEFT(H3p,H5p,Z);
        Gammatot += GetGammaHPartialEFT(H3p,H5pp,W);
        Gammatot += GetGammaHPartialEFT(H3p,H50,W);
    }
    else if (Mother == H5p)
    {
        Gammatot += GetGammaHPartialEFT(H5p,W,A);
        Gammatot += GetGammaHPartialEFT(H5p,W,Z);
        Gammatot += GetGammaHPartialEFT(H5p,H3p,Z);
        Gammatot += GetGammaHPartialEFT(H5p,H30,W);
    }
    else if (Mother == H5pp)
    {
        Gammatot += GetGammaHPartialEFT(H5pp,W,W);
        Gammatot += GetGammaHPartialEFT(H5pp,H3p,W);
    }
    return Gammatot;
}

double GMDecay::GetMulti(PID Mother, PID P1, PID P2)
{
    double multi = 1;
    if (Mother == H && ((P1 == H3p && P2 == W) || (P1 == W && P2 == H3p)))
    {
        multi = 2;
    }
    else if (Mother == H30 && ((P1 == H5p && P2 == W) || (P1 == W && P2 == H5p)))
    {
        multi = 2;
    }
    else if (Mother == H50 && ((P1 == H3p && P2 == W) || (P1 == W && P2 == H3p)))
    {
        multi = 2;
    }
    return multi;
}
























