#include "GMHVVHHVChannel.h"

double lambda(double a, double b, double c)
{
    return a*a+b*b+c*c-2*a*b-2*b*c-2*c*a;
}

double lambdaij(double ki, double kj)
{
    return -1 + 2*ki + 2*kj - pow(ki-kj,2);
}

double deltaV(PID P1)
{
    if (P1 == W)
    {
        return 3.0/2.0;
    }
    else if (P1 == Z)
    {
        return 3.0*(7.0/12.0-10.0/9.0*sw*sw+40.0/27.0*sw*sw*sw*sw);
    }
    else
    {
        return 0;
    }
}

double Gij(double MH1, double MH2, double MV)
{
    double ki = pow(MH2/MH1,2);
    double kj = pow(MV/MH1,2);
    double lamij = lambdaij(ki,kj);
    return 1.0/4.0*(2*(-1+kj-ki)*sqrt(lamij)*(PI/2.0+atan((kj*(1-kj+ki)-lamij)/((1-ki)*sqrt(lamij))))+(lamij-2*ki)*log(ki)+1.0/3.0*(1-ki)*(5*(1+ki)-4*kj+2*lamij/kj));
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

double GammaHHVOFF_DIRECT(double MH1, double MH2, double MV, PID PIDV, double Coupling)
{
    double del = deltaV(PIDV);
    double G = Gij(MH1,MH2,MV);
    return del*3*Coupling*Coupling*MV*MV*MH1*G/(16*pow(PI,3)*v2);
}

// For VEGAS Integral
double GammaHVV_Integrand_GSL(double *Q2, size_t dim, void * modparams)
{
    GMModel_GSL *fp = (GMModel_GSL *) modparams;
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
    return result/PI2;
}

double GammaHVV_Flatten_Integrand_GSL(double *rho, size_t dim, void * modparams)
{
    GMModel_GSL *fp = (GMModel_GSL *) modparams;
    // Integral under Q2, not the rho
    double result;
    double Q1,Q2;
    Q1 = Q(rho[0],fp->MV1,fp->GAV1);
    Q2 = Q(rho[1],fp->MV2,fp->GAV2);
    if (rho[1]>Rho(fp->Mmother-Q1,fp->MV2,fp->GAV2))
    {
        result = 0.0;
    }
    else
    {
        result = pow(Q1*Q2/(fp->MV1*fp->MV2),2)*GammaHVV(fp->Mmother,Q1,Q2,fp->Mod,fp->S,fp->Stilde,fp->eta);
    }
    return result;
}

double GammaHHV_Integrand_GSL(double *Q2, size_t dim, void * modparams)
{
    GMModel_HHV_GSL *fp = (GMModel_HHV_GSL *) modparams;
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
    return result/PI2;
}

double GammaHHV_Flatten_Integrand_GSL(double *rho, size_t dim, void * modparams)
{
    GMModel_HHV_GSL *fp = (GMModel_HHV_GSL *) modparams;
    // Integral over the rho, not Q2
    double result;
    double Q1,Q2;
    Q1 = Q(rho[0],fp->MH2,fp->GAH2);
    Q2 = Q(rho[1],fp->MV,fp->GAV);
    if (rho[1]>Rho(fp->MH1-Q1,fp->MV,fp->GAV))
    {
        result = 0.0;
    }
    else
    {
        result = pow(Q1*Q2/(fp->MH2*fp->MV),2)*GammaHHV(fp->MH1,Q1,Q2,fp->CHHV);
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

    GMModel_GSL fp = {MH,MV1,GA1,MV2,GA2,FS,FStilde,eta,Mod};
    gsl_monte_function G = {&GammaHVV_Integrand_GSL, 2, &fp};

    size_t calls = GSLCALLS;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);

    gsl_monte_vegas_integrate(&G,Q2L,Q2U,2,calls/5,r,s,&res,&err);

    int tries = 0;
    do
    {
        gsl_monte_vegas_integrate(&G,Q2L,Q2U,2,calls/5,r,s,&res,&err);
        tries++;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.2&&tries<10);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    return res;

}

double GammaHVVOFF_Flatten_VEGAS(double MH, double MV1, double GA1, double MV2, double GA2, GMModel &Mod, SFunc FS, SFunc FStilde, double eta)
{
    if (MH > MV1 + MV2)
    {
        return GammaHVV(MH, MV1, MV2, Mod, FS, FStilde, eta);
    }

    double res, err;

    double rhoL[2] = {Rho(0,MV1,GA1),Rho(0,MV2,GA2)};
    double rhoU[2] = {Rho(MH,MV1,GA1),Rho(MH,MV2,GA2)};

    const gsl_rng_type *T;
    gsl_rng *r;

    GMModel_GSL fp = {MH,MV1,GA1,MV2,GA2,FS,FStilde,eta,Mod};
    gsl_monte_function G = {&GammaHVV_Flatten_Integrand_GSL, 2, &fp};

    size_t calls = GSLCALLS;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);

    gsl_monte_vegas_integrate(&G,rhoL,rhoU,2,calls/5,r,s,&res,&err);

    int tries = 0;
    do
    {
        gsl_monte_vegas_integrate(&G,rhoL,rhoU,2,calls/5,r,s,&res,&err);
        tries++;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.2&&tries<10);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    return res;

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

    GMModel_GSL fp = {MH,MV1,GA1,MV2,GA2,FS,FStilde,eta,Mod};
    gsl_monte_function G = {&GammaHVV_Integrand_GSL, 2, &fp};

    size_t calls = GSLCALLS;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_miser_state *s = gsl_monte_miser_alloc(2);

    gsl_monte_miser_integrate(&G,Q2L,Q2U,2,calls,r,s,&res,&err);

    gsl_monte_miser_free(s);
    gsl_rng_free(r);

    return res;

}

double GammaHVVOFF_Flatten_MISER(double MH, double MV1, double GA1, double MV2, double GA2, GMModel &Mod, SFunc FS, SFunc FStilde, double eta)
{
    if (MH > MV1 + MV2)
    {
        return GammaHVV(MH, MV1, MV2, Mod, FS, FStilde, eta);
    }

    double res, err;

    double rhoL[2] = {Rho(0,MV1,GA1),Rho(0,MV2,GA2)};
    double rhoU[2] = {Rho(MH,MV1,GA1),Rho(MH,MV2,GA2)};

    const gsl_rng_type *T;
    gsl_rng *r;

    GMModel_GSL fp = {MH,MV1,GA1,MV2,GA2,FS,FStilde,eta,Mod};
    gsl_monte_function G = {&GammaHVV_Flatten_Integrand_GSL, 2, &fp};

    size_t calls = GSLCALLS;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_miser_state *s = gsl_monte_miser_alloc(2);

    gsl_monte_miser_integrate(&G,rhoL,rhoU,2,calls,r,s,&res,&err);

    gsl_monte_miser_free(s);
    gsl_rng_free(r);

    return res;

}

double GammaHHVOFF_VEGAS(double MH1, double MH2, double GAH2, double MV, double GAV, double CHHV)
{
    if (MH1 > MH2 + MV)
    {
        return GammaHHV(MH1, MH2, MV, CHHV);
    }
    if (MH1 <= MH2)
    {
        return 0;
    }

    double res, err;

    double Q2L[2] = {0,0};
    double Q2U[2] = {MH1*MH1,MH1*MH1};

    const gsl_rng_type *T;
    gsl_rng *r;

    GMModel_HHV_GSL fp = {MH1,MH2,GAH2,MV,GAV,CHHV};
    gsl_monte_function G = {&GammaHHV_Integrand_GSL, 2, &fp};

    size_t calls = GSLCALLS;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);

    gsl_monte_vegas_integrate(&G,Q2L,Q2U,2,calls/5,r,s,&res,&err);

    int tries = 0;
    do
    {
        gsl_monte_vegas_integrate(&G,Q2L,Q2U,2,calls/5,r,s,&res,&err);
        tries++;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.2&&tries<10);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    return res;

}

double GammaHHVOFF_MISER(double MH1, double MH2, double GAH2, double MV, double GAV, double CHHV)
{
    if (MH1 > MH2 + MV)
    {
        return GammaHHV(MH1, MH2, MV, CHHV);
    }
    if (MH1 <= MH2)
    {
        return 0;
    }

    double res, err;

    double Q2L[2] = {0,0};
    double Q2U[2] = {MH1*MH1,MH1*MH1};

    const gsl_rng_type *T;
    gsl_rng *r;

    GMModel_HHV_GSL fp = {MH1,MH2,GAH2,MV,GAV,CHHV};
    gsl_monte_function G = {&GammaHHV_Integrand_GSL, 2, &fp};

    size_t calls = GSLCALLS;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_miser_state *s = gsl_monte_miser_alloc(2);

    gsl_monte_miser_integrate(&G,Q2L,Q2U,2,calls,r,s,&res,&err);

    gsl_monte_miser_free(s);
    gsl_rng_free(r);

    return res;

}

double GammaHHVOFF_Flatten_VEGAS(double MH1, double MH2, double GAH2, double MV, double GAV, double CHHV)
{
    if (MH1 > MH2 + MV)
    {
        return GammaHHV(MH1, MH2, MV, CHHV);
    }
    if (MH1 <= MH2)
    {
        return 0;
    }

    double res, err;

    double rhoL[2] = {Rho(0,MH2,GAH2),Rho(0,MV,GAV)};
    double rhoU[2] = {Rho(MH1,MH2,GAH2),Rho(MH1,MV,GAV)};

    const gsl_rng_type *T;
    gsl_rng *r;

    GMModel_HHV_GSL fp = {MH1,MH2,GAH2,MV,GAV,CHHV};
    gsl_monte_function G = {&GammaHHV_Flatten_Integrand_GSL, 2, &fp};

    size_t calls = GSLCALLS;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(2);

    gsl_monte_vegas_integrate(&G,rhoL,rhoU,2,calls/5,r,s,&res,&err);

    int tries = 0;
    do
    {
        gsl_monte_vegas_integrate(&G,rhoL,rhoU,2,calls/5,r,s,&res,&err);
        tries++;
    }
    while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.2&&tries<10);

    gsl_monte_vegas_free(s);
    gsl_rng_free(r);

    return res;

}

double GammaHHVOFF_Flatten_MISER(double MH1, double MH2, double GAH2, double MV, double GAV, double CHHV)
{
    if (MH1 > MH2 + MV)
    {
        return GammaHHV(MH1, MH2, MV, CHHV);
    }
    if (MH1 <= MH2)
    {
        return 0;
    }

    double res, err;

    double rhoL[2] = {Rho(0,MH2,GAH2),Rho(0,MV,GAV)};
    double rhoU[2] = {Rho(MH1,MH2,GAH2),Rho(MH1,MV,GAV)};

    const gsl_rng_type *T;
    gsl_rng *r;

    GMModel_HHV_GSL fp = {MH1,MH2,GAH2,MV,GAV,CHHV};
    gsl_monte_function G = {&GammaHHV_Flatten_Integrand_GSL, 2, &fp};

    size_t calls = GSLCALLS;
    gsl_rng_env_setup();

    T = gsl_rng_default;
    r = gsl_rng_alloc(T);

    gsl_monte_miser_state *s = gsl_monte_miser_alloc(2);

    gsl_monte_miser_integrate(&G,rhoL,rhoU,2,calls,r,s,&res,&err);

    gsl_monte_miser_free(s);
    gsl_rng_free(r);

    return res;

}