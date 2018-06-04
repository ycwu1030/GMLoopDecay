#include "GMDecay.h"
#include <iostream>

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

void GMDecay::SetIntegralMethod(INTEGRALFLAGS flag)
{
    _flag = flag;
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
    if (GammaHPartial[channel] >= 0 || channel == UNKOWN)
    {
        return;
    }
    SFunc SHVV = GetSFunc(Mother, P1, P2);
    SFunc STilde = GetSTildeFunc(Mother, P1, P2);

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
        if (_flag == Q2VEGAS)
        {
            GammaHPartial[channel] = GammaHVVOFF_VEGAS(m1,m2,ga2,m3,ga3,_Mod,SHVV,STilde,eta);
        }
        else if (_flag == RHOVEGAS)
        {
            GammaHPartial[channel] = GammaHVVOFF_Flatten_VEGAS(m1,m2,ga2,m3,ga3,_Mod,SHVV,STilde,eta);
        }
        else if (_flag == Q2MISER)
        {
            GammaHPartial[channel] = GammaHVVOFF_MISER(m1,m2,ga2,m3,ga3,_Mod,SHVV,STilde,eta);
        }
        else if (_flag == RHOMISER)
        {
            GammaHPartial[channel] = GammaHVVOFF_Flatten_MISER(m1,m2,ga2,m3,ga3,_Mod,SHVV,STilde,eta);
        }
        else if (_flag == NONE)
        {
            GammaHPartial[channel] = GammaHVVOFF_VEGAS(m1,m2,ga2,m3,ga3,_Mod,SHVV,STilde,eta);
        }
        else
        {
            std::cout<<"Warning: Unknown Integral flags, suing Q2VEGAS instead"<<std::endl;
            GammaHPartial[channel] = GammaHVVOFF_VEGAS(m1,m2,ga2,m3,ga3,_Mod,SHVV,STilde,eta);
        }
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
    if (GammaHPartial[channel] >= 0 || channel == UNKOWN)
    {
        return;
    }
    double CHHV = GetHHVCoupling(Mother,P1,P2);
    double eta = GetMulti(Mother, P1, P2);
    double m1, m2, m3, ga2, ga3;
    GetScalarMass(Mother,m1);
    GetScalarMass(P1,m2);
    if (m1 <= m2)
    {
        GammaHPartial[channel] = 0.0;
        return;
    }
    ga2 = GetGammaHtot(P1);
    GetVectorMassWidth(P2,m3,ga3);
    if (_flag == Q2VEGAS)
    {
        GammaHPartial[channel] = eta*GammaHHVOFF_VEGAS(m1,m2,ga2,m3,ga3,CHHV);
    }
    else if (_flag == RHOVEGAS)
    {
        GammaHPartial[channel] = eta*GammaHHVOFF_Flatten_VEGAS(m1,m2,ga2,m3,ga3,CHHV);
    }
    else if (_flag == Q2MISER)
    {
        GammaHPartial[channel] = eta*GammaHHVOFF_MISER(m1,m2,ga2,m3,ga3,CHHV);
    }
    else if (_flag == RHOMISER)
    {
        GammaHPartial[channel] = eta*GammaHHVOFF_Flatten_MISER(m1,m2,ga2,m3,ga3,CHHV);
    }
    else if (_flag == NONE)
    {
        if (m1 > m2 + m3)
        {
            GammaHPartial[channel] = eta*GammaHHV(m1,m2,m3,CHHV);
        }
        else
        {
            GammaHPartial[channel] = eta*GammaHHVOFF_DIRECT(m1,m2,m3,P2,CHHV);
        }
    }
    else
    {
        std::cout<<"Warning: Unknown Integral flags, suing Q2VEGAS instead"<<std::endl;
        GammaHPartial[channel] = eta*GammaHHVOFF_VEGAS(m1,m2,ga2,m3,ga3,CHHV);
    }
}

void GMDecay::CalcGammaHVVEFT(PID Mother, PID P1, PID P2)
{
    int channel = ToChannelID(Mother, P1, P2);
    if (GammaHPartialEFT[channel] >= 0 || channel == UNKOWN)
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
    if (GammaHPartialEFT[channel] >= 0 || channel == UNKOWN)
    {
        return;
    }
    else if (GammaHPartial[channel] >= 0)
    {
        GammaHPartialEFT[channel] = GammaHPartial[channel];
        return;
    }
    double CHHV = GetHHVCoupling(Mother,P1,P2);
    double eta = GetMulti(Mother, P1, P2);
    double m1, m2, m3, ga2, ga3;
    GetScalarMass(Mother,m1);
    GetScalarMass(P1,m2);
    if (m1 <= m2)
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
    if (GammaHPartial[channel] >= 0 || channel == UNKOWN)
    {
        return GammaHPartial[channel];
    }
    if (P1 < 30 && P2 < 30)
    {
        CalcGammaHVV(Mother, P1, P2);
    }
    else
    {
        CalcGammaHHV(Mother, P1, P2);
    }
    return GammaHPartial[channel];
}

double GMDecay::GetGammaHPartialEFT(PID Mother, PID P1, PID P2)
{
    int channel = ToChannelID(Mother, P1, P2);
    if (GammaHPartialEFT[channel] >= 0 || channel == UNKOWN)
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
























