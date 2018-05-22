#include "GMFormFactorFunction.h"
#include <iostream>

SFunc GetSFunc(PID Mother, PID P1, PID P2)
{
  if (Mother == H)
  {
    if (P1 == A && P2 == A)
    {
      return GM::SHAA;
    }
    else if (P1 == W && P2 == W)
    {
      return GM::SHWW;
    }
    else if ( (P1 == Z && P2 == A) || (P1 == A && P2 == Z) )
    {
      return GM::SHZA;
    }
    else if ( P1 == Z && P2 == Z)
    {
      return GM::SHZZ;
    }
    else
    {
      std::cout<<"Warning: Didn't Find the desired channel, using HZZ instead"<<std::endl;
      return GM::SHZZ;
    }
  }
  else if (Mother == H30)
  {
    if (P1 == A && P2 == A)
    {
      return GM::SH30AA;
    }
    else if (P1 == W && P2 == W)
    {
      return GM::SH30WW;
    }
    else if ( (P1 == Z && P2 == A) || (P1 == A && P2 == Z) )
    {
      return GM::SH30ZA;
    }
    else if ( P1 == Z && P2 == Z)
    {
      return GM::SH30ZZ;
    }
    else
    {
      std::cout<<"Warning: Didn't Find the desired channel, using H30ZZ instead"<<std::endl;
      return GM::SH30ZZ;
    }
  }
  else if (Mother == H50)
  {
    if (P1 == A && P2 == A)
    {
      return GM::SH50AA;
    }
    else if (P1 == W && P2 == W)
    {
      return GM::SH50WW;
    }
    else if ( (P1 == Z && P2 == A) || (P1 == A && P2 == Z) )
    {
      return GM::SH50ZA;
    }
    else if ( P1 == Z && P2 == Z)
    {
      return GM::SH50ZZ;
    }
    else
    {
      std::cout<<"Warning: Didn't Find the desired channel, using H50ZZ instead"<<std::endl;
      return GM::SH50ZZ;
    }
  }
  else if (Mother == H3p)
  {
    if ((P1 == W && P2 == A) || (P1 == A && P2 == W))
    {
      return GM::SH3pWA;
    }
    else if ((P1 == W && P2 == Z) || (P1 == Z && P2 == W))
    {
      return GM::SH3pWZ;
    }
    else
    {
      std::cout<<"Warning: Didn't Find the desired channel, using H3pWZ instead"<<std::endl;
      return GM::SH3pWZ;
    }
  }
  else if (Mother == H5p)
  {
    if ((P1 == W && P2 == A) || (P1 == A && P2 == W))
    {
      return GM::SH5pWA;
    }
    else if ((P1 == W && P2 == Z) || (P1 == Z && P2 == W))
    {
      return GM::SH5pWZ;
    }
    else
    {
      std::cout<<"Warning: Didn't Find the desired channel, using H5pWZ instead"<<std::endl;
      return GM::SH5pWZ;
    }
  }
  else if (Mother == H5pp)
  {
    if (P1 == W && P2 == W)
    {
      return GM::SH5ppWW;
    }
    else
    {
      std::cout<<"Warning: Didn't Find the desired channel, using H5ppWW instead"<<std::endl;
      return GM::SH5ppWW;
    }
  }
  else
  {
    std::cout<<"Warning: Didn't Find the desired channel, using HWW instead"<<std::endl;
    return GM::SHWW;
  }
}

double GetSFuncEFTRatio(PID Mother, PID P1, PID P2)
{
  if (Mother == H)
  {
    if (P1 == A && P2 == A)
    {
      return 1.0/(sqrt(2));
    }
    else if (P1 == W && P2 == W)
    {
      return 0;
    }
    else if ( (P1 == Z && P2 == A) || (P1 == A && P2 == Z) )
    {
      return (1.0-2.0*sw*sw)/(2.0*sqrt(2)*sw*cw);
    }
    else if ( P1 == Z && P2 == Z)
    {
      return (-1.0/sqrt(2));
    }
    else
    {
      std::cout<<"Warning: Didn't Find the desired channel, using HZZ instead"<<std::endl;
      return (-1.0/sqrt(2));
    }
  }
  else if (Mother == H30)
  {
    if (P1 == A && P2 == A)
    {
      return 0;
    }
    else if (P1 == W && P2 == W)
    {
      return 0;
    }
    else if ( (P1 == Z && P2 == A) || (P1 == A && P2 == Z) )
    {
      return 0;
    }
    else if ( P1 == Z && P2 == Z)
    {
      return 0;
    }
    else
    {
      std::cout<<"Warning: Didn't Find the desired channel, using H30ZZ instead"<<std::endl;
      return 0;
    }
  }
  else if (Mother == H50)
  {
    if (P1 == A && P2 == A)
    {
      return 1.0;
    }
    else if (P1 == W && P2 == W)
    {
      return 0;
    }
    else if ( (P1 == Z && P2 == A) || (P1 == A && P2 == Z) )
    {
      return ((1.0-2.0*sw*sw)/(2*sw*cw));
    }
    else if ( P1 == Z && P2 == Z)
    {
      return -1.0;
    }
    else
    {
      std::cout<<"Warning: Didn't Find the desired channel, using H50ZZ instead"<<std::endl;
      return -1.0;
    }
  }
  else if (Mother == H3p)
  {
    if ((P1 == W && P2 == A) || (P1 == A && P2 == W))
    {
      return (sqrt(3)/(4.0*cw));
    }
    else if ((P1 == W && P2 == Z) || (P1 == Z && P2 == W))
    {
      return (-sqrt(3)/(4.0*cw));
    }
    else
    {
      std::cout<<"Warning: Didn't Find the desired channel, using H3pWZ instead"<<std::endl;
      return (-sqrt(3)/(4.0*cw));
    }
  }
  else if (Mother == H5p)
  {
    if ((P1 == W && P2 == A) || (P1 == A && P2 == W))
    {
      return (-sqrt(3)/(4.0*cw));
    }
    else if ((P1 == W && P2 == Z) || (P1 == Z && P2 == W))
    {
      return (sqrt(3)/(4.0*cw));
    }
    else
    {
      std::cout<<"Warning: Didn't Find the desired channel, using H5pWZ instead"<<std::endl;
      return (sqrt(3)/(4.0*cw));
    }
  }
  else if (Mother == H5pp)
  {
    if (P1 == W && P2 == W)
    {
      return 0;
    }
    else
    {
      std::cout<<"Warning: Didn't Find the desired channel, using H5ppWW instead"<<std::endl;
      return 0;
    }
  }
  else
  {
    std::cout<<"Warning: Didn't Find the desired channel, using HWW instead"<<std::endl;
    return 0;
  }
}

namespace GM{
ComplexType NOIMPLEMENTED(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
  return 0.0;
}

ComplexType SH30AA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return 0
;
}

ComplexType SH30WW(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return (2.*Alfa*M2coeff*(1.*C0i(cc1, m12, m22, m32, MH2, MH32, MH32) + 
   3.75*C0i(cc1, m12, m22, m32, MH32, MH52, MH52) + 
   1.*C0i(cc11, m12, m22, m32, MH2, MH32, MH32) + 
   3.75*C0i(cc11, m12, m22, m32, MH32, MH52, MH52) + 
   1.*C0i(cc12, m12, m22, m32, MH2, MH32, MH32) + 
   3.75*C0i(cc12, m12, m22, m32, MH32, MH52, MH52) - 
   1.*C0i(cc12, m22, m32, m12, MH2, MH32, MH32) + 
   1.25*C0i(cc12, m22, m32, m12, MH32, MH32, MH52) - 
   3.75*C0i(cc12, m22, m32, m12, MH32, MH52, MH52) + 
   1.25*C0i(cc12, m32, m12, m22, MH32, MH32, MH52) - 
   1.*C0i(cc2, m22, m32, m12, MH2, MH32, MH32) + 
   1.25*C0i(cc2, m22, m32, m12, MH32, MH32, MH52) - 
   3.75*C0i(cc2, m22, m32, m12, MH32, MH52, MH52) - 
   1.*C0i(cc22, m22, m32, m12, MH2, MH32, MH32) + 
   1.25*C0i(cc22, m22, m32, m12, MH32, MH32, MH52) - 
   3.75*C0i(cc22, m22, m32, m12, MH32, MH52, MH52)))/
 (PI*(pow(2*1., 1./2))*(pow(sw*1., 2*1.)))
;
}

ComplexType SH30ZA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return 0
;
}

ComplexType SH30ZZ(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return 0
;
}

ComplexType SH3pWA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return (Alfa*M2coeff*(-4.*C0i(cc12, m22, m32, m12, MH2, MH32, MH32) - 
   15.*C0i(cc12, m22, m32, m12, MH32, MH52, MH52) + 
   5.*C0i(cc12, m32, m12, m22, MH32, MH32, MH52) - 
   4.*C0i(cc2, m22, m32, m12, MH2, MH32, MH32) - 
   15.*C0i(cc2, m22, m32, m12, MH32, MH52, MH52) - 
   4.*C0i(cc22, m22, m32, m12, MH2, MH32, MH32) - 
   15.*C0i(cc22, m22, m32, m12, MH32, MH52, MH52)))/
 (PI*sw*(pow(2*1., 1./2)))
;
}

ComplexType SH3pWZ(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return -(Alfa*M2coeff*(4.*(1.*C0i(cc12, m22, m32, m12, MH2, MH32, MH32) - 
      1.25*C0i(cc12, m32, m12, m22, MH32, MH32, MH52) + 
      1.*C0i(cc2, m22, m32, m12, MH2, MH32, MH32) + 
      1.*C0i(cc22, m22, m32, m12, MH2, MH32, MH32))*
     ((pow(cw*1., 2*1.)) - 1.*(pow(sw*1., 2*1.))) + 
    15.*(1.*C0i(cc12, m22, m32, m12, MH32, MH52, MH52) + 
      1.*C0i(cc2, m22, m32, m12, MH32, MH52, MH52) + 
      1.*C0i(cc22, m22, m32, m12, MH32, MH52, MH52))*
     ((pow(cw*1., 2*1.)) - 1.*(pow(sw*1., 2*1.))) - 
    4.*(1.*C0i(cc1, m12, m22, m32, MH2, MH32, MH32) + 
      1.*C0i(cc11, m12, m22, m32, MH2, MH32, MH32) + 
      1.*C0i(cc12, m12, m22, m32, MH2, MH32, MH32))*
     ((pow(cw*1., 2*1.)) + (pow(sw*1., 2*1.))) - 
    15.*(1.*C0i(cc1, m12, m22, m32, MH32, MH52, MH52) + 
      1.*C0i(cc11, m12, m22, m32, MH32, MH52, MH52) + 
      1.*C0i(cc12, m12, m22, m32, MH32, MH52, MH52))*
     ((pow(cw*1., 2*1.)) + (pow(sw*1., 2*1.))) - 
    5.*(1.*C0i(cc12, m22, m32, m12, MH32, MH32, MH52) + 
      1.*C0i(cc2, m22, m32, m12, MH32, MH32, MH52) + 
      1.*C0i(cc22, m22, m32, m12, MH32, MH32, MH52))*
     ((pow(cw*1., 2*1.)) + (pow(sw*1., 2*1.)))))/
 (2*cw*PI*(pow(2*1., 1./2))*(pow(sw*1., 2*1.)))
;
}

ComplexType SH50AA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return (2*Alfa*M2coeff*(1.*C0i(cc12, m22, m32, m12, MH32, MH32, MH32) + 
   7.*C0i(cc12, m22, m32, m12, MH52, MH52, MH52) + 
   1.*C0i(cc2, m22, m32, m12, MH32, MH32, MH32) + 
   7.*C0i(cc2, m22, m32, m12, MH52, MH52, MH52) + 
   1.*C0i(cc22, m22, m32, m12, MH32, MH32, MH32) + 
   7.*C0i(cc22, m22, m32, m12, MH52, MH52, MH52))*(pow(6*1., 1./2)))/PI
;
}

ComplexType SH50WW(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return (0.5*Alfa*M2coeff*(1.*C0i(cc1, m12, m32, m22, MH32, MH32, MH52) + 
   1.*C0i(cc11, m12, m32, m22, MH32, MH32, MH52) + 
   1.*C0i(cc12, m12, m32, m22, MH32, MH32, MH52) - 
   8.*C0i(cc12, m22, m12, m32, MH2, MH32, MH32) - 
   7.*C0i(cc12, m22, m12, m32, MH32, MH52, MH52) + 
   4.*C0i(cc12, m22, m32, m12, MH2, MH32, MH52) - 
   3.*C0i(cc12, m22, m32, m12, MH32, MH32, MH32) - 
   21.*C0i(cc12, m22, m32, m12, MH52, MH52, MH52) + 
   4.*C0i(cc12, m32, m22, m12, MH2, MH32, MH52) + 
   4.*C0i(cc2, m22, m32, m12, MH2, MH32, MH52) - 
   3.*C0i(cc2, m22, m32, m12, MH32, MH32, MH32) - 
   21.*C0i(cc2, m22, m32, m12, MH52, MH52, MH52) + 
   4.*C0i(cc2, m32, m22, m12, MH2, MH32, MH52) + 
   4.*C0i(cc22, m22, m32, m12, MH2, MH32, MH52) - 
   3.*C0i(cc22, m22, m32, m12, MH32, MH32, MH32) - 
   21.*C0i(cc22, m22, m32, m12, MH52, MH52, MH52) + 
   4.*C0i(cc22, m32, m22, m12, MH2, MH32, MH52)))/
 (PI*(pow(6*1., 1./2))*(pow(sw*1., 2*1.)))
;
}

ComplexType SH50ZA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return (Alfa*M2coeff*(1.*C0i(cc12, m22, m32, m12, MH32, MH32, MH32) + 
   7.*C0i(cc12, m22, m32, m12, MH52, MH52, MH52) + 
   1.*C0i(cc2, m22, m32, m12, MH32, MH32, MH32) + 
   7.*C0i(cc2, m22, m32, m12, MH52, MH52, MH52) + 
   1.*C0i(cc22, m22, m32, m12, MH32, MH32, MH32) + 
   7.*C0i(cc22, m22, m32, m12, MH52, MH52, MH52))*(pow(6*1., 1./2))*
  ((pow(cw*1., 2*1.)) - (pow(sw*1., 2*1.))))/(cw*PI*sw)
;
}

ComplexType SH50ZZ(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return (Alfa*M2coeff*(3.*(1.*C0i(cc12, m22, m32, m12, MH32, MH32, MH32) + 
     1.*C0i(cc2, m22, m32, m12, MH32, MH32, MH32) + 
     1.*C0i(cc22, m22, m32, m12, MH32, MH32, MH32))*
    (pow(((pow(cw*1., 2*1.)) - (pow(sw*1., 2*1.)))*1., 2*1.)) + 
   21.*(1.*C0i(cc12, m22, m32, m12, MH52, MH52, MH52) + 
     1.*C0i(cc2, m22, m32, m12, MH52, MH52, MH52) + 
     1.*C0i(cc22, m22, m32, m12, MH52, MH52, MH52))*
    (pow(((pow(cw*1., 2*1.)) - (pow(sw*1., 2*1.)))*1., 2*1.)) - 
   1.*(1.*C0i(cc1, m12, m32, m22, MH32, MH32, MH52) + 
     1.*C0i(cc11, m12, m32, m22, MH32, MH32, MH52) + 
     1.*C0i(cc12, m12, m32, m22, MH32, MH32, MH52) - 
     8.*C0i(cc12, m22, m12, m32, MH2, MH32, MH32))*
    (pow(((pow(cw*1., 2*1.)) + (pow(sw*1., 2*1.)))*1., 2*1.)) + 
   7.*C0i(cc12, m22, m12, m32, MH32, MH52, MH52)*
    (pow(((pow(cw*1., 2*1.)) + (pow(sw*1., 2*1.)))*1., 2*1.)) - 
   4.*(1.*C0i(cc12, m22, m32, m12, MH2, MH32, MH52) + 
     1.*C0i(cc12, m32, m22, m12, MH2, MH32, MH52) + 
     1.*C0i(cc2, m22, m32, m12, MH2, MH32, MH52) + 
     1.*C0i(cc2, m32, m22, m12, MH2, MH32, MH52) + 
     1.*C0i(cc22, m22, m32, m12, MH2, MH32, MH52) + 
     1.*C0i(cc22, m32, m22, m12, MH2, MH32, MH52))*
    (pow(((pow(cw*1., 2*1.)) + (pow(sw*1., 2*1.)))*1., 2*1.))))/
 (PI*(pow(6*1., 1./2))*(pow(cw*1., 2*1.))*(pow(sw*1., 2*1.)))
;
}

ComplexType SH5pWA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return (3*Alfa*M2coeff*(1.*C0i(cc12, m22, m32, m12, MH32, MH32, MH32) + 
   7.*C0i(cc12, m22, m32, m12, MH52, MH52, MH52) + 
   1.*C0i(cc2, m22, m32, m12, MH32, MH32, MH32) + 
   7.*C0i(cc2, m22, m32, m12, MH52, MH52, MH52) + 
   1.*C0i(cc22, m22, m32, m12, MH32, MH32, MH32) + 
   7.*C0i(cc22, m22, m32, m12, MH52, MH52, MH52)))/
 (PI*sw*(pow(2*1., 1./2)))
;
}

ComplexType SH5pWZ(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return (Alfa*M2coeff*(3.*(1.*C0i(cc12, m22, m32, m12, MH32, MH32, MH32) + 
     1.*C0i(cc2, m22, m32, m12, MH32, MH32, MH32) + 
     1.*C0i(cc22, m22, m32, m12, MH32, MH32, MH32))*
    ((pow(cw*1., 2*1.)) - (pow(sw*1., 2*1.))) + 
   21.*(1.*C0i(cc12, m22, m32, m12, MH52, MH52, MH52) + 
     1.*C0i(cc2, m22, m32, m12, MH52, MH52, MH52) + 
     1.*C0i(cc22, m22, m32, m12, MH52, MH52, MH52))*
    ((pow(cw*1., 2*1.)) - 1.*(pow(sw*1., 2*1.))) - 
   1.*(1.*C0i(cc1, m12, m32, m22, MH32, MH32, MH52) + 
     1.*C0i(cc11, m12, m32, m22, MH32, MH32, MH52) + 
     1.*C0i(cc12, m12, m32, m22, MH32, MH32, MH52) - 
     8.*C0i(cc12, m22, m12, m32, MH2, MH32, MH32))*
    ((pow(cw*1., 2*1.)) + (pow(sw*1., 2*1.))) + 
   7.*C0i(cc12, m22, m12, m32, MH32, MH52, MH52)*((pow(cw*1., 2*1.)) + 
     (pow(sw*1., 2*1.))) - 
   4.*(1.*C0i(cc12, m22, m32, m12, MH2, MH32, MH52) + 
     1.*C0i(cc2, m22, m32, m12, MH2, MH32, MH52) + 
     1.*C0i(cc22, m22, m32, m12, MH2, MH32, MH52))*
    ((pow(cw*1., 2*1.)) + (pow(sw*1., 2*1.))) - 
   4.*(1.*C0i(cc12, m32, m22, m12, MH2, MH32, MH52) + 
     1.*C0i(cc2, m32, m22, m12, MH2, MH32, MH52) + 
     1.*C0i(cc22, m32, m22, m12, MH2, MH32, MH52))*
    ((pow(cw*1., 2*1.)) + (pow(sw*1., 2*1.)))))/
 (2*cw*PI*(pow(2*1., 1./2))*(pow(sw*1., 2*1.)))
;
}

ComplexType SH5ppWW(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return (0.5*Alfa*M2coeff*(1.*C0i(cc1, m12, m32, m22, MH32, MH32, MH52) + 
   1.*C0i(cc11, m12, m32, m22, MH32, MH32, MH52) + 
   1.*C0i(cc12, m12, m32, m22, MH32, MH32, MH52) - 
   8.*C0i(cc12, m22, m12, m32, MH2, MH32, MH32) - 
   7.*C0i(cc12, m22, m12, m32, MH32, MH52, MH52) + 
   4.*C0i(cc12, m22, m32, m12, MH2, MH32, MH52) - 
   3.*C0i(cc12, m22, m32, m12, MH32, MH32, MH32) - 
   21.*C0i(cc12, m22, m32, m12, MH52, MH52, MH52) + 
   4.*C0i(cc12, m32, m22, m12, MH2, MH32, MH52) + 
   4.*C0i(cc2, m22, m32, m12, MH2, MH32, MH52) - 
   3.*C0i(cc2, m22, m32, m12, MH32, MH32, MH32) - 
   21.*C0i(cc2, m22, m32, m12, MH52, MH52, MH52) + 
   4.*C0i(cc2, m32, m22, m12, MH2, MH32, MH52) + 
   4.*C0i(cc22, m22, m32, m12, MH2, MH32, MH52) - 
   3.*C0i(cc22, m22, m32, m12, MH32, MH32, MH32) - 
   21.*C0i(cc22, m22, m32, m12, MH52, MH52, MH52) + 
   4.*C0i(cc22, m32, m22, m12, MH2, MH32, MH52)))/(PI*(pow(sw*1., 2*1.)))
;
}

ComplexType SHAA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return (-4.*Alfa*M2coeff*(1.*C0i(cc12, m22, m32, m12, MH32, MH32, MH32) - 
   5.*C0i(cc12, m22, m32, m12, MH52, MH52, MH52) + 
   1.*C0i(cc2, m22, m32, m12, MH32, MH32, MH32) - 
   5.*C0i(cc2, m22, m32, m12, MH52, MH52, MH52) + 
   1.*C0i(cc22, m22, m32, m12, MH32, MH32, MH32) - 
   5.*C0i(cc22, m22, m32, m12, MH52, MH52, MH52))*(pow(3*1., 1./2)))/PI
;
}

ComplexType SHWW(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return (-8.*Alfa*M2coeff*(1.*C0i(cc1, m12, m32, m22, MH2, MH2, MH32) + 
   0.625*C0i(cc1, m12, m32, m22, MH32, MH32, MH52) + 
   1.*C0i(cc11, m12, m32, m22, MH2, MH2, MH32) + 
   0.625*C0i(cc11, m12, m32, m22, MH32, MH32, MH52) + 
   1.*C0i(cc12, m12, m32, m22, MH2, MH2, MH32) + 
   0.625*C0i(cc12, m12, m32, m22, MH32, MH32, MH52) - 
   0.5*C0i(cc12, m22, m12, m32, MH2, MH32, MH32) + 
   0.625*C0i(cc12, m22, m12, m32, MH32, MH52, MH52) + 
   0.375*C0i(cc12, m22, m32, m12, MH32, MH32, MH32) - 
   1.875*C0i(cc12, m22, m32, m12, MH52, MH52, MH52) + 
   0.375*C0i(cc2, m22, m32, m12, MH32, MH32, MH32) - 
   1.875*C0i(cc2, m22, m32, m12, MH52, MH52, MH52) + 
   0.375*C0i(cc22, m22, m32, m12, MH32, MH32, MH32) - 
   1.875*C0i(cc22, m22, m32, m12, MH52, MH52, MH52)))/
 (PI*(pow(3*1., 1./2))*(pow(sw*1., 2*1.)))
;
}

ComplexType SHZA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return (-2.*Alfa*M2coeff*(1.*C0i(cc12, m22, m32, m12, MH32, MH32, MH32) - 
   5.*C0i(cc12, m22, m32, m12, MH52, MH52, MH52) + 
   1.*C0i(cc2, m22, m32, m12, MH32, MH32, MH32) - 
   5.*C0i(cc2, m22, m32, m12, MH52, MH52, MH52) + 
   1.*C0i(cc22, m22, m32, m12, MH32, MH32, MH32) - 
   5.*C0i(cc22, m22, m32, m12, MH52, MH52, MH52))*(pow(3*1., 1./2))*
  ((pow(cw*1., 2*1.)) - 1.*(pow(sw*1., 2*1.))))/(cw*PI*sw)
;
}

ComplexType SHZZ(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32)
{
return (Alfa*M2coeff*(-3.*(1.*C0i(cc12, m22, m32, m12, MH32, MH32, MH32) + 
     1.*C0i(cc2, m22, m32, m12, MH32, MH32, MH32) + 
     1.*C0i(cc22, m22, m32, m12, MH32, MH32, MH32))*
    (pow(((pow(cw*1., 2*1.)) - (pow(sw*1., 2*1.)))*1., 2*1.)) - 
   8.*(1.*C0i(cc1, m12, m32, m22, MH2, MH2, MH32) + 
     1.*C0i(cc11, m12, m32, m22, MH2, MH2, MH32) + 
     1.*C0i(cc12, m12, m32, m22, MH2, MH2, MH32))*
    (pow(((pow(cw*1., 2*1.)) + (pow(sw*1., 2*1.)))*1., 2*1.)) - 
   5.*(1.*C0i(cc1, m12, m32, m22, MH32, MH32, MH52) + 
     1.*C0i(cc11, m12, m32, m22, MH32, MH32, MH52) + 
     1.*C0i(cc12, m12, m32, m22, MH32, MH32, MH52))*
    (pow(((pow(cw*1., 2*1.)) + (pow(sw*1., 2*1.)))*1., 2*1.)) + 
   4.*C0i(cc12, m22, m12, m32, MH2, MH32, MH32)*
    (pow(((pow(cw*1., 2*1.)) + (pow(sw*1., 2*1.)))*1., 2*1.)) + 
   5.*(3.*(1.*C0i(cc12, m22, m32, m12, MH52, MH52, MH52) + 
       1.*C0i(cc2, m22, m32, m12, MH52, MH52, MH52) + 
       1.*C0i(cc22, m22, m32, m12, MH52, MH52, MH52))*
      (pow(((pow(cw*1., 2*1.)) - (pow(sw*1., 2*1.)))*1., 2*1.)) - 
     1.*C0i(cc12, m22, m12, m32, MH32, MH52, MH52)*
      (pow(((pow(cw*1., 2*1.)) + (pow(sw*1., 2*1.)))*1., 
        2*1.)))))/(PI*(pow(3*1., 1./2))*(pow(cw*1., 2*1.))*
  (pow(sw*1., 2*1.)))
;
}

} //end namespace GM
