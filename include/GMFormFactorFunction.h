#ifndef FormFactorFunctionGM
#define FormFactorFunctionGM
#include "GMModelParameters.h"
#include "clooptools.h"

namespace GM{
ComplexType NOIMPLEMENTED(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SH30AA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SH30WW(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SH30ZA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SH30ZZ(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SH3pWA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SH3pWZ(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SH50AA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SH50WW(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SH50ZA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SH50ZZ(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SH5pWA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SH5pWZ(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SH5ppWW(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SHAA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SHWW(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SHZA(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

ComplexType SHZZ(double MH2, double MH32, double MH52, double M2coeff, double m12, double m22, double m32);

} //end namespace GM
typedef ComplexType (*SFunc)(double,double,double,double,double,double,double);
SFunc GetSFunc(PID Mother, PID P1, PID P2);
double GetSFuncEFTRatio(PID Mother, PID P1, PID P2);
double GetHHVCoupling(PID Mother, PID P1, PID P2);
#endif
