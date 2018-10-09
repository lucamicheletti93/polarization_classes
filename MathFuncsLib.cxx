// MathFuncsLib.cxx
// compile with: cl /c /EHsc MathFuncsLib.cxx
// post-build command: lib MathFuncsLib.obj

#include "MathFuncsLib.h"
#include "TMath.h"

#include <stdexcept>

using namespace std;

namespace MathFuncs
{
    double MyMathFuncs::MyFuncPolCosTheta(double *x, double *par){
      double cosTheta = x[0];
      double N = par[0];
      double lambdaTheta = par[1];
      return (N/(3 + lambdaTheta))*(1 + lambdaTheta*cosTheta*cosTheta);
    }

    double MyMathFuncs::MyFuncPolPhi(double *x, double *par){
      double phi = x[0];
      double N = par[0];
      double lambdaTheta = par[1];
      double lambdaPhi = par[2];
      return N*(1 + ((2*lambdaPhi)/(3 + lambdaTheta))*TMath::Cos(2*phi));
    }
}
