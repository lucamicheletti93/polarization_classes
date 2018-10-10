// MathFuncsLib.h
// To include in a macro this library :
//          - .L MathFuncsLib.cxx+ (in python use gROOT.ProcessLineSync(".L ../MathFuncsLib.cxx+"))
//                      -> in c++ :
//                                + #include "../MathFuncsLib.h"
//                                + TF1 *func = new TF1("func",MathFuncs::MyMathFuncs::MyFuncPolPhi,0,10,nPar)
//                      -> in python :
//                                + func = TF1("func",MathFuncs.MyMathFuncs.MyFuncPolPhi,0,10,nPar) ()

namespace MathFuncs
{
    class MyMathFuncs
    {
    public:
        static double MyFuncPolCosTheta(double *x, double *par);
        static double MyFuncPolPhi(double *x, double *par);
    };
}
