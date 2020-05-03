#ifndef SPECIALFITCALCULATOR_H
#define SPECIALFITCALCULATOR_H
#include "TObject.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TGraph.h"
#include <string>
#include <vector>
#include "TMatrixD.h"

class SpecialFitCalculator : public TObject
{

 public:
   SpecialFitCalculator();                                                      // Default constructor
   virtual ~SpecialFitCalculator();                                             // Destructor

   void SetBinning(vector <Double_t> CosThetaValues, vector <Double_t> PhiValues, vector <Double_t> PhiTildeValues);
   void SetFitRange(Int_t [], Int_t []);
   void SetFitNormalization(Double_t initNormCosTheta, Double_t initNormPhi, Double_t initNormPhiTilde);
   void SetFitPolParameters(Double_t initLambdaTheta, Double_t initLambdaPhi, Double_t initLambdaThetaPhi);
   void SimultaneousFit(TObjArray *data, Bool_t saveCanvas, string nameCanvas);
   void BarbatruccoFit(TObjArray *data, Bool_t saveCanvas, string nameCanvas);
   void DecoupledFit(TObjArray *data, Bool_t saveCanvas, string nameCanvas, double minFitrange, double maxFitRange);

   vector <Double_t> GetCosThetaParametersList();
   vector <Double_t> GetErrorCosThetaParametersList();
   vector <Double_t> GetPhiParametersList();
   vector <Double_t> GetErrorPhiParametersList();

   void SetPlotTitle(string);

   Double_t GetLambdaTheta(){return fLambdaTheta;};
   Double_t GetErrorLambdaTheta(){return fErrorLambdaTheta;};
   Double_t GetLambdaPhi(){return fLambdaPhi;};
   Double_t GetErrorLambdaPhi(){return fErrorLambdaPhi;};
   Double_t GetLambdaThetaPhi(){return fLambdaThetaPhi;};
   Double_t GetErrorLambdaThetaPhi(){return fErrorLambdaThetaPhi;};

   Double_t GetLambdaThetaHE(){return fLambdaThetaHE;};
   Double_t GetErrorLambdaThetaHE(){return fErrorLambdaThetaHE;};
   Double_t GetLambdaPhiHE(){return fLambdaPhiHE;};
   Double_t GetErrorLambdaPhiHE(){return fErrorLambdaPhiHE;};
   Double_t GetLambdaThetaPhiHE(){return fLambdaThetaPhiHE;};
   Double_t GetErrorLambdaThetaPhiHE(){return fErrorLambdaThetaPhiHE;};

   Double_t GetLambdaThetaCS(){return fLambdaThetaCS;};
   Double_t GetErrorLambdaThetaCS(){return fErrorLambdaThetaCS;};
   Double_t GetLambdaPhiCS(){return fLambdaPhiCS;};
   Double_t GetErrorLambdaPhiCS(){return fErrorLambdaPhiCS;};
   Double_t GetLambdaThetaPhiCS(){return fLambdaThetaPhiCS;};
   Double_t GetErrorLambdaThetaPhiCS(){return fErrorLambdaThetaPhiCS;};

   Double_t GetChiSquare_NDF(){return fChiSquare_NDF;}

   TGraph*  GetContour_lambdaTheta_lambdaPhi(){return fGraContour_lambdaTheta_lambdaPhi;}

   TH1D*    GetHistFit(int index){return fHistFit[index];}
   TF1*     GetFuncFit(int index){return fFuncFit[index];}

 private:
   Bool_t   fInitNorm;
   Double_t fInitNormCosTheta;
   Double_t fInitNormPhi;
   Double_t fInitNormPhiTilde;

   Bool_t   fInitPolPar;
   Double_t fInitLambdaTheta;
   Double_t fInitLambdaPhi;
   Double_t fInitLambdaThetaPhi;

   Double_t fNormCosTheta;
   Double_t fErrorNormCosTheta;
   Double_t fLambdaTheta;
   Double_t fErrorLambdaTheta;
   Double_t fNormPhi;
   Double_t fErrorNormPhi;
   Double_t fLambdaPhi;
   Double_t fErrorLambdaPhi;
   Double_t fNormPhiTilde;
   Double_t fErrorNormPhiTilde;
   Double_t fLambdaThetaPhi;
   Double_t fErrorLambdaThetaPhi;
   Double_t fErrorLambdaThetaLambdaPhi;

   Double_t fNormCosThetaHE, fErrorNormCosThetaHE;
   Double_t fLambdaThetaHE, fErrorLambdaThetaHE;
   Double_t fNormPhiHE, fErrorNormPhiHE;
   Double_t fLambdaPhiHE, fErrorLambdaPhiHE;
   Double_t fNormPhiTildeHE, fErrorNormPhiTildeHE;
   Double_t fLambdaThetaPhiHE, fErrorLambdaThetaPhiHE;
   Double_t fErrorLambdaThetaLambdaPhiHE;

   Double_t fNormCosThetaCS, fErrorNormCosThetaCS;
   Double_t fLambdaThetaCS, fErrorLambdaThetaCS;
   Double_t fNormPhiCS, fErrorNormPhiCS;
   Double_t fLambdaPhiCS, fErrorLambdaPhiCS;
   Double_t fNormPhiTildeCS, fErrorNormPhiTildeCS;
   Double_t fLambdaThetaPhiCS, fErrorLambdaThetaPhiCS;
   Double_t fErrorLambdaThetaLambdaPhiCS;

   Double_t fLambdaTilde, fErrorLambdaTilde;

   Double_t fChiSquare_NDF;

   vector <Double_t> fCosThetaParametersList;
   vector <Double_t> fErrorCosThetaParametersList;
   vector <Double_t> fPhiParametersList;
   vector <Double_t> fErrorPhiParametersList;
   vector <Double_t> fPhiTildeParametersList;
   vector <Double_t> fErrorPhiTildeParametersList;

   TGraph *fGraContour_lambdaTheta_lambdaPhi;

   TH1D *fHistCosTheta, *fHistPhi, *fHistPhiTilde;
   TF1  *fFuncCosTheta, *fFuncPhi, *fFuncPhiTilde;

   TH1D *fHistFit[3];
   TF1  *fFuncFit[3];

   // Global variables
   Double_t gPi;

ClassDef(SpecialFitCalculator,1)
};

#endif
