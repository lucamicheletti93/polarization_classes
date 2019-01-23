#ifndef SPECIALFITCALCULATOR_H
#define SPECIALFITCALCULATOR_H
#include "TObject.h"
#include "TObjArray.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <string>
#include <vector>
#include "TMatrixD.h"

class SpecialFitCalculator : public TObject
{

 public:
   SpecialFitCalculator();                                                      // Default constructor
   virtual ~SpecialFitCalculator();                                             // Destructor

   void SetBinning(vector <Double_t> CosThetaValues, vector <Double_t> PhiValues, vector <Double_t> PhiTildeValues);
   void SimultaneousFit(TObjArray *data, Bool_t saveCanvas, string nameCanvas);

   vector <Double_t> GetCosThetaParametersList();
   vector <Double_t> GetErrorCosThetaParametersList();
   vector <Double_t> GetPhiParametersList();
   vector <Double_t> GetErrorPhiParametersList();

   Double_t GetLambdaTheta(){return fLambdaTheta;};
   Double_t GetErrorLambdaTheta(){return fErrorLambdaTheta;};
   Double_t GetLambdaPhi(){return fLambdaPhi;};
   Double_t GetErrorLambdaPhi(){return fErrorLambdaPhi;};
   Double_t GetLambdaThetaPhi(){return fLambdaThetaPhi;};
   Double_t GetErrorLambdaThetaPhi(){return fErrorLambdaThetaPhi;};

 private:

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

   vector <Double_t> fCosThetaParametersList;
   vector <Double_t> fErrorCosThetaParametersList;
   vector <Double_t> fPhiParametersList;
   vector <Double_t> fErrorPhiParametersList;
   vector <Double_t> fPhiTildeParametersList;
   vector <Double_t> fErrorPhiTildeParametersList;

   // Global variables
   Double_t gPi;

ClassDef(SpecialFitCalculator,1)
};

#endif
