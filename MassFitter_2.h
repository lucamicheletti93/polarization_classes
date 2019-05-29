#ifndef MASSFITTER_2_2_H
#define MASSFITTER_2_2_H
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <string>
#include <vector>

class MassFitter_2 : public TObject
{

public:
  MassFitter_2();
  MassFitter_2(TH1D *);
  virtual ~MassFitter_2();

  void SetScalingFactor(Double_t);
  void SetBinning(vector <Double_t> , vector <Double_t>, vector <Double_t>);
  void SetFitRange(Double_t , Double_t);
  void SetSpecialFitConditions();
  void SetJpsiWidth(Double_t);
  void SetCosThetaPhiIndex(Int_t, Int_t);
  void fit_of_minv(string, string, string);

  Double_t GetNJpsi(){return fNJpsi;}
  Double_t GetStatJpsi(){return fStatJpsi;}
  Double_t GetMassJpsi(){return fMassJpsi;}
  Double_t GetErrMassJpsi(){return fErrMassJpsi;}
  Double_t GetSigmaJpsi(){return fSigmaJpsi;}
  Double_t GetErrSigmaJpsi(){return fErrSigmaJpsi;}
  Double_t GetChiSquare_NDF(){return fChiSquare_NDF;}
  string GetFitStatus(){return fFitStatus;}

private:
  void PlotStandard(Bool_t);
  void PlotPublications(Bool_t);

  TH1D*                  fHistMinv;
  TH1D*                  fHistJpsiWidth;
  Bool_t                 fSpecialFitConditions;
  Bool_t                 fJpsiWidthFixed;
  Bool_t                 fTailParametersFixed;

  Int_t                  fIndexCosTheta;
  Int_t                  fIndexPhi;
   
  Double_t               fPi;

  Double_t               fScalingFactorJpsiSigma;
  Double_t               fMinFitRange;
  Double_t               fMaxFitRange;

  string                 fFitStatus;

  Double_t               fNJpsi, fStatJpsi;
  Double_t               fMassJpsi, fErrMassJpsi;
  Double_t               fSigmaJpsi, fErrSigmaJpsi;
  Double_t               fChiSquare_NDF;

  Int_t                  fNCosThetaBins;
  vector <Double_t>      fCosThetaValues;
  Int_t fNPhiBins;
  vector <Double_t>      fPhiValues;
  Int_t                  fNPhiTildeBins;
  vector <Double_t>      fPhiTildeValues;

  TF1*                   fFuncBkgFix;
  TF1*                   fFuncSigJpsiFix;
  TF1*                   fFuncTot;

  string                 fSigShape;
  string                 fBkgShape;
  string                 fOutputDir;

  string                 fPlotType;
  Bool_t                 fSavePlot;

ClassDef(MassFitter_2,1)
};

#endif
