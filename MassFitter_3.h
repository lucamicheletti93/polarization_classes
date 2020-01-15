#ifndef MASSFITTER_3_H
#define MASSFITTER_3_H
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include <string>
#include <vector>

class MassFitter_3 : public TObject
{

public:
  MassFitter_3();
  MassFitter_3(TH1D *);
  virtual ~MassFitter_3();

  void SetScalingFactor(Double_t);
  void SetBinning(vector <Double_t> , vector <Double_t>, vector <Double_t>);
  void SetFitRange(Double_t , Double_t);
  void SetSpecialFitConditions(Int_t, Bool_t);
  void SetJpsiWidth(Double_t);
  void SetJpsiMass(Double_t);
  void SetPhilippeCorrectionFactor(Double_t);
  void SetCosThetaPhiIndex(Int_t, Int_t);
  void SetTailsParameters(Double_t [], Double_t []);
  void fit_of_minv(string, string, string, string, Bool_t, Bool_t);

  Double_t GetNJpsi(){return fNJpsi;}
  Double_t GetStatJpsi(){return fStatJpsi;}
  Double_t GetMassJpsi(){return fMassJpsi;}
  Double_t GetErrMassJpsi(){return fErrMassJpsi;}
  Double_t GetSigmaJpsi(){return fSigmaJpsi;}
  Double_t GetErrSigmaJpsi(){return fErrSigmaJpsi;}
  Double_t GetChiSquare_NDF(){return fChiSquare_NDF;}
  string   GetFitStatus(){return fFitStatus;}

  TF1*     GetFuncTot(){return fFuncTot;}

private:
  void PlotStandard(Bool_t);
  void PlotPublications(Bool_t);

  TH1D*                  fHistMinv;
  Bool_t                 fSpecialFitConditions;
  Bool_t                 fJpsiWidthFixed;
  Bool_t                 fJpsiMassFixedToPDG;
  Bool_t                 fJpsiMassFixedToIntegrated;
  Bool_t                 fTailParametersFixed;

  Int_t                  fIndexCosTheta;
  Int_t                  fIndexPhi;

  Double_t               fPi;

  Double_t               fScalingFactorJpsiSigma;
  Double_t               fMinFitRange;
  Double_t               fMaxFitRange;

  string                 fFitStatus;

  Double_t               fMassIntegrated;

  Double_t               fNJpsi, fStatJpsi;
  Double_t               fMassJpsi, fErrMassJpsi;
  Double_t               fSigmaJpsi, fErrSigmaJpsi;
  Double_t               fNPsi2S, fStatPsi2S;
  Double_t               fMassPsi2S, fErrMassPsi2S;
  Double_t               fSigmaPsi2S, fErrSigmaPsi2S;
  Double_t               fChiSquare_NDF;

  Int_t                  fNCosThetaBins;
  vector <Double_t>      fCosThetaValues;
  Int_t fNPhiBins;
  vector <Double_t>      fPhiValues;
  Int_t                  fNPhiTildeBins;
  vector <Double_t>      fPhiTildeValues;

  TF1*                   fFuncBkgFix;
  TF1*                   fFuncSigJpsiFix;
  TF1*                   fFuncSigPsi2SFix;
  TF1*                   fFuncTot;

  string                 fSigShape;
  string                 fBkgShape;
  string                 fOutputDir;

  string                 fPlotType;
  Bool_t                 fSavePlot;
  Int_t                  fRebin;

  Double_t               fPhilippeCorrectionFactor;

  Double_t               fParTailsCB2[4];
  Double_t               fParTailsNA60[8];

ClassDef(MassFitter_3,1)
};

#endif
