#ifndef MASSFITTER_H
#define MASSFITTER_H
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <string>
#include <vector>

class MassFitter : public TObject
{

 public:
   MassFitter();
   virtual ~MassFitter();

   void SetScalingFactor(Double_t);
   void SetBinning(vector <Double_t> , vector <Double_t>, vector <Double_t>);
   void fit_of_minv(TH1D *, string, string, int, int, string, string, TH2D *, bool );

   Double_t GetNJpsi(){return fNJpsi;}
   Double_t GetStatJpsi(){return fStatJpsi;}
   Double_t GetMassJpsi(){return fMassJpsi;}
   Double_t GetErrMassJpsi(){return fErrMassJpsi;}
   Double_t GetSigmaJpsi(){return fSigmaJpsi;}
   Double_t GetErrSigmaJpsi(){return fErrSigmaJpsi;}
   Double_t GetChiSquare_NDF(){return fChiSquare_NDF;}
   string GetFitStatus(){return fFitStatus;}

 private:
   Double_t fPi;

   Double_t fScalingFactorJpsiSigma;

   string fFitStatus;

   Double_t fNJpsi, fStatJpsi;
   Double_t fMassJpsi, fErrMassJpsi;
   Double_t fSigmaJpsi, fErrSigmaJpsi;
   Double_t fChiSquare_NDF;

   Int_t fNCosThetaBins;
   vector <Double_t> fCosThetaValues;
   Int_t fNPhiBins;
   vector <Double_t> fPhiValues;
   Int_t fNPhiTildeBins;
   vector <Double_t> fPhiTildeValues;

ClassDef(MassFitter,1)
};

#endif
