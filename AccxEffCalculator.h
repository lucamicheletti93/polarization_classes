#ifndef ACCXEFFCALCULATOR_H
#define ACCXEFFCALCULATOR_H
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <string>
#include <vector>

class AccxEffCalculator : public TObject
{

 public:
   AccxEffCalculator();
   AccxEffCalculator(TTree *treeAccxEff);
   virtual ~AccxEffCalculator();

   void SetPtBins(Int_t, Double_t [],Double_t []);
   void SetBinning(vector <Double_t> , vector <Double_t>, vector <Double_t>);
   void ComputeAccxEff(string strSample, string nameOutputFile);
   void ReWeightAccxEff(string refFrame, Double_t LambdaTheta, Double_t LambdaPhi, string strSample, Bool_t saveFile, string nameOutputFile);
   void ComputeTriggerResponseFunction(string strSample, string nameOutputFile);

   Int_t GetNPtBins(){return fNPtBins;};
   Int_t GetNCosThetaBins(){return fNCosThetaBins;};
   Int_t GetNPhiBins(){return fNPhiBins;};

 private:
   Double_t fPi;
   Int_t fNPtBins;
   vector <Double_t> fMinPt;
   vector <Double_t> fMaxPt;

   Int_t fNCosThetaBins;
   vector <Double_t> fCosThetaValues;
   Int_t fNPhiBins;
   vector <Double_t> fPhiValues;
   Int_t fNPhiTildeBins;
   vector <Double_t> fPhiTildeValues;

   double fPhiTilde;

   TTree *fTreeAccxEff;

   // Acc x Eff
   // cos(theta)
   TH1D *fHistGenCosThetaHE[13];
   TH1D *fHistRecCosThetaHE[13];
   TH1D *fHistAccxEffCosThetaHE[13];
   TH1D *fHistGenCosThetaCS[13];
   TH1D *fHistRecCosThetaCS[13];
   TH1D *fHistAccxEffCosThetaCS[13];
   // phi
   TH1D *fHistGenPhiHE[13];
   TH1D *fHistRecPhiHE[13];
   TH1D *fHistAccxEffPhiHE[13];
   TH1D *fHistGenPhiCS[13];
   TH1D *fHistRecPhiCS[13];
   TH1D *fHistAccxEffPhiCS[13];
   // phi Tilde
   TH1D *fHistGenPhiTildeHE[13];
   TH1D *fHistRecPhiTildeHE[13];
   TH1D *fHistAccxEffPhiTildeHE[13];
   TH1D *fHistGenPhiTildeCS[13];
   TH1D *fHistRecPhiTildeCS[13];
   TH1D *fHistAccxEffPhiTildeCS[13];
   // cos(theta),phi
   TH2D *fHistGenCosThetaPhiHE[13];
   TH2D *fHistRecCosThetaPhiHE[13];
   TH2D *fHistAccxEffCosThetaPhiHE[13];
   TH2D *fHistAccxEffCosThetaPhiStatRelHE[13];
   TH2D *fHistGenCosThetaPhiCS[13];
   TH2D *fHistRecCosThetaPhiCS[13];
   TH2D *fHistAccxEffCosThetaPhiCS[13];
   TH2D *fHistAccxEffCosThetaPhiStatRelCS[13];
   // cos(theta),pT
   TH2D *fHistGenCosThetaHEPt;
   TH2D *fHistRecCosThetaHEPt;
   TH2D *fHistAccxEffCosThetaHEPt;
   TH2D *fHistGenCosThetaCSPt;
   TH2D *fHistRecCosThetaCSPt;
   TH2D *fHistAccxEffCosThetaCSPt;
   // phi,pT
   TH2D *fHistGenPhiHEPt;
   TH2D *fHistRecPhiHEPt;
   TH2D *fHistAccxEffPhiHEPt;
   TH2D *fHistGenPhiCSPt;
   TH2D *fHistRecPhiCSPt;
   TH2D *fHistAccxEffPhiCSPt;

   // cos(theta)
   TH1D *fHistGenCosThetaReWeighted[13];
   TH1D *fHistRecCosThetaReWeighted[13];
   TH1D *fHistAccxEffCosThetaReWeighted[13];
   // phi
   TH1D *fHistGenPhiReWeighted[13];
   TH1D *fHistRecPhiReWeighted[13];
   TH1D *fHistAccxEffPhiReWeighted[13];

   // Trigger response fuction calculation
   TH1D *fHistLowPt;
   TH1D *fHistAllPt;
   TH1D *fHistTriggerResponseFunction;

   // tree variables
   Int_t fNDimuGen;
   Double_t fDimuPtGen[3000];
   Double_t fDimuYGen[3000];
   Double_t fCosThetaHEGen[3000];
   Double_t fPhiHEGen[3000];
   Double_t fCosThetaCSGen[3000];
   Double_t fPhiCSGen[3000];

   Int_t fNDimuRec;
   Double_t fDimuPtRec[3000];
   Double_t fDimuYRec[3000];
   Double_t fDimuMassRec[3000];
   Int_t fDimuMatchRec[3000];
   Double_t fCosThetaHERec[3000];
   Double_t fPhiHERec[3000];
   Double_t fCosThetaCSRec[3000];
   Double_t fPhiCSRec[3000];

   Int_t fNMuonsRec;
   Double_t fPtRec[3000];
   Int_t fMatchTrigRec[3000];

ClassDef(AccxEffCalculator,1)
};

#endif
