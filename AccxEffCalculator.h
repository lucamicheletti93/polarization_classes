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
   void SetBinning(vector <Double_t> , vector <Double_t>);
   void ComputeAccxEff(string strSample, string nameOutputFile);
   void ReWeightAccxEff(Double_t LambdaTheta, Double_t LambdaPhi, string strSample, string nameOutputFile);
   void ComputeTriggerResponseFunction(string strSample, string nameOutputFile);

   Int_t GetNPtBins(){return fNPtBins;};
   Int_t GetNCostBins(){return fNCostBins;};
   Int_t GetNPhiBins(){return fNPhiBins;};

 private:
   Double_t fPi;
   Int_t fNPtBins;
   vector <Double_t> fMinPt;
   vector <Double_t> fMaxPt;

   Int_t fNCostBins;
   vector <Double_t> fCostValues;
   Int_t fNPhiBins;
   vector <Double_t> fPhiValues;

   TTree *fTreeAccxEff;

   // Acc x Eff
   // cos(theta)
   TH1D *fHistGenCost[13];
   TH1D *fHistRecCost[13];
   TH1D *fHistAccxEffCost[13];
   // phi
   TH1D *fHistGenPhi[13];
   TH1D *fHistRecPhi[13];
   TH1D *fHistAccxEffPhi[13];
   // phi Tilde
   TH1D *fHistGenPhiTilde[13];
   TH1D *fHistRecPhiTilde[13];
   TH1D *fHistAccxEffPhiTilde[13];
   // cos(theta),phi
   TH2D *fHistGenCostPhi[13];
   TH2D *fHistRecCostPhi[13];
   TH2D *fHistAccxEffCostPhi[13];
   TH2D *fHistAccxEffCostPhiStatRel[13];
   // cos(theta),pT
   TH2D *fHistGenCostPt;
   TH2D *fHistRecCostPt;
   TH2D *fHistAccxEffCostPt;
   // phi,pT
   TH2D *fHistGenPhiPt;
   TH2D *fHistRecPhiPt;
   TH2D *fHistAccxEffPhiPt;

   // cos(theta)
   TH1D *fHistGenCostReWeighted[13];
   TH1D *fHistRecCostReWeighted[13];
   TH1D *fHistAccxEffCostReWeighted[13];
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
   Double_t fCostHEGen[3000];
   Double_t fPhiHEGen[3000];
   Double_t fCostCSGen[3000];
   Double_t fPhiCSGen[3000];

   Int_t fNDimuRec;
   Double_t fDimuPtRec[3000];
   Double_t fDimuYRec[3000];
   Double_t fDimuMassRec[3000];
   Int_t fDimuMatchRec[3000];
   Double_t fCostHERec[3000];
   Double_t fPhiHERec[3000];
   Double_t fCostCSRec[3000];
   Double_t fPhiCSRec[3000];

   Int_t fNMuonsRec;
   Double_t fPtRec[3000];
   Int_t fMatchTrigRec[3000];

ClassDef(AccxEffCalculator,1)
};

#endif
