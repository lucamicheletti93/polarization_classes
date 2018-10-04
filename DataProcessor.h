#ifndef DATAPROCESSOR_H
#define DATAPROCESSOR_H
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <string>
#include <vector>

class DataProcessor : public TObject
{

 public:
   DataProcessor();
   DataProcessor(TTree *treeData);
   virtual ~DataProcessor();

   void SetPtBins(Int_t, Double_t [],Double_t []);
   void SetBinning(vector <Double_t> , vector <Double_t>);
   void ComputeTriggerResponseFunction(string strSample, string nameOutputFile);

 private:
   Int_t fNPtBins;
   vector <Double_t> fMinPt;
   vector <Double_t> fMaxPt;

   Int_t fNCostBins;
   vector <Double_t> fCostValues;
   Int_t fNPhiBins;
   vector <Double_t> fPhiValues;

   TTree *fTreeData;

   vector <Int_t> fListMuonId;

   TH1D *fHistLowPtSM;
   TH1D *fHistAllPtSM;
   TH1D *fHistTriggerResponseFunctionSM;
   TH1D *fHistLowPtSMpDCA;
   TH1D *fHistAllPtSMpDCA;
   TH1D *fHistTriggerResponseFunctionSMpDCA;

   char fTrigClass[500];
   Float_t fPercV0M, fPercCL0, fPercCL1;
   Int_t fNMuons, fNTracklets, fNContributors;
   Double_t fVertex[3];
   Double_t fPt[300], fE[300], fPx[300], fPy[300], fPz[300], fY[300], fEta[300];
   Double_t fTrackChi2[300], fMatchTrigChi2[300], fDCA[300], fRAtAbsEnd[300];
   Int_t fCharge[300], fMatchTrig[300];
   Int_t fNDimu;
   Double_t fDimuPt[3000], fDimuPx[3000], fDimuPy[3000], fDimuPz[3000], fDimuY[3000];
   Double_t fDimuMass[3000];
   Int_t fDimuCharge[3000], fDimuMatch[3000];
   Int_t fDimuMu[3000][2];
   Double_t fCostHE[3000], fPhiHE[3000], fCostCS[3000], fPhiCS[3000];
   UInt_t fInpmask;
   Bool_t fIsPhysSelected;
   Int_t fPDCA[300];
   Int_t fMuonId[300];

ClassDef(DataProcessor,1)
};

#endif
