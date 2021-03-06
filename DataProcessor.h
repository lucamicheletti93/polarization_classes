#ifndef DATAPROCESSOR_H
#define DATAPROCESSOR_H
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TTree.h"
#include <string>
#include <vector>

class DataProcessor : public TObject
{

 public:
   DataProcessor();                                                             // default constructor
   DataProcessor(TTree *treeData);                                              // constructor for elaborating TTree from grid
   DataProcessor(THnSparse *histNVarHE, THnSparse *histNVarCS);                 // constructor for working with THnSparse
   DataProcessor(TTree *treeDataFilteredHE, TTree *treeDataFilteredCS);         // constructor for working with TTree
   DataProcessor(TTree *treeDataFilteredHE, string flagTRF);                    // constructor for working with Trigger Response Function
   virtual ~DataProcessor();

   void SetPtBins(Int_t, Double_t [],Double_t []);
   void SetBinning(vector <Int_t> CostBinsMin, vector <Int_t> CostBinsMax, vector <Double_t> CostValues, vector <Int_t> PhiBinsMin, vector <Int_t> PhiBinsMax, vector <Double_t> PhiValues, vector <Int_t> PhiTildeBinsMin, vector <Int_t> PhiTildeBinsMax, vector <Double_t> PhiTildeValues);
   void CreateTHnSparse(string strSample, Bool_t pDCAapplied, string nameOutputFile);
   void CutTHnSparse(string nameOutputFile);
   void CreateFilteredTree(string strSample, Bool_t pDCAapplied, string nameOutputFile);
   void CutFilteredTree(string nameOutputFile);
   void CreateInvMassHistograms(TFile *, string strSample);
   void ComputeTriggerResponseFunction(string strSample, string nameOutputFile);
   void ComputeTriggerResponseFunction_new(TTree *treeData, string flagTRF, string strSample, string nameOutputFile);

 private:
   Int_t fNPtBins;
   vector <Double_t> fMinPt;
   vector <Double_t> fMaxPt;
   vector <Int_t> fPtBinsMin;
   vector <Int_t> fPtBinsMax;

   Int_t fNCostBins;
   vector <Int_t> fCostBinsMin;
   vector <Int_t> fCostBinsMax;
   vector <Double_t> fCostValues;
   Int_t fNPhiBins;
   vector <Int_t> fPhiBinsMin;
   vector <Int_t> fPhiBinsMax;
   vector <Double_t> fPhiValues;
   Int_t fNPhiTildeBins;
   vector <Int_t> fPhiTildeBinsMin;
   vector <Int_t> fPhiTildeBinsMax;
   vector <Double_t> fPhiTildeValues;

   THnSparse *fHistNVarHE;
   THnSparse *fHistNVarCS;
   TTree *fTreeDataFilteredHE;
   TTree *fTreeDataFilteredCS;

   TTree *fTreeData;
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
   Double_t fPhiTildeHE[3000], fPhiTildeCS[3000];
   UInt_t fInpmask;
   Bool_t fIsPhysSelected;
   Int_t fPDCA[300];
   Int_t fMuonId[300];

   vector <Int_t> fListMuonId;

   TH1D *fHistLowPt_25eta4;
   TH1D *fHistAllPt_25eta4;
   TH1D *fHistLowPt_25eta3;
   TH1D *fHistAllPt_25eta3;
   TH1D *fHistLowPt_3eta35;
   TH1D *fHistAllPt_3eta35;
   TH1D *fHistLowPt_35eta4;
   TH1D *fHistAllPt_35eta4;

   TH1D *fHistCMUL7Triggers;

   Int_t fCMUL7Triggers;

ClassDef(DataProcessor,1)
};

#endif
