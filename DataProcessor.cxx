#include <Riostream.h>
#include <string>
#include <vector>

#include "TObject.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TF1.h"
#include "TFile.h"
#include "DataProcessor.h"

ClassImp(DataProcessor)

//______________________________________________________________________________
DataProcessor::DataProcessor(): TObject() {
   // default constructor
}
//______________________________________________________________________________
DataProcessor::DataProcessor(TTree *treeData): TObject() {
  // standard constructor
  fTreeData = (TTree*) treeData -> Clone();
  fTreeData -> SetBranchAddress("FiredTriggerClasses",fTrigClass);
  fTreeData -> SetBranchAddress("NMuons",&fNMuons);
  fTreeData -> SetBranchAddress("Vertex",fVertex);
  fTreeData -> SetBranchAddress("PercentV0M",&fPercV0M);
  fTreeData -> SetBranchAddress("Pt",fPt);
  fTreeData -> SetBranchAddress("E",fE);
  fTreeData -> SetBranchAddress("Px",fPx);
  fTreeData -> SetBranchAddress("Py",fPy);
  fTreeData -> SetBranchAddress("Pz",fPz);
  fTreeData -> SetBranchAddress("Y",fY);
  fTreeData -> SetBranchAddress("Eta",fEta);
  fTreeData -> SetBranchAddress("MatchTrig",fMatchTrig);
  fTreeData -> SetBranchAddress("MatchTrigChi2",fMatchTrigChi2);
  fTreeData -> SetBranchAddress("Charge",fCharge);
  fTreeData -> SetBranchAddress("RAtAbsEnd",fRAtAbsEnd);
  fTreeData -> SetBranchAddress("NDimu",&fNDimu);
  fTreeData -> SetBranchAddress("DimuPt",fDimuPt);
  fTreeData -> SetBranchAddress("DimuPx",fDimuPx);
  fTreeData -> SetBranchAddress("DimuPy",fDimuPy);
  fTreeData -> SetBranchAddress("DimuPz",fDimuPz);
  fTreeData -> SetBranchAddress("DimuY",fDimuY);
  fTreeData -> SetBranchAddress("DimuMass",fDimuMass);
  fTreeData -> SetBranchAddress("DimuCharge",fDimuCharge);
  fTreeData -> SetBranchAddress("DimuMatch",fDimuMatch);
  fTreeData -> SetBranchAddress("DimuMu",fDimuMu);
  fTreeData -> SetBranchAddress("CostHE",fCostHE);
  fTreeData -> SetBranchAddress("PhiHE",fPhiHE);
  fTreeData -> SetBranchAddress("CostCS",fCostCS);
  fTreeData -> SetBranchAddress("PhiCS",fPhiCS);
  fTreeData -> SetBranchAddress("IsPhysSelected",&fIsPhysSelected);
  fTreeData -> SetBranchAddress("pDCA",fPDCA);
  fTreeData -> SetBranchAddress("MuonId",fMuonId);
}
//______________________________________________________________________________
DataProcessor::~DataProcessor() {
  // destructor
}
//______________________________________________________________________________
void DataProcessor::SetPtBins(Int_t nPtBins, Double_t minPtBin[],Double_t maxPtBin[]) {
  fNPtBins = nPtBins;
  for(int i = 0;i < fNPtBins;i++){
    fMinPt.push_back(minPtBin[i]);
    fMaxPt.push_back(maxPtBin[i]);
  }
}
//______________________________________________________________________________
void DataProcessor::SetBinning(vector <Double_t> CostValues, vector <Double_t> PhiValues) {
  for(int i = 0;i < (int) CostValues.size();i++){fCostValues.push_back(CostValues[i]);}
  fNCostBins = fCostValues.size() - 1;
  for(int i = 0;i < (int) PhiValues.size();i++){fPhiValues.push_back(PhiValues[i]);}
  fNPhiBins = fPhiValues.size() - 1;
}
//______________________________________________________________________________
void DataProcessor::ComputeTriggerResponseFunction(string strSample, string nameOutputFile) {
  int nEvents = 0;

  if(strSample == "FullStat"){nEvents = fTreeData -> GetEntries();}
  if(strSample == "TestStat"){nEvents = 100000;}
  printf("N events = %i \n",nEvents);

  fHistLowPtSM = new TH1D("fHistLowPtSM","",100,0,10);
  fHistLowPtSM -> Sumw2();
  fHistAllPtSM = new TH1D("fHistAllPtSM","",100,0,10);
  fHistAllPtSM -> Sumw2();

  fHistLowPtSMpDCA = new TH1D("fHistLowPtSMpDCA","",100,0,10);
  fHistLowPtSMpDCA -> Sumw2();
  fHistAllPtSMpDCA = new TH1D("fHistAllPtSMpDCA","",100,0,10);
  fHistAllPtSMpDCA -> Sumw2();

  vector <int> fListMuonId;

  for(int i = 0;i < nEvents;i++){
    fTreeData -> GetEntry(i);
    printf("Reading : %2.1f %% \r",((double) i)/((double) nEvents)*100);
    TString Trigger = fTrigClass;
    Bool_t TriggerSelected = kFALSE;
    if(Trigger.Contains("CINT7-B-NOPF-MUFAST")){TriggerSelected = kTRUE;} // single muon trigger
    if(TriggerSelected){
      if(fIsPhysSelected){
      for(int k = 0;k < fNDimu;k++){
        if(fDimuY[k] > -4. && fDimuY[k] < -2.5){
            //if(DimuMatch[k] == 2){
              if(fDimuMass[k] > 2 && fDimuMass[k] < 5){
                fListMuonId.push_back(fDimuMu[k][0]);
                fListMuonId.push_back(fDimuMu[k][1]);
                //printf("%i - %i (Dimu Px = %f) \n",DimuMu[k][0],DimuMu[k][1],DimuPx[k]);
              }
            //}
          }
        }

        //printf("=== MUON COMPARISON === \n");
        for(int k = 0;k < fNMuons;k++){
          //printf("%i ",MuonId[k]);

          for(int j = 0;j < (int) fListMuonId.size();j++){
            if(fMuonId[k] == fListMuonId[j]){
              //printf(" -> ok! (Px = %f) \n",Px[k]);

              if(fMatchTrig[k] >= 1){fHistAllPtSM -> Fill(fPt[k]);}
              if(fMatchTrig[k] >= 2){fHistLowPtSM -> Fill(fPt[k]);}

              if(fMatchTrig[k] >= 1 && fPDCA[k] != 0){fHistAllPtSMpDCA -> Fill(fPt[k]);}
              if(fMatchTrig[k] >= 2 && fPDCA[k] != 0){fHistLowPtSMpDCA -> Fill(fPt[k]);}

              break;
            }
          }
        }
        fListMuonId.clear();
      }

    }
  }
  printf("\n");


  TH1D *fHistTriggerResponseFunctionSM = new TH1D("histTriggerResponseFunctionSM","",100,0,10);
  fHistTriggerResponseFunctionSM -> Divide(fHistLowPtSM,fHistAllPtSM,1,1,"B");
  fHistTriggerResponseFunctionSM -> Draw();

  TH1D *fHistTriggerResponseFunctionSMpDCA = new TH1D("histTriggerResponseFunctionSMpDCA","",100,0,10);
  fHistTriggerResponseFunctionSMpDCA -> Divide(fHistLowPtSMpDCA,fHistAllPtSMpDCA,1,1,"B");

}
