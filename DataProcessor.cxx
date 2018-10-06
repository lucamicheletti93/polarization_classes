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
#include "TH3.h"
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
  //fTreeData -> SetBranchAddress("pDCA",fPDCA);
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
void DataProcessor::CreateFilteredTrees(string strSample, string nameOutputFile) {
  int nEvents = 0;
  int indexPt = 0;
  double PI = TMath::Pi();
  int minPt[12] = {0,1,2,3,4,5,6,7,8,9,10,11};
  int maxPt[12] = {1,2,3,4,5,6,7,8,9,10,11,12};

  if(strSample == "FullStat"){nEvents = fTreeData -> GetEntries();}
  if(strSample == "TestStat"){nEvents = 100000;}
  printf("N events = %i \n",nEvents);

  TFile *fileTreeDataFiltered = new TFile(nameOutputFile.c_str(),"RECREATE");
  double DimuMass, CostHE, PhiHE;
  TTree *treeDataFiltered[12];
  TH3D *histDataFiltered[12];
  for(int i = 0;i < 12;i++){
    treeDataFiltered[i] = new TTree(Form("treeDataFiltered_%ipt%i",minPt[i],maxPt[i]),Form("treeDataFiltered_%ipt%i",minPt[i],maxPt[i]));
    treeDataFiltered[i] -> Branch("DimuMass",&DimuMass,"DimuMass/D");
    treeDataFiltered[i] -> Branch("CostHE",&CostHE,"CostHE/D");
    treeDataFiltered[i] -> Branch("PhiHE",&PhiHE,"PhiHE/D");

    histDataFiltered[i] = new TH3D(Form("histDataFiltered_%ipt%i",minPt[i],maxPt[i]),"",100,2,5,100,-1,1,50,0,PI);
  }

  for(int i = 0;i < nEvents;i++){
    printf("Reading : %2.1f %% \r",((double) i)/((double) nEvents)*100);
    fTreeData -> GetEntry(i);
    for(int k = 0;k < fNDimu;k++){

      if(fIsPhysSelected){
        TString Trigger = fTrigClass;
        Bool_t TriggerSelected = kFALSE;
        if(Trigger.Contains("CMUL7-B-NOPF-MUFAST")) TriggerSelected = kTRUE;
        if(fDimuY[k] > -4. && fDimuY[k] < -2.5){
          if(TriggerSelected){
            if(fDimuMatch[k] == 2){
              if(fDimuMass[k] > 2 && fDimuMass[k] < 5){
                DimuMass = fDimuMass[k];
                CostHE = fCostHE[k];
                PhiHE = TMath::Abs(fPhiHE[k]);
                if(fDimuPt[k] > 0 && fDimuPt[k] <= 12){
                  while(fDimuPt[k] < minPt[indexPt] || fDimuPt[k] > maxPt[indexPt]){indexPt++;}
                  treeDataFiltered[indexPt] -> Fill();
                  histDataFiltered[indexPt] -> Fill(DimuMass,CostHE,PhiHE);
                  indexPt = 0;
                }
              }
            }
          }
        }
      }
    }
  }
  printf("\n");

  fileTreeDataFiltered -> cd();
  for(int i = 0;i < 12;i++){
    treeDataFiltered[i] -> Write();
    histDataFiltered[i] -> Write();
  }
  fileTreeDataFiltered -> Close();
}
//______________________________________________________________________________
void DataProcessor::CreateInvMassHistograms(TFile *fileDataFiltered, string nameOutputFile) {
  int nEvents = 0;
  int indexCost = 0;
  const int nCostBins = fNCostBins;
  const int nPhiBins = fNPhiBins;

  int minPt[12] = {0,1,2,3,4,5,6,7,8,9,10,11};
  int maxPt[12] = {1,2,3,4,5,6,7,8,9,10,11,12};

  double DimuMass, CostHE, PhiHE;
  TTree *treeDataFiltered[12];
  for(int i = 0;i < 12;i++){
    treeDataFiltered[i] = (TTree*) fileDataFiltered -> Get(Form("treeDataFiltered_%ipt%i",minPt[i],maxPt[i]));
    treeDataFiltered[i] -> SetBranchAddress("DimuMass",&DimuMass);
    treeDataFiltered[i] -> SetBranchAddress("CostHE",&CostHE);
    treeDataFiltered[i] -> SetBranchAddress("PhiHE",&PhiHE);
    //treeDataFiltered[i] -> SetDirectory(0);
  }
  //fileDataFiltered -> Close();

  TH1D *histMassCost[nCostBins];

  TFile *fileHistMass = new TFile(nameOutputFile.c_str(),"RECREATE");
  for(int i = 0;i < 12;i++){
    // Inizializing the histogram for each pT bin
    for(int j = 0;j < nCostBins;j++){
      histMassCost[j] = new TH1D(Form("histMassCost_%ipt%i_%4.3fcost%4.3f",minPt[i],maxPt[i],fCostValues[j],fCostValues[j+1]),"",100,2,5);
    }

    nEvents = treeDataFiltered[i] -> GetEntries();
    for(int j = 0;j < nEvents;j++){
      printf("Reading : %2.1f %% \r",((double) j)/((double) nEvents)*100);
      treeDataFiltered[i] -> GetEntry(j);
      while(CostHE < fCostValues[indexCost] || CostHE > fCostValues[indexCost+1]){indexCost++;}
      histMassCost[indexCost] -> Fill(DimuMass);
      //printf("[%3.2f - %3.2f] - %f \n",fCostValues[indexCost],fCostValues[indexCost+1],CostHE);
      indexCost = 0;
    }

    // Saving the histograms
    fileHistMass -> cd();
    for(int j = 0;j < nCostBins;j++){
      histMassCost[j] -> Write();
    }
  }
  printf("\n");
  fileHistMass -> Close();
  fileDataFiltered -> Close();
}
//______________________________________________________________________________
void DataProcessor::ComputeTriggerResponseFunction(string strSample, string nameOutputFile) {
  fTreeData -> SetBranchAddress("pDCA",fPDCA); // enable pDDCA
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
  fHistTriggerResponseFunctionSM -> SetLineColor(kBlue);

  TH1D *fHistTriggerResponseFunctionSMpDCA = new TH1D("histTriggerResponseFunctionSMpDCA","",100,0,10);
  fHistTriggerResponseFunctionSMpDCA -> Divide(fHistLowPtSMpDCA,fHistAllPtSMpDCA,1,1,"B");
  fHistTriggerResponseFunctionSMpDCA -> SetLineColor(kRed);
  fHistTriggerResponseFunctionSMpDCA -> SetMarkerStyle(20);
  fHistTriggerResponseFunctionSMpDCA -> SetMarkerColor(kRed);

  TFile *fileTriggerResponseFunction = new TFile(nameOutputFile.c_str(),"RECREATE");
  fHistAllPtSM -> Write();
  fHistLowPtSM -> Write();
  fHistAllPtSMpDCA -> Write();
  fHistLowPtSMpDCA -> Write();
  fHistTriggerResponseFunctionSM -> Write();
  fHistTriggerResponseFunctionSMpDCA -> Write();
  fileTriggerResponseFunction -> Close();
}
