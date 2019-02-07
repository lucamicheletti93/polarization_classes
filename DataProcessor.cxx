#include <Riostream.h>
#include <string>
#include <vector>

#include "TObject.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TString.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TFile.h"
#include "THnSparse.h"
#include "DataProcessor.h"

ClassImp(DataProcessor)

//______________________________________________________________________________
DataProcessor::DataProcessor(): TObject() {
   // default constructor
}
//______________________________________________________________________________
DataProcessor::DataProcessor(THnSparse *histNVarHE, THnSparse *histNVarCS): TObject() {
  fHistNVarHE = (THnSparse*) histNVarHE -> Clone();
  fHistNVarCS = (THnSparse*) histNVarCS -> Clone();
}
//______________________________________________________________________________
DataProcessor::DataProcessor(TTree *treeDataFilteredHE, TTree *treeDataFilteredCS): TObject() {
  fTreeDataFilteredHE = (TTree*) treeDataFilteredHE -> Clone();
  fTreeDataFilteredCS = (TTree*) treeDataFilteredCS -> Clone();
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
void DataProcessor::SetBinning(vector <Int_t> CostBinsMin, vector <Int_t> CostBinsMax, vector <Double_t> CostValues, vector <Int_t> PhiBinsMin, vector <Int_t> PhiBinsMax, vector <Double_t> PhiValues, vector <Int_t> PhiTildeBinsMin, vector <Int_t> PhiTildeBinsMax, vector <Double_t> PhiTildeValues) {
  for(int i = 0;i < (int) CostValues.size();i++){fCostValues.push_back(CostValues[i]);}
  fNCostBins = (int) fCostValues.size() - 1;
  for(int i = 0;i < fNCostBins;i++){
    fCostBinsMin.push_back(CostBinsMin[i]);
    fCostBinsMax.push_back(CostBinsMax[i]);
  }
  for(int i = 0;i < (int) PhiValues.size();i++){fPhiValues.push_back(PhiValues[i]);}
  fNPhiBins = (int) fPhiValues.size() - 1;
  for(int i = 0;i < fNPhiBins;i++){
    fPhiBinsMin.push_back(PhiBinsMin[i]);
    fPhiBinsMax.push_back(PhiBinsMax[i]);
  }
  for(int i = 0;i < (int) PhiTildeValues.size();i++){fPhiTildeValues.push_back(PhiTildeValues[i]);}
  fNPhiTildeBins = (int) fPhiTildeValues.size() - 1;
  for(int i = 0;i < fNPhiTildeBins;i++){
    fPhiTildeBinsMin.push_back(PhiTildeBinsMin[i]);
    fPhiTildeBinsMax.push_back(PhiTildeBinsMax[i]);
  }
}
//______________________________________________________________________________
void DataProcessor::CreateTHnSparse(string strSample, Bool_t pDCAapplied, string nameOutputFile) {
  if(pDCAapplied){
    printf("- pDCA included \n");
    fTreeData -> SetBranchAddress("pDCA",fPDCA); // enable pDDCAs
  }
  else{printf("- pDCA not included \n");}
  printf("- Configuring Fiducial Box \n");
  printf("%f < CosTheta < %f \n",fCostValues[1],fCostValues[fNCostBins-1]);
  printf("%f < |Phi| < %f \n",fPhiValues[1],fPhiValues[fNPhiBins-1]);

  double PI = TMath::Pi();
  int nEvents = 0;
  if(strSample == "FullStat"){nEvents = fTreeData -> GetEntries();}
  if(strSample == "TestStat"){nEvents = 100000;}
  printf("N events = %i \n",nEvents);

  TFile *fileTHnSparse = new TFile(nameOutputFile.c_str(),"RECREATE");

  // Defining the THnSparse
  // varArray[4] = {DimuPt, DimuMass, DimuCost, DimuPhi, DimuPhiTilde}
  const int nVar = 5;
  int nBins[nVar] = {100,120,100,50,50};
  double minVar[nVar] = {0.,2.,-1.,0.,0.};
  double maxVar[nVar] = {10.,5.,1.,PI,2*PI};
  double varArray[nVar];
  double tmpVar = 0;    // variable used to store the value of Phi for PhiTilde calculation

  THnSparseD *histNVarHE = new THnSparseD("histNVarHE","histNVarHE",nVar,nBins,minVar,maxVar);
  THnSparseD *histNVarCS = new THnSparseD("histNVarCS","histNVarCS",nVar,nBins,minVar,maxVar);

  bool goodMuon1 = kFALSE;
  bool goodMuon2 = kFALSE;

  for(int i = 0;i < nEvents;i++){
    fTreeData -> GetEntry(i);
    printf("Reading : %2.1f %% \r",((double) i)/((double) nEvents)*100);

    for(int j = 0;j < fNDimu;j++){

      if(fIsPhysSelected){
        TString Trigger = fTrigClass;
        Bool_t TriggerSelected = kFALSE;
        if(Trigger.Contains("CMUL7-B-NOPF-MUFAST")) TriggerSelected = kTRUE;
        if(fDimuY[j] > -4. && fDimuY[j] < -2.5){
          if(TriggerSelected){
            if(fDimuMatch[j] == 2){
              if(fDimuMass[j] > 2 && fDimuMass[j] < 5){
                //printf("%i - %i (Dimu Px = %f) \n",fDimuMu[j][0],fDimuMu[j][1],fDimuPx[j]);
                for(int k = 0;k < fNMuons;k++){
                  if(fMuonId[k] == fDimuMu[j][0]){
                    if(fPDCA[k] != 0){goodMuon1 = kTRUE;}
                    //printf(" -> ok! (Px = %f) [%i] \n",Px[k],k);
                  }
                  if(fMuonId[k] == fDimuMu[j][1]){
                    if(fPDCA[k] != 0){goodMuon2 = kTRUE;}
                    //printf(" -> ok! (Px = %f) [%i] \n",Px[k],k);
                  }
                }
                if(pDCAapplied == kFALSE){
                  goodMuon1 = kTRUE;
                  goodMuon2 = kTRUE;
                }
                if(goodMuon1 == kTRUE && goodMuon2 == kTRUE){
                  if(TMath::Abs(fPhiHE[j]) > fPhiValues[1] && TMath::Abs(fPhiHE[j]) < fPhiValues[fNPhiBins-1]){
                    if(fCostHE[j] > fCostValues[1] && fCostHE[j] < fCostValues[fNCostBins-1]){
                      varArray[0] = fDimuPt[j];
                      varArray[1] = fDimuMass[j];
                      varArray[2] = fCostHE[j];
                      varArray[3] = TMath::Abs(fPhiHE[j]);
                      tmpVar = fPhiHE[j] + PI;
                      //if(fCostHE[j] < 0.){fPhiTildeHE[j] = fPhiHE[j] - (3./4.)*PI;}
                      //if(fCostHE[j] > 0.){fPhiTildeHE[j] = fPhiHE[j] - (1./4.)*PI;}
                      if(fCostHE[j] < 0.){fPhiTildeHE[j] = tmpVar - (3./4.)*PI;}
                      if(fCostHE[j] > 0.){fPhiTildeHE[j] = tmpVar - (1./4.)*PI;}
                      //if(fPhiTildeHE[j] < 0){fPhiTildeHE[j] = 2*PI + fPhiTildeHE[j];}
                      //if(fPhiTildeHE[j] > PI){fPhiTildeHE[j] = 2*PI - fPhiTildeHE[j];}
                      if(fPhiTildeHE[j] > 2*PI){fPhiTildeHE[j] = fPhiTildeHE[j] - 2*PI;}
                      if(fPhiTildeHE[j] < 0.){fPhiTildeHE[j] = 2*PI + fPhiTildeHE[j];}
                      varArray[4] = fPhiTildeHE[j];

                      histNVarHE -> Fill(varArray);                                 // Filling the THnSparse
                    }
                  }

                  if(TMath::Abs(fPhiCS[j]) > fPhiValues[1] && TMath::Abs(fPhiCS[j]) < fPhiValues[fNPhiBins-1]){
                    if(fCostCS[j] > fCostValues[1] && fCostCS[j] < fCostValues[fNCostBins-1]){
                      varArray[0] = fDimuPt[j];
                      varArray[1] = fDimuMass[j];
                      varArray[2] = fCostCS[j];
                      varArray[3] = TMath::Abs(fPhiCS[j]);
                      tmpVar = fPhiCS[j] + PI;
                      //if(fCostCS[j] < 0.){fPhiTildeCS[j] = fPhiCS[j] - (3./4.)*PI;}
                      //if(fCostCS[j] > 0.){fPhiTildeCS[j] = fPhiCS[j] - (1./4.)*PI;}
                      if(fCostCS[j] < 0.){fPhiTildeCS[j] = tmpVar - (3./4.)*PI;}
                      if(fCostCS[j] > 0.){fPhiTildeCS[j] = tmpVar - (1./4.)*PI;}
                      //if(fPhiTildeCS[j] < 0){fPhiTildeCS[j] = 2*PI + fPhiTildeCS[j];}
                      //if(fPhiTildeCS[j] > PI){fPhiTildeCS[j] = 2*PI - fPhiTildeCS[j];}
                      if(fPhiTildeCS[j] > 2*PI){fPhiTildeCS[j] = fPhiTildeCS[j] - 2*PI;}
                      if(fPhiTildeCS[j] < 0.){fPhiTildeCS[j] = 2*PI + fPhiTildeCS[j];}
                      varArray[4] = fPhiTildeCS[j];

                      histNVarCS -> Fill(varArray);                                 // Filling the THnSparse
                      //printf("Both Good Muons! \n");
                      //histMassPDCA -> Fill(DimuMass[j]);
                    }
                  }
                }
                //histMass -> Fill(DimuMass[j]);
                goodMuon1 = kFALSE;
                goodMuon2 = kFALSE;
              }
            }
          }
        }
      }
    }
  }
  printf("\n");

  fileTHnSparse -> cd();
  histNVarHE -> Write();
  histNVarCS -> Write();
  fileTHnSparse -> Close();
}
//______________________________________________________________________________
void DataProcessor::CutTHnSparse(string nameOutputFile) {

  TFile *fileOutput = new TFile(nameOutputFile.c_str(),"RECREATE");
  for(int i = 1;i < (int) fMinPt.size() - 1;i++){
    /*if(i == 1){
     fHistNVarHE -> GetAxis(0) -> SetRange(fHistNVarHE -> GetAxis(0) -> FindBin(fMinPt[i]),fHistNVarHE -> GetAxis(0) -> FindBin(fMaxPt[i])); // cut in pT
     fHistNVarCS -> GetAxis(0) -> SetRange(fHistNVarCS -> GetAxis(0) -> FindBin(fMinPt[i]),fHistNVarCS -> GetAxis(0) -> FindBin(fMaxPt[i])); // cut in pT
     cout << fHistNVarHE -> GetAxis(0) -> FindBin(fMinPt[i]) << " - " << fHistNVarHE -> GetAxis(0) -> FindBin(fMaxPt[i]) << endl;
   }
   else{
     fHistNVarHE -> GetAxis(0) -> SetRange(fHistNVarHE -> GetAxis(0) -> FindBin(fMinPt[i] + 1.01*(fHistNVarHE -> GetAxis(0) -> GetBinWidth(1))),fHistNVarHE -> GetAxis(0) -> FindBin(fMaxPt[i])); // cut in pT
     fHistNVarCS -> GetAxis(0) -> SetRange(fHistNVarCS -> GetAxis(0) -> FindBin(fMinPt[i] + 1.01*(fHistNVarCS -> GetAxis(0) -> GetBinWidth(1))),fHistNVarCS -> GetAxis(0) -> FindBin(fMaxPt[i])); // cut in pT
     cout << fHistNVarHE -> GetAxis(0) -> FindBin(fMinPt[i] + 1.01*(fHistNVarHE -> GetAxis(0) -> GetBinWidth(1))) << " - " << fHistNVarHE -> GetAxis(0) -> FindBin(fMaxPt[i]) << endl;
   }*/

    fHistNVarHE -> GetAxis(0) -> SetRange(fHistNVarHE -> GetAxis(0) -> FindBin(fMinPt[i]),fHistNVarHE -> GetAxis(0) -> FindBin(fMaxPt[i] - fHistNVarHE -> GetAxis(0) -> GetBinWidth(1))); // cut in pT
    fHistNVarCS -> GetAxis(0) -> SetRange(fHistNVarCS -> GetAxis(0) -> FindBin(fMinPt[i]),fHistNVarCS -> GetAxis(0) -> FindBin(fMaxPt[i] - fHistNVarCS -> GetAxis(0) -> GetBinWidth(1))); // cut in pT
    //cout << fHistNVarHE -> GetAxis(0) -> FindBin(fMinPt[i]) << " - " << fHistNVarCS -> GetAxis(0) -> FindBin(fMaxPt[i] - fHistNVarCS -> GetAxis(0) -> GetBinWidth(1)) << endl;

    TH1D *histMassHE = (TH1D*) fHistNVarHE -> Projection(1); histMassHE -> Write(Form("histHE_%2.1f_pT_%2.1f",fMinPt[i],fMaxPt[i])); delete histMassHE;
    TH1D *histMassCS = (TH1D*) fHistNVarCS -> Projection(1); histMassCS -> Write(Form("histCS_%2.1f_pT_%2.1f",fMinPt[i],fMaxPt[i])); delete histMassCS;

    for(int j = 0;j < fNCostBins;j++){
      //cout << fCostValues[j] << " " << fCostValues[j+1] << endl;

      THnSparse *histNVarHEClone = (THnSparse*) fHistNVarHE -> Clone("histNVarHEClone");
      histNVarHEClone -> GetAxis(2) -> SetRange(fCostBinsMin[j],fCostBinsMax[j]); // cut in CosTheta
      //cout << fCostBinsMin[j] << " - " << fCostBinsMax[j] << endl;
      TH1D *histMassCosThetaHE = (TH1D*) histNVarHEClone -> Projection(1);
      histMassCosThetaHE -> Write(Form("histHE_%2.1f_pT_%2.1f__%3.2f_CosTheta_%3.2f",fMinPt[i],fMaxPt[i],fCostValues[j],fCostValues[j+1]));
      delete histMassCosThetaHE;
      delete histNVarHEClone;

      THnSparse *histNVarCSClone = (THnSparse*) fHistNVarCS -> Clone("histNVarCSClone");
      histNVarCSClone -> GetAxis(2) -> SetRange(fCostBinsMin[j],fCostBinsMax[j]); // cut in CosTheta
      TH1D *histMassCosThetaCS = (TH1D*) histNVarCSClone -> Projection(1);
      histMassCosThetaCS -> Write(Form("histCS_%2.1f_pT_%2.1f__%3.2f_CosTheta_%3.2f",fMinPt[i],fMaxPt[i],fCostValues[j],fCostValues[j+1]));
      delete histMassCosThetaCS;
      delete histNVarCSClone;
    }

    for(int j = 0;j < fNPhiBins;j++){
      //cout << fPhiValues[j] << " " << fPhiValues[j+1] << endl;

      THnSparse *histNVarHEClone = (THnSparse*) fHistNVarHE -> Clone("histNVarHEClone");
      histNVarHEClone -> GetAxis(3) -> SetRange(fPhiBinsMin[j],fPhiBinsMax[j]); // cut in Phi
      TH1D *histMassPhiHE = (TH1D*) histNVarHEClone -> Projection(1);
      histMassPhiHE -> Write(Form("histHE_%2.1f_pT_%2.1f__%3.2f_Phi_%3.2f",fMinPt[i],fMaxPt[i],fPhiValues[j],fPhiValues[j+1]));
      delete histMassPhiHE;
      delete histNVarHEClone;

      THnSparse *histNVarCSClone = (THnSparse*) fHistNVarCS -> Clone("histNVarCSClone");
      histNVarCSClone -> GetAxis(3) -> SetRange(fPhiBinsMin[j],fPhiBinsMax[j]); // cut in Phi
      TH1D *histMassPhiCS = (TH1D*) histNVarCSClone -> Projection(1);
      histMassPhiCS -> Write(Form("histCS_%2.1f_pT_%2.1f__%3.2f_Phi_%3.2f",fMinPt[i],fMaxPt[i],fPhiValues[j],fPhiValues[j+1]));
      delete histMassPhiCS;
      delete histNVarCSClone;
    }

    for(int j = 0;j < fNCostBins;j++){
      for(int k = 0;k < fNPhiBins;k++){
        THnSparse *histNVarHEClone = (THnSparse*) fHistNVarHE -> Clone("histNVarHEClone");
        histNVarHEClone -> GetAxis(2) -> SetRange(fCostBinsMin[j],fCostBinsMax[j]); // cut in CosTheta
        histNVarHEClone -> GetAxis(3) -> SetRange(fPhiBinsMin[k],fPhiBinsMax[k]); // cut in Phi
        TH1D *histMassCosThetaPhiHE = (TH1D*) histNVarHEClone -> Projection(1);
        histMassCosThetaPhiHE -> Write(Form("histHE_%2.1f_pT_%2.1f__%3.2f_CosTheta_%3.2f__%3.2f_Phi_%3.2f",fMinPt[i],fMaxPt[i],fCostValues[j],fCostValues[j+1],fPhiValues[k],fPhiValues[k+1]));
        delete histMassCosThetaPhiHE;
        delete histNVarHEClone;

        THnSparse *histNVarCSClone = (THnSparse*) fHistNVarCS -> Clone("histNVarCSClone");
        histNVarCSClone -> GetAxis(2) -> SetRange(fCostBinsMin[j],fCostBinsMax[j]); // cut in CosTheta
        histNVarCSClone -> GetAxis(3) -> SetRange(fPhiBinsMin[k],fPhiBinsMax[k]); // cut in Phi
        TH1D *histMassCosThetaPhiCS = (TH1D*) histNVarCSClone -> Projection(1);
        histMassCosThetaPhiCS -> Write(Form("histCS_%2.1f_pT_%2.1f__%3.2f_CosTheta_%3.2f__%3.2f_Phi_%3.2f",fMinPt[i],fMaxPt[i],fCostValues[j],fCostValues[j+1],fPhiValues[k],fPhiValues[k+1]));
        delete histMassCosThetaPhiCS;
        delete histNVarCSClone;
      }
    }

    for(int j = 0;j < fNPhiTildeBins;j++){
      //cout << fPhiTildeValues[j] << " " << fPhiTildeValues[j+1] << endl;

      THnSparse *histNVarHEClone = (THnSparse*) fHistNVarHE -> Clone("histNVarHEClone");
      histNVarHEClone -> GetAxis(4) -> SetRange(fPhiTildeBinsMin[j],fPhiTildeBinsMax[j]); // cut in PhiTilde
      TH1D *histMassPhiTildeHE = (TH1D*) histNVarHEClone -> Projection(1);
      histMassPhiTildeHE -> Write(Form("histHE_%2.1f_pT_%2.1f__%3.2f_PhiTilde_%3.2f",fMinPt[i],fMaxPt[i],fPhiTildeValues[j],fPhiTildeValues[j+1]));
      delete histMassPhiTildeHE;
      delete histNVarHEClone;

      THnSparse *histNVarCSClone = (THnSparse*) fHistNVarCS -> Clone("histNVarCSClone");
      histNVarCSClone -> GetAxis(4) -> SetRange(fPhiTildeBinsMin[j],fPhiTildeBinsMax[j]); // cut in PhiTilde
      TH1D *histMassPhiTildeCS = (TH1D*) histNVarCSClone -> Projection(1);
      histMassPhiTildeCS -> Write(Form("histCS_%2.1f_pT_%2.1f__%3.2f_PhiTilde_%3.2f",fMinPt[i],fMaxPt[i],fPhiTildeValues[j],fPhiTildeValues[j+1]));
      delete histMassPhiTildeCS;
      delete histNVarCSClone;
    }
  }
  fileOutput -> Close();
}
//______________________________________________________________________________
void DataProcessor::CreateFilteredTree(string strSample, Bool_t pDCAapplied, string nameOutputFile) {
  printf("_____________Configure Fiducial Box_______________ \n");
  if(pDCAapplied){
    printf("--- pDCA included ---\n");
    fTreeData -> SetBranchAddress("pDCA",fPDCA); // enable pDDCAs
  }
  else{printf("--- pDCA not included ---\n");}
  printf("%f < CosTheta < %f \n",fCostValues[1],fCostValues[fNCostBins-1]);
  printf("%f < |Phi| < %f \n",fPhiValues[1],fPhiValues[fNPhiBins-1]);
  printf("______________________________________________________________ \n");

  double PI = TMath::Pi();
  int nEvents = 0;
  if(strSample == "FullStat"){nEvents = fTreeData -> GetEntries();}
  if(strSample == "TestStat"){nEvents = 100000;}
  printf("N events = %i \n",nEvents);

  double tmpVar = 0;    // variable used to store the value of Phi for PhiTilde calculation

  TFile *fileTreeDataFiltered = new TFile(nameOutputFile.c_str(),"RECREATE");

  double DimuMass, DimuPt, CosThetaHE, CosThetaCS, PhiHE, PhiCS,PhiTildeHE, PhiTildeCS;

  TTree *treeDataFilteredHE = new TTree("treeDataFilteredHE","treeDataFilteredHE");
  treeDataFilteredHE -> Branch("DimuPt",&DimuPt,"DimuPt/D");
  treeDataFilteredHE -> Branch("DimuMass",&DimuMass,"DimuMass/D");
  treeDataFilteredHE -> Branch("CosThetaHE",&CosThetaHE,"CosThetaHE/D");
  treeDataFilteredHE -> Branch("PhiHE",&PhiHE,"PhiHE/D");
  treeDataFilteredHE -> Branch("PhiTildeHE",&PhiTildeHE,"PhiTildeHE/D");

  TTree *treeDataFilteredCS = new TTree("treeDataFilteredCS","treeDataFilteredCS");
  treeDataFilteredCS -> Branch("DimuPt",&DimuPt,"DimuPt/D");
  treeDataFilteredCS -> Branch("DimuMass",&DimuMass,"DimuMass/D");
  treeDataFilteredCS -> Branch("CosThetaCS",&CosThetaCS,"CosThetaCS/D");
  treeDataFilteredCS -> Branch("PhiCS",&PhiCS,"PhiCS/D");
  treeDataFilteredCS -> Branch("PhiTildeCS",&PhiTildeCS,"PhiTildeCS/D");

  bool goodMuon1 = kFALSE;
  bool goodMuon2 = kFALSE;

  for(int i = 0;i < nEvents;i++){
    fTreeData -> GetEntry(i);
    printf("Reading : %2.1f %% \r",((double) i)/((double) nEvents)*100);

    for(int j = 0;j < fNDimu;j++){

      if(fIsPhysSelected){
        TString Trigger = fTrigClass;
        Bool_t TriggerSelected = kFALSE;
        if(Trigger.Contains("CMUL7-B-NOPF-MUFAST")) TriggerSelected = kTRUE;
        if(fDimuY[j] > -4. && fDimuY[j] < -2.5){
          if(TriggerSelected){
            if(fDimuMatch[j] == 2){
              if(fDimuMass[j] > 2 && fDimuMass[j] < 5){
                for(int k = 0;k < fNMuons;k++){
                  if(fMuonId[k] == fDimuMu[j][0]){
                    if(fPDCA[k] != 0){goodMuon1 = kTRUE;}
                  }
                  if(fMuonId[k] == fDimuMu[j][1]){
                    if(fPDCA[k] != 0){goodMuon2 = kTRUE;}
                  }
                }
                if(pDCAapplied == kFALSE){
                  goodMuon1 = kTRUE;
                  goodMuon2 = kTRUE;
                }
                if(goodMuon1 == kTRUE && goodMuon2 == kTRUE){
                  if(TMath::Abs(fPhiHE[j]) > fPhiValues[1] && TMath::Abs(fPhiHE[j]) < fPhiValues[fNPhiBins-1]){
                    if(fCostHE[j] > fCostValues[1] && fCostHE[j] < fCostValues[fNCostBins-1]){
                      DimuPt = fDimuPt[j];
                      DimuMass = fDimuMass[j];
                      CosThetaHE = fCostHE[j];
                      PhiHE = TMath::Abs(fPhiHE[j]);
                      tmpVar = fPhiHE[j] + PI;
                      if(fCostHE[j] < 0.){fPhiTildeHE[j] = tmpVar - (3./4.)*PI;}
                      if(fCostHE[j] > 0.){fPhiTildeHE[j] = tmpVar - (1./4.)*PI;}
                      if(fPhiTildeHE[j] > 2*PI){fPhiTildeHE[j] = fPhiTildeHE[j] - 2*PI;}
                      if(fPhiTildeHE[j] < 0.){fPhiTildeHE[j] = 2*PI + fPhiTildeHE[j];}
                      PhiTildeHE = fPhiTildeHE[j];
                      treeDataFilteredHE -> Fill();                                 // Filling the THnSparse
                    }
                  }

                  if(TMath::Abs(fPhiCS[j]) > fPhiValues[1] && TMath::Abs(fPhiCS[j]) < fPhiValues[fNPhiBins-1]){
                    if(fCostCS[j] > fCostValues[1] && fCostCS[j] < fCostValues[fNCostBins-1]){
                      DimuPt = fDimuPt[j];
                      DimuMass = fDimuMass[j];
                      CosThetaCS = fCostCS[j];
                      PhiCS = TMath::Abs(fPhiCS[j]);
                      tmpVar = fPhiCS[j] + PI;
                      if(fCostCS[j] < 0.){fPhiTildeCS[j] = tmpVar - (3./4.)*PI;}
                      if(fCostCS[j] > 0.){fPhiTildeCS[j] = tmpVar - (1./4.)*PI;}
                      if(fPhiTildeCS[j] > 2*PI){fPhiTildeCS[j] = fPhiTildeCS[j] - 2*PI;}
                      if(fPhiTildeCS[j] < 0.){fPhiTildeCS[j] = 2*PI + fPhiTildeCS[j];}
                      PhiTildeCS = fPhiTildeCS[j];
                      treeDataFilteredCS -> Fill();                                 // Filling the THnSparse
                    }
                  }
                }
                goodMuon1 = kFALSE;
                goodMuon2 = kFALSE;
              }
            }
          }
        }
      }
    }
  }
  printf("\n");

  fileTreeDataFiltered -> cd();
  treeDataFilteredHE -> Write();
  treeDataFilteredCS -> Write();
  fileTreeDataFiltered -> Close();
}
//______________________________________________________________________________
void DataProcessor::CutFilteredTree(string nameOutputFile) {
  printf("Colmpleta l'opera \n");

  double DimuMass, DimuPt, CosThetaHE, CosThetaCS, PhiHE, PhiCS,PhiTildeHE, PhiTildeCS;

  double PI = TMath::Pi();
  int nEvents = fTreeDataFilteredHE -> GetEntries();

  fTreeDataFilteredHE -> SetBranchAddress("DimuPt",&DimuPt);
  fTreeDataFilteredHE -> SetBranchAddress("DimuMass",&DimuMass);
  fTreeDataFilteredHE -> SetBranchAddress("CosThetaHE",&CosThetaHE);
  fTreeDataFilteredHE -> SetBranchAddress("PhiHE",&PhiHE);
  fTreeDataFilteredHE -> SetBranchAddress("PhiTildeHE",&PhiTildeHE);

  TH1D *histMass_2pT4 = new TH1D("histMass_2pT4","histMass_2pT4",120,2.,5.);
  TH1D *histMass_4pT6 = new TH1D("histMass_4pT6","histMass_4pT6",120,2.,5.);
  TH1D *histMass_6pT10 = new TH1D("histMass_6pT10","histMass_6pT10",120,2.,5.);

  TH1D *histMass_CosTheta = new TH1D("histMass_CosTheta","histMass_CosTheta",120,2.,5.);
  TH1D *histMass_Phi = new TH1D("histMass_Phi","histMass_Phi",120,2.,5.);

  for(int i = 0;i < nEvents;i++){
    fTreeDataFilteredHE -> GetEntry(i);
    if(DimuPt > 2. && DimuPt < 4.){histMass_2pT4 -> Fill(DimuMass);}
    if(DimuPt > 4. && DimuPt < 6.){histMass_4pT6 -> Fill(DimuMass);}
    if(DimuPt > 6. && DimuPt < 10.){
      histMass_6pT10 -> Fill(DimuMass);
      if(CosThetaHE > fCostValues[1] && CosThetaHE < fCostValues[2]){
        histMass_CosTheta -> Fill(DimuMass);
      }
      if(TMath::Abs(PhiHE) > fPhiValues[1] && TMath::Abs(PhiHE) < fPhiValues[2]){
        histMass_Phi -> Fill(DimuMass);
      }
    }
  }

  printf("N entries 2 < pT < 4 GeV/c = %i \n",(int) histMass_2pT4 -> GetEntries());
  printf("N entries 4 < pT < 6 GeV/c = %i \n",(int)histMass_4pT6 -> GetEntries());
  printf("N entries 6 < pT < 10 GeV/c = %i \n",(int)histMass_6pT10 -> GetEntries());

  TFile *fileOutput = new TFile(nameOutputFile.c_str(),"RECREATE");
  histMass_2pT4 -> Write();
  histMass_4pT6 -> Write();
  histMass_6pT10 -> Write();

  histMass_CosTheta -> Write();
  histMass_Phi -> Write();
  fileOutput -> Close();
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

  int counterCMUL7 = 0;
  vector <int> fListMuonId;

  for(int i = 0;i < nEvents;i++){
    fTreeData -> GetEntry(i);
    //printf("Reading : %2.1f %% \r",((double) i)/((double) nEvents)*100);
    TString Trigger = fTrigClass;
    Bool_t TriggerSelected = kFALSE;
    if(Trigger.Contains("CMUL7-B-NOPF-MUFAST")){counterCMUL7++;}
    if(Trigger.Contains("CINT7-B-NOPF-MUFAST")){TriggerSelected = kTRUE;} // single muon trigger
    if(TriggerSelected){
      if(fIsPhysSelected){       
      printf("------------------------------------------------------------\n");                                                          
      for(int k = 0;k < fNDimu;k++){
        if(fDimuY[k] > -4. && fDimuY[k] < -2.5){
            //if(fDimuMatch[k] == 2){                                                        //*****
              if(fDimuMass[k] > 2 && fDimuMass[k] < 5){
                fListMuonId.push_back(fDimuMu[k][0]);
                fListMuonId.push_back(fDimuMu[k][1]);
                printf("%i - %i (Dimu Px = %f) \n",fDimuMu[k][0],fDimuMu[k][1],fDimuPx[k]);
              }
            //}                                                                              //*****
          }
        }

        //printf("=== MUON COMPARISON === \n");
        for(int k = 0;k < fNMuons;k++){
          //printf("%i ",MuonId[k]);

          for(int j = 0;j < (int) fListMuonId.size();j++){
            if(fMuonId[k] == fListMuonId[j]){
              printf("%i -> ok! (Px = %f) \n",fMuonId[k],fPx[k]);

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

  //TH1D *fHistTriggerResponseFunctionSM = new TH1D("histTriggerResponseFunctionSM","",100,0,10);
  //fHistTriggerResponseFunctionSM -> Divide(fHistLowPtSM,fHistAllPtSM,1,1,"B");
  //fHistTriggerResponseFunctionSM -> SetLineColor(kBlue);

  //TH1D *fHistTriggerResponseFunctionSMpDCA = new TH1D("histTriggerResponseFunctionSMpDCA","",100,0,10);
  //fHistTriggerResponseFunctionSMpDCA -> Divide(fHistLowPtSMpDCA,fHistAllPtSMpDCA,1,1,"B");
  //fHistTriggerResponseFunctionSMpDCA -> SetLineColor(kRed);
  //fHistTriggerResponseFunctionSMpDCA -> SetMarkerStyle(20);
  //fHistTriggerResponseFunctionSMpDCA -> SetMarkerColor(kRed);

  fHistCMUL7Triggers = new TH1D("fHistCMUL7Triggers","",1,0,1);
  fHistCMUL7Triggers -> SetBinContent(1,counterCMUL7);

  TFile *fileTriggerResponseFunction = new TFile(nameOutputFile.c_str(),"RECREATE");
  fHistAllPtSM -> Write();
  fHistLowPtSM -> Write();
  fHistAllPtSMpDCA -> Write();
  fHistLowPtSMpDCA -> Write();
  fHistCMUL7Triggers -> Write();
  fileTriggerResponseFunction -> Close();
}
