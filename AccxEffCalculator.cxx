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
#include "AccxEffCalculator.h"

ClassImp(AccxEffCalculator)

//______________________________________________________________________________
AccxEffCalculator::AccxEffCalculator(): TObject() {
   // default constructor
}
//______________________________________________________________________________
AccxEffCalculator::AccxEffCalculator(TTree *treeAccxEff): TObject() {
  // standard constructor
  fTreeAccxEff = (TTree*) treeAccxEff -> Clone();
  fTreeAccxEff -> SetBranchAddress("NDimu_gen",&fNDimuGen);
  fTreeAccxEff -> SetBranchAddress("DimuPt_gen",fDimuPtGen);
  fTreeAccxEff -> SetBranchAddress("DimuY_gen",fDimuYGen);
  fTreeAccxEff -> SetBranchAddress("CostHE_gen",fCostHEGen);
  fTreeAccxEff -> SetBranchAddress("PhiHE_gen",fPhiHEGen);
  fTreeAccxEff -> SetBranchAddress("CostCS_gen",fCostCSGen);
  fTreeAccxEff -> SetBranchAddress("PhiCS_gen",fPhiCSGen);

  fTreeAccxEff -> SetBranchAddress("NDimu_rec",&fNDimuRec);
  fTreeAccxEff -> SetBranchAddress("DimuPt_rec",fDimuPtRec);
  fTreeAccxEff -> SetBranchAddress("DimuY_rec",fDimuYRec);
  fTreeAccxEff -> SetBranchAddress("DimuMass_rec",fDimuMassRec);
  fTreeAccxEff -> SetBranchAddress("DimuMatch_rec",fDimuMatchRec);
  fTreeAccxEff -> SetBranchAddress("CostHE_rec",fCostHERec);
  fTreeAccxEff -> SetBranchAddress("PhiHE_rec",fPhiHERec);
  fTreeAccxEff -> SetBranchAddress("CostCS_rec",fCostCSRec);
  fTreeAccxEff -> SetBranchAddress("PhiCS_rec",fPhiCSRec);

  fTreeAccxEff -> SetBranchAddress("NMuons_rec",&fNMuonsRec);
  fTreeAccxEff -> SetBranchAddress("Pt_rec",fPtRec);
  fTreeAccxEff -> SetBranchAddress("MatchTrig_rec",fMatchTrigRec);
}
//______________________________________________________________________________
AccxEffCalculator::~AccxEffCalculator() {
  // destructor
}
//______________________________________________________________________________
void AccxEffCalculator::SetPtBins(Int_t nPtBins, Double_t minPtBin[],Double_t maxPtBin[]) {
  fNPtBins = nPtBins;
  for(int i = 0;i < fNPtBins;i++){
    fMinPt.push_back(minPtBin[i]);
    fMaxPt.push_back(maxPtBin[i]);
  }
}
//______________________________________________________________________________
void AccxEffCalculator::SetBinning(vector <Double_t> CostValues, vector <Double_t> PhiValues) {
  for(int i = 0;i < (int) CostValues.size();i++){fCostValues.push_back(CostValues[i]);}
  fNCostBins = fCostValues.size() - 1;
  for(int i = 0;i < (int) PhiValues.size();i++){fPhiValues.push_back(PhiValues[i]);}
  fNPhiBins = fPhiValues.size() - 1;
}
//______________________________________________________________________________
void AccxEffCalculator::ComputeAccxEff(string strSample, string nameOutputFile) {
  int nEvents = 0;
  int indexPt = 0;
  double PI = TMath::Pi();

  if(strSample == "FullStat"){nEvents = fTreeAccxEff -> GetEntries();}
  if(strSample == "TestStat"){nEvents = 1000000;}
  printf("N events = %i \n",nEvents);

  for(int i = 0;i < fNPtBins;i++){
    fHistGenCost[i] = new TH1D(Form("histGenCost_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCostBins,&fCostValues[0]);
    fHistGenCost[i] -> Sumw2();

    fHistRecCost[i] = new TH1D(Form("histRecCost_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCostBins,&fCostValues[0]);
    fHistRecCost[i] -> Sumw2();

    fHistGenPhi[i] = new TH1D(Form("histGenPhi_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistGenPhi[i] -> Sumw2();

    fHistRecPhi[i] = new TH1D(Form("histRecPhi_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistRecPhi[i] -> Sumw2();

    fHistGenCostPhi[i] = new TH2D(Form("histGenCostPhi_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCostBins,&fCostValues[0],fNPhiBins,&fPhiValues[0]);
    fHistGenCostPhi[i] -> Sumw2();

    fHistRecCostPhi[i] = new TH2D(Form("histRecCostPhi_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCostBins,&fCostValues[0],fNPhiBins,&fPhiValues[0]);
    fHistRecCostPhi[i] -> Sumw2();
  }

  fHistGenCostPt = new TH2D("histGenCostPt","",100,-1,1,100,0,10);
  fHistGenCostPt -> Sumw2();

  fHistRecCostPt = new TH2D("histRecCostPt","",100,-1,1,100,0,10);
  fHistRecCostPt -> Sumw2();

  fHistGenPhiPt = new TH2D("histGenPhiPt","",100,0,PI,100,0,10);
  fHistGenPhiPt -> Sumw2();

  fHistRecPhiPt = new TH2D("histRecPhiPt","",100,0,PI,100,0,10);
  fHistRecPhiPt -> Sumw2();

  for(int i = 0;i < nEvents;i++){
    printf("Reading : %3.2f %% \r",((double) i/(double) nEvents)*100.);
    fTreeAccxEff -> GetEntry(i);

    for(int j = 0;j < fNDimuGen;j++){
      if(fDimuYGen[j] > -4. && fDimuYGen[j] < -2.5){
        fHistGenCostPt -> Fill(fCostHEGen[j],fDimuPtGen[j]);
        fHistGenPhiPt -> Fill(TMath::Abs(fPhiHEGen[j]),fDimuPtGen[j]);
        while(fDimuPtGen[j] < fMinPt[indexPt] || fDimuPtGen[j] > fMaxPt[indexPt]){indexPt++;}
        fHistGenCost[indexPt] -> Fill(fCostHEGen[j]);
        fHistGenPhi[indexPt] -> Fill(TMath::Abs(fPhiHEGen[j]));
        fHistGenCostPhi[indexPt] -> Fill(fCostHEGen[j],TMath::Abs(fPhiHEGen[j]));
        indexPt = 0;
      }
    }

    for(int j = 0;j < fNDimuRec;j++){
      if(fDimuYRec[j] > -4. && fDimuYRec[j] < -2.5){
        if(fDimuMatchRec[j] == 2){
          if(fDimuMassRec[j] > 2 && fDimuMassRec[j] < 5){
            fHistRecCostPt -> Fill(fCostHERec[j],fDimuPtRec[j]);
            fHistRecPhiPt -> Fill(TMath::Abs(fPhiHERec[j]),fDimuPtRec[j]);
            while(fDimuPtRec[j] < fMinPt[indexPt] || fDimuPtRec[j] > fMaxPt[indexPt]){indexPt++;}
            fHistRecCost[indexPt] -> Fill(fCostHERec[j]);
            fHistRecPhi[indexPt] -> Fill(TMath::Abs(fPhiHERec[j]));
            fHistRecCostPhi[indexPt] -> Fill(fCostHERec[j],TMath::Abs(fPhiHERec[j]));
            indexPt = 0;
          }
        }
      }
    }

  }
  printf("\n");

  for(int i = 0;i < fNPtBins;i++){
    fHistAccxEffCost[i] = new TH1D(Form("histAccxEffCost_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCostBins,&fCostValues[0]);
    fHistAccxEffCost[i] -> Divide(fHistRecCost[i],fHistGenCost[i],1,1,"B");

    fHistAccxEffPhi[i] = new TH1D(Form("histAccxEffPhi_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistAccxEffPhi[i] -> Divide(fHistRecPhi[i],fHistGenPhi[i],1,1,"B");

    fHistAccxEffCostPhi[i] = new TH2D(Form("histAccxEffCostPhi_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCostBins,&fCostValues[0],fNPhiBins,&fPhiValues[0]);
    fHistAccxEffCostPhi[i] -> Divide(fHistRecCostPhi[i],fHistGenCostPhi[i],1,1,"B");

    fHistAccxEffCostPhiStatRel[i] = new TH2D(Form("histAccxEffCostPhiStatRel_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCostBins,&fCostValues[0],fNPhiBins,&fPhiValues[0]);
  }

  fHistAccxEffCostPt = new TH2D("histAccxEffCostPt","",100,-1,1,100,0,10);
  fHistAccxEffCostPt -> Divide(fHistRecCostPt,fHistGenCostPt,1,1,"B");
  fHistAccxEffPhiPt = new TH2D("histAccxEffPhiPt","",100,0,PI,100,0,10);
  fHistAccxEffPhiPt -> Divide(fHistRecPhiPt,fHistGenPhiPt,1,1,"B");

  for(int i = 0;i < fNPtBins;i++){
    for(int j = 0;j < fNCostBins;j++){
      for(int k = 0;k < fNPhiBins;k++){
        fHistAccxEffCostPhiStatRel[i] -> SetBinContent(j+1,k+1,(fHistAccxEffCostPhi[i] -> GetBinError(j+1,k+1)/fHistAccxEffCostPhi[i] -> GetBinContent(j+1,k+1))*100);
      }
    }
  }

  TFile *fileAccxEff = new TFile(nameOutputFile.c_str(),"RECREATE");
  for(int i = 0;i < fNPtBins;i++){
    fHistGenCost[i] -> Write();
    fHistRecCost[i] -> Write();
    fHistAccxEffCost[i] -> Write();
    fHistGenPhi[i] -> Write();
    fHistRecPhi[i] -> Write();
    fHistAccxEffPhi[i] -> Write();
    fHistGenCostPhi[i] -> Write();
    fHistRecCostPhi[i] -> Write();
    fHistAccxEffCostPhi[i] -> Write();
    fHistAccxEffCostPhiStatRel[i] -> Write();
  }
  fHistGenCostPt -> Write();
  fHistGenPhiPt -> Write();
  fHistAccxEffCostPt -> Write();
  fHistRecCostPt -> Write();
  fHistRecPhiPt -> Write();
  fHistAccxEffPhiPt -> Write();
  fileAccxEff -> Close();
}
//______________________________________________________________________________
void AccxEffCalculator::ReWeightAccxEff(Double_t LambdaTheta, string strSample, string nameOutputFile) {
  TF1 *funcCosth = new TF1("funcCosth","(1/(3 + [0]))*(1 + [0]*x*x)",-1,1);
  funcCosth -> SetParameter(0,LambdaTheta);

  int nEvents = 0;
  int indexPt = 0;

  if(strSample == "FullStat"){nEvents = fTreeAccxEff -> GetEntries();}
  if(strSample == "TestStat"){nEvents = 100000;}
  printf("N events = %i \n",nEvents);

  double weight = 0;

  for(int i = 0;i < fNPtBins;i++){
    fHistGenCostReWeighted[i] = new TH1D(Form("histGenCostReWeighted%2.1f_%ipT%i",LambdaTheta, (int) fMinPt[i],(int) fMaxPt[i]),"",fNCostBins,&fCostValues[0]);
    fHistGenCostReWeighted[i] -> Sumw2();

    fHistRecCostReWeighted[i] = new TH1D(Form("histRecCostReWeighted%2.1f_%ipT%i",LambdaTheta, (int) fMinPt[i],(int) fMaxPt[i]),"",fNCostBins,&fCostValues[0]);
    fHistRecCostReWeighted[i] -> Sumw2();

    fHistGenPhiReWeighted[i] = new TH1D(Form("histGenPhiReWeighted%2.1f_%ipT%i",LambdaTheta, (int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistGenPhiReWeighted[i] -> Sumw2();

    fHistRecPhiReWeighted[i] = new TH1D(Form("histRecPhiReWeighted%2.1f_%ipT%i",LambdaTheta, (int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistRecPhiReWeighted[i] -> Sumw2();
  }

  for(int i = 0;i < nEvents;i++){
    printf("Reading : %3.2f %% \r",((double) i/(double) nEvents)*100.);
    fTreeAccxEff -> GetEntry(i);

    for(int j = 0;j < fNDimuGen;j++){
      if(fDimuYGen[j] > -4. && fDimuYGen[j] < -2.5){
        weight = (funcCosth -> Eval(fCostHEGen[j]))/(funcCosth -> GetMaximum());
        while(fDimuPtGen[j] < fMinPt[indexPt] || fDimuPtGen[j] > fMaxPt[indexPt]){
          indexPt++;
        }
        fHistGenCostReWeighted[indexPt] -> Fill(fCostHEGen[j],weight);
        fHistGenPhiReWeighted[indexPt] -> Fill(TMath::Abs(fPhiHEGen[j]),weight);
        indexPt = 0;
      }
    }

    for(int j = 0;j < fNDimuRec;j++){
      if(fDimuYRec[j] > -4. && fDimuYRec[j] < -2.5){
        if(fDimuMatchRec[j] == 2){
          if(fDimuMassRec[j] > 2 && fDimuMassRec[j] < 5){
            while(fDimuPtRec[j] < fMinPt[indexPt] || fDimuPtRec[j] > fMaxPt[indexPt]){
              indexPt++;
            }
            fHistRecCostReWeighted[indexPt] -> Fill(fCostHERec[j],weight);
            fHistRecPhiReWeighted[indexPt] -> Fill(TMath::Abs(fPhiHERec[j]),weight);
            indexPt = 0;
          }
        }
      }
    }

  }
  printf("\n");

  for(int i = 0;i < fNPtBins;i++){
    fHistAccxEffCostReWeighted[i] = new TH1D(Form("histAccxEffCostReWeighted%2.1f_%ipT%i",LambdaTheta,(int) fMinPt[i],(int) fMaxPt[i]),"",fNCostBins,&fCostValues[0]);
    fHistAccxEffCostReWeighted[i] -> Divide(fHistRecCostReWeighted[i],fHistGenCostReWeighted[i],1,1,"B");

    fHistAccxEffPhiReWeighted[i] = new TH1D(Form("histAccxEffPhiReWeighted%2.1f_%ipT%i",LambdaTheta,(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistAccxEffPhiReWeighted[i] -> Divide(fHistRecPhiReWeighted[i],fHistGenPhiReWeighted[i],1,1,"B");
  }

  TFile *fileAccxEffReWeighted = new TFile(nameOutputFile.c_str(),"RECREATE");
  for(int i = 0;i < fNPtBins;i++){
    fHistGenCostReWeighted[i] -> Write();
    fHistRecCostReWeighted[i] -> Write();
    fHistAccxEffCostReWeighted[i] -> Write();

    fHistGenPhiReWeighted[i] -> Write();
    fHistRecPhiReWeighted[i] -> Write();
    fHistAccxEffPhiReWeighted[i] -> Write();
  }
  fileAccxEffReWeighted -> Close();
}
//______________________________________________________________________________
void AccxEffCalculator::ComputeTriggerResponseFunction(string strSample, string nameOutputFile) {
  int nEvents = 0;
  int indexPt = 0;

  if(strSample == "FullStat"){nEvents = fTreeAccxEff -> GetEntries();}
  if(strSample == "TestStat"){nEvents = 100000;}
  printf("N events = %i \n",nEvents);

  fHistLowPt = new TH1D("fHistLowPt","",100,0,10);
  fHistLowPt -> Sumw2();
  fHistAllPt = new TH1D("fHistAllPt","",100,0,10);
  fHistAllPt -> Sumw2();

  for(int i = 0;i < nEvents;i++){
    fTreeAccxEff -> GetEntry(i);
    printf("Reading : %2.1f %% \r",((double) i/(double) nEvents)*100.);
    if(fMatchTrigRec[0] >= 1){fHistAllPt -> Fill(fPtRec[0]);}
    if(fMatchTrigRec[1] >= 1){fHistAllPt -> Fill(fPtRec[1]);}
    if(fMatchTrigRec[0] >= 2){fHistLowPt -> Fill(fPtRec[0]);}
    if(fMatchTrigRec[1] >= 2){fHistLowPt -> Fill(fPtRec[1]);}
  }
  printf("\n");

  TH1D *fHistTriggerResponseFunction = new TH1D("histTriggerResponseFunction","",100,0,10);
  fHistTriggerResponseFunction -> Divide(fHistLowPt,fHistAllPt,1,1,"B");

  TFile *fileTriggerResponseFunction = new TFile(nameOutputFile.c_str(),"RECREATE");
  fHistLowPt -> Write();
  fHistAllPt -> Write();
  fHistTriggerResponseFunction -> Write();
  fileTriggerResponseFunction -> Close();
}
