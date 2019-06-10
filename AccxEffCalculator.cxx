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
#include "TF2.h"
#include "TFile.h"
#include "THnSparse.h"
#include "AccxEffCalculator.h"

Double_t MyFuncPol(Double_t *, Double_t *);
Double_t computePhiTilde(Double_t , Double_t);

ClassImp(AccxEffCalculator)

//______________________________________________________________________________
AccxEffCalculator::AccxEffCalculator(): TObject() {
  fPi = TMath::Pi();
  // default constructor
}
//______________________________________________________________________________
AccxEffCalculator::AccxEffCalculator(TTree *treeAccxEff): TObject() {
  // standard constructor
  fPi = TMath::Pi();
  fTreeAccxEff = (TTree*) treeAccxEff -> Clone();
  fTreeAccxEff -> SetBranchAddress("NDimu_gen",&fNDimuGen);
  fTreeAccxEff -> SetBranchAddress("DimuPt_gen",fDimuPtGen);
  fTreeAccxEff -> SetBranchAddress("DimuY_gen",fDimuYGen);
  fTreeAccxEff -> SetBranchAddress("CostHE_gen",fCosThetaHEGen);
  fTreeAccxEff -> SetBranchAddress("PhiHE_gen",fPhiHEGen);
  fTreeAccxEff -> SetBranchAddress("CostCS_gen",fCosThetaCSGen);
  fTreeAccxEff -> SetBranchAddress("PhiCS_gen",fPhiCSGen);

  fTreeAccxEff -> SetBranchAddress("NDimu_rec",&fNDimuRec);
  fTreeAccxEff -> SetBranchAddress("DimuPt_rec",fDimuPtRec);
  fTreeAccxEff -> SetBranchAddress("DimuY_rec",fDimuYRec);
  fTreeAccxEff -> SetBranchAddress("DimuMass_rec",fDimuMassRec);
  fTreeAccxEff -> SetBranchAddress("DimuMatch_rec",fDimuMatchRec);
  fTreeAccxEff -> SetBranchAddress("CostHE_rec",fCosThetaHERec);
  fTreeAccxEff -> SetBranchAddress("PhiHE_rec",fPhiHERec);
  fTreeAccxEff -> SetBranchAddress("CostCS_rec",fCosThetaCSRec);
  fTreeAccxEff -> SetBranchAddress("PhiCS_rec",fPhiCSRec);

  fTreeAccxEff -> SetBranchAddress("NMuons_rec",&fNMuonsRec);
  fTreeAccxEff -> SetBranchAddress("Pt_rec",fPtRec);
  fTreeAccxEff -> SetBranchAddress("Eta_rec",fEtaRec);
  fTreeAccxEff -> SetBranchAddress("Charge_rec",fChargeRec);
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
void AccxEffCalculator::SetBinning(vector <Double_t> CosThetaValues, vector <Double_t> PhiValues, vector <Double_t> PhiTildeValues) {
  for(int i = 0;i < (int) CosThetaValues.size();i++){fCosThetaValues.push_back(CosThetaValues[i]);}
  fNCosThetaBins = fCosThetaValues.size() - 1;
  for(int i = 0;i < (int) PhiValues.size();i++){fPhiValues.push_back(PhiValues[i]);}
  fNPhiBins = fPhiValues.size() - 1;
  for(int i = 0;i < (int) PhiTildeValues.size();i++){fPhiTildeValues.push_back(PhiTildeValues[i]);}
  fNPhiTildeBins = fPhiTildeValues.size() - 1;
}
//______________________________________________________________________________
void AccxEffCalculator::ComputeResolution(string strSample, string nameOutputFile){
  int nEvents = 0;
  int indexPt = 0;

  if(strSample == "FullStat"){nEvents = fTreeAccxEff -> GetEntries();}
  if(strSample == "TestStat"){nEvents = 1000000;}
  printf("N events = %i \n",nEvents);

  for(int i = 0;i < fNPtBins;i++){
    // HELICITY
    fHistResolutionCosThetaHE[i] = new TH1D(Form("histResolutionCosThetaHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",1000,-0.5,0.5); fHistResolutionCosThetaHE[i] -> Sumw2();
    fHistResolutionPhiHE[i] = new TH1D(Form("histResolutionPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",1000,-0.5,0.5); fHistResolutionPhiHE[i] -> Sumw2();
    fHistResolutionCosThetaPhiHE[i] = new TH2D(Form("histResolutionCosThetaPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",1000,-0.5,0.5,1000,-0.5,0.5); fHistResolutionCosThetaPhiHE[i] -> Sumw2();
    // COLLINS-SOPER
    fHistResolutionCosThetaCS[i] = new TH1D(Form("histResolutionCosThetaCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",1000,-0.5,0.5); fHistResolutionCosThetaCS[i] -> Sumw2();
    fHistResolutionPhiCS[i] = new TH1D(Form("histResolutionPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",1000,-0.5,0.5); fHistResolutionPhiCS[i] -> Sumw2();
    fHistResolutionCosThetaPhiCS[i] = new TH2D(Form("histResolutionCosThetaPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",1000,-0.5,0.5,1000,-0.5,0.5); fHistResolutionCosThetaPhiCS[i] -> Sumw2();
  }
  /*
  // HELICITY
  fHistResolutionCosThetaHE = new TH1D("histResolutionCosThetaHE","",1000,-0.5,0.5); fHistResolutionCosThetaHE -> Sumw2();
  fHistResolutionPhiHE = new TH1D("histResolutionPhiHE","",1000,-fPi/2.,fPi/2.); fHistResolutionPhiHE -> Sumw2();
  fHistResolutionCosThetaPhiHE = new TH2D("histResolutionCosThetaPhiHE","",1000,-0.5,0.5,1000,-fPi/2.,fPi/2.); fHistResolutionCosThetaPhiHE -> Sumw2();
  // COLLINS-SOPER
  fHistResolutionCosThetaCS = new TH1D("histResolutionCosThetaCS","",1000,-0.5,0.5); fHistResolutionCosThetaCS -> Sumw2();
  fHistResolutionPhiCS = new TH1D("histResolutionPhiCS","",1000,-fPi/2.,fPi/2.); fHistResolutionPhiCS -> Sumw2();
  fHistResolutionCosThetaPhiCS = new TH2D("histResolutionCosThetaPhiCS","",1000,-0.5,0.5,1000,-fPi/2.,fPi/2.); fHistResolutionCosThetaPhiCS -> Sumw2();
  */

  if(strSample == "FullStat"){nEvents = fTreeAccxEff -> GetEntries();}
  if(strSample == "TestStat"){nEvents = 5000000;}
  printf("N events = %i \n",nEvents);

  double CosThetaHEGen = 0., CosThetaCSGen = 0, CosThetaHERec = 0., CosThetaCSRec = 0;
  double PhiHEGen = 0., PhiCSGen = 0, PhiHERec = 0., PhiCSRec = 0;
  double resCosThetaHE = 0., resCosThetaCS = 0., resPhiHE = 0., resPhiCS = 0.;

  for(int i = 0;i < nEvents;i++){
    printf("Reading : %3.2f %% \r",((double) i/(double) nEvents)*100.);
    fTreeAccxEff -> GetEntry(i);

    CosThetaHEGen = 0.;
    CosThetaCSGen = 0.;
    PhiHEGen = 0.;
    PhiCSGen = 0.;
    CosThetaHERec = 0.;
    CosThetaCSRec = 0.;
    PhiHERec = 0.;
    PhiCSRec = 0.;

    resCosThetaHE = 0.;
    resCosThetaCS = 0.;
    resPhiHE = 0.;
    resPhiCS = 0.;

    resCosThetaHE = 0.;

    if(fNDimuGen != 0 && fNDimuRec != 0){

      for(int j = 0;j < fNDimuGen;j++){
        if(fDimuYGen[j] > -4. && fDimuYGen[j] < -2.5){
          CosThetaHEGen = fCosThetaHEGen[j];
          CosThetaCSGen = fCosThetaCSGen[j];
          PhiHEGen = fPhiHEGen[j];
          PhiCSGen = fPhiCSGen[j];
        }
      }

      for(int j = 0;j < fNDimuRec;j++){
        if(fDimuYRec[j] > -4. && fDimuYRec[j] < -2.5){
          CosThetaHERec = fCosThetaHERec[j];
          CosThetaCSRec = fCosThetaCSRec[j];
          PhiHERec = fPhiHERec[j];
          PhiCSRec = fPhiCSRec[j];

          /*resCosThetaHE = (CosThetaHEGen - CosThetaHERec)/(CosThetaHEGen);
          resPhiHE = (PhiHEGen - PhiHERec)/(PhiHEGen);
          resCosThetaCS = (CosThetaCSGen - CosThetaCSRec)/CosThetaCSGen;
          resPhiCS = (PhiCSGen - PhiCSRec)/PhiCSGen;*/
          resCosThetaHE = (CosThetaHEGen - CosThetaHERec);
          resPhiHE = (PhiHEGen - PhiHERec);
          resCosThetaCS = (CosThetaCSGen - CosThetaCSRec);
          resPhiCS = (PhiCSGen - PhiCSRec);

          indexPt = 0;
          if(fDimuPtGen[j] < 1000.){
            while(fDimuPtGen[j] < fMinPt[indexPt] || fDimuPtGen[j] > fMaxPt[indexPt]){indexPt++;}
            fHistResolutionCosThetaHE[indexPt] -> Fill(resCosThetaHE);
            fHistResolutionPhiHE[indexPt] -> Fill(resPhiHE);
            fHistResolutionCosThetaPhiHE[indexPt] -> Fill(resCosThetaHE,resPhiHE);

            fHistResolutionCosThetaCS[indexPt] -> Fill(resCosThetaCS);
            fHistResolutionPhiCS[indexPt] -> Fill(resPhiCS);
            fHistResolutionCosThetaPhiCS[indexPt] -> Fill(resCosThetaCS,resPhiCS);
            indexPt = 0;
          }
        }
      }


    }
  }


  TFile *fileResolution = new TFile(nameOutputFile.c_str(),"RECREATE");
  for(int i = 0;i < fNPtBins;i++){
    fHistResolutionCosThetaHE[i] -> Write();
    fHistResolutionPhiHE[i] -> Write();
    fHistResolutionCosThetaPhiHE[i] -> Write();
    fHistResolutionCosThetaCS[i] -> Write();
    fHistResolutionPhiCS[i] -> Write();
    fHistResolutionCosThetaPhiCS[i] -> Write();
  }
  fileResolution -> Close();

}
//______________________________________________________________________________
void AccxEffCalculator::ComputeAccxEff(string strSample, string nameOutputFile) {
  int nEvents = 0;
  int indexPt = 0;

  const int nVar = 5;
  int nBins[nVar] = {100,120,100,50,50};
  double minVar[nVar] = {0.,2.,-1.,0.,0.};
  double maxVar[nVar] = {10.,5.,1.,fPi,2*fPi};
  double varArray[nVar];

  const int nVarTest = 4;
  int nBinsTest[nVarTest] = {100,100,100,100};
  double minVarTest[nVarTest] = {0.,-1.,-fPi,0.};
  double maxVarTest[nVarTest] = {10.,1.,fPi,2*fPi};
  double varArrayTest[nVarTest];

  THnSparseD *histNVarHE = new THnSparseD("histNVarHE","histNVarHE",nVar,nBins,minVar,maxVar);
  THnSparseD *histNVarCS = new THnSparseD("histNVarCS","histNVarCS",nVar,nBins,minVar,maxVar);
  THnSparseD *histNVarTestGenHE = new THnSparseD("histNVarTestGenHE","histNVarTestGenHE",nVarTest,nBinsTest,minVarTest,maxVarTest);
  THnSparseD *histNVarTestGenCS = new THnSparseD("histNVarTestGenCS","histNVarTestGenCS",nVarTest,nBinsTest,minVarTest,maxVarTest);
  THnSparseD *histNVarTestRecHE = new THnSparseD("histNVarTestRecHE","histNVarTestRecHE",nVarTest,nBinsTest,minVarTest,maxVarTest);
  THnSparseD *histNVarTestRecCS = new THnSparseD("histNVarTestRecCS","histNVarTestRecCS",nVarTest,nBinsTest,minVarTest,maxVarTest);

  if(strSample == "FullStat"){nEvents = fTreeAccxEff -> GetEntries();}
  if(strSample == "TestStat"){nEvents = 5000000;}
  printf("N events = %i \n",nEvents);

  for(int i = 0;i < fNPtBins;i++){
    // HELICITY
    fHistGenCosThetaHE[i] = new TH1D(Form("histGenCosThetaHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistGenCosThetaHE[i] -> Sumw2();
    fHistRecCosThetaHE[i] = new TH1D(Form("histRecCosThetaHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistRecCosThetaHE[i] -> Sumw2();
    fHistGenPhiHE[i] = new TH1D(Form("histGenPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistGenPhiHE[i] -> Sumw2();
    fHistRecPhiHE[i] = new TH1D(Form("histRecPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistRecPhiHE[i] -> Sumw2();
    fHistGenPhiTildeHE[i] = new TH1D(Form("histGenPhiTildeHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]); fHistGenPhiTildeHE[i] -> Sumw2();
    fHistRecPhiTildeHE[i] = new TH1D(Form("histRecPhiTildeHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]); fHistRecPhiTildeHE[i] -> Sumw2();
    fHistGenCosThetaPhiHE[i] = new TH2D(Form("histGenCosThetaPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistGenCosThetaPhiHE[i] -> Sumw2();
    fHistRecCosThetaPhiHE[i] = new TH2D(Form("histRecCosThetaPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistRecCosThetaPhiHE[i] -> Sumw2();

    fHistGenPhiTildeNarrowHE[i] = new TH1D(Form("histGenPhiTildeNarrowHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",200,-2*fPi,2*fPi); fHistGenPhiTildeNarrowHE[i] -> Sumw2();
    fHistRecPhiTildeNarrowHE[i] = new TH1D(Form("histRecPhiTildeNarrowHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",200,-2*fPi,2*fPi); fHistRecPhiTildeNarrowHE[i] -> Sumw2();

    fHistRecPhiTildeNarrowLeftBellHE[i] = new TH1D(Form("histRecPhiTildeNarrowLeftBellHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",200,-2*fPi,2*fPi); fHistRecPhiTildeNarrowLeftBellHE[i] -> Sumw2();
    fHistRecPhiTildeNarrowRightBellHE[i] = new TH1D(Form("histRecPhiTildeNarrowRightBellHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",200,-2*fPi,2*fPi); fHistRecPhiTildeNarrowRightBellHE[i] -> Sumw2();
    fHistRecPhiTildeNarrowAllBellHE[i] = new TH1D(Form("histRecPhiTildeNarrowAllBellHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",200,-2*fPi,2*fPi); fHistRecPhiTildeNarrowAllBellHE[i] -> Sumw2();

    fHistGenCosThetaPhiNarrowHE[i] = new TH2D(Form("histGenCosThetaPhiNarrowHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",100,-1.,1.,100,-fPi,fPi); fHistGenCosThetaPhiNarrowHE[i] -> Sumw2();
    fHistRecCosThetaPhiNarrowHE[i] = new TH2D(Form("histRecCosThetaPhiNarrowHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",100,-1.,1.,100,-fPi,fPi); fHistRecCosThetaPhiNarrowHE[i] -> Sumw2();

    fHistGenPtHE[i] = new TH1D(Form("histGenPtHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",100,0.,10.); fHistGenPtHE[i] -> Sumw2();
    fHistRecPtHE[i] = new TH1D(Form("histRecPtHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",100,0.,10.); fHistRecPtHE[i] -> Sumw2();

    // COLLINS-SOPER
    fHistGenCosThetaCS[i] = new TH1D(Form("histGenCosThetaCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistGenCosThetaCS[i] -> Sumw2();
    fHistRecCosThetaCS[i] = new TH1D(Form("histRecCosThetaCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistRecCosThetaCS[i] -> Sumw2();
    fHistGenPhiCS[i] = new TH1D(Form("histGenPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistGenPhiCS[i] -> Sumw2();
    fHistRecPhiCS[i] = new TH1D(Form("histRecPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistRecPhiCS[i] -> Sumw2();
    fHistGenPhiTildeCS[i] = new TH1D(Form("histGenPhiTildeCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]); fHistGenPhiTildeCS[i] -> Sumw2();
    fHistRecPhiTildeCS[i] = new TH1D(Form("histRecPhiTildeCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]); fHistRecPhiTildeCS[i] -> Sumw2();
    fHistGenCosThetaPhiCS[i] = new TH2D(Form("histGenCosThetaPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistGenCosThetaPhiCS[i] -> Sumw2();
    fHistRecCosThetaPhiCS[i] = new TH2D(Form("histRecCosThetaPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistRecCosThetaPhiCS[i] -> Sumw2();

    fHistGenPhiTildeNarrowCS[i] = new TH1D(Form("histGenPhiTildeNarrowCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",200,-2*fPi,2*fPi); fHistGenPhiTildeNarrowCS[i] -> Sumw2();
    fHistRecPhiTildeNarrowCS[i] = new TH1D(Form("histRecPhiTildeNarrowCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",200,-2*fPi,2*fPi); fHistRecPhiTildeNarrowCS[i] -> Sumw2();

    fHistRecPhiTildeNarrowLeftBellCS[i] = new TH1D(Form("histRecPhiTildeNarrowLeftBellCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",200,-2*fPi,2*fPi); fHistRecPhiTildeNarrowLeftBellCS[i] -> Sumw2();
    fHistRecPhiTildeNarrowRightBellCS[i] = new TH1D(Form("histRecPhiTildeNarrowRightBellCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",200,-2*fPi,2*fPi); fHistRecPhiTildeNarrowRightBellCS[i] -> Sumw2();
    fHistRecPhiTildeNarrowAllBellCS[i] = new TH1D(Form("histRecPhiTildeNarrowAllBellCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",200,-2*fPi,2*fPi); fHistRecPhiTildeNarrowAllBellCS[i] -> Sumw2();

    fHistGenCosThetaPhiNarrowCS[i] = new TH2D(Form("histGenCosThetaPhiNarrowCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",100,-1.,1.,100,-fPi,fPi); fHistGenCosThetaPhiNarrowCS[i] -> Sumw2();
    fHistRecCosThetaPhiNarrowCS[i] = new TH2D(Form("histRecCosThetaPhiNarrowCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",100,-1.,1.,100,-fPi,fPi); fHistRecCosThetaPhiNarrowCS[i] -> Sumw2();

    fHistGenPtCS[i] = new TH1D(Form("histGenPtCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",100,0.,10.); fHistGenPtCS[i] -> Sumw2();
    fHistRecPtCS[i] = new TH1D(Form("histRecPtCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",100,0.,10.); fHistRecPtCS[i] -> Sumw2();
  }

  fHistGenCosThetaHEPt = new TH2D("histGenCosThetaHEPt","",100,-1,1,100,0,10); fHistGenCosThetaHEPt -> Sumw2();
  fHistRecCosThetaHEPt = new TH2D("histRecCosThetaHEPt","",100,-1,1,100,0,10); fHistRecCosThetaHEPt -> Sumw2();
  fHistGenPhiHEPt = new TH2D("histGenPhiHEPt","",100,0,fPi,100,0,10); fHistGenPhiHEPt -> Sumw2();
  fHistRecPhiHEPt = new TH2D("histRecPhiHEPt","",100,0,fPi,100,0,10); fHistRecPhiHEPt -> Sumw2();
  fHistGenPhiTildeHEPt = new TH2D("histGenPhiTildeHEPt","",100,0.,2*fPi,100,0,15); fHistGenPhiTildeHEPt -> Sumw2();
  fHistRecPhiTildeHEPt = new TH2D("histRecPhiTildeHEPt","",100,0.,2*fPi,100,0,15); fHistRecPhiTildeHEPt -> Sumw2();

  fHistGenCosThetaCSPt = new TH2D("histGenCosThetaCSPt","",100,-1,1,100,0,10); fHistGenCosThetaCSPt -> Sumw2();
  fHistRecCosThetaCSPt = new TH2D("histRecCosThetaCSPt","",100,-1,1,100,0,10); fHistRecCosThetaCSPt -> Sumw2();
  fHistGenPhiCSPt = new TH2D("histGenPhiCSPt","",100,0,fPi,100,0,10); fHistGenPhiCSPt -> Sumw2();
  fHistRecPhiCSPt = new TH2D("histRecPhiCSPt","",100,0,fPi,100,0,10); fHistRecPhiCSPt -> Sumw2();
  fHistGenPhiTildeCSPt = new TH2D("histGenPhiTildeCSPt","",100,0.,2*fPi,100,0,15); fHistGenPhiTildeCSPt -> Sumw2();
  fHistRecPhiTildeCSPt = new TH2D("histRecPhiTildeCSPt","",100,0.,2*fPi,100,0,15); fHistRecPhiTildeCSPt -> Sumw2();

  printf("- Configuring Fiducial Box (only for Reconstructed events) \n");
  printf("%f < CosTheta < %f \n",fCosThetaValues[1],fCosThetaValues[fNCosThetaBins-1]);
  printf("%f < |Phi| < %f \n",fPhiValues[1],fPhiValues[fNPhiBins-1]);

  double tmpVar = 0;

  for(int i = 0;i < nEvents;i++){
    printf("Reading : %3.2f %% \r",((double) i/(double) nEvents)*100.);
    fTreeAccxEff -> GetEntry(i);

    for(int j = 0;j < fNDimuGen;j++){
      if(fDimuYGen[j] > -4. && fDimuYGen[j] < -2.5){
        ////////////////////////////////////////////////////////////////////////
        // HELICITY
            while(fDimuPtGen[j] < fMinPt[indexPt] || fDimuPtGen[j] > fMaxPt[indexPt]){indexPt++;}
            // 1D approach
            fHistGenCosThetaHE[indexPt] -> Fill(fCosThetaHEGen[j]);
            fHistGenPhiHE[indexPt] -> Fill(TMath::Abs(fPhiHEGen[j]));
            fHistGenPtHE[indexPt] -> Fill(fDimuPtGen[j]);

            tmpVar = fPhiHEGen[j] + fPi;
            if(fCosThetaHEGen[j] < 0.){fPhiTilde = tmpVar - (3./4.)*fPi;}
            if(fCosThetaHEGen[j] > 0.){fPhiTilde = tmpVar - (1./4.)*fPi;}
            if(fPhiTilde > 2*fPi){fPhiTilde = fPhiTilde - 2*fPi;}
            if(fPhiTilde < 0.){fPhiTilde = 2*fPi + fPhiTilde;}
            fHistGenPhiTildeHE[indexPt] -> Fill(fPhiTilde);
            fHistGenPhiTildeNarrowHE[indexPt] -> Fill(fPhiTilde);

            // CosTheta, Phi, PhiTilde vs Pt
            fHistGenCosThetaHEPt -> Fill(fCosThetaHEGen[j],fDimuPtGen[j]);
            fHistGenPhiHEPt -> Fill(TMath::Abs(fPhiHEGen[j]),fDimuPtGen[j]);
            fHistGenPhiTildeHEPt -> Fill(fPhiTilde,fDimuPtGen[j]);

            varArrayTest[0] = fDimuPtGen[j];
            varArrayTest[1] = fCosThetaHEGen[j];
            varArrayTest[2] = fPhiHEGen[j];
            varArrayTest[3] = fPhiTilde;
            histNVarTestGenHE -> Fill(varArrayTest);

            // 2D approach
            fHistGenCosThetaPhiNarrowHE[indexPt] -> Fill(fCosThetaHEGen[j],fPhiHEGen[j]);
            fHistGenCosThetaPhiHE[indexPt] -> Fill(fCosThetaHEGen[j],TMath::Abs(fPhiHEGen[j]));
            indexPt = 0;

        // COLLINS-SOPER
            fHistGenCosThetaCSPt -> Fill(fCosThetaCSGen[j],fDimuPtGen[j]);
            fHistGenPhiCSPt -> Fill(TMath::Abs(fPhiCSGen[j]),fDimuPtGen[j]);
            while(fDimuPtGen[j] < fMinPt[indexPt] || fDimuPtGen[j] > fMaxPt[indexPt]){indexPt++;}
            // 1D approach
            fHistGenCosThetaCS[indexPt] -> Fill(fCosThetaCSGen[j]);
            fHistGenPhiCS[indexPt] -> Fill(TMath::Abs(fPhiCSGen[j]));
            fHistGenPtCS[indexPt] -> Fill(fDimuPtGen[j]);

            tmpVar = fPhiCSGen[j] + fPi;

            if(fCosThetaCSGen[j] < 0.){fPhiTilde = tmpVar - (3./4.)*fPi;}
            if(fCosThetaCSGen[j] > 0.){fPhiTilde = tmpVar - (1./4.)*fPi;}
            if(fPhiTilde > 2*fPi){fPhiTilde = fPhiTilde - 2*fPi;}
            if(fPhiTilde < 0.){fPhiTilde = 2*fPi + fPhiTilde;}
            fHistGenPhiTildeCS[indexPt] -> Fill(fPhiTilde);
            fHistGenPhiTildeNarrowCS[indexPt] -> Fill(fPhiTilde);

            // CosTheta, Phi, PhiTilde vs Pt
            fHistGenCosThetaCSPt -> Fill(fCosThetaCSGen[j],fDimuPtGen[j]);
            fHistGenPhiCSPt -> Fill(TMath::Abs(fPhiCSGen[j]),fDimuPtGen[j]);
            fHistGenPhiTildeCSPt -> Fill(fPhiTilde,fDimuPtGen[j]);

            varArrayTest[0] = fDimuPtGen[j];
            varArrayTest[1] = fCosThetaCSGen[j];
            varArrayTest[2] = fPhiCSGen[j];
            varArrayTest[3] = fPhiTilde;
            histNVarTestGenCS -> Fill(varArrayTest);

            // 2D approach
            fHistGenCosThetaPhiNarrowCS[indexPt] -> Fill(fCosThetaCSGen[j],fPhiCSGen[j]);
            fHistGenCosThetaPhiCS[indexPt] -> Fill(fCosThetaCSGen[j],TMath::Abs(fPhiCSGen[j]));
            indexPt = 0;
        ////////////////////////////////////////////////////////////////////////
      }
    }

    for(int j = 0;j < fNDimuRec;j++){
      if(fDimuYRec[j] > -4. && fDimuYRec[j] < -2.5){
        if(fDimuMatchRec[j] == 2){
          if(fDimuMassRec[j] > 2 && fDimuMassRec[j] < 5){
            ////////////////////////////////////////////////////////////////////
            // HELICITY
            // CosTheta, Phi vs Pt
            fHistRecCosThetaHEPt -> Fill(fCosThetaHERec[j],fDimuPtRec[j]);
            fHistRecPhiHEPt -> Fill(TMath::Abs(fPhiHERec[j]),fDimuPtRec[j]);

            if(TMath::Abs(fPhiHERec[j]) > fPhiValues[1] && TMath::Abs(fPhiHERec[j]) < fPhiValues[fNPhiBins-1]){
              if(fCosThetaHERec[j] > fCosThetaValues[1] && fCosThetaHERec[j] < fCosThetaValues[fNCosThetaBins-1]){

                while(fDimuPtRec[j] < fMinPt[indexPt] || fDimuPtRec[j] > fMaxPt[indexPt]){indexPt++;}
                // 1D approach
                fHistRecCosThetaHE[indexPt] -> Fill(fCosThetaHERec[j]);
                fHistRecPhiHE[indexPt] -> Fill(TMath::Abs(fPhiHERec[j]));
                fHistRecPtHE[indexPt] -> Fill(fDimuPtRec[j]);

                tmpVar = fPhiHERec[j] + fPi;

                if(fCosThetaHERec[j] < 0.){fPhiTilde = tmpVar - (3./4.)*fPi;}
                if(fCosThetaHERec[j] > 0.){fPhiTilde = tmpVar - (1./4.)*fPi;}
                if(fPhiTilde > 2*fPi){fPhiTilde = fPhiTilde - 2*fPi;}
                if(fPhiTilde < 0.){fPhiTilde = 2*fPi + fPhiTilde;}
                fHistRecPhiTildeHE[indexPt] -> Fill(fPhiTilde);
                fHistRecPhiTildeNarrowHE[indexPt] -> Fill(fPhiTilde);
                if(fCosThetaHERec[j] < 0.){fHistRecPhiTildeNarrowLeftBellHE[indexPt] -> Fill(fPhiHERec[j]);}
                if(fCosThetaHERec[j] > 0.){fHistRecPhiTildeNarrowRightBellHE[indexPt] -> Fill(fPhiHERec[j]);}
                fHistRecPhiTildeNarrowAllBellHE[indexPt] -> Fill(fPhiHERec[j]);

                // PhiTilde vs Pt
                fHistRecPhiTildeHEPt -> Fill(fPhiTilde,fDimuPtRec[j]);

                varArrayTest[0] = fDimuPtRec[j];
                varArrayTest[1] = fCosThetaHERec[j];
                varArrayTest[2] = fPhiHERec[j];
                varArrayTest[3] = fPhiTilde;
                histNVarTestRecHE -> Fill(varArrayTest);

                ////////////////////////////////////////////////////////////////
                //Filling the THnSparse
                varArray[0] = fDimuPtRec[j];
                varArray[1] = fDimuMassRec[j];
                varArray[2] = fCosThetaHERec[j];
                varArray[3] = TMath::Abs(fPhiHERec[j]);
                varArray[4] = fPhiTilde;
                histNVarHE -> Fill(varArray);
                ////////////////////////////////////////////////////////////////

                // 2D approach
                fHistRecCosThetaPhiNarrowHE[indexPt] -> Fill(fCosThetaHERec[j],fPhiHERec[j]);
                fHistRecCosThetaPhiHE[indexPt] -> Fill(fCosThetaHERec[j],TMath::Abs(fPhiHERec[j]));
                indexPt = 0;
              }
            }
            // COLLINS-SOPER
            // CosTheta, Phi vs Pt
            fHistRecCosThetaCSPt -> Fill(fCosThetaCSRec[j],fDimuPtRec[j]);
            fHistRecPhiCSPt -> Fill(TMath::Abs(fPhiCSRec[j]),fDimuPtRec[j]);
            if(TMath::Abs(fPhiCSRec[j]) > fPhiValues[1] && TMath::Abs(fPhiCSRec[j]) < fPhiValues[fNPhiBins-1]){
              if(fCosThetaCSRec[j] > fCosThetaValues[1] && fCosThetaCSRec[j] < fCosThetaValues[fNCosThetaBins-1]){

                while(fDimuPtRec[j] < fMinPt[indexPt] || fDimuPtRec[j] > fMaxPt[indexPt]){indexPt++;}
                // 1D approach
                fHistRecCosThetaCS[indexPt] -> Fill(fCosThetaCSRec[j]);
                fHistRecPhiCS[indexPt] -> Fill(TMath::Abs(fPhiCSRec[j]));
                fHistRecPtCS[indexPt] -> Fill(fDimuPtRec[j]);

                tmpVar = fPhiCSRec[j] + fPi;

                if(fCosThetaCSRec[j] < 0.){fPhiTilde = tmpVar - (3./4.)*fPi;}
                if(fCosThetaCSRec[j] > 0.){fPhiTilde = tmpVar - (1./4.)*fPi;}
                if(fPhiTilde > 2*fPi){fPhiTilde = fPhiTilde - 2*fPi;}
                if(fPhiTilde < 0.){fPhiTilde = 2*fPi + fPhiTilde;}
                fHistRecPhiTildeCS[indexPt] -> Fill(fPhiTilde);
                fHistRecPhiTildeNarrowCS[indexPt] -> Fill(fPhiTilde);
                if(fCosThetaCSRec[j] < 0.){fHistRecPhiTildeNarrowLeftBellCS[indexPt] -> Fill(fPhiCSRec[j]);}
                if(fCosThetaCSRec[j] > 0.){fHistRecPhiTildeNarrowRightBellCS[indexPt] -> Fill(fPhiCSRec[j]);}
                fHistRecPhiTildeNarrowAllBellCS[indexPt] -> Fill(fPhiCSRec[j]);

                // PhiTilde vs Pt
                fHistRecPhiTildeCSPt -> Fill(fPhiTilde,fDimuPtRec[j]);

                varArrayTest[0] = fDimuPtRec[j];
                varArrayTest[1] = fCosThetaCSRec[j];
                varArrayTest[2] = fPhiCSRec[j];
                varArrayTest[3] = fPhiTilde;
                histNVarTestRecCS -> Fill(varArrayTest);

                ////////////////////////////////////////////////////////////////
                //Filling the THnSparse
                varArray[0] = fDimuPtRec[j];
                varArray[1] = fDimuMassRec[j];
                varArray[2] = fCosThetaCSRec[j];
                varArray[3] = TMath::Abs(fPhiCSRec[j]);
                varArray[4] = fPhiTilde;
                histNVarCS -> Fill(varArray);
                ////////////////////////////////////////////////////////////////

                // 2D approach
                fHistRecCosThetaPhiNarrowCS[indexPt] -> Fill(fCosThetaCSRec[j],fPhiCSRec[j]);
                fHistRecCosThetaPhiCS[indexPt] -> Fill(fCosThetaCSRec[j],TMath::Abs(fPhiCSRec[j]));
                indexPt = 0;
              }
            }
            ////////////////////////////////////////////////////////////////////
          }
        }
      }
    }

  }
  printf("\n");

  for(int i = 0;i < fNPtBins;i++){
    // HELICITY
    fHistAccxEffCosThetaHE[i] = new TH1D(Form("histAccxEffCosThetaHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]);
    fHistAccxEffCosThetaHE[i] -> Divide(fHistRecCosThetaHE[i],fHistGenCosThetaHE[i],1,1,"B");

    fHistAccxEffPhiHE[i] = new TH1D(Form("histAccxEffPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistAccxEffPhiHE[i] -> Divide(fHistRecPhiHE[i],fHistGenPhiHE[i],1,1,"B");

    fHistAccxEffPhiTildeHE[i] = new TH1D(Form("histAccxEffPhiTildeHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]);
    fHistAccxEffPhiTildeHE[i] -> Divide(fHistRecPhiTildeHE[i],fHistGenPhiTildeHE[i],1,1,"B");

    fHistAccxEffPhiTildeNarrowHE[i] = new TH1D(Form("histAccxEffPhiTildeNarrowHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",200,-2*fPi,2*fPi); //
    fHistAccxEffPhiTildeNarrowHE[i] -> Divide(fHistRecPhiTildeNarrowHE[i],fHistGenPhiTildeNarrowHE[i],1,1,"B");                     //

    fHistAccxEffCosThetaPhiHE[i] = new TH2D(Form("histAccxEffCosThetaPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]);
    fHistAccxEffCosThetaPhiHE[i] -> Divide(fHistRecCosThetaPhiHE[i],fHistGenCosThetaPhiHE[i],1,1,"B");

    fHistAccxEffCosThetaPhiStatRelHE[i] = new TH2D(Form("histAccxEffCosThetaPhiStatRelHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]);

    fHistAccxEffCosThetaPhiNarrowHE[i] = new TH2D(Form("histAccxEffCosThetaPhiNarrowHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",100,-1.,1.,100,-fPi,fPi);
    fHistAccxEffCosThetaPhiNarrowHE[i] -> Divide(fHistRecCosThetaPhiNarrowHE[i],fHistGenCosThetaPhiNarrowHE[i],1,1,"B");

    fHistAccxEffPtHE[i] = new TH1D(Form("histAccxEffPtHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",100,0.,10.);
    fHistAccxEffPtHE[i] -> Divide(fHistRecPtHE[i],fHistGenPtHE[i],1,1,"B");

    // COLLINS-SOPER
    fHistAccxEffCosThetaCS[i] = new TH1D(Form("histAccxEffCosThetaCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]);
    fHistAccxEffCosThetaCS[i] -> Divide(fHistRecCosThetaCS[i],fHistGenCosThetaCS[i],1,1,"B");

    fHistAccxEffPhiCS[i] = new TH1D(Form("histAccxEffPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistAccxEffPhiCS[i] -> Divide(fHistRecPhiCS[i],fHistGenPhiCS[i],1,1,"B");

    fHistAccxEffPhiTildeCS[i] = new TH1D(Form("histAccxEffPhiTildeCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]);
    fHistAccxEffPhiTildeCS[i] -> Divide(fHistRecPhiTildeCS[i],fHistGenPhiTildeCS[i],1,1,"B");

    fHistAccxEffPhiTildeNarrowCS[i] = new TH1D(Form("histAccxEffPhiTildeNarrowCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",200,-2*fPi,2*fPi); //
    fHistAccxEffPhiTildeNarrowCS[i] -> Divide(fHistRecPhiTildeNarrowCS[i],fHistGenPhiTildeNarrowCS[i],1,1,"B");                     //

    fHistAccxEffCosThetaPhiCS[i] = new TH2D(Form("histAccxEffCosThetaPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]);
    fHistAccxEffCosThetaPhiCS[i] -> Divide(fHistRecCosThetaPhiCS[i],fHistGenCosThetaPhiCS[i],1,1,"B");

    fHistAccxEffCosThetaPhiStatRelCS[i] = new TH2D(Form("histAccxEffCosThetaPhiStatRelCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]);

    fHistAccxEffCosThetaPhiNarrowCS[i] = new TH2D(Form("histAccxEffCosThetaPhiNarrowCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",100,-1.,1.,100,-fPi,fPi);
    fHistAccxEffCosThetaPhiNarrowCS[i] -> Divide(fHistRecCosThetaPhiNarrowCS[i],fHistGenCosThetaPhiNarrowCS[i],1,1,"B");

    fHistAccxEffPtCS[i] = new TH1D(Form("histAccxEffPtCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",100,0.,10.);
    fHistAccxEffPtCS[i] -> Divide(fHistRecPtCS[i],fHistGenPtCS[i],1,1,"B");
  }

  fHistAccxEffCosThetaHEPt = new TH2D("histAccxEffCosThetaHEPt","",100,-1,1,100,0,10);
  fHistAccxEffCosThetaHEPt -> Divide(fHistRecCosThetaHEPt,fHistGenCosThetaHEPt,1,1,"B");
  fHistAccxEffPhiHEPt = new TH2D("histAccxEffPhiHEPt","",100,0,fPi,100,0,10);
  fHistAccxEffPhiHEPt -> Divide(fHistRecPhiHEPt,fHistGenPhiHEPt,1,1,"B");
  fHistAccxEffPhiTildeHEPt = new TH2D("histAccxEffPhiTildeHEPt","",100,0.,2*fPi,100,0,15);
  fHistAccxEffPhiTildeHEPt -> Divide(fHistRecPhiTildeHEPt,fHistGenPhiTildeHEPt,1,1,"B");

  fHistAccxEffCosThetaCSPt = new TH2D("histAccxEffCosThetaCSPt","",100,-1,1,100,0,10);
  fHistAccxEffCosThetaCSPt -> Divide(fHistRecCosThetaCSPt,fHistGenCosThetaCSPt,1,1,"B");
  fHistAccxEffPhiCSPt = new TH2D("histAccxEffPhiCSPt","",100,0,fPi,100,0,10);
  fHistAccxEffPhiCSPt -> Divide(fHistRecPhiCSPt,fHistGenPhiCSPt,1,1,"B");
  fHistAccxEffPhiTildeCSPt = new TH2D("histAccxEffPhiTildeCSPt","",100,0.,2*fPi,100,0,15);
  fHistAccxEffPhiTildeCSPt -> Divide(fHistRecPhiTildeCSPt,fHistGenPhiTildeCSPt,1,1,"B");

  for(int i = 0;i < fNPtBins;i++){
    for(int j = 0;j < fNCosThetaBins;j++){
      for(int k = 0;k < fNPhiBins;k++){
        fHistAccxEffCosThetaPhiStatRelHE[i] -> SetBinContent(j+1,k+1,(fHistAccxEffCosThetaPhiHE[i] -> GetBinError(j+1,k+1)/fHistAccxEffCosThetaPhiHE[i] -> GetBinContent(j+1,k+1))*100);
        fHistAccxEffCosThetaPhiStatRelCS[i] -> SetBinContent(j+1,k+1,(fHistAccxEffCosThetaPhiCS[i] -> GetBinError(j+1,k+1)/fHistAccxEffCosThetaPhiCS[i] -> GetBinContent(j+1,k+1))*100);
      }
    }
  }

  TFile *fileAccxEff = new TFile(nameOutputFile.c_str(),"RECREATE");
  histNVarHE -> Write();
  histNVarCS -> Write();
  histNVarTestGenHE -> Write();
  histNVarTestGenCS -> Write();
  histNVarTestRecHE -> Write();
  histNVarTestRecCS -> Write();
  for(int i = 0;i < fNPtBins;i++){
    fHistGenCosThetaHE[i] -> Write(); fHistRecCosThetaHE[i] -> Write(); fHistAccxEffCosThetaHE[i] -> Write();
    fHistGenPhiHE[i] -> Write(); fHistRecPhiHE[i] -> Write(); fHistAccxEffPhiHE[i] -> Write();
    fHistGenPhiTildeHE[i] -> Write(); fHistRecPhiTildeHE[i] -> Write(); fHistAccxEffPhiTildeHE[i] -> Write();
    fHistGenCosThetaPhiHE[i] -> Write(); fHistRecCosThetaPhiHE[i] -> Write(); fHistAccxEffCosThetaPhiHE[i] -> Write();
    fHistAccxEffCosThetaPhiStatRelHE[i] -> Write();

    fHistGenCosThetaCS[i] -> Write(); fHistRecCosThetaCS[i] -> Write(); fHistAccxEffCosThetaCS[i] -> Write();
    fHistGenPhiCS[i] -> Write(); fHistRecPhiCS[i] -> Write(); fHistAccxEffPhiCS[i] -> Write();
    fHistGenPhiTildeCS[i] -> Write(); fHistRecPhiTildeCS[i] -> Write(); fHistAccxEffPhiTildeCS[i] -> Write();
    fHistGenCosThetaPhiCS[i] -> Write(); fHistRecCosThetaPhiCS[i] -> Write(); fHistAccxEffCosThetaPhiCS[i] -> Write();
    fHistAccxEffCosThetaPhiStatRelCS[i] -> Write();

    fHistGenPhiTildeNarrowHE[i] -> Write();
    fHistRecPhiTildeNarrowHE[i] -> Write();
    fHistAccxEffPhiTildeNarrowHE[i] -> Write();

    fHistGenPhiTildeNarrowCS[i] -> Write();
    fHistRecPhiTildeNarrowCS[i] -> Write();
    fHistAccxEffPhiTildeNarrowCS[i] -> Write();

    fHistRecPhiTildeNarrowLeftBellHE[i] -> Write();
    fHistRecPhiTildeNarrowRightBellHE[i] -> Write();
    fHistRecPhiTildeNarrowAllBellHE[i] -> Write();
    fHistRecPhiTildeNarrowLeftBellCS[i] -> Write();
    fHistRecPhiTildeNarrowRightBellCS[i] -> Write();
    fHistRecPhiTildeNarrowAllBellCS[i] -> Write();

    fHistAccxEffCosThetaPhiNarrowHE[i] -> Write();
    fHistAccxEffCosThetaPhiNarrowCS[i] -> Write();

    fHistAccxEffPtHE[i] -> Write();
    fHistAccxEffPtCS[i] -> Write();
  }
  fHistGenCosThetaHEPt -> Write(); fHistGenPhiHEPt -> Write(); fHistAccxEffCosThetaHEPt -> Write();
  fHistRecCosThetaHEPt -> Write(); fHistRecPhiHEPt -> Write(); fHistAccxEffPhiHEPt -> Write();

  fHistGenCosThetaCSPt -> Write(); fHistGenPhiCSPt -> Write(); fHistAccxEffCosThetaCSPt -> Write();
  fHistRecCosThetaCSPt -> Write(); fHistRecPhiCSPt -> Write(); fHistAccxEffPhiCSPt -> Write();

  fHistGenPhiTildeHEPt -> Write(); fHistRecPhiTildeHEPt -> Write(); fHistAccxEffPhiTildeHEPt -> Write();
  fHistGenPhiTildeCSPt -> Write(); fHistRecPhiTildeCSPt -> Write(); fHistAccxEffPhiTildeCSPt -> Write();

  fileAccxEff -> Close();
}
//______________________________________________________________________________
//void AccxEffCalculator::ReWeightAccxEff(string refFrame, Double_t LambdaTheta, Double_t LambdaPhi, string strSample, Bool_t saveFile, string nameOutputFile) {
//void AccxEffCalculator::ReWeightAccxEff(Double_t polParHE[], Double_t polParCS[], string strSample, Bool_t saveFile, string nameOutputFile) {
void AccxEffCalculator::ReWeightAccxEff(Double_t polParHE[3][4], Double_t polParCS[3][4], string strSample, Bool_t saveFile, string nameOutputFile) {
  int indexPt = 0;
  int nEvents = 0;

  TF1 *funcCosThetaHE[4], *funcPhiHE[4], *funcPhiTildeHE[4];
  TF1 *funcCosThetaCS[4], *funcPhiCS[4], *funcPhiTildeCS[4];

  TF1 *funcCosThetaPhiHE[4], *funcCosThetaPhiCS[4];

  for(int i = 0;i < 4;i++){
    funcCosThetaHE[i] = new TF1(Form("funcCosThetaHE_%i",i),"(1/(3 + [0]))*(1 + [0]*x*x)",-1,1); funcCosThetaHE[i] -> SetParameter(0,polParHE[0][i]);
    funcPhiHE[i] = new TF1(Form("funcPhiHE_%i",i),"(1 + ((2*[1])/(3 + [0]))*cos(2*x))",0,fPi); funcPhiHE[i] -> SetParameter(0,polParHE[0][i]); funcPhiHE[i] -> SetParameter(1,polParHE[1][i]);
    funcPhiTildeHE[i] = new TF1(Form("funcPhiTildeHE_%i",i),"(1 + ((sqrt(2)*[2])/(3 + [1]))*cos(x))",0.,2*fPi); funcPhiTildeHE[i] -> SetParameter(0,polParHE[0][i]); funcPhiTildeHE[i] -> SetParameter(2,polParHE[2][i]);

    funcCosThetaPhiHE[i] = new TF2(Form("funcCosThetaPhiHE_%i",i),MyFuncPol,-1,1,0,fPi,4);
    funcCosThetaPhiHE[i] -> SetParameter(0,1.);
    funcCosThetaPhiHE[i] -> SetParameter(1,polParHE[0][i]);
    funcCosThetaPhiHE[i] -> SetParameter(2,polParHE[1][i]);
    funcCosThetaPhiHE[i] -> SetParameter(3,polParHE[2][i]);


    funcCosThetaCS[i] = new TF1(Form("funcCosThetaCS_%i",i),"(1/(3 + [0]))*(1 + [0]*x*x)",-1,1); funcCosThetaCS[i] -> SetParameter(0,polParCS[0][i]);
    funcPhiCS[i] = new TF1(Form("funcPhiCS_%i",i),"(1 + ((2*[1])/(3 + [0]))*cos(2*x))",0,fPi); funcPhiCS[i] -> SetParameter(0,polParCS[0][i]); funcPhiCS[i] -> SetParameter(1,polParCS[1][i]);
    funcPhiTildeCS[i] = new TF1(Form("funcPhiTildeCS_%i",i),"(1 + ((sqrt(2)*[2])/(3 + [1]))*cos(x))",0.,2*fPi); funcPhiTildeCS[i] -> SetParameter(0,polParCS[0][i]); funcPhiTildeCS[i] -> SetParameter(2,polParCS[2][i]);

    funcCosThetaPhiCS[i] = new TF2(Form("funcCosThetaPhiCS_%i",i),MyFuncPol,-1,1,0,fPi,4);
    funcCosThetaPhiCS[i] -> SetParameter(0,1.);
    funcCosThetaPhiCS[i] -> SetParameter(1,polParCS[0][i]);
    funcCosThetaPhiCS[i] -> SetParameter(2,polParCS[1][i]);
    funcCosThetaPhiCS[i] -> SetParameter(3,polParCS[2][i]);
  }

  /*TF1 *funcCosThetaHE = new TF1("funcCosThetaHE","(1/(3 + [0]))*(1 + [0]*x*x)",-1,1);
  funcCosThetaHE -> SetParameter(0,polParHE[0]);

  TF1 *funcPhiHE = new TF1("funcPhiHE","(1 + ((2*[1])/(3 + [0]))*cos(2*x))",0,fPi);
  funcPhiHE -> SetParameter(0,polParHE[0]);
  funcPhiHE -> SetParameter(1,polParHE[1]);

  TF1 *funcPhiTildeHE = new TF1("funcPhiTildeHE","(1 + ((sqrt(2)*[2])/(3 + [1]))*cos(x))",0.,2*fPi);
  funcPhiTildeHE -> SetParameter(0,polParHE[0]);
  funcPhiTildeHE -> SetParameter(1,polParHE[2]);

  TF1 *funcCosThetaCS = new TF1("funcCosThetaCS","(1/(3 + [0]))*(1 + [0]*x*x)",-1,1);
  funcCosThetaCS -> SetParameter(0,polParCS[0]);

  TF1 *funcPhiCS = new TF1("funcPhiCS","(1 + ((2*[1])/(3 + [0]))*cos(2*x))",0,fPi);
  funcPhiCS -> SetParameter(0,polParCS[0]);
  funcPhiCS -> SetParameter(1,polParCS[1]);

  TF1 *funcPhiTildeCS = new TF1("funcPhiTildeCS","(1 + ((sqrt(2)*[2])/(3 + [1]))*cos(x))",0.,2*fPi);
  funcPhiTildeCS -> SetParameter(0,polParCS[0]);
  funcPhiTildeCS -> SetParameter(1,polParCS[2]);*/

  if(strSample == "FullStat"){nEvents = fTreeAccxEff -> GetEntries();}
  if(strSample == "TestStat"){nEvents = 300000;}
  printf("N events = %i \n",nEvents);

  double CosThetaHEGen = 0;
  double CosThetaHERec = 0;
  double PhiHEGen = 0;
  double PhiHERec = 0;
  double PhiTildeHEGen = 0;
  double PhiTildeHERec = 0;

  double weightCosThetaHE = 0;
  double weightPhiHE = 0;
  double weightPhiTildeHE = 0;

  double weightCosThetaPhiHE = 0;

  double CosThetaCSGen = 0;
  double CosThetaCSRec = 0;
  double PhiCSGen = 0;
  double PhiCSRec = 0;
  double PhiTildeCSGen = 0;
  double PhiTildeCSRec = 0;

  double weightCosThetaCS = 0;
  double weightPhiCS = 0;
  double weightPhiTildeCS = 0;

  double weightCosThetaPhiCS = 0;

  for(int i = 0;i < fNPtBins;i++){
    /*fHistGenCosThetaReWeighted[i] = new TH1D(Form("histGenCosThetaReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistGenCosThetaReWeighted[i] -> Sumw2();
    fHistRecCosThetaReWeighted[i] = new TH1D(Form("histRecCosThetaReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistRecCosThetaReWeighted[i] -> Sumw2();
    fHistGenPhiReWeighted[i] = new TH1D(Form("histGenPhiReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistGenPhiReWeighted[i] -> Sumw2();
    fHistRecPhiReWeighted[i] = new TH1D(Form("histRecPhiReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistRecPhiReWeighted[i] -> Sumw2();
    fHistGenPhiTildeReWeighted[i] = new TH1D(Form("histGenPhiTildeReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]); fHistGenPhiTildeReWeighted[i] -> Sumw2();
    fHistRecPhiTildeReWeighted[i] = new TH1D(Form("histRecPhiTildeReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]); fHistRecPhiTildeReWeighted[i] -> Sumw2();*/

    fHistGenCosThetaHEReWeighted[i] = new TH1D(Form("histGenCosThetaHEReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistGenCosThetaHEReWeighted[i] -> Sumw2();
    fHistRecCosThetaHEReWeighted[i] = new TH1D(Form("histRecCosThetaHEReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistRecCosThetaHEReWeighted[i] -> Sumw2();
    fHistGenPhiHEReWeighted[i] = new TH1D(Form("histGenPhiHEReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistGenPhiHEReWeighted[i] -> Sumw2();
    fHistRecPhiHEReWeighted[i] = new TH1D(Form("histRecPhiHEReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistRecPhiHEReWeighted[i] -> Sumw2();
    fHistGenPhiTildeHEReWeighted[i] = new TH1D(Form("histGenPhiTildeHEReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]); fHistGenPhiTildeHEReWeighted[i] -> Sumw2();
    fHistRecPhiTildeHEReWeighted[i] = new TH1D(Form("histRecPhiTildeHEReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]); fHistRecPhiTildeHEReWeighted[i] -> Sumw2();

    fHistGenCosThetaCSReWeighted[i] = new TH1D(Form("histGenCosThetaCSReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistGenCosThetaCSReWeighted[i] -> Sumw2();
    fHistRecCosThetaCSReWeighted[i] = new TH1D(Form("histRecCosThetaCSReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistRecCosThetaCSReWeighted[i] -> Sumw2();
    fHistGenPhiCSReWeighted[i] = new TH1D(Form("histGenPhiCSReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistGenPhiCSReWeighted[i] -> Sumw2();
    fHistRecPhiCSReWeighted[i] = new TH1D(Form("histRecPhiCSReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistRecPhiCSReWeighted[i] -> Sumw2();
    fHistGenPhiTildeCSReWeighted[i] = new TH1D(Form("histGenPhiTildeCSReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]); fHistGenPhiTildeCSReWeighted[i] -> Sumw2();
    fHistRecPhiTildeCSReWeighted[i] = new TH1D(Form("histRecPhiTildeCSReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]); fHistRecPhiTildeCSReWeighted[i] -> Sumw2();

    fHistGenCosThetaPhiHEReWeighted[i] = new TH2D(Form("histGenCosThetaPhiHEReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistGenCosThetaPhiHEReWeighted[i] -> Sumw2();
    fHistRecCosThetaPhiHEReWeighted[i] = new TH2D(Form("histRecCosThetaPhiHEReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistRecCosThetaPhiHEReWeighted[i] -> Sumw2();
    fHistGenCosThetaPhiCSReWeighted[i] = new TH2D(Form("histGenCosThetaPhiCSReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistGenCosThetaPhiCSReWeighted[i] -> Sumw2();
    fHistRecCosThetaPhiCSReWeighted[i] = new TH2D(Form("histRecCosThetaPhiCSReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistRecCosThetaPhiCSReWeighted[i] -> Sumw2();
  }

  double tmpVar = 0;

  for(int i = 0;i < nEvents;i++){
    printf("Reading : %3.2f %% \r",((double) i/(double) nEvents)*100.);
    fTreeAccxEff -> GetEntry(i);

    for(int j = 0;j < fNDimuGen;j++){
      if(fDimuPtGen[j] < 2.){continue;}
      if(fDimuYGen[j] > -4. && fDimuYGen[j] < -2.5){
        ////////////////////////////////////////////////////////////////////////
        // HELICITY
        while(fDimuPtGen[j] < fMinPt[indexPt] || fDimuPtGen[j] > fMaxPt[indexPt]){indexPt++;}
        if(indexPt >= 4){indexPt = 0; continue;}

        tmpVar = fPhiHEGen[j] + fPi;
        if(fCosThetaHEGen[j] < 0.){fPhiTilde = tmpVar - (3./4.)*fPi;}
        if(fCosThetaHEGen[j] > 0.){fPhiTilde = tmpVar - (1./4.)*fPi;}
        if(fPhiTilde > 2*fPi){fPhiTilde = fPhiTilde - 2*fPi;}
        if(fPhiTilde < 0.){fPhiTilde = 2*fPi + fPhiTilde;}

        weightCosThetaHE = (funcCosThetaHE[indexPt] -> Eval(fCosThetaHEGen[j]))/(funcCosThetaHE[indexPt] -> GetMaximum()); CosThetaHEGen = fCosThetaHEGen[j];
        weightPhiHE = (funcPhiHE[indexPt] -> Eval(fPhiHEGen[j]))/(funcPhiHE[indexPt] -> GetMaximum()); PhiHEGen = fPhiHEGen[j];
        weightPhiTildeHE = (funcPhiTildeHE[indexPt] -> Eval(fPhiTilde))/(funcPhiTildeHE[indexPt] -> GetMaximum()); PhiTildeHEGen = fPhiTilde;

        weightCosThetaPhiHE = (funcCosThetaPhiHE[indexPt] -> Eval(fCosThetaHEGen[j],fPhiHEGen[j]))/(funcCosThetaPhiHE[indexPt] -> GetMaximum());

        fHistGenCosThetaHEReWeighted[indexPt] -> Fill(CosThetaHEGen,weightCosThetaHE);
        fHistGenPhiHEReWeighted[indexPt] -> Fill(TMath::Abs(PhiHEGen),weightPhiHE);
        fHistGenPhiTildeHEReWeighted[indexPt] -> Fill(PhiTildeHEGen,weightPhiTildeHE);
        fHistGenCosThetaPhiHEReWeighted[indexPt] -> Fill(CosThetaHEGen,TMath::Abs(PhiHEGen),weightCosThetaPhiHE);;
        indexPt = 0;
        // COLLINS-SOPER
        while(fDimuPtGen[j] < fMinPt[indexPt] || fDimuPtGen[j] > fMaxPt[indexPt]){indexPt++;}
        if(indexPt >= 4){indexPt = 0; continue;}

        tmpVar = fPhiCSGen[j] + fPi;
        if(fCosThetaCSGen[j] < 0.){fPhiTilde = tmpVar - (3./4.)*fPi;}
        if(fCosThetaCSGen[j] > 0.){fPhiTilde = tmpVar - (1./4.)*fPi;}
        if(fPhiTilde > 2*fPi){fPhiTilde = fPhiTilde - 2*fPi;}
        if(fPhiTilde < 0.){fPhiTilde = 2*fPi + fPhiTilde;}

        weightCosThetaCS = (funcCosThetaCS[indexPt] -> Eval(fCosThetaCSGen[j]))/(funcCosThetaCS[indexPt] -> GetMaximum()); CosThetaCSGen = fCosThetaCSGen[j];
        weightPhiCS = (funcPhiCS[indexPt] -> Eval(fPhiCSGen[j]))/(funcPhiCS[indexPt] -> GetMaximum()); PhiCSGen = fPhiCSGen[j];
        weightPhiTildeCS = (funcPhiTildeCS[indexPt] -> Eval(fPhiTilde))/(funcPhiTildeCS[indexPt] -> GetMaximum()); PhiTildeCSGen = fPhiTilde;

        weightCosThetaPhiCS = (funcCosThetaPhiCS[indexPt] -> Eval(fCosThetaCSGen[j],fPhiCSGen[j]))/(funcCosThetaPhiCS[indexPt] -> GetMaximum());

        fHistGenCosThetaCSReWeighted[indexPt] -> Fill(CosThetaCSGen,weightCosThetaCS);
        fHistGenPhiCSReWeighted[indexPt] -> Fill(TMath::Abs(PhiCSGen),weightPhiCS);
        fHistGenPhiTildeCSReWeighted[indexPt] -> Fill(PhiTildeCSGen,weightPhiTildeCS);
        fHistGenCosThetaPhiCSReWeighted[indexPt] -> Fill(CosThetaCSGen,TMath::Abs(PhiCSGen),weightCosThetaPhiCS);
        indexPt = 0;
        ////////////////////////////////////////////////////////////////////////
      }
    }

    for(int j = 0;j < fNDimuRec;j++){
      if(fDimuPtRec[j] < 2.){continue;}
      if(fDimuYRec[j] > -4. && fDimuYRec[j] < -2.5){
        if(fDimuMatchRec[j] == 2){
          if(fDimuMassRec[j] > 2 && fDimuMassRec[j] < 5){
            ////////////////////////////////////////////////////////////////////
            // HELICITY
            if(TMath::Abs(fPhiHERec[j]) > fPhiValues[1] && TMath::Abs(fPhiHERec[j]) < fPhiValues[fNPhiBins-1]){
              if(fCosThetaHERec[j] > fCosThetaValues[1] && fCosThetaHERec[j] < fCosThetaValues[fNCosThetaBins-1]){
                while(fDimuPtRec[j] < fMinPt[indexPt] || fDimuPtRec[j] > fMaxPt[indexPt]){indexPt++;}
                if(indexPt >= 4){indexPt = 0; continue;}

                tmpVar = fPhiHERec[j] + fPi;
                if(fCosThetaHERec[j] < 0.){fPhiTilde = tmpVar - (3./4.)*fPi;}
                if(fCosThetaHERec[j] > 0.){fPhiTilde = tmpVar - (1./4.)*fPi;}
                if(fPhiTilde > 2*fPi){fPhiTilde = fPhiTilde - 2*fPi;}
                if(fPhiTilde < 0.){fPhiTilde = 2*fPi + fPhiTilde;}

                fHistRecCosThetaHEReWeighted[indexPt] -> Fill(fCosThetaHERec[j],weightCosThetaHE);
                fHistRecPhiHEReWeighted[indexPt] -> Fill(TMath::Abs(fPhiHERec[j]),weightPhiHE);
                fHistRecPhiTildeHEReWeighted[indexPt] -> Fill(fPhiTilde,weightPhiTildeHE);
                fHistRecCosThetaPhiHEReWeighted[indexPt] -> Fill(fCosThetaHERec[j],TMath::Abs(fPhiHERec[j]),weightCosThetaPhiHE);
                indexPt = 0;
              }
            }
            // COLLINS-SOPER
            if(TMath::Abs(fPhiCSRec[j]) > fPhiValues[1] && TMath::Abs(fPhiCSRec[j]) < fPhiValues[fNPhiBins-1]){
              if(fCosThetaCSRec[j] > fCosThetaValues[1] && fCosThetaCSRec[j] < fCosThetaValues[fNCosThetaBins-1]){
                while(fDimuPtRec[j] < fMinPt[indexPt] || fDimuPtRec[j] > fMaxPt[indexPt]){indexPt++;}
                if(indexPt >= 4){indexPt = 0; continue;}

                tmpVar = fPhiCSRec[j] + fPi;
                if(fCosThetaCSRec[j] < 0.){fPhiTilde = tmpVar - (3./4.)*fPi;}
                if(fCosThetaCSRec[j] > 0.){fPhiTilde = tmpVar - (1./4.)*fPi;}
                if(fPhiTilde > 2*fPi){fPhiTilde = fPhiTilde - 2*fPi;}
                if(fPhiTilde < 0.){fPhiTilde = 2*fPi + fPhiTilde;}

                fHistRecCosThetaCSReWeighted[indexPt] -> Fill(fCosThetaCSRec[j],weightCosThetaCS);
                fHistRecPhiCSReWeighted[indexPt] -> Fill(TMath::Abs(fPhiCSRec[j]),weightPhiCS);
                fHistRecPhiTildeCSReWeighted[indexPt] -> Fill(fPhiTilde,weightPhiTildeCS);
                fHistRecCosThetaPhiCSReWeighted[indexPt] -> Fill(fCosThetaCSRec[j],TMath::Abs(fPhiCSRec[j]),weightCosThetaPhiCS);
                indexPt = 0;
              }
            }
            ////////////////////////////////////////////////////////////////////
          }
        }
      }
    }

  }
  printf("\n");

  for(int i = 0;i < fNPtBins;i++){
    fHistAccxEffCosThetaHEReWeighted[i] = new TH1D(Form("histAccxEffCosThetaHEReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]);
    fHistAccxEffCosThetaHEReWeighted[i] -> Divide(fHistRecCosThetaHEReWeighted[i],fHistGenCosThetaHEReWeighted[i],1,1,"B");

    fHistAccxEffPhiHEReWeighted[i] = new TH1D(Form("histAccxEffPhiHEReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistAccxEffPhiHEReWeighted[i] -> Divide(fHistRecPhiHEReWeighted[i],fHistGenPhiHEReWeighted[i],1,1,"B");

    fHistAccxEffPhiTildeHEReWeighted[i] = new TH1D(Form("histAccxEffPhiTildeHEReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]);
    fHistAccxEffPhiTildeHEReWeighted[i] -> Divide(fHistRecPhiTildeHEReWeighted[i],fHistGenPhiTildeHEReWeighted[i],1,1,"B");

    fHistAccxEffCosThetaPhiHEReWeighted[i] = new TH2D(Form("histAccxEffCosThetaPhiHEReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]);
    fHistAccxEffCosThetaPhiHEReWeighted[i] -> Divide(fHistRecCosThetaPhiHEReWeighted[i],fHistGenCosThetaPhiHEReWeighted[i],1,1,"B");


    fHistAccxEffCosThetaCSReWeighted[i] = new TH1D(Form("histAccxEffCosThetaCSReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]);
    fHistAccxEffCosThetaCSReWeighted[i] -> Divide(fHistRecCosThetaCSReWeighted[i],fHistGenCosThetaCSReWeighted[i],1,1,"B");

    fHistAccxEffPhiCSReWeighted[i] = new TH1D(Form("histAccxEffPhiCSReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistAccxEffPhiCSReWeighted[i] -> Divide(fHistRecPhiCSReWeighted[i],fHistGenPhiCSReWeighted[i],1,1,"B");

    fHistAccxEffPhiTildeCSReWeighted[i] = new TH1D(Form("histAccxEffPhiTildeCSReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]);
    fHistAccxEffPhiTildeCSReWeighted[i] -> Divide(fHistRecPhiTildeCSReWeighted[i],fHistGenPhiTildeCSReWeighted[i],1,1,"B");

    fHistAccxEffCosThetaPhiCSReWeighted[i] = new TH2D(Form("histAccxEffCosThetaPhiCSReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]);
    fHistAccxEffCosThetaPhiCSReWeighted[i] -> Divide(fHistRecCosThetaPhiCSReWeighted[i],fHistGenCosThetaPhiCSReWeighted[i],1,1,"B");
  }

  if(saveFile){
    TFile *fileAccxEffReWeighted = new TFile(nameOutputFile.c_str(),"RECREATE");
    for(int i = 1;i < fNPtBins;i++){
      fHistGenCosThetaHEReWeighted[i] -> Write();
      fHistRecCosThetaHEReWeighted[i] -> Write();
      fHistAccxEffCosThetaHEReWeighted[i] -> Write();

      fHistGenPhiHEReWeighted[i] -> Write();
      fHistRecPhiHEReWeighted[i] -> Write();
      fHistAccxEffPhiHEReWeighted[i] -> Write();

      fHistGenPhiTildeHEReWeighted[i] -> Write();
      fHistRecPhiTildeHEReWeighted[i] -> Write();
      fHistAccxEffPhiTildeHEReWeighted[i] -> Write();

      fHistGenCosThetaCSReWeighted[i] -> Write();
      fHistRecCosThetaCSReWeighted[i] -> Write();
      fHistAccxEffCosThetaCSReWeighted[i] -> Write();

      fHistGenPhiCSReWeighted[i] -> Write();
      fHistRecPhiCSReWeighted[i] -> Write();
      fHistAccxEffPhiCSReWeighted[i] -> Write();

      fHistGenPhiTildeCSReWeighted[i] -> Write();
      fHistRecPhiTildeCSReWeighted[i] -> Write();
      fHistAccxEffPhiTildeCSReWeighted[i] -> Write();

      fHistGenCosThetaPhiHEReWeighted[i] -> Write();
      fHistRecCosThetaPhiHEReWeighted[i] -> Write();
      fHistAccxEffCosThetaPhiHEReWeighted[i] -> Write();

      fHistGenCosThetaPhiCSReWeighted[i] -> Write();
      fHistRecCosThetaPhiCSReWeighted[i] -> Write();
      fHistAccxEffCosThetaPhiCSReWeighted[i] -> Write();
    }
    fileAccxEffReWeighted -> Close();
  }
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
//______________________________________________________________________________
//void AccxEffCalculator::ComputeReweightTRFAccxEff(string strSample, string nameOutputFile, bool reweightAccxEff, TH1D *histReweightTRF) {
void AccxEffCalculator::ComputeReweightTRFAccxEff(string strSample, string nameOutputFile, bool reweightAccxEff, TObjArray *objArrayReweightTRF) {
  int nEvents = 0;
  int indexPt = 0;

  TH1D *histReweightTRF_25eta4 = (TH1D*) objArrayReweightTRF -> At(0);
  TH1D *histReweightTRF_25eta3 = (TH1D*) objArrayReweightTRF -> At(1);
  TH1D *histReweightTRF_3eta35 = (TH1D*) objArrayReweightTRF -> At(2);
  TH1D *histReweightTRF_35eta4 = (TH1D*) objArrayReweightTRF -> At(3);

  if(strSample == "FullStat"){nEvents = fTreeAccxEff -> GetEntries();}
  if(strSample == "TestStat"){nEvents = 1000;}
  printf("N events = %i \n",nEvents);

  for(int i = 0;i < fNPtBins;i++){
    // HELICITY
    fHistGenCosThetaHE[i] = new TH1D(Form("histGenCosThetaHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistGenCosThetaHE[i] -> Sumw2();
    fHistRecCosThetaHE[i] = new TH1D(Form("histRecCosThetaHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistRecCosThetaHE[i] -> Sumw2();
    fHistGenPhiHE[i] = new TH1D(Form("histGenPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistGenPhiHE[i] -> Sumw2();
    fHistRecPhiHE[i] = new TH1D(Form("histRecPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistRecPhiHE[i] -> Sumw2();
    fHistGenPhiTildeHE[i] = new TH1D(Form("histGenPhiTildeHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]); fHistGenPhiTildeHE[i] -> Sumw2();
    fHistRecPhiTildeHE[i] = new TH1D(Form("histRecPhiTildeHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]); fHistRecPhiTildeHE[i] -> Sumw2();
    fHistGenCosThetaPhiHE[i] = new TH2D(Form("histGenCosThetaPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistGenCosThetaPhiHE[i] -> Sumw2();
    fHistRecCosThetaPhiHE[i] = new TH2D(Form("histRecCosThetaPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistRecCosThetaPhiHE[i] -> Sumw2();

    // COLLINS-SOPER
    fHistGenCosThetaCS[i] = new TH1D(Form("histGenCosThetaCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistGenCosThetaCS[i] -> Sumw2();
    fHistRecCosThetaCS[i] = new TH1D(Form("histRecCosThetaCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistRecCosThetaCS[i] -> Sumw2();
    fHistGenPhiCS[i] = new TH1D(Form("histGenPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistGenPhiCS[i] -> Sumw2();
    fHistRecPhiCS[i] = new TH1D(Form("histRecPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistRecPhiCS[i] -> Sumw2();
    fHistGenPhiTildeCS[i] = new TH1D(Form("histGenPhiTildeCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]); fHistGenPhiTildeCS[i] -> Sumw2();
    fHistRecPhiTildeCS[i] = new TH1D(Form("histRecPhiTildeCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]); fHistRecPhiTildeCS[i] -> Sumw2();
    fHistGenCosThetaPhiCS[i] = new TH2D(Form("histGenCosThetaPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistGenCosThetaPhiCS[i] -> Sumw2();
    fHistRecCosThetaPhiCS[i] = new TH2D(Form("histRecCosThetaPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistRecCosThetaPhiCS[i] -> Sumw2();
  }

  fHistGenCosThetaHEPt = new TH2D("histGenCosThetaHEPt","",100,-1,1,100,0,10); fHistGenCosThetaHEPt -> Sumw2();
  fHistRecCosThetaHEPt = new TH2D("histRecCosThetaHEPt","",100,-1,1,100,0,10); fHistRecCosThetaHEPt -> Sumw2();
  fHistGenPhiHEPt = new TH2D("histGenPhiHEPt","",100,0,fPi,100,0,10); fHistGenPhiHEPt -> Sumw2();
  fHistRecPhiHEPt = new TH2D("histRecPhiHEPt","",100,0,fPi,100,0,10); fHistRecPhiHEPt -> Sumw2();
  fHistGenPhiTildeHEPt = new TH2D("histGenPhiTildeHEPt","",100,0.,2*fPi,100,0,15); fHistGenPhiTildeHEPt -> Sumw2();
  fHistRecPhiTildeHEPt = new TH2D("histRecPhiTildeHEPt","",100,0.,2*fPi,100,0,15); fHistRecPhiTildeHEPt -> Sumw2();

  fHistGenCosThetaCSPt = new TH2D("histGenCosThetaCSPt","",100,-1,1,100,0,10); fHistGenCosThetaCSPt -> Sumw2();
  fHistRecCosThetaCSPt = new TH2D("histRecCosThetaCSPt","",100,-1,1,100,0,10); fHistRecCosThetaCSPt -> Sumw2();
  fHistGenPhiCSPt = new TH2D("histGenPhiCSPt","",100,0,fPi,100,0,10); fHistGenPhiCSPt -> Sumw2();
  fHistRecPhiCSPt = new TH2D("histRecPhiCSPt","",100,0,fPi,100,0,10); fHistRecPhiCSPt -> Sumw2();
  fHistGenPhiTildeCSPt = new TH2D("histGenPhiTildeCSPt","",100,0.,2*fPi,100,0,15); fHistGenPhiTildeCSPt -> Sumw2();
  fHistRecPhiTildeCSPt = new TH2D("histRecPhiTildeCSPt","",100,0.,2*fPi,100,0,15); fHistRecPhiTildeCSPt -> Sumw2();

  fHistRecCosThetaHEEtaSM = new TH2D("histRecCosThetaHEEtaSM","",100,-1.,1.,100,-4.,-2.5);
  fHistRecCosThetaCSEtaSM = new TH2D("histRecCosThetaCSEtaSM","",100,-1.,1.,100,-4.,-2.5);
  fHistRecPhiHEEtaSM = new TH2D("histRecPhiHEEtaSM","",100,-fPi,fPi,100,-4.,-2.5);
  fHistRecPhiCSEtaSM = new TH2D("histRecPhiCSEtaSM","",100,-fPi,fPi,100,-4.,-2.5);
  fHistRecCosThetaPhiHESM = new TH2D("histRecCosThetaPhiHESM","",100,-1.,1.,100,-fPi,fPi);
  fHistRecCosThetaPhiCSSM = new TH2D("histRecCosThetaPhiCSSM","",100,-1.,1.,100,-fPi,fPi);

  printf("- Configuring Fiducial Box (only for Reconstructed events) \n");
  printf("%f < CosTheta < %f \n",fCosThetaValues[1],fCosThetaValues[fNCosThetaBins-1]);
  printf("%f < |Phi| < %f \n",fPhiValues[1],fPhiValues[fNPhiBins-1]);

  double tmpVar = 0;
  double weightMu1 = 0, weightMu2 = 0, weigthMu1Mu2 = 0;
  double etaMu;

  for(int i = 0;i < nEvents;i++){
    printf("Reading : %3.2f %% \r",((double) i/(double) nEvents)*100.);
    fTreeAccxEff -> GetEntry(i);

    for(int j = 0;j < fNDimuGen;j++){
      if(fDimuYGen[j] > -4. && fDimuYGen[j] < -2.5){
        ////////////////////////////////////////////////////////////////////////
        // HELICITY
        while(fDimuPtGen[j] < fMinPt[indexPt] || fDimuPtGen[j] > fMaxPt[indexPt]){indexPt++;}
        // 1D approach
        fPhiTilde = computePhiTilde(fCosThetaHEGen[j],fPhiHEGen[j]);

        fHistGenCosThetaHE[indexPt] -> Fill(fCosThetaHEGen[j]);
        fHistGenPhiHE[indexPt] -> Fill(TMath::Abs(fPhiHEGen[j]));
        fHistGenPhiTildeHE[indexPt] -> Fill(fPhiTilde);

        fHistGenCosThetaHEPt -> Fill(fCosThetaHEGen[j],fDimuPtGen[j]);
        fHistGenPhiHEPt -> Fill(TMath::Abs(fPhiHEGen[j]),fDimuPtGen[j]);
        fHistGenPhiTildeHEPt -> Fill(fPhiTilde,fDimuPtGen[j]);

        // 2D approach
        fHistGenCosThetaPhiHE[indexPt] -> Fill(fCosThetaHEGen[j],TMath::Abs(fPhiHEGen[j]));
        indexPt = 0;

        // COLLINS-SOPER
        while(fDimuPtGen[j] < fMinPt[indexPt] || fDimuPtGen[j] > fMaxPt[indexPt]){indexPt++;}
        // 1D approach
        fPhiTilde = computePhiTilde(fCosThetaCSGen[j],fPhiCSGen[j]);

        fHistGenCosThetaCS[indexPt] -> Fill(fCosThetaCSGen[j]);
        fHistGenPhiCS[indexPt] -> Fill(TMath::Abs(fPhiCSGen[j]));
        fHistGenPhiTildeCS[indexPt] -> Fill(fPhiTilde);

        fHistGenCosThetaCSPt -> Fill(fCosThetaCSGen[j],fDimuPtGen[j]);
        fHistGenPhiCSPt -> Fill(TMath::Abs(fPhiCSGen[j]),fDimuPtGen[j]);
        fHistGenPhiTildeCSPt -> Fill(fPhiTilde,fDimuPtGen[j]);

        // 2D approach
        fHistGenCosThetaPhiCS[indexPt] -> Fill(fCosThetaCSGen[j],TMath::Abs(fPhiCSGen[j]));
        indexPt = 0;
        ////////////////////////////////////////////////////////////////////////
      }
    }

    // The reweighting with TRF is only for Rec because Gen have no trigger!
    for(int j = 0;j < fNDimuRec;j++){
      if(fDimuYRec[j] > -4. && fDimuYRec[j] < -2.5){
        if(fDimuMatchRec[j] == 2){
          if(fDimuMassRec[j] > 2 && fDimuMassRec[j] < 5){
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(fChargeRec[0] > 0.){etaMu = fEtaRec[0];}
            else{etaMu = fEtaRec[1];}
            fHistRecCosThetaHEEtaSM -> Fill(fCosThetaHERec[j],etaMu);
            fHistRecCosThetaCSEtaSM -> Fill(fCosThetaCSRec[j],etaMu);
            fHistRecPhiHEEtaSM -> Fill(fPhiHERec[j],etaMu);
            fHistRecPhiCSEtaSM -> Fill(fPhiCSRec[j],etaMu);
            //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if(reweightAccxEff){
              if(strstr(nameOutputFile.c_str(),"INTEGRATED")){
                weightMu1 = histReweightTRF_25eta4 -> GetBinContent(histReweightTRF_25eta4 -> FindBin(fPtRec[0]));
                weightMu2 = histReweightTRF_25eta4 -> GetBinContent(histReweightTRF_25eta4 -> FindBin(fPtRec[1]));
              }

              if(strstr(nameOutputFile.c_str(),"DIFFERENTIAL")){
                if(fEtaRec[0] > -3.0 && fEtaRec[0] < -2.5){weightMu1 = histReweightTRF_25eta3 -> GetBinContent(histReweightTRF_25eta3 -> FindBin(fPtRec[0]));}
                if(fEtaRec[0] > -3.5 && fEtaRec[0] < -3.0){weightMu1 = histReweightTRF_3eta35 -> GetBinContent(histReweightTRF_3eta35 -> FindBin(fPtRec[0]));}
                if(fEtaRec[0] > -4.0 && fEtaRec[0] < -3.5){weightMu1 = histReweightTRF_35eta4 -> GetBinContent(histReweightTRF_35eta4 -> FindBin(fPtRec[0]));}

                if(fEtaRec[1] > -3.0 && fEtaRec[1] < -2.5){weightMu2 = histReweightTRF_25eta3 -> GetBinContent(histReweightTRF_25eta3 -> FindBin(fPtRec[1]));}
                if(fEtaRec[1] > -3.5 && fEtaRec[1] < -3.0){weightMu2 = histReweightTRF_3eta35 -> GetBinContent(histReweightTRF_3eta35 -> FindBin(fPtRec[1]));}
                if(fEtaRec[1] > -4.0 && fEtaRec[1] < -3.5){weightMu2 = histReweightTRF_35eta4 -> GetBinContent(histReweightTRF_35eta4 -> FindBin(fPtRec[1]));}
              }

              weigthMu1Mu2 = weightMu1*weightMu2;

              //cout << "1) Eta = " << fEtaRec[0] << " ; Pt = " << fPtRec[0] << " ; weight = " << weightMu1 << endl;
              //cout << "2) Eta = " << fEtaRec[1] << " ; Pt = " << fPtRec[1] << " ; weight = " << weightMu2 << endl;
              //cout << "====================================================================" << endl;
            }
            else{weigthMu1Mu2 = 1.;}
            ////////////////////////////////////////////////////////////////////
            // HELICITY
            if(TMath::Abs(fPhiHERec[j]) > fPhiValues[1] && TMath::Abs(fPhiHERec[j]) < fPhiValues[fNPhiBins-1]){
              if(fCosThetaHERec[j] > fCosThetaValues[1] && fCosThetaHERec[j] < fCosThetaValues[fNCosThetaBins-1]){
                while(fDimuPtRec[j] < fMinPt[indexPt] || fDimuPtRec[j] > fMaxPt[indexPt]){indexPt++;}
                // 1D approach
                fPhiTilde = computePhiTilde(fCosThetaHERec[j],fPhiHERec[j]);

                fHistRecCosThetaHE[indexPt] -> Fill(fCosThetaHERec[j],weigthMu1Mu2);
                fHistRecPhiHE[indexPt] -> Fill(TMath::Abs(fPhiHERec[j]),weigthMu1Mu2);
                fHistRecPhiTildeHE[indexPt] -> Fill(fPhiTilde,weigthMu1Mu2);

                fHistRecCosThetaHEPt -> Fill(fCosThetaHERec[j],fDimuPtRec[j],weigthMu1Mu2);
                fHistRecPhiHEPt -> Fill(TMath::Abs(fPhiHERec[j]),fDimuPtRec[j],weigthMu1Mu2);
                fHistRecPhiTildeHEPt -> Fill(fPhiTilde,fDimuPtRec[j],weigthMu1Mu2);

                // 2D approach
                fHistRecCosThetaPhiHE[indexPt] -> Fill(fCosThetaHERec[j],TMath::Abs(fPhiHERec[j]),weigthMu1Mu2);
                indexPt = 0;
              }
            }
            // COLLINS-SOPER
            if(TMath::Abs(fPhiCSRec[j]) > fPhiValues[1] && TMath::Abs(fPhiCSRec[j]) < fPhiValues[fNPhiBins-1]){
              if(fCosThetaCSRec[j] > fCosThetaValues[1] && fCosThetaCSRec[j] < fCosThetaValues[fNCosThetaBins-1]){
                while(fDimuPtRec[j] < fMinPt[indexPt] || fDimuPtRec[j] > fMaxPt[indexPt]){indexPt++;}
                // 1D approach
                fPhiTilde = computePhiTilde(fCosThetaCSRec[j],fPhiCSRec[j]);

                fHistRecCosThetaCS[indexPt] -> Fill(fCosThetaCSRec[j],weigthMu1Mu2);
                fHistRecPhiCS[indexPt] -> Fill(TMath::Abs(fPhiCSRec[j]),weigthMu1Mu2);
                fHistRecPhiTildeCS[indexPt] -> Fill(fPhiTilde,weigthMu1Mu2);

                fHistRecCosThetaCSPt -> Fill(fCosThetaCSRec[j],fDimuPtRec[j],weigthMu1Mu2);
                fHistRecPhiCSPt -> Fill(TMath::Abs(fPhiCSRec[j]),fDimuPtRec[j],weigthMu1Mu2);
                fHistRecPhiTildeCSPt -> Fill(fPhiTilde,fDimuPtRec[j],weigthMu1Mu2);

                // 2D approach
                fHistRecCosThetaPhiCS[indexPt] -> Fill(fCosThetaCSRec[j],TMath::Abs(fPhiCSRec[j]),weigthMu1Mu2);
                indexPt = 0;
              }
            }
            ////////////////////////////////////////////////////////////////////
          }
        }
      }
    }

  }
  printf("\n");

  for(int i = 0;i < fNPtBins;i++){
    // HELICITY
    fHistAccxEffCosThetaHE[i] = new TH1D(Form("histAccxEffCosThetaHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]);
    fHistAccxEffCosThetaHE[i] -> Divide(fHistRecCosThetaHE[i],fHistGenCosThetaHE[i],1,1,"B");

    fHistAccxEffPhiHE[i] = new TH1D(Form("histAccxEffPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistAccxEffPhiHE[i] -> Divide(fHistRecPhiHE[i],fHistGenPhiHE[i],1,1,"B");

    fHistAccxEffPhiTildeHE[i] = new TH1D(Form("histAccxEffPhiTildeHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]);
    fHistAccxEffPhiTildeHE[i] -> Divide(fHistRecPhiTildeHE[i],fHistGenPhiTildeHE[i],1,1,"B");

    fHistAccxEffCosThetaPhiHE[i] = new TH2D(Form("histAccxEffCosThetaPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]);
    fHistAccxEffCosThetaPhiHE[i] -> Divide(fHistRecCosThetaPhiHE[i],fHistGenCosThetaPhiHE[i],1,1,"B");

    // COLLINS-SOPER
    fHistAccxEffCosThetaCS[i] = new TH1D(Form("histAccxEffCosThetaCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]);
    fHistAccxEffCosThetaCS[i] -> Divide(fHistRecCosThetaCS[i],fHistGenCosThetaCS[i],1,1,"B");

    fHistAccxEffPhiCS[i] = new TH1D(Form("histAccxEffPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistAccxEffPhiCS[i] -> Divide(fHistRecPhiCS[i],fHistGenPhiCS[i],1,1,"B");

    fHistAccxEffPhiTildeCS[i] = new TH1D(Form("histAccxEffPhiTildeCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiTildeBins,&fPhiTildeValues[0]);
    fHistAccxEffPhiTildeCS[i] -> Divide(fHistRecPhiTildeCS[i],fHistGenPhiTildeCS[i],1,1,"B");

    fHistAccxEffCosThetaPhiCS[i] = new TH2D(Form("histAccxEffCosThetaPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]);
    fHistAccxEffCosThetaPhiCS[i] -> Divide(fHistRecCosThetaPhiCS[i],fHistGenCosThetaPhiCS[i],1,1,"B");
  }

  fHistAccxEffCosThetaHEPt = new TH2D("histAccxEffCosThetaHEPt","",100,-1,1,100,0,10);
  fHistAccxEffCosThetaHEPt -> Divide(fHistRecCosThetaHEPt,fHistGenCosThetaHEPt,1,1,"B");
  fHistAccxEffPhiHEPt = new TH2D("histAccxEffPhiHEPt","",100,0,fPi,100,0,10);
  fHistAccxEffPhiHEPt -> Divide(fHistRecPhiHEPt,fHistGenPhiHEPt,1,1,"B");
  fHistAccxEffPhiTildeHEPt = new TH2D("histAccxEffPhiTildeHEPt","",100,0.,2*fPi,100,0,15);
  fHistAccxEffPhiTildeHEPt -> Divide(fHistRecPhiTildeHEPt,fHistGenPhiTildeHEPt,1,1,"B");

  fHistAccxEffCosThetaCSPt = new TH2D("histAccxEffCosThetaCSPt","",100,-1,1,100,0,10);
  fHistAccxEffCosThetaCSPt -> Divide(fHistRecCosThetaCSPt,fHistGenCosThetaCSPt,1,1,"B");
  fHistAccxEffPhiCSPt = new TH2D("histAccxEffPhiCSPt","",100,0,fPi,100,0,10);
  fHistAccxEffPhiCSPt -> Divide(fHistRecPhiCSPt,fHistGenPhiCSPt,1,1,"B");
  fHistAccxEffPhiTildeCSPt = new TH2D("histAccxEffPhiTildeCSPt","",100,0.,2*fPi,100,0,15);
  fHistAccxEffPhiTildeCSPt -> Divide(fHistRecPhiTildeCSPt,fHistGenPhiTildeCSPt,1,1,"B");

  TFile *fileAccxEff = new TFile(nameOutputFile.c_str(),"RECREATE");
  for(int i = 0;i < fNPtBins;i++){
    fHistGenCosThetaHE[i] -> Write(); fHistRecCosThetaHE[i] -> Write(); fHistAccxEffCosThetaHE[i] -> Write();
    fHistGenPhiHE[i] -> Write(); fHistRecPhiHE[i] -> Write(); fHistAccxEffPhiHE[i] -> Write();
    fHistGenPhiTildeHE[i] -> Write(); fHistRecPhiTildeHE[i] -> Write(); fHistAccxEffPhiTildeHE[i] -> Write();

    fHistGenCosThetaPhiHE[i] -> Write(); fHistRecCosThetaPhiHE[i] -> Write(); fHistAccxEffCosThetaPhiHE[i] -> Write();

    fHistGenCosThetaCS[i] -> Write(); fHistRecCosThetaCS[i] -> Write(); fHistAccxEffCosThetaCS[i] -> Write();
    fHistGenPhiCS[i] -> Write(); fHistRecPhiCS[i] -> Write(); fHistAccxEffPhiCS[i] -> Write();
    fHistGenPhiTildeCS[i] -> Write(); fHistRecPhiTildeCS[i] -> Write(); fHistAccxEffPhiTildeCS[i] -> Write();

    fHistGenCosThetaPhiCS[i] -> Write(); fHistRecCosThetaPhiCS[i] -> Write(); fHistAccxEffCosThetaPhiCS[i] -> Write();
  }

  fHistGenCosThetaHEPt -> Write(); fHistGenPhiHEPt -> Write(); fHistAccxEffCosThetaHEPt -> Write();
  fHistRecCosThetaHEPt -> Write(); fHistRecPhiHEPt -> Write(); fHistAccxEffPhiHEPt -> Write();

  fHistGenCosThetaCSPt -> Write(); fHistGenPhiCSPt -> Write(); fHistAccxEffCosThetaCSPt -> Write();
  fHistRecCosThetaCSPt -> Write(); fHistRecPhiCSPt -> Write(); fHistAccxEffPhiCSPt -> Write();

  fHistGenPhiTildeHEPt -> Write(); fHistRecPhiTildeHEPt -> Write(); fHistAccxEffPhiTildeHEPt -> Write();
  fHistGenPhiTildeCSPt -> Write(); fHistRecPhiTildeCSPt -> Write(); fHistAccxEffPhiTildeCSPt -> Write();

  fHistRecCosThetaHEEtaSM -> Write();
  fHistRecCosThetaCSEtaSM -> Write();
  fHistRecPhiHEEtaSM -> Write();
  fHistRecPhiCSEtaSM -> Write();

  fileAccxEff -> Close();
}
//______________________________________________________________________________
////////////////////////////////////////////////////////////////////////////////
Double_t MyFuncPol(Double_t *x, Double_t *par){
  double cosTheta = x[0];
  double phi = x[1];

  double N = par[0];
  double lambdaTheta = par[1];
  double lambdaPhi = par[2];
  double lambdaThetaPhi = par[3];

  double cosPhi = TMath::Cos(phi);
  double cos2Phi = TMath::Cos(2*phi);

  // WARNING!!!!!!
  //return (N/(3 + lambdaTheta))*(1 + (lambdaTheta + lambdaPhi*cos2Phi)*cosTheta*cosTheta + 2*lambdaThetaPhi*cosTheta*cosPhi*TMath::Sqrt(1 - cosTheta*cosTheta) + lambdaPhi*cos2Phi);
  return (N/(3 + lambdaTheta))*(1 + (lambdaTheta - lambdaPhi*cos2Phi)*cosTheta*cosTheta + 2*lambdaThetaPhi*cosTheta*cosPhi*TMath::Sqrt(1 - cosTheta*cosTheta) + lambdaPhi*cos2Phi);
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Double_t computePhiTilde(Double_t CosTheta, Double_t Phi){
  double fPi = TMath::Pi();
  double PhiTilde = 0;
  double tmpVar = Phi + fPi;
  if(CosTheta < 0.){PhiTilde = tmpVar - (3./4.)*fPi;}
  if(CosTheta > 0.){PhiTilde = tmpVar - (1./4.)*fPi;}
  if(PhiTilde > 2*fPi){PhiTilde = PhiTilde - 2*fPi;}
  if(PhiTilde < 0.){PhiTilde = 2*fPi + PhiTilde;}
  return PhiTilde;
}
