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
#include "AccxEffCalculator.h"

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
void AccxEffCalculator::SetBinning(vector <Double_t> CosThetaValues, vector <Double_t> PhiValues) {
  for(int i = 0;i < (int) CosThetaValues.size();i++){fCosThetaValues.push_back(CosThetaValues[i]);}
  fNCosThetaBins = fCosThetaValues.size() - 1;
  for(int i = 0;i < (int) PhiValues.size();i++){fPhiValues.push_back(PhiValues[i]);}
  fNPhiBins = fPhiValues.size() - 1;
}
//______________________________________________________________________________
void AccxEffCalculator::ComputeAccxEff(string strSample, string nameOutputFile) {
  int nEvents = 0;
  int indexPt = 0;

  if(strSample == "FullStat"){nEvents = fTreeAccxEff -> GetEntries();}
  if(strSample == "TestStat"){nEvents = 1000000;}
  printf("N events = %i \n",nEvents);

  for(int i = 0;i < fNPtBins;i++){
    // HELICITY
    fHistGenCosThetaHE[i] = new TH1D(Form("histGenCosThetaHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistGenCosThetaHE[i] -> Sumw2();
    fHistRecCosThetaHE[i] = new TH1D(Form("histRecCosThetaHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistRecCosThetaHE[i] -> Sumw2();
    fHistGenPhiHE[i] = new TH1D(Form("histGenPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistGenPhiHE[i] -> Sumw2();
    fHistRecPhiHE[i] = new TH1D(Form("histRecPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistRecPhiHE[i] -> Sumw2();
    fHistGenPhiTildeHE[i] = new TH1D(Form("histGenPhiTildeHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",10,0.,fPi); fHistGenPhiTildeHE[i] -> Sumw2();
    fHistRecPhiTildeHE[i] = new TH1D(Form("histRecPhiTildeHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",10,0.,fPi); fHistRecPhiTildeHE[i] -> Sumw2();
    fHistGenCosThetaPhiHE[i] = new TH2D(Form("histGenCosThetaPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistGenCosThetaPhiHE[i] -> Sumw2();
    fHistRecCosThetaPhiHE[i] = new TH2D(Form("histRecCosThetaPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistRecCosThetaPhiHE[i] -> Sumw2();

    // COLLINS-SOPER
    fHistGenCosThetaCS[i] = new TH1D(Form("histGenCosThetaCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistGenCosThetaCS[i] -> Sumw2();
    fHistRecCosThetaCS[i] = new TH1D(Form("histRecCosThetaCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]); fHistRecCosThetaCS[i] -> Sumw2();
    fHistGenPhiCS[i] = new TH1D(Form("histGenPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistGenPhiCS[i] -> Sumw2();
    fHistRecPhiCS[i] = new TH1D(Form("histRecPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]); fHistRecPhiCS[i] -> Sumw2();
    fHistGenPhiTildeCS[i] = new TH1D(Form("histGenPhiTildeCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",10,0.,fPi); fHistGenPhiTildeCS[i] -> Sumw2();
    fHistRecPhiTildeCS[i] = new TH1D(Form("histRecPhiTildeCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",10,0.,fPi); fHistRecPhiTildeCS[i] -> Sumw2();
    fHistGenCosThetaPhiCS[i] = new TH2D(Form("histGenCosThetaPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistGenCosThetaPhiCS[i] -> Sumw2();
    fHistRecCosThetaPhiCS[i] = new TH2D(Form("histRecCosThetaPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]); fHistRecCosThetaPhiCS[i] -> Sumw2();
  }

  fHistGenCosThetaHEPt = new TH2D("histGenCosThetaHEPt","",100,-1,1,100,0,10); fHistGenCosThetaHEPt -> Sumw2();
  fHistRecCosThetaHEPt = new TH2D("histRecCosThetaHEPt","",100,-1,1,100,0,10); fHistRecCosThetaHEPt -> Sumw2();
  fHistGenPhiHEPt = new TH2D("histGenPhiHEPt","",100,0,fPi,100,0,10); fHistGenPhiHEPt -> Sumw2();
  fHistRecPhiHEPt = new TH2D("histRecPhiHEPt","",100,0,fPi,100,0,10); fHistRecPhiHEPt -> Sumw2();

  fHistGenCosThetaCSPt = new TH2D("histGenCosThetaCSPt","",100,-1,1,100,0,10); fHistGenCosThetaCSPt -> Sumw2();
  fHistRecCosThetaCSPt = new TH2D("histRecCosThetaCSPt","",100,-1,1,100,0,10); fHistRecCosThetaCSPt -> Sumw2();
  fHistGenPhiCSPt = new TH2D("histGenPhiCSPt","",100,0,fPi,100,0,10); fHistGenPhiCSPt -> Sumw2();
  fHistRecPhiCSPt = new TH2D("histRecPhiCSPt","",100,0,fPi,100,0,10); fHistRecPhiCSPt -> Sumw2();

  for(int i = 0;i < nEvents;i++){
    printf("Reading : %3.2f %% \r",((double) i/(double) nEvents)*100.);
    fTreeAccxEff -> GetEntry(i);

    for(int j = 0;j < fNDimuGen;j++){
      if(fDimuYGen[j] > -4. && fDimuYGen[j] < -2.5){
        fHistGenCosThetaHEPt -> Fill(fCosThetaHEGen[j],fDimuPtGen[j]);
        fHistGenPhiHEPt -> Fill(TMath::Abs(fPhiHEGen[j]),fDimuPtGen[j]);
        fHistGenCosThetaCSPt -> Fill(fCosThetaCSGen[j],fDimuPtGen[j]);
        fHistGenPhiCSPt -> Fill(TMath::Abs(fPhiCSGen[j]),fDimuPtGen[j]);
        while(fDimuPtGen[j] < fMinPt[indexPt] || fDimuPtGen[j] > fMaxPt[indexPt]){indexPt++;}
        // 1D approach
        fHistGenCosThetaHE[indexPt] -> Fill(fCosThetaHEGen[j]);
        fHistGenPhiHE[indexPt] -> Fill(TMath::Abs(fPhiHEGen[j]));
        if(fCosThetaHEGen[j] < 0.){fPhiTilde = fPhiHEGen[j] - (3./4.)*fPi;}
        if(fCosThetaHEGen[j] > 0.){fPhiTilde = fPhiHEGen[j] - (1./4.)*fPi;}
        if(fPhiTilde < 0){fPhiTilde = 2*fPi + fPhiTilde;}
        if(fPhiTilde > fPi){fPhiTilde = 2*fPi - fPhiTilde;}
        fHistGenPhiTildeHE[indexPt] -> Fill(fPhiTilde);

        fHistGenCosThetaCS[indexPt] -> Fill(fCosThetaCSGen[j]);
        fHistGenPhiCS[indexPt] -> Fill(TMath::Abs(fPhiCSGen[j]));
        if(fCosThetaCSGen[j] < 0.){fPhiTilde = fPhiCSGen[j] - (3./4.)*fPi;}
        if(fCosThetaCSGen[j] > 0.){fPhiTilde = fPhiCSGen[j] - (1./4.)*fPi;}
        if(fPhiTilde < 0){fPhiTilde = 2*fPi + fPhiTilde;}
        if(fPhiTilde > fPi){fPhiTilde = 2*fPi - fPhiTilde;}
        fHistGenPhiTildeCS[indexPt] -> Fill(fPhiTilde);
        // 2D approach
        fHistGenCosThetaPhiHE[indexPt] -> Fill(fCosThetaHEGen[j],TMath::Abs(fPhiHEGen[j]));
        fHistGenCosThetaPhiCS[indexPt] -> Fill(fCosThetaCSGen[j],TMath::Abs(fPhiCSGen[j]));
        indexPt = 0;
      }
    }

    for(int j = 0;j < fNDimuRec;j++){
      if(fDimuYRec[j] > -4. && fDimuYRec[j] < -2.5){
        if(fDimuMatchRec[j] == 2){
          if(fDimuMassRec[j] > 2 && fDimuMassRec[j] < 5){
            fHistRecCosThetaHEPt -> Fill(fCosThetaHERec[j],fDimuPtRec[j]);
            fHistRecPhiHEPt -> Fill(TMath::Abs(fPhiHERec[j]),fDimuPtRec[j]);
            fHistRecCosThetaCSPt -> Fill(fCosThetaCSRec[j],fDimuPtRec[j]);
            fHistRecPhiCSPt -> Fill(TMath::Abs(fPhiCSRec[j]),fDimuPtRec[j]);
            while(fDimuPtRec[j] < fMinPt[indexPt] || fDimuPtRec[j] > fMaxPt[indexPt]){indexPt++;}
            // 1D approach
            fHistRecCosThetaHE[indexPt] -> Fill(fCosThetaHERec[j]);
            fHistRecPhiHE[indexPt] -> Fill(TMath::Abs(fPhiHERec[j]));
            if(fCosThetaHERec[j] < 0.){fPhiTilde = fPhiHERec[j] - (3./4.)*fPi;}
            if(fCosThetaHERec[j] > 0.){fPhiTilde = fPhiHERec[j] - (1./4.)*fPi;}
            if(fPhiTilde < 0){fPhiTilde = 2*fPi + fPhiTilde;}
            if(fPhiTilde > fPi){fPhiTilde = 2*fPi - fPhiTilde;}
            fHistRecPhiTildeHE[indexPt] -> Fill(fPhiTilde);

            fHistRecCosThetaCS[indexPt] -> Fill(fCosThetaCSRec[j]);
            fHistRecPhiCS[indexPt] -> Fill(TMath::Abs(fPhiCSRec[j]));
            if(fCosThetaCSRec[j] < 0.){fPhiTilde = fPhiCSRec[j] - (3./4.)*fPi;}
            if(fCosThetaCSRec[j] > 0.){fPhiTilde = fPhiCSRec[j] - (1./4.)*fPi;}
            if(fPhiTilde < 0){fPhiTilde = 2*fPi + fPhiTilde;}
            if(fPhiTilde > fPi){fPhiTilde = 2*fPi - fPhiTilde;}
            fHistRecPhiTildeCS[indexPt] -> Fill(fPhiTilde);
            // 2D approach
            fHistRecCosThetaPhiHE[indexPt] -> Fill(fCosThetaHERec[j],TMath::Abs(fPhiHERec[j]));
            fHistRecCosThetaPhiCS[indexPt] -> Fill(fCosThetaCSRec[j],TMath::Abs(fPhiCSRec[j]));
            indexPt = 0;
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

    fHistAccxEffPhiTildeHE[i] = new TH1D(Form("histAccxEffPhiTildeHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",10,0.,fPi);
    fHistAccxEffPhiTildeHE[i] -> Divide(fHistRecPhiTildeHE[i],fHistGenPhiTildeHE[i],1,1,"B");

    fHistAccxEffCosThetaPhiHE[i] = new TH2D(Form("histAccxEffCosThetaPhiHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]);
    fHistAccxEffCosThetaPhiHE[i] -> Divide(fHistRecCosThetaPhiHE[i],fHistGenCosThetaPhiHE[i],1,1,"B");

    fHistAccxEffCosThetaPhiStatRelHE[i] = new TH2D(Form("histAccxEffCosThetaPhiStatRelHE_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]);

    // COLLINS-SOPER
    fHistAccxEffCosThetaCS[i] = new TH1D(Form("histAccxEffCosThetaCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]);
    fHistAccxEffCosThetaCS[i] -> Divide(fHistRecCosThetaCS[i],fHistGenCosThetaCS[i],1,1,"B");

    fHistAccxEffPhiCS[i] = new TH1D(Form("histAccxEffPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistAccxEffPhiCS[i] -> Divide(fHistRecPhiCS[i],fHistGenPhiCS[i],1,1,"B");

    fHistAccxEffPhiTildeCS[i] = new TH1D(Form("histAccxEffPhiTildeCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",10,0.,fPi);
    fHistAccxEffPhiTildeCS[i] -> Divide(fHistRecPhiTildeCS[i],fHistGenPhiTildeCS[i],1,1,"B");

    fHistAccxEffCosThetaPhiCS[i] = new TH2D(Form("histAccxEffCosThetaPhiCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]);
    fHistAccxEffCosThetaPhiCS[i] -> Divide(fHistRecCosThetaPhiCS[i],fHistGenCosThetaPhiCS[i],1,1,"B");

    fHistAccxEffCosThetaPhiStatRelCS[i] = new TH2D(Form("histAccxEffCosThetaPhiStatRelCS_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0],fNPhiBins,&fPhiValues[0]);
  }

  fHistAccxEffCosThetaHEPt = new TH2D("histAccxEffCosThetaHEPt","",100,-1,1,100,0,10);
  fHistAccxEffCosThetaHEPt -> Divide(fHistRecCosThetaHEPt,fHistGenCosThetaHEPt,1,1,"B");
  fHistAccxEffPhiHEPt = new TH2D("histAccxEffPhiHEPt","",100,0,fPi,100,0,10);
  fHistAccxEffPhiHEPt -> Divide(fHistRecPhiHEPt,fHistGenPhiHEPt,1,1,"B");

  fHistAccxEffCosThetaCSPt = new TH2D("histAccxEffCosThetaCSPt","",100,-1,1,100,0,10);
  fHistAccxEffCosThetaCSPt -> Divide(fHistRecCosThetaCSPt,fHistGenCosThetaCSPt,1,1,"B");
  fHistAccxEffPhiCSPt = new TH2D("histAccxEffPhiCSPt","",100,0,fPi,100,0,10);
  fHistAccxEffPhiCSPt -> Divide(fHistRecPhiCSPt,fHistGenPhiCSPt,1,1,"B");

  for(int i = 0;i < fNPtBins;i++){
    for(int j = 0;j < fNCosThetaBins;j++){
      for(int k = 0;k < fNPhiBins;k++){
        fHistAccxEffCosThetaPhiStatRelHE[i] -> SetBinContent(j+1,k+1,(fHistAccxEffCosThetaPhiHE[i] -> GetBinError(j+1,k+1)/fHistAccxEffCosThetaPhiHE[i] -> GetBinContent(j+1,k+1))*100);
        fHistAccxEffCosThetaPhiStatRelCS[i] -> SetBinContent(j+1,k+1,(fHistAccxEffCosThetaPhiCS[i] -> GetBinError(j+1,k+1)/fHistAccxEffCosThetaPhiCS[i] -> GetBinContent(j+1,k+1))*100);
      }
    }
  }

  TFile *fileAccxEff = new TFile(nameOutputFile.c_str(),"RECREATE");
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
  }
  fHistGenCosThetaHEPt -> Write(); fHistGenPhiHEPt -> Write(); fHistAccxEffCosThetaHEPt -> Write();
  fHistRecCosThetaHEPt -> Write(); fHistRecPhiHEPt -> Write(); fHistAccxEffPhiHEPt -> Write();

  fHistGenCosThetaCSPt -> Write(); fHistGenPhiCSPt -> Write(); fHistAccxEffCosThetaCSPt -> Write();
  fHistRecCosThetaCSPt -> Write(); fHistRecPhiCSPt -> Write(); fHistAccxEffPhiCSPt -> Write();
  fileAccxEff -> Close();
}
//______________________________________________________________________________
void AccxEffCalculator::ReWeightAccxEff(string refFrame, Double_t LambdaTheta,Double_t LambdaPhi, string strSample, Bool_t saveFile, string nameOutputFile) {
  TF1 *funcCosTheta = new TF1("funcCosTheta","(1/(3 + [0]))*(1 + [0]*x*x)",-1,1);
  funcCosTheta -> SetParameter(0,LambdaTheta);

  TF1 *funcPhi = new TF1("funcPhi","(1 + ((2*[1])/(3 + [0]))*cos(2*x))",0,fPi);
  funcPhi -> SetParameter(0,LambdaTheta);
  funcPhi -> SetParameter(1,LambdaPhi);

  TF2 *funcCosThetaPhi = new TF2("funcCosThetaPhi","(1/(3 + [0]))*(1 + [0]*x*x*(1 - cos(2*y)) + [1]*cos(2*y))",-1,1,0,fPi);
  funcCosThetaPhi -> SetParameter(0,LambdaTheta);
  funcCosThetaPhi -> SetParameter(1,LambdaPhi);

  int nEvents = 0;

  if(strSample == "FullStat"){nEvents = fTreeAccxEff -> GetEntries();}
  if(strSample == "TestStat"){nEvents = 100000;}
  printf("N events = %i \n",nEvents);

  double weightCosTheta = 0;
  double weightPhi = 0;
  double weightCosThetaPhi = 0;

  for(int i = 0;i < fNPtBins;i++){
    fHistGenCosThetaReWeighted[i] = new TH1D(Form("histGenCosThetaReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]);
    fHistGenCosThetaReWeighted[i] -> Sumw2();

    fHistRecCosThetaReWeighted[i] = new TH1D(Form("histRecCosThetaReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]);
    fHistRecCosThetaReWeighted[i] -> Sumw2();

    fHistGenPhiReWeighted[i] = new TH1D(Form("histGenPhiReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistGenPhiReWeighted[i] -> Sumw2();

    fHistRecPhiReWeighted[i] = new TH1D(Form("histRecPhiReWeighted_%ipT%i",(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistRecPhiReWeighted[i] -> Sumw2();
  }

  for(int i = 0;i < nEvents;i++){
    printf("Reading : %3.2f %% \r",((double) i/(double) nEvents)*100.);
    fTreeAccxEff -> GetEntry(i);

    for(int j = 0;j < fNDimuGen;j++){
      if(fDimuPtGen[j] > fMinPt[0] && fDimuPtGen[j] < fMaxPt[0]){
        if(fDimuYGen[j] > -4. && fDimuYGen[j] < -2.5){
            if(refFrame == "HE"){
              weightCosTheta = (funcCosTheta -> Eval(fCosThetaHEGen[j]))/(funcCosTheta -> GetMaximum());
              weightPhi = (funcPhi -> Eval(fPhiHEGen[j]))/(funcPhi -> GetMaximum());
            }
            if(refFrame == "CS"){
              weightCosTheta = (funcCosTheta -> Eval(fCosThetaCSGen[j]))/(funcCosTheta -> GetMaximum());
              weightPhi = (funcPhi -> Eval(fPhiCSGen[j]))/(funcPhi -> GetMaximum());
            }
            fHistGenCosThetaReWeighted[0] -> Fill(fCosThetaHEGen[j],weightCosTheta);
            fHistGenPhiReWeighted[0] -> Fill(TMath::Abs(fPhiHEGen[j]),weightPhi);
        }
      }
    }

    for(int j = 0;j < fNDimuRec;j++){
      if(fDimuPtRec[j] > fMinPt[0] && fDimuPtRec[j] < fMaxPt[0]){
        if(fDimuYRec[j] > -4. && fDimuYRec[j] < -2.5){
          if(fDimuMatchRec[j] == 2){
            if(fDimuMassRec[j] > 2 && fDimuMassRec[j] < 5){
              if(refFrame == "HE"){
                fHistRecCosThetaReWeighted[0] -> Fill(fCosThetaHERec[j],weightCosTheta);
                fHistRecPhiReWeighted[0] -> Fill(TMath::Abs(fPhiHERec[j]),weightPhi);
              }
              if(refFrame == "CS"){
                fHistRecCosThetaReWeighted[0] -> Fill(fCosThetaCSRec[j],weightCosTheta);
                fHistRecPhiReWeighted[0] -> Fill(TMath::Abs(fPhiCSRec[j]),weightPhi);
              }
            }
          }
        }
      }
    }
  }
  printf("\n");

  for(int i = 0;i < fNPtBins;i++){
    fHistAccxEffCosThetaReWeighted[i] = new TH1D(Form("histAccxEffCosThetaReWeighted%s_%ipT%i",refFrame.c_str(),(int) fMinPt[i],(int) fMaxPt[i]),"",fNCosThetaBins,&fCosThetaValues[0]);
    fHistAccxEffCosThetaReWeighted[i] -> Divide(fHistRecCosThetaReWeighted[i],fHistGenCosThetaReWeighted[i],1,1,"B");

    fHistAccxEffPhiReWeighted[i] = new TH1D(Form("histAccxEffPhiReWeighted%s_%ipT%i",refFrame.c_str(),(int) fMinPt[i],(int) fMaxPt[i]),"",fNPhiBins,&fPhiValues[0]);
    fHistAccxEffPhiReWeighted[i] -> Divide(fHistRecPhiReWeighted[i],fHistGenPhiReWeighted[i],1,1,"B");
  }

  if(saveFile){
    TFile *fileAccxEffReWeighted = new TFile(nameOutputFile.c_str(),"RECREATE");
    for(int i = 0;i < fNPtBins;i++){
      fHistGenCosThetaReWeighted[i] -> Write();
      fHistRecCosThetaReWeighted[i] -> Write();
      fHistAccxEffCosThetaReWeighted[i] -> Write();

      fHistGenPhiReWeighted[i] -> Write();
      fHistRecPhiReWeighted[i] -> Write();
      fHistAccxEffPhiReWeighted[i] -> Write();
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
