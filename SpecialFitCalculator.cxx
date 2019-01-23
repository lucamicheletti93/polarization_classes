#include <Riostream.h>
#include <string>
#include <vector>

#include "TROOT.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TMinuit.h"
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
#include "TLatex.h"
#include <vector>
#include "TMatrixD.h"

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TRandom.h"

#include "SpecialFitCalculator.h"

// Global variables
Int_t gNDistrib;
TH1D *gHistFit[3];
TF1 *gFuncFit[3];

Int_t gNCosThetaBins;
Int_t gNPhiBins;
Int_t gNPhiTildeBins;
vector <Double_t> gCosThetaValues;
vector <Double_t> gPhiValues;
vector <Double_t> gPhiTildeValues;

Double_t funcCosTheta(Double_t , Double_t *);
Double_t funcPhi(Double_t , Double_t *);
Double_t funcPhiTilde(Double_t , Double_t *);
void polarizationFCN(Int_t &, Double_t *, Double_t &, Double_t *, Int_t );

ClassImp(SpecialFitCalculator)

//______________________________________________________________________________
SpecialFitCalculator::SpecialFitCalculator(): TObject() {
  gPi = TMath::Pi();
  // default constructor
}
//______________________________________________________________________________
SpecialFitCalculator::~SpecialFitCalculator() {
  // destructor
}
void SpecialFitCalculator::SetBinning(vector <Double_t> CosThetaValues, vector <Double_t> PhiValues, vector <Double_t> PhiTildeValues) {
  for(int i = 0;i < (int) CosThetaValues.size();i++){gCosThetaValues.push_back(CosThetaValues[i]);}
  gNCosThetaBins = gCosThetaValues.size() - 1;
  for(int i = 0;i < (int) PhiValues.size();i++){gPhiValues.push_back(PhiValues[i]);}
  gNPhiBins = gPhiValues.size() - 1;
  for(int i = 0;i < (int) PhiTildeValues.size();i++){gPhiTildeValues.push_back(PhiTildeValues[i]);}
  gNPhiTildeBins = gPhiTildeValues.size() - 1;
}
//______________________________________________________________________________
void SpecialFitCalculator::SimultaneousFit(TObjArray *data) {
  Int_t ndf = 0;
  gNDistrib = data -> GetEntries();
  cout << "n distributions = " << gNDistrib << endl;

  for(int i = 0;i < gNDistrib;i++){
    gHistFit[i] = (TH1D*) data -> At(i);
    gHistFit[i] -> SetName(Form("Distrib%i",i));
    ndf += gHistFit[i] -> GetSize();
    cout << i << ") ndf = " << gHistFit[i] -> GetSize() << endl;
  }
  cout << "ndf = " << ndf << endl;

  TCanvas *canvasHistFit = new TCanvas("canvasHistFit");
  canvasHistFit -> Divide(gNDistrib,1);
  for(int i = 0;i < gNDistrib;i++){
    canvasHistFit -> cd(i+1);
    gHistFit[i] -> Draw("PE");
  }
  canvasHistFit -> Update();

  TMinuit *minuit = new TMinuit(6);
  minuit -> SetFCN(polarizationFCN);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);

  minuit -> mnparm(0,"normCosTheta",1.e7,1.,0.,0.,ierflg);
  minuit -> mnparm(1,"normPhi",1.e7,1.,0.,0.,ierflg);
  minuit -> mnparm(2,"normPhiTilde",1.e7,1.,0.,0.,ierflg);
  minuit -> mnparm(3,"lambdaTheta",0.,0.001,-1.,1.,ierflg);
  minuit -> mnparm(4,"lambdaPhi",0.,0.001,-1.,1.,ierflg);
  minuit -> mnparm(5,"lambdaThetaPhi",0.,0.001,-1.,1.,ierflg);

  arglist[0] = 500;
  arglist[1] = 1.;
  minuit -> mnexcm("MIGRAD",arglist,2,ierflg);

  arglist[0] = 500000;
  minuit -> mnexcm("IMPROVE",arglist,1,ierflg);

  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(6,amin);

  Double_t normCosTheta, errNormCosTheta;
  Double_t normPhi, errNormPhi;
  Double_t normPhiTilde, errNormPhiTilde;
  Double_t lambdaTheta, errLambdaTheta;
  Double_t lambdaPhi, errLambdaPhi;
  Double_t lambdaThetaPhi, errLambdaThetaPhi;

  minuit -> GetParameter(0,normCosTheta,errNormCosTheta);
  minuit -> GetParameter(1,normPhi,errNormPhi);
  minuit -> GetParameter(2,normPhiTilde,errNormPhiTilde);
  minuit -> GetParameter(3,lambdaTheta,errLambdaTheta);
  minuit -> GetParameter(4,lambdaPhi,errLambdaPhi);
  minuit -> GetParameter(5,lambdaThetaPhi,errLambdaThetaPhi);

  gFuncFit[0] = new TF1("gFuncFit0","([0]/(3 + [1]))*(1 + [1]*x*x)",-1.,1.);
  gFuncFit[0] -> SetParameter(0,normCosTheta);
  gFuncFit[0] -> SetParameter(1,lambdaTheta);

  gFuncFit[1] = new TF1("gFuncFit1","[0]*(1 + ((2*[2])/(3 + [1]))*cos(2*x))",0.,gPi);
  gFuncFit[1] -> SetParameter(0,normPhi);
  gFuncFit[1] -> SetParameter(1,lambdaTheta);
  gFuncFit[1] -> SetParameter(2,lambdaPhi);

  gFuncFit[2] = new TF1("gFuncFit2","[0]*(1 + ((sqrt(2)*[2])/(3 + [1]))*cos(x))",0.,2*gPi);
  gFuncFit[2] -> SetParameter(0,normPhiTilde);
  gFuncFit[2] -> SetParameter(1,lambdaTheta);
  gFuncFit[2] -> SetParameter(2,lambdaThetaPhi);

  for(int i = 0;i < gNDistrib;i++){
    canvasHistFit -> cd(i+1);
    gFuncFit[i] -> Draw("same");
  }
  canvasHistFit -> Update();
}
//______________________________________________________________________________
Double_t funcCosTheta(Double_t x, Double_t *par){
  return (par[0]/(3 + par[1]))*(1 + par[1]*x*x);
}
//______________________________________________________________________________
Double_t funcPhi(Double_t x, Double_t *par){
  return par[0]*(1 + ((2*par[2])/(3 + par[1]))*TMath::Cos(2*x));
}
//______________________________________________________________________________
Double_t funcPhiTilde(Double_t x, Double_t *par){
  return par[0]*(1 + ((TMath::Sqrt(2)*par[2])/(3 + par[1]))*TMath::Cos(x));
}
//______________________________________________________________________________
void polarizationFCN(Int_t &npar, Double_t *gin, Double_t &gChiSquare, Double_t *par, Int_t iflag){
  Double_t angleVar;
  Double_t val, errVal;

  //cout << normCosTheta << " " << lambdaTheta << endl;
  //cout << par[0] << " " << par[1] << " " << par[2] << " " << par[3] << endl;

  Double_t parCosTheta[2];
  parCosTheta[0] = par[0];
  parCosTheta[1] = par[3];

  Double_t parPhi[3];
  parPhi[0] = par[1];
  parPhi[1] = par[3];
  parPhi[2] = par[4];

  Double_t parPhiTilde[3];
  parPhiTilde[0] = par[2];
  parPhiTilde[1] = par[3];
  parPhiTilde[2] = par[5];

  Double_t func = 0, pull = 0, chiSquare = 0;

  for(int i = 0;i < gNDistrib;i++){
    for(int j = 0;j < gHistFit[i] -> GetSize() - 2;j++){
      val = gHistFit[i] -> GetBinContent(j+1);
      errVal = gHistFit[i] -> GetBinError(j+1);

      if(val == 0. && errVal == 0.){continue;}
      if(i == 0){
        angleVar = gCosThetaValues[j];
        pull = (val - funcCosTheta(angleVar,parCosTheta))/errVal;
      }
      if(i == 1){
        angleVar = gPhiValues[j];
        pull = (val - funcPhi(angleVar,parPhi))/errVal;
      }
      if(i == 2){
        angleVar = gPhiTildeValues[j];
        pull = (val - funcPhiTilde(angleVar,parPhiTilde))/errVal;
      }
      chiSquare += pull*pull;
    }
  }

  gChiSquare = chiSquare;
}
//______________________________________________________________________________
vector <Double_t> SpecialFitCalculator::GetCosThetaParametersList(){
  return fCosThetaParametersList;
}
//______________________________________________________________________________
vector <Double_t> SpecialFitCalculator::GetErrorCosThetaParametersList(){
  return fErrorCosThetaParametersList;
}
//______________________________________________________________________________
vector <Double_t> SpecialFitCalculator::GetPhiParametersList(){
  return fPhiParametersList;
}
//______________________________________________________________________________
vector <Double_t> SpecialFitCalculator::GetErrorPhiParametersList(){
  return fErrorPhiParametersList;
}
