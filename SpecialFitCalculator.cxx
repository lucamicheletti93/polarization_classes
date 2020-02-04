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
#include "TGraph.h"
#include <vector>
#include "TMatrixD.h"
#include "TPaveText.h"

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

TH1D *gHistBarbFit[6];
TF1 *gFuncBarbFit[6];

Int_t gNCosThetaBins;
Int_t gNPhiBins;
Int_t gNPhiTildeBins;
vector <Double_t> gCosThetaValues;
vector <Double_t> gPhiValues;
vector <Double_t> gPhiTildeValues;

Int_t    gMinFitRangeCosTheta;
Int_t    gMinFitRangePhi;
Int_t    gMinFitRangePhiTilde;
Int_t    gMaxFitRangeCosTheta;
Int_t    gMaxFitRangePhi;
Int_t    gMaxFitRangePhiTilde;

string   gPlotTitle;

Double_t globalChiSquare;

Double_t funcCosTheta(Double_t , Double_t *);
Double_t funcPhi(Double_t , Double_t *);
Double_t funcPhiTilde(Double_t , Double_t *);
void polarizationFCN(Int_t &, Double_t *, Double_t &, Double_t *, Int_t );
void polarizationBarbatruccoFCN(Int_t &, Double_t *, Double_t &, Double_t *, Int_t );

ClassImp(SpecialFitCalculator)

//______________________________________________________________________________
SpecialFitCalculator::SpecialFitCalculator(): TObject() {
  gPi = TMath::Pi();
  fInitNorm = kFALSE;
  // default constructor
}
//______________________________________________________________________________
SpecialFitCalculator::~SpecialFitCalculator() {
  // destructor
}
//______________________________________________________________________________
void SpecialFitCalculator::SetBinning(vector <Double_t> CosThetaValues, vector <Double_t> PhiValues, vector <Double_t> PhiTildeValues) {
  for(int i = 0;i < (int) CosThetaValues.size();i++){gCosThetaValues.push_back(CosThetaValues[i]);}
  gNCosThetaBins = gCosThetaValues.size() - 1;
  for(int i = 0;i < (int) PhiValues.size();i++){gPhiValues.push_back(PhiValues[i]);}
  gNPhiBins = gPhiValues.size() - 1;
  for(int i = 0;i < (int) PhiTildeValues.size();i++){gPhiTildeValues.push_back(PhiTildeValues[i]);}
  gNPhiTildeBins = gPhiTildeValues.size() - 1;
}
//______________________________________________________________________________
void SpecialFitCalculator::SetFitRange(Int_t minFitRange[], Int_t maxFitRange[]) {
  gMinFitRangeCosTheta = minFitRange[0];
  gMinFitRangePhi = minFitRange[1];
  gMinFitRangePhiTilde = minFitRange[2];
  gMaxFitRangeCosTheta = maxFitRange[0];
  gMaxFitRangePhi = maxFitRange[1];
  gMaxFitRangePhiTilde = maxFitRange[2];
}
//______________________________________________________________________________
void SpecialFitCalculator::SetFitNormalization(Double_t NormCosTheta, Double_t NormPhi, Double_t NormPhiTilde){
  fInitNorm = kTRUE;
  fInitNormCosTheta = NormCosTheta;
  fInitNormPhi = NormPhi;
  fInitNormPhiTilde = NormPhiTilde;
}
//______________________________________________________________________________
void SpecialFitCalculator::SetPlotTitle(string plotTitle){
   gPlotTitle = plotTitle;
}
//______________________________________________________________________________
void SpecialFitCalculator::SimultaneousFit(TObjArray *data, Bool_t saveCanvas, string nameCanvas) {
  printf("\n------------------------------------------------------------------------\n");
  printf("This object allows you to compute simultaneous fit with 3 distributions \n");
  printf("If the class doesn't work REMEMBER to use the SetBinning() method! \n");
  printf("------------------------------------------------------------------------\n");

  Int_t ndf = 0;
  gNDistrib = data -> GetEntries();
  cout << "n distributions = " << gNDistrib << endl;

  for(int i = 0;i < gNDistrib;i++){
    gHistFit[i] = (TH1D*) data -> At(i); gHistFit[i] -> SetMarkerStyle(20); gHistFit[i] -> SetMarkerSize(0.8);
    gHistFit[i] -> SetName(Form("Distrib%i",i));
    ndf += gHistFit[i] -> GetSize();
    cout << i << ") ndf = " << gHistFit[i] -> GetSize() << endl;
  }
  cout << "ndf = " << ndf << endl;
  gHistFit[0] -> SetAxisRange(-0.7,0.7);

  TMinuit *minuit = new TMinuit(6);
  minuit -> SetFCN(polarizationFCN);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);

  if(fInitNorm){
    minuit -> mnparm(0,"normCosTheta",fInitNormCosTheta,1.,0.,0.,ierflg);
    minuit -> mnparm(1,"normPhi",fInitNormPhi,1.,0.,0.,ierflg);
    minuit -> mnparm(2,"normPhiTilde",fInitNormPhiTilde,1.,0.,0.,ierflg);
  }
  else{
    minuit -> mnparm(0,"normCosTheta",1.e7,1.,0.,0.,ierflg);
    minuit -> mnparm(1,"normPhi",1.e7,1.,0.,0.,ierflg);
    minuit -> mnparm(2,"normPhiTilde",1.e7,1.,0.,0.,ierflg);
  }
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

  minuit -> GetParameter(0,normCosTheta,errNormCosTheta);
  minuit -> GetParameter(1,normPhi,errNormPhi);
  minuit -> GetParameter(2,normPhiTilde,errNormPhiTilde);
  minuit -> GetParameter(3,fLambdaTheta,fErrorLambdaTheta);
  minuit -> GetParameter(4,fLambdaPhi,fErrorLambdaPhi);
  minuit -> GetParameter(5,fLambdaThetaPhi,fErrorLambdaThetaPhi);

  fGraContour_lambdaTheta_lambdaPhi = (TGraph*) minuit -> Contour(40,3,4);

  gFuncFit[0] = new TF1("gFuncFit0","([0]/(3 + [1]))*(1 + [1]*x*x)",-1.,1.);
  gFuncFit[0] -> SetParameter(0,normCosTheta);
  gFuncFit[0] -> SetParameter(1,fLambdaTheta);

  gFuncFit[1] = new TF1("gFuncFit1","[0]*(1 + ((2*[2])/(3 + [1]))*cos(2*x))",0.,gPi);
  gFuncFit[1] -> SetParameter(0,normPhi);
  gFuncFit[1] -> SetParameter(1,fLambdaTheta);
  gFuncFit[1] -> SetParameter(2,fLambdaPhi);

  gFuncFit[2] = new TF1("gFuncFit2","[0]*(1 + ((sqrt(2)*[2])/(3 + [1]))*cos(x))",0.,2*gPi);
  gFuncFit[2] -> SetParameter(0,normPhiTilde);
  gFuncFit[2] -> SetParameter(1,fLambdaTheta);
  gFuncFit[2] -> SetParameter(2,fLambdaThetaPhi);

  /*
  TH2D *histGridContour = new TH2D("histGridContour","histGridContour",100,-1.,1.,100,-1.,1.);
  TCanvas *canvasContour = new TCanvas("canvasContour","canvasContour",600,600);
  histGridContour -> Draw();
  gr34 -> Draw("alp");
  */

  TH2D *histGridCosTheta = new TH2D("histGridCosTheta","",100,0.,1.,100.,0.,gHistFit[0] -> GetMaximum() + 0.5*gHistFit[0] -> GetMaximum());
  histGridCosTheta -> GetXaxis() -> SetTitle("cos#theta");

  TH2D *histGridPhi = new TH2D("histGridPhi","",100,0.,gPi,100.,0.,gHistFit[1] -> GetMaximum() + 0.3*gHistFit[1] -> GetMaximum());
  histGridPhi -> GetXaxis() -> SetTitle("#varphi");

  TH2D *histGridPhiTilde = new TH2D("histGridPhiTilde","",100,0.,2*gPi,100.,0.,gHistFit[2] -> GetMaximum() + 0.3*gHistFit[2] -> GetMaximum());
  histGridPhiTilde -> GetXaxis() -> SetTitle("#tilde{#varphi}");

  gHistFit[0] -> SetMarkerStyle(20); gHistFit[0] -> SetMarkerColor(kBlack); gHistFit[0] -> SetLineColor(kBlack);
  gHistFit[1] -> SetMarkerStyle(20); gHistFit[1] -> SetMarkerColor(kBlack); gHistFit[1] -> SetLineColor(kBlack);
  gHistFit[2] -> SetMarkerStyle(20); gHistFit[2] -> SetMarkerColor(kBlack); gHistFit[2] -> SetLineColor(kBlack);

  TPaveText *paveTextPolCosTheta = new TPaveText(0.2,0.85,0.8,0.98,"NDC");
  paveTextPolCosTheta -> SetFillColor(kWhite);
  paveTextPolCosTheta -> AddText(gPlotTitle.c_str());
  paveTextPolCosTheta -> AddText(Form("#lambda_{#theta} = %3.2f #pm %3.2f",fLambdaTheta,fErrorLambdaTheta));

  TPaveText *paveTextPolPhi = new TPaveText(0.2,0.85,0.8,0.98,"NDC");
  paveTextPolPhi -> SetFillColor(kWhite);
  paveTextPolPhi -> AddText(gPlotTitle.c_str());
  paveTextPolPhi -> AddText(Form("#lambda_{#varphi} = %3.2f #pm %3.2f",fLambdaPhi,fErrorLambdaPhi));

  TPaveText *paveTextPolPhiTilde = new TPaveText(0.2,0.85,0.8,0.98,"NDC");
  paveTextPolPhiTilde -> SetFillColor(kWhite);
  paveTextPolPhiTilde -> AddText(gPlotTitle.c_str());
  paveTextPolPhiTilde -> AddText(Form("#lambda_{#theta#varphi} = %3.2f #pm %3.2f",fLambdaThetaPhi,fErrorLambdaThetaPhi));

  TCanvas *canvasHistFitSim_CosTheta = new TCanvas("canvasHistFitSim_CosTheta","canvasHistFitSim_CosTheta",600,600);
  histGridCosTheta -> Draw(); gHistFit[0] -> Draw("EPsame"); gFuncFit[0] -> Draw("same");
  paveTextPolCosTheta -> Draw();
  //latexTitle -> DrawLatex(0.05,27000.,gPlotTitle.c_str());
  //latexTitle -> DrawLatex(0.3,25000.,Form("#lambda_{#theta} = %3.2f #pm %3.2f",fLambdaTheta,fErrorLambdaTheta));
  //latexTitle -> DrawLatex(0.2,3000.,Form("#chi^{2}/NDF = %3.2f",chiSquare/NDF));

  TCanvas *canvasHistFitSim_Phi = new TCanvas("canvasHistFitSim_Phi","canvasHistFitSim_Phi",600,600);
  histGridPhi -> Draw(); gHistFit[1] -> Draw("EPsame"); gFuncFit[1] -> Draw("same");
  paveTextPolPhi -> Draw();
  //latexTitle -> DrawLatex(0.2,9000.,gPlotTitle.c_str());
  //latexTitle -> DrawLatex(1.,8000.,Form("#lambda_{#varphi} = %3.2f #pm %3.2f",fLambdaPhi,fErrorLambdaPhi));
  //latexTitle -> DrawLatex(1,1000.,Form("#chi^{2}/NDF = %3.2f",chiSquare/NDF));

  TCanvas *canvasHistFitSim_PhiTilde = new TCanvas("canvasHistFitSim_PhiTilde","canvasHistFitSim_PhiTilde",600,600);
  histGridPhiTilde -> Draw(); gHistFit[2] -> Draw("EPsame"); gFuncFit[2] -> Draw("same");
  paveTextPolPhiTilde -> Draw();
  //latexTitle -> DrawLatex(0.1,4800.,gPlotTitle.c_str());
  //latexTitle -> DrawLatex(2.,4000.,Form("#lambda_{#theta#varphi} = %3.2f #pm %3.2f",fLambdaThetaPhi,fErrorLambdaThetaPhi));
  //latexTitle -> DrawLatex(2,500.,Form("#chi^{2}/NDF = %3.2f",chiSquare/NDF));

  if(saveCanvas){
    //canvasHistFit -> SaveAs(Form("%s.png",nameCanvas.c_str()));
    //canvasContour -> SaveAs(Form("%s_contour.png",nameCanvas.c_str()));
    canvasHistFitSim_CosTheta -> SaveAs(Form("%s_CosTheta.pdf",nameCanvas.c_str()));
    canvasHistFitSim_Phi -> SaveAs(Form("%s_Phi.pdf",nameCanvas.c_str()));
    canvasHistFitSim_PhiTilde -> SaveAs(Form("%s_PhiTilde.pdf",nameCanvas.c_str()));
  }

  delete canvasHistFitSim_CosTheta;
  delete canvasHistFitSim_Phi;
  delete canvasHistFitSim_PhiTilde;
}
//______________________________________________________________________________
void SpecialFitCalculator::BarbatruccoFit(TObjArray *data, Bool_t saveCanvas, string nameCanvas) {
  printf("\n------------------------------------------------------------------------\n");
  printf("This object allows you to compute barbatrucco fit with  distributions \n");
  printf("If the class doesn't work REMEMBER to use the SetBinning() method! \n");
  printf("------------------------------------------------------------------------\n");

  Int_t ndf = 0;
  gNDistrib = data -> GetEntries();
  cout << "n distributions = " << gNDistrib << endl;

  for(int i = 0;i < gNDistrib;i++){
    gHistBarbFit[i] = (TH1D*) data -> At(i); gHistBarbFit[i] -> SetMarkerStyle(20); gHistBarbFit[i] -> SetMarkerSize(0.8);
    gHistBarbFit[i] -> SetName(Form("Distrib%i",i));
    ndf += gHistBarbFit[i] -> GetSize();
    cout << i << ") ndf = " << gHistBarbFit[i] -> GetSize() << endl;
  }
  cout << "ndf = " << ndf << endl;

  TCanvas *canvasHistFit = new TCanvas("canvasHistFit");
  canvasHistFit -> Divide(gNDistrib/2,2);
  for(int i = 0;i < gNDistrib;i++){
    canvasHistFit -> cd(i+1);
    gHistBarbFit[i] -> Draw("PE");
  }
  canvasHistFit -> Update();

  TMinuit *minuit = new TMinuit(11);
  minuit -> SetFCN(polarizationBarbatruccoFCN);

  Double_t arglist[10];
  Int_t ierflg = 0;

  arglist[0] = 1;
  minuit -> mnexcm("SET ERR",arglist,1,ierflg);

  /*
  minuit -> mnparm(0,"normCosTheta",1.e7,1.,0.,0.,ierflg);
  minuit -> mnparm(1,"normPhi",1.e7,1.,0.,0.,ierflg);
  minuit -> mnparm(2,"normPhiTilde",1.e7,1.,0.,0.,ierflg);
  minuit -> mnparm(3,"lambdaThetaHE",0.,0.001,-1.,1.,ierflg);
  minuit -> mnparm(4,"lambdaThetaPhiHE",0.,0.001,-1.,1.,ierflg);
  minuit -> mnparm(8,"lambdaThetaCS",0.,0.001,-1.,1.,ierflg);
  minuit -> mnparm(9,"lambdaThetaPhiCS",0.,0.001,-1.,1.,ierflg);
  minuit -> mnparm(10,"lambdaTilde",0.,0.001,-0.5,0.5,ierflg);
  */


  minuit -> mnparm(0,"normCosThetaHE",1.e7,1.,0.,0.,ierflg);
  minuit -> mnparm(1,"normPhiHE",1.e7,1.,0.,0.,ierflg);
  minuit -> mnparm(2,"normPhiTildeHE",1.e7,1.,0.,0.,ierflg);
  minuit -> mnparm(3,"lambdaThetaHE",0.,0.001,-1.,1.,ierflg);
  minuit -> mnparm(4,"lambdaThetaPhiHE",0.,0.001,-1.,1.,ierflg);
  minuit -> mnparm(5,"normCosThetaCS",1.e7,1.,0.,0.,ierflg);
  minuit -> mnparm(6,"normPhiCS",1.e7,1.,0.,0.,ierflg);
  minuit -> mnparm(7,"normPhiTildeCS",1.e7,1.,0.,0.,ierflg);
  minuit -> mnparm(8,"lambdaThetaCS",0.,0.001,-1.,1.,ierflg);
  minuit -> mnparm(9,"lambdaThetaPhiCS",0.,0.001,-1.,1.,ierflg);
  minuit -> mnparm(10,"lambdaTilde",0.,0.001,-0.5,0.5,ierflg);


  arglist[0] = 500;
  arglist[1] = 1.;
  minuit -> mnexcm("MIGRAD",arglist,2,ierflg);

  arglist[0] = 500000;
  minuit -> mnexcm("IMPROVE",arglist,1,ierflg);

  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  minuit -> mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  minuit -> mnprin(11,amin);

  /*
  Double_t normCosTheta, errNormCosTheta;
  Double_t normPhi, errNormPhi;
  Double_t normPhiTilde, errNormPhiTilde;
  */

  Double_t normCosThetaHE, errNormCosThetaHE;
  Double_t normPhiHE, errNormPhiHE;
  Double_t normPhiTildeHE, errNormPhiTildeHE;

  Double_t normCosThetaCS, errNormCosThetaCS;
  Double_t normPhiCS, errNormPhiCS;
  Double_t normPhiTildeCS, errNormPhiTildeCS;

  /*
  minuit -> GetParameter(0,normCosTheta,errNormCosTheta);
  minuit -> GetParameter(1,normPhi,errNormPhi);
  minuit -> GetParameter(2,normPhiTilde,errNormPhiTilde);
  minuit -> GetParameter(3,fLambdaThetaHE,fErrorLambdaThetaHE);
  minuit -> GetParameter(4,fLambdaThetaPhiHE,fErrorLambdaThetaPhiHE);
  minuit -> GetParameter(5,fLambdaThetaCS,fErrorLambdaThetaCS);
  minuit -> GetParameter(6,fLambdaThetaPhiCS,fErrorLambdaThetaPhiCS);
  minuit -> GetParameter(7,fLambdaTilde,fErrorLambdaTilde);
  */

  minuit -> GetParameter(0,normCosThetaHE,errNormCosThetaHE);
  minuit -> GetParameter(1,normPhiHE,errNormPhiHE);
  minuit -> GetParameter(2,normPhiTildeHE,errNormPhiTildeHE);
  minuit -> GetParameter(3,fLambdaThetaHE,fErrorLambdaThetaHE);
  minuit -> GetParameter(4,fLambdaThetaPhiHE,fErrorLambdaThetaPhiHE);
  minuit -> GetParameter(5,normCosThetaCS,errNormCosThetaCS);
  minuit -> GetParameter(6,normPhiCS,errNormPhiCS);
  minuit -> GetParameter(7,normPhiTildeCS,errNormPhiTildeCS);
  minuit -> GetParameter(8,fLambdaThetaCS,fErrorLambdaThetaCS);
  minuit -> GetParameter(9,fLambdaThetaPhiCS,fErrorLambdaThetaPhiCS);
  minuit -> GetParameter(10,fLambdaTilde,fErrorLambdaTilde);


  fLambdaPhiHE = (fLambdaTilde - fLambdaThetaHE)/(fLambdaTilde + 3);
  fLambdaPhiCS = (fLambdaTilde - fLambdaThetaCS)/(fLambdaTilde + 3);

  /*
  gFuncBarbFit[0] = new TF1("gFuncBarbFit0","([0]/(3 + [1]))*(1 + [1]*x*x)",-1.,1.);
  gFuncBarbFit[0] -> SetParameter(0,normCosTheta);
  gFuncBarbFit[0] -> SetParameter(1,fLambdaThetaHE);

  gFuncBarbFit[1] = new TF1("gFuncBarbFit1","[0]*(1 + ((2*[2])/(3 + [1]))*cos(2*x))",0.,gPi);
  gFuncBarbFit[1] -> SetParameter(0,normPhi);
  gFuncBarbFit[1] -> SetParameter(1,fLambdaThetaHE);
  gFuncBarbFit[1] -> SetParameter(2,fLambdaPhiHE);

  gFuncBarbFit[2] = new TF1("gFuncBarbFit2","[0]*(1 + ((sqrt(2)*[2])/(3 + [1]))*cos(x))",0.,2*gPi);
  gFuncBarbFit[2] -> SetParameter(0,normPhiTilde);
  gFuncBarbFit[2] -> SetParameter(1,fLambdaThetaHE);
  gFuncBarbFit[2] -> SetParameter(2,fLambdaThetaPhiHE);

  gFuncBarbFit[3] = new TF1("gFuncBarbFit3","([0]/(3 + [1]))*(1 + [1]*x*x)",-1.,1.);
  gFuncBarbFit[3] -> SetParameter(0,normCosTheta);
  gFuncBarbFit[3] -> SetParameter(1,fLambdaThetaCS);

  gFuncBarbFit[4] = new TF1("gFuncBarbFit4","[0]*(1 + ((2*[2])/(3 + [1]))*cos(2*x))",0.,gPi);
  gFuncBarbFit[4] -> SetParameter(0,normPhi);
  gFuncBarbFit[4] -> SetParameter(1,fLambdaThetaCS);
  gFuncBarbFit[4] -> SetParameter(2,fLambdaPhiCS);

  gFuncBarbFit[5] = new TF1("gFuncBarbFit5","[0]*(1 + ((sqrt(2)*[2])/(3 + [1]))*cos(x))",0.,2*gPi);
  gFuncBarbFit[5] -> SetParameter(0,normPhiTilde);
  gFuncBarbFit[5] -> SetParameter(1,fLambdaThetaCS);
  gFuncBarbFit[5] -> SetParameter(2,fLambdaThetaPhiCS);
  */


  gFuncBarbFit[0] = new TF1("gFuncBarbFit0","([0]/(3 + [1]))*(1 + [1]*x*x)",-1.,1.);
  gFuncBarbFit[0] -> SetParameter(0,normCosThetaHE);
  gFuncBarbFit[0] -> SetParameter(1,fLambdaThetaHE);

  gFuncBarbFit[1] = new TF1("gFuncBarbFit1","[0]*(1 + ((2*[2])/(3 + [1]))*cos(2*x))",0.,gPi);
  gFuncBarbFit[1] -> SetParameter(0,normPhiHE);
  gFuncBarbFit[1] -> SetParameter(1,fLambdaThetaHE);
  gFuncBarbFit[1] -> SetParameter(2,fLambdaPhiHE);

  gFuncBarbFit[2] = new TF1("gFuncBarbFit2","[0]*(1 + ((sqrt(2)*[2])/(3 + [1]))*cos(x))",0.,2*gPi);
  gFuncBarbFit[2] -> SetParameter(0,normPhiTildeHE);
  gFuncBarbFit[2] -> SetParameter(1,fLambdaThetaHE);
  gFuncBarbFit[2] -> SetParameter(2,fLambdaThetaPhiHE);

  gFuncBarbFit[3] = new TF1("gFuncBarbFit3","([0]/(3 + [1]))*(1 + [1]*x*x)",-1.,1.);
  gFuncBarbFit[3] -> SetParameter(0,normCosThetaCS);
  gFuncBarbFit[3] -> SetParameter(1,fLambdaThetaCS);

  gFuncBarbFit[4] = new TF1("gFuncBarbFit4","[0]*(1 + ((2*[2])/(3 + [1]))*cos(2*x))",0.,gPi);
  gFuncBarbFit[4] -> SetParameter(0,normPhiCS);
  gFuncBarbFit[4] -> SetParameter(1,fLambdaThetaCS);
  gFuncBarbFit[4] -> SetParameter(2,fLambdaPhiCS);

  gFuncBarbFit[5] = new TF1("gFuncBarbFit5","[0]*(1 + ((sqrt(2)*[2])/(3 + [1]))*cos(x))",0.,2*gPi);
  gFuncBarbFit[5] -> SetParameter(0,normPhiTildeCS);
  gFuncBarbFit[5] -> SetParameter(1,fLambdaThetaCS);
  gFuncBarbFit[5] -> SetParameter(2,fLambdaThetaPhiCS);


  for(int i = 0;i < gNDistrib;i++){
    canvasHistFit -> cd(i+1);
    gFuncBarbFit[i] -> Draw("same");
  }
  canvasHistFit -> Update();


  if(saveCanvas){
    canvasHistFit -> SaveAs(Form("%s.png",nameCanvas.c_str()));
  }

  delete canvasHistFit;
}
//______________________________________________________________________________
void SpecialFitCalculator::DecoupledFit(TObjArray *data, Bool_t saveCanvas, string nameCanvas, double minFitRange, double maxFitRange) {
  printf("\n------------------------------------------------------------------------\n");
  printf("This object allows you to compute decoupled fit with 3 distributions \n");
  printf("If the class doesn't work REMEMBER to use the SetBinning() method! \n");
  printf("------------------------------------------------------------------------\n");

  Int_t ndf = 0;
  gNDistrib = data -> GetEntries();
  cout << "n distributions = " << gNDistrib << endl;

  for(int i = 0;i < gNDistrib;i++){
    gHistFit[i] = (TH1D*) data -> At(i); gHistFit[i] -> SetMarkerStyle(20); gHistFit[i] -> SetMarkerSize(0.8); gHistFit[i] -> SetMarkerColor(kBlack); gHistFit[i] -> SetLineColor(kBlack);
    gHistFit[i] -> SetName(Form("Distrib%i",i));
    ndf += gHistFit[i] -> GetSize();
    cout << i << ") ndf = " << gHistFit[i] -> GetSize() << endl;
  }
  cout << "ndf = " << ndf << endl;

  gFuncFit[0] = new TF1("gFuncFit0","([0]/(3 + [1]))*(1 + [1]*x*x)",minFitRange,maxFitRange);
  gHistFit[0] -> Fit(gFuncFit[0],"R0IQ");

  gFuncFit[1] = new TF1("gFuncFit1","[0]*(1 + ((2*[2])/(3 + [1]))*cos(2*x))",0.,gPi);
  gFuncFit[1] -> FixParameter(1,gFuncFit[0] -> GetParameter(1));
  gHistFit[1] -> Fit(gFuncFit[1],"R0IQ");

  gFuncFit[2] = new TF1("gFuncFit2","[0]*(1 + ((sqrt(2)*[2])/(3 + [1]))*cos(x))",0.,2*gPi);
  gFuncFit[2] -> FixParameter(1,gFuncFit[0] -> GetParameter(1));
  gHistFit[2] -> Fit(gFuncFit[2],"R0IQ");

  fLambdaTheta = gFuncFit[0] -> GetParameter(1);
  fErrorLambdaTheta = gFuncFit[0] -> GetParError(1);
  fLambdaPhi = gFuncFit[1] -> GetParameter(2);
  fErrorLambdaPhi = gFuncFit[1] -> GetParError(2);
  fLambdaThetaPhi = gFuncFit[2] -> GetParameter(2);
  fErrorLambdaThetaPhi = gFuncFit[2] -> GetParError(2);

  TH2D *histGridCosTheta = new TH2D("histGridCosTheta","",100,minFitRange,maxFitRange,100.,0.,gHistFit[0] -> GetMaximum() + 0.5*gHistFit[0] -> GetMaximum());
  histGridCosTheta -> GetXaxis() -> SetTitle("cos#theta");

  TH2D *histGridPhi = new TH2D("histGridPhi","",100,0.,gPi,100.,0.,gHistFit[1] -> GetMaximum() + 0.3*gHistFit[1] -> GetMaximum());
  histGridPhi -> GetXaxis() -> SetTitle("#varphi");

  TH2D *histGridPhiTilde = new TH2D("histGridPhiTilde","",100,0.,2*gPi,100.,0.,gHistFit[2] -> GetMaximum() + 0.3*gHistFit[2] -> GetMaximum());
  histGridPhiTilde -> GetXaxis() -> SetTitle("#tilde{#varphi}");

  gHistFit[0] -> SetMarkerStyle(20); gHistFit[0] -> SetMarkerColor(kBlack); gHistFit[0] -> SetLineColor(kBlack);
  gHistFit[1] -> SetMarkerStyle(20); gHistFit[1] -> SetMarkerColor(kBlack); gHistFit[1] -> SetLineColor(kBlack);
  gHistFit[2] -> SetMarkerStyle(20); gHistFit[2] -> SetMarkerColor(kBlack); gHistFit[2] -> SetLineColor(kBlack);

  TPaveText *paveTextPolCosTheta = new TPaveText(0.2,0.85,0.8,0.98,"NDC");
  paveTextPolCosTheta -> SetFillColor(kWhite);
  paveTextPolCosTheta -> AddText(gPlotTitle.c_str());
  paveTextPolCosTheta -> AddText(Form("#lambda_{#theta} = %3.2f #pm %3.2f",fLambdaTheta,fErrorLambdaTheta));
  paveTextPolCosTheta -> AddText(Form("#chi^{2}/NDF = %3.2f",gFuncFit[0] -> GetChisquare()/gFuncFit[0] -> GetNDF()));

  TPaveText *paveTextPolPhi = new TPaveText(0.2,0.85,0.8,0.98,"NDC");
  paveTextPolPhi -> SetFillColor(kWhite);
  paveTextPolPhi -> AddText(gPlotTitle.c_str());
  paveTextPolPhi -> AddText(Form("#lambda_{#varphi} = %3.2f #pm %3.2f",fLambdaPhi,fErrorLambdaPhi));
  paveTextPolPhi -> AddText(Form("#chi^{2}/NDF = %3.2f",gFuncFit[1] -> GetChisquare()/gFuncFit[1] -> GetNDF()));

  TPaveText *paveTextPolPhiTilde = new TPaveText(0.2,0.85,0.8,0.98,"NDC");
  paveTextPolPhiTilde -> SetFillColor(kWhite);
  paveTextPolPhiTilde -> AddText(gPlotTitle.c_str());
  paveTextPolPhiTilde -> AddText(Form("#lambda_{#theta#varphi} = %3.2f #pm %3.2f",fLambdaThetaPhi,fErrorLambdaThetaPhi));
  paveTextPolPhiTilde -> AddText(Form("#chi^{2}/NDF = %3.2f",gFuncFit[2] -> GetChisquare()/gFuncFit[2] -> GetNDF()));

  TCanvas *canvasHistFitDec_CosTheta = new TCanvas("canvasHistFitDec_CosTheta","canvasHistFitDec_CosTheta",600,600);
  histGridCosTheta -> Draw(); gHistFit[0] -> Draw("EPsame"); gFuncFit[0] -> Draw("same");
  paveTextPolCosTheta -> Draw();
  //latexTitle -> DrawLatex(0.05,27000.,gPlotTitle.c_str());
  //latexTitle -> DrawLatex(0.3,25000.,Form("#lambda_{#theta} = %3.2f #pm %3.2f",fLambdaTheta,fErrorLambdaTheta));
  //latexTitle -> DrawLatex(0.2,3000.,Form("#chi^{2}/NDF = %3.2f",chiSquare/NDF));

  TCanvas *canvasHistFitDec_Phi = new TCanvas("canvasHistFitDec_Phi","canvasHistFitDec_Phi",600,600);
  histGridPhi -> Draw(); gHistFit[1] -> Draw("EPsame"); gFuncFit[1] -> Draw("same");
  paveTextPolPhi -> Draw();
  //latexTitle -> DrawLatex(0.2,9000.,gPlotTitle.c_str());
  //latexTitle -> DrawLatex(1.,8000.,Form("#lambda_{#varphi} = %3.2f #pm %3.2f",fLambdaPhi,fErrorLambdaPhi));
  //latexTitle -> DrawLatex(1,1000.,Form("#chi^{2}/NDF = %3.2f",chiSquare/NDF));

  TCanvas *canvasHistFitDec_PhiTilde = new TCanvas("canvasHistFitDec_PhiTilde","canvasHistFitDec_PhiTilde",600,600);
  histGridPhiTilde -> Draw(); gHistFit[2] -> Draw("EPsame"); gFuncFit[2] -> Draw("same");
  paveTextPolPhiTilde -> Draw();
  //latexTitle -> DrawLatex(0.1,4800.,gPlotTitle.c_str());
  //latexTitle -> DrawLatex(2.,4000.,Form("#lambda_{#theta#varphi} = %3.2f #pm %3.2f",fLambdaThetaPhi,fErrorLambdaThetaPhi));
  //latexTitle -> DrawLatex(2,500.,Form("#chi^{2}/NDF = %3.2f",chiSquare/NDF));

  if(saveCanvas){
    canvasHistFitDec_CosTheta -> SaveAs(Form("%s_CosTheta.pdf",nameCanvas.c_str()));
    canvasHistFitDec_Phi -> SaveAs(Form("%s_Phi.pdf",nameCanvas.c_str()));
    canvasHistFitDec_PhiTilde -> SaveAs(Form("%s_PhiTilde.pdf",nameCanvas.c_str()));
  }

  delete canvasHistFitDec_CosTheta;
  delete canvasHistFitDec_Phi;
  delete canvasHistFitDec_PhiTilde;
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

  // for Jpsi
  //Int_t fitBinMin[3] = {2,0,0};
  //Int_t fitBinMax[3] = {4,2,2};
  // for Upsilon
  Int_t fitBinMin[3];
  Int_t fitBinMax[3];

  fitBinMin[0] = gMinFitRangeCosTheta;
  fitBinMin[1] = gMinFitRangePhi;
  fitBinMin[2] = gMinFitRangePhiTilde;

  fitBinMax[0] = gMaxFitRangeCosTheta;
  fitBinMax[1] = gMaxFitRangePhi;
  fitBinMax[2] = gMaxFitRangePhiTilde;

  //cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
  //cout << gMinFitRangeCosTheta << " " << gMinFitRangePhi << " " << gMinFitRangePhiTilde << endl;
  //cout << gMaxFitRangeCosTheta << " " << gMaxFitRangePhi << " " << gMaxFitRangePhiTilde << endl;
  //cout << gHistFit[0] -> GetBinCenter(fitBinMin[0] + 1) << " " << gHistFit[0] -> GetBinCenter(gHistFit[0] -> GetSize() - fitBinMax[0] - 1) << endl;
  //cout << gHistFit[1] -> GetBinCenter(fitBinMin[1] + 1) << " " << gHistFit[1] -> GetBinCenter(gHistFit[1] -> GetSize() - fitBinMax[1] - 1) << endl;
  //cout << gHistFit[2] -> GetBinCenter(fitBinMin[2] + 1) << " " << gHistFit[2] -> GetBinCenter(gHistFit[2] -> GetSize() - fitBinMax[2] - 1) << endl;

  //cout << "Histo range" << endl;
  //for(int i = 0;i < 7;i++){cout << gHistFit[0] -> GetBinContent(i+1) << " -- " << gHistFit[0] -> GetBinCenter(i+1) << endl;}
  //cout << "Fit range" << endl;
  //for(int j = fitBinMin[0];j < gHistFit[0] -> GetSize() - fitBinMax[0];j++){cout << gHistFit[0] -> GetBinContent(j+1) << " -- " << gHistFit[0] -> GetBinCenter(j+1) << endl;}
  //cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;

  for(int i = 0;i < gNDistrib;i++){
    for(int j = fitBinMin[i];j < gHistFit[i] -> GetSize() - fitBinMax[i];j++){
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
  globalChiSquare = chiSquare;
}
//______________________________________________________________________________
void polarizationBarbatruccoFCN(Int_t &npar, Double_t *gin, Double_t &gChiSquare, Double_t *par, Int_t iflag){
  Double_t angleVar;
  Double_t val, errVal;

  /*
  // HELICITY
  Double_t parCosThetaHE[2];
  parCosThetaHE[0] = par[0];
  parCosThetaHE[1] = par[3];

  Double_t parPhiHE[3];
  parPhiHE[0] = par[1];
  parPhiHE[1] = par[3];
  parPhiHE[2] = (par[7] - par[3])/(par[7] + 3);

  Double_t parPhiTildeHE[3];
  parPhiTildeHE[0] = par[2];
  parPhiTildeHE[1] = par[3];
  parPhiTildeHE[2] = par[4];

  // COLLINS-SOPER
  Double_t parCosThetaCS[2];
  parCosThetaCS[0] = par[0];
  parCosThetaCS[1] = par[5];

  Double_t parPhiCS[3];
  parPhiCS[0] = par[1];
  parPhiCS[1] = par[5];
  parPhiCS[2] = (par[7] - par[5])/(par[7] + 3);

  Double_t parPhiTildeCS[3];
  parPhiTildeCS[0] = par[2];
  parPhiTildeCS[1] = par[5];
  parPhiTildeCS[2] = par[6];
  */


  // HELICITY
  Double_t parCosThetaHE[2];
  parCosThetaHE[0] = par[0];
  parCosThetaHE[1] = par[3];

  Double_t parPhiHE[3];
  parPhiHE[0] = par[1];
  parPhiHE[1] = par[3];
  parPhiHE[2] = (par[10] - par[3])/(par[10] + 3);

  Double_t parPhiTildeHE[3];
  parPhiTildeHE[0] = par[2];
  parPhiTildeHE[1] = par[3];
  parPhiTildeHE[2] = par[4];

  // COLLINS-SOPER
  Double_t parCosThetaCS[2];
  parCosThetaCS[0] = par[5];
  parCosThetaCS[1] = par[8];

  Double_t parPhiCS[3];
  parPhiCS[0] = par[6];
  parPhiCS[1] = par[8];
  parPhiCS[2] = (par[10] - par[8])/(par[10] + 3);

  Double_t parPhiTildeCS[3];
  parPhiTildeCS[0] = par[7];
  parPhiTildeCS[1] = par[8];
  parPhiTildeCS[2] = par[9];


  Double_t func = 0, pull = 0, chiSquare = 0;

  // J/psi
  //Int_t fitBinMin[6] = {2,0,0,2,0,0};
  //Int_t fitBinMax[6] = {4,2,2,4,2,2};
  // Upsilon(1S)
  Int_t fitBinMin[6] = {0,0,0,0,0,0};
  Int_t fitBinMax[6] = {3,2,2,3,2,2};

  for(int i = 0;i < gNDistrib;i++){
    //for(int j = fitBinMin[i];j < gHistBarbFit[i] -> GetSize() - fitBinMax[i];j++){
      for(int j = 0;j < gHistBarbFit[i] -> GetSize();j++){
      val = gHistBarbFit[i] -> GetBinContent(j+1);
      errVal = gHistBarbFit[i] -> GetBinError(j+1);

      if(val == 0. && errVal == 0.){continue;}
      if(i == 0){
        angleVar = gCosThetaValues[j];
        pull = (val - funcCosTheta(angleVar,parCosThetaHE))/errVal;
      }
      if(i == 1){
        angleVar = gPhiValues[j];
        pull = (val - funcPhi(angleVar,parPhiHE))/errVal;
      }
      if(i == 2){
        angleVar = gPhiTildeValues[j];
        pull = (val - funcPhiTilde(angleVar,parPhiTildeHE))/errVal;
      }
      if(i == 3){
        angleVar = gCosThetaValues[j];
        pull = (val - funcCosTheta(angleVar,parCosThetaCS))/errVal;
      }
      if(i == 4){
        angleVar = gPhiValues[j];
        pull = (val - funcPhi(angleVar,parPhiCS))/errVal;
      }
      if(i == 5){
        angleVar = gPhiTildeValues[j];
        pull = (val - funcPhiTilde(angleVar,parPhiTildeCS))/errVal;
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
