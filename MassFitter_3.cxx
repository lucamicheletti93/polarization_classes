// STL includes
#include <Riostream.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>

// ROOT includes
#include <TROOT.h>
#include <TObjArray.h>
#include <TMinuit.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TPaveText.h>
#include <TCollection.h>
#include <TKey.h>
#include <TGaxis.h>

// My includes
#include "/home/luca/GITHUB/polarization_classes/FunctionsLibrary_3.C"
#include "MassFitter_3.h"

// List of fitting-functions
/*
Double_t Func_VWG(Double_t *, Double_t *);
Double_t Func_POL4EXP(Double_t *, Double_t *);
Double_t Func_Jpsi_CB2_VWG(Double_t *, Double_t *);
Double_t Func_Jpsi_CB2_POL4EXP(Double_t *, Double_t *);
Double_t Func_Jpsi_CB2_fix(Double_t *, Double_t *);
Double_t Func_Jpsi_NA60_VWG(Double_t *, Double_t *);
Double_t Func_Jpsi_NA60_POL4EXP(Double_t *, Double_t *);
Double_t Func_Jpsi_NA60_fix(Double_t *, Double_t *);
Double_t Func_Psi2S_CB2_VWG(Double_t *, Double_t *);
Double_t Func_Psi2S_CB2_POL4EXP(Double_t *, Double_t *);
Double_t Func_Psi2S_CB2_fix(Double_t *, Double_t *);
Double_t Func_Psi2S_NA60_VWG(Double_t *, Double_t *);
Double_t Func_Psi2S_NA60_POL4EXP(Double_t *, Double_t *);
Double_t Func_Psi2S_NA60_fix(Double_t *, Double_t *);
Double_t Func_tot_CB2_VWG(Double_t *, Double_t *);
Double_t Func_tot_CB2_POL4EXP(Double_t *, Double_t *);
Double_t Func_tot_NA60_VWG(Double_t *, Double_t *);
Double_t Func_tot_NA60_POL4EXP(Double_t *, Double_t *);
*/

// List of global variables
//double gPi = TMath::Pi();
//Long_t *dummy1 = 0, *dummy2 = 0, *dummy3 = 0, *dummy4 = 0;

ClassImp(MassFitter_3)

//______________________________________________________________________________
MassFitter_3::MassFitter_3(): TObject() {
  // default constructor
}
//______________________________________________________________________________
MassFitter_3::MassFitter_3(TH1D *histMinv): TObject() {
  fHistMinv = histMinv;
  fScalingFactorJpsiSigma = 1.;
  fSpecialFitConditions = kFALSE;
  fJpsiWidthFixed = kFALSE;
  fJpsiMassFixedToPDG = kFALSE;
  fJpsiMassFixedToIntegrated = kFALSE;
  fTailParametersFixed = kTRUE;
  fIndexCosTheta = 100;
  fIndexPhi = 100;
  fRebin = 1;
  double parTailsCB2[4] = {0.970014,3.9789,2.29866,3.0301};                                                 // tails parameter embedding Pb-Pb 5.02 TeV (0 < pT < 12 GeV/c)
  double parTailsNA60[8] = {0.227822,1.14006,0.0357395,0.187828,1.22477,0.0569524,-0.630315,2.36889};       // tails parameter embedding Pb-Pb 5.02 TeV (0 < pT < 12 GeV/c)

  for(int i = 0;i < 4;i++){fParTailsCB2[i] = parTailsCB2[i];}
  for(int i = 0;i < 8;i++){fParTailsNA60[i] = parTailsNA60[i];}
  // standard constructor

  fPhilippeCorrectionFactor = 0.;
}
//______________________________________________________________________________
MassFitter_3::~MassFitter_3() {
  // destructor
}
//______________________________________________________________________________
void MassFitter_3::SetScalingFactor(Double_t scalingFactor){
    fScalingFactorJpsiSigma = scalingFactor;
}
//______________________________________________________________________________
void MassFitter_3::SetFitRange(Double_t minFitRange, Double_t maxFitRange){
    fMinFitRange = minFitRange;
    fMaxFitRange = maxFitRange;
}
//______________________________________________________________________________
void MassFitter_3::SetSpecialFitConditions(Int_t rebin, Bool_t jpsiMassFixedToPDG){
    fSpecialFitConditions = kTRUE;
    fRebin = rebin;
    fJpsiMassFixedToPDG = jpsiMassFixedToPDG;
}
//______________________________________________________________________________
void MassFitter_3::SetJpsiMass(Double_t massJpsi){
    fJpsiMassFixedToIntegrated = kTRUE;
    fMassIntegrated = massJpsi;
}
//______________________________________________________________________________
void MassFitter_3::SetJpsiWidth(Double_t sigmaJpsi){
    fJpsiWidthFixed = kTRUE;
    fSigmaJpsi = sigmaJpsi;
}
//______________________________________________________________________________
void MassFitter_3::SetPhilippeCorrectionFactor(Double_t philippeCorrectionFactor){
  fPhilippeCorrectionFactor = philippeCorrectionFactor;
}
//______________________________________________________________________________
void MassFitter_3::SetCosThetaPhiIndex(Int_t indexCosTheta, Int_t indexPhi){
    fIndexCosTheta = indexCosTheta;
    fIndexPhi = indexPhi;
}
//______________________________________________________________________________
void MassFitter_3::SetTailsParameters(Double_t tmpParTailsCB2[4], Double_t tmpParTailsNA60[8]){
    for(int i = 0;i < 4;i++){fParTailsCB2[i] = tmpParTailsCB2[i];}
    for(int i = 0;i < 8;i++){fParTailsNA60[i] = tmpParTailsNA60[i];}
}
//______________________________________________________________________________
void MassFitter_3::SetBinning(vector <Double_t> CosThetaValues, vector <Double_t> PhiValues, vector <Double_t> PhiTildeValues) {
  for(int i = 0;i < (int) CosThetaValues.size();i++){fCosThetaValues.push_back(CosThetaValues[i]);}
  fNCosThetaBins = fCosThetaValues.size() - 1;
  for(int i = 0;i < (int) PhiValues.size();i++){fPhiValues.push_back(PhiValues[i]);}
  fNPhiBins = fPhiValues.size() - 1;
  for(int i = 0;i < (int) PhiTildeValues.size();i++){fPhiTildeValues.push_back(PhiTildeValues[i]);}
  fNPhiTildeBins = fPhiTildeValues.size() - 1;
}
//______________________________________________________________________________
void MassFitter_3::fit_of_minv(string sigShape, string bkgShape, string outputDir, string plotType, Bool_t showPlot, Bool_t savePlot){
  gStyle -> SetOptStat(0);
  TGaxis::SetMaxDigits(1);

  fSigShape = sigShape;
  fBkgShape = bkgShape;
  fOutputDir = outputDir;

  fPlotType = plotType;
  fSavePlot = savePlot;

  double min_fit_range = fMinFitRange;
  double max_fit_range = fMaxFitRange;

  //============================================================================
  // FIT STARTING PARAMETERS
  //============================================================================
  int nParBkg, nParSig, nParTails;

  // background function
  double parBkgVWG[4] = {500000.,0.6,0.2,0.2};
  double parBkgPOL4EXP[6] = {5000.,-10.,1e6,-1e6,1e3,1e3};
  double *parBkg;
  if(bkgShape == "VWG"){parBkg = &(parBkgVWG[0]); nParBkg = sizeof(parBkgVWG)/sizeof(double);}
  if(bkgShape == "POL4EXP"){parBkg = &(parBkgPOL4EXP[0]); nParBkg = sizeof(parBkgPOL4EXP)/sizeof(double);}

  double parSigCB2[3] = {500.,3.096,7.0e-02};
  double parSigNA60[3] = {500.,3.096,7.0e-02};
  double *parSig;
  if(sigShape == "CB2"){parSig = &(parSigCB2[0]); nParSig = sizeof(parSigCB2)/sizeof(double);}
  if(sigShape == "NA60"){parSig = &(parSigNA60[0]); nParSig = sizeof(parSigCB2)/sizeof(double);}

  //double parTailsCB2[4] = {0.970014,3.9789,2.29866,3.0301};                                                 // tails parameter embedding Pb-Pb 5.02 TeV (0 < pT < 12 GeV/c)
  //double parTailsNA60[8] = {0.227822,1.14006,0.0357395,0.187828,1.22477,0.0569524,-0.630315,2.36889};       // tails parameter embedding Pb-Pb 5.02 TeV (0 < pT < 12 GeV/c)
  double *parTails;
  if(sigShape == "CB2"){parTails = &(fParTailsCB2[0]); nParTails = sizeof(fParTailsCB2)/sizeof(double);}
  if(sigShape == "NA60"){parTails = &(fParTailsNA60[0]); nParTails = sizeof(fParTailsNA60)/sizeof(double);}

  cout << nParBkg << " " << nParSig << " " << nParTails << endl;

  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  // N.B. : NA60 tails parameters
  // As reported in the AN
  // double parTailsNA60[8] = {-0.630315,0.227822,1.14006,0.0357395,2.36889,0.187828,1.22477,0.0569524};
  // double parTailsNA60[8] = {aL,p1L,p2L,p3L,aR,p1R,p2R,p3R};
  // As requested in my FunctionsLibrary.C
  // double parTailsNA60[8] = {0.227822,1.14006,0.0357395,0.187828,1.22477,0.0569524,-0.630315,2.36889};
  // double parTailsNA60[8] = {p1L,p2L,p3L,p1R,p2R,p3R,aL,aR};
  //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  //============================================================================
  // INIZIALIZING THE HISTOGRAM
  //============================================================================
  TH1D *histMinvBkg = (TH1D*) fHistMinv -> Clone("histMinvBkg");
  double m_width = histMinvBkg -> GetBinWidth(1);
  double m_min = histMinvBkg -> FindBin(2.9);
  double m_max = histMinvBkg -> FindBin(3.3);

  cout << "-----------------------------------" << endl;
  cout << "Bin width : " << m_width << endl;
  cout << "bin min : " << m_min << "; bin max : " << m_max << endl;
  cout << "-----------------------------------" << endl;

  for(int i = m_min; i < m_max;i++){
    histMinvBkg -> SetBinContent(i+1,0);
    histMinvBkg -> SetBinError(i+1,0);
  }

  if(fSpecialFitConditions){
    fHistMinv -> Rebin(fRebin);
    //min_fit_range = 2.2;
    //max_fit_range = 4.5;
  }

  //============================================================================
  // FIT OF THE BACKGROUND
  //============================================================================
  TF1 *funcBkg;
  if(bkgShape == "VWG"){funcBkg = new TF1("funcBkg",Func_VWG,min_fit_range,max_fit_range,nParBkg);}
  if(bkgShape == "POL4EXP"){funcBkg = new TF1("funcBkg",Func_POL4EXP,min_fit_range,max_fit_range,nParBkg);}
  funcBkg -> SetLineColor(kRed);
  funcBkg -> SetParameters(parBkg);
  histMinvBkg -> Fit(funcBkg,"R0Q");

  //============================================================================
  // FIT OF THE JPSI SIGNAL
  //============================================================================
  TF1 *funcSigJpsi;
  if(sigShape == "CB2" && bkgShape == "VWG"){funcSigJpsi = new TF1("funcSigJpsi",Func_Jpsi_CB2_VWG,2.9,3.2,nParBkg + nParSig + nParTails);}
  if(sigShape == "NA60" && bkgShape == "VWG"){funcSigJpsi = new TF1("funcSigJpsi",Func_Jpsi_NA60_VWG,2.9,3.2,nParBkg + nParSig + nParTails);}
  if(sigShape == "CB2" && bkgShape == "POL4EXP"){funcSigJpsi = new TF1("funcSigJpsi",Func_Jpsi_CB2_POL4EXP,2.9,3.2,nParBkg + nParSig + nParTails);}
  if(sigShape == "NA60" && bkgShape == "POL4EXP"){funcSigJpsi = new TF1("funcSigJpsi",Func_Jpsi_NA60_POL4EXP,2.9,3.2,nParBkg + nParSig + nParTails);}
  funcSigJpsi -> SetNpx(100000);
  funcSigJpsi -> SetParameter(nParBkg + 0,parSig[0]);
  funcSigJpsi -> SetParameter(nParBkg + 1,parSig[1]);
  funcSigJpsi -> FixParameter(nParBkg + 2,parSig[2]);
  for(int i = 0;i < nParTails;i++){funcSigJpsi -> FixParameter(nParBkg + nParSig + i,parTails[i]);}
  fHistMinv -> Fit(funcSigJpsi,"R0Q");

  //============================================================================
  // FIT OF THE TOTAL SPECTRUM
  //============================================================================
  //TF1 *funcTot;
  if(sigShape == "CB2" && bkgShape == "VWG"){fFuncTot = new TF1("fFuncTot",Func_tot_CB2_VWG,min_fit_range,max_fit_range,nParBkg + nParSig + nParTails + 1);}
  if(sigShape == "NA60" && bkgShape == "VWG"){fFuncTot = new TF1("fFuncTot",Func_tot_NA60_VWG,min_fit_range,max_fit_range,nParBkg + nParSig + nParTails + 1);}
  if(sigShape == "CB2" && bkgShape == "POL4EXP"){fFuncTot = new TF1("fFuncTot",Func_tot_CB2_POL4EXP,min_fit_range,max_fit_range,nParBkg + nParSig + nParTails + 1);}
  if(sigShape == "NA60" && bkgShape == "POL4EXP"){fFuncTot = new TF1("fFuncTot",Func_tot_NA60_POL4EXP,min_fit_range,max_fit_range,nParBkg + nParSig + nParTails + 1);}

  TFitResultPtr fit_ptr;
  for(int i = 0;i < 100;i++){
    //cout << i << " " << fFuncTot -> GetParameter(nParBkg + 1) << " " << fFuncTot -> GetParameter(nParBkg + 2) << endl;
    if(i == 0){
      for(int i = 0;i < nParBkg;i++){fFuncTot -> SetParameter(i,funcBkg -> GetParameter(i));}

      fFuncTot -> SetParameter(nParBkg + 0,funcSigJpsi -> GetParameter(nParBkg));
      fFuncTot -> SetParLimits(nParBkg + 0,0,10000000);
      if(fJpsiMassFixedToPDG){fFuncTot -> FixParameter(nParBkg + 1,3.096);}
      if(fJpsiMassFixedToIntegrated){fFuncTot -> FixParameter(nParBkg + 1,fMassIntegrated);}
      if(fJpsiMassFixedToIntegrated == kFALSE && fJpsiMassFixedToPDG == kFALSE){
        fFuncTot -> SetParameter(nParBkg + 1,3.096);
        fFuncTot -> SetParLimits(nParBkg + 1,2.9,3.2);
      }
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // J/psi width conditions
      if(fJpsiWidthFixed){
        // RESTORE IN CASE OF PROBLEMS
        //fSigmaJpsi = fScalingFactorJpsiSigma*fSigmaJpsi;
        // PHILIPPE CORRECTION FACTOR
        fSigmaJpsi = fScalingFactorJpsiSigma*fSigmaJpsi + fPhilippeCorrectionFactor;
        fFuncTot -> FixParameter(nParBkg + 2,fSigmaJpsi);
      }
      else{
        fFuncTot -> SetParameter(nParBkg + 2,7.0e-02);
        fFuncTot -> SetParLimits(nParBkg + 2,2.0e-02,2.0e-01);
      }
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(fTailParametersFixed){for(int i = 0;i < nParTails;i++){fFuncTot -> FixParameter(nParBkg + nParSig + i,funcSigJpsi -> GetParameter(nParBkg + nParSig + i));}}
      else{for(int i = 0;i < nParTails;i++){fFuncTot -> SetParameter(nParBkg + nParSig + i,funcSigJpsi -> GetParameter(nParBkg + nParSig + i));}}

      // Psi(2S) : only normalization is a parameter
      fFuncTot -> SetParameter(nParBkg + nParSig + nParTails,2.37767e-02);
      fFuncTot -> SetParLimits(nParBkg + nParSig + nParTails,0.,0.1);

    }
    else{fFuncTot -> SetParameters(fFuncTot -> GetParameters());}
    fit_ptr = (TFitResultPtr) fHistMinv -> Fit(fFuncTot,"RLS0Q");
    if(gMinuit->fCstatu.Contains("CONVERGED")) break;
  }

  fChiSquare_NDF = fFuncTot -> GetChisquare()/fFuncTot -> GetNDF();

  if(gMinuit -> fCstatu.Contains("FAILED")){
    cout << "WARNING : FIT STATUS FAILED" << endl;
    fFitStatus = "FAILED";
    delete histMinvBkg, fHistMinv;
    delete funcBkg, funcSigJpsi, fFuncTot;
    fNJpsi = 0.; fStatJpsi = 0.; fMassJpsi = 0.; fErrMassJpsi = 0.; fSigmaJpsi = 0.; fErrSigmaJpsi = 0.; fChiSquare_NDF = 666.;
    string nameHistMinv = (string) fHistMinv -> GetName() + ".png";

    TCanvas *canvasMinv = new TCanvas("canvasMinv","canvasMinv",65,73,900,806);
    fHistMinv -> Draw();

    TLatex *latexTitle = new TLatex(); latexTitle -> SetTextSize(0.04); latexTitle -> SetNDC(); latexTitle -> SetTextFont(42);
    latexTitle -> DrawLatex(0.4,0.5,"FIT FAILED");

    if(strstr(nameHistMinv.c_str(),"HE")){
        if(strstr(outputDir.c_str(),"BadFit")){canvasMinv -> SaveAs(Form("%s/%s",outputDir.c_str(),nameHistMinv.c_str()));}
        else{canvasMinv -> SaveAs(Form("%s/Helicity/%s_%s/%s",outputDir.c_str(),sigShape.c_str(),bkgShape.c_str(),nameHistMinv.c_str()));}
        canvasMinv -> Close();
    }
    if(strstr(nameHistMinv.c_str(),"CS")){
        if(strstr(outputDir.c_str(),"BadFit")){canvasMinv -> SaveAs(Form("%s/%s",outputDir.c_str(),nameHistMinv.c_str()));}
        else{canvasMinv -> SaveAs(Form("%s/Collins_Soper/%s_%s/%s",outputDir.c_str(),sigShape.c_str(),bkgShape.c_str(),nameHistMinv.c_str()));}
        canvasMinv -> Close();
    }
    return;
  }
  else{fFitStatus = "SUCCESS";}

  const int dimParBkg = nParBkg;
  const int dimParSig = nParSig + nParTails;

  // Getting sub matrix for J/psi & Psi(2S)(error calculation)
  /*TMatrixDSym cov = fit_ptr -> GetCovarianceMatrix().GetSub(dimParBkg,(dimParBkg+dimParSig)-1,dimParBkg,(dimParBkg+dimParSig)-1);
  double Jpsi_par[dimParSig];
  for(int i = 0;i < dimParSig;i++){Jpsi_par[i] = fFuncTot -> GetParameter(dimParBkg+i);}*/
  TMatrixDSym covJpsi  = fit_ptr -> GetCovarianceMatrix().GetSub(dimParBkg,(dimParBkg+dimParSig)-1,dimParBkg,(dimParBkg+dimParSig)-1);
  TMatrixDSym covPsi2S = fit_ptr -> GetCovarianceMatrix().GetSub(dimParBkg,(dimParBkg+dimParSig),dimParBkg,(dimParBkg+dimParSig));
  double Jpsi_par[dimParSig];
  for(int i = 0;i < dimParSig;i++){Jpsi_par[i] = fFuncTot -> GetParameter(dimParBkg+i);}
  double Psi2S_par[dimParSig];
  for(int i = 0;i < dimParSig+1;i++){Psi2S_par[i] = fFuncTot -> GetParameter(dimParBkg+i);}

  //============================================================================
  // PLOT OF JPSI AND PSI(2S) SHAPES
  //============================================================================
  //TF1 *fFuncBkgFix;
  if(bkgShape == "VWG"){fFuncBkgFix = new TF1("fFuncBkgFix",Func_VWG,2.,5.,nParBkg);}
  if(bkgShape == "POL4EXP"){fFuncBkgFix = new TF1("fFuncBkgFix",Func_POL4EXP,2.,5.,nParBkg);}
  fFuncBkgFix -> SetNpx(1000);
  for(int i = 0;i < nParBkg;i++){fFuncBkgFix -> SetParameter(i,fFuncTot -> GetParameter(i));}
  fFuncBkgFix -> SetLineStyle(4);
  fFuncBkgFix -> SetLineColor(kGray+1);
  fFuncBkgFix -> Draw("same");

  //TF1 *funcSigJpsiFix;
  if(sigShape == "CB2"){
    fFuncSigJpsiFix = new TF1("fFuncSigJpsiFix",Func_Jpsi_CB2_fix,2.,5.,nParSig + nParTails);
    fFuncSigPsi2SFix = new TF1("fFuncSigPsi2SFix",Func_Psi2S_CB2_fix,2.,5.,nParSig + nParTails + 1);
  }
  if(sigShape == "NA60"){
    fFuncSigJpsiFix = new TF1("fFuncSigJpsiFix",Func_Jpsi_NA60_fix,2.,5.,nParSig + nParTails);
    fFuncSigPsi2SFix = new TF1("fFuncSigPsi2SFix",Func_Psi2S_NA60_fix,2.,5.,nParSig + nParTails + 1);
  }
  fFuncSigJpsiFix -> SetNpx(1000);
  for(int i = 0;i < nParSig + nParTails;i++){fFuncSigJpsiFix -> SetParameter(i,fFuncTot -> GetParameter(nParBkg + i));}
  fFuncSigJpsiFix -> SetLineColor(kBlue+1);
  fFuncSigJpsiFix -> SetFillStyle(3335);
  fFuncSigJpsiFix -> SetFillColorAlpha(kBlue,0.3);
  fFuncSigJpsiFix -> Draw("same");

  fFuncSigPsi2SFix -> SetParameter(0,fFuncTot -> GetParameter(nParBkg + 0));
  fFuncSigPsi2SFix -> SetParameter(1,(fFuncTot -> GetParameter(nParBkg + 1)) + (3.686097-3.0969));
  fFuncSigPsi2SFix -> SetParameter(2,(fFuncTot -> GetParameter(nParBkg + 2))*1.07365);
  for(int i = nParSig;i < nParSig + nParTails + 1;i++){fFuncSigPsi2SFix -> SetParameter(i,fFuncTot -> GetParameter(nParBkg + i));}
  fFuncSigPsi2SFix -> SetLineColor(kGreen+1);
  fFuncSigPsi2SFix -> SetFillStyle(3353);
  fFuncSigPsi2SFix -> SetFillColorAlpha(kGreen,0.4);

  fMassJpsi = fFuncTot -> GetParameter(nParBkg + 1);
  fErrMassJpsi = fFuncTot -> GetParError(nParBkg + 1);
  fSigmaJpsi = fFuncTot -> GetParameter(nParBkg + 2);
  fErrSigmaJpsi = fFuncTot -> GetParError(nParBkg + 2);
  double sigma_min_Jpsi = fFuncTot -> GetParameter(nParBkg + 1) - 3*(fFuncTot -> GetParameter(nParBkg + 2));
  double sigma_max_Jpsi = fFuncTot -> GetParameter(nParBkg + 1) + 3*(fFuncTot -> GetParameter(nParBkg + 2));
  double N_Jpsi_3sigma = fFuncSigJpsiFix -> Integral(sigma_min_Jpsi,sigma_max_Jpsi)/m_width;
  double N_bck_Jpsi_3sigma = fFuncBkgFix -> Integral(sigma_min_Jpsi,sigma_max_Jpsi)/m_width;
  double SB_Jpsi = N_Jpsi_3sigma/N_bck_Jpsi_3sigma;

  fNJpsi = fFuncSigJpsiFix -> Integral(0,5)/m_width;
  fStatJpsi = fFuncSigJpsiFix -> IntegralError(0.,5.,Jpsi_par,covJpsi.GetMatrixArray())/m_width;

  fMassPsi2S = (fFuncTot -> GetParameter(nParBkg + 1))+(3.686097-3.0969);
  fErrMassPsi2S = fFuncTot -> GetParError(nParBkg + 1);
  fSigmaPsi2S = (fFuncTot -> GetParameter(nParBkg + 2))*1.07365;
  fErrSigmaPsi2S = (fFuncTot -> GetParError(nParBkg + 2))*1.07365;
  fNPsi2S = fFuncSigPsi2SFix -> Integral(0,5)/m_width;
  fStatPsi2S = fFuncSigPsi2SFix -> IntegralError(0,5,Psi2S_par,covPsi2S.GetMatrixArray())/m_width;

  if(showPlot){
    if(fPlotType == "STANDARD"){PlotStandard(fSavePlot);}
    if(fPlotType == "PUBLICATION"){PlotPublications(fSavePlot);}
  }

  /*
  delete histGridMinv, histGridRatio, histFuncTot;
  delete histMinvBkg, fHistMinv;
  delete funcBkg, funcSigJpsi, funcTot, funcBkgFix, funcSigJpsiFix;
  delete letexTitle;
  //delete axisHistMinv, lineRatio, axisRatio, padHistMinv, padRatio;
  delete axisHistMinv, lineRatio, axisRatio;
  delete canvasMinv;
  */
}
//______________________________________________________________________________
void MassFitter_3::PlotStandard(bool savePlot){
  TGaxis::SetMaxDigits(1);
  gStyle -> SetCanvasPreferGL(kTRUE);
  //============================================================================
  // PREPARING THE FRAME
  //============================================================================
  double max_histo_value = fHistMinv -> GetMaximum();
  max_histo_value = max_histo_value + 0.6*max_histo_value;

  TH2D *histGridMinv = new TH2D("histGridMinv","",120,2,5,100,0,max_histo_value);
  histGridMinv -> GetYaxis() -> SetTitle(Form("Counts per %3.0f MeV/#it{c}^{2}",fHistMinv -> GetBinWidth(1)*1000));
  histGridMinv -> GetYaxis() -> CenterTitle(true);
  histGridMinv -> GetYaxis() -> SetTitleSize(0.04);
  histGridMinv -> GetYaxis() -> SetLabelSize(0);

  TGaxis *axisHistMinv = new TGaxis(2.,1.,2.,max_histo_value,1.,max_histo_value,510,"");
  axisHistMinv -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axisHistMinv -> SetLabelSize(25);

  TLine *lineRatio = new TLine(2.,1.,5.,1.0); lineRatio -> SetLineStyle(kDashed); lineRatio -> SetLineWidth(2);

  TH2D *histGridRatio = new TH2D("histGridRatio","",120,2,5,10,0.9,1.1);
  histGridRatio -> GetXaxis() -> SetTitle("#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})");
  histGridRatio -> GetXaxis() -> SetTitleSize(0.12);
  histGridRatio -> GetXaxis() -> SetTitleOffset(0.9);
  histGridRatio -> GetXaxis() -> SetLabelSize(0.12);
  histGridRatio -> GetYaxis() -> SetTitle("Data/Fit");
  histGridRatio -> GetYaxis() -> CenterTitle(true);
  histGridRatio -> GetYaxis() -> SetTitleSize(0.1);
  histGridRatio -> GetYaxis() -> SetTitleOffset(0.4);
  histGridRatio -> GetYaxis() -> SetLabelSize(0);

  TGaxis *axisRatio = new TGaxis(2.,0.92,2.,1.08,0.92,1.08,510,"");
  axisRatio -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axisRatio -> SetLabelSize(15);

  TLatex *letexTitle = new TLatex(); letexTitle -> SetTextSize(0.04); letexTitle -> SetNDC(); letexTitle -> SetTextFont(42);

  //============================================================================
  // DRAWING...
  //============================================================================
  TH1D *histFuncTot = new TH1D("histFuncTot","",120,2.,5.);
  if(fSpecialFitConditions){histFuncTot -> Rebin(fRebin);}
  for(int i = 0;i < 120;i++){
    histFuncTot -> SetBinContent(i+1,fFuncTot -> Eval(histFuncTot -> GetBinCenter(i+1)));
    histFuncTot -> SetBinError(i+1,0.);
  }

  TH1D *histRatio = (TH1D*) fHistMinv -> Clone("histRatio");
  histRatio -> SetMarkerStyle(20); histRatio -> SetMarkerSize(0.8); histRatio -> SetMarkerColor(kBlack);
  histRatio -> Divide(histFuncTot);

  fHistMinv -> SetMarkerStyle(20); fHistMinv -> SetMarkerSize(1.1); fHistMinv -> SetMarkerColor(kBlack);

  TLegend *legendMassFit1 = new TLegend(0.15,0.65,0.3,0.75); legendMassFit1 -> SetTextSize(0.04); legendMassFit1 -> SetBorderSize(0);
  legendMassFit1 -> AddEntry(fHistMinv,"Data","LP");
  legendMassFit1 -> AddEntry(fFuncSigPsi2SFix,"#psi(2S)","L");

  TLegend *legendMassFit2 = new TLegend(0.3,0.65,0.45,0.75); legendMassFit2 -> SetTextSize(0.04); legendMassFit2 -> SetBorderSize(0);
  legendMassFit2 -> AddEntry(fFuncTot,"Total fit","L");
  legendMassFit2 -> AddEntry(fFuncSigJpsiFix,"J/#psi","L");


  TCanvas *canvasMinv = new TCanvas("canvasMinv","canvasMinv",65,73,900,806);
  canvasMinv -> Range(1.825,-6776.052,5.019444,37862.12);
  canvasMinv -> SetFillColor(0);
  canvasMinv -> SetBorderMode(0);
  canvasMinv -> SetBorderSize(0);
  canvasMinv -> SetTickx(1);
  canvasMinv -> SetTicky(1);
  canvasMinv -> SetLeftMargin(0.18);
  canvasMinv -> SetBottomMargin(0.1518219);
  canvasMinv -> SetFrameBorderMode(0);
  canvasMinv -> SetFrameBorderMode(0);

  canvasMinv -> cd();
  TPad *padHistMinv = new TPad("padHistMinv","padHistMinv",0.,0.25,1.,1.); padHistMinv -> SetBottomMargin(0);
  padHistMinv -> Draw();
  padHistMinv -> cd();
  histGridMinv -> Draw();
  axisHistMinv -> Draw("same");
  fHistMinv -> Draw("sameE");
  fFuncBkgFix -> Draw("same");
  fFuncTot -> Draw("same");
  fFuncSigJpsiFix -> Draw("same");
  fFuncSigPsi2SFix -> Draw("same");
  legendMassFit1 -> Draw("same");
  legendMassFit2 -> Draw("same");

  letexTitle -> DrawLatex(0.15,0.82,Form("%2.1f < #it{m}_{#mu^{#plus}#mu^{#minus}} < %2.1f GeV/#it{c}^{2}",fMinFitRange,fMaxFitRange));
  letexTitle -> DrawLatex(0.15,0.76,Form("%s + %s",fSigShape.c_str(),fBkgShape.c_str()));
  letexTitle -> DrawLatex(0.55,0.82,Form("N_{J/#psi} = %2.0f #pm %2.0f",fNJpsi,fStatJpsi));
  letexTitle -> DrawLatex(0.55,0.76,Form("#it{m}_{J/#psi} = %4.3f #pm %4.3f GeV/#it{c}^{2}",fMassJpsi,fErrMassJpsi));
  letexTitle -> DrawLatex(0.55,0.7,Form("#it{#sigma}_{J/#psi} = %3.0f #pm %3.0f MeV/#it{c}^{2}",fSigmaJpsi*1000,fErrSigmaJpsi*1000));
  //letexTitle -> DrawLatex(0.55,0.64,Form("N_{#psi(2S)} = %2.0f #pm %2.0f",fNPsi2S,fStatPsi2S));
  //letexTitle -> DrawLatex(0.55,0.58,Form("#it{m}_{#psi(2S)} = %4.3f #pm %4.3f GeV/#it{c}^{2}",fMassPsi2S,fErrMassPsi2S));
  //letexTitle -> DrawLatex(0.55,0.52,Form("#it{#sigma}_{#psi(2S)} = %3.0f #pm %3.0f MeV/#it{c}^{2}",fSigmaPsi2S*1000,fErrSigmaPsi2S*1000));
  //if(fIndexCosTheta == 100) letexTitle -> DrawLatex(0.6,0.45,"-0.8 < cos#it{#theta} < 0.8");
  //else{letexTitle -> DrawLatex(0.6,0.45,Form("%3.2f < cos#it{#theta} < %3.2f",fCosThetaValues[fIndexCosTheta],fCosThetaValues[fIndexCosTheta+1]));}
  //if(fIndexPhi == 100) letexTitle -> DrawLatex(0.6,0.4,"0.5 < |#it{#varphi}| < 2.64 rad");
  //else{letexTitle -> DrawLatex(0.6,0.4,Form("%3.2f < |#it{#varphi}| < %3.2f rad",fPhiValues[fIndexPhi],fPhiValues[fIndexPhi+1]));}
  if(fIndexCosTheta == 100) letexTitle -> DrawLatex(0.6,0.64,"-0.8 < cos#it{#theta} < 0.8");
  else{letexTitle -> DrawLatex(0.6,0.64,Form("%3.2f < cos#it{#theta} < %3.2f",fCosThetaValues[fIndexCosTheta],fCosThetaValues[fIndexCosTheta+1]));}
  if(fIndexPhi == 100) letexTitle -> DrawLatex(0.6,0.58,"0.5 < |#it{#varphi}| < 2.64 rad");
  else{letexTitle -> DrawLatex(0.6,0.58,Form("%3.2f < |#it{#varphi}| < %3.2f rad",fPhiValues[fIndexPhi],fPhiValues[fIndexPhi+1]));}
  letexTitle -> DrawLatex(0.6,0.52,Form("#chi^{2}/ndf = %3.1f/%i",(double) fFuncTot -> GetChisquare(),(int) fFuncTot -> GetNDF()));

  if(fChiSquare_NDF > 2. || fMassJpsi < 3.){letexTitle -> DrawLatex(0.6,0.41,"#color[2]{BAD FIT}");}
  canvasMinv -> cd();
  TPad *padRatio = new TPad("padRatio","padRatio",0.,0.05,1.,0.25); padRatio -> SetTopMargin(0); padRatio -> SetBottomMargin(0.25);
  padRatio -> Draw();
  padRatio -> cd();
  histGridRatio -> Draw();
  axisRatio -> Draw("same");
  lineRatio -> Draw("same");
  histRatio -> Draw("Esame");

  string nameHistMinv = (string) fHistMinv -> GetName();

  if(savePlot){
    if(strstr(nameHistMinv.c_str(),"HE")){
      if(strstr(fOutputDir.c_str(),"BadFit")){canvasMinv -> SaveAs(Form("%s/%s.png",fOutputDir.c_str(),nameHistMinv.c_str()));}
      else{canvasMinv -> SaveAs(Form("%s/Helicity/%s_%s/%s.png",fOutputDir.c_str(),fSigShape.c_str(),fBkgShape.c_str(),nameHistMinv.c_str()));}
      canvasMinv -> Close();
    }
    if(strstr(nameHistMinv.c_str(),"CS")){
      if(strstr(fOutputDir.c_str(),"BadFit")){canvasMinv -> SaveAs(Form("%s/%s.png",fOutputDir.c_str(),nameHistMinv.c_str()));}
      else{canvasMinv -> SaveAs(Form("%s/Collins_Soper/%s_%s/%s.png",fOutputDir.c_str(),fSigShape.c_str(),fBkgShape.c_str(),nameHistMinv.c_str()));}
      canvasMinv -> Close();
    }
    canvasMinv -> SaveAs(Form("%s",fOutputDir.c_str()));

    delete histGridMinv, histGridRatio, histFuncTot;
    delete letexTitle;
    delete axisHistMinv, lineRatio, axisRatio;
    delete canvasMinv;
  }
}
//______________________________________________________________________________
void MassFitter_3::PlotPublications(bool savePlot){
  double max_histo_value = fHistMinv -> GetMaximum();
  max_histo_value = max_histo_value + 0.6*max_histo_value;
  fHistMinv -> SetMarkerStyle(20); fHistMinv -> SetMarkerSize(1.1); fHistMinv -> SetMarkerColor(kBlack);

  TH2D *histGridMinv = new TH2D("histGridMinv","",120,2,5,100,0,max_histo_value);
  histGridMinv -> GetYaxis() -> SetTitle(Form("Counts per %3.0f MeV/#it{c}^{2}",fHistMinv -> GetBinWidth(1)*1000));
  histGridMinv -> GetYaxis() -> CenterTitle(true);
  histGridMinv -> GetYaxis() -> SetTitleSize(0.04);
  histGridMinv -> GetYaxis() -> SetTitleOffset(1.2);
  histGridMinv -> GetYaxis() -> SetLabelSize(0);

  TGaxis *axisHistMinv = new TGaxis(2.,1.,2.,max_histo_value,1.,max_histo_value,510,"");
  axisHistMinv -> SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axisHistMinv -> SetLabelSize(35);

  TLegend *legendMassFit = new TLegend(0.6,0.35,0.77,0.5); legendMassFit -> SetTextSize(0.035); legendMassFit -> SetBorderSize(0);
  legendMassFit -> AddEntry(fHistMinv,"Data","LP");
  legendMassFit -> AddEntry(fFuncTot,"Total fit","L");
  legendMassFit -> AddEntry(fFuncBkgFix,"Background fit","L");
  legendMassFit -> AddEntry(fFuncSigJpsiFix,"Signal fit","L");


  TCanvas *canvasMinv = new TCanvas("canvasMinv","canvasMinv",65,73,900,806);
  canvasMinv -> Range(1.825,-6776.052,5.019444,37862.12);
  canvasMinv -> SetFillColor(0);
  canvasMinv -> SetBorderMode(0);
  canvasMinv -> SetBorderSize(0);
  canvasMinv -> SetTickx(1);
  canvasMinv -> SetTicky(1);
  canvasMinv -> SetLeftMargin(0.18);
  canvasMinv -> SetBottomMargin(0.1518219);
  canvasMinv -> SetFrameBorderMode(0);
  canvasMinv -> SetFrameBorderMode(0);

  canvasMinv -> cd();
  histGridMinv -> Draw();
  axisHistMinv -> Draw("same");
  fHistMinv -> Draw("sameE");
  fFuncBkgFix -> Draw("same");
  fFuncTot -> Draw("same");
  fFuncSigJpsiFix -> Draw("same");
  TLatex *latexTitle = new TLatex(); latexTitle -> SetTextSize(0.04); latexTitle -> SetNDC(); latexTitle -> SetTextFont(42);
  latexTitle -> DrawLatex(0.25,0.82,"ALICE Preliminary, Inclusive J/#psi #rightarrow #mu^{#plus}#mu^{#minus}");
  latexTitle -> DrawLatex(0.25,0.76,"Pb-Pb #sqrt{#it{s}_{NN}} = 5.02 TeV");
  latexTitle -> DrawLatex(0.25,0.7,"2 < #it{p}_{T} < 4 GeV/#it{c} , 2.5 < #it{y} < 4 , Helicity");
  if(fIndexCosTheta == 100) latexTitle -> DrawLatex(0.55,0.62,"-1. < cos#it{#theta} < 1.");
  else{latexTitle -> DrawLatex(0.55,0.62,Form("%3.2f < cos#it{#theta} < %3.2f",fCosThetaValues[fIndexCosTheta],fCosThetaValues[fIndexCosTheta+1]));}
  if(fIndexPhi == 100) latexTitle -> DrawLatex(0.55,0.55,"0 < |#it{#varphi}| < #pi rad");
  else{latexTitle -> DrawLatex(0.55,0.55,Form("%3.2f < |#it{#varphi}| < %3.2f rad",fPhiValues[fIndexPhi],fPhiValues[fIndexPhi+1]));}
  legendMassFit -> Draw("same");
}
