#include <Riostream.h>
/*#include <string>
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
#include "MassFitter.h"*/
#include <stdio.h>
#include <string>
#include <vector>
#include <sstream>

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
#include "/home/luca/GITHUB/polarization/2D_approach/signal_extraction/FunctionsLibrary.C"
#include "MassFitter.h"

// List of fitting-functions
Double_t Func_VWG(Double_t *, Double_t *);
Double_t Func_POL4EXP(Double_t *, Double_t *);
Double_t Func_Jpsi_CB2(Double_t *, Double_t *);
Double_t Func_Jpsi_CB2_POL4EXP(Double_t *, Double_t *);
Double_t Func_Jpsi_CB2_fix(Double_t *, Double_t *);
Double_t Func_Jpsi_NA60_VWG(Double_t *, Double_t *);
Double_t Func_Jpsi_NA60_POL4EXP(Double_t *, Double_t *);
Double_t Func_Jpsi_NA60_fix(Double_t *, Double_t *);
Double_t Func_tot(Double_t *, Double_t *);
Double_t Func_tot_CB2_POL4EXP(Double_t *, Double_t *);
Double_t Func_tot_NA60_VWG(Double_t *, Double_t *);
Double_t Func_tot_NA60_POL4EXP(Double_t *, Double_t *);

// List of global variables
double gPi = TMath::Pi();
Long_t *dummy1 = 0, *dummy2 = 0, *dummy3 = 0, *dummy4 = 0;

ClassImp(MassFitter)

//______________________________________________________________________________
MassFitter::MassFitter(): TObject() {
  //fPi = TMath::Pi();
  // default constructor
}
//______________________________________________________________________________
MassFitter::~MassFitter() {
  // destructor
}
//______________________________________________________________________________
void MassFitter::SetScalingFactor(Double_t scalingFactor){
    fScalingFactorJpsiSigma = scalingFactor;
}
//______________________________________________________________________________
void MassFitter::SetBinning(vector <Double_t> CosThetaValues, vector <Double_t> PhiValues, vector <Double_t> PhiTildeValues) {
  for(int i = 0;i < (int) CosThetaValues.size();i++){fCosThetaValues.push_back(CosThetaValues[i]);}
  fNCosThetaBins = fCosThetaValues.size() - 1;
  for(int i = 0;i < (int) PhiValues.size();i++){fPhiValues.push_back(PhiValues[i]);}
  fNPhiBins = fPhiValues.size() - 1;
  for(int i = 0;i < (int) PhiTildeValues.size();i++){fPhiTildeValues.push_back(PhiTildeValues[i]);}
  fNPhiTildeBins = fPhiTildeValues.size() - 1;
}
//______________________________________________________________________________
void MassFitter::fit_of_minv(TH1D *histMinv, string sigShape, string bkgShape, int counterCosTheta, int counterPhi, string dataset, string outputDir, TH2D *histSigmaJpsiMC, bool specialFitConditions){
  gStyle -> SetOptStat(0);
  TGaxis::SetMaxDigits(2);

  double min_fit_range = 2.1;
  double max_fit_range = 4.9;
  string tails_fix = "yes";

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

  double parTailsCB2[4] = {0.970014,3.9789,2.29866,3.0301};                                                 // tails parameter embedding Pb-Pb 5.02 TeV (0 < pT < 12 GeV/c)
  double parTailsNA60[8] = {0.227822,1.14006,0.0357395,0.187828,1.22477,0.0569524,-0.630315,2.36889};       // tails parameter embedding Pb-Pb 5.02 TeV (0 < pT < 12 GeV/c)
  double *parTails;
  if(sigShape == "CB2"){parTails = &(parTailsCB2[0]); nParTails = sizeof(parTailsCB2)/sizeof(double);}
  if(sigShape == "NA60"){parTails = &(parTailsNA60[0]); nParTails = sizeof(parTailsNA60)/sizeof(double);}

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
  TH1D *histMinvBkg = (TH1D*) histMinv -> Clone("histMinvBkg");
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

  if(specialFitConditions){
    histMinv -> Rebin(2);
    min_fit_range = 2.;
    max_fit_range = 5.;
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
  if(sigShape == "CB2" && bkgShape == "VWG"){funcSigJpsi = new TF1("funcSigJpsi",Func_Jpsi_CB2,2.9,3.2,nParBkg + nParSig + nParTails);}
  if(sigShape == "NA60" && bkgShape == "VWG"){funcSigJpsi = new TF1("funcSigJpsi",Func_Jpsi_NA60_VWG,2.9,3.2,nParBkg + nParSig + nParTails);}
  if(sigShape == "CB2" && bkgShape == "POL4EXP"){funcSigJpsi = new TF1("funcSigJpsi",Func_Jpsi_CB2_POL4EXP,2.9,3.2,nParBkg + nParSig + nParTails);}
  if(sigShape == "NA60" && bkgShape == "POL4EXP"){funcSigJpsi = new TF1("funcSigJpsi",Func_Jpsi_NA60_POL4EXP,2.9,3.2,nParBkg + nParSig + nParTails);}
  funcSigJpsi -> SetNpx(100000);
  funcSigJpsi -> SetParameter(nParBkg + 0,parSig[0]);
  funcSigJpsi -> SetParameter(nParBkg + 1,parSig[1]);
  funcSigJpsi -> FixParameter(nParBkg + 2,parSig[2]);
  for(int i = 0;i < nParTails;i++){funcSigJpsi -> FixParameter(nParBkg + nParSig + i,parTails[i]);}
  histMinv -> Fit(funcSigJpsi,"R0Q");

  //============================================================================
  // FIT OF THE TOTAL SPECTRUM
  //============================================================================
  TF1 *funcTot;
  if(sigShape == "CB2" && bkgShape == "VWG"){funcTot = new TF1("funcTot",Func_tot,min_fit_range,max_fit_range,nParBkg + nParSig + nParTails);}
  if(sigShape == "NA60" && bkgShape == "VWG"){funcTot = new TF1("funcTot",Func_tot_NA60_VWG,min_fit_range,max_fit_range,nParBkg + nParSig + nParTails);}
  if(sigShape == "CB2" && bkgShape == "POL4EXP"){funcTot = new TF1("funcTot",Func_tot_CB2_POL4EXP,min_fit_range,max_fit_range,nParBkg + nParSig + nParTails);}
  if(sigShape == "NA60" && bkgShape == "POL4EXP"){funcTot = new TF1("funcTot",Func_tot_NA60_POL4EXP,min_fit_range,max_fit_range,nParBkg + nParSig + nParTails);}

  TFitResultPtr fit_ptr;
  for(int i = 0;i < 100;i++){
    if(i == 0){
      for(int i = 0;i < nParBkg;i++){funcTot -> SetParameter(i,funcBkg -> GetParameter(i));}

      funcTot -> SetParameter(nParBkg + 0,funcSigJpsi -> GetParameter(nParBkg));
      funcTot -> SetParLimits(nParBkg + 0,0,10000000);
      if(specialFitConditions){
        funcTot -> FixParameter(nParBkg + 1,3.096);
      }
      else{
        funcTot -> SetParameter(nParBkg + 1,3.096);
        funcTot -> SetParLimits(nParBkg + 1,2.9,3.2);
      }
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // J/psi width conditions
      if(counterCosTheta == 100 && counterPhi == 100){
        funcTot -> SetParameter(nParBkg + 2,7.0e-02);
        funcTot -> SetParLimits(nParBkg + 2,2.0e-02,2.0e-01);}
      else{
        double fSigmaJpsiFixed = histSigmaJpsiMC -> GetBinContent(counterCosTheta+1,counterPhi+1);
        cout << fSigmaJpsiFixed << endl;
        fSigmaJpsiFixed = fScalingFactorJpsiSigma*fSigmaJpsiFixed;
        funcTot -> FixParameter(nParBkg + 2,fSigmaJpsiFixed);
      }
      //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if(tails_fix == "yes"){for(int i = 0;i < nParTails;i++){funcTot -> FixParameter(nParBkg + nParSig + i,funcSigJpsi -> GetParameter(nParBkg + nParSig + i));}}
      if(tails_fix == "no"){for(int i = 0;i < nParTails;i++){funcTot -> SetParameter(nParBkg + nParSig + i,funcSigJpsi -> GetParameter(nParBkg + nParSig + i));}}
    }
    else{funcTot -> SetParameters(funcTot -> GetParameters());}
    fit_ptr = (TFitResultPtr) histMinv -> Fit(funcTot,"RLS0Q");
    if(gMinuit->fCstatu.Contains("CONVERGED")) break;
  }

  fChiSquare_NDF = funcTot -> GetChisquare()/funcTot -> GetNDF();

  //printf("\n\nfit status: %s \n\n",gMinuit.fCstatu.Data());
  if(gMinuit -> fCstatu.Contains("FAILED")){
    cout << "WARNING : FIT STATUS FAILED" << endl;
    fFitStatus = "FAILED";
    delete histMinvBkg, histMinv;
    delete funcBkg, funcSigJpsi, funcTot;
    fNJpsi = 0.; fStatJpsi = 0.; fSigmaJpsi = 0.; fErrSigmaJpsi = 0.; fChiSquare_NDF = 666.;
    string nameHistMinv = (string) histMinv -> GetName() + ".png";

    TCanvas *canvasMinv = new TCanvas("canvasMinv","canvasMinv",65,73,900,806);
    histMinv -> Draw();

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

  // Getting sub matrix for Jpsi (error calculation)
  TMatrixDSym cov = fit_ptr -> GetCovarianceMatrix().GetSub(dimParBkg,(dimParBkg+dimParSig)-1,dimParBkg,(dimParBkg+dimParSig)-1);
  double Jpsi_par[dimParSig];
  for(int i = 0;i < dimParSig;i++){Jpsi_par[i] = funcTot -> GetParameter(dimParBkg+i);}

  //============================================================================
  // PLOT OF JPSI AND PSI(2S) SHAPES
  //============================================================================
  TF1 *funcBkgFix;
  if(bkgShape == "VWG"){funcBkgFix = new TF1("funcBkgFix",Func_VWG,2.,5.,nParBkg);}
  if(bkgShape == "POL4EXP"){funcBkgFix = new TF1("funcBkgFix",Func_POL4EXP,2.,5.,nParBkg);}
  funcBkgFix -> SetNpx(1000);
  for(int i = 0;i < nParBkg;i++){funcBkgFix -> SetParameter(i,funcTot -> GetParameter(i));}
  funcBkgFix -> SetLineStyle(4);
  funcBkgFix -> SetLineColor(kBlue+1);
  funcBkgFix -> Draw("same");

  TF1 *funcSigJpsiFix;
  if(sigShape == "CB2"){funcSigJpsiFix = new TF1("funcSigJpsiFix",Func_Jpsi_CB2_fix,2.,5.,nParSig + nParTails);}
  if(sigShape == "NA60"){funcSigJpsiFix = new TF1("funcSigJpsiFix",Func_Jpsi_NA60_fix,2.,5.,nParSig + nParTails);}
  funcSigJpsiFix -> SetNpx(1000);
  for(int i = 0;i < nParSig + nParTails;i++){funcSigJpsiFix -> SetParameter(i,funcTot -> GetParameter(nParBkg + i));}
  funcSigJpsiFix -> SetLineStyle(2);
  funcSigJpsiFix -> Draw("same");

  double mass_Jpsi = funcTot -> GetParameter(nParBkg + 1);
  fSigmaJpsi = funcTot -> GetParameter(nParBkg + 2);
  fErrSigmaJpsi = funcTot -> GetParError(nParBkg + 2);
  double sigma_min_Jpsi = funcTot -> GetParameter(nParBkg + 1) - 3*(funcTot -> GetParameter(nParBkg + 2));
  double sigma_max_Jpsi = funcTot -> GetParameter(nParBkg + 1) + 3*(funcTot -> GetParameter(nParBkg + 2));
  double N_Jpsi_3sigma = funcSigJpsiFix -> Integral(sigma_min_Jpsi,sigma_max_Jpsi)/m_width;
  double N_bck_Jpsi_3sigma = funcBkgFix -> Integral(sigma_min_Jpsi,sigma_max_Jpsi)/m_width;
  double SB_Jpsi = N_Jpsi_3sigma/N_bck_Jpsi_3sigma;

  fNJpsi = funcSigJpsiFix -> Integral(0,5)/m_width;
  fStatJpsi = funcSigJpsiFix -> IntegralError(0.,5.,Jpsi_par,cov.GetMatrixArray())/m_width;

  //============================================================================
  // PLOT OF THE TOTAL SPECTRUM
  //============================================================================
  char title[100];
  double max_histo_value = histMinv -> GetMaximum();
  max_histo_value = max_histo_value + 0.6*max_histo_value;

  TH2D *histGridMinv = new TH2D("histGridMinv","",120,2,5,100,0,max_histo_value);
  histGridMinv -> GetYaxis() -> SetTitle(Form("Counts per %3.0f MeV/#it{c}^{2}",histMinv -> GetBinWidth(1)*1000));
  histGridMinv -> GetYaxis() -> CenterTitle(true);
  histGridMinv -> GetYaxis() -> SetTitleSize(0.05);
  histGridMinv -> GetYaxis() -> SetTitleOffset(0.95);
  histGridMinv -> GetYaxis() -> SetLabelSize(0);

  TGaxis *axisHistMinv = new TGaxis(2.,1.,2.,max_histo_value,1.,max_histo_value,510,"");
  axisHistMinv->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axisHistMinv->SetLabelSize(25);

  TLine *lineRatio = new TLine(2.,1.,5.,1.0); lineRatio -> SetLineStyle(kDashed); lineRatio -> SetLineWidth(2);

  TH2D *histGridRatio = new TH2D("histGridRatio","",120,2,5,10,0.5,1.5);
  histGridRatio -> GetXaxis() -> SetTitle("#it{m}_{#mu^{#plus}#mu^{#minus}} (GeV/#it{c}^{2})");
  histGridRatio -> GetXaxis() -> SetTitleSize(0.12);
  histGridRatio -> GetXaxis() -> SetTitleOffset(0.9);
  histGridRatio -> GetXaxis() -> SetLabelSize(0.12);
  histGridRatio -> GetYaxis() -> SetTitle("Data/Fit");
  histGridRatio -> GetYaxis() -> CenterTitle(true);
  histGridRatio -> GetYaxis() -> SetTitleSize(0.1);
  histGridRatio -> GetYaxis() -> SetTitleOffset(0.4);
  histGridRatio -> GetYaxis() -> SetLabelSize(0);

  TGaxis *axisRatio = new TGaxis(2.,0.55,2.,1.45,0.55,1.45,510,"");
  axisRatio->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  axisRatio->SetLabelSize(15);

  sprintf(title,"N_{J/#psi} = %2.0f #pm %2.0f",fNJpsi,fStatJpsi);
  TLatex *lat0 = new TLatex(0.5,0.82,title); lat0 -> SetTextSize(0.04); lat0 -> SetNDC(); lat0 -> SetTextFont(42);

  sprintf(title,"#it{m}_{J/#psi} = %4.3f GeV/#it{c}^{2}",mass_Jpsi);
  TLatex *lat1 = new TLatex(0.5,0.76,title); lat1 -> SetTextSize(0.04); lat1 -> SetNDC(); lat1 -> SetTextFont(42);

  sprintf(title,"#it{#sigma}_{J/#psi} = %3.0f #pm %3.0f MeV/#it{c}^{2}",fSigmaJpsi*1000,fErrSigmaJpsi*1000);
  TLatex *lat2 = new TLatex(0.5,0.70,title); lat2 -> SetTextSize(0.04); lat2 -> SetNDC(); lat2 -> SetTextFont(42);

  if(counterCosTheta == 100) sprintf(title,"%2.1f < cos#it{#theta} < %2.1f",-1.,1.);
  else{sprintf(title,"%3.2f < cos#it{#theta} < %3.2f",fCosThetaValues[counterCosTheta],fCosThetaValues[counterCosTheta+1]);}
  TLatex *lat3 = new TLatex(0.55,0.62,title);
  lat3 -> SetTextSize(0.04); lat3 -> SetNDC(); lat3 -> SetTextFont(42);

  if(counterPhi == 100) sprintf(title,"0 < |#it{#varphi}| < #pi rad");
  else{sprintf(title,"%3.2f < |#it{#varphi}| < %3.2f rad",fPhiValues[counterPhi],fPhiValues[counterPhi+1]);}
  TLatex *lat4 = new TLatex(0.55,0.55,title);
  lat4 -> SetTextSize(0.04); lat4 -> SetNDC(); lat4 -> SetTextFont(42);

  sprintf(title,"#chi^{2}/ndf = %3.1f",fChiSquare_NDF);
  TLatex *lat5 = new TLatex(0.55,0.48,title);
  lat5 -> SetTextSize(0.04); lat5 -> SetNDC(); lat5 -> SetTextFont(42);

  //============================================================================
  // DRAVING THE TH2D
  //============================================================================
  TH1D *histFuncTot = new TH1D("histFuncTot","",120,2.,5.);
  if(specialFitConditions){histFuncTot -> Rebin(2);}
  for(int i = 0;i < 120;i++){
    histFuncTot -> SetBinContent(i+1,funcTot -> Eval(histFuncTot -> GetBinCenter(i+1)));
    histFuncTot -> SetBinError(i+1,0.);
  }

  TH1D *histRatio = (TH1D*) histMinv -> Clone("histRatio");
  histRatio -> Divide(histFuncTot);

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
  histMinv -> Draw("sameE");
  funcBkgFix -> Draw("same");
  funcTot -> Draw("same");
  funcSigJpsiFix -> Draw("same");
  lat0 -> Draw();
  lat1 -> Draw();
  lat2 -> Draw();
  lat3 -> Draw();
  lat4 -> Draw();
  lat5 -> Draw();

  canvasMinv -> cd();
  TPad *padRatio = new TPad("padRatio","padRatio",0.,0.05,1.,0.25); padRatio -> SetTopMargin(0); padRatio -> SetBottomMargin(0.25);
  padRatio -> Draw();
  padRatio -> cd();
  histGridRatio -> Draw();
  axisRatio -> Draw("same");
  lineRatio -> Draw("same");
  histRatio -> Draw("Esame");

  string nameHistMinv = (string) histMinv -> GetName() + ".png";

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

  delete histGridMinv, histGridRatio, histFuncTot;
  delete histMinvBkg, histMinv;
  delete funcBkg, funcSigJpsi, funcTot, funcBkgFix, funcSigJpsiFix;
  delete lat0, lat1, lat2, lat3, lat4, lat5;
  delete axisHistMinv, lineRatio, axisRatio, padHistMinv, padRatio;
  delete canvasMinv;
}