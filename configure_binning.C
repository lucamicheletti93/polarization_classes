#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>
#include <string>
#include <vector>

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
#include <TGaxis.h>

#include "Binning.h"
#include "Binning.cxx"
#endif

void configure_binning(){

  string pathOutput = "output/";

  const int nCostBins = 19;
  int minCostBin[nCostBins] = {1,11,16,21,26,31,36,41,45,49,53,57,61,66,71,76,81,86,91};
  int maxCostBin[nCostBins] = {10,15,20,25,30,35,40,44,48,52,56,60,65,70,75,80,85,90,100};

  const int nPhiBins = 10;
  int minPhiBin[nPhiBins] = {1,9,17,21,24,26,28,31,35,43};
  int maxPhiBin[nPhiBins] = {8,16,20,23,25,27,30,34,42,50};

  Binning *binning = new Binning();
  binning -> ConfigureBinValues(nCostBins,minCostBin,maxCostBin,nPhiBins,minPhiBin,maxPhiBin);

  TFile *fileBinning = new TFile(Form("%sbinning.root",pathOutput.c_str()),"RECREATE");
  binning -> Write();
  fileBinning -> Close();
}
