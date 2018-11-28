#include <Riostream.h>
#include <string>
#include <vector>

#include "TROOT.h"
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
#include "TLatex.h"
#include "TMatrixD.h"

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "Fit/Chi2FCN.h"
#include "TList.h"
#include "Math/WrappedMultiTF1.h"
#include "HFitInterface.h"
#include "TRandom.h"

#include "SpecialFitCalculator.h"

ClassImp(SpecialFitCalculator)

////////////////////////////////////////////////////////////////////////////////
// Definition of shared parameters
// Definition of FuncCosTheta parameters
int iparB1[2] = { 0, // Amplitude funcCosTheta
                  2  // lambdaTheta (common parameter) <-
};
// Definition of FuncPhi parameters
int iparB2[3] = { 1, // Amplitude funcPhi
                  2, // lambdaTheta (common parameter) <-
                  3, // lambdaPhi
};
// Definition of FuncPhiTilde parameters
int iparB3[3] = { 4, // Amplitude funcPhiTilde
                  2, // lambdaTheta (common parameter) <-
                  5, // lambdaPhi
};
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
SpecialFitCalculator::SpecialFitCalculator(): TObject() {
  fPi = TMath::Pi();
  // default constructor
}
//______________________________________________________________________________
SpecialFitCalculator::~SpecialFitCalculator() {
  // destructor
}
//______________________________________________________________________________
struct GlobalChi2 {
   GlobalChi2(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2,
                ROOT::Math::IMultiGenFunction & f3) :
      fChi2_1(&f1), fChi2_2(&f2), fChi2_3(&f3) {}

   // parameter vector is first background (in common 1 and 2)
   // and then is signal (only in 2)
   double operator() (const double *par) const {
      double p1[2];
      for (int i = 0; i < 2; ++i) p1[i] = par[iparB1[i] ];

      double p2[3];
      for (int i = 0; i < 3; ++i) p2[i] = par[iparB2[i] ];

      double p3[3];
      for (int i = 0; i < 3; ++i) p3[i] = par[iparB3[i] ];

      return (*fChi2_1)(p1) + (*fChi2_2)(p2) + (*fChi2_3)(p3);
   }

   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
   const  ROOT::Math::IMultiGenFunction * fChi2_3;
};
//______________________________________________________________________________
struct GlobalChi2_old {
   GlobalChi2_old(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2) :
      fChi2_1(&f1), fChi2_2(&f2) {}

   // parameter vector is first background (in common 1 and 2)
   // and then is signal (only in 2)
   double operator() (const double *par) const {
      double p1[2];
      for (int i = 0; i < 2; ++i) p1[i] = par[iparB1[i] ];

      double p2[3];
      for (int i = 0; i < 3; ++i) p2[i] = par[iparB2[i] ];
   }

   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
};
//______________________________________________________________________________
void SpecialFitCalculator::SimultaneousFit(TH1D *histNJpsiCost, TH1D *histNJpsiPhi, TH1D *histNJpsiPhiTilde, TH1D *histAccxEffCost, TH1D *histAccxEffPhi, TH1D *histAccxEffPhiTilde, Double_t minFitRangeCost, Double_t maxFitRangeCost, string nameOutputPlot) {
  TH1D *histNJpsiCostCorr = (TH1D*) histNJpsiCost -> Clone("histNJpsiCostCorr");
  for(int i = 0;i < 19;i++){
    histNJpsiCostCorr -> SetBinContent(i+1,(histNJpsiCost -> GetBinContent(i+1))/(histNJpsiCost -> GetBinWidth(i+1)));
    histNJpsiCostCorr -> SetBinError(i+1,(histNJpsiCost -> GetBinError(i+1))/(histNJpsiCost -> GetBinWidth(i+1)));
  }
  histNJpsiCostCorr -> Divide(histAccxEffCost);

  TH1D *histNJpsiPhiCorr = (TH1D*) histNJpsiPhi -> Clone("histNJpsiPhiCorr");
  for(int i = 0;i < 10;i++){
    histNJpsiPhiCorr -> SetBinContent(i+1,(histNJpsiPhi -> GetBinContent(i+1))/(histNJpsiPhi -> GetBinWidth(i+1)));
    histNJpsiPhiCorr -> SetBinError(i+1,(histNJpsiPhi -> GetBinError(i+1))/(histNJpsiPhi -> GetBinWidth(i+1)));
  }
  histNJpsiPhiCorr -> Divide(histAccxEffPhi);

  TH1D *histNJpsiPhiTildeCorr = (TH1D*) histNJpsiPhiTilde -> Clone("histNJpsiPhiTildeCorr");
  for(int i = 0;i < 10;i++){
    histNJpsiPhiTildeCorr -> SetBinContent(i+1,(histNJpsiPhiTilde -> GetBinContent(i+1))/(histNJpsiPhiTilde -> GetBinWidth(i+1)));
    histNJpsiPhiTildeCorr -> SetBinError(i+1,(histNJpsiPhiTilde -> GetBinError(i+1))/(histNJpsiPhiTilde -> GetBinWidth(i+1)));
  }
  histNJpsiPhiTildeCorr -> Divide(histAccxEffPhiTilde);

  TF1 *ffitB1 = new TF1("ffitB1","([0]/(3 + [1]))*(1 + [1]*x*x)",-1.,1.);
  TF1 *ffitB2 = new TF1("ffitB2","[0]*(1 + ((2*[2])/(3 + [1]))*cos(2*x))",0.,fPi);
  TF1 *ffitB3 = new TF1("ffitB3","[0]*(1 + ((sqrt(2)*[2])/(3 + [1]))*cos(x))",0.,fPi);

  ROOT::Math::WrappedMultiTF1 wfB1(*ffitB1,1);
  ROOT::Math::WrappedMultiTF1 wfB2(*ffitB2,1);
  ROOT::Math::WrappedMultiTF1 wfB3(*ffitB3,1);

  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange rangeB1;
  ROOT::Fit::DataRange rangeB2;
  ROOT::Fit::DataRange rangeB3;

  // set the data range
  rangeB1.SetRange(minFitRangeCost,maxFitRangeCost);
  rangeB2.SetRange(histNJpsiPhi -> GetBinLowEdge(2),histNJpsiPhi -> GetBinLowEdge(10));
  rangeB3.SetRange(histNJpsiPhiTilde -> GetBinLowEdge(2),histNJpsiPhiTilde -> GetBinLowEdge(10));

  ROOT::Fit::BinData dataB1(opt,rangeB1);
  ROOT::Fit::BinData dataB2(opt,rangeB2);
  ROOT::Fit::BinData dataB3(opt,rangeB3);

  ROOT::Fit::FillData(dataB1, histNJpsiCostCorr);
  ROOT::Fit::FillData(dataB2, histNJpsiPhiCorr);
  ROOT::Fit::FillData(dataB3, histNJpsiPhiTildeCorr);


  ROOT::Fit::Chi2Function chi2_B1(dataB1, wfB1);
  ROOT::Fit::Chi2Function chi2_B2(dataB2, wfB2);
  ROOT::Fit::Chi2Function chi2_B3(dataB3, wfB3);

  GlobalChi2 globalChi2(chi2_B1, chi2_B2, chi2_B3);

  ROOT::Fit::Fitter fitter;

  const int Npar = 6;
  double par0[Npar] = {10000.,10000.,-0.2,0.1,10000.,0.};

  // create before the parameter settings in order to fix or set range on them
  fitter.Config().SetParamsSettings(6,par0);

  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit","Minimize");

  // fit FCN function directly
  // (specify optionally data size and flag to indicate that is a chi2 fit)
  fitter.FitFCN(6,globalChi2,0,dataB1.Size()+dataB2.Size()+dataB3.Size(),true);
  ROOT::Fit::FitResult result = fitter.Result();
  TMatrixD matrixCov(6,6);
  result.GetCovarianceMatrix(matrixCov);
  matrixCov.Print();
  fErrorLambdaThetaLambdaPhi = matrixCov(2,3);

  result.Print(std::cout);

  fNormCosTheta = result.Parameter(0);
  fErrorNormCosTheta = result.ParError(0);
  fNormPhi = result.Parameter(1);
  fErrorNormPhi = result.ParError(1);
  fLambdaTheta = result.Parameter(2);
  fErrorLambdaTheta = result.ParError(2);
  fLambdaPhi = result.Parameter(3);
  fErrorLambdaPhi = result.ParError(3);
  fNormPhiTilde = result.Parameter(4);
  fErrorNormPhiTilde = result.ParError(4);
  fLambdaThetaPhi = result.Parameter(5);
  fErrorLambdaThetaPhi = result.ParError(5);

  fCosThetaParametersList.push_back(fNormCosTheta);
  fErrorCosThetaParametersList.push_back(fErrorNormCosTheta);
  fCosThetaParametersList.push_back(fLambdaTheta);
  fErrorCosThetaParametersList.push_back(fErrorLambdaTheta);

  fPhiParametersList.push_back(fNormPhi);
  fErrorPhiParametersList.push_back(fErrorNormPhi);
  fPhiParametersList.push_back(fLambdaTheta);
  fErrorPhiParametersList.push_back(fErrorLambdaTheta);
  fPhiParametersList.push_back(fLambdaPhi);
  fErrorPhiParametersList.push_back(fErrorLambdaPhi);

  gStyle->SetOptFit(1111);

  TCanvas * canvas = new TCanvas("Simfit","Simultaneous fit of two graphs",10,10,1400,700);
  canvas -> Divide(3,1);

  canvas -> cd(1);
  ffitB1 -> SetFitResult(result, iparB1);
  ffitB1 -> SetRange(rangeB1().first, rangeB1().second);
  ffitB1 -> SetLineColor(kBlue);
  histNJpsiCostCorr -> GetListOfFunctions() -> Add(ffitB1);
  histNJpsiCostCorr -> Draw("E");

  canvas -> cd(2);
  ffitB2 -> SetFitResult( result, iparB2);
  ffitB2 -> SetRange(rangeB2().first, rangeB2().second);
  ffitB2 -> SetLineColor(kRed);
  histNJpsiPhiCorr -> GetListOfFunctions()->Add(ffitB2);
  histNJpsiPhiCorr -> Draw("E");

  canvas -> cd(3);
  ffitB3 -> SetFitResult( result, iparB3);
  ffitB3 -> SetRange(rangeB3().first, rangeB3().second);
  ffitB3 -> SetLineColor(kRed);
  histNJpsiPhiTildeCorr -> GetListOfFunctions()->Add(ffitB3);
  histNJpsiPhiTildeCorr -> Draw("E");

  canvas -> SaveAs(nameOutputPlot.c_str());
  delete canvas;
}
//______________________________________________________________________________
struct GlobalChi2_2Var {
   GlobalChi2_2Var(  ROOT::Math::IMultiGenFunction & f1,
                ROOT::Math::IMultiGenFunction & f2) :
      fChi2_1(&f1), fChi2_2(&f2) {}

   // parameter vector is first background (in common 1 and 2)
   // and then is signal (only in 2)
   double operator() (const double *par) const {
      double p1[2];
      for (int i = 0; i < 2; ++i) p1[i] = par[iparB1[i] ];

      double p2[3];
      for (int i = 0; i < 3; ++i) p2[i] = par[iparB2[i] ];
   }

   const  ROOT::Math::IMultiGenFunction * fChi2_1;
   const  ROOT::Math::IMultiGenFunction * fChi2_2;
};
//______________________________________________________________________________
void SpecialFitCalculator::SimultaneousFit2Var(TH1D *histNJpsiCost, TH1D *histNJpsiPhi, TH1D *histAccxEffCost, TH1D *histAccxEffPhi, Double_t minFitRangeCost, Double_t maxFitRangeCost, string nameOutputPlot) {
  TH1D *histNJpsiCostCorr = (TH1D*) histNJpsiCost -> Clone("histNJpsiCostCorr");
  for(int i = 0;i < 19;i++){
    histNJpsiCostCorr -> SetBinContent(i+1,(histNJpsiCost -> GetBinContent(i+1))/(histNJpsiCost -> GetBinWidth(i+1)));
    histNJpsiCostCorr -> SetBinError(i+1,(histNJpsiCost -> GetBinError(i+1))/(histNJpsiCost -> GetBinWidth(i+1)));
  }
  histNJpsiCostCorr -> Divide(histAccxEffCost);

  TH1D *histNJpsiPhiCorr = (TH1D*) histNJpsiPhi -> Clone("histNJpsiPhiCorr");
  for(int i = 0;i < 10;i++){
    histNJpsiPhiCorr -> SetBinContent(i+1,(histNJpsiPhi -> GetBinContent(i+1))/(histNJpsiPhi -> GetBinWidth(i+1)));
    histNJpsiPhiCorr -> SetBinError(i+1,(histNJpsiPhi -> GetBinError(i+1))/(histNJpsiPhi -> GetBinWidth(i+1)));
  }
  histNJpsiPhiCorr -> Divide(histAccxEffPhi);

  TF1 *ffitB1 = new TF1("ffitB1","([0]/(3 + [1]))*(1 + [1]*x*x)",-1.,1.);
  TF1 *ffitB2 = new TF1("ffitB2","[0]*(1 + ((2*[2])/(3 + [1]))*cos(2*x))",0.,fPi);

  ROOT::Math::WrappedMultiTF1 wfB1(*ffitB1,1);
  ROOT::Math::WrappedMultiTF1 wfB2(*ffitB2,1);

  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange rangeB1;
  ROOT::Fit::DataRange rangeB2;

  // set the data range
  rangeB1.SetRange(minFitRangeCost,maxFitRangeCost);
  rangeB2.SetRange(histNJpsiPhi -> GetBinLowEdge(2),histNJpsiPhi -> GetBinLowEdge(10));

  ROOT::Fit::BinData dataB1(opt,rangeB1);
  ROOT::Fit::BinData dataB2(opt,rangeB2);

  ROOT::Fit::FillData(dataB1, histNJpsiCostCorr);
  ROOT::Fit::FillData(dataB2, histNJpsiPhiCorr);


  ROOT::Fit::Chi2Function chi2_B1(dataB1, wfB1);
  ROOT::Fit::Chi2Function chi2_B2(dataB2, wfB2);

  GlobalChi2_2Var globalChi2(chi2_B1, chi2_B2);

  ROOT::Fit::Fitter fitter;

  const int Npar = 4;
  double par0[Npar] = {10000.,10000.,-0.2,0.1};

  // create before the parameter settings in order to fix or set range on them
  fitter.Config().SetParamsSettings(4,par0);

  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit","Minimize");

  // fit FCN function directly
  // (specify optionally data size and flag to indicate that is a chi2 fit)
  fitter.FitFCN(4,globalChi2,0,dataB1.Size()+dataB2.Size(),true);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);

  fNormCosTheta = result.Parameter(0);
  fErrorNormCosTheta = result.ParError(0);
  fNormPhi = result.Parameter(1);
  fErrorNormPhi = result.ParError(1);
  fLambdaTheta = result.Parameter(2);
  fErrorLambdaTheta = result.ParError(2);
  fLambdaPhi = result.Parameter(3);
  fErrorLambdaPhi = result.ParError(3);

  fCosThetaParametersList.push_back(fNormCosTheta);
  fErrorCosThetaParametersList.push_back(fErrorNormCosTheta);
  fCosThetaParametersList.push_back(fLambdaTheta);
  fErrorCosThetaParametersList.push_back(fErrorLambdaTheta);

  fPhiParametersList.push_back(fNormPhi);
  fErrorPhiParametersList.push_back(fErrorNormPhi);
  fPhiParametersList.push_back(fLambdaTheta);
  fErrorPhiParametersList.push_back(fErrorLambdaTheta);
  fPhiParametersList.push_back(fLambdaPhi);
  fErrorPhiParametersList.push_back(fErrorLambdaPhi);

  gStyle->SetOptFit(1111);

  TCanvas * canvas = new TCanvas("Simfit","Simultaneous fit of two graphs",10,10,1400,700);
  canvas -> Divide(2,1);

  canvas -> cd(1);
  ffitB1 -> SetFitResult(result, iparB1);
  ffitB1 -> SetRange(rangeB1().first, rangeB1().second);
  ffitB1 -> SetLineColor(kBlue);
  histNJpsiCostCorr -> GetListOfFunctions() -> Add(ffitB1);
  histNJpsiCostCorr -> Draw("E");

  canvas -> cd(2);
  ffitB2 -> SetFitResult( result, iparB2);
  ffitB2 -> SetRange(rangeB2().first, rangeB2().second);
  ffitB2 -> SetLineColor(kRed);
  histNJpsiPhiCorr -> GetListOfFunctions()->Add(ffitB2);
  histNJpsiPhiCorr -> Draw("E");

  canvas -> SaveAs(nameOutputPlot.c_str());
  delete canvas;
}
//______________________________________________________________________________
void SpecialFitCalculator::SimultaneousFit2(TH1D *histNJpsiCostCorr, TH1D *histNJpsiPhiCorr, Double_t minFitRangeCost, Double_t maxFitRangeCost, string nameOutputPlot) {
  TF1 *ffitB1 = new TF1("ffitB1","([0]/(3 + [1]))*(1 + [1]*x*x)",-1.,1.);
  TF1 *ffitB2 = new TF1("ffitB2","[0]*(1 + ((2*[2])/(3 + [1]))*cos(2*x))",0.,fPi);

  ROOT::Math::WrappedMultiTF1 wfB1(*ffitB1,1);
  ROOT::Math::WrappedMultiTF1 wfB2(*ffitB2,1);

  ROOT::Fit::DataOptions opt;
  ROOT::Fit::DataRange rangeB1;
  ROOT::Fit::DataRange rangeB2;

  // set the data range
  rangeB1.SetRange(minFitRangeCost,maxFitRangeCost);
  rangeB2.SetRange(histNJpsiPhiCorr -> GetBinLowEdge(2),histNJpsiPhiCorr -> GetBinLowEdge(10));

  ROOT::Fit::BinData dataB1(opt,rangeB1);
  ROOT::Fit::BinData dataB2(opt,rangeB2);

  ROOT::Fit::FillData(dataB1, histNJpsiCostCorr);
  ROOT::Fit::FillData(dataB2, histNJpsiPhiCorr);


  ROOT::Fit::Chi2Function chi2_B1(dataB1, wfB1);
  ROOT::Fit::Chi2Function chi2_B2(dataB2, wfB2);

  GlobalChi2_old globalChi2_old(chi2_B1, chi2_B2);

  ROOT::Fit::Fitter fitter;

  const int Npar = 4;
  double par0[Npar] = {10000.,10000.,-0.2,0.1};

  // create before the parameter settings in order to fix or set range on them
  fitter.Config().SetParamsSettings(4,par0);

  fitter.Config().MinimizerOptions().SetPrintLevel(0);
  fitter.Config().SetMinimizer("Minuit","Minimize");

  // fit FCN function directly
  // (specify optionally data size and flag to indicate that is a chi2 fit)
  fitter.FitFCN(4,globalChi2_old,0,dataB1.Size()+dataB2.Size(),true);
  ROOT::Fit::FitResult result = fitter.Result();
  result.Print(std::cout);

  fNormCosTheta = result.Parameter(0);
  fErrorNormCosTheta = result.ParError(0);
  fNormPhi = result.Parameter(1);
  fErrorNormPhi = result.ParError(1);
  fLambdaTheta = result.Parameter(2);
  fErrorLambdaTheta = result.ParError(2);
  fLambdaPhi = result.Parameter(3);
  fErrorLambdaPhi = result.ParError(3);

  fCosThetaParametersList.push_back(fNormCosTheta);
  fErrorCosThetaParametersList.push_back(fErrorNormCosTheta);
  fCosThetaParametersList.push_back(fLambdaTheta);
  fErrorCosThetaParametersList.push_back(fErrorLambdaTheta);

  fPhiParametersList.push_back(fNormPhi);
  fErrorPhiParametersList.push_back(fErrorNormPhi);
  fPhiParametersList.push_back(fLambdaTheta);
  fErrorPhiParametersList.push_back(fErrorLambdaTheta);
  fPhiParametersList.push_back(fLambdaPhi);
  fErrorPhiParametersList.push_back(fErrorLambdaPhi);

  TCanvas * canvas = new TCanvas("Simfit","Simultaneous fit of two graphs",10,10,700,700);
  canvas -> Divide(1,2);

  char title[100];
  sprintf(title,"#lambda_{#theta} = %f #pm %f",fLambdaTheta,fErrorLambdaTheta);
  TLatex *lat = new TLatex(0.5,0.82,title);
  lat -> SetTextSize(0.04);
  lat -> SetNDC();
  lat -> SetTextFont(42);

  canvas -> cd(1);
  ffitB1 -> SetFitResult(result, iparB1);
  ffitB1 -> SetRange(rangeB1().first, rangeB1().second);
  ffitB1 -> SetLineColor(kBlue);
  histNJpsiCostCorr -> GetListOfFunctions() -> Add(ffitB1);
  histNJpsiCostCorr -> Draw("E");
  lat -> Draw("same");

  canvas -> cd(2);
  ffitB2 -> SetFitResult( result, iparB2);
  ffitB2 -> SetRange(rangeB2().first, rangeB2().second);
  ffitB2 -> SetLineColor(kRed);
  histNJpsiPhiCorr -> GetListOfFunctions()->Add(ffitB2);
  histNJpsiPhiCorr -> Draw("E");

  canvas -> SaveAs(nameOutputPlot.c_str());
  delete canvas;
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
