#include <Riostream.h>
#include <string>
#include <vector>

#include "TObject.h"
#include "TSystem.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "Binning.h"

ClassImp(Binning)

//______________________________________________________________________________
Binning::Binning(): TObject() {
   // default constructor
   double gPi = TMath::Pi();
   fCosThetaMin = -1.;
   fCosThetaMax = 1.;
   fPhiMin = 0.;
   fPhiMax = gPi;
   fPhiTildeMin = 0.;
   fPhiTildeMax = 2*gPi;
}
//______________________________________________________________________________
Binning::~Binning() {
  // destructor
}
//______________________________________________________________________________
void Binning::ConfigureBinValues(int NCostBins, int CostBinsMin[], int CostBinsMax[], int NPhiBins, int PhiBinsMin[], int PhiBinsMax[], int NPhiTildeBins, int PhiTildeBinsMin[], int PhiTildeBinsMax[]) {

  double gPi = TMath::Pi();
  double bin;

  printf("_______________Range Configuration________________ \n");
  printf("%f < CosTheta < %f \n",fCosThetaMin,fCosThetaMax);
  printf("%f < |Phi| < %f \n",fPhiMin,fPhiMax);
  printf("%f < PhiTilde < %f \n",fPhiTildeMin,fPhiTildeMax);
  printf("__________________________________________________ \n");
  //TH2D *hist = new TH2D("hist","hist",100,-1.,1.,50,0.,gPi);
  TH2D *hist = new TH2D("hist","hist",100,fCosThetaMin,fCosThetaMax,50,fPhiMin,fPhiMax);

  fNCostBins = NCostBins;
  for(int i = 0;i < NCostBins;i++){
    bin = hist -> GetXaxis() -> GetBinLowEdge(CostBinsMin[i]);
    fCostValues.push_back(bin);
    fCostBinsMin.push_back(CostBinsMin[i]);
    fCostBinsMax.push_back(CostBinsMax[i]);
    //printf("CosTheta = %f \n",bin);
  }
  bin = hist -> GetXaxis() -> GetBinLowEdge(101);
  fCostValues.push_back(bin);
  //printf("CosTheta = %f \n",bin);

  fNPhiBins = NPhiBins;
  for(int i = 0;i < NPhiBins;i++){
    bin = hist -> GetYaxis() -> GetBinLowEdge(PhiBinsMin[i]);
    fPhiValues.push_back(bin);
    fPhiBinsMin.push_back(PhiBinsMin[i]);
    fPhiBinsMax.push_back(PhiBinsMax[i]);
    //printf("Phi = %f \n",bin);
  }
  bin = hist -> GetYaxis() -> GetBinLowEdge(51);
  fPhiValues.push_back(bin);
  //printf("Phi = %f \n",bin);
  delete hist;

  //TH1D *histPhiTilde = new TH1D("histPhiTilde","histPhiTilde",50,0.,gPi);    // old binning
  //TH1D *histPhiTilde = new TH1D("histPhiTilde","histPhiTilde",50,0.,2*gPi);      // new binning
  TH1D *histPhiTilde = new TH1D("histPhiTilde","histPhiTilde",50,fPhiTildeMin,fPhiTildeMax);      // new binning

  fNPhiTildeBins = NPhiTildeBins;
  for(int i = 0;i < NPhiTildeBins;i++){
    bin = histPhiTilde -> GetXaxis() -> GetBinLowEdge(PhiTildeBinsMin[i]);
    fPhiTildeValues.push_back(bin);
    fPhiTildeBinsMin.push_back(PhiTildeBinsMin[i]);
    fPhiTildeBinsMax.push_back(PhiTildeBinsMax[i]);
    //printf("PhiTilde = %f \n",bin);
  }
  bin = histPhiTilde -> GetXaxis() -> GetBinLowEdge(51);
  //printf("PhiTilde = %f \n",bin);
  fPhiTildeValues.push_back(bin);
  delete histPhiTilde;
}
//______________________________________________________________________________
void Binning::PrintBinValues() {
  double gPi = TMath::Pi();
  printf("CostValues[%i] = {",fNCostBins+1);
  for(int i = 0;i < fNCostBins-1;i++) printf("%3.2f,",fCostValues[i]);
  printf("%3.2f,%3.2f}\n",fCostValues[fNCostBins-1],1.);

  printf("PhiValues[%i] = {",fNPhiBins+1);
  for(int i = 0;i < fNPhiBins-1;i++) printf("%f,",fPhiValues[i]);
  printf("%f,%f}\n",fPhiValues[fNPhiBins-1],gPi);

  printf("PhiTildeValues[%i] = {",fNPhiTildeBins+1);
  for(int i = 0;i < fNPhiTildeBins-1;i++) printf("%f,",fPhiTildeValues[i]);
  printf("%f,%f}\n",fPhiTildeValues[fNPhiTildeBins-1],2*gPi);
}
//______________________________________________________________________________
vector <Double_t> Binning::GetCostValues(){
  return fCostValues;
}
//______________________________________________________________________________
vector <Double_t> Binning::GetPhiValues(){
  return fPhiValues;
}
//______________________________________________________________________________
vector <Double_t> Binning::GetPhiTildeValues(){
  return fPhiTildeValues;
}
//______________________________________________________________________________
vector <Double_t> Binning::GetCostWidth(){
  for(int i = 0;i < fNCostBins;i++){
    fCostWidth.push_back(fCostValues[i+1] - fCostValues[i]);
  }
  return fCostWidth;
}
//______________________________________________________________________________
vector <Double_t> Binning::GetPhiWidth(){
  for(int i = 0;i < fNPhiBins;i++){
    fPhiWidth.push_back(fPhiValues[i+1] - fPhiValues[i]);
  }
  return fPhiWidth;
}
//______________________________________________________________________________
vector <Double_t> Binning::GetPhiTildeWidth(){
  for(int i = 0;i < fNPhiTildeBins;i++){
    fPhiTildeWidth.push_back(fPhiTildeValues[i+1] - fPhiTildeValues[i]);
  }
  return fPhiTildeWidth;
}
//______________________________________________________________________________
vector <Int_t> Binning::GetCostBinsMin(){
  return fCostBinsMin;
}
//______________________________________________________________________________
vector <Int_t> Binning::GetCostBinsMax(){
  return fCostBinsMax;
}
//______________________________________________________________________________
vector <Int_t> Binning::GetPhiBinsMin(){
  return fPhiBinsMin;
}
//______________________________________________________________________________
vector <Int_t> Binning::GetPhiBinsMax(){
  return fPhiBinsMax;
}
//______________________________________________________________________________
vector <Int_t> Binning::GetPhiTildeBinsMin(){
  return fPhiTildeBinsMin;
}
//______________________________________________________________________________
vector <Int_t> Binning::GetPhiTildeBinsMax(){
  return fPhiTildeBinsMax;
}
//______________________________________________________________________________
vector < vector <Double_t> > Binning::GetCellAreaMatrix(){
  vector <double> ColumnMatrix;
  for(int i = 0;i < fNCostBins;i++){
    for(int j = 0;j < fNPhiBins;j++){
      ColumnMatrix.push_back((TMath::Abs(fCostValues[i+1] - fCostValues[i]))*(fPhiValues[j+1] - fPhiValues[j]));
    }
    fCellAreaMatrix.push_back(ColumnMatrix);
    ColumnMatrix.clear();
  }
  return fCellAreaMatrix;
}
