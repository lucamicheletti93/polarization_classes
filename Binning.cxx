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
}
//______________________________________________________________________________
Binning::~Binning() {
  // destructor
}
//______________________________________________________________________________
void Binning::ConfigureBinValues(int NCostBins, int CostBinsMin[], int CostBinsMax[], int NPhiBins, int PhiBinsMin[], int PhiBinsMax[]) {

  double PI = TMath::Pi();
  TH2D *hist = new TH2D("hist","hist",100,-1.,1.,50,0.,PI);
  double bin;

    fNCostBins = NCostBins;
    for(int i = 0;i < NCostBins;i++){
      bin = hist -> GetXaxis() -> GetBinLowEdge(CostBinsMin[i]);
      fCostValues.push_back(bin);
      fCostBinsMin.push_back(CostBinsMin[i]);
      fCostBinsMax.push_back(CostBinsMax[i]);
    }
    bin = hist -> GetXaxis() -> GetBinLowEdge(101);
    fCostValues.push_back(bin);

    fNPhiBins = NPhiBins;
    for(int i = 0;i < NPhiBins;i++){
      bin = hist -> GetYaxis() -> GetBinLowEdge(PhiBinsMin[i]);
      fPhiValues.push_back(bin);
      fPhiBinsMin.push_back(PhiBinsMin[i]);
      fPhiBinsMax.push_back(PhiBinsMax[i]);
    }
    bin = hist -> GetYaxis() -> GetBinLowEdge(51);
    fPhiValues.push_back(bin);
    delete hist;
}
//______________________________________________________________________________
void Binning::PrintBinValues() {
  double PI = TMath::Pi();
  printf("CostValues[%i] = {",fNCostBins+1);
  for(int i = 0;i < fNCostBins-1;i++) printf("%3.2f,",fCostValues[i]);
  printf("%3.2f,%3.2f}\n",fCostValues[fNCostBins-1],1.);

  printf("PhiValues[%i] = {",fNPhiBins+1);
  for(int i = 0;i < fNPhiBins-1;i++) printf("%f,",fPhiValues[i]);
  printf("%f,%f}\n",fPhiValues[fNPhiBins-1],PI);
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
vector < vector <Double_t> > Binning::GetCellAreaMatrix(){
  vector <double> ColumnMatrix;
  for(int i = 0;i < fNCostBins;i++){
    for(int j = 0;j < fNPhiBins;j++){
      ColumnMatrix.push_back((fCostValues[i+1] - fCostValues[i])*(fPhiValues[j+1] - fPhiValues[j]));
    }
    fCellAreaMatrix.push_back(ColumnMatrix);
    ColumnMatrix.clear();
  }
  return fCellAreaMatrix;
}
