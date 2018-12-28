#ifndef BINNING_H
#define BINNING_H
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include <string>
#include <vector>

class Binning : public TObject
{

 public:
   Binning();
   virtual ~Binning();

   void ConfigureBinValues(int NCostBins, int CostBinsMin[], int CostBinsMax[], int NPhiBins, int PhiBinsMin[], int PhiBinsMax[], int NPhiTildeBins, int PhiTildeBinsMin[], int PhiTildeBinsMax[]);
   void PrintBinValues();

   vector <Double_t> GetCostValues();
   vector <Double_t> GetPhiValues();
   vector <Double_t> GetPhiTildeValues();
   vector <Double_t> GetCostWidth();
   vector <Double_t> GetPhiWidth();
   vector <Double_t> GetPhiTildeWidth();
   vector <Int_t> GetCostBinsMin();
   vector <Int_t> GetCostBinsMax();
   vector <Int_t> GetPhiBinsMin();
   vector <Int_t> GetPhiBinsMax();
   vector <Int_t> GetPhiTildeBinsMin();
   vector <Int_t> GetPhiTildeBinsMax();
   vector < vector <Double_t> > GetCellAreaMatrix();



 private:
   Double_t fCosThetaMin;
   Double_t fCosThetaMax;
   Double_t fPhiMin;
   Double_t fPhiMax;
   Double_t fPhiTildeMin;
   Double_t fPhiTildeMax;
   Int_t fNCostBins;
   Int_t fNPhiBins;
   Int_t fNPhiTildeBins;
   vector <Double_t> fCostValues;
   vector <Double_t> fPhiValues;
   vector <Double_t> fPhiTildeValues;
   vector <Double_t> fCostWidth;
   vector <Double_t> fPhiWidth;
   vector <Double_t> fPhiTildeWidth;
   vector <Int_t> fCostBinsMin;
   vector <Int_t> fCostBinsMax;
   vector <Int_t> fPhiBinsMin;
   vector <Int_t> fPhiBinsMax;
   vector <Int_t> fPhiTildeBinsMin;
   vector <Int_t> fPhiTildeBinsMax;
   vector < vector <Double_t> > fCellAreaMatrix;

ClassDef(Binning,1)
};

#endif
