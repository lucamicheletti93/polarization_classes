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

   void ConfigureBinValues(int NCostBins, int CostBinsMin[], int CostBinsMax[], int NPhiBins, int PhiBinsMin[], int PhiBinsMax[]);
   void PrintBinValues();

   vector <Double_t> GetCostValues();
   vector <Double_t> GetPhiValues();
   vector <Double_t> GetCostWidth();
   vector <Double_t> GetPhiWidth();
   vector <Int_t> GetCostBinsMin();
   vector <Int_t> GetCostBinsMax();
   vector <Int_t> GetPhiBinsMin();
   vector <Int_t> GetPhiBinsMax();
   vector < vector <Double_t> > GetCellAreaMatrix();



 private:
   Int_t fNCostBins;
   Int_t fNPhiBins;
   vector <Double_t> fCostValues;
   vector <Double_t> fPhiValues;
   vector <Double_t> fCostWidth;
   vector <Double_t> fPhiWidth;
   vector <Int_t> fCostBinsMin;
   vector <Int_t> fCostBinsMax;
   vector <Int_t> fPhiBinsMin;
   vector <Int_t> fPhiBinsMax;
   vector < vector <Double_t> > fCellAreaMatrix;

ClassDef(Binning,1)
};

#endif
