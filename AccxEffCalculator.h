#ifndef ACCXEFFCALCULATOR_H
#define ACCXEFFCALCULATOR_H
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TObjArray.h"
#include <string>
#include <vector>

class AccxEffCalculator : public TObject
{

 public:
   AccxEffCalculator();
   AccxEffCalculator(TTree *treeAccxEff);
   AccxEffCalculator(Double_t scaleFactor_PbPb2015, TTree *treeAccxEff_PbPb2015, Double_t  scaleFactor_PbPb2018, TTree *treeAccxEff_PbPb2018);
   virtual ~AccxEffCalculator();

   void SetPtBins(Int_t, Double_t [],Double_t []);
   void SetBinning(vector <Double_t> , vector <Double_t>, vector <Double_t>);
   void ComputeResolution(string strSample, string nameOutputFile);
   void ComputeAccxEff(string strSample, string nameOutputFile);
   //void ReWeightAccxEff(string refFrame, Double_t LambdaTheta, Double_t LambdaPhi, string strSample, Bool_t saveFile, string nameOutputFile);
   //void ReWeightAccxEff(Double_t polParHE[], Double_t polParCS[], string strSample, Bool_t saveFile, string nameOutputFile);
   void ReWeightAccxEff(Double_t polParHE[4][4], Double_t polParCS[4][4], string strSample, Bool_t saveFile, string nameOutputFile);
   void ReWeightAccxEff_PbPb2015_PbPb2018(Double_t polParHE[4][4], Double_t polParCS[4][4], string strSample, Bool_t saveFile, string nameOutputFile);
   void ReWeightAccxEff_PbPb2015_PbPb2018_Upsilon(Double_t polParHE[3][1], Double_t polParCS[3][1], string strSample, Bool_t saveFile, string nameOutputFile);
   void ComputeTriggerResponseFunction(string strSample, string nameOutputFile);
   //void ComputeReweightTRFAccxEff(string , string , bool , TH1D *);
   void ComputeReweightTRFAccxEff(string , string , bool , TObjArray *);

   Int_t GetNPtBins(){return fNPtBins;};
   Int_t GetNCosThetaBins(){return fNCosThetaBins;};
   Int_t GetNPhiBins(){return fNPhiBins;};

 private:
   Double_t fPi;
   Int_t fNPtBins;
   vector <Double_t> fMinPt;
   vector <Double_t> fMaxPt;

   Int_t fNCosThetaBins;
   vector <Double_t> fCosThetaValues;
   Int_t fNPhiBins;
   vector <Double_t> fPhiValues;
   Int_t fNPhiTildeBins;
   vector <Double_t> fPhiTildeValues;

   Double_t fPhiTilde;

   TTree *fTreeAccxEff;
   
   Double_t fScaleFactor_PbPb2015_PbPb2018[2];
   TTree *fTreeAccxEff_PbPb2015_PbPb2018[2];

   // Resolution
   /*TH1D *fHistResolutionCosThetaHE;
   TH1D *fHistResolutionPhiHE;
   TH2D *fHistResolutionCosThetaPhiHE;

   TH1D *fHistResolutionCosThetaCS;
   TH1D *fHistResolutionPhiCS;
   TH2D *fHistResolutionCosThetaPhiCS;*/

   TH1D *fHistResolutionCosThetaHE[13];
   TH1D *fHistResolutionPhiHE[13];
   TH2D *fHistResolutionCosThetaPhiHE[13];

   TH1D *fHistResolutionCosThetaCS[13];
   TH1D *fHistResolutionPhiCS[13];
   TH2D *fHistResolutionCosThetaPhiCS[13];

   // cos(theta)
   TH1D *fHistGenCosThetaHE[13];
   TH1D *fHistRecCosThetaHE[13];
   TH1D *fHistAccxEffCosThetaHE[13];
   TH1D *fHistGenCosThetaCS[13];
   TH1D *fHistRecCosThetaCS[13];
   TH1D *fHistAccxEffCosThetaCS[13];
   // phi
   TH1D *fHistGenPhiHE[13];
   TH1D *fHistRecPhiHE[13];
   TH1D *fHistAccxEffPhiHE[13];
   TH1D *fHistGenPhiCS[13];
   TH1D *fHistRecPhiCS[13];
   TH1D *fHistAccxEffPhiCS[13];
   // phi Tilde
   TH1D *fHistGenPhiTildeHE[13];
   TH1D *fHistRecPhiTildeHE[13];
   TH1D *fHistAccxEffPhiTildeHE[13];
   TH1D *fHistGenPhiTildeCS[13];
   TH1D *fHistRecPhiTildeCS[13];
   TH1D *fHistAccxEffPhiTildeCS[13];

   // pT
   TH1D *fHistGenPtHE[13];
   TH1D *fHistRecPtHE[13];
   TH1D *fHistAccxEffPtHE[13];
   TH1D *fHistGenPtCS[13];
   TH1D *fHistRecPtCS[13];
   TH1D *fHistAccxEffPtCS[13];

   /////////////////////////////////////////////////////////////////////////////
   TH1D *fHistGenPhiTildeNarrowHE[13];
   TH1D *fHistRecPhiTildeNarrowHE[13];
   TH1D *fHistAccxEffPhiTildeNarrowHE[13];
   TH1D *fHistRecPhiTildeNarrowLeftBellHE[13];
   TH1D *fHistRecPhiTildeNarrowRightBellHE[13];
   TH1D *fHistRecPhiTildeNarrowAllBellHE[13];

   TH1D *fHistGenPhiTildeNarrowCS[13];
   TH1D *fHistRecPhiTildeNarrowCS[13];
   TH1D *fHistAccxEffPhiTildeNarrowCS[13];
   TH1D *fHistRecPhiTildeNarrowLeftBellCS[13];
   TH1D *fHistRecPhiTildeNarrowRightBellCS[13];
   TH1D *fHistRecPhiTildeNarrowAllBellCS[13];
   /////////////////////////////////////////////////////////////////////////////

   // cos(theta),phi
   TH2D *fHistGenCosThetaPhiNarrowHE[13];
   TH2D *fHistRecCosThetaPhiNarrowHE[13];
   TH2D *fHistAccxEffCosThetaPhiNarrowHE[13];

   TH2D *fHistGenCosThetaPhiHE[13];
   TH2D *fHistRecCosThetaPhiHE[13];
   TH2D *fHistAccxEffCosThetaPhiHE[13];
   TH2D *fHistAccxEffCosThetaPhiStatRelHE[13];

   TH2D *fHistGenCosThetaPhiNarrowCS[13];
   TH2D *fHistRecCosThetaPhiNarrowCS[13];
   TH2D *fHistAccxEffCosThetaPhiNarrowCS[13];

   TH2D *fHistGenCosThetaPhiCS[13];
   TH2D *fHistRecCosThetaPhiCS[13];
   TH2D *fHistAccxEffCosThetaPhiCS[13];
   TH2D *fHistAccxEffCosThetaPhiStatRelCS[13];

   // cos(theta),pT
   TH2D *fHistGenCosThetaHEPt;
   TH2D *fHistRecCosThetaHEPt;
   TH2D *fHistAccxEffCosThetaHEPt;
   TH2D *fHistGenCosThetaCSPt;
   TH2D *fHistRecCosThetaCSPt;
   TH2D *fHistAccxEffCosThetaCSPt;
   // phi,pT
   TH2D *fHistGenPhiHEPt;
   TH2D *fHistRecPhiHEPt;
   TH2D *fHistAccxEffPhiHEPt;
   TH2D *fHistGenPhiCSPt;
   TH2D *fHistRecPhiCSPt;
   TH2D *fHistAccxEffPhiCSPt;
   // phiTilde,pT
   TH2D *fHistGenPhiTildeHEPt;
   TH2D *fHistRecPhiTildeHEPt;
   TH2D *fHistAccxEffPhiTildeHEPt;
   TH2D *fHistGenPhiTildeCSPt;
   TH2D *fHistRecPhiTildeCSPt;
   TH2D *fHistAccxEffPhiTildeCSPt;

   // cosTheta,phi,eta (winn test)
   TH2D *fHistRecCosThetaHEEtaSM;
   TH2D *fHistRecCosThetaCSEtaSM;
   TH2D *fHistRecPhiHEEtaSM;
   TH2D *fHistRecPhiCSEtaSM;
   TH2D *fHistRecCosThetaPhiHESM;
   TH2D *fHistRecCosThetaPhiCSSM;

   // cos(theta)
   TH1D *fHistGenCosThetaHEReWeighted[13];
   TH1D *fHistRecCosThetaHEReWeighted[13];
   TH1D *fHistAccxEffCosThetaHEReWeighted[13];
   TH1D *fHistGenCosThetaCSReWeighted[13];
   TH1D *fHistRecCosThetaCSReWeighted[13];
   TH1D *fHistAccxEffCosThetaCSReWeighted[13];
   // phi
   TH1D *fHistGenPhiHEReWeighted[13];
   TH1D *fHistRecPhiHEReWeighted[13];
   TH1D *fHistAccxEffPhiHEReWeighted[13];
   TH1D *fHistGenPhiCSReWeighted[13];
   TH1D *fHistRecPhiCSReWeighted[13];
   TH1D *fHistAccxEffPhiCSReWeighted[13];
   // phiTilde
   TH1D *fHistGenPhiTildeHEReWeighted[13];
   TH1D *fHistRecPhiTildeHEReWeighted[13];
   TH1D *fHistAccxEffPhiTildeHEReWeighted[13];
   TH1D *fHistGenPhiTildeCSReWeighted[13];
   TH1D *fHistRecPhiTildeCSReWeighted[13];
   TH1D *fHistAccxEffPhiTildeCSReWeighted[13];

   // (cos(theta),phi)
   TH2D *fHistGenCosThetaPhiHEReWeighted[13];
   TH2D *fHistRecCosThetaPhiHEReWeighted[13];
   TH2D *fHistAccxEffCosThetaPhiHEReWeighted[13];
   TH2D *fHistGenCosThetaPhiCSReWeighted[13];
   TH2D *fHistRecCosThetaPhiCSReWeighted[13];
   TH2D *fHistAccxEffCosThetaPhiCSReWeighted[13];

   // Trigger response fuction calculation
   TH1D *fHistLowPt;
   TH1D *fHistAllPt;
   TH1D *fHistTriggerResponseFunction;

   // tree variables
   Int_t fNDimuGen;
   Float_t  fPercV0M;
   Double_t fDimuPtGen[3000];
   Double_t fDimuYGen[3000];
   Double_t fCosThetaHEGen[3000];
   Double_t fPhiHEGen[3000];
   Double_t fCosThetaCSGen[3000];
   Double_t fPhiCSGen[3000];

   Int_t fNDimuRec;
   Double_t fDimuPtRec[3000];
   Double_t fDimuYRec[3000];
   Double_t fDimuMassRec[3000];
   Int_t fDimuMatchRec[3000];
   Double_t fCosThetaHERec[3000];
   Double_t fPhiHERec[3000];
   Double_t fCosThetaCSRec[3000];
   Double_t fPhiCSRec[3000];

   Int_t fNMuonsRec;
   Double_t fPtRec[3000];
   Double_t fEtaRec[3000];
   Int_t fChargeRec[3000];
   Int_t fMatchTrigRec[3000];

   Int_t    fDimuPDG_gen[3000];
   Int_t    fDimuPDG_rec[3000];

   Int_t    fDimuMu_rec[3000][2];
   Double_t fPt_rec[3000];
   Double_t fPt_gen[3000];

ClassDef(AccxEffCalculator,1)
};

#endif
