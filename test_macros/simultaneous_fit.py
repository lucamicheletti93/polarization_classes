from ROOT import *
import math

gROOT.ProcessLineSync(".x ../Binning.cxx+")

PI = math.pi

fileBinning = TFile.Open("output/binning.root")
binning = fileBinning.Get("Binning")
CostValues = binning.GetCostValues()
CostWidth = binning.GetCostWidth()
PhiValues = binning.GetPhiValues()
PhiWidth = binning.GetPhiWidth()

# Create the observables
cosTheta = RooRealVar("cosTheta","cos#it{#theta}_{HE}",-0.8,0.8)
phi = RooRealVar("phi","#it{#varphi}_{HE}",0,PI)

fileNJpsi = TFile.Open("/home/luca/GITHUB/polarization/1D_approach/signal_extraction/binned_1D_2pt4_test/2pt4.root")
histNJpsiCost = fileNJpsi.Get("histNJpsiCost")

fileAccxEff = TFile.Open("output/AccxEffFullStat.root")
histAccxEffCost = fileAccxEff.Get("histAccxEffCost_2pT4")


histNJpsiCostCorr = histNJpsiCost.Clone("histNJpsiCostCorr")
for i in range(19):
    histNJpsiCostCorr.SetBinContent(i+1,(histNJpsiCost.GetBinContent(i+1))/CostWidth[i])
    histNJpsiCostCorr.SetBinError(i+1,(histNJpsiCost.GetBinError(i+1))/CostWidth[i])
histNJpsiCostCorr.Divide(histAccxEffCost)
histNJpsiCostCorr.Sumw2(kTRUE)

rooHistNJpsiCostCorr = RooDataHist("rooHistNJpsiCostCorr","rooHistNJpsiCostCorr",RooArgList(cosTheta),histNJpsiCostCorr)

lambdaTheta = RooRealVar("lambdaTheta","#lambda_{#theta}",-1.,1.)

# Create cosTheta PDF
funcCosTheta = RooGenericPdf("funcCosTheta","funcCosTheta","(1/(3. + lambdaTheta))*(1 + lambdaTheta*cosTheta*cosTheta)",RooArgList(cosTheta,lambdaTheta))
cosTheta.setRange("R1",-0.6,0.6)
#funcCosTheta.chi2FitTo(rooHistNJpsiCostCorr,RooFit.Range("R1"))
#funcCosTheta.chi2FitTo(rooHistNJpsiCostCorr)
#funcCosTheta.fitTo(rooHistNJpsiCostCorr,RooFit.SumW2Error(kFALSE),RooFit.Range("R1"))


#chi2Fit = RooChi2Var("chi2Fit","chi2Fit",funcCosTheta,rooHistNJpsiCostCorr,RooFit.DataError(RooAbsData.SumW2Error(kTRUE)),RooFit.Range("R1"))
#chi2Fit = RooChi2Var("chi2Fit","chi2Fit",funcCosTheta,rooHistNJpsiCostCorr,RooFit.SumW2Error(kFALSE),RooFit.Range("R1"))
chi2Fit = RooChi2Var("chi2Fit","chi2Fit",funcCosTheta,rooHistNJpsiCostCorr,RooFit.SumW2Error(kFALSE))
minuit = RooMinuit(chi2Fit)
minuit.migrad()
minuit.hesse()

  #// Plot chi^2 fit result on frame as well
  #RooFitResult* r_chi2_wgt = minuit.save() ;
  #p2.plotOn(frame,LineStyle(kDashed),LineColor(kRed)) ;


# Create phi PDF
#lambdaPhi = RooRealVar("lambdaPhi","#lambda_{#varphi}",0.5,-1.,1.)
#funcPhi = RooGenericPdf("funcPhi","funcPhi","1 + ((2.*lambdaPhi)/(3. + lambdaTheta))*cos(2.*phi)",RooArgList(phi,lambdaTheta,lambdaPhi))

# Create a composite model
#norm = RooRealVar("norm","N",0.2,0.,1.)
#model = RooAddPdf("model","model",RooArgList(funcCosTheta,funcPhi),norm)


cosThetaframe = cosTheta.frame()
#funcCosTheta.plotOn(cosThetaframe)
#rooHistNJpsiCost.plotOn(cosThetaframe)
rooHistNJpsiCostCorr.plotOn(cosThetaframe)
funcCosTheta.plotOn(cosThetaframe)
cosThetaframe.Draw()

raw_input()
