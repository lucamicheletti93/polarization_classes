from ROOT import *
from array import array
import os.path
import math

gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("0.2g");
gROOT.ProcessLineSync(".x AccxEffCalculator.cxx+")
gROOT.ProcessLineSync(".x Binning.cxx+")

PI = math.pi

# crea in maniera automatica la cartella di output se questa non esiste

fileBinning = TFile.Open("output/binning.root")
binning = fileBinning.Get("Binning")
CostValues = binning.GetCostValues()
CostWidth = binning.GetCostWidth()
PhiValues = binning.GetPhiValues()
PhiWidth = binning.GetPhiWidth()

fileNJpsi = TFile.Open("/home/luca/GITHUB/polarization/1D_approach/signal_extraction/binned_1D_2pt4_test/2pt4.root")
histNJpsiCost = fileNJpsi.Get("histNJpsiCost")
histNJpsiPhi = fileNJpsi.Get("histNJpsiPhi")

canvasNJpsiCost = TCanvas("canvasNJpsiCost","canvasNJpsiCost",20,20,600,600)
histNJpsiCost.Draw("E")

canvasNJpsiPhi = TCanvas("canvasNJpsiPhi","canvasNJpsiPhi",20,20,600,600)
histNJpsiPhi.Draw("E")

fileAccxEff = TFile.Open("output/AccxEffFullStat.root")
histAccxEffCost = fileAccxEff.Get("histAccxEffCost_2pT4")
histAccxEffPhi = fileAccxEff.Get("histAccxEffPhi_2pT4")

histNJpsiCostCorr = histNJpsiCost.Clone("histNJpsiCostCorr")
for i in range(19):
    histNJpsiCostCorr.SetBinContent(i+1,(histNJpsiCost.GetBinContent(i+1))/CostWidth[i])
    histNJpsiCostCorr.SetBinError(i+1,(histNJpsiCost.GetBinError(i+1))/CostWidth[i])

histNJpsiPhiCorr = histNJpsiPhi.Clone("histNJpsiPhiCorr")
for i in range(10):
    histNJpsiPhiCorr.SetBinContent(i+1,(histNJpsiPhi.GetBinContent(i+1))/PhiWidth[i])
    histNJpsiPhiCorr.SetBinError(i+1,(histNJpsiPhi.GetBinError(i+1))/PhiWidth[i])

histNJpsiCostCorr.Divide(histAccxEffCost)
funcCost = TF1("funcCost","([0]/(3 + [1]))*(1 + [1]*x*x)",-0.8,0.8)
histNJpsiCostCorr.Fit(funcCost,"R0")

canvasNJpsiCostCorr = TCanvas("canvasNJpsiCostCorr","canvasNJpsiCostCorr",20,20,600,600)
histNJpsiCostCorr.Draw("E")
funcCost.Draw("same")

histNJpsiPhiCorr.Divide(histAccxEffPhi)
funcPhi = TF1("funcPhi","[0]*(1 + ((2*[2])/(3 + [1]))*cos(2*x))",0,PI)
funcPhi.FixParameter(1,funcCost.GetParameter(1))
histNJpsiPhiCorr.Fit(funcPhi,"R0")

canvasNJpsiPhiCorr = TCanvas("canvasNJpsiPhiCorr","canvasNJpsiPhiCorr",20,20,600,600)
histNJpsiPhiCorr.Draw("E")
funcPhi.Draw("same")

################################################################################
#if os.path.isfile("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root"):
    #fileDataMC = TFile.Open("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root")  # lxplus
#else:
    #fileDataMC = TFile.Open("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root") # local

#treeDataMC = fileDataMC.Get("MCTree")

#AccxEffReWeight1stStep = AccxEffCalculator(treeDataMC)
#AccxEffReWeight1stStep.SetPtBins(5,array('d',[0.,2.,4.,6.,10.]),array('d',[2.,4.,6.,10.,1000.]))
#AccxEffReWeight1stStep.SetBinning(CostValues,PhiValues)
#AccxEffReWeight1stStep.ReWeightAccxEff(funcCost.GetParameter(1),funcPhi.GetParameter(2),"FullStat","iterative_procedure/AccxEffReWeighted1stStep.root")

#fileAccxEffReWeight1stStep = TFile.Open("iterative_procedure/AccxEffReWeighted1stStep.root")
#histAccxEffCostReWeighted1stStep = fileAccxEffReWeight1stStep.Get("histAccxEffCostReWeighted_2pT4")
#histAccxEffCostReWeighted1stStep.SetLineColor(kRed);
#histAccxEffPhiReWeighted1stStep = fileAccxEffReWeight1stStep.Get("histAccxEffPhiReWeighted_2pT4")
#histAccxEffPhiReWeighted1stStep.SetLineColor(kRed);

#canvasAccxEffCost = TCanvas("canvasAccxEffCost","canvasAccxEffCost",20,20,600,600)
#histAccxEffCost.Draw("E")
#histAccxEffCostReWeighted1stStep.Draw("Esame")

#canvasAccxEffPhi = TCanvas("canvasAccxEffPhi","canvasAccxEffPhi",20,20,600,600)
#histAccxEffPhi.Draw("E")
#histAccxEffPhiReWeighted1stStep.Draw("Esame")

raw_input()
