from ROOT import *
from array import array
import os.path
import math

gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("0.2g");
gROOT.ProcessLineSync(".x ../AccxEffCalculator.cxx+")
gROOT.ProcessLineSync(".x ../Binning.cxx+")

# Setting main quantities
PI = math.pi
fileBinning = TFile.Open("output/binning.root")
binning = fileBinning.Get("Binning")
CostValues = binning.GetCostValues()
CostWidth = binning.GetCostWidth()
PhiValues = binning.GetPhiValues()
PhiWidth = binning.GetPhiWidth()

fileAccxEff = TFile.Open("output/AccxEffFullStat.root")
histGenCost = fileAccxEff.Get("histGenCost_2pT4")
histAccxEffCost = fileAccxEff.Get("histAccxEffCost_2pT4")
histAccxEffPhi = fileAccxEff.Get("histAccxEffPhi_2pT4")

# Fitting procedure
#fileNJpsi = TFile.Open("/home/luca/GITHUB/polarization/1D_approach/signal_extraction/binned_1D_2pt4_test/2pt4.root")
#histNJpsiCost = fileNJpsi.Get("histNJpsiCost")
#histNJpsiPhi = fileNJpsi.Get("histNJpsiPhi")

#canvasNJpsiCost = TCanvas("canvasNJpsiCost","canvasNJpsiCost",20,20,600,600)
#histNJpsiCost.Draw("E")

#canvasNJpsiPhi = TCanvas("canvasNJpsiPhi","canvasNJpsiPhi",20,20,600,600)
#histNJpsiPhi.Draw("E")

#histNJpsiCostCorr = histNJpsiCost.Clone("histNJpsiCostCorr")
#for i in range(19):
    #histNJpsiCostCorr.SetBinContent(i+1,(histNJpsiCost.GetBinContent(i+1))/CostWidth[i])
    #histNJpsiCostCorr.SetBinError(i+1,(histNJpsiCost.GetBinError(i+1))/CostWidth[i])

#histNJpsiPhiCorr = histNJpsiPhi.Clone("histNJpsiPhiCorr")
#for i in range(10):
    #histNJpsiPhiCorr.SetBinContent(i+1,(histNJpsiPhi.GetBinContent(i+1))/PhiWidth[i])
    #histNJpsiPhiCorr.SetBinError(i+1,(histNJpsiPhi.GetBinError(i+1))/PhiWidth[i])

#histNJpsiCostCorr.Divide(histAccxEffCost)
#funcCost = TF1("funcCost","([0]/(3 + [1]))*(1 + [1]*x*x)",-0.8,0.8)
#histNJpsiCostCorr.Fit(funcCost,"R0")

#canvasNJpsiCostCorr = TCanvas("canvasNJpsiCostCorr","canvasNJpsiCostCorr",20,20,600,600)
#histNJpsiCostCorr.Draw("E")
#funcCost.Draw("same")

#histNJpsiPhiCorr.Divide(histAccxEffPhi)
#funcPhi = TF1("funcPhi","[0]*(1 + ((2*[2])/(3 + [1]))*cos(2*x))",0,PI)
#funcPhi.FixParameter(1,funcCost.GetParameter(1))
#histNJpsiPhiCorr.Fit(funcPhi,"R0")

#canvasNJpsiPhiCorr = TCanvas("canvasNJpsiPhiCorr","canvasNJpsiPhiCorr",20,20,600,600)
#histNJpsiPhiCorr.Draw("E")
#funcPhi.Draw("same")

################################################################################
if os.path.isfile("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root"):
    fileDataMC = TFile.Open("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root")  # lxplus
else:
    fileDataMC = TFile.Open("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root") # local

treeDataMC = fileDataMC.Get("MCTree")

# To check if a directory exists
if not os.path.exists('iterative_procedure'):
    os.makedirs('iterative_procedure')

nameOutputFile = 'AccxEffReWeighted3rdStep.root'

#if os.path.isfile("iterative_procedure/AccxEffReWeighted1stStep.root"):
if os.path.isfile("iterative_procedure/" + nameOutputFile):
    print nameOutputFile + "has already produced"
    print "if you want to reproduce it delete" + nameOutputFile + "and re-run"
else:
    AccxEffReWeight1stStep = AccxEffCalculator(treeDataMC)
    AccxEffReWeight1stStep.SetPtBins(1,array('d',[2.]),array('d',[4.]))
    AccxEffReWeight1stStep.SetBinning(CostValues,PhiValues)
    #AccxEffReWeight1stStep.ReWeightAccxEff(-0.189948,-0.2228,"FullStat",kTRUE,"iterative_procedure/" + nameOutputFile) # 1st iteration
    #AccxEffReWeight1stStep.ReWeightAccxEff(-0.20244,-0.262088,"FullStat",kTRUE,"iterative_procedure/" + nameOutputFile)  # 2nd iteration
    AccxEffReWeight1stStep.ReWeightAccxEff(-0.203284,-0.269251,"FullStat",kTRUE,"iterative_procedure/" + nameOutputFile)  # 3rd iteration

fileAccxEffReWeight1stStep = TFile.Open("iterative_procedure/" + nameOutputFile)
histAccxEffCostReWeighted1stStep = fileAccxEffReWeight1stStep.Get("histAccxEffCostReWeighted_2pT4")
histAccxEffCostReWeighted1stStep.SetLineColor(kRed);
histAccxEffPhiReWeighted1stStep = fileAccxEffReWeight1stStep.Get("histAccxEffPhiReWeighted_2pT4")
histAccxEffPhiReWeighted1stStep.SetLineColor(kRed);

histGenCostReWeighted1stStep = fileAccxEffReWeight1stStep.Get("histGenCostReWeighted_2pT4")
histGenCostReWeighted1stStep.SetLineColor(kRed);

canvasAccxEffCost = TCanvas("canvasAccxEffCost","canvasAccxEffCost",20,20,600,600)
histAccxEffCostReWeighted1stStep.Draw("E")
histAccxEffCost.Draw("Esame")


canvasAccxEffPhi = TCanvas("canvasAccxEffPhi","canvasAccxEffPhi",20,20,600,600)
histAccxEffPhi.Draw("E")
histAccxEffPhiReWeighted1stStep.Draw("Esame")


# Check plots
for i in range(19):
    histGenCostReWeighted1stStep.SetBinContent(i+1,histGenCostReWeighted1stStep.GetBinContent(i+1)/CostWidth[i])
    histGenCostReWeighted1stStep.SetBinError(i+1,histGenCostReWeighted1stStep.GetBinError(i+1)/CostWidth[i])

canvasGenCost = TCanvas("canvasGenCost","canvasGenCost",20,20,600,600)
histGenCostReWeighted1stStep.Draw()
histGenCost.Draw("Esame")

raw_input()
