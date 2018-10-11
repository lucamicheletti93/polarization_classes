from ROOT import *
from array import array
import os.path
import math
import numpy as np

gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("0.2g");
gROOT.ProcessLineSync(".x ../Binning.cxx+")
gROOT.ProcessLineSync(".x ../AccxEffCalculator.cxx+")
gROOT.ProcessLineSync(".x ../SpecialFitCalculator.cxx+")

# Setting main quantities
PI = math.pi
fileBinning = TFile.Open("output/binning.root")
binning = fileBinning.Get("Binning")
CostValues = binning.GetCostValues()
CostWidth = binning.GetCostWidth()
PhiValues = binning.GetPhiValues()
PhiWidth = binning.GetPhiWidth()

namePtRanges = ['2pT4','4pT6','6pT10']
nameIdIterations = ['1st','2nd','3rd']
minNamePt = ['2','4','6']
maxNamePt = ['4','6','10']
minPtBin = [2.,4.,6.]
maxPtBin = [4.,6.,10.]

for i in range(len(namePtRanges)):
    if not os.path.exists('iterative_procedure/' + namePtRanges[i]):
        os.makedirs('iterative_procedure/' + namePtRanges[i])

lambdaTheta = np.zeros((len(namePtRanges), len(nameIdIterations)))
lambdaPhi = np.zeros((len(namePtRanges), len(nameIdIterations)))

for iterStep in range(len(nameIdIterations)):
    for i in range(len(namePtRanges)):
        if os.path.isfile("/afs/cern.ch/user/l/lmichele/private/Jpsi_polarization/signal_extraction_1D/" + minNamePt[i] + "pt" + maxNamePt[i] + ".root"):
            fileNJpsi = TFile.Open("/afs/cern.ch/user/l/lmichele/private/Jpsi_polarization/signal_extraction_1D/" + minNamePt[i] + "pt" + maxNamePt[i] + ".root")  # lxplus
        else:
            fileNJpsi = TFile.Open("/home/luca/GITHUB/polarization/1D_approach/signal_extraction/binned_1D_" + minNamePt[i] + "pt" + maxNamePt[i] + "_test/" + minNamePt[i] + "pt" + maxNamePt[i] + ".root") # local
        histNJpsiCost = fileNJpsi.Get("histNJpsiCost")
        histNJpsiPhi = fileNJpsi.Get("histNJpsiPhi")

        if iterStep == 0:
            print "output/AccxEffFullStat.root"
            fileAccxEff = TFile.Open("output/AccxEffFullStat.root","READ")
            histAccxEffCost = fileAccxEff.Get("histAccxEffCost_" + namePtRanges[i])
            histAccxEffPhi = fileAccxEff.Get("histAccxEffPhi_" + namePtRanges[i])
        else:
            print "iterative_procedure/" + namePtRanges[i] + "/AccxEffReWeighted" + nameIdIterations[iterStep] + "Step.root"
            fileAccxEff = TFile.Open("iterative_procedure/" + namePtRanges[i] + "/AccxEffReWeighted" + nameIdIterations[iterStep] + "Step.root")
            histAccxEffCost = fileAccxEff.Get("histAccxEffCostReWeighted_" + namePtRanges[i])
            histAccxEffPhi = fileAccxEff.Get("histAccxEffPhiReWeighted_" + namePtRanges[i])

        SimFit = SpecialFitCalculator()
        SimFit.SimultaneousFit(histNJpsiCost,histNJpsiPhi,histAccxEffCost,histAccxEffPhi)
        lambdaTheta[i, iterStep] = SimFit.GetLambdaTheta()
        lambdaPhi[i, iterStep] = SimFit.GetLambdaPhi()
        fileNJpsi.Close()
        fileAccxEff.Close()

    for i in range(len(namePtRanges)):
        if os.path.isfile("iterative_procedure/" + namePtRanges[i] + "/AccxEffReWeighted" + nameIdIterations[iterStep] + "Step.root"):
            print "--> " + "AccxEffReWeighted" + nameIdIterations[iterStep] + "Step.root" + " has already produced"
            print "if you want to reproduce it delete " + "/AccxEffReWeighted" + nameIdIterations[iterStep] + "Step.root" + " and re-run"
        else:
            print "--> iterative_procedure/" + namePtRanges[i] + "/AccxEffReWeighted" + nameIdIterations[iterStep] + "Step.root" + " does not exist"
            if os.path.isfile("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root"):
                fileDataMC = TFile.Open("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root")  # lxplus
            else:
                fileDataMC = TFile.Open("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root") # local
            treeDataMC = fileDataMC.Get("MCTree")
            AccxEffReWeight = AccxEffCalculator(treeDataMC)
            AccxEffReWeight.SetPtBins(1,array('d',[minPtBin[i]]),array('d',[maxPtBin[i]]))
            AccxEffReWeight.SetBinning(CostValues,PhiValues)
            AccxEffReWeight.ReWeightAccxEff(lambdaTheta[i, iterStep],lambdaPhi[i, iterStep],"FullStat",kTRUE,"iterative_procedure/" + namePtRanges[i] + "/AccxEffReWeighted" + nameIdIterations[iterStep] + "Step.root")
            del AccxEffReWeight
            fileDataMC.Close

raw_input()
