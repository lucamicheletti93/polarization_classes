from ROOT import *
from array import array
import os.path
import numpy as np

gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("0.2g");
gROOT.ProcessLineSync(".x ../SpecialFitCalculator.cxx+")

namePtRanges = ['2pT4','4pT6','6pT10']
nameIdIterations = ['1st','2nd','3rd']
minNamePt = ['2','4','6']
maxNamePt = ['4','6','10']

lambdaTheta = np.zeros((len(namePtRanges), len(nameIdIterations)))
lambdaPhi = np.zeros((len(namePtRanges), len(nameIdIterations)))

for iterStep in range(len(nameIdIterations)):
    for i in range(len(namePtRanges)):
        fileNJpsi = TFile.Open("/home/luca/GITHUB/polarization/1D_approach/signal_extraction/binned_1D_" + minNamePt[i] + "pt" + maxNamePt[i] + "_test/" + minNamePt[i] + "pt" + maxNamePt[i] + ".root")
        histNJpsiCost = fileNJpsi.Get("histNJpsiCost")
        histNJpsiPhi = fileNJpsi.Get("histNJpsiPhi")

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
            print "if you want to reproduce it delete " + "/AccxEffReWeighted" + nameIdIterations[i] + "Step.root" + " and re-run"
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
            AccxEffReWeight.ReWeightAccxEff(lambdaTheta[i, iterStep],lambdaPhi[i, iterStep],"FullStat",kTRUE,"iterative_procedure/" + namePtRanges[i] + "/AccxEffReWeighted" + nameIdIterations[i] + "Step.root")
            del AccxEffReWeight
            fileDataMC.Close

raw_input()