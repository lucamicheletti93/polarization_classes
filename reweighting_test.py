from ROOT import *
from array import array

gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("0.2g");
gROOT.ProcessLineSync(".x AccxEffCalculator.cxx+")
gROOT.ProcessLineSync(".x Binning.cxx+")

fileBinning = TFile.Open("/home/luca/GITHUB/polarization_classes/output/binning.root")
binning = fileBinning.Get("Binning")
CostValues = binning.GetCostValues()
PhiValues = binning.GetPhiValues()

#fileDataMC = TFile.Open("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root")  # lxplus
fileDataMC = TFile.Open("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root") # local
treeDataMC = fileDataMC.Get("MCTree")

minPt = [0,2,4,6,10]
maxPt = [2,4,6,10,1000]

# lambda theta = 0.8
AccxEffReWeight1 = AccxEffCalculator(treeDataMC)
AccxEffReWeight1.SetPtBins(5,array('d',[0.,2.,4.,6.,10.]),array('d',[2.,4.,6.,10.,1000.]))
AccxEffReWeight1.SetBinning(CostValues,PhiValues)
AccxEffReWeight1.ReWeightAccxEff(0.8,"TestStat","output/AccxEffReWeightedPol1.root")

fileAccxEffReWeightedPol1 = TFile.Open("output/AccxEffReWeightedPol1.root")

histGenCostReWeightedPol1 = []
histAccxEffCostReWeightedPol1 = []

for i in range(5):
    histGenCostReWeightedPol1.append(fileAccxEffReWeightedPol1.Get('histGenCostReWeighted' + str(0.8) + '_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histGenCostReWeightedPol1[i].SetDirectory(0)
    histAccxEffCostReWeightedPol1.append(fileAccxEffReWeightedPol1.Get('histAccxEffCostReWeighted' + str(0.8) + '_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffCostReWeightedPol1[i].SetDirectory(0)
fileDataMC.Close()
fileAccxEffReWeightedPol1.Close()


fileDataMC = TFile.Open("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root") # local
treeDataMC = fileDataMC.Get("MCTree")

# lambda theta = 0.0
AccxEffReWeight2 = AccxEffCalculator(treeDataMC)
AccxEffReWeight2.SetPtBins(5,array('d',[0.,2.,4.,6.,10.]),array('d',[2.,4.,6.,10.,1000.]))
AccxEffReWeight2.SetBinning(CostValues,PhiValues)
AccxEffReWeight2.ReWeightAccxEff(0.0,"TestStat","output/AccxEffReWeightedPol2.root")

fileAccxEffReWeightedPol2 = TFile.Open("output/AccxEffReWeightedPol2.root")

histGenCostReWeightedPol2 = []
histAccxEffCostReWeightedPol2 = []

for i in range(5):
    histGenCostReWeightedPol2.append(fileAccxEffReWeightedPol2.Get('histGenCostReWeighted' + str(0.0) + '_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histGenCostReWeightedPol2[i].SetDirectory(0)
    histGenCostReWeightedPol2[i].SetLineColor(kRed)
    histAccxEffCostReWeightedPol2.append(fileAccxEffReWeightedPol2.Get('histAccxEffCostReWeighted' + str(0.0) + '_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffCostReWeightedPol2[i].SetDirectory(0)
    histAccxEffCostReWeightedPol2[i].SetLineColor(kRed)
fileDataMC.Close()
fileAccxEffReWeightedPol2.Close()

canvasAccxEffCostReWeightedComp = TCanvas("canvasAccxEffCostReWeightedComp","canvasAccxEffCostReWeightedComp",20,20,600,600)
canvasAccxEffCostReWeightedComp.cd()
#histGenCostReWeightedPol2[2].Draw()
#histGenCostReWeightedPol1[2].Draw("same")
histAccxEffCostReWeightedPol2[3].Draw()
histAccxEffCostReWeightedPol1[3].Draw("same")


raw_input()
