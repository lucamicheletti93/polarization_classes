from ROOT import *
from array import array
import os.path

gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("0.2g");
gROOT.ProcessLineSync(".x ../AccxEffCalculator.cxx+")
gROOT.ProcessLineSync(".x ../Binning.cxx+")

fileBinning = TFile.Open("../output/binning.root")
binning = fileBinning.Get("Binning")
CostValues = binning.GetCostValues()
PhiValues = binning.GetPhiValues()

minPt = [0,2,4,6,10]
maxPt = [2,4,6,10,1000]

if os.path.isfile("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root"):
    fileDataMC = TFile.Open("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root")  # lxplus
else:
    fileDataMC = TFile.Open("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root") # local

treeDataMC = fileDataMC.Get("MCTree")

print "########################################################################"
print "Setting lambda theta = 0.8"
lambdaTheta1 = 0.8
AccxEffReWeight1 = AccxEffCalculator(treeDataMC)
AccxEffReWeight1.SetPtBins(5,array('d',[0.,2.,4.,6.,10.]),array('d',[2.,4.,6.,10.,1000.]))
AccxEffReWeight1.SetBinning(CostValues,PhiValues)

if os.path.isfile("../output/AccxEffReWeightedFullStatPol1.root"):
    pass
else:
    AccxEffReWeight1.ReWeightAccxEff(lambdaTheta1,"FullStat","../output/AccxEffReWeightedFullStatPol1.root")

fileAccxEffReWeightedPol1 = TFile.Open("../output/AccxEffReWeightedFullStatPol1.root")

histGenCostReWeightedPol1 = []
histAccxEffCostReWeightedPol1 = []
histAccxEffPhiReWeightedPol1 = []

for i in range(5):
    histGenCostReWeightedPol1.append(fileAccxEffReWeightedPol1.Get('histGenCostReWeighted' + str(lambdaTheta1) + '_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histGenCostReWeightedPol1[i].SetDirectory(0)
    histAccxEffCostReWeightedPol1.append(fileAccxEffReWeightedPol1.Get('histAccxEffCostReWeighted' + str(lambdaTheta1) + '_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffCostReWeightedPol1[i].SetDirectory(0)
    histAccxEffPhiReWeightedPol1.append(fileAccxEffReWeightedPol1.Get('histAccxEffPhiReWeighted' + str(lambdaTheta1) + '_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffPhiReWeightedPol1[i].SetDirectory(0)
fileDataMC.Close()
fileAccxEffReWeightedPol1.Close()


if os.path.isfile("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root"):
    fileDataMC = TFile.Open("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root")  # lxplus
else:
    fileDataMC = TFile.Open("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root") # local

treeDataMC = fileDataMC.Get("MCTree")

print "########################################################################"
print "Setting lambda theta = 0.0"
lambdaTheta2 = 0.0
AccxEffReWeight2 = AccxEffCalculator(treeDataMC)
AccxEffReWeight2.SetPtBins(5,array('d',[0.,2.,4.,6.,10.]),array('d',[2.,4.,6.,10.,1000.]))
AccxEffReWeight2.SetBinning(CostValues,PhiValues)

if os.path.isfile("../output/AccxEffReWeightedFullStatPol2.root"):
    pass
else:
    AccxEffReWeight2.ReWeightAccxEff(lambdaTheta2,"FullStat","../output/AccxEffReWeightedFullStatPol2.root")

fileAccxEffReWeightedPol2 = TFile.Open("../output/AccxEffReWeightedFullStatPol2.root")

histGenCostReWeightedPol2 = []
histAccxEffCostReWeightedPol2 = []
histAccxEffPhiReWeightedPol2 = []

for i in range(5):
    histGenCostReWeightedPol2.append(fileAccxEffReWeightedPol2.Get('histGenCostReWeighted' + str(lambdaTheta2) + '_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histGenCostReWeightedPol2[i].SetDirectory(0)
    histGenCostReWeightedPol2[i].SetLineColor(kRed)
    histAccxEffCostReWeightedPol2.append(fileAccxEffReWeightedPol2.Get('histAccxEffCostReWeighted' + str(lambdaTheta2) + '_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffCostReWeightedPol2[i].SetDirectory(0)
    histAccxEffCostReWeightedPol2[i].SetLineColor(kRed)
    histAccxEffPhiReWeightedPol2.append(fileAccxEffReWeightedPol2.Get('histAccxEffPhiReWeighted' + str(lambdaTheta2) + '_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffPhiReWeightedPol2[i].SetDirectory(0)
    histAccxEffPhiReWeightedPol2[i].SetLineColor(kRed)
fileDataMC.Close()
fileAccxEffReWeightedPol2.Close()

legendAccxEffComp = TLegend(0.35,0.1,0.65,0.3)
legendAccxEffComp.AddEntry(histAccxEffCostReWeightedPol2[1],"#lambda_{#theta} = 0.0","l")
legendAccxEffComp.AddEntry(histAccxEffCostReWeightedPol1[1],"#lambda_{#theta} = 0.8","l")

canvasAccxEffReWeightedComp = TCanvas("canvasAccxEffReWeightedComp","canvasAccxEffReWeightedComp",4,132,1024,768)
canvasAccxEffReWeightedComp.Divide(3,2)

canvasAccxEffReWeightedComp.cd(1)
histAccxEffCostReWeightedPol2[1].Draw()
histAccxEffCostReWeightedPol1[1].Draw("same")
legendAccxEffComp.Draw("same")

canvasAccxEffReWeightedComp.cd(2)
histAccxEffCostReWeightedPol2[2].Draw()
histAccxEffCostReWeightedPol1[2].Draw("same")
legendAccxEffComp.Draw("same")

canvasAccxEffReWeightedComp.cd(3)
histAccxEffCostReWeightedPol2[3].Draw()
histAccxEffCostReWeightedPol1[3].Draw("same")
legendAccxEffComp.Draw("same")

canvasAccxEffReWeightedComp.cd(4)
histAccxEffPhiReWeightedPol2[1].Draw()
histAccxEffPhiReWeightedPol1[1].Draw("same")
legendAccxEffComp.Draw("same")

canvasAccxEffReWeightedComp.cd(5)
histAccxEffPhiReWeightedPol2[2].Draw()
histAccxEffPhiReWeightedPol1[2].Draw("same")
legendAccxEffComp.Draw("same")

canvasAccxEffReWeightedComp.cd(6)
histAccxEffPhiReWeightedPol2[3].Draw()
histAccxEffPhiReWeightedPol1[3].Draw("same")
legendAccxEffComp.Draw("same")

raw_input()
