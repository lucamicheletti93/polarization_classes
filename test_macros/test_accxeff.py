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

if os.path.isfile("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root"):
    fileDataMC = TFile.Open("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root")  # lxplus
else:
    fileDataMC = TFile.Open("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root") # local

treeDataMC = fileDataMC.Get("MCTree")

minPt = [0,2,4,6,10]
maxPt = [2,4,6,10,1000]

AccxEff = AccxEffCalculator(treeDataMC)
AccxEff.SetPtBins(5,array('d',[0.,2.,4.,6.,10.]),array('d',[2.,4.,6.,10.,1000.]))
AccxEff.SetBinning(CostValues,PhiValues)

if os.path.isfile("../output/AccxEffFullStat.root"):
    print "The file ../output/AccxEffFullStat.root already exists, proceed with the following part of the macro"

else:
    AccxEff.ComputeAccxEff("FullStat","../output/AccxEffFullStat.root")

fileAccxEff = TFile.Open("../output/AccxEffFullStat.root")
histRecCost = []
histAccxEffCost = []
histAccxEffPhi = []
histAccxEffPhiTilde = []
histAccxEffCostPhi = []
histAccxEffCostPhiStatRel = []

for i in range(5):
    histRecCost.append(fileAccxEff.Get('histRecCost_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffCost.append(fileAccxEff.Get('histAccxEffCost_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffPhi.append(fileAccxEff.Get('histAccxEffPhi_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffPhiTilde.append(fileAccxEff.Get('histAccxEffPhiTilde_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffCostPhi.append(fileAccxEff.Get('histAccxEffCostPhi_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffCostPhiStatRel.append(fileAccxEff.Get('histAccxEffCostPhiStatRel_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histRecCost[i].SetDirectory(0)
    histRecCost[i].GetXaxis().SetTitle("cos#it{#theta}_{HE}")
    histAccxEffCost[i].SetDirectory(0)
    histAccxEffCost[i].GetXaxis().SetTitle("cos#it{#theta}_{HE}")
    histAccxEffPhi[i].SetDirectory(0)
    histAccxEffPhi[i].GetXaxis().SetTitle("#it{#varphi}_{HE}")
    histAccxEffPhiTilde[i].SetDirectory(0)
    histAccxEffPhiTilde[i].GetXaxis().SetTitle("#tilde{#it{#varphi}}_{HE}")
    histAccxEffCostPhi[i].SetDirectory(0)
    histAccxEffCostPhi[i].GetXaxis().SetTitle("#it{#varphi}_{HE}")
    histAccxEffCostPhiStatRel[i].SetDirectory(0)

histAccxEffCostPt = fileAccxEff.Get('histAccxEffCostPt')
histAccxEffCostPt.GetXaxis().SetTitle("cos#it{#theta}_{HE}")
histAccxEffCostPt.GetYaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
histAccxEffCostPt.SetDirectory(0)
histAccxEffPhiPt = fileAccxEff.Get('histAccxEffPhiPt')
histAccxEffPhiPt.GetXaxis().SetTitle("#it{#varphi}_{HE}")
histAccxEffPhiPt.GetYaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
histAccxEffPhiPt.SetDirectory(0)
histGenCostPt = fileAccxEff.Get("histGenCostPt")
histGenCostPt.SetDirectory(0)
fileAccxEff.Close()

print "########################################################################"
print "Error comparison between 1D and 2D in 6 < pT < 10 Gev/c"
for i in range(AccxEff.GetNCostBins()):
    print histRecCost[3].GetBinContent(i+1), "+-", histRecCost[3].GetBinError(i+1), "[", histRecCost[3].GetBinError(i+1)/histRecCost[3].GetBinContent(i+1)*100, "]"
canvasAccxEffCost = TCanvas("canvasAccxEffCost","canvasAccxEffCost",20,20,600,600)
histAccxEffCost[3].Draw()

canvasAccxEffPhi = TCanvas("canvasAccxEffPhi","canvasAccxEffPhi",20,20,600,600)
histAccxEffPhi[3].Draw()

canvasAccxEffCostPhi = TCanvas("canvasAccxEffCostPhi","canvasAccxEffCostPhi",20,20,600,600)
histAccxEffCostPhi[3].Draw("COLZtexterror")

print "########################################################################"
print "Acc x Eff vs (Cos(theta),pT)"
canvasAccxEffCostPt = TCanvas("canvasAccxEffCostPt","canvasAccxEffCostPt",20,20,600,600)
histAccxEffCostPt.Draw("COLZ")

canvasAccxEffPhiPt = TCanvas("canvasAccxEffPhiPt","canvasAccxEffPhiPt",20,20,600,600)
histAccxEffPhiPt.Draw("COLZ")

print "Acc x Eff vs (Cos(theta),Phi) compared at different pT"
canvasAccxEffCostComp = TCanvas("canvasAccxEffCostComp","canvasAccxEffCostComp",20,20,600,600)
histAccxEffCost[3].SetLineColor(kBlue)
histAccxEffCost[3].Draw()
histAccxEffCost[2].SetLineColor(kRed)
histAccxEffCost[2].Draw("same")
histAccxEffCost[1].SetLineColor(kGreen)
histAccxEffCost[1].Draw("same")

legendAccxEffCost = TLegend(0.35,0.1,0.65,0.3)
legendAccxEffCost.AddEntry(histAccxEffCost[1],"2 < #it{p}_{T} < 4 GeV/#it{c}","l")
legendAccxEffCost.AddEntry(histAccxEffCost[2],"4 < #it{p}_{T} < 6 GeV/#it{c}","l")
legendAccxEffCost.AddEntry(histAccxEffCost[3],"6 < #it{p}_{T} < 10 GeV/#it{c}","l")
legendAccxEffCost.Draw("same")

canvasAccxEffPhiComp = TCanvas("canvasAccxEffPhiComp","canvasAccxEffPhiComp",20,20,600,600)
histAccxEffPhi[3].SetLineColor(kBlue)
histAccxEffPhi[3].Draw()
histAccxEffPhi[2].SetLineColor(kRed)
histAccxEffPhi[2].Draw("same")
histAccxEffPhi[1].SetLineColor(kGreen)
histAccxEffPhi[1].Draw("same")

legendAccxEffPhi = TLegend(0.35,0.1,0.65,0.3)
legendAccxEffPhi.AddEntry(histAccxEffPhi[1],"2 < #it{p}_{T} < 4 GeV/#it{c}","l")
legendAccxEffPhi.AddEntry(histAccxEffPhi[2],"4 < #it{p}_{T} < 6 GeV/#it{c}","l")
legendAccxEffPhi.AddEntry(histAccxEffPhi[3],"6 < #it{p}_{T} < 10 GeV/#it{c}","l")
legendAccxEffPhi.Draw("same")

canvasAccxEffPhiTildeComp = TCanvas("canvasAccxEffPhiTildeComp","canvasAccxEffPhiTildeComp",20,20,600,600)
histAccxEffPhiTilde[3].SetLineColor(kBlue)
histAccxEffPhiTilde[3].Draw()
histAccxEffPhiTilde[2].SetLineColor(kRed)
histAccxEffPhiTilde[2].Draw("same")
histAccxEffPhiTilde[1].SetLineColor(kGreen)
histAccxEffPhiTilde[1].Draw("same")

legendAccxEffPhiTilde = TLegend(0.35,0.1,0.65,0.3)
legendAccxEffPhiTilde.AddEntry(histAccxEffPhiTilde[1],"2 < #it{p}_{T} < 4 GeV/#it{c}","l")
legendAccxEffPhiTilde.AddEntry(histAccxEffPhiTilde[2],"4 < #it{p}_{T} < 6 GeV/#it{c}","l")
legendAccxEffPhiTilde.AddEntry(histAccxEffPhiTilde[3],"6 < #it{p}_{T} < 10 GeV/#it{c}","l")
legendAccxEffPhiTilde.Draw("same")

raw_input()
