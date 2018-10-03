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

AccxEff = AccxEffCalculator(treeDataMC)
AccxEff.SetPtBins(5,array('d',[0.,2.,4.,6.,10.]),array('d',[2.,4.,6.,10.,1000.]))
AccxEff.SetBinning(CostValues,PhiValues)
AccxEff.ComputeAccxEff("TestStat","output/AccxEff.root")

fileAccxEff = TFile.Open("output/AccxEff.root")
histRecCost = []
histAccxEffCost = []
histAccxEffPhi = []
histAccxEffCostPhi = []
histAccxEffCostPhiStatRel = []

for i in range(5):
    histRecCost.append(fileAccxEff.Get('histRecCost_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffCost.append(fileAccxEff.Get('histAccxEffCost_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffPhi.append(fileAccxEff.Get('histAccxEffPhi_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffCostPhi.append(fileAccxEff.Get('histAccxEffCostPhi_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histAccxEffCostPhiStatRel.append(fileAccxEff.Get('histAccxEffCostPhiStatRel_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    histRecCost[i].SetDirectory(0)
    histAccxEffCost[i].SetDirectory(0)
    histAccxEffPhi[i].SetDirectory(0)
    histAccxEffCostPhi[i].SetDirectory(0)
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

for i in range(AccxEff.GetNCostBins()):
    #print histAccxEffCost[3].GetBinContent(i+1), "+-", histAccxEffCost[3].GetBinError(i+1), "[", histAccxEffCost[3].GetBinError(i+1)/histAccxEffCost[3].GetBinContent(i+1)*100, "]"
    print histRecCost[3].GetBinContent(i+1), "+-", histRecCost[3].GetBinError(i+1), "[", histRecCost[3].GetBinError(i+1)/histRecCost[3].GetBinContent(i+1)*100, "]"
canvasAccxEffCost = TCanvas("canvasAccxEffCost","canvasAccxEffCost",20,20,600,600)
histAccxEffCost[3].Draw()

canvasAccxEffPhi = TCanvas("canvasAccxEffPhi","canvasAccxEffPhi",20,20,600,600)
histAccxEffPhi[3].Draw()

canvasAccxEffCostPhi = TCanvas("canvasAccxEffCostPhi","canvasAccxEffCostPhi",20,20,600,600)
histAccxEffCostPhi[3].Draw("COLZtexterror")

canvasAccxEffCostPhiStatRel = TCanvas("canvasAccxEffCostPhiStatRel","canvasAccxEffCostPhiStatRel",20,20,600,600)
histAccxEffCostPhiStatRel[3].Draw("COLZtext")

canvasAccxEffCostPt = TCanvas("canvasAccxEffCostPt","canvasAccxEffCostPt",20,20,600,600)
histAccxEffCostPt.Draw("COLZ")

canvasAccxEffPhiPt = TCanvas("canvasAccxEffPhiPt","canvasAccxEffPhiPt",20,20,600,600)
histAccxEffPhiPt.Draw("COLZ")

histGenCostPt = histGenCostPt.ProjectionY("histGenPt")
canvasGenPt = TCanvas("canvasGenPt","canvasGenPt",20,20,600,600)
histGenPt.Draw()

################################################################################

#AccxEffReWeighted1 = AccxEffCalculator(treeDataMC)
#AccxEffReWeighted1.ReWeightAccxEff(0.8,"TestStat","pol08prova.root")

#fileAccxEffReWeightedpol08prova = TFile.Open("pol08prova.root")

#histAccxEffReWeightedpol08prova = []
#histGenReWeightedpol08prova = []
#histRecReWeightedpol08prova = []
#for i in range(12):
    #histAccxEffReWeightedpol08prova.append(fileAccxEffReWeightedpol08prova.Get('histAccxEffCostReWeighted0.8_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    #histAccxEffReWeightedpol08prova[i].SetDirectory(0)
    #histGenReWeightedpol08prova.append(fileAccxEffReWeightedpol08prova.Get('histGenCostReWeighted0.8_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    #histGenReWeightedpol08prova[i].SetDirectory(0)
    #histRecReWeightedpol08prova.append(fileAccxEffReWeightedpol08prova.Get('histRecCostReWeighted0.8_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    #histRecReWeightedpol08prova[i].SetDirectory(0)
#fileAccxEffReWeightedpol08prova.Close()
#fileDataMC.Close()


#fileDataMC = TFile.Open("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root")
#treeDataMC = fileDataMC.Get("MCTree")

#AccxEffReWeighted2 = AccxEffCalculator(treeDataMC)
#AccxEffReWeighted2.ReWeightAccxEff(0.,"TestStat","pol00prova.root")

#fileAccxEffReWeightedpol00prova = TFile.Open("pol00prova.root")

#histAccxEffReWeightedpol00prova = []
#histGenReWeightedpol00prova = []
#histRecReWeightedpol00prova = []
#for i in range(12):
    #histAccxEffReWeightedpol00prova.append(fileAccxEffReWeightedpol00prova.Get('histAccxEffCostReWeighted0.0_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    #histAccxEffReWeightedpol00prova[i].SetDirectory(0)
    #histGenReWeightedpol00prova.append(fileAccxEffReWeightedpol00prova.Get('histGenCostReWeighted0.0_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    #histGenReWeightedpol00prova[i].SetDirectory(0)
    #histRecReWeightedpol00prova.append(fileAccxEffReWeightedpol00prova.Get('histRecCostReWeighted0.0_' + str(minPt[i]) + 'pT' + str(maxPt[i])))
    #histRecReWeightedpol00prova[i].SetDirectory(0)
#fileAccxEffReWeightedpol00prova.Close()


#canvasGenReWeighted = TCanvas("canvasGenReWeighted","canvasGenReWeighted",20,20,600,600)
#canvasGenReWeighted.Divide(6,2)

#for i in range(12):
    #canvasGenReWeighted.cd(i+1)
    #histGenReWeightedpol08prova[i].Draw()
    #histGenReWeightedpol00prova[i].SetLineColor(2)
    #histGenReWeightedpol00prova[i].Draw("same")

#canvasRecReWeighted = TCanvas("canvasRecReWeighted","canvasRecReWeighted",20,20,600,600)
#canvasRecReWeighted.Divide(6,2)

#for i in range(12):
    #canvasRecReWeighted.cd(i+1)
    #histRecReWeightedpol08prova[i].Draw()
    #histRecReWeightedpol00prova[i].SetLineColor(2)
    #histRecReWeightedpol00prova[i].Draw("same")

#canvasAccxEffReWeighted = TCanvas("canvasAccxEffReWeighted","canvasAccxEffReWeighted",20,20,600,600)
#canvasAccxEffReWeighted.Divide(6,2)

#for i in range(12):
    #canvasAccxEffReWeighted.cd(i+1)
    #histAccxEffReWeightedpol08prova[i].Draw()
    #histAccxEffReWeightedpol00prova[i].SetLineColor(2)
    #histAccxEffReWeightedpol00prova[i].Draw("same")

################################################################################

#for i in range(12):
    #hist.append(fileAccxEff.Get("histAccxEffCost_%ipT%i",minPt[i],maxPt[i])))

#canvas = TCanvas("canvas","canvas",20,20,600,600)
#canvas.Divide(2,6)


#hist = []
#for i in range(1):
    #hist.append(TH1D("hist","hist",64,-4,4))

#histGen = AccxEff.GetHistGenCost()

#canvasGen = TCanvas("canvas","canvas",20,20,600,600)
#histGen.Draw()



#AccxEffReWeighted1 = AccxEffReWeighting(treeDataMC)
#AccxEffReWeighted1.ReWeightAccxEff(-0.5,"TestStat")

#AccxEffReWeighted1 = AccxEffReWeighting(treeDataMC)
#AccxEffReWeighted1.ReWeightAccxEff(-0.5,"TestStat")

#histGenReWeighted1 = AccxEffReWeighted1.GetHistGenCostReWeighted()
#histRecReWeighted1 = AccxEffReWeighted1.GetHistRecCostReWeighted()
#histAccxEffReWeighted1 = AccxEffReWeighted1.GetHistAccxEffCostReWeighted()
#histGenReWeighted1.SetLineColor(kBlue)
#histRecReWeighted1.SetLineColor(kBlue)
#histAccxEffReWeighted1.SetLineColor(kBlue)

#AccxEffReWeighted2 = AccxEffReWeighting(treeDataMC)
#AccxEffReWeighted2.ReWeightAccxEff(0,"TestStat")

#histGenReWeighted2 = AccxEffReWeighted2.GetHistGenCostReWeighted()
#histRecReWeighted2 = AccxEffReWeighted2.GetHistRecCostReWeighted()
#histAccxEffReWeighted2 = AccxEffReWeighted2.GetHistAccxEffCostReWeighted()
#histGenReWeighted2.SetLineColor(kRed)
#histRecReWeighted2.SetLineColor(kRed)
#histAccxEffReWeighted2.SetLineColor(kRed)

#canvasGenRec = TCanvas("canvas","canvas",20,20,600,600)
#canvasGenRec.Divide(1,2)

#canvasGenRec.cd(1)
#histGenReWeighted1.Draw()
#histGenReWeighted2.Draw("same")

#canvasGenRec.cd(2)
#histRecReWeighted1.Draw()
#histRecReWeighted2.Draw("same")

#canvasAccxEff = TCanvas("canvasAccxEff","canvasAccxEff",20,20,600,600)
#histAccxEffReWeighted1.Draw()
#histAccxEffReWeighted2.Draw("same")

raw_input()
