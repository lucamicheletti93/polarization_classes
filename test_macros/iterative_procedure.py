from ROOT import *
from array import array
import os.path
import math

gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("0.2g");
gROOT.ProcessLineSync(".L ../MathFuncsLib.cxx+")
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
#if os.path.isfile("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root"):
    #fileDataMC = TFile.Open("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root")  # lxplus
#else:
    #fileDataMC = TFile.Open("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root") # local
#treeDataMC = fileDataMC.Get("MCTree")

# To check if a directory exists
namePtRanges = ['2pt4','4pt6','6pt10']
minPtBin = [2.,4.,6.]
maxPtBin = [4.,6.,10.]
for i in range(len(namePtRanges)):
    if not os.path.exists('iterative_procedure/' + namePtRanges[i]):
        os.makedirs('iterative_procedure/' + namePtRanges[i])

#nameOutputFile = ['AccxEffReWeighted1stStep.root','AccxEffReWeighted2ndStep.root','AccxEffReWeighted3rdStep.root']
nameOutputFile = 'AccxEffReWeighted1stStep.root'
# 2 < pT < 4 GeV/c
#lambdaTheta_2pt4 = [-0.189948,-0.20244,-0.203284,-0.203346]
#lambdaPhi_2pt4 = [-0.2228,-0.262088,-0.269251,-0.270491]
# 4 < pT < 6 GeV/c
#lambdaTheta_4pt6 = [-0.0197807]
#lambdaPhi_4pt6 = [-0.00780492]
# 6 < pT < 10 GeV/c
#lambdaTheta_6pt10 = [0.0662521]
#lambdaPhi_6pt10 = [-0.0603516]

lambdaTheta = [-0.189948,-0.0197807,0.0662521]
lambdaPhi = [-0.2228,-0.00780492,-0.0603516]

for i in range(len(namePtRanges)):
    if os.path.isfile("iterative_procedure/" + namePtRanges[i] + "/" + nameOutputFile):
        print "--> " + nameOutputFile + " has already produced"
        print "if you want to reproduce it delete " + nameOutputFile + " and re-run"
    else:
        print "--> iterative_procedure/" + namePtRanges[i] + "/" + nameOutputFile + " does not exist"
        if os.path.isfile("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root"):
            fileDataMC = TFile.Open("/afs/cern.ch/user/l/lmichele/CERNBox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root")  # lxplus
        else:
            fileDataMC = TFile.Open("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/JIRA_TICKET/READ_MC/OUTPUT/MC_official_tree_Jpsi_PbPb_Nopol.root") # local
        treeDataMC = fileDataMC.Get("MCTree")
        AccxEffReWeight = AccxEffCalculator(treeDataMC)
        AccxEffReWeight.SetPtBins(1,array('d',[minPtBin[i]]),array('d',[maxPtBin[i]]))
        AccxEffReWeight.SetBinning(CostValues,PhiValues)
        AccxEffReWeight.ReWeightAccxEff(lambdaTheta[i],lambdaPhi[i],"FullStat",kTRUE,"iterative_procedure/" + namePtRanges[i] + "/" + nameOutputFile)
        del AccxEffReWeight
        fileDataMC.Close()
        #AccxEffReWeight.ReWeightAccxEff(-0.189948,-0.2228,"FullStat",kTRUE,"iterative_procedure/" + nameOutputFile)    # 1st iteration
        #AccxEffReWeight.ReWeightAccxEff(-0.20244,-0.262088,"FullStat",kTRUE,"iterative_procedure/" + nameOutputFile)   # 2nd iteration
        #AccxEffReWeight.ReWeightAccxEff(-0.203284,-0.269251,"FullStat",kTRUE,"iterative_procedure/" + nameOutputFile)  # 3rd iteration
        #AccxEffReWeight.ReWeightAccxEff(-0.203346,-0.270491,"FullStat",kTRUE,"iterative_procedure/" + nameOutputFile)  # 4th iteration (to be performed)

################################################################################
#histAccxEffCostReWeighted = []
#histAccxEffPhiReWeighted = []
#histGenCostReWeighted = []
#histGenPhiReWeighted = []
#funcCost = []
#funcPhi = []

#for i in range(3):
    #fileAccxEffReWeight = TFile.Open("iterative_procedure/" + nameOutputFile[i])
    #histAccxEffCostReWeighted.append(fileAccxEffReWeight.Get("histAccxEffCostReWeighted_2pT4"))
    #histAccxEffCostReWeighted[i].SetLineColor(i+2)
    #histAccxEffCostReWeighted[i].SetDirectory(0)
    #histAccxEffPhiReWeighted.append(fileAccxEffReWeight.Get("histAccxEffPhiReWeighted_2pT4"))
    #histAccxEffPhiReWeighted[i].SetLineColor(i+2)
    #histAccxEffPhiReWeighted[i].SetDirectory(0)
    #histGenCostReWeighted.append(fileAccxEffReWeight.Get("histGenCostReWeighted_2pT4"))
    #histGenCostReWeighted[i].SetLineColor(i+2)
    #histGenCostReWeighted[i].SetDirectory(0)
    #histGenPhiReWeighted.append(fileAccxEffReWeight.Get("histGenPhiReWeighted_2pT4"))
    #histGenPhiReWeighted[i].SetLineColor(i+2)
    #histGenPhiReWeighted[i].SetDirectory(0)
    #fileAccxEffReWeight.Close()

#canvasAccxEffCost = TCanvas("canvasAccxEffCost","canvasAccxEffCost",20,20,600,600)
#histAccxEffCostReWeighted[0].Draw("E")
#histAccxEffCostReWeighted[1].Draw("Esame")
#histAccxEffCostReWeighted[2].Draw("Esame")

#canvasAccxEffPhi = TCanvas("canvasAccxEffPhi","canvasAccxEffPhi",20,20,600,600)
#histAccxEffPhiReWeighted[0].Draw("E")
#histAccxEffPhiReWeighted[1].Draw("Esame")
#histAccxEffPhiReWeighted[2].Draw("Esame")

# Check plots
#for i in range(3):
    #for j in range(19):
        #histGenCostReWeighted[i].SetBinContent(j+1,histGenCostReWeighted[i].GetBinContent(j+1)/CostWidth[j])
        #histGenCostReWeighted[i].SetBinError(j+1,histGenCostReWeighted[i].GetBinError(j+1)/CostWidth[j])
        #histGenCostReWeighted[i].GetXaxis().SetTitle("cos#it{#theta}_{HE}")
    #for j in range(10):
        #histGenPhiReWeighted[i].SetBinContent(j+1,histGenPhiReWeighted[i].GetBinContent(j+1)/PhiWidth[j])
        #histGenPhiReWeighted[i].SetBinError(j+1,histGenPhiReWeighted[i].GetBinError(j+1)/PhiWidth[j])
        #histGenPhiReWeighted[i].GetXaxis().SetTitle("#it{#varphi}_{HE}")

#for i in range(3):
    #funcCost.append(TF1("funcCost",MathFuncs.MyMathFuncs.MyFuncPolCosTheta,-1.,1.,2))
    #funcCost[i].SetLineColor(i+2)
    #histGenCostReWeighted[i].Fit(funcCost[i],"R0")
    #funcPhi.append(TF1("funcPhi",MathFuncs.MyMathFuncs.MyFuncPolPhi,0,PI,3))
    #funcPhi[i].SetLineColor(i+2)
    #histGenPhiReWeighted[i].Fit(funcPhi[i],"R0")

#canvasGenCost = TCanvas("canvasGenCost","canvasGenCost",20,20,600,600)
#histGenCostReWeighted[0].Draw("E")
#funcCost[0].Draw("Esame")
#for i in range(2):
    #histGenCostReWeighted[i+1].Draw("Esame")
    #funcCost[i+1].Draw("Esame")
#legendGenCost = TLegend(0.35,0.1,0.65,0.3)
#legendGenCost.AddEntry(funcCost[0],"1^{st} itreration","l")
#legendGenCost.AddEntry(funcCost[1],"2^{nd} itreration","l")
#legendGenCost.AddEntry(funcCost[2],"3^{rd} itreration","l")
#legendGenCost.Draw("same")


#canvasGenPhi = TCanvas("canvasGenPhi","canvasGenPhi",20,20,600,600)
#histGenPhiReWeighted[0].Draw("E")
#funcPhi[0].Draw("Esame")
#for i in range(2):
    #histGenPhiReWeighted[i+1].Draw("Esame")
    #funcPhi[i+1].Draw("Esame")
#legendGenPhi = TLegend(0.35,0.1,0.65,0.3)
#legendGenPhi.AddEntry(funcPhi[0],"1^{st} itreration","l")
#legendGenPhi.AddEntry(funcPhi[1],"2^{nd} itreration","l")
#legendGenPhi.AddEntry(funcPhi[2],"3^{rd} itreration","l")
#legendGenPhi.Draw("same")

raw_input()
