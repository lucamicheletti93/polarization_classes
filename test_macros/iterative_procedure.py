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

namePtRanges = ['2pT4','4pT6','6pT10']

histGenCost = []; histGenPhi = [];
fileAccxEff = TFile.Open("output/AccxEffFullStat.root")
for i in range(len(namePtRanges)):
    histGenCost.append(fileAccxEff.Get("histGenCost_" + namePtRanges[i]))
    histGenPhi.append(fileAccxEff.Get("histGenPhi_" + namePtRanges[i]))

for i in range(len(namePtRanges)):
    for j in range(19):
        histGenCost[i].SetBinContent(j+1,histGenCost[i].GetBinContent(j+1)/CostWidth[j])
        histGenCost[i].SetBinError(j+1,histGenCost[i].GetBinError(j+1)/CostWidth[j])
        histGenCost[i].Scale(1./(histGenCost[i].Integral()))
    for j in range(10):
        histGenPhi[i].SetBinContent(j+1,histGenPhi[i].GetBinContent(j+1)/PhiWidth[j])
        histGenPhi[i].SetBinError(j+1,histGenPhi[i].GetBinError(j+1)/PhiWidth[j])
        histGenPhi[i].Scale(1./(histGenPhi[i].Integral()))

funcCost = []; funcPhi = [];
for i in range(len(namePtRanges)):
    funcCost.append(TF1("funcCost",MathFuncs.MyMathFuncs.MyFuncPolCosTheta,-1.,1.,2))
    funcCost[i].SetLineColor(kBlack)
    funcCost[i].SetLineStyle(kDashed)
    histGenCost[i].Fit(funcCost[i],"R0")

    funcPhi.append(TF1("funcPhi",MathFuncs.MyMathFuncs.MyFuncPolPhi,0,PI,3))
    funcPhi[i].SetLineColor(kBlack)
    funcPhi[i].SetLineStyle(kDashed)
    histGenPhi[i].Fit(funcPhi[i],"R0")

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
minPtBin = [2.,4.,6.]
maxPtBin = [4.,6.,10.]
for i in range(len(namePtRanges)):
    if not os.path.exists('iterative_procedure/' + namePtRanges[i]):
        os.makedirs('iterative_procedure/' + namePtRanges[i])

#nameOutputFile = ['AccxEffReWeighted1stStep.root','AccxEffReWeighted2ndStep.root','AccxEffReWeighted3rdStep.root']
nameOutputFile = 'AccxEffReWeighted3rdStep.root'
# 2 < pT < 4 GeV/c
lambdaTheta_2pt4 = [-0.189948,-0.20244,-0.203284,-0.203346]
lambdaPhi_2pt4 = [-0.2228,-0.262088,-0.269251,-0.270491]
# 4 < pT < 6 GeV/c
lambdaTheta_4pt6 = [-0.0197807,-0.0210581,-0.0211456,-0.0211515]
lambdaPhi_4pt6 = [-0.00780492,-0.00874564,-0.00886213,-0.00887629]
# 6 < pT < 10 GeV/c
lambdaTheta_6pt10 = [0.0662521,0.069291,0.0694284,0.0694346]
lambdaPhi_6pt10 = [-0.0603516,-0.064812,-0.0651313,-0.0651542]

lambdaTheta = [-0.203284,-0.0211456,0.0694284]
lambdaPhi = [-0.269251,-0.00886213,-0.0651313]

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
        #AccxEffReWeight.ReWeightAccxEff(lambdaTheta[i],lambdaPhi[i],"FullStat",kTRUE,"iterative_procedure/" + namePtRanges[i] + "/" + nameOutputFile)
        del AccxEffReWeight
        fileDataMC.Close()

################################################################################
histAccxEffCostReWeighted_2pT4 = []; histAccxEffCostReWeighted_4pT6 = []; histAccxEffCostReWeighted_6pT10 = [];
histAccxEffPhiReWeighted_2pT4 = []; histAccxEffPhiReWeighted_4pT6 = []; histAccxEffPhiReWeighted_6pT10 = [];
histGenCostReWeighted_2pT4 = []; histGenCostReWeighted_4pT6 = []; histGenCostReWeighted_6pT10 = [];
histGenPhiReWeighted_2pT4 = []; histGenPhiReWeighted_4pT6 = []; histGenPhiReWeighted_6pT10 = [];

listNameOutputFile = ['AccxEffReWeighted1stStep.root','AccxEffReWeighted2ndStep.root','AccxEffReWeighted3rdStep.root']
nIterations = 3
nPtBins = 3
index = 0

for i in range(nIterations):
    fileAccxEffReWeight = TFile.Open("iterative_procedure/2pT4/" +  listNameOutputFile[i])
    histAccxEffCostReWeighted_2pT4.append(fileAccxEffReWeight.Get("histAccxEffCostReWeighted_2pT4"))
    histAccxEffCostReWeighted_2pT4[i].SetLineColor(i+2)
    histAccxEffCostReWeighted_2pT4[i].SetDirectory(0)
    histAccxEffPhiReWeighted_2pT4.append(fileAccxEffReWeight.Get("histAccxEffPhiReWeighted_2pT4"))
    histAccxEffPhiReWeighted_2pT4[i].SetLineColor(i+2)
    histAccxEffPhiReWeighted_2pT4[i].SetDirectory(0)
    histGenCostReWeighted_2pT4.append(fileAccxEffReWeight.Get("histGenCostReWeighted_2pT4"))
    histGenCostReWeighted_2pT4[i].SetLineColor(i+2)
    histGenCostReWeighted_2pT4[i].SetDirectory(0)
    histGenPhiReWeighted_2pT4.append(fileAccxEffReWeight.Get("histGenPhiReWeighted_2pT4"))
    histGenPhiReWeighted_2pT4[i].SetLineColor(i+2)
    histGenPhiReWeighted_2pT4[i].SetDirectory(0)
    fileAccxEffReWeight.Close()

for i in range(nIterations):
    fileAccxEffReWeight = TFile.Open("iterative_procedure/4pT6/" +  listNameOutputFile[i])
    histAccxEffCostReWeighted_4pT6.append(fileAccxEffReWeight.Get("histAccxEffCostReWeighted_4pT6"))
    histAccxEffCostReWeighted_4pT6[i].SetLineColor(i+2)
    histAccxEffCostReWeighted_4pT6[i].SetDirectory(0)
    histAccxEffPhiReWeighted_4pT6.append(fileAccxEffReWeight.Get("histAccxEffPhiReWeighted_4pT6"))
    histAccxEffPhiReWeighted_4pT6[i].SetLineColor(i+2)
    histAccxEffPhiReWeighted_4pT6[i].SetDirectory(0)
    histGenCostReWeighted_4pT6.append(fileAccxEffReWeight.Get("histGenCostReWeighted_4pT6"))
    histGenCostReWeighted_4pT6[i].SetLineColor(i+2)
    histGenCostReWeighted_4pT6[i].SetDirectory(0)
    histGenPhiReWeighted_4pT6.append(fileAccxEffReWeight.Get("histGenPhiReWeighted_4pT6"))
    histGenPhiReWeighted_4pT6[i].SetLineColor(i+2)
    histGenPhiReWeighted_4pT6[i].SetDirectory(0)
    fileAccxEffReWeight.Close()

for i in range(nIterations):
    fileAccxEffReWeight = TFile.Open("iterative_procedure/6pT10/" +  listNameOutputFile[i])
    histAccxEffCostReWeighted_6pT10.append(fileAccxEffReWeight.Get("histAccxEffCostReWeighted_6pT10"))
    histAccxEffCostReWeighted_6pT10[i].SetLineColor(i+2)
    histAccxEffCostReWeighted_6pT10[i].SetDirectory(0)
    histAccxEffPhiReWeighted_6pT10.append(fileAccxEffReWeight.Get("histAccxEffPhiReWeighted_6pT10"))
    histAccxEffPhiReWeighted_6pT10[i].SetLineColor(i+2)
    histAccxEffPhiReWeighted_6pT10[i].SetDirectory(0)
    histGenCostReWeighted_6pT10.append(fileAccxEffReWeight.Get("histGenCostReWeighted_6pT10"))
    histGenCostReWeighted_6pT10[i].SetLineColor(i+2)
    histGenCostReWeighted_6pT10[i].SetDirectory(0)
    histGenPhiReWeighted_6pT10.append(fileAccxEffReWeight.Get("histGenPhiReWeighted_6pT10"))
    histGenPhiReWeighted_6pT10[i].SetLineColor(i+2)
    histGenPhiReWeighted_6pT10[i].SetDirectory(0)
    fileAccxEffReWeight.Close()

# Check plots
for i in range(nIterations):
    for j in range(19):
        histGenCostReWeighted_2pT4[i].SetBinContent(j+1,histGenCostReWeighted_2pT4[i].GetBinContent(j+1)/CostWidth[j])
        histGenCostReWeighted_2pT4[i].SetBinError(j+1,histGenCostReWeighted_2pT4[i].GetBinError(j+1)/CostWidth[j])
        histGenCostReWeighted_2pT4[i].Scale(1./(histGenCostReWeighted_2pT4[i].Integral()))
        histGenCostReWeighted_2pT4[i].GetXaxis().SetTitle("cos#it{#theta}_{HE}")
        histGenCostReWeighted_4pT6[i].SetBinContent(j+1,histGenCostReWeighted_4pT6[i].GetBinContent(j+1)/CostWidth[j])
        histGenCostReWeighted_4pT6[i].SetBinError(j+1,histGenCostReWeighted_4pT6[i].GetBinError(j+1)/CostWidth[j])
        histGenCostReWeighted_4pT6[i].Scale(1./(histGenCostReWeighted_4pT6[i].Integral()))
        histGenCostReWeighted_4pT6[i].GetXaxis().SetTitle("cos#it{#theta}_{HE}")
        histGenCostReWeighted_6pT10[i].SetBinContent(j+1,histGenCostReWeighted_6pT10[i].GetBinContent(j+1)/CostWidth[j])
        histGenCostReWeighted_6pT10[i].SetBinError(j+1,histGenCostReWeighted_6pT10[i].GetBinError(j+1)/CostWidth[j])
        histGenCostReWeighted_6pT10[i].Scale(1./(histGenCostReWeighted_6pT10[i].Integral()))
        histGenCostReWeighted_6pT10[i].GetXaxis().SetTitle("cos#it{#theta}_{HE}")
    for j in range(10):
        histGenPhiReWeighted_2pT4[i].SetBinContent(j+1,histGenPhiReWeighted_2pT4[i].GetBinContent(j+1)/PhiWidth[j])
        histGenPhiReWeighted_2pT4[i].SetBinError(j+1,histGenPhiReWeighted_2pT4[i].GetBinError(j+1)/PhiWidth[j])
        histGenPhiReWeighted_2pT4[i].Scale(1./(histGenPhiReWeighted_2pT4[i].Integral()))
        histGenPhiReWeighted_2pT4[i].GetXaxis().SetTitle("#it{#varphi}_{HE}")
        histGenPhiReWeighted_4pT6[i].SetBinContent(j+1,histGenPhiReWeighted_4pT6[i].GetBinContent(j+1)/PhiWidth[j])
        histGenPhiReWeighted_4pT6[i].SetBinError(j+1,histGenPhiReWeighted_4pT6[i].GetBinError(j+1)/PhiWidth[j])
        histGenPhiReWeighted_4pT6[i].Scale(1./(histGenPhiReWeighted_4pT6[i].Integral()))
        histGenPhiReWeighted_4pT6[i].GetXaxis().SetTitle("#it{#varphi}_{HE}")
        histGenPhiReWeighted_6pT10[i].SetBinContent(j+1,histGenPhiReWeighted_6pT10[i].GetBinContent(j+1)/PhiWidth[j])
        histGenPhiReWeighted_6pT10[i].SetBinError(j+1,histGenPhiReWeighted_6pT10[i].GetBinError(j+1)/PhiWidth[j])
        histGenPhiReWeighted_6pT10[i].Scale(1./(histGenPhiReWeighted_6pT10[i].Integral()))
        histGenPhiReWeighted_6pT10[i].GetXaxis().SetTitle("#it{#varphi}_{HE}")

funcCost_2pT4 = []; funcCost_4pT6 = []; funcCost_6pT10 = [];
funcPhi_2pT4 = []; funcPhi_4pT6 = []; funcPhi_6pT10 = [];

for i in range(nIterations):
    funcCost_2pT4.append(TF1("funcCost_2pT4",MathFuncs.MyMathFuncs.MyFuncPolCosTheta,-1.,1.,2))
    funcCost_2pT4[i].SetLineColor(i+2)
    histGenCostReWeighted_2pT4[i].Fit(funcCost_2pT4[i],"R0")
    funcCost_4pT6.append(TF1("funcCost_4pT6",MathFuncs.MyMathFuncs.MyFuncPolCosTheta,-1.,1.,2))
    funcCost_4pT6[i].SetLineColor(i+2)
    histGenCostReWeighted_4pT6[i].Fit(funcCost_4pT6[i],"R0")
    funcCost_6pT10.append(TF1("funcCost_6pT10",MathFuncs.MyMathFuncs.MyFuncPolCosTheta,-1.,1.,2))
    funcCost_6pT10[i].SetLineColor(i+2)
    histGenCostReWeighted_6pT10[i].Fit(funcCost_6pT10[i],"R0")

    funcPhi_2pT4.append(TF1("funcPhi_2pT4",MathFuncs.MyMathFuncs.MyFuncPolPhi,0,PI,3))
    funcPhi_2pT4[i].SetLineColor(i+2)
    histGenPhiReWeighted_2pT4[i].Fit(funcPhi_2pT4[i],"R0")
    funcPhi_4pT6.append(TF1("funcPhi_4pT6",MathFuncs.MyMathFuncs.MyFuncPolPhi,0,PI,3))
    funcPhi_4pT6[i].SetLineColor(i+2)
    histGenPhiReWeighted_4pT6[i].Fit(funcPhi_4pT6[i],"R0")
    funcPhi_6pT10.append(TF1("funcPhi_6pT10",MathFuncs.MyMathFuncs.MyFuncPolPhi,0,PI,3))
    funcPhi_6pT10[i].SetLineColor(i+2)
    histGenPhiReWeighted_6pT10[i].Fit(funcPhi_6pT10[i],"R0")

legendGenCost = TLegend(0.3,0.1,0.7,0.3)
legendGenCost.AddEntry(funcCost[0],"0 iteration","l")
legendGenCost.AddEntry(funcCost_2pT4[0],"1^{st} iteration","l")
legendGenCost.AddEntry(funcCost_2pT4[1],"2^{nd} iteration","l")
legendGenCost.AddEntry(funcCost_2pT4[2],"3^{rd} iteration","l")

canvasGenCost = TCanvas("canvasGenCost","canvasGenCost",200,10,1400,600)
canvasGenCost.Divide(3,1)

canvasGenCost.cd(1)
histGridCost_2pT4 = TH2D("histGridCost_2pT4","2 < #it{p}_{T} < 4 GeV/#it{c}",100,-1.,1.,100,0.04,0.06)
histGridCost_2pT4.Draw()
funcCost[0].Draw("same")
for i in range(nIterations):
    funcCost_2pT4[i].Draw("same")
    legendGenCost.Draw("same")

canvasGenCost.cd(2)
histGridCost_4pT6 = TH2D("histGridCost_4pT6","4 < #it{p}_{T} < 6 GeV/#it{c}",100,-1.,1.,100,0.04,0.06)
histGridCost_4pT6.Draw()
funcCost[1].Draw("same")
for i in range(nIterations):
    funcCost_4pT6[i].Draw("same")
    legendGenCost.Draw("same")

canvasGenCost.cd(3)
histGridCost_6pT10 = TH2D("histGridCost_6pT10","6 < #it{p}_{T} < 10 GeV/#it{c}",100,-1.,1.,100,0.04,0.06)
histGridCost_6pT10.Draw()
funcCost[2].Draw("same")
for i in range(nIterations):
    funcCost_6pT10[i].Draw("same")
    legendGenCost.Draw("same")

legendGenPhi = TLegend(0.3,0.1,0.7,0.3)
legendGenPhi.AddEntry(funcPhi[0],"0 iteration","l")
legendGenPhi.AddEntry(funcPhi_2pT4[0],"1^{st} iteration","l")
legendGenPhi.AddEntry(funcPhi_2pT4[1],"2^{nd} iteration","l")
legendGenPhi.AddEntry(funcPhi_2pT4[2],"3^{rd} iteration","l")

canvasGenPhi = TCanvas("canvasGenPhi","canvasGenPhi",200,10,1400,600)
canvasGenPhi.Divide(3,1)

canvasGenPhi.cd(1)
histGridPhi_2pT4 = TH2D("histGridPhi_2pT4","2 < #it{p}_{T} < 4 GeV/#it{c}",100,0.,PI,100,0.08,0.12)
histGridPhi_2pT4.Draw()
funcPhi[0].Draw("same")
for i in range(nIterations):
    funcPhi_2pT4[i].Draw("same")
    legendGenPhi.Draw("same")

canvasGenPhi.cd(2)
histGridPhi_4pT6 = TH2D("histGridPhi_4pT6","4 < #it{p}_{T} < 6 GeV/#it{c}",100,0.,PI,100,0.08,0.12)
histGridPhi_4pT6.Draw()
funcPhi[1].Draw("same")
for i in range(nIterations):
    funcPhi_4pT6[i].Draw("same")
    legendGenPhi.Draw("same")

canvasGenPhi.cd(3)
histGridPhi_6pT10 = TH2D("histGridPhi_6pT10","6 < #it{p}_{T} < 10 GeV/#it{c}",100,0.,PI,100,0.08,0.12)
histGridPhi_6pT10.Draw()
funcPhi[2].Draw("same")
for i in range(nIterations):
    funcPhi_6pT10[i].Draw("same")
    legendGenPhi.Draw("same")

################################################################################
# Results from the 2D approach
lambdaTheta2DSigmaFixedRange1 = [0.090071,0.042878,-0.187015]
errLambdaTheta2DSigmaFixedRange1 = [0.109741,0.120239,0.148915]
lambdaPhi2DSigmaFixedRange1 = [0.034499,0.002710,-0.023982]
errLambdaPhi2DSigmaFixedRange1 = [0.027053,0.030219,0.037875]
lambdaThetaPhi2DSigmaFixedRange1 = [-0.144842,-0.060507,-0.029900]
errLambdaThetaPhi2DSigmaFixedRange1 = [0.045375,0.051753,0.056858]

lambdaTheta2DSigmaFreeRange1 = [-0.012614,0.116041,-0.364637]
errLambdaTheta2DSigmaFreeRange1 = [0.144894,0.174576,0.159683]
lambdaPhi2DSigmaFreeRange1 = [-0.105262,-0.013791,-0.020319]
errLambdaPhi2DSigmaFreeRange1 = [0.038864,0.038414,0.042074]
lambdaThetaPhi2DSigmaFreeRange1 = [-0.393252,-0.090899,-0.023839]
errLambdaThetaPhi2DSigmaFreeRange1 = [0.059848,0.064584,0.063189]
#==============================

lambdaTheta2DSigmaFixedRange2 = [0.151902,0.114259,-0.211925]
errLambdaTheta2DSigmaFixedRange2 = [0.086331,0.098645,0.118930]
lambdaPhi2DSigmaFixedRange2 = [0.029227,0.009894,-0.020011]
errLambdaPhi2DSigmaFixedRange2 = [0.024582,0.029836,0.039160]
lambdaThetaPhi2DSigmaFixedRange2 = [-0.116712,-0.084138,-0.038373]
errLambdaThetaPhi2DSigmaFixedRange2 = [0.041448,0.045351,0.054891]

lambdaTheta2DSigmaFreeRange2 = [-0.228718,-0.018856,-0.355285]
errLambdaTheta2DSigmaFreeRange2 = [0.113497,0.122214,0.123781]
lambdaPhi2DSigmaFreeRange2 = [-0.130421,-0.005728,-0.015552]
errLambdaPhi2DSigmaFreeRange2 = [0.037736,0.037244,0.041935]
lambdaThetaPhi2DSigmaFreeRange2 = [-0.328887,-0.101709,-0.016528]
errLambdaThetaPhi2DSigmaFreeRange2 = [0.053661,0.057556,0.058968]
#===============================

lambdaTheta2DSigmaFixedRange3 = [-1.475312,0.112109,-0.102163]
errLambdaTheta2DSigmaFixedRange3 = [0.028088,0.091883,0.104001]
lambdaPhi2DSigmaFixedRange3 = [0.075501,0.009568,-0.021483]
errLambdaPhi2DSigmaFixedRange3 = [0.020608,0.029891,0.037627]
lambdaThetaPhi2DSigmaFixedRange3 = [-0.240544,-0.081499,-0.081782]
errLambdaThetaPhi2DSigmaFixedRange3 = [0.036157,0.044702,0.051858]

lambdaTheta2DSigmaFreeRange3 = [-0.342158,-0.136331,-0.407555]
errLambdaTheta2DSigmaFreeRange3 = [0.098063,0.102047,0.101571]
lambdaPhi2DSigmaFreeRange3 = [-0.117459,-0.004814,-0.023401]
errLambdaPhi2DSigmaFreeRange3 = [0.036713,0.036359,0.039777]
lambdaThetaPhi2DSigmaFreeRange3 = [-0.331377,-0.078762,-0.070618]
errLambdaThetaPhi2DSigmaFreeRange3 = [0.052152,0.055155,0.055898]

filePolarization2D = TFile.Open("/home/luca/GITHUB/polarization_classes/test_macros/output/PolarizationParameters2D.root")
histLambdaTheta2DSigmaFixed1 = filePolarization2D.Get("histLambdaThetaSigmaFixed1"); histLambdaTheta2DSigmaFixed1.SetDirectory(0);
histLambdaTheta2DSigmaFixed2 = filePolarization2D.Get("histLambdaThetaSigmaFixed2"); histLambdaTheta2DSigmaFixed2.SetDirectory(0);
histLambdaTheta2DSigmaFixed3 = filePolarization2D.Get("histLambdaThetaSigmaFixed3"); histLambdaTheta2DSigmaFixed3.SetDirectory(0);
histLambdaTheta2DSigmaFree1 = filePolarization2D.Get("histLambdaThetaSigmaFree1"); histLambdaTheta2DSigmaFree1.SetDirectory(0);
histLambdaTheta2DSigmaFree2 = filePolarization2D.Get("histLambdaThetaSigmaFree2"); histLambdaTheta2DSigmaFree2.SetDirectory(0);
histLambdaTheta2DSigmaFree3 = filePolarization2D.Get("histLambdaThetaSigmaFree3"); histLambdaTheta2DSigmaFree3.SetDirectory(0);
histLambdaPhi2DSigmaFixed1 = filePolarization2D.Get("histLambdaPhiSigmaFixed1"); histLambdaPhi2DSigmaFixed1.SetDirectory(0);
histLambdaPhi2DSigmaFixed2 = filePolarization2D.Get("histLambdaPhiSigmaFixed2"); histLambdaPhi2DSigmaFixed2.SetDirectory(0);
histLambdaPhi2DSigmaFixed3 = filePolarization2D.Get("histLambdaPhiSigmaFixed3"); histLambdaPhi2DSigmaFixed3.SetDirectory(0);
histLambdaPhi2DSigmaFree1 = filePolarization2D.Get("histLambdaPhiSigmaFree1"); histLambdaPhi2DSigmaFree1.SetDirectory(0);
histLambdaPhi2DSigmaFree2 = filePolarization2D.Get("histLambdaPhiSigmaFree2"); histLambdaPhi2DSigmaFree2.SetDirectory(0);
histLambdaPhi2DSigmaFree3 = filePolarization2D.Get("histLambdaPhiSigmaFree3"); histLambdaPhi2DSigmaFree3.SetDirectory(0);
histLambdaThetaPhi2DSigmaFixed1 = filePolarization2D.Get("histLambdaThetaPhiSigmaFixed1"); histLambdaThetaPhi2DSigmaFixed1.SetDirectory(0);
histLambdaThetaPhi2DSigmaFixed2 = filePolarization2D.Get("histLambdaThetaPhiSigmaFixed2"); histLambdaThetaPhi2DSigmaFixed2.SetDirectory(0);
histLambdaThetaPhi2DSigmaFixed3 = filePolarization2D.Get("histLambdaThetaPhiSigmaFixed3"); histLambdaThetaPhi2DSigmaFixed3.SetDirectory(0);
histLambdaThetaPhi2DSigmaFree1 = filePolarization2D.Get("histLambdaThetaPhiSigmaFree1"); histLambdaThetaPhi2DSigmaFree1.SetDirectory(0);
histLambdaThetaPhi2DSigmaFree2 = filePolarization2D.Get("histLambdaThetaPhiSigmaFree2"); histLambdaThetaPhi2DSigmaFree2.SetDirectory(0);
histLambdaThetaPhi2DSigmaFree3 = filePolarization2D.Get("histLambdaThetaPhiSigmaFree3"); histLambdaThetaPhi2DSigmaFree3.SetDirectory(0);

histProjCosThetaLambdaThetaSigmaFixed1 = filePolarization2D.Get("histProjCosThetaLambdaThetaSigmaFixed1"); histProjCosThetaLambdaThetaSigmaFixed1.SetDirectory(0);
histProjCosThetaLambdaThetaSigmaFixed2 = filePolarization2D.Get("histProjCosThetaLambdaThetaSigmaFixed2"); histProjCosThetaLambdaThetaSigmaFixed2.SetDirectory(0);
histProjCosThetaLambdaThetaSigmaFixed3 = filePolarization2D.Get("histProjCosThetaLambdaThetaSigmaFixed3"); histProjCosThetaLambdaThetaSigmaFixed3.SetDirectory(0);
histProjCosThetaLambdaThetaSigmaFree1 = filePolarization2D.Get("histProjCosThetaLambdaThetaSigmaFree1"); histProjCosThetaLambdaThetaSigmaFree1.SetDirectory(0);
histProjCosThetaLambdaThetaSigmaFree2 = filePolarization2D.Get("histProjCosThetaLambdaThetaSigmaFree2"); histProjCosThetaLambdaThetaSigmaFree2.SetDirectory(0);
histProjCosThetaLambdaThetaSigmaFree3 = filePolarization2D.Get("histProjCosThetaLambdaThetaSigmaFree3"); histProjCosThetaLambdaThetaSigmaFree3.SetDirectory(0);
filePolarization2D.Close()

################################################################################
lambdaThetaFin = [-0.203346,-0.0211515,0.0694346]
errLambdaThetaFin = [0.0965229,0.121205,0.150435]

lineZero = TLine(0.,0.,15.,0.); lineZero.SetLineStyle(kDashed);

histLambdaThetaFin = TH1D("histLambdaThetaFin","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaThetaFin.SetLineColor(kBlack); histLambdaThetaFin.SetLineWidth(2);
histLambdaTheta2DSigmaFixedRange1 = TH1D("histLambdaTheta2DSigmaFixedRange1","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaTheta2DSigmaFixedRange1.SetLineColor(kRed);
histLambdaTheta2DSigmaFixedRange2 = TH1D("histLambdaTheta2DSigmaFixedRange2","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaTheta2DSigmaFixedRange2.SetLineColor(kGreen);
histLambdaTheta2DSigmaFixedRange3 = TH1D("histLambdaTheta2DSigmaFixedRange3","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaTheta2DSigmaFixedRange3.SetLineColor(kMagenta);
histLambdaTheta2DSigmaFreeRange1 = TH1D("histLambdaTheta2DSigmaFreeRange1","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaTheta2DSigmaFreeRange1.SetLineColor(kRed+2);
histLambdaTheta2DSigmaFreeRange2 = TH1D("histLambdaTheta2DSigmaFreeRange2","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaTheta2DSigmaFreeRange2.SetLineColor(kGreen+2);
histLambdaTheta2DSigmaFreeRange3 = TH1D("histLambdaTheta2DSigmaFreeRange3","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaTheta2DSigmaFreeRange3.SetLineColor(kMagenta+2);
for i in range(1,4):
    histLambdaThetaFin.SetBinContent(i+1,lambdaThetaFin[i-1]); histLambdaThetaFin.SetBinError(i+1,errLambdaThetaFin[i-1]);
    histLambdaTheta2DSigmaFixedRange1.SetBinContent(i+1,lambdaTheta2DSigmaFixedRange1[i-1]); histLambdaTheta2DSigmaFixedRange1.SetBinError(i+1,errLambdaTheta2DSigmaFixedRange1[i-1]);
    histLambdaTheta2DSigmaFixedRange2.SetBinContent(i+1,lambdaTheta2DSigmaFixedRange2[i-1]); histLambdaTheta2DSigmaFixedRange2.SetBinError(i+1,errLambdaTheta2DSigmaFixedRange2[i-1]);
    histLambdaTheta2DSigmaFixedRange3.SetBinContent(i+1,lambdaTheta2DSigmaFixedRange3[i-1]); histLambdaTheta2DSigmaFixedRange3.SetBinError(i+1,errLambdaTheta2DSigmaFixedRange3[i-1]);
    histLambdaTheta2DSigmaFreeRange1.SetBinContent(i+1,lambdaTheta2DSigmaFreeRange1[i-1]); histLambdaTheta2DSigmaFreeRange1.SetBinError(i+1,errLambdaTheta2DSigmaFreeRange1[i-1]);
    histLambdaTheta2DSigmaFreeRange2.SetBinContent(i+1,lambdaTheta2DSigmaFreeRange2[i-1]); histLambdaTheta2DSigmaFreeRange2.SetBinError(i+1,errLambdaTheta2DSigmaFreeRange2[i-1]);
    histLambdaTheta2DSigmaFreeRange3.SetBinContent(i+1,lambdaTheta2DSigmaFreeRange3[i-1]); histLambdaTheta2DSigmaFreeRange3.SetBinError(i+1,errLambdaTheta2DSigmaFreeRange3[i-1]);

histGridLambdaTheta = TH2D("histGridLambdaTheta","#lambda_{#it{#theta}} vs #it{p}_{T}",100,0,15,100,-0.5,0.5)
histGridLambdaTheta.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
histGridLambdaTheta.GetYaxis().SetTitle("#lambda_{#it{#theta}}"); histGridLambdaTheta.GetYaxis().SetTitleOffset(1.2);
canvasLambdaThetaFin = TCanvas("canvasLambdaThetaFin","canvasLambdaThetaFin",20,20,600,600)
histGridLambdaTheta.Draw()
lineZero.Draw("same")
histLambdaThetaFin.Draw("Esame")
histLambdaTheta2DSigmaFixed2.Draw("Esame")
histLambdaTheta2DSigmaFree2.Draw("Esame")
histProjCosThetaLambdaThetaSigmaFixed2.Draw("Esame")
histProjCosThetaLambdaThetaSigmaFree2.Draw("Esame")
#histLambdaTheta2DSigmaFixedRange1.Draw("Esame")
#histLambdaTheta2DSigmaFixedRange2.Draw("Esame")
#histLambdaTheta2DSigmaFixedRange3.Draw("Esame")
#histLambdaTheta2DSigmaFreeRange1.Draw("Esame")
#histLambdaTheta2DSigmaFreeRange2.Draw("Esame")
#histLambdaTheta2DSigmaFreeRange3.Draw("Esame")
legendLambdaTheta = TLegend(0.6,0.7,0.9,0.9)
legendLambdaTheta.SetTextSize(0.04)
legendLambdaTheta.AddEntry(histLambdaThetaFin,"#lambda_{#theta}^{1D}","l")
#legendLambdaTheta.AddEntry(histLambdaTheta2DSigmaFixedRange2,"#lambda_{#theta}^{2D} #sigma_{J/#psi}^{Fixed}","l")
#legendLambdaTheta.AddEntry(histLambdaTheta2DSigmaFreeRange2,"#lambda_{#theta}^{2D} #sigma_{J/#psi}^{Free}","l")
legendLambdaTheta.AddEntry(histLambdaTheta2DSigmaFixed2,"#lambda_{#theta}^{2D} #sigma_{J/#psi}^{Fixed}","lp")
legendLambdaTheta.AddEntry(histLambdaTheta2DSigmaFree2,"#lambda_{#theta}^{2D} #sigma_{J/#psi}^{Free}","lp")
legendLambdaTheta.AddEntry(histProjCosThetaLambdaThetaSigmaFixed2,"#lambda_{#theta}^{2D projection} #sigma_{J/#psi}^{Fixed}","lp")
legendLambdaTheta.AddEntry(histProjCosThetaLambdaThetaSigmaFree2,"#lambda_{#theta}^{2D projection} #sigma_{J/#psi}^{Free}","lp")
legendLambdaTheta.Draw("same")

lambdaPhiFin = [-0.270491,-0.00887629,-0.0651542]
errLambdaPhiFin = [0.0483929,0.0541244,0.0678619]

histLambdaPhiFin = TH1D("histLambdaPhiFin","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaPhiFin.SetLineColor(kBlack); histLambdaPhiFin.SetLineWidth(2);
histLambdaPhi2DSigmaFixedRange1 = TH1D("histLambdaPhi2DSigmaFixedRange1","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaPhi2DSigmaFixedRange1.SetLineColor(kRed);
histLambdaPhi2DSigmaFixedRange2 = TH1D("histLambdaPhi2DSigmaFixedRange2","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaPhi2DSigmaFixedRange2.SetLineColor(kGreen);
histLambdaPhi2DSigmaFixedRange3 = TH1D("histLambdaPhi2DSigmaFixedRange3","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaPhi2DSigmaFixedRange3.SetLineColor(kMagenta);
histLambdaPhi2DSigmaFreeRange1 = TH1D("histLambdaPhi2DSigmaFreeRange1","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaPhi2DSigmaFreeRange1.SetLineColor(kRed+2);
histLambdaPhi2DSigmaFreeRange2 = TH1D("histLambdaPhi2DSigmaFreeRange2","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaPhi2DSigmaFreeRange2.SetLineColor(kGreen+2);
histLambdaPhi2DSigmaFreeRange3 = TH1D("histLambdaPhi2DSigmaFreeRange3","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaPhi2DSigmaFreeRange3.SetLineColor(kMagenta+2);
for i in range(1,4):
    histLambdaPhiFin.SetBinContent(i+1,lambdaPhiFin[i-1])
    histLambdaPhiFin.SetBinError(i+1,errLambdaPhiFin[i-1])
    histLambdaPhi2DSigmaFixedRange1.SetBinContent(i+1,lambdaPhi2DSigmaFixedRange1[i-1]); histLambdaPhi2DSigmaFixedRange1.SetBinError(i+1,errLambdaPhi2DSigmaFixedRange1[i-1]);
    histLambdaPhi2DSigmaFixedRange2.SetBinContent(i+1,lambdaPhi2DSigmaFixedRange2[i-1]); histLambdaPhi2DSigmaFixedRange2.SetBinError(i+1,errLambdaPhi2DSigmaFixedRange2[i-1]);
    histLambdaPhi2DSigmaFixedRange3.SetBinContent(i+1,lambdaPhi2DSigmaFixedRange3[i-1]); histLambdaPhi2DSigmaFixedRange3.SetBinError(i+1,errLambdaPhi2DSigmaFixedRange3[i-1]);
    histLambdaPhi2DSigmaFreeRange1.SetBinContent(i+1,lambdaPhi2DSigmaFreeRange1[i-1]); histLambdaPhi2DSigmaFreeRange1.SetBinError(i+1,errLambdaPhi2DSigmaFreeRange1[i-1]);
    histLambdaPhi2DSigmaFreeRange2.SetBinContent(i+1,lambdaPhi2DSigmaFreeRange2[i-1]); histLambdaPhi2DSigmaFreeRange2.SetBinError(i+1,errLambdaPhi2DSigmaFreeRange2[i-1]);
    histLambdaPhi2DSigmaFreeRange3.SetBinContent(i+1,lambdaPhi2DSigmaFreeRange3[i-1]); histLambdaPhi2DSigmaFreeRange3.SetBinError(i+1,errLambdaPhi2DSigmaFreeRange3[i-1]);

histGridLambdaPhi = TH2D("histGridLambdaPhi","#lambda_{#it{#varphi}} vs #it{p}_{T}",100,0,15,100,-0.5,0.5)
histGridLambdaPhi.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
histGridLambdaPhi.GetYaxis().SetTitle("#lambda_{#it{#varphi}}"); histGridLambdaPhi.GetYaxis().SetTitleOffset(1.2);
canvasLambdaPhiFin = TCanvas("canvasLambdaPhiFin","canvasLambdaPhiFin",20,20,600,600)
histGridLambdaPhi.Draw()
lineZero.Draw("same")
histLambdaPhiFin.Draw("Esame")
histLambdaPhi2DSigmaFixed2.Draw("Esame")
histLambdaPhi2DSigmaFree2.Draw("Esame")
#histLambdaPhi2DSigmaFixedRange1.Draw("Esame")
#histLambdaPhi2DSigmaFixedRange2.Draw("Esame")
#histLambdaPhi2DSigmaFixedRange3.Draw("Esame")
#histLambdaPhi2DSigmaFreeRange1.Draw("Esame")
#histLambdaPhi2DSigmaFreeRange2.Draw("Esame")
#histLambdaPhi2DSigmaFreeRange3.Draw("Esame")
legendLambdaPhi = TLegend(0.6,0.7,0.9,0.9)
legendLambdaPhi.SetTextSize(0.04)
legendLambdaPhi.AddEntry(histLambdaPhiFin,"#lambda_{#phi}^{1D}","l")
#legendLambdaPhi.AddEntry(histLambdaPhi2DSigmaFixedRange2,"#lambda_{#phi}^{2D} #sigma_{J/#psi}^{Fixed}","l")
#legendLambdaPhi.AddEntry(histLambdaPhi2DSigmaFreeRange2,"#lambda_{#phi}^{2D} #sigma_{J/#psi}^{Free}","l")
legendLambdaPhi.AddEntry(histLambdaTheta2DSigmaFixed2,"#lambda_{#theta}^{2D} #sigma_{J/#psi}^{Fixed}","lp")
legendLambdaPhi.AddEntry(histLambdaTheta2DSigmaFree2,"#lambda_{#theta}^{2D} #sigma_{J/#psi}^{Free}","lp")
legendLambdaPhi.Draw("same")


histLambdaThetaPhi2DSigmaFixedRange1 = TH1D("histLambdaThetaPhi2DSigmaFixedRange1","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaThetaPhi2DSigmaFixedRange1.SetLineColor(kRed);
histLambdaThetaPhi2DSigmaFixedRange2 = TH1D("histLambdaThetaPhi2DSigmaFixedRange2","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaThetaPhi2DSigmaFixedRange2.SetLineColor(kGreen);
histLambdaThetaPhi2DSigmaFixedRange3 = TH1D("histLambdaThetaPhi2DSigmaFixedRange3","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaThetaPhi2DSigmaFixedRange3.SetLineColor(kMagenta);
histLambdaThetaPhi2DSigmaFreeRange1 = TH1D("histLambdaThetaPhi2DSigmaFreeRange1","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaThetaPhi2DSigmaFreeRange1.SetLineColor(kRed+2);
histLambdaThetaPhi2DSigmaFreeRange2 = TH1D("histLambdaThetaPhi2DSigmaFreeRange2","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaThetaPhi2DSigmaFreeRange2.SetLineColor(kGreen+2);
histLambdaThetaPhi2DSigmaFreeRange3 = TH1D("histLambdaThetaPhi2DSigmaFreeRange3","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaThetaPhi2DSigmaFreeRange3.SetLineColor(kMagenta+2);
for i in range(1,4):
    histLambdaThetaPhi2DSigmaFixedRange1.SetBinContent(i+1,lambdaThetaPhi2DSigmaFixedRange1[i-1]); histLambdaThetaPhi2DSigmaFixedRange1.SetBinError(i+1,errLambdaThetaPhi2DSigmaFixedRange1[i-1]);
    histLambdaThetaPhi2DSigmaFixedRange2.SetBinContent(i+1,lambdaThetaPhi2DSigmaFixedRange2[i-1]); histLambdaThetaPhi2DSigmaFixedRange2.SetBinError(i+1,errLambdaThetaPhi2DSigmaFixedRange2[i-1]);
    histLambdaThetaPhi2DSigmaFixedRange3.SetBinContent(i+1,lambdaThetaPhi2DSigmaFixedRange3[i-1]); histLambdaThetaPhi2DSigmaFixedRange3.SetBinError(i+1,errLambdaThetaPhi2DSigmaFixedRange3[i-1]);
    histLambdaThetaPhi2DSigmaFreeRange1.SetBinContent(i+1,lambdaThetaPhi2DSigmaFreeRange1[i-1]); histLambdaThetaPhi2DSigmaFreeRange1.SetBinError(i+1,errLambdaThetaPhi2DSigmaFreeRange1[i-1]);
    histLambdaThetaPhi2DSigmaFreeRange2.SetBinContent(i+1,lambdaThetaPhi2DSigmaFreeRange2[i-1]); histLambdaThetaPhi2DSigmaFreeRange2.SetBinError(i+1,errLambdaThetaPhi2DSigmaFreeRange2[i-1]);
    histLambdaThetaPhi2DSigmaFreeRange3.SetBinContent(i+1,lambdaThetaPhi2DSigmaFreeRange3[i-1]); histLambdaThetaPhi2DSigmaFreeRange3.SetBinError(i+1,errLambdaThetaPhi2DSigmaFreeRange3[i-1]);

histGridLambdaThetaPhi = TH2D("histGridLambdaThetaPhi","#lambda_{#it{#theta#varphi}} vs #it{p}_{T}",100,0,15,100,-0.5,0.5)
histGridLambdaThetaPhi.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
histGridLambdaThetaPhi.GetYaxis().SetTitle("#lambda_{#it{#theta#varphi}}"); histGridLambdaThetaPhi.GetYaxis().SetTitleOffset(1.2);
canvasLambdaThetaPhiFin = TCanvas("canvasLambdaThetaPhiFin","canvasLambdaThetaPhiFin",20,20,600,600)
histGridLambdaThetaPhi.Draw()
lineZero.Draw("same")
histLambdaThetaPhi2DSigmaFixed2.Draw("Esame")
histLambdaThetaPhi2DSigmaFree2.Draw("Esame")
#histLambdaThetaPhi2DSigmaFixedRange1.Draw("Esame")
#histLambdaThetaPhi2DSigmaFixedRange2.Draw("Esame")
#histLambdaThetaPhi2DSigmaFixedRange3.Draw("Esame")
#histLambdaThetaPhi2DSigmaFreeRange1.Draw("Esame")
#histLambdaThetaPhi2DSigmaFreeRange2.Draw("Esame")
#histLambdaThetaPhi2DSigmaFreeRange3.Draw("Esame")
legendLambdaThetaPhi = TLegend(0.6,0.7,0.9,0.9)
legendLambdaThetaPhi.SetTextSize(0.05)
#legendLambdaThetaPhi.AddEntry(histLambdaThetaPhi2DSigmaFixedRange2,"#lambda_{#theta#varphi}^{2D} #sigma_{J/#psi}^{Fixed}","l")
#legendLambdaThetaPhi.AddEntry(histLambdaThetaPhi2DSigmaFreeRange2,"#lambda_{#theta#varphi}^{2D} #sigma_{J/#psi}^{Free}","l")
legendLambdaThetaPhi.AddEntry(histLambdaThetaPhi2DSigmaFixed2,"#lambda_{#theta#varphi}^{2D} #sigma_{J/#psi}^{Fixed}","lp")
legendLambdaThetaPhi.AddEntry(histLambdaThetaPhi2DSigmaFree2,"#lambda_{#theta#varphi}^{2D} #sigma_{J/#psi}^{Free}","lp")
legendLambdaThetaPhi.Draw("same")

raw_input()
