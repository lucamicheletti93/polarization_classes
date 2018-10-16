from ROOT import *
from array import array
import os.path
import math

gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("0.2g");
gROOT.ProcessLineSync(".L /home/luca/GITHUB/polarization_classes/MathFuncsLib.cxx+")
gROOT.ProcessLineSync(".x /home/luca/GITHUB/polarization_classes/Binning.cxx+")
gROOT.ProcessLineSync(".x /home/luca/GITHUB/polarization_classes/SpecialFitCalculator.cxx+")

# inizializing main quantities
PI = math.pi
fileBinning = TFile.Open("/home/luca/GITHUB/polarization_classes/test_macros/output/binning.root")
binning = fileBinning.Get("Binning")
CostValues = binning.GetCostValues()
CostWidth = binning.GetCostWidth()
NCostBins = CostValues.size() - 1
PhiValues = binning.GetPhiValues()
PhiWidth = binning.GetPhiWidth()
NPhiBins = PhiValues.size() - 1

CellAreaMatrix = binning.GetCellAreaMatrix();

nPtRanges = 3
titleHistPtRanges = ['2 < #it{p}_{T} < 4 GeV/#it{c}','4 < #it{p}_{T} < 6 GeV/#it{c}','6 < #it{p}_{T} < 10 GeV/#it{c}']
namePtRanges = ['2pt4','4pt6','6pt10']
namePtRangesNew = ['2pT4','4pT6','6pT10']

#minFitRangeCost = [-0.8,-0.7,-0.6]
#maxFitRangeCost = [0.8,0.7,0.6]
minFitRangeCost = -0.6
maxFitRangeCost = 0.6
minFitRangePhi = 0.502655
maxFitRangePhi = 2.63894

################
# 2D approach  #
################
histNJpsiSigmaFixed = []; histNJpsiSigmaFixedCorr = [];
histNJpsiSigmaFree = []; histNJpsiSigmaFreeCorr = [];
histAccxEff = []

histLambdaThetaSigmaFixed = TH1D("histLambdaThetaSigmaFixed","",4,array('d',[0.,2.,4.,6.,10.]))
histLambdaThetaSigmaFree = TH1D("histLambdaThetaSigmaFree","",4,array('d',[0.,2.,4.,6.,10.]))
histLambdaPhiSigmaFixed = TH1D("histLambdaPhiSigmaFixed","",4,array('d',[0.,2.,4.,6.,10.]))
histLambdaPhiSigmaFree = TH1D("histLambdaPhiSigmaFree","",4,array('d',[0.,2.,4.,6.,10.]))
histLambdaThetaPhiSigmaFixed = TH1D("histLambdaThetaPhiSigmaFixed","",4,array('d',[0.,2.,4.,6.,10.]))
histLambdaThetaPhiSigmaFree = TH1D("histLambdaThetaPhiSigmaFree","",4,array('d',[0.,2.,4.,6.,10.]))

for i in range(nPtRanges):
    fileNJpsiSigmaFixed = TFile.Open("/home/luca/GITHUB/polarization/2D_approach/signal_extraction/" + namePtRanges[i] + "_fixed_sigma_test/N_Jpsi_" + namePtRanges[i] + "_test.root")
    histNJpsiSigmaFixed.append(fileNJpsiSigmaFixed.Get("histNJpsi")); histNJpsiSigmaFixed[i].SetDirectory(0);
    fileNJpsiSigmaFree = TFile.Open("/home/luca/GITHUB/polarization/2D_approach/signal_extraction/" + namePtRanges[i] + "_free_sigma_test/N_Jpsi_" + namePtRanges[i] + "_test.root")
    histNJpsiSigmaFree.append(fileNJpsiSigmaFree.Get("histNJpsi")); histNJpsiSigmaFixed[i].SetDirectory(0);
    fileAccxEff = TFile.Open("/home/luca/GITHUB/polarization/2D_approach/accxeff/" + namePtRanges[i] + "_test/accxeff_" + namePtRanges[i] + "_test.root")
    histAccxEff.append(fileAccxEff.Get("histCostPhiAccxEffRebin")); histAccxEff[i].SetDirectory(0);

    for j in range(NCostBins):
        for k in range(NPhiBins):
            histNJpsiSigmaFixed[i].SetBinContent(j+1,k+1,(histNJpsiSigmaFixed[i].GetBinContent(j+1,k+1))/CellAreaMatrix[j][k])
            histNJpsiSigmaFixed[i].SetBinError(j+1,k+1,(histNJpsiSigmaFixed[i].GetBinError(j+1,k+1))/CellAreaMatrix[j][k])
            histNJpsiSigmaFree[i].SetBinContent(j+1,k+1,(histNJpsiSigmaFree[i].GetBinContent(j+1,k+1))/CellAreaMatrix[j][k])
            histNJpsiSigmaFree[i].SetBinError(j+1,k+1,(histNJpsiSigmaFree[i].GetBinError(j+1,k+1))/CellAreaMatrix[j][k])

    # Deleting bad fit points
    if namePtRanges[i] == '2pt4':
        histNJpsiSigmaFixed[i].SetBinContent(2,2,0.); histNJpsiSigmaFixed[i].SetBinError(2,2,0.);
        histNJpsiSigmaFree[i].SetBinContent(2,2,0.); histNJpsiSigmaFree[i].SetBinError(2,2,0.);

        histNJpsiSigmaFixed[i].SetBinContent(2,8,0.); histNJpsiSigmaFixed[i].SetBinError(2,8,0.);
        histNJpsiSigmaFixed[i].SetBinContent(2,9,0.); histNJpsiSigmaFixed[i].SetBinError(2,9,0.);
        histNJpsiSigmaFixed[i].SetBinContent(18,5,0.); histNJpsiSigmaFixed[i].SetBinError(18,5,0.);
        histNJpsiSigmaFixed[i].SetBinContent(18,9,0.); histNJpsiSigmaFixed[i].SetBinError(18,9,0.);

        histNJpsiSigmaFree[i].SetBinContent(2,9,0.); histNJpsiSigmaFree[i].SetBinError(2,9,0.);
        histNJpsiSigmaFree[i].SetBinContent(18,9,0.); histNJpsiSigmaFree[i].SetBinError(18,9,0.);

    if namePtRanges[i] == '6pt10':
        histNJpsiSigmaFixed[i].SetBinContent(2,3,0.); histNJpsiSigmaFixed[i].SetBinError(2,3,0.);
        histNJpsiSigmaFree[i].SetBinContent(2,3,0.); histNJpsiSigmaFree[i].SetBinError(2,3,0.);

        histNJpsiSigmaFixed[i].SetBinContent(18,9,0.); histNJpsiSigmaFixed[i].SetBinError(18,9,0.);
        histNJpsiSigmaFree[i].SetBinContent(18,9,0.); histNJpsiSigmaFree[i].SetBinError(18,9,0.);

    histNJpsiSigmaFixedCorr.append(histNJpsiSigmaFixed[i].Clone("histNJpsiSigmaFixedCorr")); histNJpsiSigmaFixedCorr[i].SetTitle(titleHistPtRanges[i] + " #sigma_{J/#psi}^{Fixed}"); histNJpsiSigmaFixedCorr[i].SetDirectory(0);
    histNJpsiSigmaFixedCorr[i].Divide(histAccxEff[i])

    funcPolSigmaFixed = TF2("funcPolSigmaFixed",MathFuncs.MyMathFuncs.MyFuncPol,-0.6,0.6,minFitRangePhi,maxFitRangePhi,4)
    funcPolSigmaFixed.SetParameter(0,1000.)
    funcPolSigmaFixed.SetParameter(1,0.)
    funcPolSigmaFixed.SetParameter(2,0.)
    funcPolSigmaFixed.FixParameter(3,0.)
    histNJpsiSigmaFixedCorr[i].Fit(funcPolSigmaFixed,"RSI0")
    histLambdaThetaSigmaFixed.SetBinContent(i+2,funcPolSigmaFixed.GetParameter(1)); histLambdaThetaSigmaFixed.SetBinError(i+2,funcPolSigmaFixed.GetParError(1));
    histLambdaPhiSigmaFixed.SetBinContent(i+2,funcPolSigmaFixed.GetParameter(2)); histLambdaPhiSigmaFixed.SetBinError(i+2,funcPolSigmaFixed.GetParError(2));
    histLambdaThetaPhiSigmaFixed.SetBinContent(i+2,funcPolSigmaFixed.GetParameter(3)); histLambdaThetaPhiSigmaFixed.SetBinError(i+2,funcPolSigmaFixed.GetParError(3));

    histNJpsiSigmaFreeCorr.append(histNJpsiSigmaFree[i].Clone("histNJpsiSigmaFreeCorr")); histNJpsiSigmaFreeCorr[i].SetTitle(titleHistPtRanges[i] + " #sigma_{J/#psi}^{Free}"); histNJpsiSigmaFreeCorr[i].SetDirectory(0);
    histNJpsiSigmaFreeCorr[i].Divide(histAccxEff[i])

    funcPolSigmaFree = TF2("funcPolSigmaFree",MathFuncs.MyMathFuncs.MyFuncPol,-0.6,0.6,minFitRangePhi,maxFitRangePhi,4)
    funcPolSigmaFree.SetParameter(0,1000.)
    funcPolSigmaFree.SetParameter(1,0.)
    funcPolSigmaFree.SetParameter(2,0.)
    funcPolSigmaFree.FixParameter(3,0.)
    histNJpsiSigmaFreeCorr[i].Fit(funcPolSigmaFree,"RSI0")
    histLambdaThetaSigmaFree.SetBinContent(i+2,funcPolSigmaFree.GetParameter(1)); histLambdaThetaSigmaFree.SetBinError(i+2,funcPolSigmaFree.GetParError(1));
    histLambdaPhiSigmaFree.SetBinContent(i+2,funcPolSigmaFree.GetParameter(2)); histLambdaPhiSigmaFree.SetBinError(i+2,funcPolSigmaFree.GetParError(2));
    histLambdaThetaPhiSigmaFree.SetBinContent(i+2,funcPolSigmaFree.GetParameter(3)); histLambdaThetaPhiSigmaFree.SetBinError(i+2,funcPolSigmaFree.GetParError(3));

    fileNJpsiSigmaFixed.Close()
    fileNJpsiSigmaFree.Close()
    fileAccxEff.Close()

# Control plots
canvasNJpsiCorr = TCanvas("canvasNJpsiCorr","canvasNJpsiCorr",20,20,1400,1200)
canvasNJpsiCorr.Divide(3,2)
canvasNJpsiCorr.cd(1); histNJpsiSigmaFixedCorr[0].Draw("COLZ");
canvasNJpsiCorr.cd(2); histNJpsiSigmaFixedCorr[1].Draw("COLZ");
canvasNJpsiCorr.cd(3); histNJpsiSigmaFixedCorr[2].Draw("COLZ");
canvasNJpsiCorr.cd(4); histNJpsiSigmaFreeCorr[0].Draw("COLZ");
canvasNJpsiCorr.cd(5); histNJpsiSigmaFreeCorr[1].Draw("COLZ");
canvasNJpsiCorr.cd(6); histNJpsiSigmaFreeCorr[2].Draw("COLZ");

# Projecting the 2D histogram in 1 dimension
histProjCosThetaNJpsiSigmaFixedCorr = []; histProjPhiNJpsiSigmaFixedCorr = [];
histProjCosThetaNJpsiSigmaFreeCorr = []; histProjPhiNJpsiSigmaFreeCorr = [];
histProjCosThetaNJpsiSigmaFixedCorrClone = []; histProjPhiNJpsiSigmaFixedCorrClone = [];
histProjCosThetaNJpsiSigmaFreeCorrClone = []; histProjPhiNJpsiSigmaFreeCorrClone = [];
histGridProjCosThetaNJpsi = []; histGridProjPhiNJpsi = [];

for i in range(nPtRanges):
    histProjCosThetaNJpsiSigmaFixedCorr.append(histNJpsiSigmaFixedCorr[i].ProjectionX("histProjCosThetaNJpsiSigmaFixedCorr" + namePtRangesNew[i])); histProjCosThetaNJpsiSigmaFixedCorr[i].SetLineColor(kRed); histProjCosThetaNJpsiSigmaFixedCorr[i].SetTitle("");
    histProjPhiNJpsiSigmaFixedCorr.append(histNJpsiSigmaFixedCorr[i].ProjectionY("histProjPhiNJpsiSigmaFixedCorr" + namePtRangesNew[i])); histProjPhiNJpsiSigmaFixedCorr[i].SetLineColor(kRed); histProjPhiNJpsiSigmaFixedCorr[i].SetTitle("");
    histProjCosThetaNJpsiSigmaFreeCorr.append(histNJpsiSigmaFreeCorr[i].ProjectionX("histProjCosThetaNJpsiSigmaFreeCorr" + namePtRangesNew[i])); histProjCosThetaNJpsiSigmaFreeCorr[i].SetLineColor(kBlue); histProjCosThetaNJpsiSigmaFreeCorr[i].SetTitle("");
    histProjPhiNJpsiSigmaFreeCorr.append(histNJpsiSigmaFreeCorr[i].ProjectionY("histProjPhiNJpsiSigmaFreeCorr" + namePtRangesNew[i])); histProjPhiNJpsiSigmaFreeCorr[i].SetLineColor(kBlue); histProjPhiNJpsiSigmaFreeCorr[i].SetTitle("");

    histProjCosThetaNJpsiSigmaFixedCorrClone.append(histProjCosThetaNJpsiSigmaFixedCorr[i].Clone()); histProjCosThetaNJpsiSigmaFixedCorrClone[i].Scale(1./(histProjCosThetaNJpsiSigmaFixedCorrClone[i].Integral()))
    histProjPhiNJpsiSigmaFixedCorrClone.append(histProjPhiNJpsiSigmaFixedCorr[i].Clone()); histProjPhiNJpsiSigmaFixedCorrClone[i].Scale(1./(histProjPhiNJpsiSigmaFixedCorrClone[i].Integral()))
    histProjCosThetaNJpsiSigmaFreeCorrClone.append(histProjCosThetaNJpsiSigmaFreeCorr[i].Clone()); histProjCosThetaNJpsiSigmaFreeCorrClone[i].Scale(1./(histProjCosThetaNJpsiSigmaFreeCorrClone[i].Integral()))
    histProjPhiNJpsiSigmaFreeCorrClone.append(histProjPhiNJpsiSigmaFreeCorr[i].Clone()); histProjPhiNJpsiSigmaFreeCorrClone[i].Scale(1./(histProjPhiNJpsiSigmaFreeCorrClone[i].Integral()))

    histGridProjCosThetaNJpsi.append(TH2D("cos#it{#theta}_{HE}" + titleHistPtRanges[i],"",100,-1.,1.,100,0.,0.1)); histGridProjCosThetaNJpsi[i].GetXaxis().SetTitle("cos#it{#theta}_{HE}"); histGridProjCosThetaNJpsi[i].SetTitle(titleHistPtRanges[i]);
    histGridProjPhiNJpsi.append(TH2D("#it{#varphi}_{HE}" + titleHistPtRanges[i],"",100,0.,PI,100,0.,0.2)); histGridProjPhiNJpsi[i].GetXaxis().SetTitle("#it{#varphi}_{HE}"); histGridProjPhiNJpsi[i].SetTitle(titleHistPtRanges[i]);

funcFit2DProjectionCosThetaSigmaFixed = []; funcFit2DProjectionPhiSigmaFixed = [];
funcFit2DProjectionCosThetaSigmaFree = []; funcFit2DProjectionPhiSigmaFree = [];

# Performing simultaneous fit
histLambdaThetaFit2DProjectionSigmaFixed = TH1D("histLambdaThetaFit2DProjectionSigmaFixed" + str(i+1),"",4,array('d',[0.,2.,4.,6.,10.])); histLambdaThetaFit2DProjectionSigmaFixed.SetLineColor(kRed);
histLambdaPhiFit2DProjectionSigmaFixed = TH1D("histLambdaPhiFit2DProjectionSigmaFixed" + str(i+1),"",4,array('d',[0.,2.,4.,6.,10.])); histLambdaPhiFit2DProjectionSigmaFixed.SetLineColor(kRed);

histLambdaThetaFit2DProjectionSigmaFree = TH1D("histLambdaThetaFit2DProjectionSigmaFree" + str(i+1),"",4,array('d',[0.,2.,4.,6.,10.])); histLambdaThetaFit2DProjectionSigmaFree.SetLineColor(kBlue);
histLambdaPhiFit2DProjectionSigmaFree = TH1D("histLambdaPhiFit2DProjectionSigmaFree" + str(i+1),"",4,array('d',[0.,2.,4.,6.,10.])); histLambdaPhiFit2DProjectionSigmaFree.SetLineColor(kBlue);

for i in range(nPtRanges):
    SimFit2DProjectionSigmaFixed = SpecialFitCalculator()
    SimFit2DProjectionSigmaFixed.SimultaneousFit2(histProjCosThetaNJpsiSigmaFixedCorrClone[i],histProjPhiNJpsiSigmaFixedCorrClone[i],minFitRangeCost,maxFitRangeCost,"")

    parameters = SimFit2DProjectionSigmaFixed.GetCosThetaParametersList(); errorParameters = SimFit2DProjectionSigmaFixed.GetErrorCosThetaParametersList();
    funcFit2DProjectionCosThetaSigmaFixed.append(TF1("funcFit2DProjectionCosThetaSigmaFixed_" + namePtRangesNew[i],MathFuncs.MyMathFuncs.MyFuncPolCosTheta,minFitRangeCost,maxFitRangeCost,2))
    funcFit2DProjectionCosThetaSigmaFixed[i].SetLineColor(kRed)
    funcFit2DProjectionCosThetaSigmaFixed[i].SetParameter(0,parameters[0]); funcFit2DProjectionCosThetaSigmaFixed[i].SetParameter(1,parameters[1]);
    histLambdaThetaFit2DProjectionSigmaFixed.SetBinContent(i+2,parameters[1]); histLambdaThetaFit2DProjectionSigmaFixed.SetBinError(i+2,errorParameters[1]);

    parameters = SimFit2DProjectionSigmaFixed.GetPhiParametersList(); errorParameters = SimFit2DProjectionSigmaFixed.GetErrorPhiParametersList();
    funcFit2DProjectionPhiSigmaFixed.append(TF1("funcFit2DProjectionPhiSigmaFixed_" + namePtRangesNew[i],MathFuncs.MyMathFuncs.MyFuncPolPhi,0.,PI,3))
    funcFit2DProjectionPhiSigmaFixed[i].SetLineColor(kRed)
    funcFit2DProjectionPhiSigmaFixed[i].SetParameter(0,parameters[0]); funcFit2DProjectionPhiSigmaFixed[i].SetParameter(1,parameters[1]); funcFit2DProjectionPhiSigmaFixed[i].SetParameter(2,parameters[2])
    histLambdaPhiFit2DProjectionSigmaFixed.SetBinContent(i+2,parameters[2]); histLambdaPhiFit2DProjectionSigmaFixed.SetBinError(i+2,errorParameters[2]);

    SimFit2DProjectionSigmaFree = SpecialFitCalculator()
    SimFit2DProjectionSigmaFree.SimultaneousFit2(histProjCosThetaNJpsiSigmaFreeCorrClone[i],histProjPhiNJpsiSigmaFreeCorrClone[i],minFitRangeCost,maxFitRangeCost,"")

    parameters = SimFit2DProjectionSigmaFree.GetCosThetaParametersList(); errorParameters = SimFit2DProjectionSigmaFree.GetErrorCosThetaParametersList();
    funcFit2DProjectionCosThetaSigmaFree.append(TF1("funcFit2DProjectionCosThetaSigmaFree_" + namePtRangesNew[i],MathFuncs.MyMathFuncs.MyFuncPolCosTheta,minFitRangeCost,maxFitRangeCost,2))
    funcFit2DProjectionCosThetaSigmaFree[i].SetLineColor(kBlue)
    funcFit2DProjectionCosThetaSigmaFree[i].SetParameter(0,parameters[0]); funcFit2DProjectionCosThetaSigmaFree[i].SetParameter(1,parameters[1]);
    histLambdaThetaFit2DProjectionSigmaFree.SetBinContent(i+2,parameters[1]); histLambdaThetaFit2DProjectionSigmaFree.SetBinError(i+2,errorParameters[1]);

    parameters = SimFit2DProjectionSigmaFree.GetPhiParametersList(); errorParameters = SimFit2DProjectionSigmaFree.GetErrorPhiParametersList();
    funcFit2DProjectionPhiSigmaFree.append(TF1("funcFit2DProjectionPhiSigmaFree_" + namePtRangesNew[i],MathFuncs.MyMathFuncs.MyFuncPolPhi,0.,PI,3))
    funcFit2DProjectionPhiSigmaFree[i].SetLineColor(kBlue)
    funcFit2DProjectionPhiSigmaFree[i].SetParameter(0,parameters[0]); funcFit2DProjectionPhiSigmaFree[i].SetParameter(1,parameters[1]); funcFit2DProjectionPhiSigmaFree[i].SetParameter(2,parameters[2])
    histLambdaPhiFit2DProjectionSigmaFree.SetBinContent(i+2,parameters[2]); histLambdaPhiFit2DProjectionSigmaFree.SetBinError(i+2,errorParameters[2]);

    del SimFit2DProjectionSigmaFixed
    del SimFit2DProjectionSigmaFree

# Control plots
canvasProjNJpsiCorr = TCanvas("canvasProjNJpsiCorr","canvasProjNJpsiCorr",20,20,1400,1200)
canvasProjNJpsiCorr.Divide(3,2)
canvasProjNJpsiCorr.cd(1);
histGridProjCosThetaNJpsi[0].Draw(); histProjCosThetaNJpsiSigmaFixedCorrClone[0].Draw("Esame"); histProjCosThetaNJpsiSigmaFreeCorrClone[0].Draw("Esame");
funcFit2DProjectionCosThetaSigmaFixed[0].Draw("Esame"); funcFit2DProjectionCosThetaSigmaFree[0].Draw("Esame");
canvasProjNJpsiCorr.cd(2);
histGridProjCosThetaNJpsi[1].Draw(); histProjCosThetaNJpsiSigmaFixedCorrClone[1].Draw("Esame"); histProjCosThetaNJpsiSigmaFreeCorrClone[1].Draw("Esame");
funcFit2DProjectionCosThetaSigmaFixed[1].Draw("Esame"); funcFit2DProjectionCosThetaSigmaFree[1].Draw("Esame");
canvasProjNJpsiCorr.cd(3);
histGridProjCosThetaNJpsi[2].Draw(); histProjCosThetaNJpsiSigmaFixedCorrClone[2].Draw("Esame"); histProjCosThetaNJpsiSigmaFreeCorrClone[2].Draw("Esame");
funcFit2DProjectionCosThetaSigmaFixed[2].Draw("Esame"); funcFit2DProjectionCosThetaSigmaFree[2].Draw("Esame");
canvasProjNJpsiCorr.cd(4);
histGridProjPhiNJpsi[0].Draw(); histProjPhiNJpsiSigmaFixedCorrClone[0].Draw("Esame"); histProjPhiNJpsiSigmaFreeCorrClone[0].Draw("Esame");
funcFit2DProjectionPhiSigmaFixed[0].Draw("Esame"); funcFit2DProjectionPhiSigmaFree[0].Draw("Esame");
canvasProjNJpsiCorr.cd(5);
histGridProjPhiNJpsi[1].Draw(); histProjPhiNJpsiSigmaFixedCorrClone[1].Draw("Esame"); histProjPhiNJpsiSigmaFreeCorrClone[1].Draw("Esame");
funcFit2DProjectionPhiSigmaFixed[1].Draw("Esame"); funcFit2DProjectionPhiSigmaFree[1].Draw("Esame");
canvasProjNJpsiCorr.cd(6);
histGridProjPhiNJpsi[2].Draw(); histProjPhiNJpsiSigmaFixedCorrClone[2].Draw("Esame"); histProjPhiNJpsiSigmaFreeCorrClone[2].Draw("Esame");
funcFit2DProjectionPhiSigmaFixed[2].Draw("Esame"); funcFit2DProjectionPhiSigmaFree[2].Draw("Esame");

################
# 1D approach  #
################
histNJpsiCost = []; histNJpsiCostCorr = [];
histNJpsiPhi = []; histNJpsiPhiCorr = [];
histAccxEffCost = [];
histAccxEffPhi = [];

funcFit1DCosTheta = []; funcFit1DPhi = [];

histLambdaThetaPure1D = TH1D("histLambdaThetaPure1D" + str(i+1),"",4,array('d',[0.,2.,4.,6.,10.])); histLambdaThetaPure1D.SetLineColor(kBlack); histLambdaThetaPure1D.SetLineWidth(2); histLambdaThetaPure1D.SetMarkerStyle(24); histLambdaThetaPure1D.SetMarkerColor(kBlack);
histLambdaPhiPure1D = TH1D("histLambdaPhiPure1D" + str(i+1),"",4,array('d',[0.,2.,4.,6.,10.])); histLambdaPhiPure1D.SetLineColor(kBlack); histLambdaPhiPure1D.SetLineWidth(2); histLambdaPhiPure1D.SetMarkerStyle(24); histLambdaPhiPure1D.SetMarkerColor(kBlack);

for i in range(nPtRanges):
    fileNJpsi1D = TFile.Open("/home/luca/GITHUB/polarization/1D_approach/signal_extraction/binned_1D_"  + namePtRanges[i] +  "_test/"   + namePtRanges[i] +   ".root")
    histNJpsiCost.append(fileNJpsi1D.Get("histNJpsiCost")); histNJpsiCost[i].SetDirectory(0);
    histNJpsiPhi.append(fileNJpsi1D.Get("histNJpsiPhi")); histNJpsiPhi[i].SetDirectory(0);
    fileAccxEff1D = TFile.Open("/home/luca/GITHUB/polarization_classes/test_macros/iterative_procedure/" + namePtRangesNew[i] + "/AccxEffReWeighted3rdStep.root")
    histAccxEffCost.append(fileAccxEff1D.Get("histAccxEffCostReWeighted_" + namePtRangesNew[i])); histAccxEffCost[i].SetDirectory(0);
    histAccxEffPhi.append(fileAccxEff1D.Get("histAccxEffPhiReWeighted_" + namePtRangesNew[i])); histAccxEffPhi[i].SetDirectory(0);

    histNJpsiCostCorr.append(histNJpsiCost[i].Clone("histNJpsiCostCorr_" + namePtRangesNew[i]))
    for j in range(19):
        histNJpsiCostCorr[i].SetBinContent(j+1,(histNJpsiCost[i].GetBinContent(j+1))/(histNJpsiCost[i].GetBinWidth(j+1)))
        histNJpsiCostCorr[i].SetBinError(j+1,(histNJpsiCost[i].GetBinError(j+1))/(histNJpsiCost[i].GetBinWidth(j+1)))
    histNJpsiCostCorr[i].Divide(histAccxEffCost[i]); histNJpsiCostCorr[i].Scale(1./(histNJpsiCostCorr[i].Integral())); histNJpsiCostCorr[i].SetLineColor(kBlack); histNJpsiCostCorr[i].SetDirectory(0);

    histNJpsiPhiCorr.append(histNJpsiPhi[i].Clone("histNJpsiPhiCorr_" + namePtRangesNew[i]))
    for j in range(10):
        histNJpsiPhiCorr[i].SetBinContent(j+1,(histNJpsiPhi[i].GetBinContent(j+1))/(histNJpsiPhi[i].GetBinWidth(j+1)));
        histNJpsiPhiCorr[i].SetBinError(j+1,(histNJpsiPhi[i].GetBinError(j+1))/(histNJpsiPhi[i].GetBinWidth(j+1)));
    histNJpsiPhiCorr[i].Divide(histAccxEffPhi[i]); histNJpsiPhiCorr[i].Scale(1./(histNJpsiPhiCorr[i].Integral())); histNJpsiPhiCorr[i].SetLineColor(kBlack); histNJpsiPhiCorr[i].SetDirectory(0);

    SimFitPure1D = SpecialFitCalculator()
    SimFitPure1D.SimultaneousFit2(histNJpsiCostCorr[i],histNJpsiPhiCorr[i],-0.8,0.8,"")

    parameters = SimFitPure1D.GetCosThetaParametersList(); errorParameters = SimFitPure1D.GetErrorCosThetaParametersList()
    funcFit1DCosTheta.append(TF1("funcFit1DCosTheta_" + namePtRangesNew[i],MathFuncs.MyMathFuncs.MyFuncPolCosTheta,minFitRangeCost,maxFitRangeCost,2))
    funcFit1DCosTheta[i].SetLineColor(kBlack); funcFit1DCosTheta[i].SetLineStyle(kDashed);
    funcFit1DCosTheta[i].SetParameter(0,parameters[0]); funcFit1DCosTheta[i].SetParameter(1,parameters[1]);
    histLambdaThetaPure1D.SetBinContent(i+2,parameters[1]); histLambdaThetaPure1D.SetBinError(i+2,errorParameters[1]);

    parameters = SimFitPure1D.GetPhiParametersList(); errorParameters = SimFitPure1D.GetErrorPhiParametersList()
    funcFit1DPhi.append(TF1("funcFit1DPhi_" + namePtRangesNew[i],MathFuncs.MyMathFuncs.MyFuncPolPhi,0.,PI,3))
    funcFit1DPhi[i].SetLineColor(kBlack); funcFit1DPhi[i].SetLineStyle(kDashed);
    funcFit1DPhi[i].SetParameter(0,parameters[0]); funcFit1DPhi[i].SetParameter(1,parameters[1]); funcFit1DPhi[i].SetParameter(2,parameters[2])
    histLambdaPhiPure1D.SetBinContent(i+2,parameters[2]); histLambdaPhiPure1D.SetBinError(i+2,errorParameters[2]);

    del SimFitPure1D
    fileNJpsi1D.Close()
    fileAccxEff.Close()

canvasPure1DNJpsiCorr = TCanvas("canvasPure1DNJpsiCorr","canvasPure1DNJpsiCorr",20,20,1400,1200)
canvasPure1DNJpsiCorr.Divide(3,2)
canvasPure1DNJpsiCorr.cd(1);
histGridProjCosThetaNJpsi[0].Draw(); histNJpsiCostCorr[0].Draw("Esame");
funcFit1DCosTheta[0].Draw("Esame")
canvasPure1DNJpsiCorr.cd(2);
histGridProjCosThetaNJpsi[1].Draw(); histNJpsiCostCorr[1].Draw("Esame");
funcFit1DCosTheta[1].Draw("Esame")
canvasPure1DNJpsiCorr.cd(3);
histGridProjCosThetaNJpsi[2].Draw(); histNJpsiCostCorr[2].Draw("Esame");
funcFit1DCosTheta[2].Draw("Esame")
canvasPure1DNJpsiCorr.cd(4);
histGridProjPhiNJpsi[0].Draw(); histNJpsiPhiCorr[0].Draw("Esame");
funcFit1DPhi[0].Draw("Esame")
canvasPure1DNJpsiCorr.cd(5);
histGridProjPhiNJpsi[1].Draw(); histNJpsiPhiCorr[1].Draw("Esame");
funcFit1DPhi[1].Draw("Esame")
canvasPure1DNJpsiCorr.cd(6);
histGridProjPhiNJpsi[2].Draw(); histNJpsiPhiCorr[2].Draw("Esame");
funcFit1DPhi[2].Draw("Esame")

################################################################################
# Polarization parameters comparison
lineZero = TLine(0.,0.,15.,0.); lineZero.SetLineStyle(kDashed);

histGridLambdaTheta = TH2D("histGridLambdaTheta","#lambda_{#it{#theta}} vs #it{p}_{T}",100,0,15,100,-0.5,0.5)
histGridLambdaTheta.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
histGridLambdaTheta.GetYaxis().SetTitle("#lambda_{#it{#theta}}"); histGridLambdaTheta.GetYaxis().SetTitleOffset(1.2);

canvasLambdaTheta = TCanvas("canvasLambdaTheta","canvasLambdaTheta",20,20,1200,600)
canvasLambdaTheta.Divide(2,1)

canvasLambdaTheta.cd(1)
histGridLambdaTheta.Draw()
lineZero.Draw("same")
histLambdaThetaPure1D.Draw("EPsame")
histLambdaThetaFit2DProjectionSigmaFixed.Draw("EPsame")
histLambdaThetaSigmaFixed.Draw("EPsame")
#for i in range(len(minFitRangeCost)):
    #histLambdaThetaSigmaFixed[i].Draw("Esame")
    #histProjLambdaThetaSigmaFixed[i].Draw("EPsame")
#legendSigmaFixed.Draw("same")

canvasLambdaTheta.cd(2)
histGridLambdaTheta.Draw()
lineZero.Draw("same")
histLambdaThetaPure1D.Draw("EPsame")
histLambdaThetaFit2DProjectionSigmaFree.Draw("EPsame")
histLambdaThetaSigmaFree.Draw("EPsame")
#for i in range(len(minFitRangeCost)):
    #histLambdaThetaSigmaFree[i].Draw("EPsame")
    #histProjLambdaThetaSigmaFree[i].Draw("EPsame")
#legendSigmaFree.Draw("same")

histGridLambdaPhi = TH2D("histGridLambdaPhi","#lambda_{#it{#varphi}} vs #it{p}_{T}",100,0,15,100,-0.5,0.5)
histGridLambdaPhi.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
histGridLambdaPhi.GetYaxis().SetTitle("#lambda_{#it{#varphi}}"); histGridLambdaPhi.GetYaxis().SetTitleOffset(1.2);

canvasLambdaPhi = TCanvas("canvasLambdaPhi","canvasLambdaPhi",20,20,1200,600)
canvasLambdaPhi.Divide(2,1)

canvasLambdaPhi.cd(1)
histGridLambdaPhi.Draw()
lineZero.Draw("same")
histLambdaPhiPure1D.Draw("EPsame")
histLambdaPhiFit2DProjectionSigmaFixed.Draw("EPsame")
histLambdaPhiSigmaFixed.Draw("EPsame")
#for i in range(len(minFitRangeCost)):
    #histLambdaPhiSigmaFixed[i].Draw("Esame")
    #histProjLambdaPhiSigmaFixed[i].Draw("EPsame")
#legendSigmaFixed.Draw("same")

canvasLambdaPhi.cd(2)
histGridLambdaPhi.Draw()
lineZero.Draw("same")
histLambdaPhiPure1D.Draw("EPsame")
histLambdaPhiFit2DProjectionSigmaFree.Draw("EPsame")
histLambdaPhiSigmaFree.Draw("EPsame")
#for i in range(len(minFitRangeCost)):
    #histLambdaPhiSigmaFree[i].Draw("EPsame")
    #histProjLambdaPhiSigmaFree[i].Draw("EPsame")
#legendSigmaFree.Draw("same")



raw_input()
