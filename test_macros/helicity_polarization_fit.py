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
#fileBinning = TFile.Open("/home/luca/GITHUB/polarization_classes/test_macros/output/binning2D.root")
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
minFitRangeCost = -0.7
maxFitRangeCost = 0.7
minFitRangePhi = 0.502655
maxFitRangePhi = 2.63894

################
# 2D approach  #
################
histNJpsiSigmaFixed = []; histNJpsiSigmaFixedCorr = [];
histNJpsiSigmaFree = []; histNJpsiSigmaFreeCorr = [];
histAccxEff = []
histAccxEffProjectionCosTheta = []; histAccxEffProjectionPhi = [];

funcPolSigmaFixed = []; funcPolSigmaFree = [];
funcPolProjectionCosThetaSigmaFixed = []; funcPolProjectionPhiSigmaFixed = [];                # Projection on CosTheta and Phi of the 2D fit

histLambdaThetaSigmaFixed = TH1D("histLambdaThetaSigmaFixed","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaThetaSigmaFixed.SetLineColor(kRed+1); histLambdaThetaSigmaFixed.SetLineStyle(kDashed); histLambdaThetaSigmaFixed.SetLineWidth(2);
histLambdaThetaSigmaFree = TH1D("histLambdaThetaSigmaFree","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaThetaSigmaFree.SetLineColor(kBlue+1); histLambdaThetaSigmaFree.SetLineStyle(kDashed); histLambdaThetaSigmaFree.SetLineWidth(2);
histLambdaPhiSigmaFixed = TH1D("histLambdaPhiSigmaFixed","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaPhiSigmaFixed.SetLineColor(kRed+1); histLambdaPhiSigmaFixed.SetLineStyle(kDashed); histLambdaPhiSigmaFixed.SetLineWidth(2);
histLambdaPhiSigmaFree = TH1D("histLambdaPhiSigmaFree","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaPhiSigmaFree.SetLineColor(kBlue+1); histLambdaPhiSigmaFree.SetLineStyle(kDashed); histLambdaPhiSigmaFree.SetLineWidth(2);
histLambdaThetaPhiSigmaFixed = TH1D("histLambdaThetaPhiSigmaFixed","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaThetaPhiSigmaFixed.SetLineColor(kRed+1); histLambdaThetaPhiSigmaFixed.SetLineStyle(kDashed); histLambdaThetaPhiSigmaFixed.SetLineWidth(2);
histLambdaThetaPhiSigmaFree = TH1D("histLambdaThetaPhiSigmaFree","",4,array('d',[0.,2.,4.,6.,10.])); histLambdaThetaPhiSigmaFree.SetLineColor(kBlue+1); histLambdaThetaPhiSigmaFree.SetLineStyle(kDashed); histLambdaThetaPhiSigmaFree.SetLineWidth(2);

for i in range(nPtRanges):
    fileNJpsiSigmaFixed = TFile.Open("/home/luca/GITHUB/polarization/2D_approach/signal_extraction/" + namePtRanges[i] + "_fixed_sigma_test/N_Jpsi_" + namePtRanges[i] + "_test.root")
    #fileNJpsiSigmaFixed = TFile.Open("/home/luca/GITHUB/polarization/2D_approach/signal_extraction/" + namePtRanges[i] + "_fixed_sigma_test2DBinning/N_Jpsi_" + namePtRanges[i] + "_test2DBinning.root")
    histNJpsiSigmaFixed.append(fileNJpsiSigmaFixed.Get("histNJpsi")); histNJpsiSigmaFixed[i].SetDirectory(0);
    fileNJpsiSigmaFree = TFile.Open("/home/luca/GITHUB/polarization/2D_approach/signal_extraction/" + namePtRanges[i] + "_free_sigma_test/N_Jpsi_" + namePtRanges[i] + "_test.root")
    #fileNJpsiSigmaFree = TFile.Open("/home/luca/GITHUB/polarization/2D_approach/signal_extraction/" + namePtRanges[i] + "_free_sigma_test2DBinning/N_Jpsi_" + namePtRanges[i] + "_test2DBinning.root")
    histNJpsiSigmaFree.append(fileNJpsiSigmaFree.Get("histNJpsi")); histNJpsiSigmaFixed[i].SetDirectory(0);
    fileAccxEff = TFile.Open("/home/luca/GITHUB/polarization/2D_approach/accxeff/" + namePtRanges[i] + "_test/accxeff_" + namePtRanges[i] + "_test.root")
    #fileAccxEff = TFile.Open("/home/luca/GITHUB/polarization/2D_approach/accxeff/" + namePtRanges[i] + "_test2DBinning/accxeff_" + namePtRanges[i] + "_test2DBinning.root")
    histAccxEff.append(fileAccxEff.Get("histCostPhiAccxEffRebin")); histAccxEff[i].SetDirectory(0);
    #histAccxEffProjectionCosTheta.append(histAccxEff[i].ProjectionX("histAccxEffProjectionCosTheta" + namePtRangesNew[i])); histAccxEffProjectionCosTheta[i].SetLineColor(kRed); histAccxEffProjectionCosTheta[i].SetDirectory(0);
    #histAccxEffProjectionCosTheta[i].Scale(1./(histAccxEffProjectionCosTheta[i].Integral()))
    #histAccxEffProjectionPhi.append(histAccxEff[i].ProjectionY("histAccxEffProjectionPhi" + namePtRangesNew[i])); histAccxEffProjectionPhi[i].SetLineColor(kRed); histAccxEffProjectionPhi[i].SetDirectory(0);
    #histAccxEffProjectionPhi[i].Scale(1./(histAccxEffProjectionPhi[i].Integral()))

    for j in range(NCostBins):
        for k in range(NPhiBins):
            histNJpsiSigmaFixed[i].SetBinContent(j+1,k+1,(histNJpsiSigmaFixed[i].GetBinContent(j+1,k+1))/CellAreaMatrix[j][k])
            histNJpsiSigmaFixed[i].SetBinError(j+1,k+1,(histNJpsiSigmaFixed[i].GetBinError(j+1,k+1))/CellAreaMatrix[j][k])
            histNJpsiSigmaFree[i].SetBinContent(j+1,k+1,(histNJpsiSigmaFree[i].GetBinContent(j+1,k+1))/CellAreaMatrix[j][k])
            histNJpsiSigmaFree[i].SetBinError(j+1,k+1,(histNJpsiSigmaFree[i].GetBinError(j+1,k+1))/CellAreaMatrix[j][k])

    # Deleting bad fit points
    if namePtRanges[i] == '2pt4':
        #histNJpsiSigmaFixed[i].SetBinContent(2,2,0.); histNJpsiSigmaFixed[i].SetBinError(2,2,0.);
        histNJpsiSigmaFree[i].SetBinContent(2,2,0.); histNJpsiSigmaFree[i].SetBinError(2,2,0.);

        histNJpsiSigmaFixed[i].SetBinContent(2,8,0.); histNJpsiSigmaFixed[i].SetBinError(2,8,0.);
        histNJpsiSigmaFixed[i].SetBinContent(2,9,0.); histNJpsiSigmaFixed[i].SetBinError(2,9,0.);
        #histNJpsiSigmaFixed[i].SetBinContent(18,5,0.); histNJpsiSigmaFixed[i].SetBinError(18,5,0.);
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

    funcPolSigmaFixed.append(TF2("funcPolSigmaFixed",MathFuncs.MyMathFuncs.MyFuncPol,minFitRangeCost,maxFitRangeCost,minFitRangePhi,maxFitRangePhi,4))
    funcPolSigmaFixed[i].SetParameter(0,1000.)
    funcPolSigmaFixed[i].SetParameter(1,0.)
    funcPolSigmaFixed[i].SetParameter(2,0.)
    funcPolSigmaFixed[i].FixParameter(3,0.)
    #histNJpsiSigmaFixedCorr[i].Fit(funcPolSigmaFixed[i],"RSI0")
    histNJpsiSigmaFixedCorr[i].Fit(funcPolSigmaFixed[i],"RSIWL0")               # Barbatrucco
    histLambdaThetaSigmaFixed.SetBinContent(i+2,funcPolSigmaFixed[i].GetParameter(1)); histLambdaThetaSigmaFixed.SetBinError(i+2,funcPolSigmaFixed[i].GetParError(1));
    histLambdaPhiSigmaFixed.SetBinContent(i+2,funcPolSigmaFixed[i].GetParameter(2)); histLambdaPhiSigmaFixed.SetBinError(i+2,funcPolSigmaFixed[i].GetParError(2));
    histLambdaThetaPhiSigmaFixed.SetBinContent(i+2,funcPolSigmaFixed[i].GetParameter(3)); histLambdaThetaPhiSigmaFixed.SetBinError(i+2,funcPolSigmaFixed[i].GetParError(3));

    funcPolProjectionCosThetaSigmaFixed.append(TF1("funcPolProjectionCosThetaSigmaFixed_" + namePtRangesNew[i],MathFuncs.MyMathFuncs.MyFuncPolCosTheta,minFitRangeCost,maxFitRangeCost,2))
    funcPolProjectionCosThetaSigmaFixed[i].FixParameter(1,funcPolSigmaFixed[i].GetParameter(1))

    funcPolProjectionPhiSigmaFixed.append(TF1("funcPolProjectionPhiSigmaFixed_" + namePtRangesNew[i],MathFuncs.MyMathFuncs.MyFuncPolPhi,0.,PI,3))
    funcPolProjectionPhiSigmaFixed[i].FixParameter(1,funcPolSigmaFixed[i].GetParameter(1))
    funcPolProjectionPhiSigmaFixed[i].FixParameter(2,funcPolSigmaFixed[i].GetParameter(2))

    histNJpsiSigmaFreeCorr.append(histNJpsiSigmaFree[i].Clone("histNJpsiSigmaFreeCorr")); histNJpsiSigmaFreeCorr[i].SetTitle(titleHistPtRanges[i] + " #sigma_{J/#psi}^{Free}"); histNJpsiSigmaFreeCorr[i].SetDirectory(0);
    histNJpsiSigmaFreeCorr[i].Divide(histAccxEff[i])

    #funcPolSigmaFree.append(TF2("funcPolSigmaFree",MathFuncs.MyMathFuncs.MyFuncPol,minFitRangeCost,maxFitRangeCost,minFitRangePhi,maxFitRangePhi,4))
    #funcPolSigmaFree[i].SetParameter(0,1000.)
    #funcPolSigmaFree[i].SetParameter(1,0.)
    #funcPolSigmaFree[i].SetParameter(2,0.)
    #funcPolSigmaFree[i].FixParameter(3,0.)
    #histNJpsiSigmaFreeCorr[i].Fit(funcPolSigmaFree[i],"RSI0")
    #histLambdaThetaSigmaFree.SetBinContent(i+2,funcPolSigmaFree[i].GetParameter(1)); histLambdaThetaSigmaFree.SetBinError(i+2,funcPolSigmaFree[i].GetParError(1));
    #histLambdaPhiSigmaFree.SetBinContent(i+2,funcPolSigmaFree[i].GetParameter(2)); histLambdaPhiSigmaFree.SetBinError(i+2,funcPolSigmaFree[i].GetParError(2));
    #histLambdaThetaPhiSigmaFree.SetBinContent(i+2,funcPolSigmaFree[i].GetParameter(3)); histLambdaThetaPhiSigmaFree.SetBinError(i+2,funcPolSigmaFree[i].GetParError(3));

    fileNJpsiSigmaFixed.Close()
    fileNJpsiSigmaFree.Close()
    fileAccxEff.Close()

# Control plots
canvasNJpsiCorr = TCanvas("canvasNJpsiCorr","canvasNJpsiCorr",20,20,1400,1200)
canvasNJpsiCorr.Divide(3,2)
canvasNJpsiCorr.cd(1); histNJpsiSigmaFixedCorr[0].Draw("COLZ"); funcPolSigmaFixed[0].Draw("same");
canvasNJpsiCorr.cd(2); histNJpsiSigmaFixedCorr[1].Draw("COLZ"); funcPolSigmaFixed[1].Draw("same");
canvasNJpsiCorr.cd(3); histNJpsiSigmaFixedCorr[2].Draw("COLZ"); funcPolSigmaFixed[2].Draw("same");
canvasNJpsiCorr.cd(4); histNJpsiSigmaFreeCorr[0].Draw("COLZ"); #funcPolSigmaFree[0].Draw("same");
canvasNJpsiCorr.cd(5); histNJpsiSigmaFreeCorr[1].Draw("COLZ"); #funcPolSigmaFree[1].Draw("same");
canvasNJpsiCorr.cd(6); histNJpsiSigmaFreeCorr[2].Draw("COLZ"); #funcPolSigmaFree[2].Draw("same");

# Projecting the 2D histogram in 1 dimension
histProjCosThetaNJpsiSigmaFixedCorr = []; histProjPhiNJpsiSigmaFixedCorr = [];
histProjCosThetaNJpsiSigmaFreeCorr = []; histProjPhiNJpsiSigmaFreeCorr = [];
histProjCosThetaNJpsiSigmaFixedCorrClone = []; histProjPhiNJpsiSigmaFixedCorrClone = [];
histProjCosThetaNJpsiSigmaFreeCorrClone = []; histProjPhiNJpsiSigmaFreeCorrClone = [];
histGridProjCosThetaNJpsi = []; histGridProjPhiNJpsi = [];

for i in range(nPtRanges):
    histProjCosThetaNJpsiSigmaFixedCorr.append(histNJpsiSigmaFixedCorr[i].ProjectionX("histProjCosThetaNJpsiSigmaFixedCorr" + namePtRangesNew[i])); histProjCosThetaNJpsiSigmaFixedCorr[i].SetLineColor(kRed); histProjCosThetaNJpsiSigmaFixedCorr[i].SetTitle("");
    #histProjPhiNJpsiSigmaFixedCorr.append(histNJpsiSigmaFixedCorr[i].ProjectionY("histProjPhiNJpsiSigmaFixedCorr" + namePtRangesNew[i])); histProjPhiNJpsiSigmaFixedCorr[i].SetLineColor(kRed); histProjPhiNJpsiSigmaFixedCorr[i].SetTitle("");
    histProjPhiNJpsiSigmaFixedCorr.append(histNJpsiSigmaFixedCorr[i].ProjectionY("histProjPhiNJpsiSigmaFixedCorr" + namePtRangesNew[i],3,17)); histProjPhiNJpsiSigmaFixedCorr[i].SetLineColor(kRed); histProjPhiNJpsiSigmaFixedCorr[i].SetTitle("");   #Barbatrucco
    histProjCosThetaNJpsiSigmaFreeCorr.append(histNJpsiSigmaFreeCorr[i].ProjectionX("histProjCosThetaNJpsiSigmaFreeCorr" + namePtRangesNew[i])); histProjCosThetaNJpsiSigmaFreeCorr[i].SetLineColor(kBlue); histProjCosThetaNJpsiSigmaFreeCorr[i].SetTitle("");
    #histProjPhiNJpsiSigmaFreeCorr.append(histNJpsiSigmaFreeCorr[i].ProjectionY("histProjPhiNJpsiSigmaFreeCorr" + namePtRangesNew[i])); histProjPhiNJpsiSigmaFreeCorr[i].SetLineColor(kBlue); histProjPhiNJpsiSigmaFreeCorr[i].SetTitle("");
    histProjPhiNJpsiSigmaFreeCorr.append(histNJpsiSigmaFreeCorr[i].ProjectionY("histProjPhiNJpsiSigmaFreeCorr" + namePtRangesNew[i],3,17)); histProjPhiNJpsiSigmaFreeCorr[i].SetLineColor(kBlue); histProjPhiNJpsiSigmaFreeCorr[i].SetTitle("");       #Barbatrucco

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

    histProjCosThetaNJpsiSigmaFixedCorrClone[i].Fit(funcPolProjectionCosThetaSigmaFixed[i],"RSI0")
    funcPolProjectionCosThetaSigmaFixed[i].SetLineColor(kRed); funcPolProjectionCosThetaSigmaFixed[i].SetLineStyle(kDashed);

    parameters = SimFit2DProjectionSigmaFixed.GetPhiParametersList(); errorParameters = SimFit2DProjectionSigmaFixed.GetErrorPhiParametersList();
    funcFit2DProjectionPhiSigmaFixed.append(TF1("funcFit2DProjectionPhiSigmaFixed_" + namePtRangesNew[i],MathFuncs.MyMathFuncs.MyFuncPolPhi,0.,PI,3))
    funcFit2DProjectionPhiSigmaFixed[i].SetLineColor(kRed)
    funcFit2DProjectionPhiSigmaFixed[i].SetParameter(0,parameters[0]); funcFit2DProjectionPhiSigmaFixed[i].SetParameter(1,parameters[1]); funcFit2DProjectionPhiSigmaFixed[i].SetParameter(2,parameters[2])
    histLambdaPhiFit2DProjectionSigmaFixed.SetBinContent(i+2,parameters[2]); histLambdaPhiFit2DProjectionSigmaFixed.SetBinError(i+2,errorParameters[2]);

    histProjCosThetaNJpsiSigmaFixedCorrClone[i].Fit(funcPolProjectionPhiSigmaFixed[i],"RSI0")
    funcPolProjectionPhiSigmaFixed[i].SetLineColor(kRed); funcPolProjectionPhiSigmaFixed[i].SetLineStyle(kDashed);

    #SimFit2DProjectionSigmaFree = SpecialFitCalculator()
    #SimFit2DProjectionSigmaFree.SimultaneousFit2(histProjCosThetaNJpsiSigmaFreeCorrClone[i],histProjPhiNJpsiSigmaFreeCorrClone[i],minFitRangeCost,maxFitRangeCost,"")

    #parameters = SimFit2DProjectionSigmaFree.GetCosThetaParametersList(); errorParameters = SimFit2DProjectionSigmaFree.GetErrorCosThetaParametersList();
    #funcFit2DProjectionCosThetaSigmaFree.append(TF1("funcFit2DProjectionCosThetaSigmaFree_" + namePtRangesNew[i],MathFuncs.MyMathFuncs.MyFuncPolCosTheta,minFitRangeCost,maxFitRangeCost,2))
    #funcFit2DProjectionCosThetaSigmaFree[i].SetLineColor(kBlue)
    #funcFit2DProjectionCosThetaSigmaFree[i].SetParameter(0,parameters[0]); funcFit2DProjectionCosThetaSigmaFree[i].SetParameter(1,parameters[1]);
    #histLambdaThetaFit2DProjectionSigmaFree.SetBinContent(i+2,parameters[1]); histLambdaThetaFit2DProjectionSigmaFree.SetBinError(i+2,errorParameters[1]);

    #parameters = SimFit2DProjectionSigmaFree.GetPhiParametersList(); errorParameters = SimFit2DProjectionSigmaFree.GetErrorPhiParametersList();
    #funcFit2DProjectionPhiSigmaFree.append(TF1("funcFit2DProjectionPhiSigmaFree_" + namePtRangesNew[i],MathFuncs.MyMathFuncs.MyFuncPolPhi,0.,PI,3))
    #funcFit2DProjectionPhiSigmaFree[i].SetLineColor(kBlue)
    #funcFit2DProjectionPhiSigmaFree[i].SetParameter(0,parameters[0]); funcFit2DProjectionPhiSigmaFree[i].SetParameter(1,parameters[1]); funcFit2DProjectionPhiSigmaFree[i].SetParameter(2,parameters[2])
    #histLambdaPhiFit2DProjectionSigmaFree.SetBinContent(i+2,parameters[2]); histLambdaPhiFit2DProjectionSigmaFree.SetBinError(i+2,errorParameters[2]);

    del SimFit2DProjectionSigmaFixed
    #del SimFit2DProjectionSigmaFree

# Control plots
canvasProjNJpsiCorr = TCanvas("canvasProjNJpsiCorr","canvasProjNJpsiCorr",20,20,1400,1200)
canvasProjNJpsiCorr.Divide(3,2)
canvasProjNJpsiCorr.cd(1);
histGridProjCosThetaNJpsi[0].Draw(); histProjCosThetaNJpsiSigmaFixedCorrClone[0].Draw("Esame"); #histProjCosThetaNJpsiSigmaFreeCorrClone[0].Draw("Esame");
funcFit2DProjectionCosThetaSigmaFixed[0].Draw("Esame"); funcPolProjectionCosThetaSigmaFixed[0].Draw("Esame"); #funcFit2DProjectionCosThetaSigmaFree[0].Draw("Esame");
canvasProjNJpsiCorr.cd(2);
histGridProjCosThetaNJpsi[1].Draw(); histProjCosThetaNJpsiSigmaFixedCorrClone[1].Draw("Esame"); #histProjCosThetaNJpsiSigmaFreeCorrClone[1].Draw("Esame");
funcFit2DProjectionCosThetaSigmaFixed[1].Draw("Esame"); funcPolProjectionCosThetaSigmaFixed[1].Draw("Esame"); #funcFit2DProjectionCosThetaSigmaFree[1].Draw("Esame");
canvasProjNJpsiCorr.cd(3);
histGridProjCosThetaNJpsi[2].Draw(); histProjCosThetaNJpsiSigmaFixedCorrClone[2].Draw("Esame"); #histProjCosThetaNJpsiSigmaFreeCorrClone[2].Draw("Esame");
funcFit2DProjectionCosThetaSigmaFixed[2].Draw("Esame"); funcPolProjectionCosThetaSigmaFixed[2].Draw("Esame"); #funcFit2DProjectionCosThetaSigmaFree[2].Draw("Esame");
canvasProjNJpsiCorr.cd(4);
histGridProjPhiNJpsi[0].Draw(); histProjPhiNJpsiSigmaFixedCorrClone[0].Draw("Esame"); #histProjPhiNJpsiSigmaFreeCorrClone[0].Draw("Esame");
funcFit2DProjectionPhiSigmaFixed[0].Draw("Esame"); funcPolProjectionPhiSigmaFixed[0].Draw("Esame"); #funcFit2DProjectionPhiSigmaFree[0].Draw("Esame");
canvasProjNJpsiCorr.cd(5);
histGridProjPhiNJpsi[1].Draw(); histProjPhiNJpsiSigmaFixedCorrClone[1].Draw("Esame"); #histProjPhiNJpsiSigmaFreeCorrClone[1].Draw("Esame");
funcFit2DProjectionPhiSigmaFixed[1].Draw("Esame"); funcPolProjectionPhiSigmaFixed[1].Draw("Esame"); #funcFit2DProjectionPhiSigmaFree[1].Draw("Esame");
canvasProjNJpsiCorr.cd(6);
histGridProjPhiNJpsi[2].Draw(); histProjPhiNJpsiSigmaFixedCorrClone[2].Draw("Esame"); #histProjPhiNJpsiSigmaFreeCorrClone[2].Draw("Esame");
funcFit2DProjectionPhiSigmaFixed[2].Draw("Esame"); funcPolProjectionPhiSigmaFixed[2].Draw("Esame"); #funcFit2DProjectionPhiSigmaFree[2].Draw("Esame");

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
    histAccxEffCost.append(fileAccxEff1D.Get("histAccxEffCostReWeighted_" + namePtRangesNew[i])); histAccxEffCost[i].SetLineColor(kBlack); histAccxEffCost[i].SetDirectory(0);
    histAccxEffPhi.append(fileAccxEff1D.Get("histAccxEffPhiReWeighted_" + namePtRangesNew[i])); histAccxEffPhi[i].SetLineColor(kBlack); histAccxEffPhi[i].SetDirectory(0);

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
    SimFitPure1D.SimultaneousFit2(histNJpsiCostCorr[i],histNJpsiPhiCorr[i],minFitRangeCost,maxFitRangeCost,"")

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

#for i in range(3):
    #histAccxEffCost[i].Scale(1./(histAccxEffCost[i].Integral()))
    #histAccxEffPhi[i].Scale(1./(histAccxEffPhi[i].Integral()))

#canvasAccxEffComp = TCanvas("canvasAccxEffComp","canvasAccxEffComp",20,20,1400,1200)
#canvasAccxEffComp.Divide(3,2)
#canvasAccxEffComp.cd(1); histAccxEffProjectionCosTheta[0].Draw("E"); histAccxEffCost[0].Draw("Esame");
#canvasAccxEffComp.cd(2); histAccxEffProjectionCosTheta[1].Draw("E"); histAccxEffCost[1].Draw("Esame");
#canvasAccxEffComp.cd(3); histAccxEffProjectionCosTheta[2].Draw("E"); histAccxEffCost[2].Draw("Esame");
#canvasAccxEffComp.cd(4); histAccxEffProjectionPhi[0].Draw("E"); histAccxEffPhi[0].Draw("Esame");
#canvasAccxEffComp.cd(5); histAccxEffProjectionPhi[1].Draw("E"); histAccxEffPhi[1].Draw("Esame");
#canvasAccxEffComp.cd(6); histAccxEffProjectionPhi[2].Draw("E"); histAccxEffPhi[2].Draw("Esame");

canvasPure1DNJpsiCorr = TCanvas("canvasPure1DNJpsiCorr","canvasPure1DNJpsiCorr",20,20,1400,1200)
canvasPure1DNJpsiCorr.Divide(3,2)
canvasPure1DNJpsiCorr.cd(1);
histGridProjCosThetaNJpsi[0].Draw(); histNJpsiCostCorr[0].Draw("Esame"); histProjCosThetaNJpsiSigmaFixedCorrClone[0].Draw("Esame");
funcFit1DCosTheta[0].Draw("Esame")
canvasPure1DNJpsiCorr.cd(2);
histGridProjCosThetaNJpsi[1].Draw(); histNJpsiCostCorr[1].Draw("Esame"); histProjCosThetaNJpsiSigmaFixedCorrClone[1].Draw("Esame");
funcFit1DCosTheta[1].Draw("Esame")
canvasPure1DNJpsiCorr.cd(3);
histGridProjCosThetaNJpsi[2].Draw(); histNJpsiCostCorr[2].Draw("Esame"); histProjCosThetaNJpsiSigmaFixedCorrClone[2].Draw("Esame");
funcFit1DCosTheta[2].Draw("Esame")
canvasPure1DNJpsiCorr.cd(4);
histGridProjPhiNJpsi[0].Draw(); histNJpsiPhiCorr[0].Draw("Esame"); histProjPhiNJpsiSigmaFixedCorrClone[0].Draw("Esame");
funcFit1DPhi[0].Draw("Esame")
canvasPure1DNJpsiCorr.cd(5);
histGridProjPhiNJpsi[1].Draw(); histNJpsiPhiCorr[1].Draw("Esame"); histProjPhiNJpsiSigmaFixedCorrClone[1].Draw("Esame");
funcFit1DPhi[1].Draw("Esame")
canvasPure1DNJpsiCorr.cd(6);
histGridProjPhiNJpsi[2].Draw(); histNJpsiPhiCorr[2].Draw("Esame"); histProjPhiNJpsiSigmaFixedCorrClone[2].Draw("Esame");
funcFit1DPhi[2].Draw("Esame")

################################################################################
# Polarization parameters comparison
lineZero = TLine(0.,0.,15.,0.); lineZero.SetLineStyle(kDashed);

histGridLambdaTheta = TH2D("histGridLambdaTheta","#lambda_{#it{#theta}} vs #it{p}_{T}",100,0,15,100,-1.,1.)
histGridLambdaTheta.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
histGridLambdaTheta.GetYaxis().SetTitle("#lambda_{#it{#theta}}"); histGridLambdaTheta.GetYaxis().SetTitleOffset(1.2);

histGridLambdaPhi = TH2D("histGridLambdaPhi","#lambda_{#it{#varphi}} vs #it{p}_{T}",100,0,15,100,-1.,1.)
histGridLambdaPhi.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
histGridLambdaPhi.GetYaxis().SetTitle("#lambda_{#it{#varphi}}"); histGridLambdaPhi.GetYaxis().SetTitleOffset(1.2);

histGridLambdaThetaPhi = TH2D("histGridLambdaThetaPhi","#lambda_{#it{#theta#varphi}} vs #it{p}_{T}",100,0,15,100,-1.,1.)
histGridLambdaThetaPhi.GetXaxis().SetTitle("#it{p}_{T} (GeV/#it{c})")
histGridLambdaThetaPhi.GetYaxis().SetTitle("#lambda_{#it{#theta#varphi}}"); histGridLambdaPhi.GetYaxis().SetTitleOffset(1.2);

canvasPolarizationParameters = TCanvas("canvasPolarizationParameters","canvasPolarizationParameters",20,20,1400,600)
canvasPolarizationParameters.Divide(3,1)

canvasPolarizationParameters.cd(1)
histGridLambdaTheta.Draw()
lineZero.Draw("same")
histLambdaThetaPure1D.Draw("EPsame")
histLambdaThetaFit2DProjectionSigmaFixed.Draw("EPsame")
histLambdaThetaSigmaFixed.Draw("EPsame")
legendLambdaThetaSigmaFixed = TLegend(0.55,0.75,0.9,0.9); legendLambdaThetaSigmaFixed.SetTextSize(0.03);
legendLambdaThetaSigmaFixed.AddEntry(histLambdaPhiPure1D,"1D approach","lp")
legendLambdaThetaSigmaFixed.AddEntry(histLambdaPhiFit2DProjectionSigmaFixed,"2D projection #sigma_{J/#psi}^{Fixed}","l")
legendLambdaThetaSigmaFixed.AddEntry(histLambdaPhiSigmaFixed,"2D approach #sigma_{J/#psi}^{Fixed}","l")
legendLambdaThetaSigmaFixed.Draw("same")

canvasPolarizationParameters.cd(2)
histGridLambdaPhi.Draw()
lineZero.Draw("same")
histLambdaPhiPure1D.Draw("EPsame")
histLambdaPhiFit2DProjectionSigmaFixed.Draw("EPsame")
histLambdaPhiSigmaFixed.Draw("EPsame")
legendLambdaPhiSigmaFixed = TLegend(0.55,0.75,0.9,0.9); legendLambdaPhiSigmaFixed.SetTextSize(0.03);
legendLambdaPhiSigmaFixed.AddEntry(histLambdaPhiPure1D,"1D approach","lp")
legendLambdaPhiSigmaFixed.AddEntry(histLambdaPhiFit2DProjectionSigmaFixed,"2D projection #sigma_{J/#psi}^{Fixed}","l")
legendLambdaPhiSigmaFixed.AddEntry(histLambdaPhiSigmaFixed,"2D approach #sigma_{J/#psi}^{Fixed}","l")
legendLambdaPhiSigmaFixed.Draw("same")

canvasPolarizationParameters.cd(3)
histGridLambdaThetaPhi.Draw()
lineZero.Draw("same")
histLambdaThetaPhiSigmaFixed.Draw("EPsame")
legendLambdaThetaPhiSigmaFixed = TLegend(0.55,0.75,0.9,0.9); legendLambdaThetaPhiSigmaFixed.SetTextSize(0.03);
legendLambdaThetaPhiSigmaFixed.AddEntry(histLambdaPhiSigmaFixed,"2D approach #sigma_{J/#psi}^{Fixed}","l")
legendLambdaThetaPhiSigmaFixed.Draw("same")

#canvasPolarizationParameters.cd(4)
#histGridLambdaTheta.Draw()
#lineZero.Draw("same")
#histLambdaThetaPure1D.Draw("EPsame")
#histLambdaThetaFit2DProjectionSigmaFree.Draw("EPsame")
#histLambdaThetaSigmaFree.Draw("EPsame")
#legendLambdaThetaSigmaFree = TLegend(0.55,0.75,0.9,0.9); legendLambdaThetaSigmaFree.SetTextSize(0.03);
#legendLambdaThetaSigmaFree.AddEntry(histLambdaPhiPure1D,"1D approach","lp")
#legendLambdaThetaSigmaFree.AddEntry(histLambdaPhiFit2DProjectionSigmaFree,"2D projection #sigma_{J/#psi}^{Free}","l")
#legendLambdaThetaSigmaFree.AddEntry(histLambdaPhiSigmaFree,"2D approach #sigma_{J/#psi}^{Free}","l")
#legendLambdaThetaSigmaFree.Draw("same")

#canvasPolarizationParameters.cd(5)
#histGridLambdaPhi.Draw()
#lineZero.Draw("same")
#histLambdaPhiPure1D.Draw("EPsame")
#histLambdaPhiFit2DProjectionSigmaFree.Draw("EPsame")
#histLambdaPhiSigmaFree.Draw("EPsame")
#legendLambdaPhiSigmaFree = TLegend(0.55,0.75,0.9,0.9); legendLambdaPhiSigmaFree.SetTextSize(0.03);
#legendLambdaPhiSigmaFree.AddEntry(histLambdaPhiPure1D,"1D approach","lp")
#legendLambdaPhiSigmaFree.AddEntry(histLambdaPhiFit2DProjectionSigmaFree,"2D projection #sigma_{J/#psi}^{Free}","l")
#legendLambdaPhiSigmaFree.AddEntry(histLambdaPhiSigmaFree,"2D approach #sigma_{J/#psi}^{Free}","l")
#legendLambdaPhiSigmaFree.Draw("same")



raw_input()
