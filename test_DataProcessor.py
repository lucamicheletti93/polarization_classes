from ROOT import *
from array import array
import os.path

gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("0.2g");
gROOT.ProcessLineSync(".x DataProcessor.cxx+")

if os.path.isfile("/media/luca/488AE2208AE20A70/PbPb_2015_Trees_pDCA/new_tree/Tree_full_stat.root"):
    fileData = TFile.Open("/media/luca/488AE2208AE20A70/PbPb_2015_Trees_pDCA/new_tree/Tree_full_stat.root")  # local
else:
    print "This class can run only in local on the hard disk!!!"

treeData = fileData.Get("PbPbTree")

test = DataProcessor(treeData)
#test.ComputeTriggerResponseFunction("FullStat","output/TriggerRespondeFunctionData.root")

fileTriggerResponseFunctionData = TFile.Open("output/TriggerRespondeFunctionData.root")
histTriggerResponseFunctionSM = fileTriggerResponseFunctionData.Get("histTriggerResponseFunctionSM")
histTriggerResponseFunctionSMpDCA = fileTriggerResponseFunctionData.Get("histTriggerResponseFunctionSMpDCA")

canvasTriggerResponseFunctionComp = TCanvas("canvasTriggerResponseFunctionComp","canvasTriggerResponseFunctionComp",20,20,600,600)
histTriggerResponseFunctionSM.Draw()
histTriggerResponseFunctionSMpDCA.Draw("same")

raw_input()
