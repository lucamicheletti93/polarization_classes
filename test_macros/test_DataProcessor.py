from ROOT import *
from array import array
import os.path
import string

gStyle.SetOptStat(0)
gStyle.SetPaintTextFormat("0.2g");
gROOT.ProcessLineSync(".x ../DataProcessor.cxx+")
gROOT.ProcessLineSync(".x ../Binning.cxx+")

fileBinning = TFile.Open("output/binning.root")
binning = fileBinning.Get("Binning")
CostValues = binning.GetCostValues()
PhiValues = binning.GetPhiValues()

#print "#######################################################################"
#print "Create filtered histograms from the original trees"
#fileRunList = open("run_list.txt","r")
#for eachLine in fileRunList:
    #runNumber = ''
    #for char in eachLine:
        #if char.isdigit():
            #runNumber += char
    #if runNumber.isdigit():
        #if os.path.isfile('/media/luca/488AE2208AE20A70/PbPb_2015_Trees/Tree_' + runNumber + '.root'):
            #print os.path.join('/media/luca/488AE2208AE20A70/PbPb_2015_Trees/Tree_' + runNumber + '.root')
            #fileData = TFile.Open('/media/luca/488AE2208AE20A70/PbPb_2015_Trees/Tree_' + runNumber + '.root')  # local
            #treeData = fileData.Get("PbPbTree")
            #test = DataProcessor(treeData)
            #test.CreateFilteredTrees("FullStat",'output/filtered_trees/TreeDataFiltered_' + runNumber + '.root')
            #fileData.Close()
            #del test
    #else:
        #print os.path.join('File /media/luca/488AE2208AE20A70/PbPb_2015_Trees/Tree_' + runNumber + '.root not found')

#fileRunList.close()

#print "#######################################################################"
print "Create histograms from the filtered trees"
fileHistMass = TFile.Open('output/filtered_trees/TreeDataFiltered_246994.root')
test = DataProcessor()
test.SetBinning(CostValues,PhiValues)
test.CreateInvMassHistograms(fileHistMass,'output/HistMass_246994.root')

#print "#######################################################################"
#print "Compute the trigger response function from data"
#if os.path.isfile("/media/luca/488AE2208AE20A70/PbPb_2015_Trees_pDCA/new_tree/Tree_full_stat.root"):
    #fileData = TFile.Open("/media/luca/488AE2208AE20A70/PbPb_2015_Trees_pDCA/new_tree/Tree_full_stat.root")  # local
    #treeData = fileData.Get("PbPbTree")
    #test = DataProcessor(treeData)
    #test.ComputeTriggerResponseFunction("FullStat","output/TriggerRespondeFunctionData.root")

#else:
    #print "This class can be run only in local on the hard disk!!!"

#if os.path.isfile("output/TriggerRespondeFunctionData.root"):
    #fileTriggerResponseFunctionData = TFile.Open("output/TriggerRespondeFunctionData.root")
    #histTriggerResponseFunctionSM = fileTriggerResponseFunctionData.Get("histTriggerResponseFunctionSM")
    #histTriggerResponseFunctionSMpDCA = fileTriggerResponseFunctionData.Get("histTriggerResponseFunctionSMpDCA")

    #canvasTriggerResponseFunctionComp = TCanvas("canvasTriggerResponseFunctionComp","canvasTriggerResponseFunctionComp",20,20,600,600)
    #histTriggerResponseFunctionSM.Draw()
    #histTriggerResponseFunctionSMpDCA.Draw("same")

#else:
    #print "Connect the hard disk to produce the file output/TriggerRespondeFunctionData.root!!!"

raw_input()
