import ROOT
import sys,os
import argparse
sys.path.append('../GPBackgrounds')
import RooGP
import RooGPBkg
from SignalModel import *
"""Check roofit for option called 'offset' in RooAbsPdf. This caches some offset
value and might be causing the discontinuity.

Another option is that it doesnt know that the pdf is binned. Look for an option
to make either the pdf or the variable binned.
"""


def run(args, trainHisto, dataHisto):
    myy = ROOT.RooRealVar('myy','myy',args.lowRange,args.hiRange)
    nSigGP = ROOT.RooRealVar('nSigGP','nSigGP',-500,500)

    #The PDF for the GP bkg + signal DSCB
    GPpdf = RooGP.RooGP("mypdf", "CustomPDF",myy ,nSigGP, args.sigMass, trainHisto, dataHisto, [args.lowRange,args.hiRange])

    #get data from histogram to fit to
    data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(myy), dataHisto)

    #Fit the pdf to the data
    fitResult = GPpdf.fitTo(data, ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.Save())
    print "-------RooRealVar min at:", nSigGP.getVal(), nSigGP.getError()
    print "-------FitResult min at:", fitResult.floatParsFinal().Print("nSigGP")
    gphisto = GPpdf.getGPHisto(fitResult.floatParsFinal().find('nSigGP').getVal())
    #gphisto125 = GPpdf.getGPHisto(-1300)


    # Get negative log likelihoods
    nllGP = GPpdf.createNLL(data, ROOT.RooCmdArg("OffsetLikelihood"))
    profileGP = nllGP.createProfile(ROOT.RooArgSet(nSigGP))

    # Plot the background fit and data
    c1 = ROOT.TCanvas('BkgFit','Background Fit')
    frame = myy.frame()
    data.plotOn(frame)
    GPpdf.plotOn(frame)
    frame.Draw()
    #nSigGP.Print()
    #dataHisto.SetMarkerStyle(20)
    gphisto.SetLineColor(ROOT.kRed)
    #gphisto100.Draw()
    gphisto.Draw('same')
    #dataHisto.Draw('samep')
    c1.Print(args.outDir+'/test_GP.pdf')

    c2 = ROOT.TCanvas('nll', 'nll')
    nSigFrame = nSigGP.frame()
    nSigFrame.SetTitle("Min: "+str(nSigGP.getVal()))
    profileGP.plotOn(nSigFrame)
    nSigFrame.Draw()
    c2.Print(args.outDir+'/NLL.pdf')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="../data/CouplingsTemplatesMoriond2017.root",help="Input root file with background histogram")
    parser.add_argument("--outDir", help="Name tag to add to begining of matplotlib plot")
    parser.add_argument("--doSig","-s", action='store_true', help="Inject signal into background and make some signal plots")

    args = parser.parse_args()
    #args.CL = '95'
    args.CL = '68'
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    args.sigMass = 90
    args.lowRange = 70
    args.hiRange = 200
    stats = 100000

    f = ROOT.TFile(args.input)

    #trainHisto = f.Get('hmgg_c0')
    trainHisto = f.Get('m_yy')
    #trainHisto = f.Get('m_yy_c1_M17_ggH_0J_Cen_BkgTemplate')
    norm = 333000/ trainHisto.Integral()
    #norm = 1100000/ trainHisto.Integral()
    #norm = 333000/ trainHisto.Integral()
    trainHisto.Scale(norm)
    #trainHisto.Sumw2()
    #trainHisto.Rebin(8)
    #dataHisto = trainHisto.Clone('dataHisto')
    #dataHisto.Reset()
    dataHisto = ROOT.TH1F('dataHisto', 'Background data', 27, 105, 159)

    ROOT.gRandom.SetSeed(4)

    dataHisto.FillRandom(trainHisto, stats)
    if args.doSig:
        #Make a signal toy and add it to the
        datasighist = buildSignal(125,1000, dataHisto.GetNbinsX())
        dataHisto.Add(datasighist)

    #run(args, trainHisto, dataHisto)
    run(args, trainHisto, trainHisto)
