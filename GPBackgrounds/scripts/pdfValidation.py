import ROOT
import sys,os
import argparse
sys.path.append('../GPBackgrounds')
import RooGP
import RooGPBkg
from SignalModel import *


def run(args, trainHisto, dataHisto):
    myy = ROOT.RooRealVar('myy','myy',105,159)
    nSigGP = ROOT.RooRealVar('nSigGP','nSigGP',-100)

    #The PDF for the GP bkg + signal DSCB
    GPpdf = RooGP.RooGP("mypdf", "CustomPDF",myy ,nSigGP, args.sigMass, trainHisto, dataHisto)

    #get data from histogram to fit to
    data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(myy), dataHisto)

    #Fit the pdf to the data
    #fitResult = GPpdf.fitTo(data, ROOT.RooFit.Save())
    gphisto100 = GPpdf.getGPHisto(-1100)
    gphisto125 = GPpdf.getGPHisto(-1300)


    # Get negative log likelihoods
    #nllGP = GPpdf.createNLL(data)
    #profileGP = nllGP.createProfile(ROOT.RooArgSet(nSigGP))

    # Plot the background fit and data
    c1 = ROOT.TCanvas('BkgFit','Background Fit')
    #frame = myy.frame()
    #data.plotOn(frame)
    #GPpdf.plotOn(frame)
    #frame.Draw()
    #nSigGP.Print()
    dataHisto.SetMarkerStyle(20)
    gphisto125.SetLineColor(ROOT.kRed)
    gphisto100.Draw()
    gphisto125.Draw('same')
    dataHisto.Draw('samep')
    c1.Print(args.outDir+'/test_GP.pdf')

    """
    c2 = ROOT.TCanvas('nll', 'nll')
    nSigFrame = nSigGP.frame()
    profileGP.plotOn(nSigFrame)
    nSigFrame.Draw()
    c2.Print(args.outDir+'/NLL.pdf')
    """

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input root file with background histogram")
    parser.add_argument("--outDir", help="Name tag to add to begining of matplotlib plot")
    parser.add_argument("--doSig","-s", action='store_true', help="Inject signal into background and make some signal plots")

    args = parser.parse_args()
    #args.CL = '95'
    args.CL = '68'
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    args.sigMass = 125
    stats = 100000

    f = ROOT.TFile(args.input)

    trainHisto = f.Get('hmgg_c0')
    trainHisto.Rebin(8)
    #dataHisto = trainHisto.Clone('dataHisto')
    #dataHisto.Reset()
    dataHisto = ROOT.TH1F('dataHisto', 'Background data', 27, 105, 159)

    ROOT.gRandom.SetSeed(2)

    dataHisto.FillRandom(trainHisto, stats)
    if args.doSig:
        #Make a signal toy and add it to the
        datasighist = buildSignal(125,1000, dataHisto.GetNbinsX())
        dataHisto.Add(datasighist)

    run(args, trainHisto, dataHisto)
