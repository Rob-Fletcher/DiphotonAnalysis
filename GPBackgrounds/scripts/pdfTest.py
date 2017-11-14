import ROOT
import sys,os
import argparse
sys.path.append('../GPBackgrounds')
import RooGP
import RooGPBkg
from SignalModel import *


def run(args, trainHisto, dataHisto):
    myy = ROOT.RooRealVar('myy','myy',105,160)
    nSig = ROOT.RooRealVar('nSig','nSig',-200,200)
    sigMass = 125
    pdf = RooGP.RooGP("mypdf", "CustomPDF",myy ,nSig, sigMass, trainHisto, dataHisto)
    #pdf = RooGPBkg.RooGPBkg("bkgPDF", "bkg only PDF",myy, trainHisto, dataHisto)

    data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(myy), dataHisto)

    c1 = ROOT.TCanvas('c1','c1')
    frame = myy.frame()
    data.plotOn(frame)
    fitResult = pdf.fitTo(data, ROOT.RooFit.Save())
    pdf.plotOn(frame)
    frame.Draw()
    nSig.Print()
    c1.Print(args.outDir+'/test_GP.pdf')


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

    stats = 100000

    f = ROOT.TFile(args.input)

    trainHisto = f.Get('hmgg_c0')
    trainHisto.Rebin(8)
    dataHisto = trainHisto.Clone('dataHisto')
    dataHisto.Reset()

    ROOT.gRandom.SetSeed(2)

    dataHisto.FillRandom(trainHisto, stats)
    if args.doSig:
        #Make a signal toy and add it to the
        datasighist = buildSignal(125,1000, dataHisto.GetNbinsX())
        dataHisto.Add(datasighist)

    run(args, trainHisto, dataHisto)
