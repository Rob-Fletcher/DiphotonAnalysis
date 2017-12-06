from ROOT import *
import sys,os
import argparse
sys.path.append('../GPBackgrounds')
import RooGP
import RooGPBkg
from SignalModel import *


def run(args, trainHisto, dataHisto):
    sigMass = 125
    myy =  RooRealVar('myy','myy',105,160)
    #nSigGP =  RooRealVar('nSigGP','nSigGP',-500,500)
    nSigGP =  RooRealVar('nSigGP','nSigGP',-100, 500)
    #nSig =  RooRealVar('nSig','nSig',-500,500)
    nSig =  RooRealVar('nSig','nSig',-100, 500)
    nBkg =  RooRealVar('nBkg', 'nBkg', 1,90000000 )

    m0 = RooRealVar('m0','m0', sigMass)
    sigma = RooRealVar('sigma','sigma', 1.86)
    alphalo = RooRealVar('  alphalo','  alphalo',   1.475)
    nlo = RooRealVar('nlo','nlo', 12.1)
    alphahi = RooRealVar('alphahi','alphahi', 1.902)
    nHi = RooRealVar('nHi','nHi', 11.6)

    #Gaussian process PDF
    pdfGP = RooGP.RooGP("mypdf", "CustomPDF",myy ,nSigGP, sigMass, trainHisto, dataHisto)
    #pdf = RooGPBkg.RooGPBkg("bkgPDF", "bkg only PDF",myy, trainHisto, dataHisto)


    #The data to fit to
    data =  RooDataHist("dh", "dh",  RooArgList(myy), dataHisto)

    #GP fit result to the data histogram
    fitResultGP = pdfGP.fitTo(data,  RooFit.Save())

    #Get the best fit nSig and use that to get the corresponding GP background prediction.
    #print "Fit Result: ",fitResultGP.Print()

    #optBkgGPHisto = pdfGP.getGPHisto(nSigGP.getValV())

    """
    ########## Fixed background ######################
    #Fixed Background PDFs
    sigPdf = HggTwoSidedCBPdf("sigPdf", "DSCB sig", myy, m0, sigma, alphalo,nlo,alphahi,nHi)

    dh =  RooDataHist("mc","mc",  RooArgList(myy), trainHisto)
    #dh =  RooDataHist("mc","mc",  RooArgList(myy), optBkgGPHisto)
    BkgPdf =  RooHistPdf("noSmooth", "HS MC PDF",  RooArgSet(myy), dh)

    pdf =  RooAddPdf("sigBkg", "Signal + Background",  RooArgList(sigPdf, BkgPdf),  RooArgList(nSig, nBkg))
    fitResultNorm = pdf.fitTo(data,  RooFit.Save())
"""

    c1 =  TCanvas('c1','c1')
    #c2 =  TCanvas('c2','c2')
    frame = myy.frame()
    #frame2 = nSig.frame()
    frame3 = nSigGP.frame()
    nllGP = pdfGP.createNLL(data)
    pllGP = nllGP.createProfile( RooArgSet(nSigGP))
    #nll = pdf.createNLL(data)
    #pll = nll.createProfile( RooArgSet(nSig))
    data.plotOn(frame)
    pdfGP.plotOn(frame, RooFit.LineColor(kRed))
    #optBkgGPHisto.Draw()
#    pdf.plotOn(frame)
    pllGP.plotOn(frame3, RooFit.LineColor(kRed))
    #pll.plotOn(frame2)
    frame.Draw()
    #dataHisto.Draw()
    c1.Print(args.outDir+'/test_GP_both.pdf')
    #c2.cd()
    frame3.Draw()
    frame2.Draw('same')
    c1.Print(args.outDir+'/NLL_GP_both.pdf')
    print "Evaluate result:",pdfGP.evaluate(130, 20.0)
#    nllGP.Print()
#    nSig.Print()
#    nSigGP.Print()

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input  file with background histogram")
    parser.add_argument("--outDir", help="Name tag to add to begining of matplotlib plot")
    parser.add_argument("--doSig","-s", action='store_true', help="Inject signal into background and make some signal plots")

    args = parser.parse_args()
    #args.CL = '95'
    args.CL = '68'
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    stats = 100000

    f =  TFile(args.input)

    trainHisto = f.Get('hmgg_c0')
    trainHisto.Rebin(8)
    dataHisto = trainHisto.Clone('dataHisto')
    dataHisto.Reset()
    #dataHisto = TH1F("dataHisto", "dataHisto", 27, 105, 160)

    gRandom.SetSeed(4)

    dataHisto.FillRandom(trainHisto, stats)
    print "nEntries in dataHisto", dataHisto.GetEntries()
    if args.doSig:
        #Make a signal toy and add it to the
        datasighist = buildSignal(125,200, dataHisto.GetNbinsX())
        dataHisto.Add(datasighist)

    run(args, trainHisto, dataHisto)
