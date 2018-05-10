import ROOT
import sys,os
import argparse
sys.path.append('../GPBackgrounds')
import RooGP
import RooGPBkg
from SignalModel import *


def run(args, gpConfig, dataHisto):
    myy = ROOT.RooRealVar('myy','myy',55,120)
    nSigGP = ROOT.RooRealVar('nSigGP','nSigGP',-2000,2000)

    fixedMyy = ROOT.RooRealVar('fixMyy', 'fixMyy', 55, 120)
    nSigGP = ROOT.RooRealVar('fixnSigGP','fixnSigGP',-2000,2000)

    #The PDF for the GP bkg + signal DSCB
    #GPpdf = RooGP.RooGP("mypdf", "CustomPDF",myy ,nSigGP, args.sigMass, trainHisto, dataHisto)
    GPpdf = RooGP.RooGP("mypdf", "CustomPDF",myy ,nSigGP, args.sigMass, gpConfig)

    #get data from histogram to fit to
    data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(myy), dataHisto)

    #Fit the pdf to the data
    fitResult = GPpdf.fitTo(data, ROOT.RooFit.Save())

    # Get negative log likelihoods
    nllGP = GPpdf.createNLL(data)
    profileGP = nllGP.createProfile(ROOT.RooArgSet(nSigGP))

    #Get the best fit nSig and use that to get the corresponding GP background prediction.
    #print "Fit Result: ",fitResultGP.Print()
    fixedGPpdf = RooGP.RooGP("fixedPdf", "fixedPdf", fixedMyy, fixnSigGP, args.sigMass, gpConfig)

    optBkgGPHisto = fixedGPpdf.getGPHisto(nSigGP.getValV())

    ########## Fixed background ######################
    #Fixed Background PDFs
    sigPdf = ROOT.HggTwoSidedCBPdf("sigPdf", "DSCB sig", myy, m0, sigma, alphalo,nlo,alphahi,nHi)

    #dh =  RooDataHist("mc","mc",  RooArgList(myy), trainHisto)
    dh =  ROOT.RooDataHist("mc","mc",  RooArgList(myy), optBkgGPHisto)
    BkgPdf =  ROOT.RooHistPdf("noSmooth", "HS MC PDF",  RooArgSet(myy), dh)

    pdf =  ROOT.RooAddPdf("sigBkg", "Signal + Background",  RooArgList(sigPdf, BkgPdf),  RooArgList(nSig, nBkg))
    fitResultNorm = pdf.fitTo(data,  RooFit.Save())


    # Plot the background fit and data
    c1 = ROOT.TCanvas('BkgFit','Background Fit')
    frame = myy.frame()
    data.plotOn(frame)
    GPpdf.plotOn(frame)
    frame.Draw()
    nSigGP.Print()
    c1.Print(args.outDir+'/test_GP.pdf')

    c2 = ROOT.TCanvas('nll', 'nll')
    nSigFrame = nSigGP.frame()
    profileGP.plotOn(nSigFrame)
    nSigFrame.Draw()
    c2.Print(args.outDir+'/NLL.pdf')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default='../data/bkg_comb_Ioannis_h021.root',help="Input root file with background histogram")
    parser.add_argument("--outDir", help="Name tag to add to begining of matplotlib plot")
    parser.add_argument("--doSig","-s", action='store_true', help="Inject signal into background and make some signal plots")

    args = parser.parse_args()
    #args.CL = '95'
    args.CL = '68'
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    args.sigMass = 85
    #stats = 1000000

    f = ROOT.TFile(args.input)

    trainHisto = f.Get('hbkg_tmp_nominal_UU_incl')
    #norm = 333000/ trainHisto.Integral()
    #trainHisto.Scale(norm)
    #trainHisto.Rebin(8)
    dataHisto = trainHisto.Clone('dataHisto')
    #dataHisto.Reset()
    #dataHisto = ROOT.TH1F('dataHisto', 'Background data', trainHisto.GetNbinsX(), trainHisto.GetBinLowEdge(1), trainHisto.GetBinLowEdge(trainHisto.GetNbinsX())+trainHisto.GetBinWidth(1))

    ROOT.gRandom.SetSeed(2)

    #dataHisto.FillRandom(trainHisto, stats)
    if args.doSig:
        #Make a signal toy and add it to the
        datasighist = buildSignal(125,125, dataHisto.GetNbinsX())
        dataHisto.Add(datasighist)


    gpConfig = {
                'train_range'    : [55,120]
                ,'trainHisto' : trainHisto
                #'lengthScale' : 60
                #,'lengthScale_min': 1
                #,'lengthScale_max': 200
                #,'amplitude'      : 100
                #,'amplitude_min'      : 1e-3
                #,'amplitude_max'      : 1e15
                #,'priorMean' : priorMean
    }

    run(args, trainHisto, dataHisto)
    print "Signal Mass Histo Bin:", trainHisto.GetBinContent(trainHisto.FindBin(115)),"+=",trainHisto.GetBinError(trainHisto.FindBin(115))
