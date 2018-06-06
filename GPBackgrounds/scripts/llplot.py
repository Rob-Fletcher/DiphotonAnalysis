import ROOT
import sys,os
import argparse
sys.path.append('../GPBackgrounds')
import RooGP
#import RooGPBkg
from SignalModel import *


def run(args, gpConfig, dataHisto):
    myy = ROOT.RooRealVar('myy','myy',60,120)
    nSigGP = ROOT.RooRealVar('nSigGP','nSigGP',-2000,2000)

    #The PDF for the GP bkg + signal DSCB
    #GPpdf = RooGP.RooGP("mypdf", "CustomPDF",myy ,nSigGP, args.sigMass, trainHisto, dataHisto)
    GPpdf = RooGP.RooGP("mypdf", "CustomPDF",myy ,nSigGP, args.sigMass, gpConfig)

    #get data from histogram to fit to
    data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(myy), dataHisto)

    #Fit the pdf to the data
    fitResult = GPpdf.fitTo(data,ROOT.RooFit.SumW2Error(ROOT.kTRUE), ROOT.RooFit.Save())
    #GPpdf.fz(True)
    GPHisto = GPpdf.getGPHisto(nSigGP.getValV())[0]

    # Get negative log likelihoods
    nllGP = GPpdf.createNLL(data)
    profileGP = nllGP.createProfile(ROOT.RooArgSet(nSigGP))

    ########## Fixed background ######################
    nSig = ROOT.RooRealVar('nSig','nSig',-2000,2000)
    nBkg = ROOT.RooRealVar('nBkg','nBkg', 0, 9e15)

    fixedMyy = ROOT.RooRealVar('fixMyy', 'fixMyy', 60, 120)
    fixnSigGP = ROOT.RooRealVar('fixnSigGP','fixnSigGP',-2000,2000)


    mass = args.sigMass

    mCBSignal0 =  0.048803
    mCBSignal1 = -0.00512848
    sCBSignal0 =  1.3612
    sCBSignal1 =  0.78698
    aLoSignal0 =  1.6738
    aLoSignal1 = -0.0159122
    aHiSignal0 =  1.4824
    aHiSignal1 =  0.12803
    nLoSignal0 =  9.5595
    nHiSignal0 =  72085

    mnX = (mass - 100)/100
    sigma  = ROOT.RooRealVar("sigma","sigma",(sCBSignal0   + sCBSignal1*mnX))
    deltaMSignal  = (mCBSignal0 + mCBSignal1*mnX)
    alpha_low = ROOT.RooRealVar("alpha_low","alpha_low",(aLoSignal0 + aLoSignal1*mnX))
    n_low = ROOT.RooRealVar("n_low","n_low", (nLoSignal0))
    alpha_high = ROOT.RooRealVar("alpha_high","alpha_high",(aHiSignal0 + aHiSignal1*mnX))
    n_high = ROOT.RooRealVar("n_high","n_high",(nHiSignal0))

    Mx  = ROOT.RooRealVar("Mx","Mx",(mass + deltaMSignal))

    #Get the best fit nSig and use that to get the corresponding GP background prediction.
    #print "Fit Result: ",fitResultGP.Print()

    #print "Getting background histogram with best fit nSig value:", nSigGP.getValV()
    #optBkgGPHisto = GPpdf.getGPHisto(nSigGP.getValV())[0]

    #Fixed Background PDFs
    sigPdf = ROOT.HggTwoSidedCBPdf("sigPdf", "DSCB sig", myy, Mx, sigma, alpha_low, n_low, alpha_high, n_high)

    #dh =  RooDataHist("mc","mc",  RooArgList(myy), trainHisto)
    #dh =  ROOT.RooDataHist("mc","mc",  ROOT.RooArgList(myy), optBkgGPHisto)
    #BkgPdf =  ROOT.RooHistPdf("noSmooth", "HS MC PDF",  ROOT.RooArgSet(myy), dh)

    pdf =  ROOT.RooAddPdf("sigBkg", "Signal + Background",  ROOT.RooArgList(sigPdf, GPpdf),  ROOT.RooArgList(nSig, nBkg))
    fitResultNorm = pdf.fitTo(data,  ROOT.RooFit.Save())

    fixnllGP = pdf.createNLL(data)
    fixprofileGP = fixnllGP.createProfile(ROOT.RooArgSet(nSig))

    ########## End fixed background

    # Plot the background fit and data
    c1 = ROOT.TCanvas('BkgFit','Background Fit')
    frame = myy.frame()
    data.plotOn(frame)
    GPpdf.plotOn(frame)
    frame.Draw()
    #GPHisto.Draw('same')
    #nSigGP.Print()
    c1.Print(args.outDir+'/test_GP.pdf')

    c2 = ROOT.TCanvas('nll', 'nll')
    nSigFrame = nSigGP.frame()
    profileGP.plotOn(nSigFrame)
    fixprofileGP.plotOn(nSigFrame)
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

    """
    ROOT.gRandom.SetSeed(2)

    #dataHisto.FillRandom(trainHisto, stats)
    if args.doSig:
        #Make a signal toy and add it to the
        datasighist = buildSignal(125,125, dataHisto.GetNbinsX())
        dataHisto.Add(datasighist)
    """

    gpConfig = {
                'train_range'    : [60,120]
                ,'trainHisto' : trainHisto
                #'lengthScale' : 60
                #,'lengthScale_min': 1
                #,'lengthScale_max': 200
                #,'amplitude'      : 100
                #,'amplitude_min'      : 1e-3
                #,'amplitude_max'      : 1e15
                #,'priorMean' : priorMean
    }

    run(args, gpConfig, dataHisto)
