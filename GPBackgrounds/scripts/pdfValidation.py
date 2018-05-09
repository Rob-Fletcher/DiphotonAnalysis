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

def fitData(args, histo, fitfunc):
    """Fit the data with a simple function to use as the mean of the GP.
    This is just for the purposes of demonstrating that as lengthScale gets long
    the GP approaches whatever function you use as the mean, thus the SS number
    also approaches the value found with this mean.

    """
    temp = histo.Clone('tmp')
    temp.Fit(fitfunc, 'L')
    func = temp.GetFunction(fitfunc).Clone('myfunc')
    return func

def run(args, trainHisto, dataHisto, gpConfig):
    myy = ROOT.RooRealVar('myy','myy',args.lowRange,args.hiRange)
    nSigGP = ROOT.RooRealVar('nSigGP','nSigGP',-20000,20000)

    #The PDF for the GP bkg + signal DSCB
    GPpdf = RooGP.RooGP("mypdf", "CustomPDF",myy ,nSigGP, args.sigMass, gpConfig)
    GPpdf.setTrainData(trainHisto, dataHisto)

    #get data from histogram to fit to
    data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(myy), dataHisto)

    #Fit the pdf to the data
    fitResult = GPpdf.fitTo(data, ROOT.RooFit.Save())
    #fitResult = GPpdf.fitTo(data)


    print "-------RooRealVar min at:", nSigGP.getVal(), nSigGP.getError()
    print "-------FitResult min at:", fitResult.floatParsFinal().Print("nSigGP")
    gphisto, sigSub, gpSigSub = GPpdf.getGPHisto(nSigGP.getVal())
    trainResult = GPpdf.trainResult
    #gphisto, sigSub, gpSigSub = GPpdf.getGPHisto(-19000)
    #gphisto125 = GPpdf.getGPHisto(-1300)


    # Get negative log likelihoods
    #nllGP = GPpdf.createNLL(data)
    #profileGP = nllGP.createProfile(ROOT.RooArgSet(nSigGP))

    # Plot the background fit and data
    c1 = ROOT.TCanvas('BkgFit','Background Fit')
    #frame = myy.frame(ROOT.RooFit.Range(120,130))
    frame = myy.frame()
    frame.SetTitle("GP Background, Kernel: {0}        SS: {1:.2f}".format( GPpdf.gp.kernel_ ,nSigGP.getVal() ))
    data.plotOn(frame)
    GPpdf.plotOn(frame, ROOT.RooFit.Normalization(1.0/data.binVolume()))
    frame.Draw()
    #gpConfig['priorMean'].Draw('same')
    #nSigGP.Print()
    #dataHisto.SetMarkerStyle(20)
    #dataHisto.SetLineColor(ROOT.kRed)
    #dataHisto.Draw('same')
    #gphisto.SetLineColor(ROOT.kRed)
    #gphisto.SetLineWidth(3)
    #gphisto100.Draw()
    #gphisto.Draw('same')

    #trainResult.SetLineColor(ROOT.kRed)
    #trainResult.SetLineWidth(3)
    #trainResult.Draw('same')
    #sigSub.SetLineWidth(3)
    #sigSub.SetLineColor(ROOT.kBlue)
    #sigSub.Draw('samehist')
    #gpSigSub.SetLineWidth(3)
    #gpSigSub.SetLineColor(ROOT.kGreen)
    #gpSigSub.Draw('samehist')
    #dataHisto.Draw('samep')
    c1.Print(args.outDir+'/test_GP_'+args.tag+'.pdf')

    """
    c2 = ROOT.TCanvas('nll', 'nll')
    nSigFrame = nSigGP.frame()
    nSigFrame.SetTitle("Min: "+str(nSigGP.getVal()))
    #nllGP.plotOn(nSigFrame,ROOT.RooFit.LineColor(ROOT.kRed))
    profileGP.plotOn(nSigFrame)
    nSigFrame.Draw()
    c2.Print(args.outDir+'/NLL_'+args.tag+'.pdf')
    """



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="../data/CouplingsTemplatesMoriond2017.root",help="Input root file with background histogram")
    parser.add_argument("--outDir", help="Name tag to add to begining of matplotlib plot")
    parser.add_argument("--doSig","-s", action='store_true', help="Inject signal into background and make some signal plots")
    parser.add_argument("--tag", help="Name tag to add to begining of matplotlib plot")

    args = parser.parse_args()
    #args.CL = '95'
    args.CL = '68'
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    #priorMeanFile = ROOT.TFile('../data/KDEOutput.root')
    #priorMean = priorMeanFile.Get('myy_high')

    gpConfig = {
                'lengthScale' : None
                #,'lengthScale_min': 1
                #,'lengthScale_max': 200
                #,'amplitude'      : 693**2
                ,'amplitude'      : None
                ,'train_range'    : [105,160]
    }
    args.sigMass = 125
    args.lowRange = 105
    args.hiRange = 160
    stats = 330000

    f = ROOT.TFile(args.input)

    #trainHisto = f.Get('myy_high')
    #trainHisto.Rebin(30)
    #trainHisto = f.Get('hmgg_c0')
    #trainHisto = f.Get('m_yy_fine')
    #trainHisto = f.Get('m_yy_c1_M17_ggH_0J_Cen_BkgTemplate') # plot 1,2 and 3
    trainHisto = f.Get('m_yy_c0_Inclusive_BkgTemplate') # plot 1,2 and 3
    #trainHisto = f.Get('m_yy_c23_M17_VHdilep_BkgTemplate') # plot 4
    #Normalize the trainig histo to the correct number of expected background events
    norm = 333000/ trainHisto.Integral() # plot 1
    #norm = 1100000/ trainHisto.Integral() # plot 2
    #norm = 10/ trainHisto.Integral() # plot 3 and 4
    trainHisto.Scale(norm)


    #gpConfig['lengthScale'] = float(1000)
    #dataHisto = ROOT.TH1F('dataHisto', 'Background data', 27, 105, 159)
    dataHisto = trainHisto.Clone()
    dataHisto.Reset()

    ROOT.gRandom.SetSeed(5)

    dataHisto.FillRandom(trainHisto, stats)

    priorMean = fitData(args, dataHisto, "pol1")
    #print "Function returns:",priorMean
    gpConfig['priorMean'] = priorMean


    #run(args, trainHisto, dataHisto)
    run(args, trainHisto, dataHisto, gpConfig)
