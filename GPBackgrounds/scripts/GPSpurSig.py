import ROOT
import sys,os
import argparse
sys.path.append('../GPBackgrounds')
import RooGP
import RooGPBkg
from SignalModel import *

ROOT.gROOT.SetBatch(ROOT.kTRUE)

def run(args, gpConfig, bkgOnly=False):
    myy = ROOT.RooRealVar('myy','myy',args.lowRange,args.hiRange)
    nSigGP = ROOT.RooRealVar('nSigGP','nSigGP',-2000,2000)

    #The PDF for the GP bkg + signal DSCB
    GPpdf = RooGP.RooGP("mypdf", "CustomPDF",myy ,nSigGP, args.sigMass, gpConfig)

    #get data from histogram to fit to
    data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(myy), gpConfig['trainHisto'])

    #Fit the pdf to the data
    fitResult = GPpdf.fitTo(data, ROOT.RooFit.SumW2Error(ROOT.kTRUE), ROOT.RooFit.Save())

    if bkgOnly:
        return GPpdf.getGPHisto(0)


    return nSigGP.getVal(), nSigGP.getError()

if __name__ == '__main__':
    """
    You can find the background templates, provided by Kurt, in the attached
    file CouplingsTemplatesMoriond2017.root. Use the histogram
    m_yy_c1_M17_ggH_0J_Cen_BkgTemplate for the test 1.,  2 and 3..
    Use m_yy_c23_M17_VHdilep_BkgTemplate for test 4. You have to normalize the
    histogram to the number of expected background event: I suggest to use 333k
    for 1., 1.1M for 2 and 10 for 3. and 4.

    For the value of S to be used to normalize the SS in the first plot I suggest
    to use 335 for 1., 1115 for 2. and 0.9 for 3 and 4.

    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--input",default='../data/CouplingsTemplatesMoriond2017.root', help="Input root file with background histogram")
    parser.add_argument("--outDir", help="Directory to write the results out to")
    parser.add_argument("--tag", help="Name tag to add to begining of plot")

    args = parser.parse_args()
    #args.CL = '95'
    args.CL = '68'
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)


    args.lowRange = 55
    args.hiRange = 120
    args.lowMass = 65
    args.hiMass  = 110
    S = 335 # plot 1
    #S = 1115 # plot 2
    #S = 0.9 # plot 3 and 4
    f = ROOT.TFile(args.input)

    trainHisto = f.Get('m_yy_conv0_fine_98perW') # plot 1,2 and 3
    #trainHisto = f.Get('m_yy_c1_M17_ggH_0J_Cen_BkgTemplate') # plot 1,2 and 3
    #trainHisto = f.Get('m_yy_c23_M17_VHdilep_BkgTemplate') # plot 4
    #trainHisto = f.Get('m_yy')
    #Normalize the trainig histo to the correct number of expected background events
    #norm = 333000/ trainHisto.Integral() # plot 1
    #norm = 1100000/ trainHisto.Integral() # plot 2
    #norm = 10/ trainHisto.Integral() # plot 3 and 4
    norm = 79.801/1.484
    trainHisto.Scale(norm)

    outfile = ROOT.TFile(args.outDir+"/output.root", "RECREATE")


    SSHisto = ROOT.TH1F("nSS", "N Spurious Signal",args.hiMass-args.lowMass, args.lowMass, args.hiMass )
    SSHisto.SetDirectory(outfile)
    SSRelHisto = ROOT.TH1F("ssS","SS / S", 110-65, 65,110)
    SSRelHisto.SetDirectory(outfile)
    SSErrHisto = ROOT.TH1F("ssErr", "SS / Err", 110-65, 65, 110)
    SSErrHisto.SetDirectory(outfile)

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
    #SSHisto = trainHisto.Clone('dataHisto')
    #SSHisto.Reset()
    #SSHisto.Rebin(10)
    #SSRelHisto = SSHisto.Clone('SSRelHisto')
    #SSErrHisto = SSHisto.Clone('SSErrHisto')
    loopN = 1
    ssMax = 0
    for sigMass in range(args.lowMass, args.hiMass):

        print "====================================================="
        print "({0}/{1})Getting SS for signal mass: {2}".format(loopN, args.hiMass-args.lowMass, sigMass)
        args.sigMass = sigMass
        recoNsig, recoNsigError =    run(args, gpConfig)
        SSHisto.SetBinContent(SSHisto.FindBin(sigMass), recoNsig )
        SSHisto.SetBinError(SSHisto.FindBin(sigMass), recoNsigError)
        SSRelHisto.SetBinContent(SSRelHisto.FindBin(sigMass), recoNsig/S)
        SSErrHisto.SetBinContent(SSErrHisto.FindBin(sigMass), recoNsig/recoNsigError)
        print "----  Found SS:", recoNsig, "+-", recoNsigError
        ssMax = max(ssMax, abs(recoNsig))
        loopN = loopN + 1

    cBkgOnly = ROOT.TCanvas('bkgOnly', 'Bkg Only Fit')
    bkgOnlyHisto = run(args, gpConfig, True)[0]
    trainHisto.SetStats(0)
    trainHisto.SetMarkerStyle(20)
    trainHisto.Draw()
    trainHisto.GetXaxis().SetRangeUser(55, 120)
    bkgOnlyHisto.Draw("samehist")
    cBkgOnly.Print(args.outDir+'/bkgOnly_'+args.tag+'.pdf')

    c = ROOT.TCanvas('fitNsig', 'Fit nSig')
    SSHisto.SetMarkerStyle(20)
    SSHisto.SetStats(0)
    SSHisto.Draw("PCE")
    c.Print(args.outDir+'/fitNsig_'+args.tag+'.pdf')

    c_ssRel = ROOT.TCanvas('SSRel', 'SS / S')
    SSRelHisto.SetMarkerStyle(20)
    SSRelHisto.SetStats(0)
    SSRelHisto.Draw("PC")
    c_ssRel.Print(args.outDir+'/ssRel_'+args.tag+'.pdf')

    c_ssErr = ROOT.TCanvas()
    SSErrHisto.SetMarkerStyle(20)
    #SSErrHisto.GetYaxis().SetRangeUser(-2,2)
    SSErrHisto.SetStats(0)
    SSErrHisto.Draw("PC")
    c_ssErr.Print(args.outDir+'/ssErr_'+args.tag+'.pdf')

    print "Max SS:", ssMax #SSHisto.GetMaximum()

    outfile.Write()
    outfile.Close()
    f.Close()
