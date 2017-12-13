import ROOT
import sys,os
import argparse
sys.path.append('../GPBackgrounds')
import RooGP
import RooGPBkg
from SignalModel import *

ROOT.gROOT.SetBatch(ROOT.kTRUE)

def run(args, trainHisto):
    myy = ROOT.RooRealVar('myy','myy',105,160)
    nSigGP = ROOT.RooRealVar('nSigGP','nSigGP',-1000,1000)

    #The PDF for the GP bkg + signal DSCB
    GPpdf = RooGP.RooGP("mypdf", "CustomPDF",myy ,nSigGP, args.sigMass, trainHisto, trainHisto)

    #get data from histogram to fit to
    data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(myy), trainHisto)

    #Fit the pdf to the data
    fitResult = GPpdf.fitTo(data, ROOT.RooFit.SumW2Error(ROOT.kTRUE), ROOT.RooFit.Save())


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
    parser.add_argument("--outDir", help="Name tag to add to begining of matplotlib plot")
    parser.add_argument("--tag", help="Name tag to add to begining of matplotlib plot")

    args = parser.parse_args()
    #args.CL = '95'
    args.CL = '68'
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)


    args.lowMass = 110
    args.hiMass  = 158
    f = ROOT.TFile(args.input)

    trainHisto = f.Get('m_yy_c1_M17_ggH_0J_Cen_BkgTemplate')
    cont_unscaled = trainHisto.GetBinContent(trainHisto.FindBin(125))
    err_unscaled  = trainHisto.GetBinError(trainHisto.FindBin(125))
    print "Unscaled - Content: {}    Error: {}   Relative Error: {}".format(str(cont_unscaled), str(err_unscaled), str(err_unscaled/cont_unscaled))
    #Normalize the trainig histo to the correct number of expected background events
    norm = 333000/ trainHisto.Integral()
    #norm = 1100000/ trainHisto.Integral()
    #norm = 333000/ trainHisto.Integral()
    trainHisto.Scale(norm)
    cont_scaled = trainHisto.GetBinContent(trainHisto.FindBin(125))
    err_scaled  = trainHisto.GetBinError(trainHisto.FindBin(125))
    print "Scaled - Content: {}    Error: {}   Relative Error: {}".format(str(cont_scaled), str(err_scaled), str(err_scaled/cont_scaled))

    outfile = ROOT.TFile("output.root", "RECREATE")


    SSHisto = ROOT.TH1F("nSS", "N Spurious Signal",158-110, 110, 158 )
    SSHisto.SetDirectory(outfile)
    SSRelHisto = ROOT.TH1F("ssS","SS / S", 158-110, 110,158)
    SSRelHisto.SetDirectory(outfile)
    SSErrHisto = ROOT.TH1F("ssErr", "SS / Err", 158-110, 110, 158)
    SSErrHisto.SetDirectory(outfile)

    #SSHisto = trainHisto.Clone('dataHisto')
    #SSHisto.Reset()
    #SSHisto.Rebin(10)
    #SSRelHisto = SSHisto.Clone('SSRelHisto')
    #SSErrHisto = SSHisto.Clone('SSErrHisto')
    loopN = 1
    ssMax = 0
    for sigMass in range(110, 158):

        print "====================================================="
        print "({0}/{1})Getting SS for signal mass: {2}".format(loopN, 158-110, sigMass)
        args.sigMass = sigMass
        recoNsig, recoNsigError =    run(args, trainHisto)
        SSHisto.SetBinContent(SSHisto.FindBin(sigMass), recoNsig )
        SSHisto.SetBinError(SSHisto.FindBin(sigMass), recoNsigError)
        SSRelHisto.SetBinContent(SSRelHisto.FindBin(sigMass), recoNsig/335)
        SSErrHisto.SetBinContent(SSErrHisto.FindBin(sigMass), recoNsig/recoNsigError)
        print "----  Found SS:", recoNsig, "+-", recoNsigError
        loopN = loopN + 1


    c = ROOT.TCanvas('fitNsig', 'Fit nSig')
    SSHisto.SetMarkerStyle(20)
    SSHisto.Draw("PC")
    c.Print(args.outDir+'/firNsig_'+args.tag+'.pdf')

    c_ssRel = ROOT.TCanvas('SSRel', 'SS / S')
    SSRelHisto.SetMarkerStyle(20)
    SSRelHisto.Draw("PC")
    c_ssRel.Print(args.outDir+'/ssRel_'+args.tag+'.pdf')

    c_ssErr = ROOT.TCanvas()
    SSErrHisto.SetMarkerStyle(20)
    SSErrHisto.GetYaxis().SetRangeUser(-2,2)
    SSErrHisto.Draw("PC")
    c_ssErr.Print(args.outDir+'/ssErr_'+args.tag+'.pdf')

    print "Max SS:", SSHisto.GetMaximum()

    outfile.Write()
    outfile.Close()
    f.Close()
