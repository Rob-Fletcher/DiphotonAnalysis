import ROOT
import sys,os
import argparse
sys.path.append('../GPBackgrounds')
import RooGP
import RooGPBkg
from SignalModel import *
import json

ROOT.gROOT.SetBatch(ROOT.kTRUE)
ROOT.TH1.AddDirectory(ROOT.kFALSE)

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


def run(args, gpConfig):
    myy = ROOT.RooRealVar('myy','myy',args.lowRange,args.hiRange)
    nSigGP = ROOT.RooRealVar('nSigGP','nSigGP',-2000,2000)
    dataHisto = gpConfig['dataHisto']

    #The PDF for the GP bkg + signal DSCB
    GPpdf = RooGP.RooGP("mypdf", "CustomPDF",myy ,nSigGP, args.sigMass, gpConfig)
    #GPpdf.setTrainData(trainHisto, dataHisto)

    #get data from histogram to fit to
    data = ROOT.RooDataHist("dh", "dh", ROOT.RooArgList(myy), dataHisto)

    #Fit the pdf to the data
    fitResult = GPpdf.fitTo(data, ROOT.RooFit.SumW2Error(ROOT.kFALSE), ROOT.RooFit.Save())


    return nSigGP.getVal(), nSigGP.getError()

if __name__ == '__main__':
    """Run a spurious signal test for varying length scale.

    Make a plot of the maximum (maybe average would be better for low stats templates) SS as a function of length scale

    """

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default='../data/bkg_comb_Ioannis_h021.root',help="Input root file with background histogram")
    parser.add_argument("--outDir", help="Name tag to add to begining of matplotlib plot")
    parser.add_argument("--tag", help="Name tag to add to begining of matplotlib plot")

    args = parser.parse_args()
    #args.CL = '95'
    args.CL = '68'
    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    logFile = open(args.outDir+'/output.log', 'w')

    #priorMeanFile = ROOT.TFile('../data/KDEOutput.root')
    #priorMean = priorMeanFile.Get('myy_high')

    gpConfig = {
                'lengthScale' : 60
                #,'lengthScale_min': 1
                #,'lengthScale_max': 200
                #,'amplitude'      : 100
                ,'amplitude_min'      : 1e-3
                ,'amplitude_max'      : 1e15
                ,'train_range'    : [60,120]
                #,'priorMean' : priorMean
    }
    args.lowRange = 60
    args.hiRange = 120
    args.sigMass = 85
    nToys = 100
    stats = 333000
    f = ROOT.TFile(args.input)

    #trainHisto = f.Get('m_yy_c1_M17_ggH_0J_Cen_BkgTemplate') # plot 1,2 and 3
    #trainHisto = f.Get('m_yy_c0_Inclusive_BkgTemplate') # plot 1,2 and 3
    trainHisto = f.Get('hbkg_tmp_nominal_UU_incl')
    #trainHisto = f.Get('m_yy_c23_M17_VHdilep_BkgTemplate') # plot 4
    #trainHisto = f.Get('m_yy')
    #Normalize the trainig histo to the correct number of expected background events
    #norm = 333000/ trainHisto.Integral() # plot 1
    #norm = 1100000/ trainHisto.Integral() # plot 2
    #norm = 10/ trainHisto.Integral() # plot 3 and 4
    #trainHisto.Scale(norm)

    outfile = ROOT.TFile(args.outDir+"/output.root", "RECREATE")
    outfile.cd()

    ls = range(1,1001,50)
    totSSNom = ROOT.TH1F("SS_Nom","zero", len(ls), 0, 1000)
    totSSNom.SetDirectory(outfile)
    totSSMeanFunc = ROOT.TH1F("SS_MeanFunc","poly1", len(ls), 0, 1000)
    totSSMeanFunc.SetDirectory(outfile)
    totSSMeanFunc1 = ROOT.TH1F("SS_MeanFunc1","exp", len(ls), 0, 1000)
    totSSMeanFunc.SetDirectory(outfile)
    SSdata = {}
    SSdata['nom'] = []
    SSdata['mean'] = []
    SSdata['mean1'] = []

    for lengthScale in ls:

        gpConfig['lengthScale'] = float(lengthScale)
        #logFile.write("LengthScale: {}\n".format(str(lengthScale)))
        fitNSigNominal  = ROOT.TH1F("fitNSigNominal_len"+str(lengthScale), "Fitted nSig", 5000, -10000,10000)
        fitNSigMeanFunc  = ROOT.TH1F("fitNSigMeanFunc_len"+str(lengthScale), "Fitted nSig", 5000, -10000,10000)
        fitNSigMeanFunc1  = ROOT.TH1F("fitNSigMeanFunc1_len"+str(lengthScale), "Fitted nSig", 5000, -10000,10000)
        print "++++ Setting lengthScale to {}".format(gpConfig['lengthScale'])

        for i in range(nToys):
            #dataHisto = ROOT.TH1F('dataHisto', 'Background data', 27, 105, 159)
            print "======== Toy Seed: {}".format(i)
            trainHisto1 = trainHisto.Clone()

            dataHisto1 = trainHisto.Clone()
            dataHisto1.Reset()

            #dataHisto2 = trainHisto1.Clone()
            #dataHisto2.Reset()

            ROOT.gRandom.SetSeed(i)
            dataHisto1.FillRandom(trainHisto, stats)
            dataHisto2.FillRandom(trainHisto1, stats)

            gpConfig['amplitude'] = 6.05e+03**2
            gpConfig['trainHisto'] = trainHisto
            gpConfig['dataHisto'] = dataHisto1
            recoNsig, recoNsigError =    run(args, gpConfig)
            fitNSigNominal.Fill( recoNsig )

            gpConfig['amplitude'] = 992**2
            gpConfig['priorMean'] = fitData(args, dataHisto2, "pol1")
            gpConfig['trainHisto'] = trainHisto
            gpConfig['dataHisto'] = dataHisto1
            recoNsig, recoNsigError =    run(args, gpConfig)
            fitNSigMeanFunc.Fill( recoNsig )

            gpConfig['amplitude'] = 4.19e+03**2
            gpConfig['priorMean'] = fitData(args, dataHisto2, "expo")
            gpConfig['trainHisto'] = trainHisto
            gpConfig['dataHisto'] = dataHisto1
            recoNsig, recoNsigError =    run(args, gpConfig)
            fitNSigMeanFunc1.Fill( recoNsig )
            # end Toys loop

        meanNom = fitNSigNominal.GetMean()
        errNom = fitNSigNominal.GetRMS()
        meanMeanFunc = fitNSigMeanFunc.GetMean()
        errMeanFunc = fitNSigMeanFunc.GetRMS()
        meanMeanFunc1 = fitNSigMeanFunc1.GetMean()
        errMeanFunc1 = fitNSigMeanFunc1.GetRMS()
        #logFile.write("   Nominal   mean: {}   err: {}\n".format(str(meanNom), str(errNom)))
        #logFile.write("   MeanFunc  mean: {}   err: {}\n".format(str(meanMeanFunc), str(errMeanFunc)))
        totSSNom.SetBinContent(totSSNom.FindBin(lengthScale), meanNom)
        totSSNom.SetBinError(totSSNom.FindBin(lengthScale), errNom)
        totSSMeanFunc.SetBinContent(totSSMeanFunc.FindBin(lengthScale), meanMeanFunc)
        totSSMeanFunc.SetBinError(totSSMeanFunc.FindBin(lengthScale), errMeanFunc)
        totSSMeanFunc1.SetBinContent(totSSMeanFunc1.FindBin(lengthScale), meanMeanFunc1)
        totSSMeanFunc1.SetBinError(totSSMeanFunc1.FindBin(lengthScale), errMeanFunc1)
        SSdata['nom'].append( (lengthScale, meanNom) )
        SSdata['mean'].append( (lengthScale, meanMeanFunc) )
        SSdata['mean1'].append( (lengthScale, meanMeanFunc1) )

    c = ROOT.TCanvas('fitNsig', 'Spurious Signal')
    totSSNom.SetMarkerStyle(20)
    totSSNom.SetMarkerColor(ROOT.kBlack)
    totSSNom.SetStats(0)
    totSSNom.Draw("PCE")
    totSSMeanFunc.SetMarkerStyle(20)
    totSSMeanFunc.SetMarkerColor(ROOT.kRed)
    totSSMeanFunc.Draw("PCE same")
    totSSMeanFunc1.SetMarkerStyle(20)
    totSSMeanFunc1.SetMarkerColor(ROOT.kBlue)
    totSSMeanFunc1.Draw("PCE same")
    c.BuildLegend()
    c.Print(args.outDir+'/totSS_'+args.tag+'.pdf')

    totSSNom.Write()
    totSSMeanFunc.Write()
    totSSMeanFunc1.Write()
    json.dump(SSdata, logFile)
    outfile.Close()
    logFile.close()
    f.Close()
