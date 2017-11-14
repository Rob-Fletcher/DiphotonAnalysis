from ROOT import *
import os
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF,Matern, ConstantKernel as C
#from root_numpy import *
#from matplotlib import pyplot as plt
import argparse
import sys
sys.path.append('../GPBackgrounds')
import RooGP
import RooGPBkg
from RootToNp import *
from SignalModel import *
from MYYKernel import *
import re

gROOT.SetBatch(True)

def run(args, bkghist, trainHisto, optKernel):


    GPh = GPHisto(bkghist)

    # The distributions with no window removed.
    X_t = GPh.getXArr()
    y_t = GPh.getYArr()
    dy_t = GPh.getErrArr()
    if args.noWindow:
        X = X_t
        y = y_t
        dy = dy_t


    x = GPh.getXArr()
    gp = None

    kernel = optKernel
    gp = GaussianProcessRegressor(kernel=kernel
                                    ,optimizer=None
                                    ,alpha=dy**2
                                    )


    gp.fit(X,y)
    print gp.kernel_

    length = float(re.search('length_scale=(\d+(\.\d+)?)', gp.kernel_.__repr__()).group(1))

    y_pred, sigma = gp.predict(x, return_std=True)

    outhist = GPh.getHisto(y_pred, sigma, 'GP Fit')
    #bkg = GPh.getHisto(y, dy, 'Full Background')


    ###  RooFit part
    myy = RooRealVar('myy','myy',105,160)
    #nSig = RooRealVar('nSig','nSig',-200,1000)
    sigMass = 125
    #pdf = RooGP.RooGP("mypdf", "CustomPDF",myy ,nSig, sigMass, trainHisto, dataHisto)
    pdf = RooGPBkg.RooGPBkg("bkgPDF", "bkg only PDF",myy, trainHisto, bkghist)

    data = RooDataHist("dh", "dh", RooArgList(myy), bkghist)


    c1 = TCanvas('c1','c1')
    frame = myy.frame()
    data.plotOn(frame, RooFit.MarkerColor(kRed))
    fitResult = pdf.fitTo(data, RooFit.Save())
#    pdf.gpHisto.Draw()
#    data.Draw('same')
#    outhist.Draw('samehist')
#    outhist.SetFillColorAlpha(kWhite, 0)
#    bkghist.Draw('same')
#    bkghist.SetMarkerColor(kBlack)
#    outhist.Divide(bkghist)
#    outhist.Draw()
    pdf.plotOn(frame)
    #fitResult.plotOn(frame)
    frame.Draw()
    #nSig.Print()
    c1.Print(args.outDir+'/test_GP.pdf')


    pass  #Run

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input root file with background histogram")
    parser.add_argument("--outDir", help="Name outDir to add to begining of matplotlib plot")
    parser.add_argument("--doSig","-s", action='store_true', help="Inject signal into background and make some signal plots")

    args = parser.parse_args()

    args.noWindow = True

    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    f = TFile(args.input)

    stats = 10000000
    seed = 1

    bkghist_template = f.Get('hmgg_c0')
    bkghist_template.Rebin(8)
    bkghist_template.Scale(stats/bkghist_template.Integral())

    trainHisto = bkghist_template.Clone('trainHisto')

    optKernel = None

    GPh = GPHisto(bkghist_template)

    X = GPh.getXArr()
    y = GPh.getYArr()
    dy = GPh.getErrArr()

    x = GPh.getXArr()

    kernel = C(10.0, (1e-3, 1e15)) * RBF(50.0, (1e-2, 1e5)) #squared exponential kernel
    #kernel = C(1000.0, (1e-3, 1e15)) * FallExp(1.0, (1e-5, 1e2), 1.0, (1e-3,1e15)) * Gibbs(0.01, (1e-7, 1e5), 0.01, (1e-7,1e5))


    gp = GaussianProcessRegressor(kernel=kernel
                                    #,optimizer='fmin'
                                    ,alpha=dy**2
                                    ,n_restarts_optimizer=15
                                    )

    gp.fit(X,y)
    print gp.kernel_
    lengthOpt = float(re.search('length_scale=(\d+(\.\d+)?)', gp.kernel_.__repr__()).group(1))
    #raise()
    optKernel = gp.kernel_

    bkghist = toyModel(bkghist_template, stats, seed)

    run(args, bkghist, trainHisto, optKernel)
