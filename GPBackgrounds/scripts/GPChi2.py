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
from RootToNp import *
from SignalModel import *
from MYYKernel import *
import re

gROOT.SetBatch(True)

def run(args, mass, winLow, winHigh, bkghist, template, optKernel):


    GPh = GPHisto(bkghist)
    GPh.setWindow(winLow,winHigh)
    X = GPh.getXWindowArr()
    y = GPh.getYWindowArr()
    dy = GPh.getErrWindowArr()

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

    if optKernel is None:
        kernel = C(10.0, (1e-3, 1e15)) * RBF(22.0, (1e-2, 1e6)) #squared exponential kernel
        #kernel = C(1000.0, (1e-3, 1e15)) * FallExp(1.0, (1e-5, 1e2), 1.0, (1e-3,1e15)) * Gibbs(1.0, (1e-3, 1e5), 1.0, (1e-3,1e5))
        gp = GaussianProcessRegressor(kernel=kernel
                                        #,optimizer='fmin'
                                        ,alpha=dy**2
                                        ,n_restarts_optimizer=10
                                        )
    else:
        kernel = optKernel
        gp = GaussianProcessRegressor(kernel=kernel
                                        ,optimizer=None
                                        ,alpha=dy**2
                                        ,n_restarts_optimizer=10
                                        )


    gp.fit(X,y)
    print gp.kernel_

    length = float(re.search('length_scale=(\d+(\.\d+)?)', gp.kernel_.__repr__()).group(1))

    y_pred, sigma = gp.predict(x, return_std=True)

    outhist = GPh.getHisto(y_pred, sigma, 'GP Fit')
    bkgWindow = GPh.getHisto(y, dy, 'Full Background')

    ####### Poly2 fit
    #canv4 = TCanvas('c4','c4')
    expPol_func = TF1("expPol","[0]*exp((x-100)/100 * ([1] + [2]*(x-100)/100))",105,160)
    expPol_func.SetParameters(0,0,0)
    expPol_func.SetParLimits(1,-10.,10.)
    expPol_func.SetParLimits(2,-10.,10.)
    bkgWindow.Fit("expPol","","",105,160)
    expFitResult = bkgWindow.GetFunction("expPol")
    #expPolHist = expFitResult.GetHistogram()


    #expChi2 = expFitResult.GetChisquare()
    #GPChi2  = bkgWindow.Chi2Test(outhist,"CHI2")

    expChi2 = 0
    GPChi2  = 0

    for nbin in range(1,bkgWindow.GetNbinsX()):
        binCenter = bkgWindow.GetBinCenter(nbin)
        binContent = bkgWindow.GetBinContent(nbin)
        fitval = expFitResult.Eval(binCenter)
        GPVal  = outhist.GetBinContent(nbin)
        expChi2 += (binContent - fitval)**2 / fitval
        GPChi2 += (binContent - GPVal)**2 / GPVal





    print "expChi2:", expChi2
    print "GPChi2:", GPChi2
    print "length:", length


    return expChi2, GPChi2, length

    pass #run



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input root file with background histogram")
    parser.add_argument("--outDir", help="Name outDir to add to begining of matplotlib plot")
    parser.add_argument("--mplot","-m", action='store_true', help="Plot with matplotlib")
    parser.add_argument("--doSig","-s", action='store_true', help="Inject signal into background and make some signal plots")
    parser.add_argument("--noWindow","-w", action='store_true', help="Dont do the window cut. Train on whole sample.")

    args = parser.parse_args()

    args.noWindow = True

    if not os.path.exists(args.outDir):
        os.makedirs(args.outDir)

    f = TFile(args.input)

    stats = 10000000
    nToys = 1000

    bkghist_template = f.Get('hmgg_c0')
    template = f.Get('hmgg_c0')
    bkghist_template.Rebin(8)
    bkghist_template.Scale(stats/bkghist_template.Integral())

    optKernel = None

    GPh = GPHisto(bkghist_template)
    #GPh.setWindow(120,130)
    #X = GPh.getXWindowArr()
    #y = GPh.getYWindowArr()
    #dy = GPh.getErrWindowArr()

    # The distributions with no window removed.
    X = GPh.getXArr()
    y = GPh.getYArr()
    dy = GPh.getErrArr()
    #if args.noWindow:
    #    X = X_t
    #    y = y_t
    #    dy = dy_t

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
    #optKernel = gp.kernel_

    h_expChi2 = TH1F("expChi2","expChi2",15,10,40)
    h_GPChi2  = TH1F("GPChi2","GPChi2",15,10,40)
    h_hyperparams = TH1F("length", "length", 40, 20, 100)


    args.windowSize = 10 #GeV
    for seed in range(1,nToys):
        print "Throwing Toy ({0},{1})".format(seed, nToys)
        bkghist = toyModel(bkghist_template, stats, seed)



        expChi2, GPChi2, length = run(args, 125,120,130, bkghist, template, optKernel)

        h_expChi2.Fill(expChi2)
        h_GPChi2.Fill(GPChi2)
        h_hyperparams.Fill(length)



    canv = TCanvas('c','c')
    h_expChi2.SetLineColor(kRed)
    h_GPChi2.SetLineColor(kBlack)
    h_expChi2.Draw('hist')
    h_GPChi2.Draw('samehist')
    canv.Print(args.outDir+'/Chi2Test.pdf')

    canv1 = TCanvas('c1', 'c1')
    h_hyperparams.Draw('hist')
    canv1.Update()
    ymin = canv1.GetFrame().GetY1()
    ymax = canv1.GetFrame().GetY2()
    line = TLine(lengthOpt,ymin,lengthOpt,ymax)
    line.SetLineColor(kRed)
    line.SetLineWidth(3)
    line.Draw("same")
    canv1.Print(args.outDir+'/lengthDistribution.pdf')
    #print gp.kernel_
