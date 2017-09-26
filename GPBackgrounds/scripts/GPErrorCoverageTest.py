from ROOT import *
import os
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF,Matern, ConstantKernel as C
#from root_numpy import *
#from matplotlib import pyplot as plt
import argparse
from RootToNp import *
from SignalModel import *

def run(args, mass, winLow, winHigh, bkghist, template):


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

    kernel = C(10.0, (1e-3, 1e15)) * RBF(22.0, (np.sqrt(2)*(8**2), 1e6)) #squared exponential kernel


    gp = GaussianProcessRegressor(kernel=kernel
                                    #,optimizer=None
                                    ,alpha=dy**2
                                    ,n_restarts_optimizer=15
                                    )

    gp.fit(X,y)
    y_pred, sigma = gp.predict(x, return_std=True)

    outhist = GPh.getHisto(y_pred, 1.96*sigma, 'GP Fit')
    print "outhist nbins:", outhist.GetNbinsX()
    if args.noWindow:
        bkgWindow = GPh.getHisto(y, dy, 'Full Background')
    else:
        bkgWindow = GPh.getWinHisto(y, dy, 'Full Background')

    errCov = outhist.Clone('errorCoverage')
    errCov.Reset()

    norm = (bkghist.Integral()/template.Integral())

    for nbin in range(1,outhist.GetNbinsX()):
        print nbin
        Err = outhist.GetBinError(nbin)
        cont = outhist.GetBinContent(nbin)
        binCenter = outhist.GetBinCenter(nbin)
        temp = norm*template.GetBinContent(template.FindBin(binCenter))
        if temp > (cont + Err) or temp < (cont - Err):
            print "found uncovered bin: ", nbin
            print "LowErr : template : HighErr ----  {0} : {1} : {2}".format((cont-Err), temp,(cont+Err))
            print "----------"
            errCov.Fill(binCenter)

    return errCov

    pass #run



if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument("--input", help="Input root file with background histogram")
    parser.add_argument("--tag", help="Name tag to add to begining of matplotlib plot")
    parser.add_argument("--mplot","-m", action='store_true', help="Plot with matplotlib")
    parser.add_argument("--doSig","-s", action='store_true', help="Inject signal into background and make some signal plots")
    parser.add_argument("--noWindow","-w", action='store_true', help="Dont do the window cut. Train on whole sample.")

    args = parser.parse_args()

    if not os.path.exists(args.tag):
        os.makedirs(args.tag)

    f = TFile(args.input)

    bkghist_template = f.Get('hmgg_c0')
    template = f.Get('hmgg_c0')
    bkghist_template.Rebin(8)
    stats = 60000



    errCov = None


    args.windowSize = 10 #GeV
    for seed in range(1,100):
        print "Throwing Toy ({0},{1})".format(seed, 40)
        bkghist = toyModel(bkghist_template, stats, seed)
        print "Bins in bkghist:", bkghist.GetNbinsX()

        args.noWindow = True

        if not errCov:
            errCov = run(args, 125,120,130, bkghist, template)
        else:
            errCov.Add(run(args, 125,120,130, bkghist, template))


    #errCov.Scale(100/40)
    canv = TCanvas('c','c')
    errCov.Draw()
    canv.Print(args.tag+'/ErrorCoverage.pdf')
