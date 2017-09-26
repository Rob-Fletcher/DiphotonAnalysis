from ROOT import *
import os
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF,Matern, ConstantKernel as C
from root_numpy import *
from matplotlib import pyplot as plt
import argparse
from RootToNp import *
from SignalModel import *

def run(args, mass, winLow, winHigh, seed, bkghist):


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

    kernel = C(10.0, (1e-3, 1e15)) * RBF(22.0, (1, 1e6)) #squared exponential kernel


    gp = GaussianProcessRegressor(kernel=kernel
                                    #,optimizer=None
                                    ,alpha=dy**2
                                    ,n_restarts_optimizer=15
                                    )

    gp.fit(X,y)
    y_pred, sigma = gp.predict(x, return_std=True)

    outhist = GPh.getHisto(y_pred, sigma, 'GP Fit')
    if args.noWindow:
        bkgWindow = GPh.getHisto(y, dy, 'Full Background')
    else:
        bkgWindow = GPh.getWinHisto(y, dy, 'Full Background')
    bkgSubtracted = bkghist.Clone('bkgSubtracted')
    bkgSubtracted.Add(outhist,-1)  #Subtract background prediction from background with injected signal.

    dscb_func = TF1("dscb", DSCB, 105, 160, 7)
    dscb_func.SetParameters(1  # Normalization
                        ,mass  # mu
                        ,1.475   # alpha_low
                        ,1.902   # alpha_high
                        ,12.1   # n_low
                        ,11.6   # n_high
                        ,1.86  ) # sigma

    #dscb_func.FixParameter(0,1)  #Normalization Dont want to fix this
    dscb_func.FixParameter(1,mass) #Mass Fixed to middle of window
    dscb_func.FixParameter(2, 1.475) #alpha_low
    dscb_func.FixParameter(3, 1.902) #alpha_high
    dscb_func.FixParameter(4, 12.1) #n_low
    dscb_func.FixParameter(5, 11.6) #n_high
    dscb_func.FixParameter(6, 1.68) # sigma
    bkgSubtracted.Fit("dscb","","", winLow, winHigh)
    fitResult = bkgSubtracted.GetFunction("dscb")
    norm = fitResult.GetParameter(0)
    ss = fitResult.Integral(winLow,winHigh)

    bkgSub = None
    fitRes = None
    GPHisto = None
    ToyHisto = None

    results = {
                "bkgSub"    : bkgSub,
                "fitResult" : fitRes
                "ss"        : ss,
                "norm"      : norm
    }
    if args.doGPPlot:
        results["GPHisto"] = outhist
        results["ToyHisto"] = bkgWindow

    return ss, norm

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
    bkghist_template.Rebin(8)
    stats = 60000



    ssHisto = TH1F("ss","ss",20, -200,200)

    index = 1

    args.windowSize = 10 #GeV
    for seed in range(40):
        print "Throwing Toy ({0},{1})".format(index, 40)
        bkghist = toyModel(bkghist_template, stats, seed)
        index+=1

        #loop over all windows in the range
        for mass in np.arange(114,150,2):
            winlow = mass - args.windowSize/2
            winHigh = mass + args.windowSize/2
            ss, norm = run(args, mass, winlow, winHigh, seed, bkghist)
            print "====== mass:", mass,"=========="
            print "SS is:", ss
            print "Norm is:", norm
            ssHisto.Fill(ss)

    canv = TCanvas('c','c')
    ssHisto.Draw()
    canv.Print(args.tag+'/spuriousSignal.pdf')
