import ROOT
import numpy as np
from math import exp
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
from RootToNp import *

def DSCB(myy, par):
    norm = par[0]
    mu = par[1]

    # Need to set all these parameters from the paper
    #alpha_low = par[2]
    #alpha_high = par[3]
    #n_low = par[4]
    #n_high = par[5]
    #sigma = par[6]
    alpha_low = 1.475   # alpha_low
    alpha_high = 1.902   # alpha_high
    n_low = 12.1   # n_low
    n_high = 11.6   # n_high
    sigma = 1.86  # sigma


    t = (myy[0] - mu)/sigma
    if ((-1*alpha_low) <= t) or (t < alpha_high):
        return norm * exp((-1*t**2)/2)
    if (t < -1*alpha_low):
        num = exp(-0.5*alpha_low**2)
        A = alpha_low/n_low
        B = (n_low/alpha_low) - alpha_low - t
        return norm * num/((A*B)**n_low)
    if (t > alpha_high):
        num = exp(-0.5*alpha_low**2)
        A = alpha_high/n_high
        B = (n_high/alpha_high) - alpha_high + t
        return norm * num/((A*B)**n_high)
    else:
        raise("Inputs caused t to not fall in valid range")


class RooGPBkg(ROOT.TPyRooGPPdf):
    def __init__(self, name, title, Myy, trainHisto, dataHisto):
        super(RooGPBkg, self).__init__(self, name, title, Myy)
        self.name = name
        self.title = title
        self.Myy = Myy
        self.trainHisto = trainHisto
        self.dataHisto = dataHisto

        self.kernel = C((2.98e4)**2, (1e-3, 1e15)) * RBF(60, (1,1e5 )) #squared exponential kernel
        trainHisto.Scale(dataHisto.Integral()/trainHisto.Integral())
        self.opt_kernel = self.setTrainMC(trainHisto)

        #Need to think a bit more about how to set the range correctly here.
        self.sigFunction = ROOT.TF1("dscb", DSCB, 105, 160, 2)

        self.currentNSig = None
        self.currentSBHist = None

        GPh = GPHisto(dataHisto)
        X = GPh.getXArr()
        Y = GPh.getYArr()
        dataErrs = GPh.getErrArr()
        self.gp = GaussianProcessRegressor(kernel=self.opt_kernel,  #previously optimized kernel
                                            optimizer=None,         # Dont reoptimize hyperparamters.
                                            alpha=dataErrs**2)

        self.gp.fit(X,Y)
        y_pred, sigma = self.gp.predict(X, return_std=True)
        #self.gpHisto = arrayToHisto("GPhisto", 105, 160, y_pred)
        self.gpHisto = GPh.getHisto(y_pred, sigma, "GPhisto")

    def setTrainMC(self, MCHisto):
        """Train a GP on a histogram to get hyperparamters.

        Use a high stats MC sample to optimize the hyperparamters then
        return a kernel object to be used in the data fit.
        """
        print "===== Optimizing hyperparamters on the training sample."
        GPh = GPHisto(MCHisto)
        X_t = GPh.getXArr()
        Y_t = GPh.getYArr()
        dy_t = GPh.getErrArr()

        gp = GaussianProcessRegressor(kernel=self.kernel
                                        ,alpha=dy_t**2
                                        ,n_restarts_optimizer=10)

        # Fit for the hyperparameters.
        gp.fit(X_t, Y_t)
        print "Optimized hyperparameters:"
        print gp.kernel_

        # return a kernel object with hyperparameters optimized
        return gp.kernel_

    def evaluate(self, myy):
        """Overide the base class evaluate function.

        """
        return self.gpHisto.GetBinContent(self.gpHisto.FindBin(myy))
