import ROOT
import numpy as np
from math import exp
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
from RootToNp import *

def DSCB(myy, par):
    norm = par[0]
    mu = par[1]

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


class RooGP(ROOT.TPyRooGPSigPdf):
    def __init__(self, name, title, Myy, nSig, sigMass, trainHisto, dataHisto):
        super(RooGP, self).__init__(self, name, title, Myy, nSig)
        self.name = name
        self.title = title
        self.Myy = Myy
        self.nSig = nSig
        self.sigMass = sigMass
        self.trainHisto = trainHisto
        self.dataHisto = dataHisto

        self.kernel = C((2.98e4)**2, (1e-3, 1e15)) * RBF(60, (1,1e5 )) #squared exponential kernel
        # Scale the training histogram to the datahisto stats. This will ensure that the amplitude
        # of the optimized kernel will correspond correctly to the y-range of the data.
        trainHisto.Scale(dataHisto.Integral()/trainHisto.Integral())
        self.opt_kernel = self.setTrainMC(trainHisto)

        #Need to think a bit more about how to set the range correctly here.
        self.sigFunction = ROOT.TF1("dscb", DSCB, 105, 160, 2)

        self.currentNSig = None
        self.currentSBHist = None

        GPh = GPHisto(dataHisto)
        dataErrs = GPh.getErrArr()
        self.gp = GaussianProcessRegressor(kernel=self.opt_kernel,  #previously optimized kernel
                                            optimizer=None,         # Dont reoptimize hyperparamters.
                                            alpha=dataErrs**2)

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

    def evaluate(self, myy, nSig):
        """Overide the base class evaluate function.

        Get the signal distribution for the nSig that the PDF is trying and
        convert it to a histogram. Then subtract that histogram from the data
        histogram. Then fit the GP to the remaining histo. Return the resulting
        signal histo + the GP fit of the background for the myy bin in question.
        """
        if not (nSig == self.currentNSig):
            self.currentNSig = nSig
            print "==== New current nSig:", self.currentNSig
            #print "=======  Python - In Evaluate()."
            #print "Myy:", myy, "nSig", nSig
            #Get DSCB histogram
            sig_func = self.sigFunction
            sig_func.SetParameters(nSig, self.sigMass)
            #sig_Histo = sig_func.CreateHistogram()
            #subtract DSCB distribution from background
            sigSubtractedHist = self.dataHisto.Clone("sigSubtractedHist")
            sigSubtractedHist.Add(sig_func, -1) #Subtrac the signal from the data
            #Get the Array corresponding to the subtracted histogram
            X, y_subtracted, errArr = histoToArray(sigSubtractedHist)
            #print "X array size", X.size
            # fit the GP to what is left. Then get a GP prediction on the result.
            self.gp.fit(np.atleast_2d(X).T,y_subtracted)
            y_pred, sigma = self.gp.predict(np.atleast_2d(X).T, return_std=True)

            gpHisto = arrayToHisto("GPhisto", 105, 160, y_pred)
            gpHisto.Add(sig_func, 1) #Add the signal back into the data
            self.currentSBHist = gpHisto
            return gpHisto.GetBinContent(gpHisto.FindBin(myy))
        else:  #If we are still on the same nSig number just return the stored histogram.
            return self.currentSBHist.GetBinContent(self.currentSBHist.FindBin(myy))
