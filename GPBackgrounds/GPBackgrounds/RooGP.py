import ROOT
import numpy as np
from math import exp, sqrt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
#from RootToNp import *

ROOT.TH1.AddDirectory(ROOT.kFALSE)


def DSCB(myy, par):
    mu = par[0]
    alpha_low = 1.475   # alpha_low
    alpha_high = 1.902   # alpha_high
    n_low = 12.1   # n_low
    n_high = 11.6   # n_high
    sigma = 1.86  # sigma


    t = (myy[0] - mu)/sigma
    if ((-1*alpha_low) <= t) or (t < alpha_high):
        return exp((-1*t**2)/2)
    if (t < -1*alpha_low):
        num = exp(-0.5*alpha_low**2)
        A = alpha_low/n_low
        B = (n_low/alpha_low) - alpha_low - t
        return num/((A*B)**n_low)
    if (t > alpha_high):
        num = exp(-0.5*alpha_low**2)
        A = alpha_high/n_high
        B = (n_high/alpha_high) - alpha_high + t
        return num/((A*B)**n_high)
    else:
        raise("Inputs caused t to not fall in valid range")

def histoToArray(histo, xmin=False, xmax=False):
    """Take a histogram and convert the data to a numpy array for use in
    GaussianProcessRegressor

    \return arr,errarr  retuns a tuple containing the array of data points and the array of errors.
    """
    xaxis = np.array([])
    yaxis = np.array([])
    err = np.array([])
    minBin = 1
    maxBin = histo.GetNbinsX()+1
    if xmin and xmax:
        #print "Training range [{0},{1}]".format(str(xmin), str(xmax))
        minBin = histo.FindBin(xmin)
        maxBin = histo.FindBin(xmax)

    for nbin in range(minBin,maxBin):
        xaxis = np.append(xaxis, [histo.GetBinCenter(nbin)])
        yaxis = np.append(yaxis, [histo.GetBinContent(nbin)])
        err   = np.append(err, [histo.GetBinError(nbin)])
    return xaxis,yaxis,err

def arrayToHisto(title, xmin, xmax, arr, errarr=None):
    """Convert a numpy array into a histogram to plot with root.

    """
    #print "Length of arr:", len(arr)
    hist = ROOT.TH1F(title, title, len(arr),xmin,xmax)
    for nbin in range(1, len(arr)+1):
        hist.SetBinContent(nbin, arr[nbin-1]) # arr starts at 0 but histograms start at 1
        if False:
            hist.SetBinError(nbin, errarr[nbin-1]) # -1 because root starts at 1 but array index starts at 0

    return hist


class RooGP(ROOT.TPyRooGPSigPdf):
    def __init__(self, name, title, Myy, nSig, sigMass, trainHisto, dataHisto, trainRange=None):
        super(RooGP, self).__init__(self, name, title, Myy, nSig)
        self.name = name
        self.title = title
        self.Myy = Myy
        self.nSig = nSig
        self.sigMass = sigMass
        self.trainHisto = trainHisto
        self.dataHisto = dataHisto
        self.dataXmin = self.dataHisto.GetBinLowEdge(1)
        self.dataXmax = self.dataHisto.GetBinLowEdge(self.dataHisto.GetNbinsX())+self.dataHisto.GetBinWidth(1)
        if trainRange != None:
            self.trainMin = trainRange[0]  # The range to consider when training the GP
            self.trainMax = trainRange[1]
        else:
            self.trainMin = None  # The range to consider when training the GP
            self.trainMax = None 
        print "dataXmin: {0}   dataXmax: {1}".format(str(self.dataXmin), str(self.dataXmax))

        self.kernel = C((2.98e4)**2, (1e-3, 1e15)) * RBF(60, (1,200 )) #squared exponential kernel
        # Scale the training histogram to the datahisto stats. This will ensure that the amplitude
        # of the optimized kernel will correspond correctly to the y-range of the data.
        print "Scaling training histo to data:", dataHisto.Integral()/trainHisto.Integral()
        self.trainHisto.Scale(dataHisto.Integral()/trainHisto.Integral())
        self.trainHisto.Sumw2()
        self.opt_kernel = self.setTrainMC(self.trainHisto)
        #print "WARNING!!!!!!  Kernel is not being optimized. FOR TESTING ONLY. (RooGP.py)"
        #self.opt_kernel = C((2.1e3)**2) * RBF(40.8)

        #Need to think a bit more about how to set the range correctly here.
        #self.sigFunction = ROOT.TF1("dscb", DSCB, 105, 159, 1)
        self.sigFunction = ROOT.TF1("dscb", DSCB, self.dataXmin, self.dataXmax, 1)
        self.sigFunction.SetParameters(self.sigMass,0)
        self.sigNorm = self.sigFunction.Integral(self.dataXmin, self.dataXmax )
        print "Signal Normalization is", self.sigNorm

        self.currentNSig = None
        self.currentSBHist = None

        X, y, dataErrs = histoToArray(dataHisto, self.trainMin, self.trainMax)
        #GPh = GPHisto(dataHisto)
        #dataErrs = GPh.getErrArr()
        self.gp = GaussianProcessRegressor(kernel=self.opt_kernel,  #previously optimized kernel
                                            optimizer=None,         # Dont reoptimize hyperparamters.
                                            alpha=dataErrs**2)

    def setTrainMC(self, MCHisto):
        """Train a GP on a histogram to get hyperparamters.

        Use a high stats MC sample to optimize the hyperparamters then
        return a kernel object to be used in the data fit.
        """
        print "===== Optimizing hyperparamters on the training sample."

        #GPh = GPHisto(MCHisto)
        #X_t = GPh.getXArr()
        #Y_t = GPh.getYArr()
        #dy_t = GPh.getErrArr()

        X_t, Y_t, dy_t = histoToArray(MCHisto, self.trainMin, self.trainMax)

        gp = GaussianProcessRegressor(kernel=self.kernel
                                        ,alpha=dy_t**2
                                        ,n_restarts_optimizer=10)

        # Fit for the hyperparameters.
        gp.fit(np.atleast_2d(X_t).T, Y_t)
        print "Optimized hyperparameters:"
        print gp.kernel_

        # return a kernel object with hyperparameters optimized
        return gp.kernel_

    def getSigHisto(self, nSig):
        sig_func = self.sigFunction
        nSigScale = nSig*self.dataHisto.GetBinWidth(1) / self.sigNorm
        sigHisto = self.dataHisto.Clone("sigHisto")
        sigHisto.Reset()
        sigHisto.Add(sig_func, nSigScale)
        return sigHisto



    def getGPHisto(self, nSig):
        if self.currentNSig == nSig:
            return self.currentSBHist
        else:
            sig_func = self.sigFunction
            nSigScale = nSig*self.dataHisto.GetBinWidth(1) / self.sigNorm
            sigSubtractedHist = self.dataHisto.Clone("sigSubtractedHist")
            sigSubtractedHist.Add(sig_func, -1*nSigScale) #Subtrac the signal from the data
            #Get the Array corresponding to the subtracted histogram
            X, y_subtracted, errArr = histoToArray(sigSubtractedHist)
            #print "X array size", X.size
            # fit the GP to what is left. Then get a GP prediction on the result.
            self.gp.fit(np.atleast_2d(X).T,y_subtracted)
            y_pred, sigma = self.gp.predict(np.atleast_2d(X).T, return_std=True)

            gpHisto = arrayToHisto("GPhisto", self.dataXmin, self.dataXmax, y_pred)
            gpHisto.Add(sig_func, nSigScale) #Add the signal back into the data
            return gpHisto

    def evaluate(self, myy, nSig):
        """Overide the base class evaluate function.

        Get the signal distribution for the nSig that the PDF is trying and
        convert it to a histogram. Then subtract that histogram from the data
        histogram. Then fit the GP to the remaining histo. Return the resulting
        signal histo + the GP fit of the background for the myy bin in question.
        """
        #print "=== Myy:", myy, "   nSig:", nSig
        if not (nSig == self.currentNSig):
            self.currentNSig = nSig
            #print "==== New current nSig:", self.currentNSig
            sig_func = self.sigFunction
            nSigScale = nSig*self.dataHisto.GetBinWidth(1) / self.sigNorm
            #print "==== nSigScale:", nSigScale
            #sig_func.SetParameters (self.sigMass)
            #sig_Histo = sig_func.CreateHistogram()
            #subtract DSCB distribution from background
            sigSubtractedHist = self.dataHisto.Clone("sigSubtractedHist")
            sigSubtractedHist.Add(sig_func, -1*nSigScale) #Subtrac the signal from the data
            #Get the Array corresponding to the subtracted histogram
            X, y_subtracted, errArr = histoToArray(sigSubtractedHist, self.trainMin, self.trainMax)
            #print "X array size", X.size
            # fit the GP to what is left. Then get a GP prediction on the result.
            self.gp.fit(np.atleast_2d(X).T,y_subtracted)
            y_pred, sigma = self.gp.predict(np.atleast_2d(X).T, return_std=True)

            #gpHisto = arrayToHisto("GPhisto", self.dataXmin, self.dataXmax, y_pred)
            gpHisto = arrayToHisto("GPhisto", self.trainMin, self.trainMax, y_pred)
            gpHisto.Add(sig_func, nSigScale) #Add the signal back into the data
            if self.currentSBHist:
                self.currentSBHist.Delete()
            self.currentSBHist = gpHisto
            #print "=== evaluate return value:", gpHisto.GetBinContent(gpHisto.FindBin(myy))
            return gpHisto.GetBinContent(gpHisto.FindBin(myy))
        else:  #If we are still on the same nSig number just return the stored histogram.
            #print "=== evaluate return value:", self.currentSBHist.GetBinContent(self.currentSBHist.FindBin(myy))
            return self.currentSBHist.GetBinContent(self.currentSBHist.FindBin(myy))

    def analyticalIntegral(self, myyMin, myyMax):
        """Tell RooFit how to integrate the GP histo.

        """
        minBin = self.currentSBHist.FindBin(myyMin)
        maxBin = self.currentSBHist.FindBin(myyMax)
        result = self.currentSBHist.Integral(minBin, maxBin)
        return result
