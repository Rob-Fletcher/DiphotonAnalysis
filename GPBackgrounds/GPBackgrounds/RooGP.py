import ROOT
import numpy as np
from math import exp, sqrt
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C

ROOT.TH1.AddDirectory(ROOT.kFALSE)


def DSCB(myy, par):
    """Signal double sided crystal ball function

    """
    mass = par[0]

    mCBSignal0 =  0.048803
    mCBSignal1 = -0.00512848
    sCBSignal0 =  1.3612
    sCBSignal1 =  0.78698
    aLoSignal0 =  1.6738
    aLoSignal1 = -0.0159122
    aHiSignal0 =  1.4824
    aHiSignal1 =  0.12803
    nLoSignal0 =  9.5595
    nHiSignal0 =  72085

    mnX = (mass - 100)/100
    sigma  = (sCBSignal0   + sCBSignal1*mnX)
    deltaMSignal  = (mCBSignal0 + mCBSignal1*mnX)
    alpha_low = (aLoSignal0 + aLoSignal1*mnX)
    n_low = (nLoSignal0)
    alpha_high = (aHiSignal0 + aHiSignal1*mnX)
    n_high = (nHiSignal0)

    Mx  = (mass + deltaMSignal)

    #alpha_low = 1.475   # alpha_low
    #alpha_high = 1.902   # alpha_high
    #n_low = 12.1   # n_low
    #n_high = 11.6   # n_high
    #sigma = 1.86  # sigma

    """
    MnX = (Mx-100.0)/100.0
    deltaMx = -0.322 - 0.14 * MnX
    sigma = 1.418 + 1.03 * MnX
    alpha_low = 1.69 - 0.02 * MnX
    n_low = 13
    alpha_hi = 2.24 + 0 * MnX
    n_high = 7.9
    """

    t = (myy[0] - Mx)/sigma
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

################################################################################
#           Helper Functions   ROOT / numpy conversions
################################################################################

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
        #err   = np.append(err, [sqrt( abs(histo.GetBinContent(nbin)) )])
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
################################################################################


class RooGP(ROOT.TPyRooGPSigPdf):
    """Roofit implememntation of the Gaussian Process

    """
    def __init__(self, name, title, Myy, nSig, sigMass, gpConfig=None):
        super(RooGP, self).__init__(self, name, title, Myy, nSig)
        self.name = name
        self.title = title
        self.Myy = Myy
        self.nSig = nSig
        self.sigMass = sigMass
        self.trainResult = None

        # Setup the prior Mean function if passed through config
        self.priorMean = None
        if 'priorMean' in gpConfig:
            print "Using prior Mean function from config."
            self.priorMean = gpConfig['priorMean']
            if not self.priorMean.InheritsFrom(ROOT.TF1.Class()):
                raise "Prior mean must be a TF1 for now..."
        # Set the training range. If none was given will default to the entire histogram.
        self.trainMin = None  # The range to consider when training the GP
        self.trainMax = None
        if 'train_range' in gpConfig:
            self.trainMin = gpConfig['train_range'][0]  # The range to consider when training the GP
            self.trainMax = gpConfig['train_range'][1]

        # Setup the Kernel
        self.lengthScale = None
        self.amplitude = None
        self.k1 = RBF(60, (1, 200)) #squared exponential kernel
        self.k2 = C((1e3)**2, (1e-10, 1e15))

        if 'lengthScale' in gpConfig:
            self.lengthScale = gpConfig['lengthScale']
            print "Found lengthScale: {} in config".format(str(self.lengthScale))
            if self.lengthScale: # make sure its not set to None in the config
                print "******** LengthScale: {}".format(self.lengthScale)
                self.k1 = RBF(self.lengthScale, (self.lengthScale, self.lengthScale)) #squared exponential kernel

        if 'amplitude' in gpConfig:
            self.amplitude = gpConfig['amplitude']
            if self.amplitude: #make sure its not set to None in the config
                print "******** Amplitude: {}".format(self.amplitude)
                self.k2 = C(self.amplitude, (self.amplitude, self.amplitude)) #Constant kernel

        # asseble the Kernel
        self.kernel = self.k1 * self.k2
        print "Setting GP kernel: {}".format(self.kernel)

        # Internal vars to keep track of things
        self.currentNSig = None
        self.currentSBHist = None

        self.setTrainData(gpConfig['trainHisto'], gpConfig['trainHisto'])
        ####  End of __init__()


    def setTrainData(self, trainHisto, dataHisto, binErrors=False):
        """Train a GP on a histogram to get hyperparamters.

        Use a high stats MC sample to optimize the hyperparamters then
        return a kernel object to be used in the data fit.
        """
        print "===== Optimizing hyperparamters on the training sample."
        self.dataHisto = dataHisto
        self.dataBinErrors = binErrors
        self.dataXmin = self.dataHisto.GetBinLowEdge(1)
        self.dataXmax = self.dataHisto.GetBinLowEdge(self.dataHisto.GetNbinsX())+self.dataHisto.GetBinWidth(1)
        
        # Get the signal function
        self.sigFunction = ROOT.TF1("dscb", DSCB, self.dataXmin, self.dataXmax, 1)
        self.sigFunction.SetParameters(self.sigMass,0)
        self.sigNorm = self.sigFunction.Integral(self.dataXmin, self.dataXmax )



        if not self.dataHisto:
            raise "Must set fit data before training data. Use setFitData( dataHisto, binErrors)"

        self.trainHisto = trainHisto.Clone()
        print "Scaling training histo to data:", self.dataHisto.Integral()/self.trainHisto.Integral()
        self.trainHisto.Scale(self.dataHisto.Integral()/self.trainHisto.Integral())
        self.trainHisto.Sumw2()

        X_t, Y_t, dy_t = histoToArray(self.trainHisto, self.trainMin, self.trainMax)
        binErr = Y_t

        #Get the prior mean shape. Subtract if from the training histo.
        if self.priorMean:
            self.trainHisto.Add(self.priorMean, -1.0)
            X_t, Y_t, dy_t = histoToArray(self.trainHisto, self.trainMin, self.trainMax)

        #Need Y values from original histogram to get sqrt(N) errors correctly
        print "Running GP with sqrt(N) errors."
        gp = GaussianProcessRegressor(kernel=self.kernel
                                    ,alpha=binErr    # is really sqrt(Y_t_nom)^2  Input should be error bars squared
                                    ,n_restarts_optimizer=10)

        # Fit for the hyperparameters.
        gp.fit(np.atleast_2d(X_t).T, Y_t)
        y_pred, sigma = gp.predict(np.atleast_2d(X_t).T, return_std=True)

        self.trainResult = arrayToHisto("GPhisto", self.trainMin, self.trainMax, y_pred)
        self.opt_kernel = gp.kernel_
        print "Optimized hyperparameters:"
        print gp.kernel_

        X, y, dataErrs = histoToArray(dataHisto, self.trainMin, self.trainMax)
        #GPh = GPHisto(dataHisto)
        #dataErrs = GPh.getErrArr()
        if self.dataBinErrors:
            self.gp = GaussianProcessRegressor(kernel=self.opt_kernel  #previously optimized kernel
                                                ,optimizer=None         # Dont reoptimize hyperparamters.
                                                ,alpha=dataErrs**2
                                                )
        else:
            self.gp = GaussianProcessRegressor(kernel=self.opt_kernel  #previously optimized kernel
                                                ,optimizer=None         # Dont reoptimize hyperparamters.
                                                ,alpha=y                # sqrt(y)^2. Seems that the error should be squared here.
                                                )

    def getSigHisto(self, nSig):
        """Get histogram used for signal subtration. For testing.

        """
        sig_func = self.sigFunction
        nSigScale = nSig*self.dataHisto.GetBinWidth(1) / self.sigNorm
        sigHisto = self.dataHisto.Clone("sigHisto")
        sigHisto.Reset()
        sigHisto.Add(sig_func, nSigScale)
        return sigHisto



    def getGPHisto(self, nSig):
        """Get the histogram result of fitting the S+B GP. For testing.

        """
        #if self.currentNSig == nSig:
        #    return self.currentSBHist
        #else:
        sig_func = self.sigFunction
        nSigScale = nSig*self.dataHisto.GetBinWidth(1) / self.sigNorm
        sigSubtractedHist = self.dataHisto.Clone("sigSubtractedHist")
        sigSubtractedHist.Add(sig_func, -1*nSigScale) #Subtrac the signal from the data
        #Get the Array corresponding to the subtracted histogram
        X, y_subtracted, errArr = histoToArray(sigSubtractedHist, self.trainMin, self.trainMax)
        #print "X array size", X.size
        # fit the GP to what is left. Then get a GP prediction on the result.
        self.gp.fit(np.atleast_2d(X).T,y_subtracted)
        y_pred, sigma = self.gp.predict(np.atleast_2d(X).T, return_std=True)

        gpHisto = arrayToHisto("GPhisto", self.trainMin, self.trainMax, y_pred)
        gpHistoSigSub = gpHisto.Clone('GPSigSub')
        gpHisto.Add(sig_func, nSigScale) #Add the signal back into the data
        return gpHisto, sigSubtractedHist,gpHistoSigSub

    def evaluate(self, myy, nSig):
        """Overide the base class evaluate function.

        Get the signal distribution for the nSig that the PDF is trying and
        convert it to a histogram. Then subtract that histogram from the data
        histogram. Then fit the GP to the remaining histo. Return the resulting
        signal histo + the GP fit of the background for the myy bin in question.
        """
        if not (nSig == self.currentNSig): #If we have moved on to a different nSig, redo the GP
            self.currentNSig = nSig
            sig_func = self.sigFunction
            nSigScale = nSig*self.dataHisto.GetBinWidth(1) / self.sigNorm

            #subtract prior mean and DSCB distribution from background
            sigSubtractedHist = self.dataHisto.Clone("sigSubtractedHist")
            if self.priorMean:
                sigSubtractedHist.Add(self.priorMean, -1.0)

            sigSubtractedHist.Add(sig_func, -1*nSigScale) #Subtrac the signal from the data

            #Get the Array corresponding to the subtracted histogram
            X, y_subtracted, errArr = histoToArray(sigSubtractedHist, self.trainMin, self.trainMax)

            # fit the GP to what is left. Then get a GP prediction on the result.
            self.gp.fit(np.atleast_2d(X).T,y_subtracted)
            y_pred, sigma = self.gp.predict(np.atleast_2d(X).T, return_std=True)

            #gpHisto = arrayToHisto("GPhisto", self.dataXmin, self.dataXmax, y_pred)
            gpHisto = arrayToHisto("GPhisto", self.trainMin, self.trainMax, y_pred)
            gpHisto.Add(sig_func, nSigScale) #Add the signal back into the data
            if self.priorMean:
                gpHisto.Add(self.priorMean,1.0)
            if self.currentSBHist:
                self.currentSBHist.Delete()
            self.currentSBHist = gpHisto

            #if gpHisto.GetBinContent(gpHisto.FindBin(myy)) <= 0.0:
            #    return 0
            return gpHisto.GetBinContent(gpHisto.FindBin(myy))

        else:  #If we are still on the same nSig number just return the stored histogram.
            #print "=== evaluate return value:", self.currentSBHist.GetBinContent(self.currentSBHist.FindBin(myy))
            #if self.currentSBHist.GetBinContent(self.currentSBHist.FindBin(myy)) <= 0.0:
            #    return 0
            return self.currentSBHist.GetBinContent(self.currentSBHist.FindBin(myy))

    def analyticalIntegral(self, myyMin, myyMax):
        """Tell RooFit how to integrate the GP histo.

        Do a simple histogram integral.
        """
        minBin = self.currentSBHist.FindBin(myyMin)
        maxBin = self.currentSBHist.FindBin(myyMax)
        return self.currentSBHist.Integral(minBin, maxBin)
