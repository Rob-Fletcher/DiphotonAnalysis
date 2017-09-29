from ROOT import TH1F
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF,Matern, ConstantKernel as C
from MYYKernel import *


class GPFitter():
    """Class to hold all histograms and ndarrays used to do GP fitting and
    return the actual fits. Converts histograms to ndarrays internally and
    returns histograms.

    Each instance of this class should define one GP. This means it has one kernel,
    and one training/fitting histogram. If a window fit needs to be done it will
    also be defined in here. If a separate non-window fit needs to be done it should
    be done with a second instance of the class.

    """
    def __init__(self, fitHisto=None, trainHisto=None, kernel='RBF', hParams={}):
        """Initialize all internal class variables.

        @param fitHisto The histogram containing the data points you would like to fit to.
        @param trainHisto Histogram that is only used to determine hyperparameters. If nothing is passed the fitHisto will be used to determine hyperparameters.
        @param kernel Kernel to use as the covariance function. Can only use ones implemented in here:
                        Options:
                          - 'RBF': Radial basis function. Has 2 parameters 'A' and 'length_scale'.
                          - 'Gibbs': Exponential decay * Gibbs kernel. Has 5 parameters 'A', 'a', 'b', 'c', 'd'
                          - 'FallingRBF': Exponential decay * RBF. Has 4 Parameters 'A', 'a', 'b', 'length_scale'

        """
        if fitHisto is None:
            raise ValueError("Must pass a fit histogram to GPFitter()!")
        self.fitHisto = fitHisto
        self.trainHisto = trainHisto
        self.kernelFunc = kernel #internally the self.kernel variable will hold the actual kernel object.
        self.hParams = hParams

    def predict(self, predPoints=None):
        """Get the GP prediction.

        Will run the full gaussian process on the fitHisto given when initializing.

        """

    def getParams(self):
