from ROOT import TH1F
import numpy as np

class GPHisto():
    """Class to hold np arrays with all needed information and also the histogram
    for easy converting between the two.

    """
    def __init__(self, histo):

        self._histo = histo
        self._binLowEdges = np.array([])
        self._binContent = np.array([])
        self._binErrors = np.array([])

        self._windowLowEdge = None
        self._windowUpEdge = None

        self._binLowEdgesWindow = np.array([])
        self._binContentWindow = np.array([])
        self._binErrorsWindow = np.array([])

        for nbin in range(1, self._histo.GetNbinsX()+1):
            self._binLowEdges = np.append(self._binLowEdges, [self._histo.GetBinCenter(nbin)])
            self._binContent = np.append(self._binContent, [self._histo.GetBinContent(nbin)])
            self._binErrors = np.append(self._binErrors, [self._histo.GetBinError(nbin)])

        pass # init

    def getXArr(self):
        """Get the np array defining the x axis

        """
        return np.atleast_2d(self._binLowEdges).T

    def getYArr(self):
        """Get the np array defining the bin content

        """
        return self._binContent

    def getErrArr(self):
        """Get the np array containing the bin errors

        """
        return self._binErrors

    def setWindow(self,winLowEdge, winUpEdge):
        """Set the window to remove from the background spectrum. Must be called
        before you try to get the arrays.

        Sets the window to remove from the background and fills arrays that
        contain the info with the window removed. Once this function is called
        you can return each array individually with their getters.

        """
        self._windowUpEdge = winUpEdge
        self._windowLowEdge = winLowEdge

        for nbin in range(1, self._histo.GetNbinsX()+1):
            if not (self._histo.GetBinLowEdge(nbin) < self._windowLowEdge or self._histo.GetBinLowEdge(nbin) > self._windowUpEdge):
                continue
            self._binLowEdgesWindow = np.append(self._binLowEdgesWindow, [self._histo.GetBinLowEdge(nbin)])
            self._binContentWindow = np.append(self._binContentWindow, [self._histo.GetBinContent(nbin)])
            self._binErrorsWindow = np.append(self._binErrorsWindow, [self._histo.GetBinError(nbin)])


    def getXWindowArr(self):
        """Get the X axis array with the window removed. Returns the 'at least 2d'
        array with the transpose taken so that the return value can be plugged
        directly into the GP.

        """
        if not self._binLowEdgesWindow.any():
            raise("You must call setWindow() before you can get the arrays with windows removed")
        return np.atleast_2d(self._binLowEdgesWindow).T

    def getYWindowArr(self):
        """Get the y-axis array with the window removed

        """
        if not self._binContentWindow.any():
            raise("You must call setWindow() before you can get the arrays with windows removed")
        return self._binContentWindow.ravel()

    def getErrWindowArr(self):
        """Get the Error array with the window removed

        """
        if not self._binErrorsWindow.any():
            raise("You must call setWindow() before you can get the arrays with windows removed")
        return self._binErrorsWindow.ravel()

    def getWinHisto(self, yArr, errArr, name='histo'):
        """Get the histogram corresponding to the base histo with the window
        removed.

        """
        if not self._binLowEdgesWindow.any():
            raise("You must call setWindow() before you can get the histo with windows removed")
        histo = self._histo.Clone(name)
        histo.Reset()
        for index,value in enumerate(self._binLowEdgesWindow):
            histo.SetBinContent(histo.FindBin(value), yArr[index])
            histo.SetBinError(histo.FindBin(value), errArr[index])

        return histo

    def getHisto(self, yArr, errArr, name='histo'):
        """Get the histogram with no window removed. This can be used to get the
        histogram for the prediction

        """
        histo = self._histo.Clone(name)
        histo.Reset()
        for index,value in enumerate(self._binLowEdges):
            histo.SetBinContent(histo.FindBin(value), yArr[index])
            histo.SetBinError(histo.FindBin(value), errArr[index])
        return histo

    pass # NpHisto



############### Older Methods

def histoToArray(histo):
    """Take a histogram and convert the data to a numpy array for use in
    GaussianProcessRegressor

    \return arr,errarr  retuns a tuple containing the array of data points and the array of errors.
    """
    xaxis = np.array([])
    yaxis = np.array([])
    err = np.array([])
    for nbin in range(1,histo.GetNbinsX()):
        xaxis = np.append(xaxis, [histo.GetBinLowEdge(nbin)])
        yaxis = np.append(yaxis, [histo.GetBinContent(nbin)])
        err   = np.append(err, [histo.GetBinError(nbin)])
    return xaxis,yaxis,err
    pass

def histoToArrayCut(histo, cutlow, cuthigh):
    """Take a histogram and convert the data to a numpy array for use in
    GaussianProcessRegressor


    \return xaxis,yaxis,errarr  retuns a tuple containing the array of x-axis values, data points and the array of errors.
    """
    xaxis = np.array([])
    yaxis = np.array([])
    err = np.array([])
    for nbin in range(1,histo.GetNbinsX()):
        if histo.GetBinCenter(nbin) > cutlow and histo.GetBinCenter(nbin) < cuthigh:
            continue
        xaxis = np.append(xaxis, [histo.GetBinLowEdge(nbin)])
        yaxis = np.append(yaxis, [histo.GetBinContent(nbin)])
        err   = np.append(err, [histo.GetBinError(nbin)])
    return xaxis,yaxis,err
    pass

def histoToArrayTest(histo, cutlow, cuthigh):
    """Take a histogram and convert the data to a numpy array for use in
    GaussianProcessRegressor


    \return xaxis,yaxis,errarr  retuns a tuple containing the array of x-axis values, data points and the array of errors.
    """
    xaxis = np.array([])
    yaxis = np.array([])
    err = np.array([])
    for nbin in range(1,histo.GetNbinsX()):
        if not nbin%22 == 0:
            continue
        if histo.GetBinCenter(nbin) > cutlow and histo.GetBinCenter(nbin) < cuthigh:
            continue
        xaxis = np.append(xaxis, [histo.GetBinLowEdge(nbin)])
        yaxis = np.append(yaxis, [histo.GetBinContent(nbin)])
        err   = np.append(err, [histo.GetBinError(nbin)])
    return xaxis,yaxis,err
    pass

def arrayToHisto(title, xmin, xmax, arr, errarr=None):
    """Convert a numpy array into a histogram to plot with root.

    """
    hist = TH1F(title, title, len(arr)+1,xmin,xmax)
    for nbin, binContent in enumerate(arr,start=1):
        hist.SetBinContent(nbin, binContent)
        hist.SetBinError(nbin, errarr[nbin-1]) # -1 because root starts at 1 but array index starts at 0

    return hist
    pass
