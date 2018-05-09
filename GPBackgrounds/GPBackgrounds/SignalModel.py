from math import exp
from ROOT import *

def DSCB(myy, par):
    norm = par[0]
    mu = par[1]
    alpha_low = par[2]
    alpha_high = par[3]
    n_low = par[4]
    n_high = par[5]
    sigma = par[6]

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


def buildSignal(mass, nevents, nbins):
    """Retuns a histogram with signal that can be injected into background.

    """
    histo = TH1F('sig', 'sig', nbins, 60, 120)
    dscb_func = TF1("dscb", DSCB, 60, 120, 7)

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
    aLo = (aLoSignal0 + aLoSignal1*mnX)
    nLo = (nLoSignal0)
    aHi = (aHiSignal0 + aHiSignal1*mnX)
    nHi = (nHiSignal0)

    mu  = (mass + deltaMSignal)

    dscb_func.SetParameters(1  # Normalization
                            ,mu  # mu
                            ,aLo   # alpha_low
                            ,aHi   # alpha_high
                            ,nLo   # n_low
                            ,nHi   # n_high
                            ,sigma  ) # sigma

    histo.FillRandom("dscb",nevents)
    return histo

def toyModel(PDF, nEvents, seed):
    """Make a toy from a PDF.

    """
    gRandom.SetSeed(seed)
    histo = PDF.Clone('toyModel')
    histo.Reset()
    histo.FillRandom(PDF, nEvents)

    return histo
