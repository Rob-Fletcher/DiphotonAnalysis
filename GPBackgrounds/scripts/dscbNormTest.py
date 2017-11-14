import ROOT
import sys
sys.path.append('../GPBackgrounds')
from SignalModel import *

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

f = ROOT.TFile("../data/templateHists_inclusive.root")
stats = 1000
bkgTemp = f.Get('hmgg_c0')
bkgTemp.Rebin(8)

datasighist = buildSignal(125,stats, bkgTemp.GetNbinsX())

integralHisto = datasighist.Integral()

dscb_func = ROOT.TF1("dscb", DSCB, 105, 160, 7)
dscb_func.SetParameters(401  # Normalization
                        ,125  # mu
                        ,1.475   # alpha_low
                        ,1.902   # alpha_high
                        ,12.1   # n_low
                        ,11.6   # n_high
                        ,1.86  ) # sigma

integralFunc = dscb_func.Integral(105, 160)

print "Histo Integral:", integralHisto, "Function Integral:", integralFunc

canv1 = ROOT.TCanvas('c1','c1')
datasighist.Draw()
dscb_func.Draw('same')
canv1.Print('NormTest.pdf')
