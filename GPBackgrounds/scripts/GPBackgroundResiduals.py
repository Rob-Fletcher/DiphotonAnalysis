from ROOT import *
import os
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF,Matern, ConstantKernel as C
#from root_numpy import *
#from matplotlib import pyplot as plt
import argparse
import sys
sys.path.append('../GPBackgrounds')
from RootToNp import *
from SignalModel import *
from MYYKernel import *


def run(args, mass, winLow, winHigh):
    f = TFile(args.input)
    bkghist_template = f.Get('hmgg_c0')
    bkghist_template.Rebin(8)

    stats = 100000
    seed = 10

    bkghist = toyModel(bkghist_template, stats, seed)

    if args.doSig:
        #get signal hist
        sighist = buildSignal(125,1000, bkghist.GetNbinsX())

        #inject signal into background
        bkghist.Add(sighist)

    GPh = GPHisto(bkghist)
    GPh.setWindow(winLow,winHigh)
    X = GPh.getXWindowArr()
    y = GPh.getYWindowArr()
    dy = GPh.getErrWindowArr()

    X_t = GPh.getXArr()
    y_t = GPh.getYArr()
    dy_t = GPh.getErrArr()

    if args.noWindow:
        X = X_t
        y = y_t
        dy = dy_t

    #X, y, dy = histoToArrayTest(bkghist, 120, 140)
    #X, y, dy = histoToArrayCut(bkghist, 120, 125)
    #X_t, y_t, dy_t = histoToArray(bkghist)

    #X = np.atleast_2d(X).T
    #y = y.ravel()
    #dy = dy.ravel()

    #x = np.atleast_2d(np.linspace(start=105, stop=160, num=1000)).T  # Predict a relatively smooth function
    #x = np.atleast_2d(np.linspace(start=105, stop=160, num=219)).T  # Predict at each data point

    x = GPh.getXArr()

    #kernel = C(800.0, (1e-3, 1e3)) * RBF(100.0, (1e-3, 1e3)) #squared exponential kernel
    #kernel = C(10.0, (1e-3, 1e15)) * RBF(np.sqrt(2)*(7**2), (1e-3,1e5 )) #squared exponential kernel
    kernel = C(1000.0, (1e-3, 1e15)) * FallExp(1.0, (1e-5, 1e2), 1.0, (1e-3,1e15)) * Gibbs(1.0, (1e-3, 1e5), 1.0, (1e-3,1e5))
    #kernel = C(10.0, (1e-3, 1e6)) * Gibbs(1.0, (1e-3, 1e5), 1.0, (1e-3,1e5))


    print "dy[5] =",dy[5]
    print "err =", bkghist.GetBinError(5), "Original =", bkghist_template.GetBinError(5)
    gp = GaussianProcessRegressor(kernel=kernel
                                    ,optimizer='fmin'
                                    ,alpha=dy**2
                                    ,n_restarts_optimizer=15
                                    )

    gp.fit(X,y)
    print gp.kernel_
    y_pred, sigma = gp.predict(x, return_std=True)


    if args.mplot:
        fig = plt.figure()
        #plt.plot(X, y, 'r.', markersize=10, label=u'Background')
        plt.errorbar(X.ravel(), y, dy, fmt='r.', markersize=8, label=u'Training Points', zorder=2)
        plt.errorbar(X_t.ravel(), y_t, dy_t, fmt='k.', markersize=7, label=u'Background', zorder=1)
        plt.plot(x, y_pred, 'b-', label=u'Prediction', zorder=3)
        plt.fill(np.concatenate([x, x[::-1]]),
                 np.concatenate([y_pred - 1.9600 * sigma,
                                (y_pred + 1.9600 * sigma)[::-1]]),
                 alpha=.5, fc='b', ec='None', label='95% confidence interval', zorder=3)
        plt.xlabel('$M_{\gamma \gamma}$')
        plt.ylabel('$events$')
        plt.title('Optimized Kernel: {}'.format(gp.kernel_))
        #plt.yscale('log')
        #plt.ylim(-10, 20)
        plt.legend(loc='upper right')
        plt.savefig(args.tag+'GPFit.pdf')
        #plt.show()
    else:
        outfile = TFile('out.root','RECREATE')
        #outhist = arrayToHisto('GP Fit', 105, 160, y_pred, sigma)
        outhist = GPh.getHisto(y_pred, 1.96*sigma, 'GP Fit')
        if args.noWindow:
            bkgWindow = GPh.getHisto(y, dy, 'Full Background')
        else:
            bkgWindow = GPh.getWinHisto(y, dy, 'Full Background')
        bkgSubtracted = bkghist.Clone('bkgSubtracted')
        bkgSubtracted.Add(outhist,-1)  #Subtract background prediction from background with injected signal.

        canv = TCanvas('canv', 'canv')
        pad1 = TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0)
        pad1.SetGridx()
        pad1.Draw()
        pad1.cd()
        outhist.SetStats(0)
        bkghist_template.SetStats(0)
        bkgWindow.SetStats(0)
        bkgWindow.SetMarkerColor(kBlue)
        bkgWindow.SetLineColor(kBlue)
        outhist.SetMarkerColor(kBlack)
        outhist.SetLineColor(kBlack)
        print outhist.GetBinError(10)

        #bkgNorm = bkgWindow.Integral(1, bkgWindow.FindBin(winLow))
        #tmpNorm = bkghist_template.Integral(1,bkghist_template.FindBin(winLow))
        bkgNorm = bkgWindow.Integral()
        tmpNorm = bkghist_template.Integral()
        bkghist_template.Scale(bkgNorm/tmpNorm)
        bkghist_template.SetTitle(str(gp.kernel_)+" nToys: "+str(stats))
        print "Bin 24:  {0} : {1} : {2}".format((outhist.GetBinContent(24)-outhist.GetBinError(24)), bkghist_template.GetBinContent(24), (outhist.GetBinContent(24)+outhist.GetBinError(24))    )

        ####### Poly2 fit
        #canv4 = TCanvas('c4','c4')
        expPol_func = TF1("expPol","[0]*exp((x-100)/100 * ([1] + [2]*(x-100)/100))",105,160)
        expPol_func.SetParameters(0,0,0)
        expPol_func.SetParLimits(1,-10.,10.)
        expPol_func.SetParLimits(2,-10.,10.)
        bkgWindow.Fit("expPol","","",105,160)
        expFitResult = bkgWindow.GetFunction("expPol")
        expPolHist = expFitResult.GetHistogram()
        print expPolHist.GetNbinsX()
        #expPolHist.Divide(outhist)
        #expPolHist.Draw()
        #canv4.Print(args.tag+'/expPol_GP_ratio.pdf')

        bkghist_template.Draw('')
        bkgWindow.Draw('same')
        outhist.Draw('histsame')

        #outhist.GetYaxis().SetLabelSize(0.)
        axis = TGaxis( -5, 20, -5, 220, 20,220,510,"")
        axis.SetLabelFont(43)
        axis.SetLabelSize(15)
        axis.Draw()

        canv.cd()
        pad2 = TPad("pad2", "pad2", 0, 0.02, 1, 0.3)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.28)
        pad2.SetGridx()
        pad2.Draw()
        pad2.cd()


        h3 = bkghist_template.Clone("h3")
        h3.SetLineColor(kBlack)
        h3.SetMinimum(0.95)
        h3.SetMaximum(1.05)
        h3.Sumw2()
        h3.SetStats(0)
        h3.Divide(outhist)
        h3.SetMarkerColor(kBlack)
        h3.SetMarkerStyle(20)
        h3.SetMarkerSize(0.5)
        h3.Draw("ep")

        h4 = bkghist_template.Clone("h4")
        h4.SetLineColor(kRed)
        h4.SetMinimum(0.95)
        h4.SetMaximum(1.05)
        h4.Sumw2()
        h4.SetStats(0)
        h4.Divide(expFitResult)
        h4.SetMarkerColor(kRed)
        h4.SetMarkerStyle(20)
        h4.SetMarkerSize(0.5)
        h4.Draw("epsame")

        line = TLine(105,1,160,1)
        line.Draw('same')


        # outhist settings
        outhist.SetLineColor(kBlack);
        outhist.SetFillColorAlpha(33, 0.5)
        outhist.SetLineWidth(2);

        # Y axis outhist plot settings
        outhist.GetYaxis().SetTitleSize(20);
        outhist.GetYaxis().SetTitleFont(43);
        outhist.GetYaxis().SetTitleOffset(1.55);

        # bkghist settings
        bkghist.SetLineColor(kBlack);
        bkghist.SetMarkerSize(0.7)
        bkghist.SetLineWidth(2);

        # Ratio plot (h3) settings
        h3.SetTitle(""); # Remove the ratio title

        # Y axis ratio plot settings
        h3.GetYaxis().SetTitle("data/fit ");
        h3.GetYaxis().SetNdivisions(505);
        h3.GetYaxis().SetTitleSize(20);
        h3.GetYaxis().SetTitleFont(43);
        h3.GetYaxis().SetTitleOffset(1.);
        h3.GetYaxis().SetLabelFont(43); # Absolute font size in pixel (precision 3)
        h3.GetYaxis().SetLabelSize(15);

        # X axis ratio plot settings
        h3.GetXaxis().SetTitleSize(20);
        h3.GetXaxis().SetTitleFont(43);
        h3.GetXaxis().SetTitleOffset(3.);
        h3.GetXaxis().SetLabelFont(43); # Absolute font size in pixel (precision 3)
        h3.GetXaxis().SetLabelSize(15)

        canv.SetBottomMargin(0)

        canv.Write()
        #canv.Print(winLow+'_'+winHigh+'_GPFit.pdf')
        #canv.Print(args.tag+'/GPFit_'+str(winLow)+'_'+str(winHigh)+'.pdf')
        canv.Print(args.tag+'/GPFit_'+str(seed)+'.pdf')

        if args.doSig:
            ###  Plot signal stuff
            canv2 = TCanvas('c2','c2')
            canv2.cd()
            sighist.SetMarkerColor(kBlack)
            sighist.SetMarkerStyle(20)
            bkgSubtracted.GetXaxis().SetRangeUser(105,158)
            bkgSubtracted.Draw('hist')
            sighist.Draw('samep')
            #canv2.Write()
            canv2.Print(args.tag+'/SigYield_root.pdf')

            canv3 = TCanvas('c3', 'c3')
            canv3.cd()
            ratio = sighist.Clone('ratio')
            ratio.Divide(bkgSubtracted)
            ratio.GetYaxis().SetRangeUser(-5,5)
            ratio.Draw()
            #canv3.Write()
            canv3.Print(args.tag+'/SigYield_Ratio_root.pdf')

        """
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
        """

        #canv.cd()
        #bkgSubtracted.GetXaxis().SetRangeUser(120,130)
        #bkgSubtracted.Draw()
        #fitResult.Draw('same')
        #canv.Print(args.tag+'/fitResult.pdf')
        #print fitResult.Integral(120,130)


        """
        canv4 = TCanvas('c4','c4')
        gp_pred_full = GPh.getHisto(y_pred_full, sigma_full, 'GP Fit full')
        gp_pred_full.Divide(outhist)
        gp_pred_full.GetYaxis().SetRangeUser(0.95,1.05)
        gp_pred_full.Draw()
        canv4.Print(args.tag+'/Full_window_ratio.pdf')
        """
        f.Close()
        #outfile.Close()
        #return ss, norm

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

    ssHisto = TH1F("ss","ss",20, -20,20)

    args.windowSize = 10 #GeV
    run(args, 125,120,130)
    """
    for mass in np.arange(115,150,0.25):
        winlow = mass - args.windowSize/2
        winHigh = mass + args.windowSize/2
        ss, norm = run(args, mass, winlow, winHigh)
        print "====== mass:", mass,"=========="
        print "SS is:", ss
        print "Norm is:", norm
        ssHisto.Fill(ss)

    canv = TCanvas('c','c')
    ssHisto.Draw()
    canv.Print(args.tag+'/spuriousSignal.pdf')
    """
