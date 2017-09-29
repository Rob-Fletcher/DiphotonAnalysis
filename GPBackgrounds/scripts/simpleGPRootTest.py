from ROOT import *
import numpy as np
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF,Matern, ConstantKernel as C
import sys
sys.path.append('../GPBackgrounds')
from RootToNp import *
from SignalModel import *
from MYYKernel import *
import sys

sigmaMin = 1e-1
lMin = 1e-1

gRandom.SetSeed(int(sys.argv[1]))
func = TF1('xsinx', 'x*sin(x)+10', 0, 10)

hist = TH1F('noisy', 'noisy', 20, 0, 10)
hist.FillRandom('xsinx', int(sys.argv[2]))


GPh = GPHisto(hist)

X = GPh.getXArr()
y = GPh.getYArr()
dy = GPh.getErrArr()

x = np.atleast_2d(np.linspace(start=0., stop=10, num=1000)).T  # Predict at each data point

#kernel = C(1.0, (sigmaMin, 1e5)) * RBF(2.0, (lMin, 1e3)) #squared exponential kernel
kernel = C(1.0, 1e-3,1e5) * FallExp() * Gibbs()

gp = GaussianProcessRegressor(kernel=kernel
                                #,optimizer=None
                                ,alpha=(dy**2)
                                ,n_restarts_optimizer=15
                                )

gp.fit(X,y)
print gp.kernel_
y_pred, sigma = gp.predict(x, return_std=True)

outhist = TH1F('GP','GP', 1000,0,10)
for index,cont in enumerate(y_pred):
    outhist.SetBinContent(index+1, cont)
    outhist.SetBinError(index+1, 1.96*sigma[index])


canv = TCanvas('plot')
canv.cd()

funcHist = func.GetHistogram()
funcHist.Scale(hist.Integral()*hist.GetBinWidth(2)/(funcHist.Integral()*funcHist.GetBinWidth(2)))

hist.SetLineColor(kBlue)
hist.SetMarkerStyle(20)
hist.SetMarkerColor(kBlue)

outhist.SetLineColor(kGray)
outhist.SetMarkerColor(kBlack)
#outhist.SetMarkerStyle(20)
#outhist.SetFillColor(kGrey)
outhist.GetYaxis().SetRangeUser(-20, int(sys.argv[2])/10)

meanhist = outhist.Clone('mean')
meanhist.SetMarkerSize(0.5)

outhist.Draw('pe4')
meanhist.Draw('same')
funcHist.Draw('histsame')
hist.Draw('samepe1')

canv.Print('SimpleGPTest.pdf')
