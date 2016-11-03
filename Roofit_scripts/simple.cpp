//Get the background histograms from the template file

#ifndef __CINT__
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooCategory.h"
#include "RooGenericPdf.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
using namespace RooFit ;
#endif
#include "AtlasStyle.h"

void simple(){
    TFile* background_file = new TFile("../histos/BkgEstimation_NoCategory_HKHI_Lin_all_20.7_L4.root");
    TFile* signal_600M1200 = new TFile("../histos/mc15a.Sherpa_ADDyy_MS3500_600M1200.MxAOD.p2610.h013x.root");
    TFile* signal_1200M1800 = new TFile("../histos/mc15a.Sherpa_ADDyy_MS3500_1200M1800.MxAOD.p2610.h013x.root");
    TFile* signal_1800M = new TFile("../histos/mc15a.Sherpa_ADDyy_MS3500_1800M.MxAOD.p2610.h013x.root");

    TH1 * h_total_background = (TH1*) background_file->Get("bkg_total_gg_full");
    TTree * tree_600M1200 = (TTree*) signal_600M1200->Get("CollectionTree");
    TTree * tree_1200M1800 = (TTree*) signal_1200M1800->Get("CollectionTree");
    TTree * tree_1800M = (TTree*) signal_1800M->Get("CollectionTree");
    //TH1 * h_reducible_background = (TH1*) gDirectory->Get("bkg_reducible_gg_full");
    //TH1 * h_gammajet_background = (TH1*) gDirectory->Get("bkg_gammajet_gg_full");
    //TH1 * h_jetgamma_background = (TH1*) gDirectory->Get("bkg_jetgamma_gg_full");
    //TH1 * h_jetjet_background = (TH1*) gDirectory->Get("bkg_jetjet_gg_full");

    RooRealVar mgg("HGamEventInfoAuxDyn.m_yy","HGamEventInfoAuxDyn.m_yy",200,3500);
    //RooRealVar sig_mgg("HGamEventInfoAuxDyn.m_yy", "HGamEventInfoAuxDyn.m_yy",200, 3500);
    RooCategory isPassed("HGamEventInfoAuxDyn.isPassed", "HGamEventInfoAuxDyn.isPassed");

    RooDataHist total_background("totalBackground", "total background", mgg, Import(*h_total_background));

    RooDataSet signal_total("signal_600M1200","signal_600M1200", RooArgList(mgg,isPassed), Import(*tree_600M1200), Cut("isPassed == 1"));
    RooDataSet signal1("signal_1200M1800","signal_1200M1800", RooArgList(mgg,isPassed), Import(*tree_1200M1800), Cut("isPassed ==1"));
    RooDataSet signal2("signal_1800M","signal_1800M", RooArgList(mgg,isPassed), Import(*tree_1800M), Cut("isPassed == 1"));

    signal_total.append(signal1);
    signal_total.append(signal2);
    //RooDataHist reducible_background("reducible_background", "reducible background", x, h_reducible_background);
    //RooDataHist gammajet_background("gammajet_background", "gammajet background", x, h_gammajet_background);
    //RooDataHist jetgamma_background("jetgamma_background", "jetgamma background", x, h_jetgamma_background);
    //RooDataHist jetjet_background("jetjet_background", "jetjet background", x, h_jetjet_background);

//    RooRealVar mass_var("mass_gev", "mass_gev", 200, 3500);
    RooRealVar par0("par0", "par0", 0.3, 0, 10);
    RooRealVar par1("par1", "par1", 65, 0, 100);
    RooRealVar par2("par2", "par2", 3, 0, 20);
    RooRealVar par3("par3", "par3", 0.1, 0, 5);
    RooArgList bkg_list(mgg, par0, par1, par2, par3);
    RooGenericPdf bkg_model( "bkg_model", "bkg_model",
                                    "par0*(1-mgg/13000)^(par1)*(mgg/13000)^(-par2-par3*log(mgg/13000))", bkg_list);

    bkg_model.fitTo(total_background); //chi2FitTo(data, ROOT.RooFit.NumCPU(20), ROOT.RooCmdArg())

    TCanvas * c = new TCanvas("c");
    RooPlot* xframe = mgg.frame();
    total_background.plotOn(xframe,DrawOption("hist"));
    bkg_model.plotOn(xframe);
    signal_total.plotOn(xframe);
    c->SetLogy();
    xframe->Draw();
};
