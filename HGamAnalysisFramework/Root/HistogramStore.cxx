#include "HGamAnalysisFramework/HistogramStore.h"
#include "HGamAnalysisFramework/HgammaUtils.h"
#include "TH1.h"
#include "TH1F.h"
#include "TString.h"
#include <iostream>
#include <vector>

// Create and store TH1F histogram
void HistogramStore::createTH1F(TString name, int Nbins, double xmin, double xmax, TString title) 
{
  if (hasTH1F(name)) HG::fatal("HistogramStore::createHistoTH1F: Attempt to create second histogram named " + name);
  m_histoTH1F[name] = new TH1F(name,title,Nbins,xmin,xmax);
  m_histoTH1F[name]->Sumw2();
}

// Create and Store TH1F histogram
void HistogramStore::createTH1F(TString name, const std::vector<double> &bins, TString title) 
{
  if (hasTH1F(name)) HG::fatal("HistogramStore::createHistoTH1F: Attempt to create second histogram named " + name);
  m_histoTH1F[name] = new TH1F(name,title,-1+bins.size(),&bins[0]);
  m_histoTH1F[name]->Sumw2();
}

// Create and store TH2F histogram
void HistogramStore::createTH2F(TString name, int NbinsX, double xmin, double xmax, int NBinsY, double ymin, double ymax, TString title) 
{
  if (hasTH2F(name)) HG::fatal("HistogramStore::createHistoTH2F: Attempt to create second histogram named " + name);
  m_histoTH2F[name] = new TH2F(name, title, NbinsX, xmin, xmax, NBinsY, ymin, ymax);
  m_histoTH2F[name]->Sumw2();
}

// Create and store TH2F histogram
void HistogramStore::createTH2F(TString name, const std::vector<double> &xbins, const std::vector<double> &ybins, TString title) 
{
  if (hasTH2F(name)) HG::fatal("HistogramStore::createHistoTH2F: Attempt to create second histogram named " + name);
  m_histoTH2F[name] = new TH2F(name, title, xbins.size() - 1, &xbins[0],  ybins.size() - 1, &ybins[0]);
  m_histoTH2F[name]->Sumw2();
}

// Create and store TH3F histogram
void HistogramStore::createTH3F(TString name, int NbinsX, double xmin, double xmax, int NBinsY, double ymin, double ymax, int NBinsZ, double zmin, double zmax, TString title)
{
  if (hasTH3F(name)) HG::fatal("HistogramStore::createHistoTH3F: Attempt to create second histogram named " + name);
  m_histoTH3F[name] = new TH3F(name, title, NbinsX, xmin, xmax, NBinsY, ymin, ymax, NBinsZ, zmin, zmax);
  m_histoTH3F[name]->Sumw2();
}

// Create and store TH3F histogram
void HistogramStore::createTH3F(TString name, const std::vector<double> &xbins, const std::vector<double> &ybins, const std::vector<double> &zbins, TString title)
{
  if (hasTH3F(name)) HG::fatal("HistogramStore::createHistoTH3F: Attempt to create second histogram named " + name);
  m_histoTH3F[name] = new TH3F(name, title, xbins.size() - 1, &xbins[0], ybins.size() - 1, &ybins[0], zbins.size() - 1, &zbins[0]);
  m_histoTH3F[name]->Sumw2();
}

// Create and store TProfile histogram
void HistogramStore::createTProfile(TString name, int NbinsX, double xmin, double xmax, TString title) 
{
  if (hasTProfile(name)) HG::fatal("HistogramStore::createHistoTProfile: Attempt to create second histogram named " + name);
  m_histoTProfile[name] = new TProfile(name, title, NbinsX, xmin, xmax);
  m_histoTProfile[name]->Sumw2();
}

// Create and store TProfile histogram
void HistogramStore::createTProfile(TString name, const std::vector<double> &xbins, TString title) 
{
  if (hasTProfile(name)) HG::fatal("HistogramStore::createHistoTProfile: Attempt to create second histogram named " + name);
  m_histoTProfile[name] = new TProfile(name, title, xbins.size() - 1, &xbins[0]);
  m_histoTProfile[name]->Sumw2();
}




// Retrieve a TH1F histogram from the internal store
TH1F* HistogramStore::getTH1F(TString name)
{
  if (!hasTH1F(name)) HG::fatal("HistoStore::getTH1F requested histogram " + name + " cannot be accessed (did you forget to declare it?)");
  return m_histoTH1F[name];
}

// Retrieve a TH2F histogram from the internal store
TH2F* HistogramStore::getTH2F(TString name)
{
  if (!hasTH2F(name)) HG::fatal("HistoStore::getTH2F requested histogram " + name + " cannot be accessed (did you forget to declare it?)");
  return m_histoTH2F[name];
}

// Retrieve a TH3F histogram from the internal store
TH3F* HistogramStore::getTH3F(TString name)
{
  if (!hasTH3F(name)) HG::fatal("HistoStore::getTH3F requested histogram " + name + " cannot be accessed (did you forget to declare it?)");
  return m_histoTH3F[name];
}

// Retrieve a TH2F histogram from the internal store
TProfile* HistogramStore::getTProfile(TString name)
{
  if (!hasTProfile(name)) HG::fatal("HistoStore::getTProfile " + name + " cannot be accessed (did you forget to declare it?)");
  return m_histoTProfile[name];
}



// Return vector of all histograms in internal store
std::vector<TH1*> HistogramStore::getListOfHistograms()
{
  std::vector<TH1*> allHistos;
  // an iterator of a map is pair<keyType,valueType>
  for (auto h : m_histoTH1F ) allHistos.push_back(h.second);
  for (auto h : m_histoTH2F ) allHistos.push_back(h.second);
  for (auto h : m_histoTH3F ) allHistos.push_back(h.second);
  for (auto h : m_histoTProfile ) allHistos.push_back(h.second);
  return allHistos;
}
