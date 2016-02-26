#ifndef HISTOGRAM_STORE_H
#define HISTOGRAM_STORE_H

/********************
 * HistogramStore: 
 *   class to create, store, fill and retrieve TH1 and TH2 histograms
 *
 * Usage:
 *   HistogramStore HistoStore;
 *   HistoStore.createTH1F("Nphotons",40,-0.5,39.5,";#it{N}_{photon-clusters}");
 *   vector<TH1*> AllHistos = HistoStore.getListOfHistograms();
 *
 */

#include <map>
#include <vector>

#include <TString.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TH3.h>
#include <TH3F.h>
#include <TProfile.h>

/*! \brief Class that handles creating, filling, and retrieving histograms
 *        through use of an internal map
 *  \author Nathan Readioff
 */
class HistogramStore
{
 private:
  std::map<TString,TH1F*>     m_histoTH1F;
  std::map<TString,TH2F*>     m_histoTH2F;
  std::map<TString,TH3F*>     m_histoTH3F;
  std::map<TString,TProfile*> m_histoTProfile;

 public:
  //! \brief Create and store TH1F histogram
  void createTH1F(TString name, int Nbins, double xmin, double xmax, TString title = "");

  //! \brief Create and store TH1F histogram
  void createTH1F(TString name, const std::vector<double> &bins, TString title = "");

  //! \brief Create and store TH2F histogram
  void createTH2F(TString name, int NbinsX, double xmin, double xmax, int NBinsY, double ymin, double ymax, TString title = ""); 

  //! \brief Create and store TH2F histogram
  void createTH2F(TString name, const std::vector<double> &xbins, const std::vector<double> &ybins, TString title = ""); 

  //! \brief Create and store TH3F histogram
  void createTH3F(TString name, int NbinsX, double xmin, double xmax, int NBinsY, double ymin, double ymax, int NBinsZ, double zmin, double zmax, TString title = "");

  //! \brief Create and store TH3F histogram
  void createTH3F(TString name, const std::vector<double> &xbins, const std::vector<double> &ybins, const std::vector<double> &zbins, TString title = "");

  //! \brief Create and store TProfile histogram
  void createTProfile(TString name, int NbinsX, double xmin, double xmax, TString title = ""); 

  //! \brief Create and store TProfile histogram
  void createTProfile(TString name, const std::vector<double> &xbins, TString title = ""); 

  
  
  //! \brief Fill existing TH1F histogram
  inline void fillTH1F(TString name, double x, double w=1.0) {getTH1F(name)->Fill(x,w);}
  
  //! \brief Fill existing TH2F histogram
  inline void fillTH2F(TString name, double x, double y, double w=1.0) {getTH2F(name)->Fill(x, y, w);}
  
  //! \brief Fill existing TH3F histogram
  inline void fillTH3F(TString name, double x, double y, double z, double w=1.0) {getTH3F(name)->Fill(x, y, z, w);}
  
  //! \brief Fill existing TProfile histogram
  inline void fillTProfile(TString name, double x, double y, double w=1.0) {getTProfile(name)->Fill(x, y, w);}
  
  
  //! \brief check whether a given TH1F exist in the store
  inline bool hasTH1F(TString name) {return m_histoTH1F.count(name)>0;}

  //! \brief check whether a given TH2F exist in the store
  inline bool hasTH2F(TString name) {return m_histoTH2F.count(name)>0;}

  //! \brief check whether a given TH3F exists in the store
  inline bool hasTH3F(TString name) {return m_histoTH3F.count(name)>0;}

  //! \brief check whether a given TH2F exist in the store
  inline bool hasTProfile(TString name) {return m_histoTProfile.count(name)>0;}
  


  //! \brief Retrieve TH1F histogram from internal store
  TH1F* getTH1F(TString name);

  //! \brief Retrieve TH2F histogram from internal store
  TH2F* getTH2F(TString name);

  //! \brief Retrieve TH3F histogram from internal store
  TH3F* getTH3F(TString name);

  //! \briefRetrieve TProfile histogram from internal store
  TProfile* getTProfile(TString name);

  //! \brief Retrieve List of all histograms in internal store
  std::vector<TH1*> getListOfHistograms();
};

#endif
