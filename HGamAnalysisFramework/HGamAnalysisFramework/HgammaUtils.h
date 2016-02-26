#ifndef HGamAnalysisFramework_HgammaUtils_H
#define HGamAnalysisFramework_HgammaUtils_H

#include "HGamAnalysisFramework/HgammaIncludes.h"


//! \brief Hgamma namespace
namespace HG {

  
  //! \name   A few general helper methods and definitions
  //! \author Dag Gillberg
  //@{
  
  //! \brief typedef for a vector of doubles (to save some typing)
  typedef std::vector<double>  NumV;
  //! \brief typedef for a vector of ints (to save some typing)
  typedef std::vector<int>     IntV;
  //! \brief typedef for a vector of strings (to save some typing)
  typedef std::vector<TString> StrV;

  //! \brief Converts a text line to a vector of words
  //  \param str input string with words
  //  \param sep separator to define where a word ends or starts
  StrV vectorize(TString str, TString sep=" ");
  
  //! \brief Converts string of separated numbers to vector<double>
  //  \param str input string with numbers
  //  \param sep separator to define where a number ends or starts
  NumV vectorizeNum(TString str, TString sep=" ");
  
  //! \brief method to abort program with error message
  void fatal(TString msg);

  //! \brief returns true if a given file or directory exist
  bool fileExist(TString fn);

  //! \brief calculates DeltaR in (y,phi)-space instead of (eta,phi) given by p4().DeltaR()
  inline double DRrap(const TLorentzVector &p1, const TLorentzVector &p2) {
    double dy=p1.Rapidity()-p2.Rapidity(), dphi=p1.DeltaPhi(p2);
    return sqrt(dy*dy+dphi*dphi);
  }

  TH1* getHistogramFromFile(TString fname, TString hname);
  
#ifndef __MAKECINT__

  //! \brief print 4-vector as string
  TString fourVecAsText(const TLorentzVector &p4);
  TString fourVecAsText(const xAOD::IParticle *p);
  //const char* fourVecAsCharStar(const xAOD::IParticle *p);

  //! \brief calculates DeltaR in (y,phi)-space instead of (eta,phi) given by p4().DeltaR()
  inline double DRrap(const xAOD::IParticle *p1, const xAOD::IParticle *p2) {
    return DRrap(p1->p4(),p2->p4());
  }

  //! \brief calculates DeltaR in (eta,phi)-space
  inline double DR(const xAOD::IParticle *p1, const xAOD::IParticle *p2) {
    return p1->p4().DeltaR(p2->p4());
  }
  
  //! \brief returns smallest DR between ptcl and any of the objects in ptcls
  //!        if ptcl occurs in the list of particles, it is ignored
  template <class T> double minDR(const xAOD::IParticle *ptcl, T ptcls) {
    double mindr = 99;
    for ( auto p : ptcls ) if ( p!=ptcl && DR(ptcl,p)<mindr ) mindr=DR(ptcl,p);
    return mindr;
  }
  
  //! \brief returns smallest DR between ptcl and any of the objects in ptcls in (y,phi)-space
  //!        if ptcl occurs in the list of particles, it is ignored
  template <class T> double minDRrap(const xAOD::IParticle *ptcl, T ptcls) {
    double mindr = 99;
    for ( auto p : ptcls ) if ( p!=ptcl && DRrap(ptcl,p)<mindr ) mindr=DRrap(ptcl,p);
    return mindr;
  }
  
#endif
  
  //@}
}

#endif
