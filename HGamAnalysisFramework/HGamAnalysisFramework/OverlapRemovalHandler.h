#ifndef HGamAnalysisFramework_OverlapRemoval
#define HGamAnalysisFramework_OverlapRemoval

#include "HGamAnalysisFramework/HgammaIncludes.h"

// \brief Hgamma namespace
namespace HG {
  
  /*! \brief Class that removes overlapping objects to avoid double counting
   *
   *  \details
   *  Class that removes overlapping objects in order to avoid
   *  double counting based on the HGam strategy: photons first!
   *  Currently the official tool is not used as it does not
   *  support this strategy, but instead follows the
   *  <a href="https://cds.cern.ch/record/1700874/files/ATL-COM-PHYS-2014-451.pdf">overlap
   *  harmonization recommendation</a>.
   *
   *  Run 1 HGam strategy<ul>
   *  <li>The two leading photons are always kept;
   *  <li>Electrons with &Delta;R(e,&gamma;) < 0.4 are removed;
   *  <li>Jets such as &Delta;R(jet,e) < 0.2 or &Delta;R(jet,&gamma;) < 0.4 are removed
   *  <li>Muons with &Delta;R(&mu;,jet) < 0.4 or &Delta;R(&mu;,&gamma;) < 0.4 are removed
   *  </ul>
   *
   *  \author Kaicheng Li
   *  \author Dag Gillberg
   */
  class OverlapRemovalHandler {
    
  public:
    
    //! \brief constructor
    OverlapRemovalHandler(TString name="OverlapRemoval");

    //! \brief destructor
    ~OverlapRemovalHandler();
    
    //! \brief initalizaiton
    EL::StatusCode initialize(Config &config);
    
    //! \brief removes overlap between objects according to configuration
    void removeOverlap(xAOD::PhotonContainer &photons,
                       xAOD::JetContainer &jets,
                       xAOD::ElectronContainer &electrons,
                       xAOD::MuonContainer &muons);
    
    //! \brief removes overlap between objects. A nullptr can be passed for
    //!        electrons, jets and/or muons which means no overlap will be
    //!        performed for these objects
    void removeOverlap(xAOD::PhotonContainer *photons,
                       xAOD::JetContainer *jets,
                       xAOD::ElectronContainer *electrons,
                       xAOD::MuonContainer *muons);
    
  private:

    TString m_name;

    enum MatchingMode { eta_phi=0, y_phi=1 };
    MatchingMode m_matchMode;

    double m_e_DR_y, m_jet_DR_y, m_jet_DR_e, m_e_DR_jet;
    double m_mu_DR_y, m_mu_DR_jet;

    //! \brief returns true if the ptcl matches any of the ptcls within DeltaR < DRcut
    template <class T> bool overlap(const xAOD::IParticle *ptcl, T ptcls, double DRcut) {

      // negative DRcut values means that overlaps are not considered
      if (DRcut<0) return false;      
      if      ( m_matchMode == eta_phi ) return HG::minDR(ptcl,ptcls)    < DRcut;
      else if ( m_matchMode == y_phi   ) return HG::minDRrap(ptcl,ptcls) < DRcut;
      else HG::fatal("Unsupported DR matching mode");
      
      return false;
    }

    //! \brief returns true if the ptcl doesn't match any of the ptcls within DeltaR < DRcut
    template <class T> bool noOverlap(const xAOD::IParticle *ptcl, T ptcls, double DRcut) {
      return !overlap(ptcl,ptcls,DRcut);
    }

    //! \brief returns true if the ptcl matches any of the ptcls within DeltaEta < Deta_cut and DeltaPhi < Dphi_cut
    template <class T> bool overlap(const xAOD::IParticle *ptcl, T ptcls, double Deta_cut, double Dphi_cut) {
      for ( auto p : ptcls ) {
	if (p==ptcl) continue;
	if ( fabs(p->eta()-ptcl->eta())<Deta_cut && fabs(p->p4().DeltaPhi(ptcl->p4()))<Dphi_cut) return true;
      }
      return false;
    }
    
    

  public:
    
    //! \brief Removes any particle from ptcls that overlaps with any of particle in refPtcls
    template <class T1, class T2> void removeOverlap(T1 &ptcls, T2 refPtcls, double DRcut) {
      
      // Reverse for loop over probe-particles
      for ( auto p=ptcls.rbegin(); p!=ptcls.rend(); ++p) {

        // If the current probe-particle overlaps with any reference particle remove it!
        // Need to convert reverse_iterator to a normal iterator, see:
        // http://en.cppreference.com/w/cpp/iterator/reverse_iterator
        if ( overlap(*p,refPtcls,DRcut))
          ptcls.erase(p.base()-1);
      }
    }

    //! \brief Removes any particle from ptcls that overlaps with any of particles with higher pT in the same container
    template <class T1> void removeOverlap(T1 &ptcls, double Deta_cut,double Dphi_cut) {
      
      // Reverse for loop over probe-particles
      for ( auto p=ptcls.rbegin(); p!=ptcls.rend(); ++p) {

        // If the current probe-particle overlaps with any reference particle remove it!
        // Need to convert reverse_iterator to a normal iterator, see:
        // http://en.cppreference.com/w/cpp/iterator/reverse_iterator
        if ( overlap(*p,ptcls,Deta_cut,Dphi_cut))
          ptcls.erase(p.base()-1);
      }
    }

    //! \brief returns the subset of muons in jets.
    //  By default the same DRcut is used as for the muon-jet removal
    xAOD::MuonContainer muonsInJets(xAOD::MuonContainer muons,
                                    xAOD::JetContainer jets,
                                    double DRcut = -1);
    
    //! \brief returns the subset of ptcls that overlap with any particle in refPtcls
    template <class T1, class T2> T1 getOverlaps(T1 ptcls, T2 refPtcls, double DRcut) {
      T1 overlaps(SG::VIEW_ELEMENTS);
      for (auto p:ptcls)
        if (overlap(p,refPtcls,DRcut)) overlaps.push_back(p);
      return overlaps;
    }

    
  };
}

#endif // HGamAnalysisFramework_OverlapRemoval
