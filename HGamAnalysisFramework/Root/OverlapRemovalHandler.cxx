#include "HGamAnalysisFramework/OverlapRemovalHandler.h"

namespace HG {
  
  //______________________________________________________________________________
  OverlapRemovalHandler::OverlapRemovalHandler(TString name)
  : m_name(name) { }
  
  //______________________________________________________________________________
  OverlapRemovalHandler::~OverlapRemovalHandler() { }
  
  //______________________________________________________________________________
  EL::StatusCode OverlapRemovalHandler::initialize(Config &conf)
  {

    // Read the matching method (y,phi) or (eta,phi) DR matching
    TString matchStr = conf.getStr(m_name+".MatchingMode","RapidityPhi");
    if      (matchStr=="RapidityPhi") m_matchMode=y_phi;
    else if (matchStr=="EtaPhi") m_matchMode=eta_phi;
    else HG::fatal("Can not interpret MatchingMode "+matchStr);
    
    // Read in the DeltaR distances used in the overlap removal
    // Default values are from the Run 1 analysis.
    m_e_DR_y    = conf.getNum(m_name+".Electron_DR_Photon",0.4);
    m_jet_DR_y  = conf.getNum(m_name+".Jet_DR_Photon",0.4);
    m_jet_DR_e  = conf.getNum(m_name+".Jet_DR_Electron",0.2);
    m_e_DR_jet  = conf.getNum(m_name+".Electron_DR_Jet",0.4);
    m_mu_DR_y   = conf.getNum(m_name+".Muon_DR_Photon",0.4);
    m_mu_DR_jet = conf.getNum(m_name+".Muon_DR_Jet",0.4);
    
    return EL::StatusCode::SUCCESS;
  }
  
  // Remove overlap. The input containers are modified: overlapping elements are removed
  void OverlapRemovalHandler::removeOverlap(xAOD::PhotonContainer &photons,
                                            xAOD::JetContainer &jets,
                                            xAOD::ElectronContainer &elecs,
                                            xAOD::MuonContainer &muons)
  {
    removeOverlap(&photons,&jets,&elecs,&muons);
  }
  
  // Remove overlap. The input containers are modified: overlapping elements are removed
  void OverlapRemovalHandler::removeOverlap(xAOD::PhotonContainer *photons,
                                            xAOD::JetContainer *jets,
                                            xAOD::ElectronContainer *elecs,
                                            xAOD::MuonContainer *muons)
  {
    if (photons==nullptr) HG::fatal("removeOverlap cannot be done without photons!");
    
    // 1. remove electrons overlapping with photons
    if (elecs!=nullptr) removeOverlap(*elecs,*photons,m_e_DR_y);
    
    // 2. jets
    if (jets!=nullptr) {
      
      // 2.a remove jets overlapping with photons
      removeOverlap(*jets,*photons,m_jet_DR_y);

      // 2.b remove jets overlapping with electrons
      if (elecs!=nullptr) removeOverlap(*jets,*elecs,m_jet_DR_e);
    }
    
    // 3. remove electrons too close to jets (usually 0.4)
    if ( jets!=nullptr && elecs!=nullptr )
      removeOverlap(*elecs,*jets,m_e_DR_jet);
    
    // 4. remove muons overlapping photons and jets
    if (muons!=nullptr) {
      
      removeOverlap(*muons,*photons,m_mu_DR_y);
      if (jets!=nullptr) removeOverlap(*muons,*jets,m_mu_DR_jet);
    }
    
  }
  
  xAOD::MuonContainer OverlapRemovalHandler::muonsInJets(xAOD::MuonContainer muons,
                                                         xAOD::JetContainer jets,
                                                         double DRcut) {
    // if no DR cut is specified, use the same distance as ued for muons-jet removal
    if (DRcut<0) DRcut=m_mu_DR_jet;
    
    static int nwarn=0; // if a crazy value is used, print a warning the first 20 times
    if (DRcut<0&&nwarn++<20)
      Warning("muonsInJets","The DeltaR value used dosn't make sense: %.3f",DRcut);

    return getOverlaps(muons,jets,DRcut);
  }
  
  
}// namespace HG
