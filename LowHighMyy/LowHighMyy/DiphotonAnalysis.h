//***********
// Analysis for diphoton used in the Exotic diphoton analysis and Higgs low-high mass search
// Supported by: Simone Mazza <simone.mazza@mi.infn.it>
//               Kirill Grevtsov <kirill.grevtsov@cern.ch>
//***********


#ifndef DiphotonAnalysis_DiphotonAnalysis_H
#define DiphotonAnalysis_DiphotonAnalysis_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"
#include "LowHighMyy/CutFlowHisto.h"

#include <EventLoop/Algorithm.h>
#include <EventLoop/OutputStream.h>
#include <EventLoopAlgs/NTupleSvc.h>

// SUSY tools includes
#include "SUSYTools/SUSYCrossSection.h"

#include <TParameter.h>

// header for systematics:
#include "PATInterfaces/SystematicRegistry.h"

// Forward declarations
class GoodRunsListSelectionTool;
class AsgElectronIsEMSelector;
class AsgPhotonIsEMSelector;
class AsgFudgeMCTool;
namespace CP {
  class EgammaCalibrationAndSmearingTool;
  class IsolationCorrectionTool;
  class PhotonVertexSelectionTool;
}

namespace CP{
  class IsolationSelectionTool;
}

namespace Trig {
    class TrigDecisionTool;
}
namespace TrigConf {
    class xAODConfigTool;
    class ITrigConfigTool;
}

class DiphotonAnalysis : public HgammaAnalysis
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;

  #ifndef __CINT__
    AsgPhotonIsEMSelector* photonID_tight_tool_SS; //!
    AsgFudgeMCTool* m_fudgeMCTool; //!
    CP::IsolationCorrectionTool* isoCorr_tool; //!
    CP::PhotonVertexSelectionTool* m_phVerSel_tool; //!
    std::map<std::string,CP::IsolationSelectionTool*> m_isoTool; //!

    // For triggers
    TrigConf::xAODConfigTool* configTool; //!
    Trig::TrigDecisionTool* trigDecTool; //!
    ToolHandle<TrigConf::ITrigConfigTool>* configHandle; //!

  #endif // not __CINT__
  // configurable tools

  //Type of Analysis
  int AnalysisBranch;

  // SUSYTool Cross-section
  SUSY::CrossSectionDB *my_XsecDB;  //!

  // configuration variables
  unsigned int min_nphotons;
  unsigned int min_nvertex;
  unsigned int min_ntracks;


  /*          Isolation cuts                */
  float isolation_cut;       // both for Exotics and Hohhs (different values)
  float isolation_track_cut; // for Higgs only

  /*          E_T cuts                */
  float leading_min_pt;      // p_T cut for Exotics
  float subleading_min_pt;   // p_T cut for Exotics

  float leading_rel_cut_pt;      // E_T/m_gg cut for Higgs
  float subleading_rel_cut_pt;   // E_T/m_gg cut for Higgs

  int fudge_tune;
  int loose_tune;
  int tight_tune;

  bool save_all_photons;
  bool do_systematics;
  bool correct_isolation;
  bool save_all_events;
  bool save_all_preselection;

  bool test_bkg;
  bool no_xaod;

  // output stream names
  std::string outputStreamName;

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
private:
  // Tree *myTree; //!
  // TH1 *myHist; //!

  TH1 *nph_hist; //!
  TH1 *pileup_hist; //!
  CutFlowHisto* cutflow; //!

  // weights
  float weight_pileup; //!
  TParameter<float>* weight_sum_before; //!
  float weight_sum_before_f; //!
  TParameter<float>* weight_sum_selected; //!

  // other variables
  float averageIntPerXing; //!
  float actualIntPerXing; //!
  int run_number; //!
  unsigned long long event_number; //!

  int npvs; //!
  float xs; //!
  float xs_ami; //!
  float filter_eff; //!
  float filter_eff_ami; //!
  float pileup_weight; //!
  float MC_weight; //!
  float xs_weight; //!
  float event_weight; //!
  float prel_weight; //!

  int total_events; //!
  int non_derived_total_events; //!

  double ED_central; //!
  double ED_forward; //!

  // pass flags
  std::vector< std::pair<std::string, bool> > pass_flags; //!

  // isEM flags
  std::vector< std::pair<std::string, int> > isEM_flags_lead; //!
  std::vector< std::pair<std::string, int> > isEM_flags_sublead; //!

  std::vector< std::vector< std::pair<std::string, int> > > ph_isEM_flags; //!

  // is MonteCarlo
  bool is_mc; //!

  // list of systematics
  std::vector<CP::SystematicSet> sysList; //!

  // cinematic variables
  float mass; //!
  float mass_gev; //!
  float costhetastar; //!
  double deltaphi; //!

  TLorentzVector LV_leading; //!
  TLorentzVector LV_subleading; //!
  TLorentzVector LV_diphoton; //!

  float rho_median; //!
  float pt_subleading; //!
  float pt_leading; //!
  float phi_subleading; //!
  float phi_leading; //!
  float eta_subleading; //!
  float eta_leading; //!
  int conv_subleading; //!
  int conv_leading; //!
  unsigned int isEM_subleading; //!
  unsigned int isEM_leading; //!

  float topoetcone40_leading; //!
  float topoetcone40_rel17_leading; //!

  float topoetcone40_rel17_electron_leading; //!
  float topoetcone40_trouble_electron_leading; //!
  float topoetcone40_electron_leading; //!
  int author_electron_leading; //!

  float topoetcone40_subleading; //!
  float topoetcone20_leading; //!
  float topoetcone20_subleading; //!

  float truth_etcone40_leading; //!
  float truth_etcone40_subleading; //!
  float truth_etcone20_leading; //!
  float truth_etcone20_subleading; //!

  float truth_etcone40_PUcorr_leading; //!
  float truth_etcone40_PUcorr_subleading; //!
  float truth_etcone20_PUcorr_leading; //!
  float truth_etcone20_PUcorr_subleading; //!

  float truth_ptcone40_leading; //!
  float truth_ptcone40_subleading; //!
  float truth_ptcone20_leading; //!
  float truth_ptcone20_subleading; //!

  float truth_ptcone40_PUcorr_leading; //!
  float truth_ptcone40_PUcorr_subleading; //!
  float truth_ptcone20_PUcorr_leading; //!
  float truth_ptcone20_PUcorr_subleading; //!

  float truth_local_etcone40_leading; //!
  float truth_local_etcone40_subleading; //!
  float truth_local_etcone20_leading; //!
  float truth_local_etcone20_subleading; //!

  float truth_rho_central; //!
  float truth_rho_forward; //!

  int loose_leading; //!
  int loose_prime_leading; //!
  int tight_leading; //!
  int my_tight_leading; //!
  int tight_subleading; //!
  int loose_subleading; //!
  int loose_prime_subleading; //!
  bool match_leading; //!
  bool match_subleading; //!
  int parent_pdgID_ld; //!
  int parent_pdgID_subld; //!

  bool bg_truth_match_leading; //!
  bool bg_truth_match_origin_leading; //!

  float Rhad_leading; //!
  float e277_leading; //!
  float Reta_leading; //!
  float Rphi_leading; //!
  float weta2_leading; //!
  float f1_leading; //!
  float DeltaE_leading; //!
  float wtots1_leading; //!
  float weta1_leading; //!
  float fracs1_leading; //!
  float Eratio_leading; //!

  float E1_E2_leading; //!
  float etaS2_leading; //!

  float ptvarcone20_leading; //!
  float ptcone20_leading; //!
  float ptvarcone40_leading; //!
  float ptcone40_leading; //!
  float topoetcone20_Pt_leading; //!
  float ptvarcone20_Pt_leading; //!

  float ptcone20_subleading; //!
  float ptcone40_subleading; //!

  int origin_leading; //!
  int type_leading; //!

  bool two_truth_photons; //!
  bool pass_eta_truth_analysis; //!
  bool truth_leading_matched_leading_ph; //!
  bool truth_subleading_matched_subleading_ph; //!
  TLorentzVector truth_leading_LV; //!
  TLorentzVector truth_subleading_LV; //!
  TLorentzVector truth_diphoton_LV; //!
  const xAOD::TruthParticle* truth_leading_photon; //!
  const xAOD::TruthParticle* truth_subleading_photon; //!

  std::vector<float> ph_pt; //!
  std::vector<float> ph_eta; //!
  std::vector<float> ph_etas2; //!
  std::vector<float> ph_phi; //!
  std::vector<float> ph_cl_phi; //!

  std::vector<bool> ph_tight; //!
  std::vector<bool> ph_matched; //!

  // systematics vectors
  std::vector<float> mass_gev_syst; //!
  std::vector<int> accepted_syst; //!
  std::vector<TLorentzVector> leading_LV_syst; //!
  std::vector<TLorentzVector> subleading_LV_syst; //!
  std::vector<TLorentzVector> diphoton_LV_syst; //!

  // Isolation
  std::vector<float> topoetcone40; //!
  std::vector<float> topoetcone30; //!
  std::vector<float> topoetcone20; //!

  std::vector<float> ptcone40; //!
  std::vector<float> ptcone30; //!
  std::vector<float> ptcone20; //!

  std::vector<float> ptvarcone40; //!
  std::vector<float> ptvarcone30; //!
  std::vector<float> ptvarcone20; //!

  std::vector<int> ph_parent_pdgID; //!

  bool pass_truth_match; //!

  bool FixedCutTightCaloOnly_ld; //!
  bool FixedCutTight_ld; //!
  bool FixedCutLoose_ld; //!

  bool FixedCutTightCaloOnly_subld; //!
  bool FixedCutTight_subld; //!
  bool FixedCutLoose_subld; //!

  std::vector< std::string > ph_ss_names; //!
  std::vector< std::vector<float> > ph_ss; //!

  // output NTUPLES
  EL::NTupleSvc *sample_NTUP; //!

public:
  // this is a standard constructor
  DiphotonAnalysis() { }
  DiphotonAnalysis(const char *name);
  virtual ~DiphotonAnalysis();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode execute();

  // these are all the other functions inherited from the Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(DiphotonAnalysis, 1);

private:

  void clear_variables();

  // weights
  //float compute_pileup_weight();

  bool pass_cut(std::string cut_name, bool pass);

  // custom pass selections
  bool cut_PV();

  template <class T> void make_syst_branches(std::string tag, std::vector<T>& vect);
  bool attach_variables_egamma(xAOD::IParticle* particle) const;
  bool attach_variables_eventinfo(xAOD::EventInfo* event_info);
  bool attach_truth_particle(xAOD::IParticle* particle) const;
  std::vector<double> get_truth_iso(const xAOD::TruthParticle *ptcl, HG::TruthPtcls &stblPtcls);

  EL::StatusCode run_truth_analysis();

};

#endif // DiphotonAnalysis_DiphotonAnalysis_H
