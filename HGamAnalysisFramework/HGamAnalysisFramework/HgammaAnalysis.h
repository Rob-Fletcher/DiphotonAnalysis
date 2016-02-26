#ifndef HGamAnalysisFramework_HgammaAnalysis_H
#define HGamAnalysisFramework_HgammaAnalysis_H

#include <EventLoop/Algorithm.h>

#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamAnalysisFramework/HgammaUtils.h"
#include "HGamAnalysisFramework/HistogramStore.h"
#include "HGamAnalysisFramework/Config.h"

#ifndef __CINT__
#include "PhotonVertexSelection/PhotonVertexSelectionTool.h"
#include "HGamAnalysisFramework/PhotonHandler.h"
#include "HGamAnalysisFramework/ElectronHandler.h"
#include "HGamAnalysisFramework/JetHandler.h"
#include "HGamAnalysisFramework/MuonHandler.h"
#include "HGamAnalysisFramework/EventHandler.h"
#include "HGamAnalysisFramework/TruthHandler.h"
#include "HGamAnalysisFramework/OverlapRemovalHandler.h"
#include "HGamAnalysisFramework/ETmissHandler.h"
#endif

// Forward declarations

class HgammaAnalysis : public EL::Algorithm
{
  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  
  // Init
  //  EL::StatusCode initTools();
  
  //! \brief event counter
  int m_eventCounter;//!
  //! \brief start time (used to calculate the processing rate)
  long m_startTime; //!

  enum class TrigType : char {
      SinglePhoton,
      DiPhoton,
      SingleMuon,
      DiMuon,
      SingleElectron,
      DiElectron,
      Undefined
  };


  
private:
  HG::Config m_config;

  xAOD::TEvent *m_event; //!
  xAOD::TStore *m_store; //!
  //! \brief histogram store
  HistogramStore *m_histoStore; //!
  //! \brief configuration file
  HG::SystematicList m_sysList; //!
  TString m_name;
#ifndef __CINT__
  CP::PhotonVertexSelectionTool *m_vertexTool; //!
  HG::PhotonHandler *m_photonHandler; //!
  HG::ElectronHandler *m_electronHandler; //!
  HG::JetHandler *m_jetHandler; //!
  HG::MuonHandler *m_muonHandler; //!
  HG::EventHandler *m_eventHandler; //!
  HG::TruthHandler *m_truthHandler; //!
  HG::OverlapRemovalHandler *m_overlapHandler; //!
  HG::ETmissHandler *m_etmissHandler; //!
#endif // __CINT__
  bool   m_isInit;
  bool   m_isAOD;
  bool   m_isDAOD;
  bool   m_isMAOD;
  bool   m_isMC;
  bool   m_isData;
  bool   m_doRelPtCut, m_doMyyCut;
  bool   m_doTwoGoodPhotonsCut;
  double m_relPtCut1;
  double m_relPtCut2;
  double m_myyLow, m_myyHigh;
  bool   m_doVertex;
  bool   m_doHardPV;
  bool   m_doJetClean;
  double m_jetCleanPt;
  double m_jetCleanJvt;
  bool   m_doTrigMatch;
  StrV   m_requiredTriggers;
  std::map<TString, TrigType> m_trigMatch;
  std::map<int, double>  m_crossSections, m_genEffs, m_kFactors, m_NevtsInitial;
  std::map<int, double>  m_skimmingEff, m_weightXsec, m_higgsMass;
  std::map<int, TString> m_mcNames;
  std::map<int, int> m_nTotEvents;

protected:
  virtual const HG::SystematicList& getSystematics() { return m_sysList; }
  virtual CP::SystematicCode applySystematicVariation(const CP::SystematicSet &sys);

  // access to useful pointers
  inline virtual xAOD::TEvent*      event() { return m_event; }
  inline virtual xAOD::TStore*      store() { return m_store; }
  inline virtual const xAOD::EventInfo*  eventInfo();

  inline virtual HistogramStore*      histoStore()      { return m_histoStore; }
  inline virtual HG::Config*          config()          { return &m_config; }
  inline virtual HG::PhotonHandler*   photonHandler()   { return m_photonHandler; }
  inline virtual HG::ElectronHandler* electronHandler() { return m_electronHandler; }
  inline virtual HG::JetHandler*      jetHandler()      { return m_jetHandler; }
  inline virtual HG::MuonHandler*     muonHandler()     { return m_muonHandler; }
  inline virtual HG::EventHandler*    eventHandler()    { return m_eventHandler; }
  inline virtual HG::TruthHandler*    truthHandler()    { return m_truthHandler; }
  inline virtual HG::OverlapRemovalHandler* overlapHandler() { return m_overlapHandler; }
  inline virtual HG::ETmissHandler*   etmissHandler()   { return m_etmissHandler; }

  virtual void setWeightInitial();

  virtual double weight();
  virtual double weightInitial();
  inline virtual double weightFinal();
  virtual double weightCategory();

  TString getMCSampleName(int mcChannelNumber = -1);
  int getNtotalEvents(int mcChannelNumber = -1);
  TH1F *getCutFlowHistogram(int mcID, TString suffix="");
  virtual double getCrossSection(int mcChannelNumber = -1);
  virtual double getGeneratorEfficiency(int mcChannelNumber = -1);
  virtual double getGeneratorHiggsMass(int mcChannelNumber = -1);
  virtual double getKFactor(int mcChannelNumber = -1);
  virtual double getIntialSumOfWeights(int mcChannelNumber = -1);

  /// \brief get intL * sigma / Ninitial
  double lumiXsecWeight(double intLumiPbInv = -1, int mcChannelNumber = -1, bool printFirst = true);

  virtual void selectVertex();

  virtual void setSelectedTruthObjects(const xAOD::TruthParticleContainer *photons   = nullptr,
                                       const xAOD::TruthParticleContainer *electrons = nullptr,
                                       const xAOD::TruthParticleContainer *muons     = nullptr,
                                       const xAOD::JetContainer           *jets      = nullptr);
  virtual void setSelectedObjects(const xAOD::PhotonContainer   *photons   = nullptr,
                                  const xAOD::ElectronContainer *electrons = nullptr,
                                  const xAOD::MuonContainer     *muons     = nullptr,
                                  const xAOD::JetContainer      *jets      = nullptr);
  virtual bool pass(const xAOD::PhotonContainer *photons,
                    const xAOD::ElectronContainer *electrons = nullptr,
                    const xAOD::MuonContainer *muons = nullptr,
                    const xAOD::JetContainer *jets = nullptr);
  virtual bool passTwoGoodPhotonsCut(const xAOD::PhotonContainer &photons);
  virtual bool passRelativePtCuts(const xAOD::PhotonContainer &photons);
  virtual bool passMyyWindowCut(const xAOD::PhotonContainer &photons);
  virtual bool passJetCleaning();
  virtual bool passTriggerMatch(const xAOD::PhotonContainer *photons,
                                const xAOD::ElectronContainer *electrons = nullptr,
                                const xAOD::MuonContainer *muons = nullptr,
                                const xAOD::JetContainer *jets = nullptr);
  virtual bool passTriggerMatch(const TString &trig,
                                const xAOD::PhotonContainer *photons,
                                const xAOD::ElectronContainer *electrons = nullptr,
                                const xAOD::MuonContainer *muons = nullptr,
                                const xAOD::JetContainer *jets = nullptr);
  virtual TrigType getTriggerType(TString Trigger);

  inline bool isAOD();
  inline bool isDAOD();
  inline bool isMAOD();
  inline bool isMC();
  inline bool isData();

public:
  
  // this is a standard constructor/destructor
  HgammaAnalysis() : m_name() { }
  HgammaAnalysis(const char *name);
  virtual ~HgammaAnalysis() { }

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  virtual EL::StatusCode createOutput();

  inline virtual void setConfig(const HG::Config &config);
  
  // this is needed to distribute the algorithm to the workers
  ClassDef(HgammaAnalysis, 1);
};

inline
void HgammaAnalysis::setConfig(const HG::Config &config)
{
  m_config = config;
}

inline
bool HgammaAnalysis::isAOD()
{ return m_isAOD; }

inline
bool HgammaAnalysis::isDAOD()
{ return m_isDAOD; }

inline
bool HgammaAnalysis::isMAOD()
{ return m_isMAOD; }

inline
bool HgammaAnalysis::isMC()
{ return m_isMC; }

inline
bool HgammaAnalysis::isData()
{ return m_isData; }

inline double HgammaAnalysis::weightFinal() { return lumiXsecWeight() * weight(); }

inline
const xAOD::EventInfo* HgammaAnalysis::eventInfo()
{
  const xAOD::EventInfo *eventinfo = 0;
  if (m_event->retrieve(eventinfo, "EventInfo").isFailure()) {
    HG::fatal("Cannot access EventInfo");
  }
  return eventinfo;
}
#endif
