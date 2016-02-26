#ifndef HGamAnalysisFramework_MuonHandler
#define HGamAnalysisFramework_MuonHandler

#include "HGamAnalysisFramework/HgammaHandler.h"

namespace CP {
  class IsolationSelectionTool;
}

namespace HG {
  class MuonHandler : public HgammaHandler<xAOD::Muon, xAOD::MuonContainer, xAOD::MuonAuxContainer> {
  private:
    CP::MuonEfficiencyScaleFactors     *m_muonEffScaleFactors;
    CP::MuonTriggerScaleFactors        *m_muonTrigScaleFactors;
    CP::MuonCalibrationAndSmearingTool *m_muonCalibTool;
    CP::MuonSelectionTool              *m_muonSelectTool;
    std::map<HG::Iso::IsolationType, CP::IsolationSelectionTool*> m_isoTools;
    std::map<HG::Iso::IsolationType, SG::AuxElement::Accessor<char>* > m_isoDecorators;


    TString m_pidCut;
    
    bool   m_ApplyPtCut;
    double m_PtCut;
    
    bool   m_ApplyIPCuts;
    double m_d0Cut;
    double m_z0Cut;
    
    bool m_doIsoCut;
    StrV m_isoCuts;
    HG::Iso::IsolationType m_defaultIso;
    
    double m_MaxEta;
    
  public:
    static SG::AuxElement::Decorator<float> effSF;
    static SG::AuxElement::Decorator<float> scaleFactor;
    static SG::AuxElement::Decorator<char> passIPCut;
    static SG::AuxElement::Decorator<char> isAccepted;
    
  public:
    MuonHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store);
    virtual ~MuonHandler();
    
    virtual EL::StatusCode initialize(Config &config);
    
    virtual xAOD::MuonContainer getCorrectedContainer();
    virtual xAOD::MuonContainer applySelection(xAOD::MuonContainer &container);
    virtual CP::SystematicCode  applySystematicVariation(const CP::SystematicSet &sys);
    
    /// applies kinematic preselection cuts: not-in-crack + pT cut
    bool passSelection(const xAOD::Muon *muon);
    
    /// applies PID cut specified in config file
    bool passIsoCut(const xAOD::Muon *muon, HG::Iso::IsolationType iso = HG::Iso::Undefined);
    void decorateIso(xAOD::Muon &muon);

    /// Decorates IP cuts result to muon
    void decorateIPCut(const xAOD::Muon *muon);
    
    /// calibrates and smears a photon
    static void calibrateAndSmearMuon      (xAOD::Muon *muon, CP::MuonCalibrationAndSmearingTool *muonCalibTool);

    /// decorate photon with efficiency scale factor and uncertainty
    void applyScaleFactor(xAOD::Muon *muon);
    
    /// print details about photon to the screen
    static void printMuon(const xAOD::Muon *muon, TString comment="");
    
    
  };
}

#endif // HgammaSandbox_PhotonHandler
