#ifndef HGamAnalysisFramework_ElectronHandler
#define HGamAnalysisFramework_ElectronHandler

#include "HGamAnalysisFramework/HgammaHandler.h"

namespace CP {
  class IsolationSelectionTool;
}

namespace HG {
  class ElectronHandler : public HgammaHandler<xAOD::Electron, xAOD::ElectronContainer, xAOD::ElectronAuxContainer> {
  private:
    CP::EgammaCalibrationAndSmearingTool *m_electronCalibTool;
    AsgElectronLikelihoodTool            *m_electronLHIDTool;
    std::map<TString, AsgElectronLikelihoodTool*> m_electronSelectors;
    std::map<TString, SG::AuxElement::Accessor<char>* > m_electronSelDecorators;

    std::map<HG::Iso::IsolationType, CP::IsolationSelectionTool*> m_isoTools;
    std::map<HG::Iso::IsolationType, SG::AuxElement::Accessor<char>* > m_isoDecorators;

    AsgElectronEfficiencyCorrectionTool  *m_electronIDSF;
    AsgElectronEfficiencyCorrectionTool  *m_electronRecoSF;

    bool    m_doPidCut;
    StrV    m_pidCuts;
    TString m_defaultPid;
    StrV    m_pidConfigs;    
    bool    m_doIsoCut;
    StrV    m_isoCuts;
    HG::Iso::IsolationType m_defaultIso;
    double  m_etaCut;
    double  m_ptCut;

    bool    m_crackReject;
    double  m_barrelMax;
    double  m_endcapMin;

    bool   m_applyIPCuts;
    double m_d0BySigd0Cut;
    double m_z0Cut;


  public:
    static SG::AuxElement::Decorator<float> effIDSF;
    static SG::AuxElement::Decorator<float> effRecoSF;
    static SG::AuxElement::Decorator<float> scaleFactor;
    static SG::AuxElement::Decorator<float> Ecalib_ratio, Ereso;
    static SG::AuxElement::Decorator<char>  passIPCut;
    static SG::AuxElement::Decorator<char>  isTight, isMedium, isLoose;

    
  public:
    /// constructor
    ElectronHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store);

    /// destructor
    virtual ~ElectronHandler();

    virtual EL::StatusCode initialize(Config &config);

    virtual xAOD::ElectronContainer getCorrectedContainer();
    virtual xAOD::ElectronContainer applySelection(xAOD::ElectronContainer &container);
    virtual CP::SystematicCode    applySystematicVariation(const CP::SystematicSet &sys);

    /// applies kinematic preselection cuts: not-in-crack + pT cut
    bool passPtEtaCuts(const xAOD::Electron *ele);

    /// IP cuts
    void decorateIPCut(xAOD::Electron &ele);
    bool passIPCuts(const xAOD::Electron *ele);
    
    /// applies PID cut
    bool passPIDCut(const xAOD::Electron *ele, TString pid = "Default");
    void decoratePID(xAOD::Electron &ele);

    /// applies Iso cut specified in config file
    bool passIsoCut(const xAOD::Electron *ele, HG::Iso::IsolationType iso = HG::Iso::Undefined);
    void decorateIso(xAOD::Electron &ele);
    
    /// calibrates and smears an electron
    static void calibrateAndSmearElectron(xAOD::Electron *ele,
                                          const xAOD::EventInfo *evtInfo,
                                          CP::EgammaCalibrationAndSmearingTool *electronCalibTool);
    /// access the isolation types needed to initalize the tools
    HG::Iso::IsolationType getIsoType(TString isoName);

    /// decorate electron with efficiency scale factor and uncertainty
    void applyScaleFactor(xAOD::Electron *ele, const xAOD::EventInfo *evtInfo);

    
  };
}

#endif // HGamAnalysisFramework_ElectronHandler
