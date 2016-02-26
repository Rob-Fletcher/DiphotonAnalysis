#ifndef HGamAnalysisFramework_PhotonHandler
#define HGamAnalysisFramework_PhotonHandler

#include "HGamAnalysisFramework/HgammaHandler.h"

namespace CP {
  class PhotonPointingTool;
  class IsolationSelectionTool;
}

class ElectronPhotonShowerShapeFudgeTool;
class EGammaAmbiguityTool;

namespace HG {
  class PhotonHandler : public HgammaHandler<xAOD::Photon, xAOD::PhotonContainer, xAOD::PhotonAuxContainer> {
  private:
    CP::EgammaCalibrationAndSmearingTool *m_photonCalibTool;
    CP::PhotonPointingTool               *m_pointTool;
    std::map<egammaPID::PID, AsgPhotonIsEMSelector*> m_photonSelectors;
    std::map<egammaPID::PID, SG::AuxElement::Accessor<char>* > m_pidAcc;
    
    std::map<HG::Iso::IsolationType, CP::IsolationSelectionTool*> m_isoTools;
    std::map<HG::Iso::IsolationType, SG::AuxElement::Accessor<char>* > m_isoAcc;

    ElectronPhotonShowerShapeFudgeTool   *m_fudgeMC;
    EGammaAmbiguityTool                  *m_fakeTool;
    AsgPhotonEfficiencyCorrectionTool    *m_photonSF;

    bool           m_doPidCut;
    StrV           m_pidCuts;
    egammaPID::PID m_defaultPid;
    bool           m_doQuality;
    bool           m_doIsoCut;
    bool           m_doAmbCut;
    StrV           m_isoCuts;
    HG::Iso::IsolationType m_defaultIso;
    double         m_etaCut;
    double         m_ptCut;
    int            m_fudgeSet;

    bool    m_crackReject;
    double  m_barrelMax;
    double  m_endcapMin;

    bool    m_doCalib;
    bool    m_doFudge;
    bool    m_doScale;
    bool    m_doAuthor;
    bool    m_correctIso;
    bool    m_isAFII;



  private:
    void removeAuthorCut(xAOD::PhotonContainer &photons);



  public:
    static SG::AuxElement::Decorator<float> effSF, scaleFactor;
    static SG::AuxElement::Decorator<float> effSFunc;
    static SG::AuxElement::Decorator<float> Ecalib_ratio;
    static SG::AuxElement::Decorator<float> r_SL1, z_SL1;
    static SG::AuxElement::Decorator<float> etaS1, etaS2, cl_eta, cl_phi;
    static SG::AuxElement::Decorator<float> relEreso;
    static SG::AuxElement::Decorator<char> isLoose, isTight;
    static SG::AuxElement::Decorator<char> isConv, passOQ;


  public:
    PhotonHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store);
    virtual ~PhotonHandler();

    virtual EL::StatusCode initialize(Config &config);

    virtual xAOD::PhotonContainer getCorrectedContainer();
    virtual xAOD::PhotonContainer applySelection(xAOD::PhotonContainer &container);
    virtual xAOD::PhotonContainer applyPreSelection(xAOD::PhotonContainer &container);
    virtual CP::SystematicCode    applySystematicVariation(const CP::SystematicSet &sys);

    /// applies kinematic selection cuts: not-in-crack + pT cut
    bool passPtEtaCuts(const xAOD::Photon *gam);

    /// applies OQ cut, if specified
    bool passOQCut(const xAOD::Photon *gam);
    void decorateOQ(xAOD::Photon &gam);

    /// applies ambiguity cut, if specified
    bool passAmbCut(const xAOD::Photon *gam);
    void decorateAmbCut(xAOD::Photon &gam);

    /// applies Iso cut specified in config file
    bool passIsoCut(const xAOD::Photon *gam, HG::Iso::IsolationType iso = HG::Iso::Undefined);
    void decorateIso(xAOD::Photon &gam);

    /// applies PID cut specified in config file
    bool passPIDCut(const xAOD::Photon *gam, egammaPID::PID pid = egammaPID::LastEgammaPID);
    void decoratePID(xAOD::Photon &gam);

    /// Requires photon to pass preselection, defined by
    /// a) OQ, b) pT>25, |eta_s2| selection, and c) Loose PID
    bool passPreSelection(const xAOD::Photon *gam);
    
    /// applies author cut, to remove topocluster seeded photons
    bool passAuthorCut(const xAOD::Photon *gam);

    /// calibrates and smears a photon
    static void calibrateAndSmearPhoton(xAOD::Photon *gam, 
                                 const xAOD::EventInfo *evtInfo,
                                 CP::EgammaCalibrationAndSmearingTool *photonCalibTool);

    /// access the isolation and PID types needed to initalize the tools
    HG::Iso::IsolationType getIsoType(TString isoName);
    unsigned int getPIDmask(TString pidName);
    egammaPID::PID getPID(TString pidName);
    
    /// Use NN to find PV vertex of two leading photons
    double getPVz(xAOD::PhotonContainer &container, const xAOD::EventInfo *eventInfo);

    /// applies PVz correction
    virtual void correctPrimaryVertex(const xAOD::Vertex *vertex, xAOD::Photon &gam);
    
    /// if MC, the photon is smeared using randSeed = 100*evtNum+photonIndex
    /// and the shower shape variables are "fudged"
    void applyFudgeFactor(xAOD::Photon *gam, const xAOD::EventInfo *evtInfo);
    
    /// decorate photon with efficiency scale factor and uncertainty
    void applyScaleFactor(xAOD::Photon *gam, const xAOD::EventInfo *evtInfo);
    
    /// access of the PVz corrected eta of a photon based on the eta in the first sampling
    static double PVz_corrected_eta(double eta_s1, double PVz);

    /// calculates transverse distance [mm] in first sampling layer given the eta coordinate
    static double BarrelR_s1(double eta_s1);

    /// calculates endcap z position [mm] in the first sampling layer given the eta coordinate
    static double EndcapZ_s1(double eta_s1);
    
    /// calculates transverse distance [mm] in first sampling layer given the eta coordinate
    static double r_s1(double eta_s1);
    
    /// calculates z position [mm] in the first sampling layer
    static double z_s1(double eta_s1);

    /// Fernando's method
    static const xAOD::Photon *findTrigMatchedPhoton(const xAOD::Photon *photon,
                                              const xAOD::PhotonContainer *trigPhotons,
                                              bool debug=false);

    /// print details about photon to the screen
    static void printPhoton(const xAOD::Photon *gam, TString comment="");
    
  };
}

#endif // HGamAnalysisFramework_PhotonHandler
