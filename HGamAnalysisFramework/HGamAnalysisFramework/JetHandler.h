#ifndef HGamAnalysisFramework_JetHandler_H
#define HGamAnalysisFramework_JetHandler_H

#include "HGamAnalysisFramework/HgammaHandler.h"

class JetUncertaintiesTool;
class JERSmearingTool;
class JERTool;
class JetVertexTaggerTool;
namespace InDet {
  class InDetTrackSelectionTool;
}
class JetOriginCorrectionTool;
class BTaggingEfficiencyTool;

namespace HG {
  class JetHandler : public HgammaHandler<xAOD::Jet, xAOD::JetContainer, xAOD::JetAuxContainer> {
  private:
    JetCalibrationTool     *m_jetCalibTool; //!
    JetCleaningTool        *m_jetCleaning; //!
    JetUncertaintiesTool   *m_jesProvider; //!
    JERTool                *m_jerTool; //!
    JERSmearingTool        *m_jerSmear; //!
    InDet::InDetTrackSelectionTool *m_trackTool; //!
    TH2F                   *m_jvtLikelihood; //!
    JetVertexTaggerTool    *m_jvtTool; //!
    JetOriginCorrectionTool *m_jetOriginTool; //!
    StrV                    m_bTagNames; //!
    std::map<TString, StrV> m_bTagEffs; //!
    std::map<TString, StrV> m_bTagOPs; //!
    std::map<TString, NumV> m_bTagCuts; //!
    std::map<TString, std::vector<BTaggingEfficiencyTool*> > m_bTagEffTools; //!

    double  m_rapidityCut;
    double  m_ptCut;
    double  m_jvf;
    double  m_jvt;
    double  m_bTagRapidityCut;
    double  m_bTagJvfCut;
    bool    m_enableBTagging;
    bool    m_isAOD;
    bool    m_isAFII;
    bool    m_correctVertex;



  public:
    static SG::AuxElement::Accessor<std::vector<float> > JVF;
    static SG::AuxElement::Accessor<float>  Jvt;
    static SG::AuxElement::Decorator<float> Jvf;
    static SG::AuxElement::Decorator<float> scaleFactor;
    static SG::AuxElement::Decorator<char>  isClean;
    static SG::AuxElement::Accessor<float>  DetectorEta;



  public:
    JetHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store);
    virtual ~JetHandler();

    virtual EL::StatusCode initialize(Config &config);

    /// Set if it's an xAOD (true) or a DxAOD (false)
    virtual void setAOD(bool flag) { m_isAOD = flag; }

    virtual xAOD::JetContainer getCorrectedContainer();
    virtual xAOD::JetContainer applySelection(xAOD::JetContainer &container);
    virtual CP::SystematicCode applySystematicVariation(const CP::SystematicSet &sys);

    /// applies kinematic preselection cuts: not-in-crack + pT cut
    bool passPtEtaCuts(const xAOD::Jet *jet);

    /// corrects JVT value
    void rescaleJVT(xAOD::Jet &jet);
    void recalculateJVT(xAOD::Jet &jet);

    /// applies kinematic preselection cuts: not-in-crack + pT cut
    static bool passJVTCut(const xAOD::Jet *jet, float jvtCut);

    /// applies kinematic preselection cuts: not-in-crack + pT cut
    bool passJVFCut(const xAOD::Jet *jet, bool useBTagCut = false);

    /// calibrates and smears a jet
    void calibrateJet(xAOD::Jet *jet);

    /// apply btags to jets
    void decorateBJet(xAOD::Jet *jet);

    /// apply btags to jets
    void decorateJVF(xAOD::Jet *jet);

    /// print details about photon to the screen
    static void printJet(const xAOD::Jet *gam, TString comment="");
    
  };
}

#endif // HGamAnalysisFramework_JetHandler_H
