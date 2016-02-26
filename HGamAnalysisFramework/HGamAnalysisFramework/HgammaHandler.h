#ifndef HGamAnalysisFramework_HgammaHandler_H
#define HGamAnalysisFramework_HgammaHandler_H

#include "HGamAnalysisFramework/HgammaIncludes.h"

#include "IsolationTool/TrackIsolationTool.h"

// Used HSG5 ObjectHandler as example, then simplified

namespace HG {

  template <class partType, class partContainer, class auxContainer>
  class HgammaHandler {
  protected:
    TString       m_name;
    xAOD::TEvent *m_event;
    xAOD::TStore *m_store;

    xAOD::TrackIsolationTool              *m_trackIsoTool;
    std::vector<xAOD::Iso::IsolationType>  m_isoT;
    xAOD::TrackCorrection                  m_corrList;

    TString       m_sysName;
    TString       m_containerName;
    TString       m_truthName;
    TString       m_MxAODname;
    bool          m_isMxAOD;
    bool          m_isData;
    bool          m_isMC;
    bool          m_isVtxCorr;



  protected:
    // Old method
    partContainer getShallowContainer(std::string name = ""); 
    // New methods
    partContainer getShallowContainer(bool &calibStatus, bool makeShallowCopy = true); 
    partContainer getStoreContainer(std::string name); 
    partContainer getEventContainer(std::string name, std::string post = "", bool makeShallowCopy = true); 

    bool          isMC();

    virtual void correctContainerPV(partContainer &cont);
    virtual void correctPrimaryVertex(const xAOD::Vertex *vertex, partType &part);
    


  public:
    HgammaHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store);
    virtual ~HgammaHandler();
     
    virtual EL::StatusCode initialize(Config &config);

    partContainer              getContainer(std::string name);
    virtual partContainer      getCorrectedContainer() = 0;
    virtual partContainer      applySelection(partContainer &container) = 0;
    virtual CP::SystematicCode applySystematicVariation(const CP::SystematicSet &sys) = 0;
    virtual EL::StatusCode     writeContainer(partContainer &container, TString name = "");
    virtual EL::StatusCode     writeTruthContainer(partContainer &container, TString name = "");

    static bool comparePt(const partType *a, const partType *b);

    virtual void setVertexCorrected(bool flag = true) { m_isVtxCorr = flag; }



  };

}

#include "HGamAnalysisFramework/HgammaHandler.hpp"

#endif // HGamAnalysisFramework_HgammaHandler_H
