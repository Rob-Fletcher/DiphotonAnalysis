#ifndef HGamAnalysisFramework_ETmissHandler_H
#define HGamAnalysisFramework_ETmissHandler_H

#define __HGamMET__
#include "HGamAnalysisFramework/HgammaHandler.h"
#undef __HGamMET__

//#include "HGamAnalysisFramework/HgammaHandler.h"

// MET EDM
#include "xAODMissingET/MissingET.h"
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"

#include "METUtilities/METMaker.h"
#include "METUtilities/METSystematicsTool.h"

namespace HG {

  class ETmissHandler : public HgammaHandler<xAOD::MissingET, xAOD::MissingETContainer, xAOD::MissingETAuxContainer> {

  private:
    met::METMaker           *m_metMakerTool; //!
    met::METSystematicsTool *m_metSysTool; //!

    std::string m_assocMapName;
    std::string m_coreName;
    std::vector<TString> m_metTypes;

  private:
    // Hide the function which doesn't inclue objects
    virtual xAOD::MissingETContainer getCorrectedContainer();

  public:
    ETmissHandler(const char *name, xAOD::TEvent *event, xAOD::TStore *store);
    virtual ~ETmissHandler();

    virtual EL::StatusCode initialize(Config &config);

    virtual xAOD::MissingETContainer getCorrectedContainer(const xAOD::PhotonContainer   *photons  ,
                                                           const xAOD::JetContainer      *jets     ,
                                                           const xAOD::ElectronContainer *electrons,
                                                           const xAOD::MuonContainer     *muons    );
    virtual xAOD::MissingETContainer applySelection(xAOD::MissingETContainer &container);
    virtual CP::SystematicCode applySystematicVariation(const CP::SystematicSet &sys);

  };
}

#endif // HGamAnalysisFramework_ETmissHandler_H
