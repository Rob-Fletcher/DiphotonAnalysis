#include "HGamAnalysisFramework/HGamCategories.h"

#include "HGamAnalysisFramework/VarHandler.h"

namespace HG {

  std::pair<int, float> getCategoryAndWeight(const xAOD::PhotonContainer   *photons  ,
                                             const xAOD::ElectronContainer *electrons,
                                             const xAOD::MuonContainer     *muons    ,
                                             const xAOD::JetContainer      *jets     )
  {
    // Not yet implemented...
    return std::make_pair(1, 1.0);
  }

} // namespace HG
