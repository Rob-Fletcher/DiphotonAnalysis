#ifndef HGamAnalysisFramework_HGamCategories_H
#define HGamAnalysisFramework_HGamCategories_H

#include "HGamAnalysisFramework/HgammaIncludes.h"
#include <utility>

namespace HG {
  // Get category for default coupling analysis
  std::pair<int, float> getCategoryAndWeight(const xAOD::PhotonContainer   *photons   = nullptr,
                                             const xAOD::ElectronContainer *electrons = nullptr,
                                             const xAOD::MuonContainer     *muons     = nullptr,
                                             const xAOD::JetContainer      *jets      = nullptr);
}
#endif // HGamAnalysisFramework_HGamCategories_H
