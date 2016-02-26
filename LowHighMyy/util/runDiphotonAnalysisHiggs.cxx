#include "LowHighMyy/DiphotonAnalysis.h"
#include "HGamAnalysisFramework/RunUtils.h"

int main(int argc, char *argv[])
{
  // Set up the job for xAOD access
  xAOD::Init().ignore();

  // Create our algorithm
  DiphotonAnalysis *alg = new DiphotonAnalysis("DiphotonAnalysis");
  
  //Higgs
  alg->AnalysisBranch=1;

  //Hardcode selection parameters
  alg->isolation_cut = 6.*1E3;
  alg->isolation_track_cut = 2.6*1E3;
  alg->leading_rel_cut_pt = 0.4; //Relative cut E_T/m_gg
  alg->subleading_rel_cut_pt = 0.3;
    
  // Use helper to start the job
  HG::runJob(alg, argc, argv);

  return 0;
}
