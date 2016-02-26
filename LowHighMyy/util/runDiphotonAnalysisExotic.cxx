#include "LowHighMyy/DiphotonAnalysis.h"
#include "HGamAnalysisFramework/RunUtils.h"

int main(int argc, char *argv[])
{
  // Set up the job for xAOD access
  xAOD::Init().ignore();

  // Create our algorithm
  DiphotonAnalysis *alg = new DiphotonAnalysis("DiphotonAnalysis");

  //Exotics
  //alg->AnalysisBranch=11;
  alg->AnalysisBranch=11;

  //Hardcode selection parameters
  alg->isolation_cut = 8.*1E3;
  alg->leading_min_pt = 50.*1E3;  
  alg->subleading_min_pt = 50.*1E3; 
    
  // Use helper to start the job
  HG::runJob(alg, argc, argv);

  return 0;
}
