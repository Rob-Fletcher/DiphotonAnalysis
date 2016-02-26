//#define DEBUG_BUILD true

// Internal includes
#include "LowHighMyy/DiphotonAnalysis.h"
#include "LowHighMyy/debug.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "xAODEventShape/EventShape.h"
#include "xAODCutFlow/CutBookkeeper.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"

// other includes
#include <EventLoop/Worker.h>
#include <EventLoop/Job.h>

// Tools includes
#include "IsolationCorrections/IsolationCorrectionTool.h"
//#include "ElectronIsolationSelection/IsolationSelectionTool.h"
#include "IsolationSelection/IsolationSelectionTool.h"
#include "PhotonVertexSelection/PhotonVertexSelectionTool.h"

// Egamma includes
#include "xAODEgamma/PhotonxAODHelpers.h"
#include "xAODEgamma/EgammaDefs.h"
#include "xAODEgamma/EgammaEnums.h"
#include "xAODEgamma/EgammaxAODHelpers.h"
//#include "xAODEgamma/​EgammaTruthxAODHelpers.h"

#include "xAODTruth/xAODTruthHelpers.h"

// SG includes
#include <xAODTruth/TruthParticle.h>
#include <xAODTruth/TruthParticleContainer.h>

// Trigger includes
//#include "xAODTrigger/"
#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"

// this is needed to distribute the algorithm to the workers
ClassImp(DiphotonAnalysis)

template <class T, class U>
class ShallowContainerDeleter{
  public:
    ShallowContainerDeleter(std::pair< T*, U* > shallowpair){
      m_shallowpair = shallowpair;
    }

    void release(){
      m_shallowpair.first.release();
      m_shallowpair.second.release();
    }

    ~ShallowContainerDeleter(){
      if(m_shallowpair.first) delete m_shallowpair.first;
      if(m_shallowpair.second) delete m_shallowpair.second;
    }

  private:
    std::pair< T*, U* > m_shallowpair;
};

DiphotonAnalysis::DiphotonAnalysis(const char *name)
: HgammaAnalysis(name)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

  // Type of Analysis
    AnalysisBranch = 11; // 0 - not specified (should give error)
                        // 1 - Higgs, 11 - Exotics;

    // cut values

    /*   different for Exotic and Higgs  */
    //should be defind in steering macro corresponding to analysis type
    isolation_cut = -999;
    isolation_track_cut = -999;

    leading_min_pt = -999;
    subleading_min_pt = -999;
    leading_rel_cut_pt = -999;
    subleading_rel_cut_pt = -999;

    /*   same for Exotic and Higgs  */
    // hardcoded here
    min_nvertex = 1;
    min_ntracks = 3;

    // tunes
    fudge_tune = 14;
    loose_tune = 4;
    tight_tune = 2012;

    // configurations
    save_all_photons = false;
    do_systematics = false;
    save_all_events = false;
    save_all_preselection = false;
    correct_isolation = true;
    test_bkg = false;
    no_xaod = false;

}



DiphotonAnalysis::~DiphotonAnalysis()
{
  // Here you delete any memory you allocated during your analysis.
}



EL::StatusCode DiphotonAnalysis::createOutput()
{
  // Here you setup the histograms needed for you analysis. This method
  // gets called after the Handlers are initialized, so that the systematic
  // registry is already filled.

  //histoStore()->createTH1F("m_yy", 60, 110, 140);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode DiphotonAnalysis::execute()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  clear_variables();
  pass_cut("all", true);

  // Important to keep this, so that internal tools / event variables
  // are filled properly.
  HgammaAnalysis::execute();

  //----------------------------
  // Get Event information
  //---------------------------
  const xAOD::EventInfo* eventInfo = 0;
  if( ! event()->retrieve( eventInfo, "EventInfo").isSuccess() ){
    Error("execute()", "Failed to retrieve event info collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // check if the event is data or MC
  // (many tools are applied either to data or MC)
  if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ){
       is_mc = true; // can do something with this later
  }else is_mc = false;

  // ----------------------------------------------------------------
  // run truth analysis
  // ----------------------------------------------------------------

  if(is_mc) run_truth_analysis();

  // ----------------------------------------------------------------
  // weights for MC
  // ----------------------------------------------------------------

  if(is_mc){
    xs = 1E-3 * my_XsecDB->rawxsect(eventInfo->mcChannelNumber());
    filter_eff = my_XsecDB->efficiency(eventInfo->mcChannelNumber());
    xs_ami = wk()->metaData()->castDouble ("xs",1.);
    filter_eff_ami = wk()->metaData()->castDouble ("filter_eff",1.);
    total_events = wk()->metaData()->castDouble ("total_events",1);
    non_derived_total_events = wk()->metaData()->castDouble ("non_derived_total_events",1);

    double t_events = 1;
    double t_xs = 1;
    double t_filter_eff = 1;

    if(non_derived_total_events!=1.) t_events = non_derived_total_events;
    else t_events = total_events;

    if(xs > 0.){
      t_xs = xs;
      t_filter_eff = filter_eff;
    }else{
      t_xs = xs_ami;
      t_filter_eff = filter_eff_ami;
    }

    MC_weight = eventHandler()->mcWeight();
    // Pileup weight
    pileup_weight = eventHandler()->pileupWeight();

    xs_weight = t_xs*t_filter_eff;
    event_weight = MC_weight*pileup_weight;

    //prel_weight = event_weight*xs_weight*1./weight_sum_before_f;
    prel_weight = event_weight*xs_weight*1./t_events;

    /*
    std::cout<<"xsec: "<<xs<<" from ami: "<<xs_ami<<std::endl;
    std::cout<<"GenFiltEff: "<<filter_eff<<" from ami "<<filter_eff_ami<<std::endl;
    std::cout<<"Total Events: "<<total_events<<" non derived "<<non_derived_total_events<<std::endl;
    std::cout<<"MC weight: "<<eventHandler()->mcWeight()<<std::endl;
    std::cout<<"Pileup weight: "<<eventHandler()->pileupWeight()<<std::endl;
    std::cout<<"weight: "<<weight<<std::endl;
    */
  }

  // Create output collections
  xAOD::PhotonContainer* sel_photons       = new xAOD::PhotonContainer();
  xAOD::PhotonAuxContainer* sel_photonsAux = new xAOD::PhotonAuxContainer();
  sel_photons->setStore( sel_photonsAux ); //< Connect the two

  xAOD::EventInfo *my_store       = new xAOD::EventInfo();
  xAOD::AuxInfoBase * my_storeAux = new xAOD::AuxInfoBase();
  my_store->setStore(my_storeAux);

  // Get ED
  const xAOD::EventShape* eventShape_central_xaod = 0;
  if( ! event()->retrieve( eventShape_central_xaod, "TopoClusterIsoCentralEventShape").isSuccess() ){
      Error("execute()", "Failed to retrieve event TopoClusterIsoCentralEventShape collection. Exiting." );
      return EL::StatusCode::FAILURE;
  }
  eventShape_central_xaod->getDensity(xAOD::EventShape::Density, ED_central);

  const xAOD::EventShape* eventShape_forward_xaod = 0;
  if( ! event()->retrieve( eventShape_forward_xaod, "TopoClusterIsoForwardEventShape").isSuccess() ){
    Error("execute()", "Failed to retrieve event TopoClusterIsoForwardEventShape collection. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  eventShape_forward_xaod->getDensity(xAOD::EventShape::Density, ED_forward);

  averageIntPerXing = eventHandler()->mu(); //?? doesn't work with data
  //std::cout << "averageIntPerXing from evthandler: " << averageIntPerXing << "\n";
  averageIntPerXing = eventInfo->averageInteractionsPerCrossing();
  //std::cout << "averageIntPerXing: " << averageIntPerXing << "\n";
  run_number = eventInfo->runNumber();
  event_number = eventInfo->eventNumber();

  // ----------------------------------------------------------------
  // Pass Dalitz
  // ----------------------------------------------------------------
  bool pass_dalitz = true;
  if(is_mc) pass_dalitz = !eventHandler()->isDalitz();
  if( !pass_cut("pass_dalitz", pass_dalitz) ) return EL::StatusCode::SUCCESS;

  // ----------------------------------------------------------------
  // Pass trigger: HLT_2g50_loose (for Graviton analysis)
  // Pass trigger: HLT_g35_medium_g25_medium (for Higgs analysis)
  // ----------------------------------------------------------------
  bool pass_trigger = true;
  if( AnalysisBranch==11 ) pass_trigger = eventHandler()->passTrigger("HLT_2g50_loose");
  if( AnalysisBranch==1 ) pass_trigger = eventHandler()->passTrigger("HLT_g35_medium_g25_medium");
  if( !pass_cut("pass_trigger", pass_trigger) ) return EL::StatusCode::SUCCESS;

  // ----------------------------------------------------------------
  // Pass event selection: GRL+LAr+passTile+passCore
  // ----------------------------------------------------------------
  bool pass_grl = true;
  if(!is_mc) pass_grl = eventHandler()->passGRL(eventInfo);
  if( !pass_cut("pass_grl", pass_grl) ) return EL::StatusCode::SUCCESS;

  bool pass_detector_DQ = true;
  if(!is_mc) pass_detector_DQ =  ( eventHandler()->passLAr (eventInfo) &&
                                   eventHandler()->passTile(eventInfo) &&
                                   eventHandler()->passCore(eventInfo) );
  if( !pass_cut("pass_detector_DQ", pass_detector_DQ) ) return EL::StatusCode::SUCCESS;

  // ----------------------------------------------------------------
  // Pass PV selection
  // ----------------------------------------------------------------

  bool pass_PV = eventHandler()->passVertex(eventInfo);
  //if( !pass_cut("pass_PV", cut_PV() )) return EL::StatusCode::SUCCESS;
  cut_PV(); // Called to get nPV
  if( !pass_cut("pass_PV", pass_PV )) return EL::StatusCode::SUCCESS;


  // *** Loop over recommended systematics
  int sysList_i = 0;
  for (auto sys : sysList){ // auto sysListItr = sysList.begin(); sysListItr != sysList.end(); ++sysListItr){
      bool is_nominal = false;
      if(sys.name()==""){
        DEBUG("Nominal (no syst)");
        is_nominal = true;
      }else{
        if(!do_systematics) break;
        DEBUG("Systematic: " << (sys).name());
      }

      // ************ APPLY SYSTEMATICS HERE **************
      applySystematicVariation(sys);

      // ----------------------------------------------------------------
      // Get photons
      // ----------------------------------------------------------------

      auto Allphotons = photonHandler()->getCorrectedContainer();

      // ----------------------------------------------------------------
      // Apply preselection cuts on photons (using framework)
      // Loose + Et > 25 GeV + eta < 2.37 + no crack
      // ----------------------------------------------------------------

      auto photons = photonHandler()->applySelection(Allphotons);
      bool pass_preselection = (photons.size() >= 2);
      if( !pass_cut( "pass_preselection", pass_preselection) ) return EL::StatusCode::SUCCESS;

      // If there is at least one photon save info about it (working only with --all option)
      if(photons.size() >= 1){
          xAOD::Photon* leading_photon = photons[0];
          //xAOD::Photon* subleading_photon = photons[1];
          xAOD::Photon* subleading_photon = 0;
          if(photons.size() >= 2) subleading_photon = photons[1];
          else  subleading_photon = photons[0];

          origin_leading = xAOD::TruthHelpers::getParticleTruthOrigin(*leading_photon);
          type_leading = xAOD::TruthHelpers::getParticleTruthType(*leading_photon);

          if(test_bkg){
            bg_truth_match_leading = false;
            bg_truth_match_origin_leading = false;
            int nphotons = 0;
            for(auto photon : photons){
              const xAOD::TruthParticle* true_photon = xAOD::TruthHelpers::getTruthParticle(*photon);
              if(true_photon){
                nphotons++;
                if(nphotons > 2) break;
                //if ((origin_leading != 3 && type_leading != 14) &&
                //   (origin_leading != 0 && type_leading != 0) &&
                //   (origin_leading != 0 && type_leading != 13)) bg_truth_match_origin_leading = true;
                origin_leading = xAOD::TruthHelpers::getParticleTruthOrigin(*photon);
				        type_leading = xAOD::TruthHelpers::getParticleTruthType(*photon);
                if (origin_leading != 14 && origin_leading != 38 && origin_leading != 39 && origin_leading != 40){
                   bg_truth_match_origin_leading = true;
                }

                int ph_truth_type = true_photon->pdgId();
                int ph_truth_mothertype = true_photon->parent(0)->pdgId();
                bool from_meson = ( (abs(ph_truth_mothertype) > 100) && (abs(ph_truth_mothertype) != 5100039) );
                bool not_photon = (ph_truth_type != 22);
                bool is_photon = (ph_truth_type == 22);
                bool not_electron = (abs(ph_truth_type) != 11);
                if(not_electron && (not_photon || (is_photon && from_meson)) ){
                  bg_truth_match_leading = true;
                  leading_photon = photon;
                  break;
                }
              }
            }
            if(bg_truth_match_leading != bg_truth_match_origin_leading) std::cout<<"WARNING: no match in pdgID check and truth origin\n";
            if(!(bg_truth_match_leading || bg_truth_match_origin_leading)) return EL::StatusCode::SUCCESS;
            //if(!(bg_truth_match_leading)) return EL::StatusCode::SUCCESS;
          }

          float Et_leading = leading_photon->pt();
          float Et_subleading = subleading_photon->pt();

	        //Define mass (required for Higgs E_T selection)
          TLorentzVector Lead_lv = leading_photon->p4();
          TLorentzVector SubLead_lv = subleading_photon->p4();
          TLorentzVector gamgam_lv = Lead_lv + SubLead_lv;

          LV_diphoton = gamgam_lv;
          LV_leading = Lead_lv;
          LV_subleading = SubLead_lv;

          mass = gamgam_lv.M();

          // ----------------------------------------------------------------
          // Apply leading/subleading tight ID
          // ----------------------------------------------------------------

          bool leading_is_tight(false), subleading_is_tight(false);
          leading_is_tight = photonHandler()->passPIDCut(leading_photon, egammaPID::IsEMTight);
          subleading_is_tight = photonHandler()->passPIDCut(subleading_photon, egammaPID::IsEMTight);

          if( !pass_cut( "pass_ld_subld_id", leading_is_tight && subleading_is_tight) ) return EL::StatusCode::SUCCESS;

          loose_leading = photonHandler()->passPIDCut(leading_photon, egammaPID::IsEMLoose);
          loose_subleading = photonHandler()->passPIDCut(subleading_photon, egammaPID::IsEMLoose);

          tight_leading = leading_is_tight;
          tight_subleading = subleading_is_tight;

          ph_tight.push_back(leading_is_tight);
          ph_tight.push_back(subleading_is_tight);

          // ----------------------------------------------------------------
          // Apply leading/subleading isolation cut
          // AnalysisBranch - 11=Exotics, 1=Higgs, 0-default, give error
          // ----------------------------------------------------------------

      	  //Isolation correction
          // Reverting xAOD corrections and applying new corrections
          if(correct_isolation){
            topoetcone40_rel17_leading = leading_photon->isolationValue(xAOD::Iso::topoetcone40);
            if(!isoCorr_tool->CorrectLeakage(*leading_photon)) return EL::StatusCode::FAILURE;
            if(!isoCorr_tool->CorrectLeakage(*subleading_photon)) return EL::StatusCode::FAILURE;

            /*
      	    auto electrons = electronHandler()->getCorrectedContainer();
      	    xAOD::Electron* electron = 0;
      	    if(electrons.size() > 0) electron = electrons[0];
      	    //std::cout << "After correction: " << leading_photon->isolationValue(xAOD::Iso::topoetcone40) << std::endl;
      	    if(electron){
      	      //std::cout<<"trying to correct electron\n";
              topoetcone40_rel17_electron_leading = electron->isolationValue(xAOD::Iso::topoetcone40);
              auto electron_tmp = new xAOD::Electron(*electron);
              if(!isoCorr_tool->CorrectLeakage(electron_tmp)) return EL::StatusCode::FAILURE;
              topoetcone40_trouble_electron_leading = electron_tmp->isolationValue(xAOD::Iso::topoetcone40);
              delete electron_tmp;
      	      isoCorr_tool->CorrectLeakage(electron);
              topoetcone40_electron_leading = electron->isolationValue(xAOD::Iso::topoetcone40);
              author_electron_leading = electron->author();
      	      //std::cout  << "Electron correction: " << electron->isolationValue(xAOD::Iso::topoetcone40) << std::endl;
            }
            */

          }

      	  if(AnalysisBranch==11){//Exotics
            bool pass_iso_lead = m_isoTool["FixedCutTightCaloOnly"]->accept(*leading_photon); //photonHandler()->passIsoCut(leading_photon, HG::Iso::FixedCutTightCaloOnly);
            bool pass_iso_sublead = m_isoTool["FixedCutTightCaloOnly"]->accept(*subleading_photon); //photonHandler()->passIsoCut(subleading_photon, HG::Iso::FixedCutTightCaloOnly);
            if( !pass_cut( "pass_ld_subld_isol", pass_iso_lead && pass_iso_sublead )) return EL::StatusCode::SUCCESS;
      	  }

      	  if(AnalysisBranch==1){//Higgs
            bool pass_iso_lead = m_isoTool["FixedCutTight"]->accept(*leading_photon);//photonHandler()->passIsoCut(leading_photon, HG::Iso::FixedCutTight);
            bool pass_iso_sublead = m_isoTool["FixedCutTight"]->accept(*subleading_photon);//photonHandler()->passIsoCut(subleading_photon, HG::Iso::FixedCutTight);
            if( !pass_cut( "pass_ld_subld_isol", pass_iso_lead && pass_iso_sublead ) ) return EL::StatusCode::SUCCESS;
      	  }

          if(AnalysisBranch==11)//Exotics pt cut
          {
            // ----------------------------------------------------------------
            // Apply leading/subleading Et selection
            // ----------------------------------------------------------------

              if( !pass_cut( "pass_ld_subld_Et", (Et_leading >= leading_min_pt) &&
                                                 (Et_subleading >= subleading_min_pt) )
            ) return EL::StatusCode::SUCCESS;
          }else pass_cut( "pass_ld_subld_Et", true);


          if(AnalysisBranch==1)//Higgs relative Et cut
          {
             // ----------------------------------------------------------------
             // Apply leading/subleading Et/m_gg selection
             // ----------------------------------------------------------------

             if( !pass_cut( "pass_ld_subld_rel_cuts", (Et_leading/mass >= leading_rel_cut_pt) &&
                                                      (Et_subleading/mass >= subleading_rel_cut_pt) )
           ) return EL::StatusCode::SUCCESS;
         }else pass_cut( "pass_ld_subld_rel_cuts",true);


          // ----------------------------------------------------------------
          // ********************* EVENT IS ACCEPTED ************************
          // ********************* KINEMATIC VARIABLES **********************
          // ----------------------------------------------------------------

          pass_cut( "accepted", true);
          if(cutflow->event_accepted()) DEBUG("----> Event accepted");

          // *** Computing cinematic variables

          //---->  (costhetastar)
          double l1_plus = Lead_lv.E()+Lead_lv.Pz();
          double l1_minus = Lead_lv.E()-Lead_lv.Pz();
          double l2_plus  = SubLead_lv.E()+SubLead_lv.Pz();
          double l2_minus = SubLead_lv.E()-SubLead_lv.Pz();
          double num = -(l2_plus*l1_minus-l1_plus*l2_minus);
          double denom =  gamgam_lv.M()*sqrt(pow(gamgam_lv.M(),2) + pow(gamgam_lv.Pt(),2));
          costhetastar = num/denom;

          if(is_nominal){
            //*** Filling variables for nominal
            mass_gev = mass*0.001;

	          // **** BLINDING! ****
            // Only for data and gravitons
            if(AnalysisBranch == 11 && !is_mc && mass_gev > 1000) return  EL::StatusCode::SUCCESS;

            weight_sum_selected->SetVal(weight_sum_selected->GetVal() + event_weight);

            pt_subleading = subleading_photon->pt();
            pt_leading = leading_photon->pt();
            phi_subleading = subleading_photon->phi();
            phi_leading = leading_photon->phi();
            eta_subleading = subleading_photon->eta();
            eta_leading = leading_photon->eta();
            topoetcone40_leading = leading_photon->isolationValue(xAOD::Iso::topoetcone40);
            topoetcone40_subleading = subleading_photon->isolationValue(xAOD::Iso::topoetcone40);

            topoetcone20_leading = leading_photon->isolationValue(xAOD::Iso::topoetcone20);
            topoetcone20_subleading = subleading_photon->isolationValue(xAOD::Iso::topoetcone20);

            conv_leading = xAOD::EgammaHelpers::conversionType(leading_photon);
            conv_subleading = xAOD::EgammaHelpers::conversionType(subleading_photon);

            my_tight_leading = photonID_tight_tool_SS->accept(leading_photon);
            isEM_leading = photonID_tight_tool_SS->IsemValue();
            bool passed_all = true;
            //std::cout<<"Now checking isEM mask\n";
            for( unsigned int i = 0; i < isEM_flags_lead.size(); i++ ){
              isEM_flags_lead[i].second = ((isEM_leading >> i) & 1);
              if(i>=10) if(isEM_flags_lead[i].second == 1) passed_all = false;
              //std::cout<<"isEM "<<i<<" leading: "<<((isEM_leading >> i) & 1)<<" "<<isEM_flags_lead[i].second<<"\n";
            }

            photonID_tight_tool_SS->accept(subleading_photon);
            isEM_subleading = photonID_tight_tool_SS->IsemValue();
            for( unsigned int i = 0; i < isEM_flags_sublead.size(); i++ ){
              isEM_flags_sublead[i].second = ((isEM_subleading >> i) & 1);
              //std::cout<<"isEM>>"<<i<<" subleading: "<<((isEM_leading >> i) & 1)<<"\n";
            }

            loose_prime_leading = !(isEM_flags_lead[10].second || isEM_flags_lead[11].second ||\
                                  isEM_flags_lead[12].second || isEM_flags_lead[13].second ||\
                                  isEM_flags_lead[14].second || isEM_flags_lead[15].second ||\
                                  isEM_flags_lead[18].second || isEM_flags_lead[0].second );
            loose_prime_subleading = !(isEM_flags_sublead[10].second || isEM_flags_sublead[11].second ||\
                                     isEM_flags_sublead[12].second || isEM_flags_sublead[13].second ||\
                                     isEM_flags_sublead[14].second || isEM_flags_sublead[15].second ||\
                                     isEM_flags_sublead[18].second || isEM_flags_lead[0].second );

            ph_isEM_flags.push_back(isEM_flags_lead);
            ph_isEM_flags.push_back(isEM_flags_sublead);

            if( my_tight_leading != tight_leading ) std::cout<<"DANGER: handler tool and pID tool not with the same output \n";
            if( my_tight_leading != (isEM_leading == 0) ) std::cout<<"DANGER: pID tool and isEM not with the same output \n";
            if( my_tight_leading != passed_all ) std::cout<<"DANGER: pID tool and single isEM bits not with the same output \n";

            leading_photon->showerShapeValue( Rhad_leading, xAOD::EgammaParameters::ShowerShapeType::Rhad);
            leading_photon->showerShapeValue( e277_leading, xAOD::EgammaParameters::ShowerShapeType::e277);
            leading_photon->showerShapeValue( Reta_leading, xAOD::EgammaParameters::ShowerShapeType::Reta);
            leading_photon->showerShapeValue( Rphi_leading, xAOD::EgammaParameters::ShowerShapeType::Rphi);
            leading_photon->showerShapeValue( weta2_leading, xAOD::EgammaParameters::ShowerShapeType::weta2);
            leading_photon->showerShapeValue( f1_leading, xAOD::EgammaParameters::ShowerShapeType::f1);
            leading_photon->showerShapeValue( DeltaE_leading, xAOD::EgammaParameters::ShowerShapeType::DeltaE);
            leading_photon->showerShapeValue( wtots1_leading, xAOD::EgammaParameters::ShowerShapeType::wtots1);
            leading_photon->showerShapeValue( fracs1_leading, xAOD::EgammaParameters::ShowerShapeType::fracs1);
            leading_photon->showerShapeValue( weta1_leading, xAOD::EgammaParameters::ShowerShapeType::weta1);
            leading_photon->showerShapeValue( Eratio_leading, xAOD::EgammaParameters::ShowerShapeType::Eratio);

            E1_E2_leading = leading_photon->caloCluster()->energyBE(1)/leading_photon->caloCluster()->energyBE(2);
            etaS2_leading = leading_photon->caloCluster()->etaBE(2);

            ptvarcone20_leading = leading_photon->isolationValue(xAOD::Iso::ptvarcone20);
            ptcone20_leading = leading_photon->isolationValue(xAOD::Iso::ptcone20);
            ptvarcone40_leading = leading_photon->isolationValue(xAOD::Iso::ptvarcone40);
            ptcone40_leading = leading_photon->isolationValue(xAOD::Iso::ptcone40);

            ptcone40_subleading = subleading_photon->isolationValue(xAOD::Iso::ptcone40);
            ptcone20_subleading = subleading_photon->isolationValue(xAOD::Iso::ptcone20);

            topoetcone20_Pt_leading = topoetcone20_leading / pt_leading;
            ptvarcone20_Pt_leading = ptvarcone20_leading / pt_leading;

            ph_pt.push_back(pt_leading);
            ph_eta.push_back(eta_leading);
            ph_phi.push_back(phi_leading);
            ph_pt.push_back(pt_subleading);
            ph_eta.push_back(eta_subleading);
            ph_phi.push_back(phi_subleading);

            ph_etas2.push_back(leading_photon->caloCluster()->etaBE(2));
            ph_etas2.push_back(subleading_photon->caloCluster()->etaBE(2));

            ph_cl_phi.push_back(leading_photon->caloCluster()->phi());
            ph_cl_phi.push_back(subleading_photon->caloCluster()->phi());

            topoetcone40.push_back(leading_photon->isolationValue(xAOD::Iso::topoetcone40));
            topoetcone40.push_back(subleading_photon->isolationValue(xAOD::Iso::topoetcone40));

            topoetcone30.push_back(leading_photon->isolationValue(xAOD::Iso::topoetcone30));
            topoetcone30.push_back(subleading_photon->isolationValue(xAOD::Iso::topoetcone30));

            topoetcone20.push_back(leading_photon->isolationValue(xAOD::Iso::topoetcone20));
            topoetcone20.push_back(subleading_photon->isolationValue(xAOD::Iso::topoetcone20));

            ptvarcone40.push_back(leading_photon->isolationValue(xAOD::Iso::ptvarcone40));
            ptvarcone40.push_back(subleading_photon->isolationValue(xAOD::Iso::ptvarcone40));

            ptvarcone30.push_back(leading_photon->isolationValue(xAOD::Iso::ptvarcone30));
            ptvarcone30.push_back(subleading_photon->isolationValue(xAOD::Iso::ptvarcone30));

            ptvarcone20.push_back(leading_photon->isolationValue(xAOD::Iso::ptvarcone20));
            ptvarcone20.push_back(subleading_photon->isolationValue(xAOD::Iso::ptvarcone20));

            ptcone40.push_back(leading_photon->isolationValue(xAOD::Iso::ptcone40));
            ptcone40.push_back(subleading_photon->isolationValue(xAOD::Iso::ptcone40));

            ptcone30.push_back(leading_photon->isolationValue(xAOD::Iso::ptcone30));
            ptcone30.push_back(subleading_photon->isolationValue(xAOD::Iso::ptcone30));

            ptcone20.push_back(leading_photon->isolationValue(xAOD::Iso::ptcone20));
            ptcone20.push_back(subleading_photon->isolationValue(xAOD::Iso::ptcone20));

            // Isolation pre-reccomendation
            FixedCutTightCaloOnly_ld = m_isoTool["FixedCutTightCaloOnly"]->accept(*leading_photon);
            FixedCutTight_ld = m_isoTool["FixedCutTight"]->accept(*leading_photon);
            FixedCutLoose_ld = m_isoTool["FixedCutLoose"]->accept(*leading_photon);

            FixedCutTightCaloOnly_subld = m_isoTool["FixedCutTightCaloOnly"]->accept(*leading_photon);
            FixedCutTight_subld = m_isoTool["FixedCutTight"]->accept(*leading_photon);
            FixedCutLoose_subld = m_isoTool["FixedCutLoose"]->accept(*leading_photon);

            // doing truth matching
            if(is_mc){
              origin_leading = xAOD::TruthHelpers::getParticleTruthOrigin(*leading_photon);
              type_leading = xAOD::TruthHelpers::getParticleTruthType(*leading_photon);

              const xAOD::TruthParticle* truth_leading_photon_from_reco = xAOD::TruthHelpers::getTruthParticle(*leading_photon);
              const xAOD::TruthParticle* truth_subleading_photon_from_reco = xAOD::TruthHelpers::getTruthParticle(*subleading_photon);

              if(truth_leading_photon_from_reco){
                match_leading = true;
                if(truth_leading_photon_from_reco == truth_leading_photon || truth_leading_photon_from_reco == truth_subleading_photon) truth_leading_matched_leading_ph = true;
                else truth_leading_matched_leading_ph = false;
              }
              if(truth_subleading_photon_from_reco){
                match_subleading = true;
                if(truth_subleading_photon_from_reco == truth_subleading_photon || truth_subleading_photon_from_reco == truth_leading_photon) truth_subleading_matched_subleading_ph = true;
                else truth_subleading_matched_subleading_ph = false;
              }

              // truth isolation only in > p2361 datasets
              bool is_p2361 = false;

              double trho[2] = { 0, 0 };

              std::string tesN[2] = { "TruthIsoCentralEventShape", "TruthIsoForwardEventShape"};
              const xAOD::EventShape *tes[2] = {0,0};

              for (int is = 0; is < 2; is++) {
                if ( !event()->retrieve( tes[is], tesN[is] ).isSuccess() ){
                   Warning("execute()", "Failed to retrieve EventShape %s. Skipping.", tesN[is].c_str());
                   is_p2361 = false;
                } else {
                   double rr(0);
                   is_p2361 = true;
                   if (!tes[is]->getDensity(xAOD::EventShape::Density,rr)) Info("","not able to get density");
                   else trho[is] = rr;
                 }
              }

              const xAOD::TruthParticleContainer* egtC = 0;
              if( !event()->retrieve( egtC, "egammaTruthParticles" ).isSuccess() ){
                Warning("execute()", "Failed to retrieve egamma truth particle container. Skipping.");
                is_p2361 = false;
              }

              if(AnalysisBranch==11) is_p2361 = false;

              if(is_p2361){
                truth_rho_central = trho[0];
                truth_rho_forward = trho[1];

                // Write out leading/sub-leading photons regardless of reco match CM
                // For now, leave the rest (including ptcone) how it was, since it's not needed for Cx
                if (truth_leading_photon) {
                  truth_etcone40_leading = truth_leading_photon->auxdataConst<float>("etcone40");
                  truth_etcone20_leading = truth_leading_photon->auxdataConst<float>("etcone20");
                  truth_ptcone40_leading = truth_leading_photon->auxdataConst<float>("ptcone40");
                  truth_ptcone20_leading = truth_leading_photon->auxdataConst<float>("ptcone20");

                  float lead_PUcorr(0);
                  if(fabs(truth_leading_photon->eta()) < 1.5) lead_PUcorr = truth_rho_central;
                  else lead_PUcorr = truth_rho_forward;

                  truth_etcone40_PUcorr_leading = truth_etcone40_leading - lead_PUcorr * 3.1415 * (0.16 - 0.875/128);
                  truth_etcone20_PUcorr_leading = truth_etcone20_leading - lead_PUcorr * 3.1415 * (0.04 - 0.875/128);
                  truth_ptcone40_PUcorr_leading = truth_ptcone40_leading - lead_PUcorr * 3.1415 * (0.16 - 0.875/128);
                  truth_ptcone20_PUcorr_leading = truth_ptcone20_leading - lead_PUcorr * 3.1415 * (0.04 - 0.875/128);

                } else {
                  truth_etcone40_leading = -99;
                  truth_etcone20_leading = -99;
                  truth_ptcone40_leading = -99;
                  truth_ptcone20_leading = -99;
                  truth_etcone40_PUcorr_leading = -99;
                  truth_etcone20_PUcorr_leading = -99;
                  truth_ptcone40_PUcorr_leading = -99;
                  truth_ptcone20_PUcorr_leading = -99;
                }
                if (truth_subleading_photon) {
                  truth_etcone40_subleading = truth_subleading_photon->auxdataConst<float>("etcone40");
                  truth_etcone20_subleading = truth_subleading_photon->auxdataConst<float>("etcone20");
                  truth_ptcone40_subleading = truth_subleading_photon->auxdataConst<float>("ptcone40");
                  truth_ptcone20_subleading = truth_subleading_photon->auxdataConst<float>("ptcone20");

                  float sublead_PUcorr(0);
                  if(fabs(truth_subleading_photon->eta()) < 1.5) sublead_PUcorr = truth_rho_central;
                  else sublead_PUcorr = truth_rho_forward;

                  truth_etcone40_PUcorr_subleading = truth_etcone40_subleading - sublead_PUcorr * 3.1415 * (0.16 - 0.875/128);
                  truth_etcone20_PUcorr_subleading = truth_etcone20_subleading - sublead_PUcorr * 3.1415 * (0.04 - 0.875/128);
                  truth_ptcone40_PUcorr_subleading = truth_ptcone40_subleading - sublead_PUcorr * 3.1415 * (0.16 - 0.875/128);
                  truth_ptcone20_PUcorr_subleading = truth_ptcone20_subleading - sublead_PUcorr * 3.1415 * (0.04 - 0.875/128);
                } else {
                  truth_etcone40_subleading = -99;
                  truth_etcone20_subleading = -99;
                  truth_ptcone40_subleading = -99;
                  truth_ptcone20_subleading = -99;
                  truth_etcone40_PUcorr_subleading = -99;
                  truth_etcone20_PUcorr_subleading = -99;
                  truth_ptcone40_PUcorr_subleading = -99;
                  truth_ptcone20_PUcorr_subleading = -99;
                }
              }
            }

            // ----------------------------------------------------------------
            // Copy selected photons to output xAOD
            // ----------------------------------------------------------------

            xAOD::Photon* leading_photon_store = new xAOD::Photon();
            leading_photon_store->makePrivateStore(*leading_photon);
            leading_photon_store->auxdata< int >( "leading" ) = 1;
            if( ! attach_variables_egamma(leading_photon_store) ) return EL::StatusCode::FAILURE;
            if( ! attach_truth_particle(leading_photon_store) ) return EL::StatusCode::FAILURE;

            xAOD::Photon* subleading_photon_store = new xAOD::Photon();
            subleading_photon_store->makePrivateStore(*subleading_photon);
            subleading_photon_store->auxdata< int >( "leading" ) = 0;
			      if( ! attach_variables_egamma(subleading_photon_store) ) return EL::StatusCode::FAILURE;
			      if( ! attach_truth_particle(subleading_photon_store) ) return EL::StatusCode::FAILURE;

            sel_photons->push_back( leading_photon_store );
            sel_photons->push_back( subleading_photon_store );

            // ----------------------------------------------------------------
            // Use vertex selector tool to get the MVA vertex for the leading/subleading photons
            // ----------------------------------------------------------------

            //const xAOD::Vertex* vertex;
            //m_phVerSel_tool->getVertex(*sel_photons, vertex);
            //std::cout << "\nIs the vertex: " << vertex << "\n";

          }else{
      			//*** Filling systematic variables
      			mass_gev_syst.push_back(mass*0.001);
      			accepted_syst.push_back(cutflow->event_accepted());
            leading_LV_syst.push_back(LV_leading);
            subleading_LV_syst.push_back(LV_subleading);
            diphoton_LV_syst.push_back(LV_diphoton);
      			DEBUG("mass_gev_syst["<<sysList_i<<"] = "<<mass_gev_syst[sysList_i]<<"\n");
      			++sysList_i;
    		  }
      }
      /*
      else{
        pass_cut( "pass_preselection", false);
        return EL::StatusCode::SUCCESS;
      }
      */
      // ***********************************************
      // End of systematics
      // ***********************************************
  }

  //if( ! attach_variables_eventinfo(my_store) ) return EL::StatusCode::FAILURE;

  // TODO: save a SelectedPhotons collection for each systematic -> BEST: use multiple shallow copies

  // Save my_store in the output
  if( ! event()->record( my_store, "my_store" ) ) { return EL::StatusCode::FAILURE; }
  if( ! event()->record( my_storeAux, "my_storeAux." ) ) { return EL::StatusCode::FAILURE; }

  // Record the photons into the output xAOD:
  if( ! event()->record( sel_photons, "SelectedPhotons" ) ) { return EL::StatusCode::FAILURE; }
  if( ! event()->record( sel_photonsAux, "SelectedPhotonsAux." ) ) { return EL::StatusCode::FAILURE; }

  cutflow->finalize_cutflow();

  // Copy EventInfo
  event()->copy("EventInfo");

  // Filling branches of NTUP
  if (sample_NTUP) sample_NTUP->setFilterPassed();
  // Filling output xAOD
  if(!no_xaod) event()->fill();

  return EL::StatusCode::SUCCESS;
}
/*
bool DiphotonAnalysis::pass_cut(std::string cut_name, bool pass){
  bool go_on = true;
  if(!pass){
    DEBUG("	Event not passing "+cut_name);
    if(!save_all_events){
      go_on = false;
      cutflow->fill_cutflow(cut_name, pass);
      cutflow->finalize_cutflow();
    }
  }else{
    DEBUG("	Event passing "+cut_name);
  }
  cutflow->fill_cutflow(cut_name, pass);
  return go_on;
}
*/
bool DiphotonAnalysis::pass_cut(std::string cut_name, bool pass){
  bool go_on = true;
  if(!pass){
    DEBUG("	Event not passing "+cut_name);
    if(save_all_events){
      go_on = true;
    }else if(save_all_preselection &&
            (cut_name == "accepted" || cut_name == "pass_ld_subld_Et" ||
            cut_name == "pass_ld_subld_id" ||  cut_name == "pass_ld_subld_isol" ||
            cut_name == "pass_ld_subld_rel_cuts")
            ){
      go_on = true;
    }else{
      go_on = false;
      cutflow->fill_cutflow(cut_name, pass);
      cutflow->finalize_cutflow();
    }
  }else{
    DEBUG("	Event passing "+cut_name);
    go_on = true;
    cutflow->fill_cutflow(cut_name, pass);
  }
  return go_on;
}

// ****************************************************
//
// Here we add desired variables to photons/electrons
//
// ****************************************************

bool DiphotonAnalysis::attach_variables_egamma(xAOD::IParticle* particle) const{

  xAOD::Egamma* egamma = dynamic_cast<xAOD::Egamma*>(particle);

  // For egammas
  if(egamma){
	  // kinematic variables
	  float etas1 = egamma->caloCluster()->etaBE(1);
	  particle->auxdata< float >( "etas1" ) = etas1;
	  float etas2 = egamma->caloCluster()->etaBE(2);
	  particle->auxdata< float >( "etas2" ) = etas2;

	  // photon exclusive part
	  xAOD::Photon* photon = dynamic_cast<xAOD::Photon*>(egamma);
	  if(photon){
		  const xAOD::Vertex* phVertex = photon->vertex();

		  float convtrk1_nPixHits(0), convtrk2_nPixHits(0);
		  float convtrk1_nSCTHits(0), convtrk2_nSCTHits(0);
		  float convTrk1_nTRTHits(0), convTrk2_nTRTHits(0);
		  float pt_conv(0), pt1_conv(0), pt2_conv(0);

		  // if there is conversion
		  if(phVertex){
			  const Amg::Vector3D pos = phVertex->position();
			  float Rconv = static_cast<float>(hypot (pos.x(), pos.y()));
			  particle->auxdata< float >( "Rconv" ) = Rconv;

			  int convFlag = xAOD::EgammaHelpers::conversionType(photon);
			  particle->auxdata< int >( "convFlag" ) = convFlag;

			  const xAOD::TrackParticle* tp0 = phVertex->trackParticle(0);
			  const xAOD::TrackParticle* tp1 = phVertex->trackParticle(1);

			  bool singleTrackConv = false;
			  if(!tp0 || !tp1) singleTrackConv = true;
			  particle->auxdata< bool >( "singleTrackConv" ) = singleTrackConv;

			  TLorentzVector sum;
			  if(tp0){
				  sum += tp0->p4();

				  uint8_t hits;
				  tp0->summaryValue(hits,xAOD::numberOfPixelHits);
				  convtrk1_nPixHits = hits;
				  tp0->summaryValue(hits,xAOD::numberOfSCTHits);
				  convtrk1_nSCTHits = hits;
				  tp0->summaryValue(hits,xAOD::numberOfTRTHits);
				  convTrk1_nTRTHits = hits;

				  pt1_conv = static_cast<float>(tp0->pt());
			  }

			  if(tp1){
				  sum += tp1->p4();

				  uint8_t hits;
				  tp1->summaryValue(hits,xAOD::numberOfPixelHits);
				  convtrk2_nPixHits = hits;
				  tp1->summaryValue(hits,xAOD::numberOfSCTHits);
				  convtrk2_nSCTHits = hits;
				  tp1->summaryValue(hits,xAOD::numberOfTRTHits);
				  convTrk2_nTRTHits = hits;

				  pt2_conv = static_cast<float>(tp1->pt());
			  }
			  pt_conv = sum.Perp();
		  }

		  particle->auxdata< float >( "pt_conv" ) = pt_conv;
		  particle->auxdata< float >( "pt1_conv" ) = pt1_conv;
		  particle->auxdata< float >( "pt2_conv" ) = pt2_conv;

		  particle->auxdata< int >( "convtrk1_nPixHits" ) = convtrk1_nPixHits;
		  particle->auxdata< int >( "convtrk1_nSCTHits" ) = convtrk1_nSCTHits;
		  particle->auxdata< int >( "convTrk1_nTRTHits" ) = convTrk1_nTRTHits;
		  particle->auxdata< int >( "convtrk2_nPixHits" ) = convtrk2_nPixHits;
		  particle->auxdata< int >( "convtrk2_nSCTHits" ) = convtrk2_nSCTHits;
      particle->auxdata< int >( "convTrk2_nTRTHits" ) = convTrk2_nTRTHits;
	  }

	  // electron exclusive part
	  xAOD::Electron* electron = dynamic_cast<xAOD::Electron*>(egamma);
	  if(electron){
		// ...
	  }
  }

  return true;
}

// ****************************************************
//
// Here we add desired variables to eventInfo
//
// ****************************************************

bool DiphotonAnalysis::attach_variables_eventinfo(xAOD::EventInfo* event_info) {
  // Save the invariant mass
  event_info->auxdata< float > ("mass_gev") = mass_gev;
  event_info->auxdata< float > ("mass") = mass;
  event_info->auxdata< float > ("costhetastar") = costhetastar;

  // Add pass flags to the event
  for( unsigned int i = 0; i < pass_flags.size(); i++ ){
    event_info->auxdata< bool >(pass_flags[i].first.c_str()) = pass_flags[i].second;
  }

  if(do_systematics){
    int sysList_i = 0;
    for(auto sys : sysList){
      std::string br_name_mass = "mass_" + (sys).name();
      event_info->auxdata< float > (br_name_mass) = mass_gev_syst[sysList_i];
      std::string br_name_pass = "pass_sel_" + (sys).name();
      event_info->auxdata< float > (br_name_pass) = accepted_syst[sysList_i];
      sysList_i++;
    }
  }

  // TODO: only for rel20
  const xAOD::EventShape* eventShape_central_xaod = 0;
  if( ! event()->retrieve( eventShape_central_xaod, "TopoClusterIsoCentralEventShape").isSuccess() ){
    Error("execute()", "Failed to retrieve event TopoClusterIsoCentralEventShape collection. Exiting." );
    return false;
  }
  double density_central_xaod = 0;
  eventShape_central_xaod->getDensity(xAOD::EventShape::Density, density_central_xaod);
  event_info->auxdata< float > ("ED_central") = density_central_xaod;

  const xAOD::EventShape* eventShape_forward_xaod = 0;
  if( ! event()->retrieve( eventShape_forward_xaod, "TopoClusterIsoForwardEventShape").isSuccess() ){
    Error("execute()", "Failed to retrieve event TopoClusterIsoForwardEventShape collection. Exiting." );
    return false;
  }
  double density_forward_xaod = 0;
  eventShape_forward_xaod->getDensity(xAOD::EventShape::Density, density_forward_xaod);
  event_info->auxdata< float > ("ED_forward") = density_forward_xaod;

  return true;
}

// ****************************************************
//
// Here we attach truth informations
//
// ****************************************************

bool DiphotonAnalysis::attach_truth_particle(xAOD::IParticle* particle) const{
  if(is_mc){
	  const xAOD::TruthParticle* true_particle = xAOD::TruthHelpers::getTruthParticle(*particle);
	  bool truth_matched(false);
	  float truth_pt(0), truth_phi(0), truth_eta(0), truth_E(0);
	  int truth_pdgId(0);

	  if(true_particle){
		truth_matched = true;
		truth_pt = true_particle->pt();
		truth_phi = true_particle->phi();
		truth_eta = true_particle->eta();
		truth_E = true_particle->e();
		truth_pdgId = true_particle->pdgId();
	  }

	  particle->auxdata< bool >( "truth_matched" ) = truth_matched;
	  particle->auxdata< float >( "truth_pt" ) = truth_pt;
	  particle->auxdata< float >( "truth_phi" ) = truth_phi;
	  particle->auxdata< float >( "truth_eta" ) = truth_eta;
	  particle->auxdata< float >( "truth_E" ) = truth_E;
	  particle->auxdata< int >( "truth_pdgId" ) = truth_pdgId;
  }
  return true;
}

// ****************************************************
//
// Here we compute truth isolation
//
// ****************************************************

std::vector<double> DiphotonAnalysis::get_truth_iso(const xAOD::TruthParticle *ptcl, HG::TruthPtcls &stblPtcls) {
  std::vector<double> isoETs;
  HG::TruthPtcls decay = HG::getStableDecayProducts(ptcl);
  TLorentzVector iso20(0,0,0,0), iso40(0,0,0,0), iso30(0,0,0,0);
  // Calculate the truth isolation!!
  for (auto p:stblPtcls) {
    if (p->p4().Pt()<1e-3) continue;
    double dr=ptcl->p4().DeltaR(p->p4());
    if (dr>0.4) continue;
    // now, check if the particle p in question is part of the decay
    bool partOfDecay=false, print=false;
    for (auto d:decay) if (p==d) partOfDecay=true;
    if (partOfDecay) {
      if (ptcl->p4().DeltaR(p->p4())>0.1) print=true;
      else continue;
    }
    if (0&&print) {
      HG::printTruthPtcl(ptcl,"Origin");
      for (auto d:decay)
        HG::printTruthPtcl(d,Form("Decay ptcl DR(orign)=%.2f",
                                  d->p4().DeltaR(ptcl->p4())));
    }
    iso40 += p->p4();
    if (dr<0.3) iso30+=p->p4();
    if (dr<0.2) iso20+=p->p4();
  }
  isoETs.push_back(iso40.Et()*HG::invGeV);
  isoETs.push_back(iso30.Et()*HG::invGeV);
  isoETs.push_back(iso20.Et()*HG::invGeV);
  return isoETs;
}

// ****************************************************
//
// Here we do a truth analysis on photons
//
// ****************************************************

EL::StatusCode DiphotonAnalysis::run_truth_analysis(){
/*
	const xAOD::TruthParticleContainer* truth_particles = 0;
	if( !event()->retrieve( truth_particles, "TruthParticles" ).isSuccess() ){
		Error("execute()", "Failed to retrieve truth particle container. Exiting." );
		return EL::StatusCode::FAILURE;
	}
*/
  auto truth_particles = truthHandler()->getPhotons();
	truth_leading_photon = 0;
  truth_subleading_photon = 0;

	float truth_pt_leading = 0.;
  float truth_pt_subleading = 0.;

	for(auto truth_particle : truth_particles){
    if(truth_particle->pdgId() != 22) continue;
		if(truth_particle->pt() > truth_pt_leading){
      truth_subleading_photon = truth_leading_photon;
      truth_pt_subleading = truth_pt_leading;
      truth_pt_leading = truth_particle->pt();
			truth_leading_photon = truth_particle;
		}else if(truth_particle->pt() > truth_pt_subleading){
      truth_subleading_photon = truth_particle;
      truth_pt_subleading = truth_particle->pt();
    }
	}

  if(truth_leading_photon && truth_subleading_photon){
    two_truth_photons = true;
    //Info("run_truth_analysis()", "Truth_leading pt = %f", truth_pt_leading);
    //Info("run_truth_analysis()", "Truth_subleading pt = %f", truth_pt_subleading);
    truth_leading_LV = truth_leading_photon->p4();
    truth_subleading_LV = truth_subleading_photon->p4();
    truth_diphoton_LV = truth_leading_LV + truth_subleading_LV;
    float eta_l = fabs(truth_leading_photon->eta());
    float eta_sl = fabs(truth_subleading_photon->eta());
    if( ( eta_l < 1.37 || (eta_l < 2.37 && eta_l > 1.52) ) &&
        ( eta_sl < 1.37 || (eta_sl < 2.37 && eta_sl > 1.52) )
    ){
      pass_eta_truth_analysis = true;
      //Info("run_truth_analysis()", "Truth_leading eta = %f", eta_l);
      //Info("run_truth_analysis()", "Truth_subleading eta = %f", eta_sl);
    }else pass_eta_truth_analysis = false;
  }else{
    two_truth_photons = false;
    pass_eta_truth_analysis = false;
    truth_leading_LV.SetPtEtaPhiE(0.,0.,0.,0.);
    truth_subleading_LV.SetPtEtaPhiE(0.,0.,0.,0.);
  }

	return EL::StatusCode::SUCCESS;
}

EL::StatusCode DiphotonAnalysis::setupJob (EL::Job& job)
{
  HgammaAnalysis::setupJob(job); // keep this line!

  // let's initialize the algorithm to use the xAODRootAccess package
  xAOD::Init( "MyxAODAnalysis" ).ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode DiphotonAnalysis::fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  HgammaAnalysis::fileExecute(); // keep this line!

  xAOD::TEvent* tmp_event = wk()->xaodEvent();

  const char *APP_NAME = "fileExecute()";
  // Check file's metadata:
  TTree *MetaData = dynamic_cast<TTree*>(wk()->inputFile()->Get("MetaData"));
  if (!MetaData) {
    Error("fileExecute()", "MetaData not found! Exiting.");
    return EL::StatusCode::FAILURE;
  }
  MetaData->LoadTree(0);
  bool m_isDerivation = !MetaData->GetBranch("StreamAOD");

  if(m_isDerivation ){

    // check for corruption
    const xAOD::CutBookkeeperContainer* incompleteCBC = nullptr;
    if(!tmp_event->retrieveMetaInput(incompleteCBC, "IncompleteCutBookkeepers").isSuccess()){
      Error("initializeEvent()","Failed to retrieve IncompleteCutBookkeepers from MetaData! Exiting.");
      return EL::StatusCode::FAILURE;
    }
    if ( incompleteCBC->size() != 0 ) {
      Warning("initializeEvent()","Found incomplete Bookkeepers! Check file for corruption.");
      //return EL::StatusCode::FAILURE;
    }


    //Read the CutBookkeeper container
    const xAOD::CutBookkeeperContainer* completeCBC = 0;
    if (!tmp_event->retrieveMetaInput(completeCBC, "CutBookkeepers").isSuccess()) {
      Error( APP_NAME, "Failed to retrieve CutBookkeepers from MetaData! Exiting.");
    }

    // First, let's find the smallest cycle number,
    // i.e., the original first processing step/cycle
    int minCycle = 10000;
    for ( auto cbk : *completeCBC ) {
      if ( ! cbk->name().empty()  && minCycle > cbk->cycle() ) { minCycle = cbk->cycle(); }
    }
    // Now, let's actually find the right one that contains all the needed info...
    const xAOD::CutBookkeeper* allEventsCBK = 0;
    //const xAOD::CutBookkeeper* DxAODEventsCBK=0;
    std::string derivationName = "EXOT10Kernel"; //need to replace by appropriate name
    for ( auto cbk :  *completeCBC ) {
      if ( minCycle == cbk->cycle() && cbk->name() == "AllExecutedEvents" ) {
        allEventsCBK = cbk;
            }

      //if ( cbk->name() == derivationName){
      //  DxAODEventsCBK = cbk;
      //}
    }

    int m_initialEvents = allEventsCBK->nAcceptedEvents();
    float m_initialSumOfWeights = allEventsCBK->sumOfEventWeights();
    float m_initialSumOfWeightsSquared = allEventsCBK->sumOfEventWeightsSquared();
    Info( APP_NAME, "CutBookkeepers Accepted %d SumWei %f sumWei2 %f ",m_initialEvents, m_initialSumOfWeights, m_initialSumOfWeightsSquared);

    //int m_initialEventsDxAOD = DxAODEventsCBK->nAcceptedEvents();
    //float m_initialSumOfWeightsDxAOD = DxAODEventsCBK->sumOfEventWeights();
    //float m_initialSumOfWeightsSquaredDxAOD = DxAODEventsCBK->sumOfEventWeightsSquared();
    //std::cout<<m_initialEvents<<" "<<m_initialSumOfWeights<<" "<<m_initialSumOfWeightsSquared<<std::endl;
    //std::cout<<m_initialEventsDxAOD<<" "<<m_initialSumOfWeightsDxAOD<<" "<<m_initialSumOfWeightsSquaredDxAOD<<std::endl;

    Info( "fileExecute()", "Total initial sum of weight before sum is: %f", weight_sum_before->GetVal());
    weight_sum_before_f = m_initialSumOfWeights;
    weight_sum_before->SetVal(weight_sum_before->GetVal() + m_initialSumOfWeights);
    Info( "fileExecute()", "  After sum: %f", weight_sum_before->GetVal());
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode DiphotonAnalysis::histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  HgammaAnalysis::histInitialize(); // keep this line!

  std::cout<<"\n\nIn my histinitialize\n\n";

  nph_hist = new TH1F("nph_hist", "nph_hist", 15, 0, 15);
  pileup_hist = new TH1F("pileup_hist", "pileup_hist", 100, 0, 5);

  pass_flags   = { std::make_pair("all", false),
                   std::make_pair("pass_dalitz", false),
                   std::make_pair("pass_trigger", false),
                   std::make_pair("pass_grl", false),
                   std::make_pair("pass_detector_DQ", false),
                   std::make_pair("pass_PV", false),
                   std::make_pair("pass_preselection", false),
                   std::make_pair("pass_ld_subld_id", false),
                   std::make_pair("pass_ld_subld_isol", false),
                   std::make_pair("pass_ld_subld_Et", false),
                   std::make_pair("pass_ld_subld_rel_cuts", false),
                   std::make_pair("accepted", false)
                 };

isEM_flags_lead  = { std::make_pair("isEM_0_ld", 0),
                     std::make_pair("isEM_1_ld", 0),
                     std::make_pair("isEM_2_ld", 0),
                     std::make_pair("isEM_3_ld", 0),
                     std::make_pair("isEM_4_ld", 0),
                     std::make_pair("isEM_5_ld", 0),
                     std::make_pair("isEM_6_ld", 0),
                     std::make_pair("isEM_7_ld", 0),
                     std::make_pair("isEM_8_ld", 0),
                     std::make_pair("isEM_9_ld", 0),
                     std::make_pair("isEM_10_ld", 0),
                     std::make_pair("isEM_11_ld", 0),
                     std::make_pair("isEM_12_ld", 0),
                     std::make_pair("isEM_13_ld", 0),
                     std::make_pair("isEM_14_ld", 0),
                     std::make_pair("isEM_15_ld", 0),
                     std::make_pair("isEM_16_ld", 0),
                     std::make_pair("isEM_17_ld", 0),
                     std::make_pair("isEM_18_ld", 0),
                     std::make_pair("isEM_19_ld", 0),
                     std::make_pair("isEM_20_ld", 0),
                     std::make_pair("isEM_21_ld", 0)
                };
  isEM_flags_sublead  = { std::make_pair("isEM_0_subld", 0),
                       std::make_pair("isEM_1_subld", 0),
                       std::make_pair("isEM_2_subld", 0),
                       std::make_pair("isEM_3_subld", 0),
                       std::make_pair("isEM_4_subld", 0),
                       std::make_pair("isEM_5_subld", 0),
                       std::make_pair("isEM_6_subld", 0),
                       std::make_pair("isEM_7_subld", 0),
                       std::make_pair("isEM_8_subld", 0),
                       std::make_pair("isEM_9_subld", 0),
                       std::make_pair("isEM_10_subld", 0),
                       std::make_pair("isEM_11_subld", 0),
                       std::make_pair("isEM_12_subld", 0),
                       std::make_pair("isEM_13_subld", 0),
                       std::make_pair("isEM_14_subld", 0),
                       std::make_pair("isEM_15_subld", 0),
                       std::make_pair("isEM_16_subld", 0),
                       std::make_pair("isEM_17_subld", 0),
                       std::make_pair("isEM_18_subld", 0),
                       std::make_pair("isEM_19_subld", 0),
                       std::make_pair("isEM_20_subld", 0),
                       std::make_pair("isEM_21_subld", 0)
                  };
  cutflow = new CutFlowHisto("cutflow", "cutflow");
  cutflow->add_cuts(pass_flags);

  wk()->addOutput(cutflow);
  wk()->addOutput(nph_hist);
  wk()->addOutput(pileup_hist);

  weight_sum_before = new TParameter<float>("weight_sum_before", 0., '+');
  weight_sum_before->SetVal(0.);
  weight_sum_selected = new TParameter<float>("weight_sum_selected", 0., '+');
  weight_sum_selected->SetVal(0.);

  wk()->addOutput(weight_sum_before);
  wk()->addOutput(weight_sum_selected);
  Info( "histInitialize()", "Total initial sum of weight is: %f", weight_sum_before->GetVal());
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode DiphotonAnalysis::changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  HgammaAnalysis::changeInput(firstFile); // keep this line!
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode DiphotonAnalysis::initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // do not get called if no events are processed.  So any objects
  // you create here will not be available in the output if you have no
  // input events.

  // TEnv uses value from first file it's specified in.
  // Fill unspecified values from default config, specified here.
  HgammaAnalysis::initialize(); // keep this line!

  my_XsecDB = new SUSY::CrossSectionDB(gSystem->ExpandPathName("$ROOTCOREBIN/data/LowHighMyy/cross_sections_13TeV.txt"));

  bool isMC_tmp = eventInfo()->eventType(xAOD::EventInfo::IS_SIMULATION);

  // initialize isolation correction tool
  isoCorr_tool = new CP::IsolationCorrectionTool("IsoCorr_tool");
  isoCorr_tool->setProperty( "CorrFile", "IsolationCorrections/isolation_ptcorrections_rel20_2.root");
  isoCorr_tool->setProperty( "Apply_datadriven", false);
  isoCorr_tool->setProperty( "Correct_etcone", false);
  isoCorr_tool->setProperty( "Trouble_categories", true);
  isoCorr_tool->setProperty( "IsMC", isMC_tmp);
  isoCorr_tool->setProperty( "ToolVer", "REL20_2");
  if (! isoCorr_tool->initialize().isSuccess() ){
    Error("initialize()", "Failed to properly initialize the isoCorr Tool. Exiting." );
    return EL::StatusCode::FAILURE;
  }

  // create the tight selector
  photonID_tight_tool_SS = new AsgPhotonIsEMSelector ( "PhotonTightIsEMSelector_SS" );
  // decide which kind of selection (loose/medium/tight) you want to use
  photonID_tight_tool_SS->setProperty("isEMMask",egammaPID::PhotonTight);
  // set the file that contains the cuts on the shower shapes (stored in http://atlas.web.cern.ch/Atlas/GROUPS/DATABASE/GroupData/)
  photonID_tight_tool_SS->setProperty("ConfigFile","ElectronPhotonSelectorTools/offline/mc15_20150712/PhotonIsEMTightSelectorCutDefs.conf");
  // Use corrected shower shapes
  if(isMC_tmp) photonID_tight_tool_SS->setProperty("UseFudgedShowers",true);
  // initialise the tool
  if (!photonID_tight_tool_SS->initialize().isSuccess()) {
    Fatal("MyFunction", "Failed to initialize PhotonTightIsEMSelector_SS");
  }

  // Photon vertex selector tool
  m_phVerSel_tool = new CP::PhotonVertexSelectionTool("PhotonVertexSelector");
  m_phVerSel_tool->initialize().ignore();

  // Isolation selection, different tunes
  CP::IsolationSelectionTool* isoFixedCutTightCaloOnly = new CP::IsolationSelectionTool("FixedCutTightCaloOnlyWPIso");
  isoFixedCutTightCaloOnly->setProperty("photonWP","FixedCutTightCaloOnly");
  if (!isoFixedCutTightCaloOnly->initialize().isSuccess()) {
    Fatal("MyFunction", "Failed to initialize IsolationSelectionTool");
  }
  m_isoTool["FixedCutTightCaloOnly"] = isoFixedCutTightCaloOnly;

  CP::IsolationSelectionTool* isoFixedCutTight = new CP::IsolationSelectionTool("FixedCutTight");
  isoFixedCutTight->setProperty("photonWP","FixedCutTight");
  if (!isoFixedCutTight->initialize().isSuccess()) {
    Fatal("MyFunction", "Failed to initialize IsolationSelectionTool");
  }
  m_isoTool["FixedCutTight"] = isoFixedCutTight;

  CP::IsolationSelectionTool* isoFixedCutLoose = new CP::IsolationSelectionTool("FixedCutLoose");
  isoFixedCutLoose->setProperty("photonWP","FixedCutLoose");
  if (!isoFixedCutLoose->initialize().isSuccess()) {
    Fatal("MyFunction", "Failed to initialize IsolationSelectionTool");
  }
  m_isoTool["FixedCutLoose"] = isoFixedCutLoose;

  // loop over systematics registry:
  const CP::SystematicRegistry& registry = CP::SystematicRegistry::getInstance();
  const CP::SystematicSet& recommendedSystematics = registry.recommendedSystematics(); // get list of recommended systematics

  // this is the nominal set, no systematic
  sysList.push_back(CP::SystematicSet());

  // loop over recommended systematics
  for(auto sys : recommendedSystematics){ //CP::SystematicSet::const_iterator sysItr = recommendedSystematics.begin(); sysItr != recommendedSystematics.end(); ++sysItr){
    std::string br_name = sys.name();
    if( (br_name.find("EG") == std::string::npos) && (br_name.find("PH") == std::string::npos) ) continue;
    sysList.push_back( CP::SystematicSet() );
    sysList.back().insert( sys );
    //std::cout<<"\nSystematic: "<<(*sysItr).name()<<"\n";
  } // end for loop over recommended systematics

  if (!outputStreamName.empty()) {
      sample_NTUP = EL::getNTupleSvc(wk(), outputStreamName);
  } else {
      sample_NTUP = 0;
  }

  if (sample_NTUP) {
      // Pass flags
      for( unsigned int i = 0; i < pass_flags.size(); i++ ){
         sample_NTUP->tree()->Branch(pass_flags[i].first.c_str(), &pass_flags[i].second);
      }
      for( unsigned int i = 0; i < isEM_flags_lead.size(); i++ ){
         sample_NTUP->tree()->Branch(isEM_flags_lead[i].first.c_str(), &isEM_flags_lead[i].second);
         sample_NTUP->tree()->Branch(isEM_flags_sublead[i].first.c_str(), &isEM_flags_sublead[i].second);
      }
      sample_NTUP->tree()->Branch("ph_isEM_flags", &ph_isEM_flags);

      // Event variables
      sample_NTUP->tree()->Branch("averageIntPerXing", &averageIntPerXing);
      sample_NTUP->tree()->Branch("ED_central", &ED_central);
      sample_NTUP->tree()->Branch("ED_forward", &ED_forward);

      sample_NTUP->tree()->Branch("run_number", &run_number);
      sample_NTUP->tree()->Branch("event_number", &event_number);

      sample_NTUP->tree()->Branch("nPV", &npvs);
      sample_NTUP->tree()->Branch("xs", &xs);
      sample_NTUP->tree()->Branch("xs_ami", &xs_ami);
      sample_NTUP->tree()->Branch("filter_eff", &filter_eff);
      sample_NTUP->tree()->Branch("filter_eff_ami", &filter_eff_ami);

      // Weights
      sample_NTUP->tree()->Branch("pileup_weight", &pileup_weight);
      sample_NTUP->tree()->Branch("MC_weight", &MC_weight);
      sample_NTUP->tree()->Branch("xs_weight", &xs_weight);
      sample_NTUP->tree()->Branch("event_weight", &event_weight);
      sample_NTUP->tree()->Branch("prel_weight", &prel_weight);
      sample_NTUP->tree()->Branch("total_events", &total_events);
      sample_NTUP->tree()->Branch("non_derived_total_events", &non_derived_total_events);
      sample_NTUP->tree()->Branch("weight_sum_before", &weight_sum_before_f);

      // Cinematic variables
      sample_NTUP->tree()->Branch("mass_gev", &mass_gev)->SetTitle("Invariant mass in GeV");
      sample_NTUP->tree()->Branch("mass", &mass)->SetTitle("Invariant mass in MeV");
      sample_NTUP->tree()->Branch("costhetastar", &costhetastar);
      sample_NTUP->tree()->Branch("deltaphi", &deltaphi);

      sample_NTUP->tree()->Branch("LV_leading", &LV_leading);
      sample_NTUP->tree()->Branch("LV_subleading", &LV_subleading);
      sample_NTUP->tree()->Branch("LV_diphoton", &LV_diphoton);

      sample_NTUP->tree()->Branch("pt_subleading", &pt_subleading);
      sample_NTUP->tree()->Branch("pt_leading", &pt_leading);
      sample_NTUP->tree()->Branch("phi_subleading", &phi_subleading);
      sample_NTUP->tree()->Branch("phi_leading", &phi_leading);
      sample_NTUP->tree()->Branch("eta_subleading", &eta_subleading);
      sample_NTUP->tree()->Branch("eta_leading", &eta_leading);

      sample_NTUP->tree()->Branch("topoetcone40_leading", &topoetcone40_leading);
      sample_NTUP->tree()->Branch("topoetcone40_rel17_leading", &topoetcone40_rel17_leading);

      sample_NTUP->tree()->Branch("topoetcone40_electron_leading", &topoetcone40_electron_leading);
      sample_NTUP->tree()->Branch("topoetcone40_trouble_electron_leading", &topoetcone40_trouble_electron_leading);
      sample_NTUP->tree()->Branch("topoetcone40_rel17_electron_leading", &topoetcone40_rel17_electron_leading);
      sample_NTUP->tree()->Branch("author_electron_leading", &author_electron_leading);

      sample_NTUP->tree()->Branch("topoetcone40_subleading", &topoetcone40_subleading);
      sample_NTUP->tree()->Branch("topoetcone20_leading", &topoetcone20_leading);
      sample_NTUP->tree()->Branch("topoetcone20_subleading", &topoetcone20_subleading);

      sample_NTUP->tree()->Branch("truth_etcone40_leading", &truth_etcone40_leading);
      sample_NTUP->tree()->Branch("truth_etcone40_subleading", &truth_etcone40_subleading);
      sample_NTUP->tree()->Branch("truth_etcone20_leading", &truth_etcone20_leading);
      sample_NTUP->tree()->Branch("truth_etcone20_subleading", &truth_etcone20_subleading);

      sample_NTUP->tree()->Branch("truth_local_etcone40_leading", &truth_local_etcone40_leading);
      sample_NTUP->tree()->Branch("truth_local_etcone40_subleading", &truth_local_etcone40_subleading);
      sample_NTUP->tree()->Branch("truth_local_etcone20_leading", &truth_local_etcone20_leading);
      sample_NTUP->tree()->Branch("truth_local_etcone20_subleading", &truth_local_etcone20_subleading);

      sample_NTUP->tree()->Branch("truth_etcone40_PUcorr_leading", &truth_etcone40_PUcorr_leading);
      sample_NTUP->tree()->Branch("truth_etcone40_PUcorr_subleading", &truth_etcone40_PUcorr_subleading);
      sample_NTUP->tree()->Branch("truth_etcone20_PUcorr_leading", &truth_etcone20_PUcorr_leading);
      sample_NTUP->tree()->Branch("truth_etcone20_PUcorr_subleading", &truth_etcone20_PUcorr_subleading);

      sample_NTUP->tree()->Branch("truth_ptcone40_leading", &truth_ptcone40_leading);
      sample_NTUP->tree()->Branch("truth_ptcone40_subleading", &truth_ptcone40_subleading);
      sample_NTUP->tree()->Branch("truth_ptcone20_leading", &truth_ptcone20_leading);
      sample_NTUP->tree()->Branch("truth_ptcone20_subleading", &truth_ptcone20_subleading);

      sample_NTUP->tree()->Branch("truth_ptcone40_PUcorr_leading", &truth_ptcone40_PUcorr_leading);
      sample_NTUP->tree()->Branch("truth_ptcone40_PUcorr_subleading", &truth_ptcone40_PUcorr_subleading);
      sample_NTUP->tree()->Branch("truth_ptcone20_PUcorr_leading", &truth_ptcone20_PUcorr_leading);
      sample_NTUP->tree()->Branch("truth_ptcone20_PUcorr_subleading", &truth_ptcone20_PUcorr_subleading);

      sample_NTUP->tree()->Branch("truth_rho_central", &truth_rho_central);
      sample_NTUP->tree()->Branch("truth_rho_forward", &truth_rho_forward);

      sample_NTUP->tree()->Branch("tight_leading", &tight_leading);
      sample_NTUP->tree()->Branch("my_tight_leading", &my_tight_leading);
      sample_NTUP->tree()->Branch("loose_leading", &loose_leading);
      sample_NTUP->tree()->Branch("loose_prime_leading", &loose_prime_leading);
      sample_NTUP->tree()->Branch("tight_subleading", &tight_subleading);
      sample_NTUP->tree()->Branch("loose_subleading", &loose_subleading);
      sample_NTUP->tree()->Branch("loose_prime_subleading", &loose_prime_subleading);

      sample_NTUP->tree()->Branch("match_leading", &match_leading);
      sample_NTUP->tree()->Branch("match_subleading", &match_subleading);

      sample_NTUP->tree()->Branch("bg_truth_match_leading", &bg_truth_match_leading);
      sample_NTUP->tree()->Branch("bg_truth_match_origin_leading", &bg_truth_match_origin_leading);

      sample_NTUP->tree()->Branch("origin_leading", &origin_leading);
      sample_NTUP->tree()->Branch("type_leading", &type_leading);

      sample_NTUP->tree()->Branch("conv_leading", &conv_leading);
      sample_NTUP->tree()->Branch("conv_subleading", &conv_subleading);

      sample_NTUP->tree()->Branch("isEM_leading", &isEM_leading);
      sample_NTUP->tree()->Branch("isEM_subleading", &isEM_subleading);

      sample_NTUP->tree()->Branch("Rhad_leading", &Rhad_leading);
      sample_NTUP->tree()->Branch("e277_leading", &e277_leading);
      sample_NTUP->tree()->Branch("Reta_leading", &Reta_leading);
      sample_NTUP->tree()->Branch("Rphi_leading", &Rphi_leading);
      sample_NTUP->tree()->Branch("weta2_leading", &weta2_leading);
      sample_NTUP->tree()->Branch("f1_leading", &f1_leading);
      sample_NTUP->tree()->Branch("DeltaE_leading", &DeltaE_leading);
      sample_NTUP->tree()->Branch("wtots1_leading", &wtots1_leading);
      sample_NTUP->tree()->Branch("fracs1_leading", &fracs1_leading);
      sample_NTUP->tree()->Branch("weta1_leading", &weta1_leading);
      sample_NTUP->tree()->Branch("Eratio_leading", &Eratio_leading);

      sample_NTUP->tree()->Branch("E1_E2_leading", &E1_E2_leading);
      sample_NTUP->tree()->Branch("etaS2_leading", &etaS2_leading);

      sample_NTUP->tree()->Branch("ptvarcone20_leading", &ptvarcone20_leading);
      sample_NTUP->tree()->Branch("ptcone20_leading", &ptcone20_leading);
      sample_NTUP->tree()->Branch("ptvarcone40_leading", &ptvarcone40_leading);
      sample_NTUP->tree()->Branch("ptcone40_leading", &ptcone40_leading);
      sample_NTUP->tree()->Branch("topoetcone20_Pt_leading", &topoetcone20_Pt_leading);
      sample_NTUP->tree()->Branch("ptvarcone20_Pt_leading", &ptvarcone20_Pt_leading);

      sample_NTUP->tree()->Branch("ptcone20_subleading", &ptcone20_subleading);
      sample_NTUP->tree()->Branch("ptcone40_subleading", &ptcone40_subleading);


      // Isolation variables
      sample_NTUP->tree()->Branch("topoetcone40", &topoetcone40);
      sample_NTUP->tree()->Branch("topoetcone30", &topoetcone30);
      sample_NTUP->tree()->Branch("topoetcone20", &topoetcone20);

      sample_NTUP->tree()->Branch("ptcone40", &ptcone40);
      sample_NTUP->tree()->Branch("ptcone30", &ptcone30);
      sample_NTUP->tree()->Branch("ptcone20", &ptcone20);

      sample_NTUP->tree()->Branch("ptvarcone40", &ptvarcone40);
      sample_NTUP->tree()->Branch("ptvarcone30", &ptvarcone30);
      sample_NTUP->tree()->Branch("ptvarcone20", &ptvarcone20);

      sample_NTUP->tree()->Branch("FixedCutTightCaloOnly_ld", &FixedCutTightCaloOnly_ld);
      sample_NTUP->tree()->Branch("FixedCutTight_ld", &FixedCutTight_ld);
      sample_NTUP->tree()->Branch("FixedCutLoose_ld", &FixedCutLoose_ld);

      sample_NTUP->tree()->Branch("FixedCutTightCaloOnly_subld", &FixedCutTightCaloOnly_subld);
      sample_NTUP->tree()->Branch("FixedCutTight_subld", &FixedCutTight_subld);
      sample_NTUP->tree()->Branch("FixedCutLoose_subld", &FixedCutLoose_subld);

      sample_NTUP->tree()->Branch("pass_truth_match", &pass_truth_match);
      sample_NTUP->tree()->Branch("parent_pdgID_ld", &parent_pdgID_ld);
      sample_NTUP->tree()->Branch("parent_pdgID_subld", &parent_pdgID_subld);

      sample_NTUP->tree()->Branch("ph_pt", &ph_pt);
      sample_NTUP->tree()->Branch("ph_eta", &ph_eta);
      sample_NTUP->tree()->Branch("ph_etas2", &ph_etas2);
      sample_NTUP->tree()->Branch("ph_phi", &ph_phi);
      sample_NTUP->tree()->Branch("ph_cl_phi", &ph_cl_phi);
      sample_NTUP->tree()->Branch("ph_parent_pdgID", &ph_parent_pdgID);
      sample_NTUP->tree()->Branch("ph_tight", &ph_tight);
      sample_NTUP->tree()->Branch("ph_matched", &ph_matched);


      // truth analysis variables
      sample_NTUP->tree()->Branch("two_truth_photons", &two_truth_photons);
      sample_NTUP->tree()->Branch("pass_eta_truth_analysis", &pass_eta_truth_analysis);
      sample_NTUP->tree()->Branch("truth_leading_LV", &truth_leading_LV);
      sample_NTUP->tree()->Branch("truth_subleading_LV", &truth_subleading_LV);
      sample_NTUP->tree()->Branch("truth_diphoton_LV", &truth_diphoton_LV);
      sample_NTUP->tree()->Branch("truth_leading_matched_leading_ph", &truth_leading_matched_leading_ph);
      sample_NTUP->tree()->Branch("truth_subleading_matched_subleading_ph", &truth_subleading_matched_subleading_ph);

      // Systematics variables
      if(do_systematics){
          make_syst_branches("mass", mass_gev_syst);
          make_syst_branches("accepted", accepted_syst);
          make_syst_branches("leading_LV", leading_LV_syst);
          make_syst_branches("subleading_LV", subleading_LV_syst);
          make_syst_branches("diphoton_LV", diphoton_LV_syst);
      }
  }

  //std::cout<<"\n\n\n\n\n\nIn my initialize\n\n\n\n\n\n";
  /* GK 23.03.15 */
  //in initialize
  std::cout<<"--------- option from steering macro: "<< AnalysisBranch<< "----------------\n";
  if(AnalysisBranch==1)
    {
      std::cout<<"            HIGGS                 \n";
    }
  if(AnalysisBranch==11)
    {
      std::cout<<"            EXOTIC                 \n";
    }

  if(AnalysisBranch!=1 && AnalysisBranch!=11 )
    {
      std::cout<<"            you need to specify Analisys                 \n";
      return EL::StatusCode::FAILURE;
    }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode DiphotonAnalysis::postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  HgammaAnalysis::postExecute(); // keep this line!
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode DiphotonAnalysis::finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.
  HgammaAnalysis::finalize(); // keep this line!


  if(photonID_tight_tool_SS){
    delete photonID_tight_tool_SS;
    photonID_tight_tool_SS = 0;
  }
  if(isoCorr_tool){
    delete isoCorr_tool;
    isoCorr_tool = 0;
  }

  cutflow->print_cutflow();

  Info( "Finalize()", "Total initial sum of weight is: %f", weight_sum_before->GetVal());
  Info( "Finalize()", "Total sum of weight is: %f", weight_sum_selected->GetVal());
/*
  Float_t new_weight_sum_selected;
  TBranch *newBranch = sample_NTUP->tree()->Branch("weight_sum_selected", &new_weight_sum_selected);

  Long64_t nentries = sample_NTUP->tree()->GetEntries();
  for (Long64_t i = 0; i < nentries; i++){
    new_weight_sum_selected = weight_sum_selected->GetVal();
    newBranch->Fill();
  }
*/
  sample_NTUP->setFilterPassed();

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode DiphotonAnalysis::histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  HgammaAnalysis::histFinalize(); // keep this line!
  return EL::StatusCode::SUCCESS;
}

template <class T> void DiphotonAnalysis::make_syst_branches(std::string tag, std::vector<T>& vect){
  int sysList_i = 0;
  vect.reserve(sysList.size());
  for(auto sys : sysList){ //= sysList.begin(); sysListItr != sysList.end(); ++sysListItr)
    std::string br_name = sys.name();
    if(sys.name()==""){
      DEBUG("Nominal (no syst)");
      continue;
    } else {
      br_name = tag + "_" + sys.name();
      DEBUG("Var with systematic: " << br_name);
    }
    vect.push_back(T());
    sample_NTUP->tree()->Branch(br_name.c_str(), &(vect[sysList_i]));
    sysList_i++;
  }
}

bool DiphotonAnalysis::cut_PV()
{
  const xAOD::VertexContainer* pvs = 0;
  if( !event()->retrieve( pvs, "PrimaryVertices" ).isSuccess() ){
    Error("execute()", "Failed to retrieve PV container. Exiting." );
    return EL::StatusCode::FAILURE;
  }
  // At least one PV
  unsigned int npv = 0;
  auto vtx_itr = pvs->begin();
  auto vtx_end = pvs->end();
  for( ; vtx_itr != vtx_end; ++vtx_itr ) {
    if( (*vtx_itr)->vertexType() == 1 || (*vtx_itr)->vertexType() == 3 ){
      if( (*vtx_itr)->nTrackParticles() >= 2 ){
      npv++;
      }
    }
  }

  npvs = npv;

  bool pass_PV = (npv >= min_nvertex);

  return pass_PV;
}

void DiphotonAnalysis::clear_variables() {
    DEBUG("\nHSG1Analysis::clear_variables");

    ph_pt.clear();
    ph_eta.clear();
    ph_etas2.clear();
    ph_phi.clear();
    ph_cl_phi.clear();
    ph_parent_pdgID.clear();

    ph_isEM_flags.clear();

    for( unsigned int i = 0; i < isEM_flags_lead.size(); i++ ){
       isEM_flags_lead[i].second = 0;
       isEM_flags_sublead[i].second = 0;
    }

    pass_truth_match = false;
    parent_pdgID_ld = 0;
    parent_pdgID_subld = 0;
    match_leading = false;
    match_subleading = false;
    tight_leading = false;
    my_tight_leading = false;
    tight_subleading = false;

    mass = 0.;
    pt_leading = 0.;

    ph_tight.clear();
    ph_matched.clear();

    topoetcone40.clear();
    topoetcone30.clear();
    topoetcone20.clear();

    ptcone40.clear();
    ptcone30.clear();
    ptcone20.clear();

    ptvarcone40.clear();
    ptvarcone30.clear();
    ptvarcone20.clear();

    mass_gev_syst.clear();
    accepted_syst.clear();

    leading_LV_syst.clear();
    subleading_LV_syst.clear();
    diphoton_LV_syst.clear();

    FixedCutTightCaloOnly_ld = false;
    FixedCutTight_ld = false;
    FixedCutLoose_ld = false;

    FixedCutTightCaloOnly_subld = false;
    FixedCutTight_subld = false;
    FixedCutLoose_subld = false;

    cutflow->reset();
}
