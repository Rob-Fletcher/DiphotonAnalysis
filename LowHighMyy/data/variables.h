// This file contains photon, jets, electrons, muons and ETmiss variables
// In the MC sample all events are kept (in order to compute CX correction factors) whereas for data only events after the preselection is applied are recorded. In both cases the reconstructed variables are computed only after preselection.
// Different flags allow to distinguish different steps in the cut-flow (preselection, tight, isolation, etc). The most important one is probably isPassingRelative: all selection cuts for the final analysis are applied.
// different weights are to be applied to the MC samples


///////////////////////////////////////////////////////////////////////
//                        Photon variables                           //
///////////////////////////////////////////////////////////////////////

// information for the two leading photons after pre-selection is recorded: ph1_* and ph2_*

x Float_t         ph1_Et; // transverse energy in GeV
x Float_t         ph1_eta; // eta corrected to the diphoton vertex
x Float_t         ph1_phi; // phi 
x Float_t         ph1_etas2; // eta in second layer of EM calo
Float_t         ph1_time; // timing
x Float_t         ph1_EDmedian; // energy density for pile-up correction
Int_t           ph1_index; // index in the Photon container
// isolation
x Float_t         ph1_isol_calo; // calorimeter isolation in GeV
Float_t         ph1_isol_track; // track isolation in GeV
// conversion
Float_t         ph1_zvertex; // zvertex
x Float_t         ph1_Rconv; // R position of conversion vertex (if any)
Float_t         ph1_zconv; // z position of conversion vertex (if any)
Float_t         ph1_materialTraversed; // ? (probably related to conversion)
x Int_t           ph1_convFlag; // conversion flag (ph1_convFlag%10==0 for unconverted)
x Int_t           ph1_singleTrackConv; // boolean if single track conversion
x Int_t           ph1_convTrk1_nPixelHits; // conversion stuff
x Int_t           ph1_convTrk2_nPixelHits;
x Int_t           ph1_convTrk1_nSCTHits;
x Int_t           ph1_convTrk2_nSCTHits;
x Int_t           ph1_convTrk1_nTRTHits;
x Int_t           ph1_convTrk2_nTRTHits;
// identification
UInt_t          ph1_isEM; // isEM bit from PhotonID tool
UInt_t          ph1_isEM_noFudge; // isEM bit from PhotonID tool without FF applied
Int_t           ph1_isLoose; // loose ID from PhotonID tool
Int_t           ph1_isTight; // tight ID from PhotonID tool
Int_t           ph1_isLoose_noFudge; // loose ID from PhotonID tool without FF applied
Int_t           ph1_isTight_noFudge;; // tight ID from PhotonID tool without FF applied
Float_t         ph1_fside; // shower-shapes
Float_t         ph1_ws3;
Float_t         ph1_weta2;
Float_t         ph1_wstot;
Float_t         ph1_rHad1;
Float_t         ph1_rHad;
Float_t         ph1_deltaE;
Float_t         ph1_eRatio;
Float_t         ph1_e277;
Float_t         ph1_reta;
Float_t         ph1_rphi;
Bool_t          ph1_passfSide; // single-bit efficiencies (after loose): is passing a given cut
Bool_t          ph1_passws3;
Bool_t          ph1_passweta2;
Bool_t          ph1_passwstot;
Bool_t          ph1_passrHad;
Bool_t          ph1_passdeltaE;
Bool_t          ph1_passeRatio;
Bool_t          ph1_passreta;
Bool_t          ph1_passrphi;
Float_t         ph1_Et_true; // true transverse energy in GeV
Float_t         ph1_eta_true; // true eta
Float_t         ph1_phi_true; // true phi
Int_t           ph1_type; // truth type
Int_t           ph1_mother; // truth mother type
Int_t           ph1_truthConv; // is true conversion
Int_t           ph1_isPhotonFromHardProc; // truth: is the photon from a hard process
Int_t           ph1_isBrem; // truth: is brem photon
Int_t           ph1_isPromptPhoton; // truth: is a prompt photon

///////////////////////////////////////////////////////////////////////
//                   Truth photon variables                          //
///////////////////////////////////////////////////////////////////////

// information on the two leading photons from the egammaTruth container. They are used to compute the CX correction factors. 
// !! Because they came from a different container than the reco photon variables there is no reason that truePhoton1_Et and ph1_Et_true are the same

x Float_t         truePhoton1_Et; // true ET in GeV
x Float_t         truePhoton1_eta; // eta
x Float_t         truePhoton1_phi; // phi
Int_t           truePhoton1_motherType; // mother type
Float_t         truePhoton1_particleIsol20_ED_corr; // particle isolation in a ΔR=0.x cone with pile-up correction
Float_t         truePhoton1_particleIsol30_ED_corr;
Float_t         truePhoton1_particleIsol40_ED_corr;
Float_t         truePhoton1_particleIsol20; // particle isolation in a ΔR=0.x cone without pile-up correction
Float_t         truePhoton1_particleIsol30;
Float_t         truePhoton1_particleIsol40;
Float_t         truePhoton1_particleIsol_fromTruth10; // ?
Float_t         truePhoton1_particleIsol_fromTruth20;
Float_t         truePhoton1_particleIsol_fromTruth40;
Float_t         truePhoton1_particleEDmedian; // energy density for pile-up correction

///////////////////////////////////////////////////////////////////////
//                         Jet variables                             //
///////////////////////////////////////////////////////////////////////

// information for the three leading jets after selection (pt, cleaning, overlap removal, etc): jet1_* and jet2_* and jet3_*

Float_t         jet1_pt; // jet pt in GeV
Float_t         jet1_eta; // jet eta
Float_t         jet1_phi; // jet phi
Float_t         jet1_E; // jet energy
Float_t         jet1_btag; // b-tagging weight (flavor_weight_MV1)


///////////////////////////////////////////////////////////////////////
//                         Jet variables                             //
///////////////////////////////////////////////////////////////////////

// information on the two leading truth jets from the TruthJet variables (same warning as for truePhoton*): truthJet1_* and truthJet2_*

Float_t         truthJet1_pt;
Float_t         truthJet1_eta;
Float_t         truthJet1_phi;
Float_t         truthJet1_E;
Float_t         truthJet1_btag;


///////////////////////////////////////////////////////////////////////
//                       Electron variables                          //
///////////////////////////////////////////////////////////////////////

// information for the two leading electrons after selection (pt, medium ID, overlap removal, etc): el1_* and el2_* 

Float_t         el1_pt; // electron pt in GeV
Float_t         el1_eta; // electron eta
Float_t         el1_phi; // electron phi
Float_t         el1_sf; // electron SF (reconstruction*medium)
Float_t         el1_sf_errreco; // error on electron reco SF
Float_t         el1_sf_errid; // error on electron medium SF
Float_t         el1_sf_tight; // electron SF for tight
Int_t           el1_isTightPP; // isPassing tight ID
Float_t         el1_et_up; // energy variations with different energy scale and resolution uncertainties
Float_t         el1_et_down;
Float_t         el1_et_ZeeStatUp;
Float_t         el1_et_ZeeMethUp;
Float_t         el1_et_ZeeGenUp;
Float_t         el1_et_MATUp;
Float_t         el1_et_PSUp;
Float_t         el1_et_LowUp;
Int_t           el1_mother; // truth : electron mother type
Int_t           el1_type; // truth : type

///////////////////////////////////////////////////////////////////////
//                         Muon variables                            //
///////////////////////////////////////////////////////////////////////

// information for the two leading (staco) muons after selection (pt, number of hits, ID, overlap removal, etc): mu1_* and mu2_* 

Float_t         mu1_pt; // pt in GeV
Float_t         mu1_eta; // eta
Float_t         mu1_phi; // phi
Float_t         mu1_sf; // SF
Float_t         mu1_sf_err; //  error on SF
Float_t         mu1_d0; // d0
Float_t         mu1_pt_msup; // variations for uncertainties
Float_t         mu1_pt_msdown;
Float_t         mu1_pt_idup;
Float_t         mu1_pt_iddown;
Int_t           mu1_isLoose; // loose ID
Int_t           mu1_isTight; // tight ID
Int_t           mu1_veto; // veto if overlap with electron
Int_t           mu1_mother; // truth: mother type
Int_t           mu1_type; // truth:  type

///////////////////////////////////////////////////////////////////////
//                        ETmiss variables                           //
///////////////////////////////////////////////////////////////////////

// I won't go through all of them, each time you have the Ex and Ey projections in GeV, so you can compute ETmiss = sqrt(EXmiss^2+EYmiss^2) and phimiss = atan(EYmiss/EXmiss)

Float_t         MET_RefFinal_etx;
Float_t         MET_RefFinal_ety;
Float_t         MET_RefFinal_sumet;
Float_t         MET_Topo_etx;
Float_t         MET_Topo_ety;
Float_t         MET_Topo_sumet;
Float_t         MET_LocHadTopo_etx;
Float_t         MET_LocHadTopo_ety;
Float_t         MET_LocHadTopo_sumet;
Float_t         MET_PhotonTight_Calib_OR_stdvert_RefFinal_etx; // the one that was used for the H->γγ analysis
Float_t         MET_PhotonTight_Calib_OR_stdvert_RefFinal_ety;
Float_t         MET_PhotonTight_Calib_OR_stdvert_RefFinal_sumet;
Float_t         MET_Truth_NonInt_etx;
Float_t         MET_Truth_NonInt_ety;
Float_t         MET_Truth_NonInt_sumet;
Float_t         MET_Truth_Int_etx;
Float_t         MET_Truth_Int_ety;
Float_t         MET_Truth_Int_sumet;
Float_t         MET_Truth_IntMuons_etx;
Float_t         MET_Truth_IntMuons_ety;
Float_t         MET_Truth_IntMuons_sumet;


///////////////////////////////////////////////////////////////////////
//              General information on the event                     //
///////////////////////////////////////////////////////////////////////

Int_t           RunNumber;
Int_t           LB;
Int_t           EventNumber;
Bool_t          collcand_passCaloTime; // it is so 2009
Bool_t          collcand_passMBTSTime;
Int_t           npv;
Int_t           trk_n; // size of track container (do we need this???)
Int_t           mc_channel_number;
Float_t         averageIntPerXing;
Float_t         trueZvertex;
Float_t         zvertex;

///////////////////////////////////////////////////////////////////////
//                diphoton and dijet variables                       //
///////////////////////////////////////////////////////////////////////

Float_t         mgg; // diphoton invariant mass
Float_t         cosThStar; // cosθ*
Float_t         cosThStarCS; // cosθ* in Collin-Sopper frame
Float_t         pT; // pt of the diphoton system
Float_t         pTt; // pT_t of the diphoton system
Float_t         pTl; // pT_l of the diphoton system
Float_t         pT_true; // true pt of the diphoton system
Float_t         pTt_true; // true pT_t of the diphoton system
Float_t         eta; // eta of the diphoton system
Float_t         phi; // phi of the diphoton system
Float_t         dpT; // ΔpT between the two photons
Float_t         dEta; // Δη between the two photons
Float_t         dPhi; // Δφ between the two photons
Float_t         dR; // ΔR between the two photons
Float_t         mggj; // invariant mass of the diphoton+1 jet system
Float_t         pTggj; // pT of the diphoton+1 jet system
Float_t         pTjj; // pT of the dijet system
Float_t         mjj; // dijet invariant mass
Float_t         delta_phi_jj; // Δφ between the two jets
Float_t         delta_phi_gg_jj; // Δφ between the diphoton and the dijet systems
Float_t         delta_eta_gg_jj; // Δη between the diphoton and the dijet systems
Float_t         cosThStar_gg_jj; // cosθ* with the diphoton and the dijet systems
Float_t         jet_DeltaR15; // minimum ΔR between any jet with pT> x GeV and any of the two photons
Float_t         jet_DeltaR20;
Float_t         jet_DeltaR25;
Float_t         jet_DeltaR30;
Float_t         jet_DeltaR50;
Float_t         jet1_DeltaR30; // minimum ΔR between the leading jet with pT> 20 GeV and any of the two photons
Float_t         jet12_DeltaR30; // minimum ΔR any of the two leading jets with pT> 20 GeV and any of the two photons
Float_t         truthJet1_DeltaR30; // same with truth variables
Float_t         truthJet12_DeltaR30;
Float_t         MET_Signif; // ETmiss significance
Float_t         HT; // HT computed with the two leading photons, all the jets, the two leading electrons and muons (after selection cuts)
Float_t         HT_woPh; // the same but without photons
Float_t         sphericity; // sphericity (same objects as HT)
Float_t         sphericity_T; // transverse sphericity (same objects as HT)

Float_t         mee; // dielectron invariant mass
Float_t         mTe; // electron+Etmiss transverse mass
Float_t         mmumu; // same for muons
Float_t         mTmu;



///////////////////////////////////////////////////////////////////////
//                               weights                             //
///////////////////////////////////////////////////////////////////////


Float_t         weight_allButLumi; // product of all weights (MC, interference, zvertex, pileup. pt)
Float_t         weight_lumi; // lumi weight ( = XS*FiltEff/Lumi/Ntot) to normalised to 1 fb-1
Double_t        weight; // total weight:  weight_allButLumi*weight_lumi
Float_t         weight_mc; // MC weight
Float_t         weight_interference; // interference
Float_t         weight_zvertex; // correction to diphoton vertex
Float_t         weight_pileup; // pileup
Float_t         weight_pt; // correction of Higgs pT for ggF
Float_t         weight_electron; // electron SF
Float_t         weight_muon; // muon SF


///////////////////////////////////////////////////////////////////////
//                             Flags                                 //
///////////////////////////////////////////////////////////////////////

// allow to distinguish different steps in the cut-flow (preselection, tight, isolation, etc).

Bool_t          isPassing; // all selection cuts except the ET/mγγ ones, and with 100 < mγγ < 160 GeV
Bool_t          isPreselected; // preselection
Bool_t          isAllButIsol; // all selection cuts except the isolation cuts
Bool_t          isAllButMass; // all selection cuts except the ET/mγγ ones, no cut on mγγ
Bool_t          isPassingRelative; // all selection cuts with the ET/mγγ ones, no cut on mγγ
Bool_t          isVBFCutsLoose; // true if the event passes the criteria to enter the loose cut-based category
Bool_t          isVBFCutsTight;; // true if the event passes the criteria to enter the tight cut-based category
Bool_t          isVBFMVALoose;; // true if the event passes the criteria to enter the loose MVA category
Bool_t          isVBFMVATight;; // true if the event passes the criteria to enter the tight MVA category
Bool_t          isVHad; // true if the event passes the criteria to enter the hadronic VH category
Bool_t          isMET; // true if the event passes the criteria to enter the ETmiss category
Bool_t          isElectron; // true if the event has an electron passing the selection cuts
Bool_t          isMuon; // true if the event has a muon passing the selection cuts (lepton category = isElectron&&isMuon)

///////////////////////////////////////////////////////////////////////
//                  Number of selected objects                       //
///////////////////////////////////////////////////////////////////////

// Number of objects passing the selection cuts even if we keep detailed information on the leading ones only

Int_t           nJets; // number of jets passing the default cuts
Int_t           nJets30; // number of jets passing the default cuts + having pt > 30 GeV
Int_t           nJetsEleLoose; // number of loose electrons
Int_t           nJetsEleMedium; // number of medium electrons
Int_t           nJetsMu; //  number of muons

///////////////////////////////////////////////////////////////////////
//                              Higgs boson                          //
///////////////////////////////////////////////////////////////////////

// true kinematic variables for the Higgs boson

Float_t         trueHiggs_m;
Float_t         trueHiggs_pt;
Float_t         trueHiggs_eta;
Float_t         trueHiggs_phi;

///////////////////////////////////////////////////////////////////////
//                            Categories                             //
///////////////////////////////////////////////////////////////////////

// category number for different categories used in the past (ask me if you want the detail of it)

Int_t           cat14MoriondMVA_idx; // categories for Moriond2013 with MVA VBF categories
Int_t           cat14MoriondCuts_idx; // categories for Moriond2013 with cut-based VBF categories
Int_t           cat10Moriond_idx; // categories for Moriond2013 without the VH categories
Int_t           catConvEtapTt_idx; // older categories
Int_t           catConvEtapTt2Jets_idx;





///////////////////////////////////////////////////////////////////////
//                             Misc                                  //
///////////////////////////////////////////////////////////////////////

// I don't really know ;)

Int_t           is_true_Zqq;
Int_t           is_true_Zee;
Int_t           is_true_Zmm;
Int_t           is_true_Ztt;
Int_t           is_true_Znn;
Int_t           is_true_Wqq;
Int_t           is_true_Wen;
Int_t           is_true_Wmn;
Int_t           is_true_Wtn;
