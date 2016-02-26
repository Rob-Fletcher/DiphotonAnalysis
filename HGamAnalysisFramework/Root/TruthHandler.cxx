#include "HGamAnalysisFramework/TruthHandler.h"

#include "HGamAnalysisFramework/PhotonHandler.h"
#include "HGamAnalysisFramework/ElectronHandler.h"
#include "HGamAnalysisFramework/MuonHandler.h"
#include "HGamAnalysisFramework/JetHandler.h"
#include "HGamAnalysisFramework/ETmissHandler.h"
#include "HGamAnalysisFramework/HGamVariables.h"

#include "HGamAnalysisFramework/HgammaIncludes.h"

#include "MCTruthClassifier/MCTruthClassifier.h"

typedef std::pair<MCTruthPartClassifier::ParticleType, MCTruthPartClassifier::ParticleOrigin> ClassifierResult;

namespace HG {

  //______________________________________________________________________________
  SG::AuxElement::Decorator<char> TruthHandler::isIsolated("isIsolated");
  SG::AuxElement::Decorator<float> TruthHandler::etcone20("etcone20");
  SG::AuxElement::Decorator<float> TruthHandler::etcone40("etcone40");
  SG::AuxElement::Decorator<float> TruthHandler::pt("pt");
  SG::AuxElement::Decorator<float> TruthHandler::eta("eta");

  //______________________________________________________________________________
  TruthHandler::TruthHandler(xAOD::TEvent *event, xAOD::TStore *store)
  : m_event(event)
  , m_store(store)
  , m_truthClass(nullptr)
  { }

  //______________________________________________________________________________
  TruthHandler::~TruthHandler()
  {
    SafeDelete(m_truthClass);
  }

  //______________________________________________________________________________
  EL::StatusCode TruthHandler::initialize(Config &config)
  {
    m_MxAODName     = "HGam";

    // Container names
    m_particleName = config.getStr("TruthHandler.ParticleContainerName"   , "TruthParticles"  );
    m_eventName    = config.getStr("TruthHandler.EventContainerName"      , "TruthEvents"     );
    m_photonName   = config.getStr("TruthHandler.PhotonContainerName"     , "TruthPhotons"    );
    m_electronName = config.getStr("TruthHandler.ElectronContainerName"   , "TruthElectrons"  );
    m_muonName     = config.getStr("TruthHandler.MuonContainerName"       , "TruthMuons"      );
    m_jetName      = config.getStr("TruthHandler.JetContainerName"        , "AntiKt4TruthJets");
    m_metName      = config.getStr("TruthHandler.MissingETContainerName"  , "MET_Truth"       );
    m_higgsName    = config.getStr("TruthHandler.HiggsBosonContainerName" , "TruthHiggsBosons");

    // Photon selections
    m_phMaxEta      = config.getNum ("TruthHandler.Photons.MaxAbsEta"          , 2.37);
    m_phMinPt       = config.getNum ("TruthHandler.Photons.PtPreCutGeV"        , 25.0)*HG::GeV;
    m_phRejectCrack = config.getBool("TruthHandler.Photons.ApplyCrackRejection", true);
    m_phMinCrack    = config.getNum ("TruthHandler.Photons.BarrelMaxAbsEta"    , 1.37);
    m_phMaxCrack    = config.getNum ("TruthHandler.Photons.EndcapMinAbsEta"    , 1.52);
    m_phIsoCone     = config.getNum ("TruthHandler.Photons.IsolationCone"      , 0.40);
    m_phIsoCut      = config.getNum ("TruthHandler.Photons.IsolationCutGeV"    , 14.0)*HG::GeV;

    // Electron selections
    m_elMaxEta      = config.getNum ("TruthHandler.Electrons.MaxAbsEta"          , 2.47);
    m_elMinPt       = config.getNum ("TruthHandler.Electrons.PtPreCutGeV"        , 25.0)*HG::GeV;
    m_elRejectCrack = config.getBool("TruthHandler.Electrons.ApplyCrackRejection", true);
    m_elMinCrack    = config.getNum ("TruthHandler.Electrons.BarrelMaxAbsEta"    , 1.37);
    m_elMaxCrack    = config.getNum ("TruthHandler.Electrons.EndcapMinAbsEta"    , 1.52);
    m_elIsoCone     = config.getNum ("TruthHandler.Electrons.IsolationCone"      , 0.40);
    m_elIsoCut      = config.getNum ("TruthHandler.Electrons.IsolationCutGeV"    , 14.0)*HG::GeV;

    // Muon selections
    m_muMaxEta  = config.getNum("TruthHandler.Muons.MaxAbsEta"      , 2.50);
    m_muMinPt   = config.getNum("TruthHandler.Muons.PtPreCutGeV"    , 10.0)*HG::GeV;
    m_muIsoCone = config.getNum("TruthHandler.Muons.IsolationCone"  , 0.40);
    m_muIsoCut  = config.getNum("TruthHandler.Muons.IsolationCutGeV", 14.0)*HG::GeV;

    // Jet selections
    m_jetMaxRapidity = config.getNum("TruthHandler.Jets.MaxAbsRapidity", 4.40);
    m_jetMinPt       = config.getNum("TruthHandler.Jets.PtPreCutGeV"   , 25.0)*HG::GeV;

    // MET selections
    m_metTypes       = config.getStrV("TruthHandler.MissingET.METTypes");

    // MC truth classifier
    m_truthClass = new MCTruthClassifier("MCTruthClassifier");
    if (m_truthClass->initialize().isFailure())
      HG::fatal("Couldn't initialize MCTruthClassifier");

    return EL::StatusCode::SUCCESS;
  }

  //______________________________________________________________________________
  template <>
  void
  TruthHandler::setP4<xAOD::Muon>(const xAOD::IParticle *p, xAOD::Muon &copy)
  { copy.setP4(p->pt(), p->eta(), p->phi()); }

  //______________________________________________________________________________
  template <>
  void
  TruthHandler::setP4<xAOD::TruthParticle>(const xAOD::IParticle *p, xAOD::TruthParticle &copy)
  {
    const xAOD::TruthParticle *truth = static_cast<const xAOD::TruthParticle*>(p);
    if (truth == nullptr)
      fatal("Couldn't cast passed IParticle to TruthParticle, which is not supported. Exiting.");

    copy.setPx(truth->px());
    copy.setPy(truth->py());
    copy.setPz(truth->pz());
    copy.setE (truth->e ());
    copy.setM (truth->m ());
  }

  //______________________________________________________________________________
  template <>
  void
  TruthHandler::setOriginalObjectLink<xAOD::MissingET>(const xAOD::MissingET *orig, xAOD::MissingET *copy)
  {
    // ASG function only works for IParticle, do nothing for MissingET
  }

  //______________________________________________________________________________
  template <>
  void
  TruthHandler::setOriginalObjectLink<xAOD::MissingET>(const DataVector<xAOD::MissingET> *orig, DataVector<xAOD::MissingET> *copy)
  {
    // ASG function only works for DataVector<IParticle>, do nothing for MissingET
  }

  //______________________________________________________________________________
  void TruthHandler::decorateClassification(xAOD::TruthParticle &part)
  {
    // Try to retrieve pointer to original particle, if possible
    const xAOD::TruthParticle *_part = &part;

    static SG::AuxElement::ConstAccessor<ElementLink<xAOD::IParticleContainer> > link("originalObjectLink");
    if (link.isAvailable(part))
      _part = dynamic_cast<const xAOD::TruthParticle*>(xAOD::getOriginalObject(part));

    // Truth type/origin decorations
    static SG::AuxElement::Accessor<int> truthType("truthType");
    static SG::AuxElement::Accessor<int> truthOrigin("truthOrigin");

    // Use MCTruthClassifier to retrieve type/origin info
    ClassifierResult result = m_truthClass->particleTruthClassifier(_part);
    truthType(part) = result.first;
    truthOrigin(part) = result.second;
  }


  //______________________________________________________________________________
  TruthParticles* TruthHandler::getTruthParticles()
  {
    if (!m_event->contains<xAOD::TruthParticleContainer>(m_particleName))
      return nullptr;

    const xAOD::TruthParticleContainer *truthParticles = nullptr;
    if (m_event->retrieve(truthParticles, m_particleName).isFailure())
      HG::fatal("Can't access TruthParticleContainer");

    return truthParticles;
  }

  //______________________________________________________________________________
  xAOD::TruthParticleContainer TruthHandler::getPhotons()
  {
    // Check if HGam truth photons are in TEvent/TStore already
    xAOD::TruthParticleContainer cont;
    if (checkEventAndStore<xAOD::TruthParticle>(cont, m_MxAODName + m_photonName))
      return cont;

    // Get all truth particles, then good photons, from truth record
    TruthParticles *truthParticles = getTruthParticles();
    if (truthParticles == nullptr)
      HG::fatal("No "+m_particleName+" and no "+m_MxAODName+m_photonName+", exiting!");

    TruthContainer  truthcont      = HG::getGoodTruthPhotons(truthParticles);

    // Make a shallow copy, sort by pT
    cont = getDeepCopy<xAOD::TruthParticle>(truthcont.asDataVector(), m_photonName);
    cont.sort(comparePt);

    // Add decorations
    static std::vector<int> ignorePdgIds = {13, 12, 14, 16, 18}; // mu, nus
    for (auto part: cont) {
      etcone20(*part) = HG::getTruthIsolation(part, truthParticles, 0.2, ignorePdgIds);
      etcone40(*part) = HG::getTruthIsolation(part, truthParticles, 0.4, ignorePdgIds);

      isIsolated(*part) = true;
      if (m_phIsoCone > 0 &&
          HG::getTruthIsolation(part, truthParticles, m_phIsoCone) > m_phIsoCut)
        isIsolated(*part) = false;

      pt (*part) = part->pt();
      eta(*part) = part->eta();

      decorateClassification(*part);
    }

    return cont;
  }

  //______________________________________________________________________________
  xAOD::TruthParticleContainer TruthHandler::getElectrons()
  {
    // Check if HGam truth electrons are in TEvent/TStore already
    xAOD::TruthParticleContainer cont;
    if (checkEventAndStore<xAOD::TruthParticle>(cont, m_MxAODName + m_electronName))
      return cont;

    // Get all truth particles, then good electrons, from truth record
    TruthParticles *truthParticles = getTruthParticles();
    if (truthParticles == nullptr)
      HG::fatal("No "+m_particleName+" and no "+m_MxAODName+m_electronName+", exiting!");

    TruthContainer  truthcont      = HG::getGoodTruthElectrons(truthParticles);

    // Make a shallow copy, sort by pT
    cont = getDeepCopy<xAOD::TruthParticle>(truthcont.asDataVector(), m_electronName);
    cont.sort(comparePt);

    // Add decorations
    for (auto part: cont) {
      isIsolated(*part) = true;
      if (m_elIsoCone > 0 &&
          HG::getTruthIsolation(part, truthParticles, m_elIsoCone) > m_elIsoCut)
        isIsolated(*part) = false;

      pt (*part) = part->pt();
      eta(*part) = part->eta();
    }

    return cont;
  }

  //______________________________________________________________________________
  xAOD::TruthParticleContainer TruthHandler::getMuons()
  {
    // Check if HGam truth muons are in TEvent/TStore already
    xAOD::TruthParticleContainer cont;
    if (checkEventAndStore<xAOD::TruthParticle>(cont, m_MxAODName + m_muonName))
      return cont;

    // Get all truth particles, then good muons, from truth record
    TruthParticles *truthParticles = getTruthParticles();
    if (truthParticles == nullptr)
      HG::fatal("No "+m_particleName+" and no "+m_MxAODName+m_muonName+", exiting!");

    TruthContainer  truthcont      = HG::getGoodTruthMuons(truthParticles);

    // Make a shallow copy, sort by pT
    cont = getDeepCopy<xAOD::TruthParticle>(truthcont.asDataVector(), m_muonName);
    cont.sort(comparePt);

    // Add decorations
    for (auto part: cont) {
      isIsolated(*part) = true;
      if (m_muIsoCone > 0 &&
          HG::getTruthIsolation(part, truthParticles, m_muIsoCone) > m_muIsoCut)
        isIsolated(*part) = false;

      pt (*part) = part->pt();
      eta(*part) = part->eta();
    }

    return cont;
  }

  //______________________________________________________________________________
  xAOD::JetContainer TruthHandler::getJets()
  {
    // Check if HGam truth jets are in TEvent/TStore already
    xAOD::JetContainer jets;
    if (checkEventAndStore<xAOD::Jet>(jets, m_MxAODName + m_jetName))
      return jets;

    // Check if raw truth jets are in TEvent/TStore already
    if (!checkEventAndStore<xAOD::Jet>(jets, m_jetName))
      fatal("Couldn't retrieve raw truth jet container, exiting.");

    jets.sort(comparePt);
    return jets;
  }

  //______________________________________________________________________________
  xAOD::MissingETContainer TruthHandler::getMissingET()
  {
    // Check if HGam truth MissingETs are in TEvent/TStore already
    xAOD::MissingETContainer met;
    if (checkEventAndStore<xAOD::MissingET>(met, m_MxAODName + m_metName))
      return met;

    // Check if raw truth MissingETs are in TEvent/TStore already
    if (!checkEventAndStore<xAOD::MissingET>(met, m_metName))
      fatal("Couldn't retrieve raw truth MissingET container, exiting.");

    return met;
  }

  //______________________________________________________________________________
  xAOD::TruthParticleContainer TruthHandler::getHiggsBosons()
  {
    // Check if HGam Higgs Bosons are in TEvent/TStore already
    xAOD::TruthParticleContainer cont;
    if (checkEventAndStore<xAOD::TruthParticle>(cont, m_MxAODName + m_higgsName))
      return cont;

    // Get all truth particles, then good muons, from truth record
    TruthParticles *truthParticles = getTruthParticles();
    if (truthParticles == nullptr)
      HG::fatal("No "+m_particleName+" and no "+m_MxAODName+m_higgsName+", exiting!");

    TruthContainer truthcont = HG::getFinalHiggsBosons(truthParticles);

    // Make a deep copy, sort by pT
    cont = getDeepCopy<xAOD::TruthParticle>(truthcont.asDataVector(), m_higgsName);
    cont.sort(comparePt);

    // Add decorations
    for (auto part: cont) {
      pt (*part) = part->pt();
      eta(*part) = part->eta();
    }

    return cont;
  }

  //______________________________________________________________________________
  xAOD::TruthParticleContainer TruthHandler::applyPhotonSelection(xAOD::TruthParticleContainer &photons)
  {
    xAOD::TruthParticleContainer selected(SG::VIEW_ELEMENTS);
    for (auto ph: photons) {
      // Pt cuts
      if (ph->pt() < m_phMinPt) continue;

      // Eta cuts
      double aeta = fabs(ph->eta());
      if (aeta > m_phMaxEta) continue;
      if (m_phRejectCrack &&
          (aeta > m_phMinCrack && aeta < m_phMaxCrack))
        continue;

      // Isolation cuts
      if (isIsolated.isAvailable(*ph) && !isIsolated(*ph))
        continue;

      // Passed cuts, add to selected container
      selected.push_back(ph);
    }

    return selected;
  }

  //______________________________________________________________________________
  xAOD::TruthParticleContainer TruthHandler::applyElectronSelection(xAOD::TruthParticleContainer &electrons)
  {
    xAOD::TruthParticleContainer selected(SG::VIEW_ELEMENTS);
    for (auto el: electrons) {
      // Pt cuts
      if (el->pt() < m_elMinPt) continue;

      // Eta cuts
      double aeta = fabs(el->eta());
      if (aeta > m_elMaxEta) continue;
      if (m_elRejectCrack &&
          (aeta > m_elMinCrack && aeta < m_elMaxCrack))
        continue;

      // Isolation cuts
      if (isIsolated.isAvailable(*el) && !isIsolated(*el))
        continue;

      // Passed cuts, add to selected container
      selected.push_back(el);
    }

    return selected;
  }

  //______________________________________________________________________________
  xAOD::TruthParticleContainer TruthHandler::applyMuonSelection(xAOD::TruthParticleContainer &muons)
  {
    xAOD::TruthParticleContainer selected(SG::VIEW_ELEMENTS);
    for (auto mu: muons) {
      // Pt cuts
      if (mu->pt() < m_muMinPt) continue;

      // Eta cuts
      if (fabs(mu->eta()) > m_muMaxEta) continue;

      // Isolation cuts
      if (isIsolated.isAvailable(*mu) && !isIsolated(*mu))
        continue;

      // Passed cuts, add to selected container
      selected.push_back(mu);
    }

    return selected;
  }

  //______________________________________________________________________________
  xAOD::JetContainer TruthHandler::applyJetSelection(xAOD::JetContainer &jets)
  {
    xAOD::JetContainer selected(SG::VIEW_ELEMENTS);
    for (auto jet: jets) {
      // Pt cuts
      if (jet->pt() < m_jetMinPt) continue;

      // Eta cuts
      if (fabs(jet->rapidity()) > m_jetMaxRapidity) continue;

      // Passed cuts, add to selected container
      selected.push_back(jet);
    }

    return selected;
  }

  //______________________________________________________________________________
  xAOD::MissingETContainer TruthHandler::applyMissingETSelection(xAOD::MissingETContainer &mets)
  {
    xAOD::MissingETContainer selected(SG::VIEW_ELEMENTS);
    for (auto met: mets) {
      // Limit MET to the types specified in config (TST, CST, ...)
      if (std::find(m_metTypes.begin(), m_metTypes.end(), met->name().c_str()) == m_metTypes.end())
        continue;

      // Passed cuts, add to selected container
      selected.push_back(met);
    }

    return selected;
  }

  //______________________________________________________________________________
  void TruthHandler::removeOverlap(xAOD::TruthParticleContainer &photons  ,
                                   xAOD::JetContainer           &jets     ,
                                   xAOD::TruthParticleContainer &electrons,
                                   xAOD::TruthParticleContainer &muons    )
  {
    // jets overlapping a pT>15 GeV photon or electron are removed
    for ( auto jet=jets.rbegin(); jet!=jets.rend(); ++jet) {
      bool overlap=false;
      for (auto gam:photons)
        if (gam->pt()>15.0*HG::GeV && HG::DRrap(gam,*jet)<0.4)
          overlap=true;
      for (auto e:electrons)
        if (e->pt()>15.0*HG::GeV && HG::DRrap(e,*jet)<0.4)
          overlap=true;
      if (overlap) jets.erase(jet.base()-1);
    }
  }

  //______________________________________________________________________________
  double TruthHandler::truthVertexZ()
  {
    if (var::truthVertexZ.exists())
      return var::truthVertexZ();

    const xAOD::TruthEventContainer *truthEvents = nullptr;
    if (m_event->retrieve(truthEvents, m_eventName).isFailure())
      HG::fatal("Can't access TruthEvents");

    if (truthEvents->size() < 1) {
      Warning("TruthHandler::truthVertexZ()","No TruthEvents?");
      return -999;
    }

    static int nNullVtx = 0;
    if (truthEvents->at(0)->signalProcessVertex() == nullptr) {
      nNullVtx++;
      if (nNullVtx < 5)
        Warning("TruthHandler::truthVertexZ()","No signalProcessVertex for event");
      if (nNullVtx == 5)
        Warning("TruthHandler::truthVertexZ()","Supressing WARNING for: No signalProcessVertex for event");
      return -999;
    }

    var::truthVertexZ.setValue(truthEvents->at(0)->signalProcessVertex()->z());
    return var::truthVertexZ();
  }

  //______________________________________________________________________________
  bool TruthHandler::passFiducial(const xAOD::TruthParticleContainer *allPhotons,
                                  const xAOD::TruthParticleContainer *electrons ,
                                  const xAOD::TruthParticleContainer *muons     ,
                                  const xAOD::JetContainer           *jets      )
  {
    if (var::isFiducial.exists())
      return var::isFiducial();

    var::isFiducial.setTruthValue(false);

    // Kinematic cuts
    if (not passFiducialKinOnly(allPhotons, electrons, muons, jets))
      return false;

    const xAOD::TruthParticle *gam1 = (*allPhotons)[0], *gam2 = (*allPhotons)[1];

    // Isolation cuts
    if (not isIsolated(*gam1) || not isIsolated(*gam2))
      return false;

    var::isFiducial.setTruthValue(true);

    // All cuts passed, this event is in fiducial volume
    return true;
  }
  //______________________________________________________________________________
  bool TruthHandler::passFiducialKinOnly(const xAOD::TruthParticleContainer *allPhotons,
                                         const xAOD::TruthParticleContainer *electrons ,
                                         const xAOD::TruthParticleContainer *muons     ,
                                         const xAOD::JetContainer           *jets      )
  {
    if (var::isFiducialKinOnly.exists())
      return var::isFiducialKinOnly();

    var::isFiducialKinOnly.setTruthValue(false);

    // Safety check
    if (allPhotons == nullptr)
      return false;

    // At least two photons
    if (allPhotons->size() < 2)
      return false;

    const xAOD::TruthParticle *gam1 = (*allPhotons)[0], *gam2 = (*allPhotons)[1];

    // Pt cuts
    if (gam1->pt() < m_phMinPt || gam2->pt() < m_phMinPt)
      return false;

    // Eta cuts
    double aeta1 = fabs(gam1->eta()), aeta2 = fabs(gam2->eta());
    if (aeta1 > m_phMaxEta || aeta2 > m_phMaxEta)
      return false;

    if (m_phRejectCrack                                  &&
        ((aeta1 > m_phMinCrack && aeta1 < m_phMaxCrack)  ||
         (aeta2 > m_phMinCrack && aeta2 < m_phMaxCrack)) )
      return false;

    // Relative pT cuts
    TLorentzVector yy = gam1->p4() + gam2->p4();
    if (gam1->pt()/yy.M() < 0.35 || gam2->pt()/yy.M() < 0.25)
      return false;

    var::isFiducialKinOnly.setTruthValue(true);

    // All cuts passed, this event is in kinematic fiducial volume
    return true;
  }

  //______________________________________________________________________________
  EL::StatusCode TruthHandler::writePhotons(xAOD::TruthParticleContainer &container)
  { return writeContainer<xAOD::TruthParticle>(container, m_MxAODName + m_photonName); }

  //______________________________________________________________________________
  EL::StatusCode TruthHandler::writeElectrons(xAOD::TruthParticleContainer &container)
  { return writeContainer<xAOD::TruthParticle>(container, m_MxAODName + m_electronName); }

  //______________________________________________________________________________
  EL::StatusCode TruthHandler::writeMuons(xAOD::TruthParticleContainer &container)
  { return writeContainer<xAOD::TruthParticle>(container, m_MxAODName + m_muonName); }

  //______________________________________________________________________________
  EL::StatusCode TruthHandler::writeJets(xAOD::JetContainer &container)
  { return writeContainer<xAOD::Jet>(container, m_MxAODName + m_jetName); }

  //______________________________________________________________________________
  EL::StatusCode TruthHandler::writeMissingET(xAOD::MissingETContainer &container)
  { return writeContainer<xAOD::MissingET>(container, m_MxAODName + m_metName); }

  //______________________________________________________________________________
  EL::StatusCode TruthHandler::writeHiggsBosons(xAOD::TruthParticleContainer &container)
  { return writeContainer<xAOD::TruthParticle>(container, m_MxAODName + m_higgsName); }

  //______________________________________________________________________________
  bool TruthHandler::comparePt(const xAOD::IParticle *a, const xAOD::IParticle *b)
  { return a->pt() > b->pt(); }

} // namespace HG
