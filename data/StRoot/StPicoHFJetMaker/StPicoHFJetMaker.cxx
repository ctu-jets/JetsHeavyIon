#include "StPicoHFJetMaker.h"
#include "JetInfo.h"

#include "BemcNewCalib.h"
#include "StEmcADCtoEMaker/StBemcData.h"
#include "StEmcADCtoEMaker/StEmcADCtoEMaker.h"
#include "StEmcRawMaker/StBemcRaw.h"
#include "StEmcRawMaker/StBemcTables.h"
#include "StEmcRawMaker/defines.h"

#include "TRandom3.h"
#include "TVector2.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/config.h"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "MyJet.h"

using namespace std;

vector<MatchedJetPair> MatchJetsEtaPhi(const vector<MyJet> &McJets,
                                       const vector<MyJet> &RecoJets,
                                       const double &R);

ClassImp(StPicoHFJetMaker)

    StPicoHFJetMaker::StPicoHFJetMaker(TString name, StPicoDstMaker *picoMaker,
                                       TString outputBaseFileName)
    : StPicoJetMaker(name, picoMaker, outputBaseFileName),
      mRefmultCorrUtil(NULL) {

  // constructor
}

// _________________________________________________________
StPicoHFJetMaker::~StPicoHFJetMaker() {
  // destructor
}

// _________________________________________________________
int StPicoHFJetMaker::InitJets() {

  mADCtoEMaker = dynamic_cast<StEmcADCtoEMaker *>(GetMaker("Eread"));
  assert(mADCtoEMaker);
  mTables = mADCtoEMaker->getBemcData()->getTables();

  // -- INITIALIZE USER HISTOGRAMS ETC HERE -------------------
  //    add them to the output list mOutList which is automatically written

  // EXAMPLE //  mOutList->Add(new TH1F(...));
  // EXAMPLE //  TH1F* hist = static_cast<TH1F*>(mOutList->Last());
  TH1::SetDefaultSumw2();

  mOutList->Add(new TH1D("hcent", "centrality", 10, -1, 9));
  TDirectory *currentDir = gDirectory;
  // create directory for each radius
  for (unsigned int i = 0; i < fR.size(); i++) {
    TString dirName = Form("R%.1f", fR[i]);
    TDirectory *dir = currentDir->mkdir(dirName);
    dir->cd();
    TTree *jetTree = new TTree("JetTree", "JetTree");
    jetTree->Branch("runId", &fRunNumber, "runId/I");
    jetTree->Branch("centrality", &fCentrality, "centrality/I");
    jetTree->Branch("centralityWeight", &fCentralityWeight,
                    "centralityWeight/F");


    if(mIsEmbedding){
    jetTree->Branch("xsecWeight", &fXsecWeight, "xsecWeight/F");
    jetTree->Branch("deltaR", &fDeltaR, "deltaR/F");
    jetTree->Branch("mc_pt", &fMcJet.pt, "mc_pt/F");
    jetTree->Branch("mc_eta", &fMcJet.eta, "mc_eta/F");
    jetTree->Branch("mc_phi", &fMcJet.phi, "mc_phi/F");
    jetTree->Branch("mc_area", &fMcJet.area, "mc_area/F");
    jetTree->Branch("mc_pt_lead", &fMcJet.pt_lead, "mc_pt_lead/F");
    jetTree->Branch("mc_n_constituents", &fMcJet.n_constituents,
                    "mc_n_constituents/I");
    jetTree->Branch("mc_neutral_fraction", &fMcJet.neutral_fraction,
                    "mc_neutral_fraction/F");
    }
    jetTree->Branch("reco_pt", &fRecoJet.pt, "reco_pt/F");
    jetTree->Branch("reco_pt_corr", &fRecoJet.pt_corr, "reco_pt_corr/F");
    jetTree->Branch("reco_eta", &fRecoJet.eta, "reco_eta/F");
    jetTree->Branch("reco_phi", &fRecoJet.phi, "reco_phi/F");
    jetTree->Branch("reco_area", &fRecoJet.area, "reco_area/F");
    jetTree->Branch("reco_rho", &fRecoJet.rho, "reco_rho/F");
    jetTree->Branch("reco_pt_lead", &fRecoJet.pt_lead, "reco_pt_lead/F");
    jetTree->Branch("reco_n_constituents", &fRecoJet.n_constituents,
                    "reco_n_constituents/I");
    jetTree->Branch("reco_neutral_fraction", &fRecoJet.neutral_fraction,
                    "reco_neutral_fraction/F");
    jetTree->Branch("reco_trigger_match", &fRecoJet.trigger_match,
                    "reco_trigger_match/O");

    fTree.push_back(jetTree);

    currentDir->cd();
  }

  return kStOK;
}

// _________________________________________________________
void StPicoHFJetMaker::ClearJets(Option_t *opt = "") { return; }

// _________________________________________________________
int StPicoHFJetMaker::FinishJets() {
  TDirectory *currentDir = gDirectory;
  for (unsigned int i = 0; i < fTree.size(); i++) {
    TString dirName = Form("R%.1f", fR[i]);
    TDirectory *dir = (TDirectory *)currentDir->Get(dirName);
    dir->cd();
    if (fTree[i]) {
      fTree[i]->Write();
    }
    currentDir->cd();
  }

  return kStOK;
}

// _________________________________________________________
int StPicoHFJetMaker::MakeJets() {
  TH1F *hcent = static_cast<TH1F *>(mOutList->FindObject("hcent"));

  vector<fastjet::PseudoJet> jetTracks;
  vector<fastjet::PseudoJet> neutraljetTracks; // from bemc towers only
  vector<fastjet::PseudoJet> fullTracks;
  vector<fastjet::PseudoJet> MCjetTracks;
  // vector<fastjet::PseudoJet> MCjetTowers;

  fRunNumber = mPicoDst->event()->runId();
  int eventId = mPicoDst->event()->eventId(); // eventID
  int refMult = mPicoDst->event()->refMult();
  // get centrality
  double vz = mPrimVtx.z();
  mRefmultCorrUtil->setEvent(fRunNumber, refMult, mPicoDst->event()->ZDCx(),
                             vz);
  fCentrality = mRefmultCorrUtil->centrality9(); // 0 = 0-5 %,..., 8 = 70-80
                                                 // %
  if (fCentrality == -1)
    return kStOk; // no fCentrality
  if (fCentrality == 0)
    fCentrality = 1; // merge 0-5% and 5-10% into 0-10%
  if (fCentrality == 8)
    fCentrality = 7; // merge 60-70% and 70-80% into 60-80%

  fCentralityWeight = mRefmultCorrUtil->weight();

  hcent->Fill(fCentrality, fCentralityWeight);

  // MC tracks
  int noMCtracks = mPicoDst->numberOfMcTracks();
  for (int i = 0; i < noMCtracks; i++) {
    StPicoMcTrack *mctrk = (StPicoMcTrack *)mPicoDst->mcTrack(i);
    if (mctrk->idVtxStart() > 1)
      continue; // only primary tracks
    int geantId = mctrk->geantId();
    double mcpt = mctrk->pt();
    double mceta = mctrk->eta();
    if ((geantId > 3 && geantId < 7) || fabs(mceta) > 1.0 || mcpt < 0.2)
      continue;
    TVector3 mcmom = mctrk->p();
    double mcphi = mcmom.Phi();
    if (mcphi < 0.0)
      mcphi += 2.0 * TMath::Pi();
    if (mcphi > 2.0 * TMath::Pi())
      mcphi -= 2.0 * TMath::Pi();

    double mcpx, mcpy, mcpz;
    mcpx = mcmom.x();
    mcpy = mcmom.y();
    mcpz = mcmom.z();

    double mcE = mctrk->energy();

    fastjet::PseudoJet inputMcParticle(mcpx, mcpy, mcpz, mcE);
    if (mctrk->charge() == 0) {
      inputMcParticle.set_user_index(0);
    } else
      inputMcParticle.set_user_index(i);
    MCjetTracks.push_back(inputMcParticle);
  }

  // RC part
  GetCaloTrackMomentum(mPicoDst, mPrimVtx); // fill array Sump with momenta of
                                            // tracks which are matched to BEMC

  StEmcPosition *mEmcPosition = new StEmcPosition();

  double TOWE = 0;
  for (int iTow = 0; iTow < 4800; iTow++) { // get btow info
    StPicoBTowHit *towHit = mPicoDst->btowHit(iTow);
    if (!towHit || towHit->isBad())
      continue; // if the tower is marked as bad or missing info
    int realtowID = towHit->numericIndex2SoftId(iTow);
    if (BadTowerMap[realtowID])
      continue; // exclude bad towers (map in JetInfo.h)

    double towE = GetTowerCalibEnergy(iTow + 1); // get tower energy
    TOWE = towE; // just keep track of the original energy for trigger
                 // approximation

    if (doTowErrPlus == true) {
      towE = towE + 0.038 * towE;
    }
    if (doTowErrMinus == true) {
      towE = towE - 0.038 * towE;
    }

    towE -= fHadronCorr * Sump[iTow]; // subtract hadronic energy deposition
    if (towE < 0)
      towE = 0;

    StEmcGeom *mEmcGeom;
    mEmcGeom = StEmcGeom::getEmcGeom("bemc");
    float Toweta_tmp = 0, Towphi = 0;
    mEmcGeom->getEtaPhi(realtowID, Toweta_tmp, Towphi);
    StThreeVectorF towerPosition = mEmcPosition->getPosFromVertex(
        StThreeVectorF(mPrimVtx.x(), mPrimVtx.y(), mPrimVtx.z()), realtowID);
    //    float Toweta2 = vertexCorrectedEta(Toweta_tmp, vz); //max eta 1.05258
    //    max difference: ET = 0.124452 for E = 0.2, if we cut on |Vz| < 30 cm
    float Toweta = towerPosition.pseudoRapidity();
    double ET = towE / cosh(Toweta);
    if (ET > 30) {
      continue;
    } // ignore E > 30 GeV towers
    // no clustering
    double px, py, pz;

    px = ET * cos(Towphi);
    py = ET * sin(Towphi);
    pz = towE * tanh(Toweta);

    fastjet::PseudoJet inputTower(px, py, pz, towE);
    if (inputTower.perp() > fETmincut) {
      inputTower.set_user_index(
          0); // default index is -1, 0 means neutral particle

      int ADC = towHit->adc() >> 4;
      if (ADC > fTrgthresh)
        inputTower.set_user_index(
            9999); // mark trigger towers with user_index 9999
      neutraljetTracks.push_back(inputTower);
    }
  } // end get btow info

  TRandom3 randGen;
  // loop over primary tracks
  for (unsigned int i = 0; i < mIdxPicoParticles.size(); i++) {
    StPicoTrack *trk = mPicoDst->track(mIdxPicoParticles[i]);
    if (doTrackErr == true) {
      double randomNumber = randGen.Rndm();
      if (randomNumber > 0.96) {
        continue;
      }
    }
    double pT = trk->pMom().Perp(); // using primary tracks
    if (pT != pT)
      continue; // NaN test.
    float eta = trk->pMom().PseudoRapidity();
    float phi = trk->pMom().Phi();
    float dca = (mPrimVtx - trk->origin()).Mag();
    float charged = trk->charge();

    fastjet::PseudoJet inputParticle(trk->pMom().x(), trk->pMom().y(),
                                     trk->pMom().z(), trk->pMom().Mag());
    if (trk->qaTruth() > 95) {
      inputParticle.set_user_index(trk->idTruth() - 1);
      jetTracks.push_back(inputParticle);
    }
  } // end loop over primary tracks

  fullTracks = neutraljetTracks;
  fullTracks.insert(
      fullTracks.end(), jetTracks.begin(),
      jetTracks.end()); // commenting this line will cause only neutral jets,
  // MAX NEUTRAL FRACTION HAS TO BE TURNED OFF

//==================================================================================//
// Jet part
//==================================================================================//
fastjet::AreaDefinition area_def(
    fastjet::active_area_explicit_ghosts,
    fastjet::GhostedAreaSpec(fGhostMaxrap, 1, 0.01));

//====================background estimate=======================//
fastjet::JetDefinition jet_def_for_rho(fastjet::kt_algorithm, fRBg);
if (fCentrality == 1)
  nJetsRemove = 2; // remove 2 hardest jets in central

fastjet::Selector selector = (!fastjet::SelectorNHardest(nJetsRemove)) *
                             fastjet::SelectorAbsEtaMax(1.0) *
                             fastjet::SelectorPtMin(0.01);

fastjet::JetMedianBackgroundEstimator bkgd_estimator(
    selector, jet_def_for_rho, area_def);
bkgd_estimator.set_particles(fullTracks);
float rho = bkgd_estimator.rho();
//======================================================================//

for (unsigned int i = 0; i < fR.size(); i++) {
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, fR[i]);
  float maxRapJet = 1 - fR[i];

  //==============================Reco jets===============================//
  fastjet::ClusterSequenceArea reco_cluster_seq(fullTracks, jet_def, area_def);
  vector<fastjet::PseudoJet> fjets_all =
      sorted_by_pt(reco_cluster_seq.inclusive_jets(fJetPtMin));

  fastjet::Selector fiducial_cut_selector = fastjet::SelectorAbsEtaMax(maxRapJet);
  vector<fastjet::PseudoJet> RecoJets = fiducial_cut_selector(fjets_all);
  vector<MyJet> myRecoJets;
  for (auto &rcJet : RecoJets)
    myRecoJets.push_back(MyJet(rcJet, rho));
  //======================================================================//

  if (mIsEmbedding) {
    //==============================MC jets===============================//
    fastjet::ClusterSequenceArea mc_cluster_seq(MCjetTracks, jet_def, area_def);
    vector<fastjet::PseudoJet> Mcjets_all =
        sorted_by_pt(mc_cluster_seq.inclusive_jets(1.));

    fastjet::Selector McFiducial_cut_selector =
        fastjet::SelectorAbsEtaMax(maxRapJet) * fastjet::SelectorPtMin(0.01) *
        fastjet::SelectorPtMax(1.5 * fpThatmax); // reject high-weight jets

    vector<fastjet::PseudoJet> McJets = McFiducial_cut_selector(Mcjets_all);
    vector<MyJet> myMcJets;
    for (auto &mcJet : McJets) {
      myMcJets.push_back(MyJet(mcJet, rho));
    }
    //======================================================================//

    vector<MatchedJetPair> MatchedJets;
    MatchedJets = MatchJetsEtaPhi(myMcJets, myRecoJets, fR[i]);

    for (unsigned int j = 0; j < MatchedJets.size(); j++) {
      fMcJet = MatchedJets[j].first;
      fRecoJet = MatchedJets[j].second;
      fDeltaR = fMcJet.deltaR(fRecoJet);
      fTree[i]->Fill();
    } // end loop over matched jets

  } else {
    //==================== Data mode: fill only reco jets ==================//
    for (const auto &recoJet : myRecoJets) {
      fRecoJet = recoJet;
      fMcJet = MyJet(); // empty
      fDeltaR = -1.0;
      fTree[i]->Fill();
    } // end loop over reco jets (data)
  } // end embedding/data condition
} // end loop over R

  for (int i = 0; i < 4800; i++) {
    Sump[i] = 0.0; // reset Sump array
  }
  return kStOK;
}



Double_t StPicoHFJetMaker::GetTowerCalibEnergy(Int_t TowerId) {
  StPicoBTowHit *tower =
      static_cast<StPicoBTowHit *>(mPicoDst->btowHit(TowerId - 1));
  Float_t pedestal, rms;
  Int_t status;
  mTables->getPedestal(BTOW, TowerId, 0, pedestal, rms);
  mTables->getStatus(BTOW, TowerId, status);
  Double_t *TowerCoeff;
  if (fRunNumber <= 15094020)
    TowerCoeff = CPre;
  else
    TowerCoeff = CLowMidHigh;
  Double_t calibEnergy = TowerCoeff[TowerId - 1] * (tower->adc() - pedestal);
  return calibEnergy;
}

//-----------------------------------------------------------------------------
////Correct tower eta for Vz position //// Not used anymore
//-----------------------------------------------------------------------------
Double_t StPicoHFJetMaker::vertexCorrectedEta(double eta, double vz) {
  double tower_theta = 2.0 * atan(exp(-eta));
  double z = 0.0;
  if (eta != 0.0)
    z = mBarrelRadius / tan(tower_theta);
  double z_diff = z - vz;
  double theta_corr = atan2(mBarrelRadius, z_diff);
  double eta_corr = -log(tan(theta_corr / 2.0));
  return eta_corr;
}

//-----------------------------------------------------------------------------
// Fill array with momentum of BEMC-matched tracks
//-----------------------------------------------------------------------------
Bool_t StPicoHFJetMaker::GetCaloTrackMomentum(StPicoDst *mPicoDst,
                                              TVector3 mPrimVtx) {
  // loop over global tracks  - towers
  UInt_t nTracks = mPicoDst->numberOfTracks();
  for (unsigned int itrack = 0; itrack < nTracks; itrack++) {
    StPicoTrack *trk = mPicoDst->track(itrack);
    TVector3 gMom = trk->gMom();
    // using global tracks
    double pT = gMom.Perp();
    if (pT != pT || pT < 0.2)
      continue;
    float eta = gMom.PseudoRapidity();
    if (fabs(eta) > 1)
      continue;
    float phi = gMom.Phi();

    float nHitsFit = trk->nHitsFit();
    float nHitsMax = trk->nHitsMax();
    if (nHitsFit < 15 || nHitsFit / nHitsMax < 0.52)
      continue; // some basic QA cuts
    double Bfield = mPicoDst->event()->bField();

    StPicoPhysicalHelix trkhelix = trk->helix(Bfield);
    float vtx_x = mPrimVtx.x();
    float vtx_y = mPrimVtx.y();
    float vtx_z = mPrimVtx.z();

    float dca_z = abs(trk->gDCAz(mPicoDst->event()->primaryVertex().z()));
    if (fabs(dca_z) > maxdcazhadroncorr)
      continue;
    int TowIndex = -99999;
    TowIndex = trk->bemcTowerIndex();
    float p = 0;
    if (TowIndex >= 0) {
      p = gMom.Mag();
      Sump[TowIndex] += p;
    }
  } // END global track loop
  return true;
}

//-----------------------------------------------------------------------------
// Jet matching function using only eta-phi criteria
//-----------------------------------------------------------------------------
vector<MatchedJetPair> MatchJetsEtaPhi(const vector<MyJet> &McJets,
                                       const vector<MyJet> &RecoJets,
                                       const double &R) {
  vector<MyJet> recoJetsCopy = RecoJets; // copy to avoid modifying the original
  vector<MatchedJetPair> matchedJets;
  MyJet dummy;

  for (const auto &mcJet : McJets) {
    bool isMatched = false;
    double minDeltaR = 10000; // Initialize to a large value
    auto bestMatch = recoJetsCopy.end();

    // Find the closest reco jet within the threshold
    for (auto rcit = recoJetsCopy.begin(); rcit != recoJetsCopy.end(); ++rcit) {
      MyJet recoJet = *rcit;
      double deltaR = mcJet.deltaR(recoJet);

      if (deltaR < 0.6 * R && deltaR < minDeltaR) {
        minDeltaR = deltaR;
        bestMatch = rcit;
        isMatched = true;
      }
    }

    // If a match was found, add it and remove from available jets
    if (isMatched) {
      matchedJets.push_back(make_pair(mcJet, *bestMatch));
      recoJetsCopy.erase(bestMatch);
    } else {
      // If no match was found for this MC jet, record it as unmatched
      matchedJets.push_back(make_pair(mcJet, dummy));
    }
  }

  // Add the remaining unmatched reco jets
  for (const auto &recoJet : recoJetsCopy)
    matchedJets.push_back(make_pair(dummy, recoJet));

  return matchedJets;
}
