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
#include <algorithm>

#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/Selector.hh"
#include "fastjet/config.h"
#include "fastjet/tools/JetMedianBackgroundEstimator.hh"

#include "MyJet.h"

using namespace std;

const char* kCentTag[4] = {
  "",            // index 0 unused
  "CENT_0_10",   // c3 = 1
  "MID_20_40",   // c3 = 2
  "PERI_60_80"   // c3 = 3
};


const double CUT_AREA_02 = 0.07; // R = 0.2
const double CUT_AREA_03 = 0.20; // R = 0.3
const double CUT_AREA_04 = 0.40; // R = 0.4

const double CUT_NEUTRAL_FRACTION = 0.95;

vector<MatchedJetPair> MatchJetsEtaPhi(const vector<MyJet> &McJets,
                                       const vector<MyJet> &RecoJets,
                                       const double &R);

inline bool passHistoCuts(const MyJet &j, double R);

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
  mADCtoEMaker = dynamic_cast<StEmcADCtoEMaker*>(GetMaker("Eread"));
  assert(mADCtoEMaker);
  mTables = mADCtoEMaker->getBemcData()->getTables();

  TH1::SetDefaultSumw2();

  mOutList->SetName("QA_histograms"); 

  mOutList->Add(new TH1D("hcent9", "centrality9;bin;events", 10, -1, 9));

  TDirectory* fileDir = gDirectory;
  fTreeRC.clear();
  fTreeRC.reserve(fR.size());

  for (size_t iR = 0; iR < fR.size(); ++iR) {
    const TString rName = Form("R%.1f", fR[iR]);
    TDirectory* rdir = fileDir->mkdir(rName);
    if (!rdir) rdir = (TDirectory*)fileDir->Get(rName);
    rdir->cd();

    std::vector<TTree*> treesC;  // 3 classes: 1..3 (we'll index 0..2)
    treesC.reserve(3);

    // book three classes (central, midcentral, peripheral)
    for (int c3 = 1; c3 <= 3; ++c3) {
      TDirectory* cdir = rdir->mkdir(kCentTag[c3]);
      if (!cdir) cdir = (TDirectory*)rdir->Get(kCentTag[c3]);
      cdir->cd();
      const int ci = c3; 

      if (iR == 0 && c3 == 1) {
        const size_t nR = fR.size();
        fH2_den.assign(nR, std::vector<TH2D*>(4, nullptr));
        fH2_num.assign(nR, std::vector<TH2D*>(4, nullptr));
        fH1_reco.assign(nR, std::vector<TH1D*>(4, nullptr));
        fH1_mc.assign(nR, std::vector<TH1D*>(4, nullptr));
        fH2_reco_mc.assign(nR, std::vector<TH2D*>(4, nullptr));
        fH2_reco_matched.assign(nR, std::vector<TH2D*>(4, nullptr));

      }

      
      const int nb_pt   = 1000;  const double pt_min   = -40.0, pt_max   = 60.0;
      const int nb_lead = 200;  const double lead_min = 0.0, lead_max = 30.0;

    fH2_den[iR][ci] = new TH2D("den_ptcorr_vs_ptlead",
      Form("Den: p_{T}^{corr} vs p_{T}^{lead} (R=%.1f, %s);p_{T}^{corr} [GeV];p_{T}^{lead} [GeV]",
       fR[iR], kCentTag[c3]), nb_pt, pt_min, pt_max, nb_lead, lead_min, lead_max);
    fH2_den[iR][ci]->SetDirectory(cdir);

    fH2_num[iR][ci] = new TH2D("num_ptcorr_vs_ptlead",
      Form("Num: (trg) p_{T}^{corr} vs p_{T}^{lead} (R=%.1f, %s);p_{T}^{corr} [GeV];p_{T}^{lead} [GeV]",
       fR[iR], kCentTag[c3]), nb_pt, pt_min, pt_max, nb_lead, lead_min, lead_max);
    fH2_num[iR][ci]->SetDirectory(cdir);

    fH1_reco[iR][ci] = new TH1D("reco_ptcorr",
      Form("Reco jet p_{T}^{corr} (R=%.1f, %s);p_{T}^{corr} [GeV];Jets",
       fR[iR], kCentTag[c3]), nb_pt, pt_min, pt_max);
    fH1_reco[iR][ci]->SetDirectory(cdir);


    if (mIsEmbedding) {
      fH1_mc[iR][ci] = new TH1D("mc_pt",
        Form("MC jet p_{T} (R=%.1f, %s);p_{T}^{MC} [GeV];Jets",
         fR[iR], kCentTag[c3]), nb_pt, pt_min, pt_max);
      fH1_mc[iR][ci]->SetDirectory(cdir);

      fH2_reco_mc[iR][ci] = new TH2D("recoptcorr_vs_mcpt",
        Form("Reco p_{T}^{corr} vs MC p_{T} (R=%.1f, %s);p_{T}^{MC} [GeV];p_{T}^{corr} [GeV]",
         fR[iR], kCentTag[c3]), nb_pt, pt_min, pt_max, nb_pt, pt_min, pt_max);
      fH2_reco_mc[iR][ci]->SetDirectory(cdir);

      fH2_reco_matched[iR][ci] = new TH2D("reco_matched_ptcorr_vs_ptlead",
        Form("Matched reco jets p_{T}^{corr} vs p_{T}^{lead} (R=%.1f, %s);p_{T}^{corr} [GeV];p_{T}^{lead} [GeV]",
         fR[iR], kCentTag[c3]), nb_pt, pt_min, pt_max, nb_lead, lead_min, lead_max);
      fH2_reco_matched[iR][ci]->SetDirectory(cdir);
    }
      
      // ---- TTree per (R,class); NO centrality branch
      TTree* jetTree = new TTree("JetTree", "JetTree");
      jetTree->Branch("runId", &fRunNumber, "runId/I");
      jetTree->Branch("centralityWeight", &fCentralityWeight, "centralityWeight/F"); // keep weight
      if (mIsEmbedding) {
        jetTree->Branch("xsecWeight", &fXsecWeight, "xsecWeight/F");
        jetTree->Branch("deltaR", &fDeltaR, "deltaR/F");
        jetTree->Branch("mc_pt", &fMcJet.pt, "mc_pt/F");
        jetTree->Branch("mc_eta", &fMcJet.eta, "mc_eta/F");
        jetTree->Branch("mc_phi", &fMcJet.phi, "mc_phi/F");
        jetTree->Branch("mc_area", &fMcJet.area, "mc_area/F");
        jetTree->Branch("mc_pt_lead", &fMcJet.pt_lead, "mc_pt_lead/F");
        jetTree->Branch("mc_n_constituents", &fMcJet.n_constituents, "mc_n_constituents/I");
        jetTree->Branch("mc_neutral_fraction", &fMcJet.neutral_fraction, "mc_neutral_fraction/F");
      }
      jetTree->Branch("reco_pt", &fRecoJet.pt, "reco_pt/F");
      jetTree->Branch("reco_pt_corr", &fRecoJet.pt_corr, "reco_pt_corr/F");
      jetTree->Branch("reco_eta", &fRecoJet.eta, "reco_eta/F");
      jetTree->Branch("reco_phi", &fRecoJet.phi, "reco_phi/F");
      jetTree->Branch("reco_area", &fRecoJet.area, "reco_area/F");
      jetTree->Branch("reco_rho", &fRecoJet.rho, "reco_rho/F");
      jetTree->Branch("reco_pt_lead", &fRecoJet.pt_lead, "reco_pt_lead/F");
      jetTree->Branch("reco_n_constituents", &fRecoJet.n_constituents, "reco_n_constituents/I");
      jetTree->Branch("reco_neutral_fraction", &fRecoJet.neutral_fraction, "reco_neutral_fraction/F");
      jetTree->Branch("reco_trigger_match", &fRecoJet.trigger_match, "reco_trigger_match/O");

      treesC.push_back(jetTree);

      rdir->cd(); // back up for next class
    }

    fTreeRC.push_back(treesC);
    fileDir->cd();
  }

  return kStOK;
}

// _________________________________________________________
void StPicoHFJetMaker::ClearJets(Option_t *opt = "") { return; }

// _________________________________________________________
int StPicoHFJetMaker::FinishJets() {
  TDirectory* fileDir = gDirectory;
  if (!fileDir) return kStOK;

  const size_t nR = fR.size();

  for (size_t iR = 0; iR < nR; ++iR) {
    TDirectory* rdir = dynamic_cast<TDirectory*>(fileDir->Get(Form("R%.1f", fR[iR])));
    if (!rdir) continue;

    for (int c3 = 1; c3 <= 3; ++c3) {
      const int ciTree = c3 - 1;  // 0..2

      TDirectory* cdir = dynamic_cast<TDirectory*>(rdir->Get(kCentTag[c3]));
      if (!cdir) continue;
      cdir->cd();

      // --- write the tree
       if (iR < fTreeRC.size() && ciTree >= 0 && ciTree < (int)fTreeRC[iR].size() && fTreeRC[iR][ciTree]) {
        fTreeRC[iR][ciTree]->Write();
      }

      // --- write the histograms (if they exist)
      const size_t ci = size_t(c3);
      if (iR < fH2_den.size()) {
        if (ci < fH2_den[iR].size()     && fH2_den[iR][ci])     fH2_den[iR][ci]->Write();
        if (ci < fH2_num[iR].size()     && fH2_num[iR][ci])     fH2_num[iR][ci]->Write();
        if (ci < fH1_reco[iR].size()    && fH1_reco[iR][ci])    fH1_reco[iR][ci]->Write();
      if (mIsEmbedding) {
        if (ci < fH1_mc[iR].size()          && fH1_mc[iR][ci])          fH1_mc[iR][ci]->Write();
        if (ci < fH2_reco_mc[iR].size()     && fH2_reco_mc[iR][ci])     fH2_reco_mc[iR][ci]->Write();
        if (ci < fH2_reco_matched[iR].size()&& fH2_reco_matched[iR][ci])fH2_reco_matched[iR][ci]->Write(); // NEW
      }
      }
    } // c3
  }   // iR

  fileDir->cd();
  return kStOK;
}


// _________________________________________________________
int StPicoHFJetMaker::MakeJets() {
  
  TH1D *hcent9 = static_cast<TH1D *>(mOutList->FindObject("hcent9"));
  TAxis* ax = hcent9->GetXaxis();
  const char* lab9[10] = {
    "undef (-1)",
    "0: 0-5%",
    "1: 5-10%",
    "2: 10-20%",
    "3: 20-30%",
    "4: 30-40%",
    "5: 40-50%",
    "6: 50-60%",
    "7: 60-70%",
    "8: 70-80%"
  };
  for (int b = 1; b <= 10; ++b) ax->SetBinLabel(b, lab9[b-1]);
  ax->CenterLabels(true);

  vector<fastjet::PseudoJet> jetTracks;
  vector<fastjet::PseudoJet> neutraljetTracks; // from bemc towers only
  vector<fastjet::PseudoJet> fullTracks;
  vector<fastjet::PseudoJet> MCjetTracks;



  fRunNumber = mPicoDst->event()->runId();
  int eventId = mPicoDst->event()->eventId(); // eventID
  (void)eventId;
  int refMult = mPicoDst->event()->refMult();
  double vz = mPrimVtx.z();
  mRefmultCorrUtil->setEvent(fRunNumber, refMult, mPicoDst->event()->ZDCx(),
                             vz);
  fCentrality = mRefmultCorrUtil->centrality9(); // 0 = 0-5 %,..., 8 = 70-80
                                                 // %
  if (fCentrality == -1)
    return kStOK; // no fCentrality

  fCentralityWeight = mRefmultCorrUtil->weight();

  if (hcent9) hcent9->Fill(fCentrality, fCentralityWeight);

  // map to 3 classes: 1=0-10%, 2=20-40%, 3=60-80%, else 0 (=skip)
  int c3 = 0;
  if (fCentrality == 0 || fCentrality == 1)       c3 = 1; // 0-10%
  else if (fCentrality == 3 || fCentrality == 4)  c3 = 2; // 20-40%
  else if (fCentrality == 7 || fCentrality == 8)  c3 = 3; // 60-80%
  if (c3 == 0) {
  // Reset Sump[] before returning
  for (int i = 0; i < 4800; i++) Sump[i] = 0.0;
  return kStOK;
  }

  const double w_event = fCentralityWeight * (mIsEmbedding ? fXsecWeight : 1.0);
  const int ci = c3;

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

if (mIsEmbedding && fpThatmax > 0.0 && !MCjetTracks.empty()) {
  float Rcheck = fR.empty() ? 0.4f : *std::max_element(fR.begin(), fR.end());
  fastjet::JetDefinition mc_jet_def_veto(fastjet::antikt_algorithm, Rcheck);
  fastjet::ClusterSequence mc_cs_veto(MCjetTracks, mc_jet_def_veto);
  std::vector<fastjet::PseudoJet> mcjets_veto =
      sorted_by_pt(mc_cs_veto.inclusive_jets(1.0)); // pT > 1 GeV

  const double ptMaxVeto = 1.5 * fpThatmax;

  bool vetoEvent = false;
  for (const auto &j : mcjets_veto) {
    if (j.perp() > ptMaxVeto) {
      vetoEvent = true;
      break;
    }
  }

  if (vetoEvent) {
    for (int i = 0; i < 4800; i++) Sump[i] = 0.0;
    return kStOK;
  }
}
  // RC part
  GetCaloTrackMomentum(mPicoDst, mPrimVtx); // fill array Sump with momenta of
                                            // tracks which are matched to BEMC

  StEmcPosition *mEmcPosition = new StEmcPosition();

//  double TOWE = 0;
  for (int iTow = 0; iTow < 4800; iTow++) { // get btow info
    StPicoBTowHit *towHit = mPicoDst->btowHit(iTow);
    if (!towHit || towHit->isBad())
      continue; // if the tower is marked as bad or missing info
    int realtowID = towHit->numericIndex2SoftId(iTow);
    if (BadTowerMap[realtowID])
      continue; // exclude bad towers (map in JetInfo.h)


    double towE = GetTowerCalibEnergy(iTow + 1); // get tower energy
//    TOWE = towE; // just keep track of the original energy for trigger approximation

    if (doTowErrPlus == true) {towE = towE + 0.038 * towE;}
    if (doTowErrMinus == true) {towE = towE - 0.038 * towE;}

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
    if (Towphi < 0) Towphi += 2.0*TMath::Pi();
    if (Towphi >= 2.0*TMath::Pi()) Towphi -= 2.0*TMath::Pi();

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
      if (ADC > fTrgthresh){
     // triggerTowersEtaPhi.emplace_back(Toweta, Towphi);
      inputTower.set_user_index(9999); // mark trigger towers with user_index 9999
      }
      neutraljetTracks.push_back(inputTower);
    }
  } // end get btow info

  delete mEmcPosition;

  // loop over primary tracks
  for (unsigned int i = 0; i < mIdxPicoParticles.size(); i++) {
    StPicoTrack *trk = mPicoDst->track(mIdxPicoParticles[i]);
   if (doTrackErr) {
    static TRandom3 randGen;
    if (randGen.Rndm() > 0.96) continue;
  }

    const TVector3 p = trk->pMom();
    const double pT  = p.Perp();
    if (!(pT > 0)) continue;                // NaN/zero guard
    const float eta = p.PseudoRapidity();
    if (fabs(eta) > 1.0) continue;          // your fiducial cut
    float phi = trk->pMom().Phi();
    float dca = (mPrimVtx - trk->origin()).Mag();
    float charged = trk->charge();


  (void)phi; (void)dca; (void)charged;
  
  fastjet::PseudoJet pj(p.x(), p.y(), p.z(), p.Mag());

      if (mIsEmbedding) {
       if (trk->qaTruth() > 95) pj.set_user_index(trk->idTruth() - 1);
       else                     pj.set_user_index(trk->charge() ? 1 : 0);
      } else {
      pj.set_user_index(trk->charge() ? 1 : 0);
      }

      jetTracks.push_back(pj);  // <-- add this line

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
nJetsRemove = (c3 == 1 ? 2 : 1); // remove 2 hardest jets in central, 1 otherwise

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
std::vector<fastjet::PseudoJet> fjets_all =
    sorted_by_pt(reco_cluster_seq.inclusive_jets(fJetPtMin));

fastjet::Selector fiducial_cut_selector = fastjet::SelectorAbsEtaMax(maxRapJet);
std::vector<fastjet::PseudoJet> RecoJets = fiducial_cut_selector(fjets_all);

std::vector<MyJet> myRecoJets;
myRecoJets.reserve(RecoJets.size());
for (auto &rcJet : RecoJets) {
  myRecoJets.push_back(MyJet(rcJet, rho));
}
  //======================================================================//

  // pick tree for this (R, class)
  TTree* jetTree = nullptr;
  const int ciTree = c3 - 1;  // 0..2
  if (i < fTreeRC.size() && ciTree >= 0 && ciTree < (int)fTreeRC[i].size())
    jetTree = fTreeRC[i][ciTree];

  // histogram pointers for this (R, class)
  TH2D* hDen    = (i < fH2_den.size()     && ci >= 0 && fH2_den[i][ci])     ? fH2_den[i][ci]     : nullptr;
  TH2D* hNum    = (i < fH2_num.size()     && ci >= 0 && fH2_num[i][ci])     ? fH2_num[i][ci]     : nullptr;
  TH1D* hReco   = (i < fH1_reco.size()    && ci >= 0 && fH1_reco[i][ci])    ? fH1_reco[i][ci]    : nullptr;
  TH1D* hMc     = (mIsEmbedding && i < fH1_mc.size() && ci >= 0 && fH1_mc[i][ci]) ? fH1_mc[i][ci] : nullptr;
  TH2D* hRecoMc = (mIsEmbedding && i < fH2_reco_mc.size() && ci >= 0 && fH2_reco_mc[i][ci]) ? fH2_reco_mc[i][ci] : nullptr;
  TH2D* hRecoMatched = (mIsEmbedding && i < fH2_reco_matched.size() && ci >= 0 && fH2_reco_matched[i][ci]) ? fH2_reco_matched[i][ci] : nullptr;

//==================== Embedding mode ==================//
if (mIsEmbedding) {
  //============================== MC jets ===============================//
  fastjet::ClusterSequenceArea mc_cluster_seq(MCjetTracks, jet_def, area_def);
  vector<fastjet::PseudoJet> Mcjets_all =
      sorted_by_pt(mc_cluster_seq.inclusive_jets(1.0));

  fastjet::Selector McFiducial_cut_selector =
      fastjet::SelectorAbsEtaMax(maxRapJet) *
      fastjet::SelectorPtMin(0.01);

vector<fastjet::PseudoJet> McJets = McFiducial_cut_selector(Mcjets_all);
vector<MyJet> myMcJets;
myMcJets.reserve(McJets.size());
for (auto &mcJet : McJets) {
  // truth jet: no background subtraction, rho stored as 0
  myMcJets.push_back(MyJet(mcJet, 0.0f));
}


  //========================= MC–Reco matching ===========================//
  vector<MatchedJetPair> MatchedJets = MatchJetsEtaPhi(myMcJets, myRecoJets, fR[i]);

  for (const auto& mp : MatchedJets) {
    fMcJet   = mp.first;
    fRecoJet = mp.second;
    fDeltaR  = fMcJet.deltaR(fRecoJet);

    const bool haveReco = (fRecoJet.pt >= 0);
    const bool haveMC   = (fMcJet.pt   >= 0);

    // --- MC-only histograms (all MC jets) ---
    if (haveMC && hMc) {
      hMc->Fill(fMcJet.pt, w_event);
    }

    // ================= Reco-part: trigger-efficiency histos =================
    if (haveReco) {
      // Always apply area & neutral-fraction cuts for trigger-efficiency histos
      const bool passCuts = passHistoCuts(fRecoJet, fR[i]);

      // Den: all reco jets that pass cuts
      if (passCuts && hDen) {
        hDen->Fill(fRecoJet.pt_corr, fRecoJet.pt_lead, w_event);
      }

      // Reco-only, non-trigger jets: only Den (if above) and nothing else,
      if (!fRecoJet.trigger_match && !haveMC) {
        continue; // skip Num/Reco/Response/Tree for these
      }

      // Num: triggered jets only, with cuts
      if (fRecoJet.trigger_match && passCuts && hNum) {
        hNum->Fill(fRecoJet.pt_corr, fRecoJet.pt_lead, w_event);
      }

      // Reco spectrum (for jets in the trigger/MC sample), with cuts
      if (passCuts && hReco) {
        hReco->Fill(fRecoJet.pt_corr, w_event);
      }

      // Response matrix: only matched jets, with cuts
      if (haveMC && passCuts && hRecoMc) {
        hRecoMc->Fill(fMcJet.pt, fRecoJet.pt_corr, w_event);
      }

      // Matched reco jets for trigger/MC sample, with cuts
      if (haveMC && passCuts && hRecoMatched) {
        hRecoMatched->Fill(fRecoJet.pt_corr, fRecoJet.pt_lead, w_event);
      }

    } // end if (haveReco)

    // ================= Tree filling =================
    // Now we want:
    //  - all MC jets (even if no reco) in the tree
    //  - plus reco jets with MC or trigger (same as before) 
    if (jetTree) {
      // MC-only jets: haveMC == true, haveReco == false -> included here
      // Reco-only triggered jets: haveMC == false, trigger_match == true -> also included
      // Matched jets: haveMC == true, haveReco == true -> included
      if (haveMC || (haveReco && fRecoJet.trigger_match)) {
        jetTree->Fill();
      }
    }
  } // end loop over MatchedJets

} else {
  //==================== Data mode ===========================//
  for (const auto& rj : myRecoJets) {
    fRecoJet = rj;
    fMcJet   = MyJet();   // dummy
    fDeltaR  = -1.0;

    // Always apply area & neutral-fraction cuts for trigger-efficiency histos
    const bool passCuts = passHistoCuts(fRecoJet, fR[i]);

    // Den: all reco jets that pass cuts
    if (passCuts && hDen) {
      hDen->Fill(fRecoJet.pt_corr, fRecoJet.pt_lead, w_event);
    }

    // Num: triggered jets only, with cuts
    if (fRecoJet.trigger_match && passCuts && hNum) {
      hNum->Fill(fRecoJet.pt_corr, fRecoJet.pt_lead, w_event);
    }

    // Reco spectrum of triggered jets, with cuts
    if (fRecoJet.trigger_match && passCuts && hReco) {
      hReco->Fill(fRecoJet.pt_corr, w_event);
    }

    // Tree: only triggered jets, independent of cuts
    if (!fRecoJet.trigger_match) continue;
    if (jetTree) jetTree->Fill();
  } // end loop over reco jets (data)
} // end embedding/data
} // end loop over R

  for (int i = 0; i < 4800; i++) {
    Sump[i] = 0.0; // reset Sump array
  }
  return kStOK;
}


//-----------------------------------------------------------------------------
////Correct tower energy
//-----------------------------------------------------------------------------
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
  //  float phi = gMom.Phi();

    float nHitsFit = trk->nHitsFit();
    float nHitsMax = trk->nHitsMax();
    if (nHitsFit < 15 || nHitsFit / nHitsMax < 0.52)
      continue; // some basic QA cuts
    double Bfield = mPicoDst->event()->bField();

    StPicoPhysicalHelix trkhelix = trk->helix(Bfield);
  //  float vtx_x = mPrimVtx.x();
  //  float vtx_y = mPrimVtx.y();
  //  float vtx_z = mPrimVtx.z();

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
  const double matchRadius = 0.6*R; 
  vector<char> recoUsed(RecoJets.size(), 0);
  vector<MatchedJetPair> matchedJets;
  matchedJets.reserve(McJets.size() + RecoJets.size());

  // For each MC jet, find the closest unused reco jet within matchRadius
  for (const auto &mcJet : McJets) {
    int bestIdx = -1;
    double bestDr = matchRadius;
    for (size_t j = 0; j < RecoJets.size(); ++j) {
      if (recoUsed[j]) continue;
      double dr = mcJet.deltaR(RecoJets[j]); 
      if (dr < bestDr) {
        bestDr = dr;
        bestIdx = (int)j;
      }
    }
    if (bestIdx >= 0) {
      matchedJets.emplace_back(mcJet, RecoJets[bestIdx]);
      recoUsed[bestIdx] = 1;
    } else {
      matchedJets.emplace_back(mcJet, MyJet()); // unmatched MC -> dummy reco
    }
  }

  // Append unmatched reco jets
  for (size_t j = 0; j < RecoJets.size(); ++j) {
    if (!recoUsed[j]) matchedJets.emplace_back(MyJet(), RecoJets[j]);
  }

  return matchedJets;
}

//-----------------------------------------------------------------------------
// Apply histogram cuts
//-----------------------------------------------------------------------------
inline bool passHistoCuts(const MyJet &j, double R) {
    // area cut per radius
    double acut = (R < 0.25 ? CUT_AREA_02 :
                  (R < 0.35 ? CUT_AREA_03 :
                              CUT_AREA_04));

    if (j.area < acut) return false;
    if (j.neutral_fraction > CUT_NEUTRAL_FRACTION) return false;

    return true;
}