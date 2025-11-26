#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TSystem.h"

#include <vector>
#include <string>
#include <iostream>

using std::string;
using std::vector;
using std::cout;
using std::endl;

// ------------------------- config -----------------------------

static const int kBayesIters[] = {1, 4, 5, 7};
static const int kNBayesIters = sizeof(kBayesIters)/sizeof(kBayesIters[0]);

static const double kTestFrac = 0.50;     // 50/50 split
static const uint32_t kSeed    = 12345;   // deterministic split

static const double kPtLeadCuts[] = {0.0, 5.0, 7.0};
static const int    kNPtLeadCuts  = sizeof(kPtLeadCuts)/sizeof(kPtLeadCuts[0]);

// measured & truth binning
static const int nbins_meas = 24;
static const double bin_meas_edges[nbins_meas+1] = {
  -100,-80,-60,-40,-20,-10,-5,-2.5,0,2.5,5,7.5,10,12.5,15,17.5,20,22.5,25,27.5,30,35,40,50,60
};
static const int nbins_truth = 10;
static const double bin_truth_edges[nbins_truth+1] = {
  0,5,10,15,20,25,30,35,40,50,60
};

static const vector<string> kCentralities = {"CENT_0_10", "MID_20_40", "PERI_60_80"};
static const vector<string> kRadii        = {"R0.2", "R0.3", "R0.4"};

// ------------------ jet-quality cuts -------------------------

// Area cuts per jet radius
static const double CUT_AREA_02 = 0.07;  // R = 0.2
static const double CUT_AREA_03 = 0.20;  // R = 0.3
static const double CUT_AREA_04 = 0.40;  // R = 0.4

// Max neutral energy fraction
static const double CUT_NEUTRAL_FRACTION = 0.95;

// --------------------------------------------------------------

static void EnsureDir(const string& path){
  if (gSystem->AccessPathName(path.c_str()))
    gSystem->mkdir(path.c_str(), /*recursive=*/true);
}

void unfold(const char* inputFile,
            const char* outDir)
{
  gStyle->SetOptStat(0);
  EnsureDir(outDir);

  TFile* fin = TFile::Open(inputFile, "READ");
  if (!fin || fin->IsZombie()) {
    cout << "[error] Cannot open input file: " << inputFile << endl;
    return;
  }

  // loop R, centrality, ptlead
  for (const auto& R : kRadii) {
    // parse numeric R to choose area cut
    double Rval = 0.0;
    if (sscanf(R.c_str(), "R%lf", &Rval) != 1) {
      cout << "[warn] Could not parse radius from " << R << ", skipping.\n";
      continue;
    }
    double areaMin = 0.0;
    if      (Rval < 0.25) areaMin = CUT_AREA_02; // ~0.2
    else if (Rval < 0.35) areaMin = CUT_AREA_03; // ~0.3
    else                  areaMin = CUT_AREA_04; // ~0.4

    for (const auto& C : kCentralities) {

      const string treePath = R + "/" + C + "/JetTree";
      TTree* tr = dynamic_cast<TTree*>(fin->Get(treePath.c_str()));
      if (!tr) {
        cout << "[note] missing tree: " << treePath << " (skip)\n";
        continue;
      }

      // speed up I/O: enable only needed branches
      tr->SetBranchStatus("*", 0);
      tr->SetBranchStatus("mc_pt", 1);
      tr->SetBranchStatus("mc_pt_lead", 1);
      tr->SetBranchStatus("reco_pt_corr", 1);
      tr->SetBranchStatus("reco_pt_lead", 1);
      tr->SetBranchStatus("reco_trigger_match", 1);
      tr->SetBranchStatus("centralityWeight", 1);
      tr->SetBranchStatus("xsecWeight", 1);
      tr->SetBranchStatus("reco_area", 1);
      tr->SetBranchStatus("reco_neutral_fraction", 1);

      // branch addresses
      float mc_pt=0, mc_pt_lead=0;
      float reco_pt_corr=0, reco_pt_lead=0;
      bool  reco_trigger_match=false;
      float centralityWeight=1.0f, xsecWeight=1.0f;
      float reco_area=0.0f, reco_neutral_fraction=0.0f;

      tr->SetBranchAddress("mc_pt", &mc_pt);
      tr->SetBranchAddress("mc_pt_lead", &mc_pt_lead);
      tr->SetBranchAddress("reco_pt_corr", &reco_pt_corr);
      tr->SetBranchAddress("reco_pt_lead", &reco_pt_lead);
      tr->SetBranchAddress("reco_trigger_match", &reco_trigger_match);
      tr->SetBranchAddress("centralityWeight", &centralityWeight);
      tr->SetBranchAddress("xsecWeight", &xsecWeight);
      tr->SetBranchAddress("reco_area", &reco_area);
      tr->SetBranchAddress("reco_neutral_fraction", &reco_neutral_fraction);

      const Long64_t n = tr->GetEntries();
      if (n <= 0) {
        cout << "[note] empty tree: " << treePath << " (skip)\n";
        continue;
      }

      for (int ic = 0; ic < kNPtLeadCuts; ++ic) {
        const double cut = kPtLeadCuts[ic];
        const string tag = R + "_" + C + Form("_ptlead%.0f", cut);

        TH1D* hMeasTrain = new TH1D(("hMeas_"+tag).c_str(),
            ";reco p_{T}^{corr} [GeV];dN/dp_{T}",
            nbins_meas, bin_meas_edges);
        TH1D* hTrueTrain = new TH1D(("hTrue_"+tag).c_str(),
            ";mc p_{T} [GeV];dN/dp_{T}",
            nbins_truth, bin_truth_edges);
        TH1D* hMeasTest  = (TH1D*)hMeasTrain->Clone(("hMeasTest_"+tag).c_str());
        TH1D* hTrueTest  = (TH1D*)hTrueTrain->Clone(("hTrueTest_"+tag).c_str());

        TH2D* hRespRecoVsTruth = new TH2D(("hResp_"+tag).c_str(),
            ";p_{T}^{reco,corr};p_{T}^{mc}",
            nbins_meas, bin_meas_edges,
            nbins_truth, bin_truth_edges);

        TRandom3 rng(kSeed);

        // event loop
        for (Long64_t i = 0; i < n; ++i) {
          if ((i % 200000) == 0) cout << "  ["<< tag <<"] " << i << "/" << n << "\r" << std::flush;
          tr->GetEntry(i);

          if (!reco_trigger_match) continue;
          if (reco_area < areaMin) continue;
          if (reco_neutral_fraction > CUT_NEUTRAL_FRACTION) continue;

          // dual ptlead cut (both reco & MC)
          if (!(reco_pt_lead >= cut && mc_pt_lead >= cut)) continue;

          const double w = (double)centralityWeight * (double)xsecWeight;

          const bool train = (rng.Uniform() > kTestFrac);
          if (train) {
            hRespRecoVsTruth->Fill(reco_pt_corr, mc_pt, w);
            hMeasTrain->Fill(reco_pt_corr, w);
            hTrueTrain->Fill(mc_pt, w);
          } else {
            hMeasTest->Fill(reco_pt_corr, w);
            hTrueTest->Fill(mc_pt, w);
          }
        }
        cout << endl;

        RooUnfoldResponse response(hMeasTrain, hTrueTrain, hRespRecoVsTruth);
        response.SetName(("response_"+tag).c_str());

        vector<TH1D*> unfolded(kNBayesIters, nullptr);
        for (int ib = 0; ib < kNBayesIters; ++ib) {
          RooUnfoldBayes u(&response, hMeasTest, kBayesIters[ib]);
          TH1D* hunf = (TH1D*)u.Hunfold();
          hunf->SetDirectory(nullptr);
          hunf->SetName(Form("Unfolded_%s_iter%d", tag.c_str(), kBayesIters[ib]));
          unfolded[ib] = hunf;
        }

        // plot
        TCanvas* c = new TCanvas(("c_"+tag).c_str(), "", 800, 1200);
        c->Divide(1,2);

        c->cd(1); gPad->SetLogy();
        TH1D* hM = (TH1D*)hMeasTest->Clone(("hMeasW_"+tag).c_str());
        TH1D* hT = (TH1D*)hTrueTest->Clone(("hTrueW_"+tag).c_str());
        hM->Scale(1.0, "width");
        hT->Scale(1.0, "width");
        hM->SetMarkerStyle(24); hM->SetLineColor(kBlue+1);
        hT->SetMarkerStyle(20); hT->SetLineColor(kBlack);
        hT->Draw("E1");
        hM->Draw("E1 SAME");

        TLegend* leg = new TLegend(0.58,0.58,0.88,0.88);
        leg->AddEntry(hT, "Truth (test)", "lp");
        leg->AddEntry(hM, "Measured (test)", "lp");

        for (int ib = 0; ib < kNBayesIters; ++ib) {
          TH1D* w = (TH1D*)unfolded[ib]->Clone(Form("UnfW_%d_%s", kBayesIters[ib], tag.c_str()));
          w->Scale(1.0, "width");
          w->SetMarkerStyle(20);
          w->SetMarkerColor(kMagenta + ib);
          w->SetLineColor(kMagenta + ib);
          w->Draw("E1 SAME");
          leg->AddEntry(w, Form("Bayes %d it.", kBayesIters[ib]), "lp");
        }
        leg->Draw();

        c->cd(2);
        TH1D* frame = (TH1D*)hT->Clone(("ratioFrame_"+tag).c_str());
        frame->Reset();
        frame->GetYaxis()->SetTitle("Unfolded / Truth");
        frame->GetYaxis()->SetRangeUser(0.5, 1.5);
        frame->Draw();
        for (int ib = 0; ib < kNBayesIters; ++ib) {
          TH1D* r = (TH1D*)unfolded[ib]->Clone(Form("ratio_%d_%s", kBayesIters[ib], tag.c_str()));
          r->Divide(hT);
          r->Draw(ib==0 ? "E1" : "E1 SAME");
        }
        TLine* ln = new TLine(frame->GetXaxis()->GetXmin(), 1.05,
                              frame->GetXaxis()->GetXmax(), 1.05);
        ln->SetLineStyle(2); ln->SetLineColor(kGray+1); ln->Draw();
        ln = new TLine(frame->GetXaxis()->GetXmin(), 0.95,
                       frame->GetXaxis()->GetXmax(), 0.95);
        ln->SetLineStyle(2); ln->SetLineColor(kGray+1); ln->Draw();

        const string tagfile = R + "_" + C + Form("_ptlead%.0f", cut);
        const string rootPath = string(outDir) + "/unfold_response_" + tagfile + ".root";
        const string pdfPath  = string(outDir) + "/closure_" + tagfile + ".pdf";

        TFile* outf = TFile::Open(rootPath.c_str(), "RECREATE");
        hRespRecoVsTruth->Write();
        response.Write();
        hMeasTrain->Write(); hTrueTrain->Write();
        hMeasTest->Write();  hTrueTest->Write();
        for (auto* u : unfolded) if (u) u->Write();
        outf->Close();

        c->SaveAs(pdfPath.c_str());

        // cleanup
        delete c;
        delete hM; delete hT;
        delete hRespRecoVsTruth;
        delete hMeasTrain; delete hTrueTrain;
        delete hMeasTest;  delete hTrueTest;
        for (auto* u: unfolded) delete u;
      } // ptlead cuts
    } // centralities
  } // radii

  cout << "Done. Outputs in: " << outDir << endl;
  fin->Close();
  delete fin;
}
