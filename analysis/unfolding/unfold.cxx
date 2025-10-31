// ROOT macro for a Bayes closure check with RooUnfold — same logic, clearer
// names. Usage (ROOT): .L unfold_simple.C++ ;
// unfold_simple(".../merged_matching_JP2_R0.2.root")

#include "RooUnfoldBayes.h"
#include "RooUnfoldResponse.h"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TRandom3.h"
#include "TStyle.h"
#include "TTree.h"

#include <iostream>
#include <vector>

using std::cout;
using std::endl;

// Bin edges
static const std::vector<double> recoPtBins = {
    5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
    22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,
    39, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 68, 80, 90};
static const std::vector<double> truthPtBins = {4.0,  6.9,  8.2,  9.7,  11.5,
                                                13.6, 16.1, 19.0, 22.5, 26.6,
                                                31.4, 37.2, 44.0, 52.0, 80.0};

// Filename parsers
static TString inferTriggerFromPath(const TString &path) {
  return path.Contains("JP2") ? "JP2" : "HT2";
}
static double inferJetRadiusFromPath(const TString &path) {
  if (path.Contains("R0.2"))
    return 0.2;
  if (path.Contains("R0.3"))
    return 0.3;
  if (path.Contains("R0.4"))
    return 0.4;
  if (path.Contains("R0.5"))
    return 0.5;
  if (path.Contains("R0.6"))
    return 0.6;
  return -1.0;
}

void unfold_simple(TString inputFile = "/home/prozorov/dev/star/jets_pp_2012/"
                                       "output/merged_matching_HT2_R0.2.root") {
  gStyle->SetOptStat(0);

  const TString triggerTag = inferTriggerFromPath(inputFile);
  const double jetRadius = inferJetRadiusFromPath(inputFile);
  if (jetRadius < 0) {
    cout << "Radius not found in filename. Abort.\n";
    return;
  }

  // Input
  TFile *input = TFile::Open(inputFile, "READ");
  if (!input || input->IsZombie()) {
    cout << "Cannot open input.\n";
    return;
  }

  TTree *matchedTree = (TTree *)input->Get("MatchedTree");
  if (!matchedTree) {
    cout << "Tree 'MatchedTree' missing.\n";
    return;
  }

  // Branches
  double mcJetPt = 0.0, recoJetPt = 0.0, eventWeight = 1.0;
  bool recoMatchedHT2 = false; // reserved for potential selection
  matchedTree->SetBranchAddress("mc_pt", &mcJetPt);
  matchedTree->SetBranchAddress("reco_pt", &recoJetPt);
  matchedTree->SetBranchAddress("mc_weight", &eventWeight);
  matchedTree->SetBranchAddress("reco_trigger_match_HT2", &recoMatchedHT2);

  // Histograms (train/test split)
  TH1D *recoSpectrumTrain =
      new TH1D("Measured", ";reco p_{T} [GeV/c];dN/dp_{T}",
               (int)recoPtBins.size() - 1, recoPtBins.data());
  TH1D *truthSpectrumTrain =
      new TH1D("Truth", ";mc p_{T} [GeV/c];dN/dp_{T}",
               (int)truthPtBins.size() - 1, truthPtBins.data());
  TH1D *recoSpectrumTest = (TH1D *)recoSpectrumTrain->Clone("MeasuredTest");
  TH1D *truthSpectrumTest = (TH1D *)truthSpectrumTrain->Clone("TruthTest");

  TH2D *responseMatrixRecoVsTruth = new TH2D(
      "hResponseMatrix", ";p_{T}^{reco};p_{T}^{mc}", (int)recoPtBins.size() - 1,
      recoPtBins.data(), (int)truthPtBins.size() - 1, truthPtBins.data());

  // Train/test split
  TRandom3 splitRng(0);
  const double testFraction = 0.20;

  const Long64_t nEntries = matchedTree->GetEntries();
  for (Long64_t entry = 0; entry < nEntries; ++entry) {
    if ((entry % 100000) == 0)
      cout << "Entry " << entry << "/" << nEntries << "\r" << std::flush;

    matchedTree->GetEntry(entry);
    if (mcJetPt <= 0 || recoJetPt <= 0)
      continue;

    const bool useForTraining = (splitRng.Uniform() > testFraction);
    if (useForTraining) {
      responseMatrixRecoVsTruth->Fill(recoJetPt, mcJetPt, eventWeight);
      recoSpectrumTrain->Fill(recoJetPt, eventWeight);
      truthSpectrumTrain->Fill(mcJetPt, eventWeight);
    } else {
      recoSpectrumTest->Fill(recoJetPt, eventWeight);
      truthSpectrumTest->Fill(mcJetPt, eventWeight);
    }
  }

  // Response
  RooUnfoldResponse bayesResponse(recoSpectrumTrain, truthSpectrumTrain,
                                  responseMatrixRecoVsTruth);
  bayesResponse.SetName("response");

  // Unfold test sample
  const int bayesIterations[] = {1, 2, 3, 4};
  const int nBayesIterations =
      sizeof(bayesIterations) / sizeof(bayesIterations[0]);

  // Canvas: spectra + ratio
  TCanvas *canvas = new TCanvas("c", "", 800, 1200);
  canvas->Divide(1, 2);

  // Top panel: spectra (bin-width normalized)
  canvas->cd(1);
  gPad->SetLogy();

  TH1D *recoTestWidthNorm =
      (TH1D *)recoSpectrumTest->Clone("recoTestWidthNorm");
  TH1D *truthTestWidthNorm =
      (TH1D *)truthSpectrumTest->Clone("truthTestWidthNorm");
  recoTestWidthNorm->Scale(1.0, "width");
  truthTestWidthNorm->Scale(1.0, "width");

  recoTestWidthNorm->SetMarkerStyle(24);
  recoTestWidthNorm->SetLineColor(kBlue + 1);
  truthTestWidthNorm->SetMarkerStyle(20);
  truthTestWidthNorm->SetLineColor(kBlack);

  truthTestWidthNorm->Draw("E1");
  recoTestWidthNorm->Draw("E1 SAME");

  TLegend *legend = new TLegend(0.60, 0.60, 0.88, 0.88);
  legend->AddEntry(truthTestWidthNorm, "Truth (test)", "lp");
  legend->AddEntry(recoTestWidthNorm, "Measured (test)", "lp");

  // Unfold and overlay
  std::vector<TH1D *> unfoldedSpectra(nBayesIterations, nullptr);
  for (int i = 0; i < nBayesIterations; ++i) {
    RooUnfoldBayes unfoldBayes(&bayesResponse, recoSpectrumTest,
                               bayesIterations[i]);
    TH1D *unfoldedWidthNorm = (TH1D *)unfoldBayes.Hunfold();
    unfoldedWidthNorm->SetDirectory(0);
    unfoldedWidthNorm->SetName(Form("Unfolded_iter%d", bayesIterations[i]));
    unfoldedWidthNorm->Scale(1.0, "width");
    unfoldedWidthNorm->SetMarkerStyle(20);
    unfoldedWidthNorm->SetMarkerColor(kMagenta + i);
    unfoldedWidthNorm->SetLineColor(kMagenta + i);
    unfoldedWidthNorm->Draw("E1 SAME");
    legend->AddEntry(unfoldedWidthNorm,
                     Form("Bayes %d it.", bayesIterations[i]), "lp");
    unfoldedSpectra[i] = unfoldedWidthNorm;
  }
  legend->Draw();

  // Bottom panel: Unfolded/Truth ratio
  canvas->cd(2);
  TH1D *ratioFrame = (TH1D *)truthTestWidthNorm->Clone("ratioFrame");
  ratioFrame->Reset();
  ratioFrame->GetYaxis()->SetTitle("Unfolded / Truth");
  ratioFrame->Draw("");
  for (int i = 0; i < nBayesIterations; ++i) {
    TH1D *ratioHist =
        (TH1D *)unfoldedSpectra[i]->Clone(Form("ratio_%d", bayesIterations[i]));
    ratioHist->Divide(truthTestWidthNorm);
    ratioHist->GetYaxis()->SetRangeUser(0.8, 1.2);
    ratioHist->Draw(i == 0 ? "E1" : "E1 SAME");
  }
  // ±5% guide lines
  TLine *bandLine = new TLine(ratioFrame->GetXaxis()->GetXmin(), 1.05,
                              ratioFrame->GetXaxis()->GetXmax(), 1.05);
  bandLine->SetLineStyle(2);
  bandLine->SetLineColor(kGray + 1);
  bandLine->Draw();
  bandLine = new TLine(ratioFrame->GetXaxis()->GetXmin(), 0.95,
                       ratioFrame->GetXaxis()->GetXmax(), 0.95);
  bandLine->SetLineStyle(2);
  bandLine->SetLineColor(kGray + 1);
  bandLine->Draw();

  // Outputs
  TString tag = Form("%s_R%.1f", triggerTag.Data(), jetRadius);
  TFile *out =
      TFile::Open(Form("simple_response_%s.root", tag.Data()), "RECREATE");
  responseMatrixRecoVsTruth->Write();
  bayesResponse.Write();
  recoSpectrumTrain->Write();
  truthSpectrumTrain->Write();
  recoSpectrumTest->Write();
  truthSpectrumTest->Write();
  for (auto *h : unfoldedSpectra)
    if (h)
      h->Write();
  out->Close();

  canvas->SaveAs(Form("simple_closure_%s.pdf", tag.Data()));
}