#include "TFile.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TKey.h"
#include "TCollection.h"
#include "TROOT.h"

#include <iostream>
#include <string>

// =========================================================
//  Jet-quality cuts
// =========================================================

// Area cuts per jet radius
const double CUT_AREA_02 = 0.07;  // R = 0.2
const double CUT_AREA_03 = 0.20;  // R = 0.3
const double CUT_AREA_04 = 0.40;  // R = 0.4

// Max neutral energy fraction
const double CUT_NEUTRAL_FRACTION = 0.95;  
// =========================================================

void make_hists(const char *infile  = "embedding_merged.root",
                const char *outfile = "weighted_hists_matched.root")
{
    TFile *fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "Error: cannot open input file " << infile << std::endl;
        return;
    }

    TFile *fout = TFile::Open(outfile, "RECREATE");
    if (!fout || fout->IsZombie()) {
        std::cerr << "Error: cannot create output file " << outfile << std::endl;
        fin->Close();
        return;
    }

    std::cout << "Input : " << infile  << std::endl;
    std::cout << "Output: " << outfile << std::endl;

    // Loop over top-level keys (radii: R0.2, R0.3, R0.4, ...)
    TIter nextR(fin->GetListOfKeys());
    TKey *keyR = 0;

    while ((keyR = (TKey *) nextR())) {
        if (std::string(keyR->GetClassName()) != "TDirectoryFile")
            continue;

        std::string rname = keyR->GetName();

        if (rname.size() < 2 || rname[0] != 'R') continue;

        double R = 0.0;
        if (sscanf(rname.c_str(), "R%lf", &R) != 1) {
            std::cerr << "  Warning: could not parse radius from '"
                      << rname << "'. Skipping this directory." << std::endl;
            continue;
        }
        const double dRmax = 0.6 * R;

        double areaMin = 0.0;
        if      (R < 0.25) areaMin = CUT_AREA_02;
        else if (R < 0.35) areaMin = CUT_AREA_03;
        else               areaMin = CUT_AREA_04;

        TDirectoryFile *dirR = (TDirectoryFile *) keyR->ReadObj();
        if (!dirR) continue;

        std::cout << "Radius directory: " << rname
                  << " (R = " << R << ", dRmax = " << dRmax
                  << ", areaMin = " << areaMin << ")" << std::endl;

        fout->cd();
        TDirectory *outR = fout->mkdir(rname.c_str());
        if (!outR) outR = fout->GetDirectory(rname.c_str());

        TIter nextC(dirR->GetListOfKeys());
        TKey *keyC = 0;

        while ((keyC = (TKey *) nextC())) {
            if (std::string(keyC->GetClassName()) != "TDirectoryFile")
                continue;

            std::string cname = keyC->GetName();
            TDirectoryFile *dirC = (TDirectoryFile *) keyC->ReadObj();
            if (!dirC) continue;

            TTree *tree = (TTree *) dirC->Get("JetTree");
            if (!tree) {
                std::cerr << "  Warning: JetTree not found in " << rname
                          << "/" << cname << std::endl;
                continue;
            }

            Long64_t nentries = tree->GetEntries();
            std::cout << "  Centrality: " << cname << " (entries: "
                      << nentries << ")" << std::endl;

            outR->cd();
            TDirectory *outC = outR->mkdir(cname.c_str());
            if (!outC) outC = outR->GetDirectory(cname.c_str());
            outC->cd();

            // ================== selection & weights ==================
            std::string wexpr = "xsecWeight*centralityWeight";

            char bufMatch[256];
            sprintf(bufMatch, "(mc_pt>0 && deltaR<%g)", dRmax);
            std::string matchCond = bufMatch;

            char bufArea[256], bufNef[256];
            sprintf(bufArea, "reco_area>=%g", areaMin);
            sprintf(bufNef,  "reco_neutral_fraction<=%g", CUT_NEUTRAL_FRACTION);

            std::string qualCond  = "(";
            qualCond += bufArea;
            qualCond += ")&&(";
            qualCond += bufNef;
            qualCond += ")";

            std::string fullCut = "(" + wexpr + ")*(" + matchCond + "&&" + qualCond + ")";
            const char *wmatch = fullCut.c_str();

            // ================= histogram booking =====================
            const Int_t    nbins_common = 600;
            const Double_t xmin_common  = 0.0;
            const Double_t xmax_common  = 60.0;

            const Int_t    nbins_corr = 1000;
            const Double_t xmin_corr  = -40.0;
            const Double_t xmax_corr  =  60.0;

            // reco_pt
            std::string h_reco_pt_name  = "h_reco_pt_"  + rname + "_" + cname;
            std::string draw_reco_pt =
                "reco_pt>>" + h_reco_pt_name +
                Form("(%d,%g,%g)", nbins_common, xmin_common, xmax_common);
            tree->Draw(draw_reco_pt.c_str(), wmatch, "goff");

            // reco_pt_corr
            std::string h_reco_ptcorr_name  = "h_reco_ptcorr_"  + rname + "_" + cname;
            std::string draw_reco_ptcorr =
                "reco_pt_corr>>" + h_reco_ptcorr_name +
                Form("(%d,%g,%g)", nbins_corr, xmin_corr, xmax_corr);
            tree->Draw(draw_reco_ptcorr.c_str(), wmatch, "goff");

            // reco_pt_lead
            std::string h_reco_pt_lead_name  = "h_reco_pt_lead_"  + rname + "_" + cname;
            std::string draw_reco_pt_lead =
                "reco_pt_lead>>" + h_reco_pt_lead_name +
                Form("(%d,%g,%g)", nbins_common, xmin_common, xmax_common);
            tree->Draw(draw_reco_pt_lead.c_str(), wmatch, "goff");

            // mc_pt_lead
            std::string h_mc_pt_lead_name  = "h_mc_pt_lead_"  + rname + "_" + cname;
            std::string draw_mc_pt_lead =
                "mc_pt_lead>>" + h_mc_pt_lead_name +
                Form("(%d,%g,%g)", nbins_common, xmin_common, xmax_common);
            tree->Draw(draw_mc_pt_lead.c_str(), wmatch, "goff");
        }
    }

    fout->Write();
    fout->Close();
    fin->Close();

    std::cout << "Done. Histograms (matched jets only, with area+NEF cuts)"
              << " written to " << outfile << std::endl;
}
