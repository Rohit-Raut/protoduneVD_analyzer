#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TPad.h"
#include <iostream>
#include <vector>
#include <algorithm>

// Returns a survival-fraction efficiency graph:
//   for each threshold, eff = (# TAs above threshold) / (total # TAs collected)
TGraph* calculateEfficiency(const std::vector<Double_t>& data, double thr_min, int nThreshold) {
    double Ntotal = static_cast<double>(data.size());
    double thr_max = data.back();

    TGraph* eff = new TGraph();
    for (int j = 0; j < nThreshold; j++) {
        double thr      = thr_min + j * (thr_max - thr_min) / (nThreshold - 1);
        double countAbove = (double)std::distance(
            std::lower_bound(data.begin(), data.end(), thr), data.end());
        double efficiency = (Ntotal > 0) ? 100.0 * (countAbove / Ntotal) : 0.0;
        eff->SetPoint(j, thr, efficiency);
    }
    return eff;
}

// Collect all BDE TAs (channels 0–6143) from a TATree into out_vec.
// Uses ALL TAs per event (not just the best one).
void collectBDE(TTree* tree, std::vector<Double_t>& out_vec) {
    std::vector<double>*    taADCSumVec      = nullptr;
    std::vector<ULong64_t>* taChannelPeakVec = nullptr;

    tree->SetBranchAddress("adc_integral", &taADCSumVec);
    tree->SetBranchAddress("channel_peak", &taChannelPeakVec);

    auto isBDE = [](ULong64_t ch) { return ch <= 6143; };

    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        if (!taADCSumVec || !taChannelPeakVec) continue;
        for (size_t k = 0; k < taADCSumVec->size(); ++k) {
            if (isBDE((*taChannelPeakVec)[k]))
                out_vec.push_back((*taADCSumVec)[k]);
        }
    }
}

void efficiency_BDE_data_only() {
    // ---- Data files (real runs) ----------------------------------------
    // Add or remove entries here to compare different runs.
    std::vector<TString> filePaths = {
        "../Run39350_BDE_Reco_stage1_100.root",
        "../PDVD_Cosmic_Run039693_BDE_Reco_stage1.root",
    };
    std::vector<TString> labels = {
        "Cosmic Run 039350",
        "Cosmic Run 039693",
    };
    std::vector<int> colors = { kBlue, kRed+1 };
    std::vector<int> markers = { 20, 21 };

    const double thr_min    = 1e6;
    const int    nThreshold = 200;

    // ---- Collect BDE TAs per run ----------------------------------------
    std::vector<std::vector<Double_t>> bdeSets(filePaths.size());

    for (size_t r = 0; r < filePaths.size(); ++r) {
        TFile* f = TFile::Open(filePaths[r], "READ");
        if (!f || f->IsZombie()) {
            std::cout << "Error opening: " << filePaths[r] << std::endl;
            return;
        }
        TTree* taTree = (TTree*)f->Get("hitdQ/TATree");
        if (!taTree) {
            std::cout << "Could not find hitdQ/TATree in " << filePaths[r] << std::endl;
            f->Close();
            return;
        }
        collectBDE(taTree, bdeSets[r]);
        std::sort(bdeSets[r].begin(), bdeSets[r].end());

        auto countAbove8M = [](const std::vector<Double_t>& v, double thr) {
            return std::distance(std::lower_bound(v.begin(), v.end(), thr), v.end());
        };
        std::cout << labels[r] << " — total BDE TAs: " << bdeSets[r].size()
                  << "  above 8M: " << countAbove8M(bdeSets[r], 8e6) << "\n";

        // Keep file open until graphs are built; close it afterwards.
        // (TATree lives in memory through the vector)
        f->Close();
    }

    // ---- Build efficiency graphs -----------------------------------------
    std::vector<TGraph*> effGraphs(filePaths.size());
    for (size_t r = 0; r < filePaths.size(); ++r)
        effGraphs[r] = calculateEfficiency(bdeSets[r], thr_min, nThreshold);

    // ---- Style --------------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelFont(42, "xyz");

    // ---- Canvas layout (title pad + plot pad) ------------------------------
    TCanvas* c = new TCanvas("c", "BDE Efficiency — Data Only", 800, 600);

    TPad* titlePad = new TPad("titlePad", "", 0, 0.90, 1, 1);
    titlePad->SetTopMargin(0);
    titlePad->SetBottomMargin(0);
    titlePad->Draw();
    titlePad->cd();
    TLatex* titletex = new TLatex(0.5, 0.5,
        "#splitline{                ADCSimpleWindow TA Algorithm Efficiency - BDE (Data Only)}"
        "{NP02 Trigger Study: Cosmic Data Runs — All TAs in BDE (ch 0#font[122]{-}6143)}");
    titletex->SetTextAlign(22);
    titletex->SetTextSize(0.35);
    titletex->SetTextFont(2);
    titletex->Draw();

    c->cd();
    TPad* plotPad = new TPad("plotPad", "", 0, 0, 1, 0.90);
    plotPad->Draw();
    plotPad->cd();
    plotPad->SetRightMargin(0.08);
    plotPad->SetLeftMargin(0.12);
    plotPad->SetLogy();

    // Style and draw graphs
    for (size_t r = 0; r < effGraphs.size(); ++r) {
        effGraphs[r]->SetLineColor(colors[r]);
        effGraphs[r]->SetLineWidth(3);
        effGraphs[r]->SetMarkerStyle(markers[r]);
        effGraphs[r]->SetMarkerSize(0.7);
        effGraphs[r]->SetMarkerColorAlpha(colors[r], 0.6);
    }

    effGraphs[0]->SetTitle("BDE;ADC Integral Sum Cut (ADC);Trigger Efficiency [%]");
    effGraphs[0]->GetHistogram()->SetMinimum(0);
    effGraphs[0]->GetHistogram()->SetMaximum(105);
    effGraphs[0]->GetXaxis()->SetLimits(1e5, 40e6);

    effGraphs[0]->Draw("ALP");
    for (size_t r = 1; r < effGraphs.size(); ++r)
        effGraphs[r]->Draw("LP SAME");

    // ---- Legend -------------------------------------------------------------
    TLegend* leg = new TLegend(0.5, 0.6, 0.87, 0.89);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    leg->SetHeader("BDE (CRP 4+5) — Data", "C");
    for (size_t r = 0; r < filePaths.size(); ++r)
        leg->AddEntry(effGraphs[r], labels[r], "lp");
    leg->Draw();

    c->SaveAs("efficiency_BDE_data_only.png");
    std::cout << "Saved efficiency_BDE_data_only.png\n";
}
