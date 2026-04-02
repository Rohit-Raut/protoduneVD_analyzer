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

TGraph* calculateEfficiency(const std::vector<Double_t>& primaryData, double thr_min, int nThreshold) {
    auto it_primary = std::lower_bound(primaryData.begin(), primaryData.end(), thr_min);
    double Nprimary = static_cast<double>(primaryData.size());
    double thr_max  = primaryData.back();

    auto count_above = [](const std::vector<Double_t>& v, auto from, double thr) {
        return (double)std::distance(std::lower_bound(from, v.end(), thr), v.end());
    };

    TGraph* eff = new TGraph();
    for (int j = 0; j < nThreshold; j++) {
        double thr        = thr_min + j * (thr_max - thr_min) / (nThreshold - 1);
        double countA     = count_above(primaryData, it_primary, thr);
        double efficiency = (Nprimary > 0) ? 100.0 * (countA / Nprimary) : 0.0;
        eff->SetPoint(j, thr, efficiency);
    }
    return eff;
}

void Efficiency_cosmic_only() {
    TFile* cosmic = TFile::Open("../yk_cosmic_sample.root", "READ");

    if (!cosmic || cosmic->IsZombie()) {
        std::cout << "Error reading cosmic file" << std::endl;
        return;
    }

    TTree* taTreeCosmic = (TTree*) cosmic->Get("hitdQ/TATree");
    if (!taTreeCosmic) {
        std::cout << "Could not find TA tree" << std::endl;
        return;
    }

    std::vector<double>*    taADCSumVec      = nullptr;
    std::vector<ULong64_t>* taChannelPeakVec = nullptr;

    std::vector<Double_t> cosmicBDE;

    auto isBDE = [](int ch) {
        return (ch >= 0 && ch <= 6143);
    };

    // ===== Cosmic =====
    taTreeCosmic->SetBranchAddress("adc_integral", &taADCSumVec);
    taTreeCosmic->SetBranchAddress("channel_peak", &taChannelPeakVec);

    for (Long64_t i = 0; i < taTreeCosmic->GetEntries(); ++i) {
        taTreeCosmic->GetEntry(i);
        if (!taADCSumVec || !taChannelPeakVec) continue;

        // find single best TA across ALL TAs in event
        int    bestIdx = -1;
        double bestSum = -1.0;
        for (size_t k = 0; k < taADCSumVec->size(); ++k) {
            if ((*taADCSumVec)[k] > bestSum) {
                bestSum = (*taADCSumVec)[k];
                bestIdx = (int)k;
            }
        }

        // only keep if best TA is in BDE
        if (bestIdx >= 0 && isBDE((int)(*taChannelPeakVec)[bestIdx]))
            cosmicBDE.push_back(bestSum);
    }

    std::sort(cosmicBDE.begin(), cosmicBDE.end());

    auto countAbove = [](const std::vector<Double_t>& v, double thr) {
        return std::distance(std::lower_bound(v.begin(), v.end(), thr), v.end());
    };

    std::cout << "\n======== BDE Event Count ========\n";
    std::cout << "Cosmic events total TA in BDE: " << cosmicBDE.size() << "\n";
    std::cout << "Cosmic above 1M ADC: " << countAbove(cosmicBDE, 1e6) << "\n";
    std::cout << "=================================\n";

    const double thr_min    = 1e6;
    const int    nThreshold = 200;

    TGraph* effCos = calculateEfficiency(cosmicBDE, thr_min, nThreshold);

    double cosmicEff8M = 100.0 * countAbove(cosmicBDE, 8e6) / (double)cosmicBDE.size();
    std::cout << "Cosmic Efficiency at 8M ADC: " << cosmicEff8M << "%\n" << std::endl;

    gStyle->SetOptStat(0);
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelFont(42, "xyz");

    TCanvas* c = new TCanvas("c", "BDE Cosmic Efficiency", 800, 600);

    TPad* titlePad = new TPad("titlePad", "", 0, 0.90, 1, 1);
    titlePad->SetTopMargin(0);
    titlePad->SetBottomMargin(0);
    titlePad->Draw();
    titlePad->cd();
    TLatex* titletex = new TLatex(0.5, 0.5,
        "#splitline{       ADCSimpleWindow TA Algorithm Efficiency - BDE}"
        "{NP02 Cosmic Data Only}");
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

    effCos->SetLineColor(kBlue);
    effCos->SetLineWidth(3);
    effCos->SetMarkerStyle(20);
    effCos->SetMarkerSize(0.7);
    effCos->SetMarkerColorAlpha(kBlue, 0.6);

    effCos->SetTitle("BDE;ADC Integral Sum Cut (ADC);Trigger Efficiency [%]");
    effCos->GetHistogram()->SetMinimum(0);
    effCos->GetHistogram()->SetMaximum(105);
    effCos->GetXaxis()->SetLimits(1e5, 40e6);

    effCos->Draw("ALP");

    TLegend* leg = new TLegend(0.5, 0.75, 0.87, 0.89);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    leg->SetHeader("BDE (CRP 4+5)", "C");
    leg->AddEntry(effCos, "Cosmic", "lp");
    leg->Draw();

    //c->SaveAs("efficiency_cosmic_only.png");
}
