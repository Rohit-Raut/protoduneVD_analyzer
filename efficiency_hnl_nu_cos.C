#include "utils_NP02.hpp"
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

TGraph* calculateEfficiency(const std::vector<Double_t>& primaryData,
                            const std::vector<Double_t>& otherData,
                            double thr_min, int nThreshold,
                            double Ntotal) {

    auto it_primary = std::lower_bound(primaryData.begin(), primaryData.end(), thr_min);

    double thr_max = std::max(primaryData.back(), otherData.back());
    auto count_above = [](const std::vector<Double_t>& v, auto from, double thr) {
        return (double)std::distance(std::lower_bound(from, v.end(), thr), v.end());
    };

    TGraph* eff = new TGraph();
    for (int j = 0; j < nThreshold; j++) {
        double thr        = thr_min + j * (thr_max - thr_min) / (nThreshold - 1);
        double countA     = count_above(primaryData, it_primary, thr);
        double efficiency = (Ntotal > 0) ? 100.0 * (countA / Ntotal) : 0.0;
        eff->SetPoint(j, thr, efficiency);
    }
    return eff;
}

// Extract best BDE TA per event from a TTree
void fillBestBDE(TTree* tree, std::vector<Double_t>& out) {
    std::vector<double>*    taADCSumVec      = nullptr;
    std::vector<ULong64_t>* taChannelPeakVec = nullptr;

    tree->SetBranchAddress("adc_integral",  &taADCSumVec);
    tree->SetBranchAddress("channel_peak",  &taChannelPeakVec);

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (!taADCSumVec || !taChannelPeakVec) continue;

        int    bestIdx = -1;
        double bestSum = -1.0;
        for (size_t k = 0; k < taADCSumVec->size(); ++k) {
            if (NP02::isBDE((int)(*taChannelPeakVec)[k]) && (*taADCSumVec)[k] > bestSum) {
                bestSum = (*taADCSumVec)[k];
                bestIdx = (int)k;
            }
        }
        if (bestIdx >= 0)
            out.push_back(bestSum);
    }
}

void efficiencyBDE() {
    TFile* cosmic   = TFile::Open("../Cosmic_Pdune_1MADC_20kticks_10k.root", "READ");
    TFile* neutrino = TFile::Open("../Nu+Cosmic_Pdune_1MADC_20kticks.root",  "READ");
    TFile* hnl = TFile::Open("/home/rraut/trigger_overlay/ta_result.root", "READ");
    //TFile* hnl      = TFile::Open("/project/dune/pdvd_hnls_bsm_flux/HNL_Data/run_hnl_1k_sample.root", "READ");

    if (!cosmic || cosmic->IsZombie() ||
        !neutrino || neutrino->IsZombie() ||
        !hnl || hnl->IsZombie()) {
        std::cout << "Error reading files" << std::endl;
        return;
    }

    TTree* taTreeCosmic   = (TTree*) cosmic->Get("hitdQ/TATree");
    TTree* taTreeNeutrino = (TTree*) neutrino->Get("hitdQ/TATree");
    //TTree* taTreeHNL      = (TTree*) hnl->Get("hitdQ/TATree");
    TTree* taTreeHNL 	= (TTree*)hnl->Get("TATree");
    if (!taTreeCosmic || !taTreeNeutrino || !taTreeHNL) {
        std::cout << "Could not find TA trees" << std::endl;
        return;
    }

    std::vector<Double_t> cosmicBDE;
    std::vector<Double_t> neutrinoBDE;
    std::vector<Double_t> hnlBDE;

    fillBestBDE(taTreeCosmic,   cosmicBDE);
    fillBestBDE(taTreeNeutrino, neutrinoBDE);
    fillBestBDE(taTreeHNL,      hnlBDE);

    std::sort(cosmicBDE.begin(),   cosmicBDE.end());
    std::sort(neutrinoBDE.begin(), neutrinoBDE.end());
    std::sort(hnlBDE.begin(),      hnlBDE.end());

    auto countAbove = [](const std::vector<Double_t>& v, double thr) {
        return std::distance(std::lower_bound(v.begin(), v.end(), thr), v.end());
    };

    std::cout << "\n======== BDE Event Count ========\n";
    std::cout << "Cosmic   events with TA in BDE: " << cosmicBDE.size()   << "\n";
    std::cout << "Neutrino events with TA in BDE: " << neutrinoBDE.size() << "\n";
    std::cout << "HNL      events with TA in BDE: " << hnlBDE.size()      << "\n";
    std::cout << "=================================\n";

    std::cout << "Cosmic   above 1M: " << countAbove(cosmicBDE,   1e6) << "\n";
    std::cout << "Neutrino above 1M: " << countAbove(neutrinoBDE, 1e6) << "\n";
    std::cout << "HNL      above 1M: " << countAbove(hnlBDE,      1e6) << "\n";

    // --- Efficiency curves ---
    const double thr_min = 1e6;
    const double thr_max = 40e6;
    const double steps   = 0.1e6;
    const int nThreshold = (int)((thr_max - thr_min) / steps) + 1;

    // Global trigger efficiency: denominator = total simulated events
    const double Ntotal_cosmic   = 10000.0;
    const double Ntotal_neutrino = 10000.0;
    const double Ntotal_hnl      = 1000.0;

    // For thr_max in calculateEfficiency we need a common upper bound
    // Combine all three into a dummy "other" so thr_max is consistent
    std::vector<Double_t> allCombined;
    allCombined.insert(allCombined.end(), cosmicBDE.begin(),   cosmicBDE.end());
    allCombined.insert(allCombined.end(), neutrinoBDE.begin(), neutrinoBDE.end());
    allCombined.insert(allCombined.end(), hnlBDE.begin(),      hnlBDE.end());
    std::sort(allCombined.begin(), allCombined.end());

    TGraph* effCos = calculateEfficiency(cosmicBDE,   allCombined, thr_min, nThreshold, Ntotal_cosmic);
    TGraph* effNeu = calculateEfficiency(neutrinoBDE, allCombined, thr_min, nThreshold, Ntotal_neutrino);
    TGraph* effHNL = calculateEfficiency(hnlBDE,      allCombined, thr_min, nThreshold, Ntotal_hnl);

    // Spot-check at 8M
    double cosmicEff8M = 100.0 * countAbove(cosmicBDE,   8e6) / Ntotal_cosmic;
    double nuEff8M     = 100.0 * countAbove(neutrinoBDE, 8e6) / Ntotal_neutrino;
    double hnlEff8M    = 100.0 * countAbove(hnlBDE,      8e6) / Ntotal_hnl;

    std::cout << "\nEfficiency at 8M ADC:\n";
    std::cout << "  Cosmic:       " << cosmicEff8M << " %\n";
    std::cout << "  Nu + Cosmic:  " << nuEff8M     << " %\n";
    std::cout << "  HNL:          " << hnlEff8M    << " %\n\n";

    // --- Plot ---
    gStyle->SetOptStat(0);
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelFont(42, "xyz");

    TCanvas* c = new TCanvas("c", "BDE Efficiency", 800, 600);

    TPad* titlePad = new TPad("titlePad", "", 0, 0.90, 1, 1);
    titlePad->SetTopMargin(0);
    titlePad->SetBottomMargin(0);
    titlePad->Draw();
    titlePad->cd();
    TLatex* titletex = new TLatex(0.5, 0.5,
        "#splitline{                ADCSimpleWindow TA Algorithm Efficiency - BDE}"
        "{NP02 Neutrino Trigger Study: Cosmic vs Neutrino + Cosmic vs HNL 20k ticks}");
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

    effCos->SetLineColor(kBlue);  effCos->SetLineWidth(3);
    effCos->SetMarkerStyle(20);   effCos->SetMarkerSize(0.7);
    effCos->SetMarkerColorAlpha(kBlue, 0.6);

    effNeu->SetLineColor(kBlack); effNeu->SetLineWidth(3);
    effNeu->SetMarkerStyle(21);   effNeu->SetMarkerSize(0.7);
    effNeu->SetMarkerColorAlpha(kBlack, 0.6);

    effHNL->SetLineColor(kRed);   effHNL->SetLineWidth(3);
    effHNL->SetMarkerStyle(22);   effHNL->SetMarkerSize(0.7);
    effHNL->SetMarkerColorAlpha(kRed, 0.6);

    effCos->SetTitle("BDE;ADC Integral Sum Cut (ADC);Trigger Efficiency [%]");
    effCos->GetHistogram()->SetMinimum(0);
    effCos->GetHistogram()->SetMaximum(105);
    effCos->GetXaxis()->SetLimits(1e5, 40e6);

    effCos->Draw("ALP");
    effNeu->Draw("LP SAME");
    effHNL->Draw("LP SAME");

    TLegend* leg = new TLegend(0.5, 0.6, 0.87, 0.89);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    leg->SetHeader("BDE (CRP 4+5)", "C");
    leg->AddEntry(effCos, "Cosmic",         "lp");
    leg->AddEntry(effNeu, "#nu + Cosmic",   "lp");
    leg->AddEntry(effHNL, "HNL",            "lp");
    leg->Draw();

    //c->SaveAs("efficiency_BDE_all.png");
}
