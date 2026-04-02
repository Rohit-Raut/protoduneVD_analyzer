#include "utils_NP02.hpp"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLine.h"
#include <iostream>
#include <vector>
#include <algorithm>
TGraph* calculateEfficiency(const std::vector<Double_t>& primaryData, const std::vector<Double_t>& otherData, double thr_min, int nThreshold, double Nprimary) {
    auto it_primary = std::lower_bound(primaryData.begin(), primaryData.end(), thr_min);

    double thr_max = std::max(primaryData.back(), otherData.back());

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

void fillBDEVector(TTree* tree, std::vector<Double_t>& bdeVec) {
    std::vector<double>*    taADCSumVec      = nullptr;
    std::vector<ULong64_t>* taChannelPeakVec = nullptr;

    tree->SetBranchAddress("adc_integral", &taADCSumVec);
    tree->SetBranchAddress("channel_peak", &taChannelPeakVec);

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
            bdeVec.push_back(bestSum);
    }
    std::sort(bdeVec.begin(), bdeVec.end());
}

void efficiency_hnl_nu_cos() {

    TFile* cosmic   = TFile::Open("../Cosmic_Pdune_1MADC_20kticks_10k.root", "READ");
    TFile* neutrino = TFile::Open("../Nu+Cosmic_Pdune_1MADC_20kticks.root", "READ");
    TFile* hnl      = TFile::Open("../hnl_cosmic_nu_mc.root", "READ");

    if (!cosmic || cosmic->IsZombie() || !neutrino || neutrino->IsZombie() || !hnl || hnl->IsZombie()) {
        std::cout << "Error reading files" << std::endl;
        return;
    }

    TTree* taTreeCosmic   = (TTree*) cosmic->Get("hitdQ/TATree");
    TTree* taTreeNeutrino = (TTree*) neutrino->Get("hitdQ/TATree");
    TTree* taTreeHNL      = (TTree*) hnl->Get("hitdQ/TATree");

    if (!taTreeCosmic || !taTreeNeutrino || !taTreeHNL) {
        std::cout << "Could not find TA trees" << std::endl;
        return;
    }

    std::vector<Double_t> cosmicBDE;
    std::vector<Double_t> neutrinoBDE;
    std::vector<Double_t> hnlBDE;

    fillBDEVector(taTreeCosmic,   cosmicBDE);
    fillBDEVector(taTreeNeutrino, neutrinoBDE);
    fillBDEVector(taTreeHNL,      hnlBDE);

    auto countAbove = [](const std::vector<Double_t>& v, double thr) {
        return std::distance(std::lower_bound(v.begin(), v.end(), thr), v.end());
    };

    std::cout << "\n======== BDE Event Count ========\n";
    std::cout << "Cosmic   events with TA in BDE: " << cosmicBDE.size()   << "\n";
    std::cout << "Neutrino events with TA in BDE: " << neutrinoBDE.size() << "\n";
    std::cout << "HNL      events with TA in BDE: " << hnlBDE.size()     << "\n";
    std::cout << "=================================\n\n";

    std::cout << "Cosmic   above 8M: " << countAbove(cosmicBDE,   8e6) << "\n";
    std::cout << "Neutrino above 8M: " << countAbove(neutrinoBDE, 8e6) << "\n";
    std::cout << "HNL      above 8M: " << countAbove(hnlBDE,      8e6) << "\n\n";

    double cosmicEff8M = 100.0 * countAbove(cosmicBDE,   8e6) / 10000.0;
    double nuEff8M     = 100.0 * countAbove(neutrinoBDE, 8e6) / 10000.0;
    double hnlEff8M    = 100.0 * countAbove(hnlBDE,      8e6) / 935.0;

    std::cout << "Cosmic   efficiency at 8M: " << cosmicEff8M << " %\n";
    std::cout << "Neutrino efficiency at 8M: " << nuEff8M     << " %\n";
    std::cout << "HNL      efficiency at 8M: " << hnlEff8M    << " %\n\n";

    const double thr_min = 1e6;
    const double thr_max = 40e6;
    const double steps   = 0.1e6;
    const int nThreshold = (int)((thr_max - thr_min) / steps) + 1;

    TGraph* effCos = calculateEfficiency(cosmicBDE,   neutrinoBDE, thr_min, nThreshold, 10000.0);
    TGraph* effNeu = calculateEfficiency(neutrinoBDE, cosmicBDE,   thr_min, nThreshold, 10000.0);
    TGraph* effHNL = calculateEfficiency(hnlBDE,      cosmicBDE,   thr_min, nThreshold, 935.0);

    gStyle->SetOptStat(0);
    gStyle->SetTextFont(42);
    gStyle->SetTitleFont(42, "");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetLabelFont(42, "xyz");

    TCanvas* c = new TCanvas("c", "BDE Efficiency", 800, 600);

//    TPad* titlePad = new TPad("titlePad", "", 0, 0.90, 1, 1);
//    titlePad->SetTopMargin(0);
//    titlePad->SetBottomMargin(0);
//    titlePad->Draw();
//    titlePad->cd();
//    TLatex* titletex = new TLatex(0.5, 0.5,
//        "#splitline{                ADCSimpleWindow TA Algorithm Efficiency - BDE}"
//        "{NP02 Trigger Study: Cosmic vs #nu + Cosmic vs HNL + #nu + Cosmic | 20k ticks}");
//    titletex->SetTextAlign(22);
//    titletex->SetTextSize(0.35);
//    titletex->SetTextFont(2);
//    titletex->Draw();

//    c->cd();
    TPad* plotPad = new TPad("plotPad", "", 0, 0, 1, 1);
    plotPad->Draw();
    plotPad->cd();
    plotPad->SetTopMargin(0.06);
    plotPad->SetRightMargin(0.08);
    plotPad->SetLeftMargin(0.12);
    plotPad->SetGridy();

    effCos->SetLineColor(kBlue);  effCos->SetLineWidth(3);
    effCos->SetMarkerStyle(20);   effCos->SetMarkerSize(0.7);
    effCos->SetMarkerColorAlpha(kBlue, 0.6);

    effNeu->SetLineColor(kBlack); effNeu->SetLineWidth(3);
    effNeu->SetMarkerStyle(21);   effNeu->SetMarkerSize(0.7);
    effNeu->SetMarkerColorAlpha(kBlack, 0.6);

    effHNL->SetLineColor(kRed);   effHNL->SetLineWidth(3);
    effHNL->SetMarkerStyle(22);   effHNL->SetMarkerSize(0.7);
    effHNL->SetMarkerColorAlpha(kRed, 0.6);

    effCos->SetTitle(";ADC Integral Sum Cut (ADC);Trigger Efficiency [%]");
    effCos->GetHistogram()->SetMinimum(0);
    effCos->GetHistogram()->SetMaximum(105);
    effCos->GetXaxis()->SetLimits(1e5, 40e6);
    auto getEff = [](TGraph* g, double x){
	Double_t *xp = g->GetX(), *yp = g->GetY();
	for(int i = 0; i<g->GetN()-1; i++){
	    if(xp[i]<=x && xp[i+1]>=x){
		double f =  (x-xp[i])/(xp[i+1] - xp[i]);
		return yp[i]+f*(yp[i+1]-yp[i]);
	    }
	}
	return -1.0;
    };

    double thrLine[] = {5e6, 6e6, 7e6, 8e6};
    for (int t = 0; t<4; t++){
	TLine* l =new TLine(thrLine[t], 0, thrLine[t], 105);
	l->SetLineColor(kGray+2);
	l->SetLineStyle(7);
	l->SetLineWidth(1);
	l->Draw();

	double eCos = getEff(effCos, thrLine[t]);
	double eNU = getEff(effNeu, thrLine[t]);
	double eHNL = getEff(effHNL, thrLine[t]);

	TLatex lt;
	lt.SetTextSize(0.025);
	lt.SetTextFont(42);
	lt.SetTextAlign(12);
        double xoff = thrLine[t] + 0.3e6;

        lt.SetTextColor(kBlue);
        lt.DrawLatex(xoff, eCos, Form("%.1f%%", eCos));
        lt.SetTextColor(kBlack);
        lt.DrawLatex(xoff, eNU, Form("%.1f%%", eNU));
        lt.SetTextColor(kRed);
        lt.DrawLatex(xoff, eHNL, Form("%.1f%%", eHNL));

    }
        

    effCos->Draw("ALP");
    effNeu->Draw("LP SAME");
    effHNL->Draw("LP SAME");



    TLegend* leg = new TLegend(0.5, 0.55, 0.87, 0.89);
    leg->SetTextFont(42);
    leg->SetTextSize(0.035);
    leg->SetBorderSize(0);
    leg->SetFillStyle(4000);
    //leg->SetHeader("BDE (CRP 4+5)", "C");
    leg->AddEntry(effCos, "Cosmic",              "lp");
    leg->AddEntry(effNeu, "#nu + Cosmic",        "lp");
    leg->AddEntry(effHNL, "HNL + #nu + Cosmic",  "lp");
    leg->Draw();

    //TPad
    c->cd();
    TPad* headerPad = new TPad("heaaderPad", "", 0, 0.95, 1, 1);
    //headerPad->SetTopMargin(0);
    //headerPad->SetBottomMargin(0);
    //headerPad->SetLeftMargin(0.12);
    //headerPad->SetRightMargin(0.15);
    headerPad->SetFillStyle(4000);
    headerPad->SetFillColor(0);
    headerPad->Draw();
    headerPad->cd();

    TLatex tex;
    tex.SetNDC();
    tex.SetTextFont(42);
    tex.SetTextSize(0.65);

    tex.SetTextAlign(11);
    tex.DrawLatex(0.12, 0.15, "#bf{DUNE} #it{Simulation}");
    tex.SetTextAlign(21);
    tex.DrawLatex(0.50, 0.15, "BDE (CRP 4 + CRP 5)");
    tex.SetTextAlign(31);
    tex.DrawLatex(0.88, 0.14, "#bf{NP02}");


    c->SaveAs("efficiency_BDE_hnl_nu_cos.png");
}
