// efficiency_all_DUNE_style.C
// Trigger efficiency for BDE (CRP 4+5): Cosmic, Nu+Cosmic, HNL+Nu+Cosmic
// ADCSimpleWindow TA algorithm, NP02 Neutrino Trigger Study
// Professional DUNE publication-style plot

#include "utils_NP02.hpp"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLine.h"
#include <iostream>
#include <vector>
#include <algorithm>

// Fill vector with best BDE TA ADC integral per event
void fillBDEVector(TTree* tree, std::vector<Double_t>& vec) {
    std::vector<double>*    adcVec = nullptr;
    std::vector<ULong64_t>* chVec  = nullptr;
    tree->SetBranchAddress("adc_integral", &adcVec);
    tree->SetBranchAddress("channel_peak", &chVec);
    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        if (!adcVec || !chVec) continue;
        int    bestIdx = -1;
        double bestVal = -1.0;
        for (size_t k = 0; k < adcVec->size(); ++k) {
            if (NP02::isBDE((int)(*chVec)[k]) && (*adcVec)[k] > bestVal) {
                bestVal = (*adcVec)[k];
                bestIdx = (int)k;
            }
        }
        if (bestIdx >= 0) vec.push_back(bestVal);
    }
    std::sort(vec.begin(), vec.end());
}

// Build efficiency-vs-threshold TGraph
// Ntot: total number of simulated events (denominator)
TGraph* makeEffCurve(const std::vector<Double_t>& data,
                     const std::vector<Double_t>& ref,
                     double thr_min, int nPts, double Ntot) {
    double thr_max = std::max(data.back(), ref.back());
    TGraph* g = new TGraph();
    for (int j = 0; j < nPts; ++j) {
        double thr = thr_min + j * (thr_max - thr_min) / (nPts - 1);
        double n   = (double)std::distance(
                         std::lower_bound(data.begin(), data.end(), thr),
                         data.end());
        g->SetPoint(j, thr, 100.0 * n / Ntot);
    }
    return g;
}

// Linear interpolation on TGraph at a given x
double interpEff(TGraph* g, double x) {
    Double_t *xp = g->GetX(), *yp = g->GetY();
    for (int i = 0; i < g->GetN() - 1; ++i) {
        if (xp[i] <= x && xp[i+1] >= x) {
            double f = (x - xp[i]) / (xp[i+1] - xp[i]);
            return yp[i] + f * (yp[i+1] - yp[i]);
        }
    }
    return -1.0;
}

void efficiency_all_DUNE() {

    // ----------------------------------------------------------------
    //  Open files
    // ----------------------------------------------------------------
    TFile* fCos = TFile::Open("../Cosmic_Pdune_1MADC_20kticks_10k.root", "READ");
    TFile* fNu  = TFile::Open("../Nu+Cosmic_Pdune_1MADC_20kticks.root",  "READ");
    TFile* fHNL = TFile::Open("../hnl_cosmic_nu_mc.root",                "READ");

    if (!fCos || fCos->IsZombie() ||
        !fNu  || fNu ->IsZombie() ||
        !fHNL || fHNL->IsZombie()) {
        std::cerr << "Error: cannot open one or more input files.\n";
        return;
    }

    TTree* tCos = (TTree*) fCos->Get("hitdQ/TATree");
    TTree* tNu  = (TTree*) fNu ->Get("hitdQ/TATree");
    TTree* tHNL = (TTree*) fHNL->Get("hitdQ/TATree");

    if (!tCos || !tNu || !tHNL) {
        std::cerr << "Error: could not retrieve TA trees.\n";
        return;
    }

    // ----------------------------------------------------------------
    //  Fill data vectors
    // ----------------------------------------------------------------
    std::vector<Double_t> vCos, vNu, vHNL;
    fillBDEVector(tCos, vCos);
    fillBDEVector(tNu,  vNu);
    fillBDEVector(tHNL, vHNL);

    // Total simulated events per sample (denominator for efficiency)
    const double N_cos = 10000.0;
    const double N_nu  = 10000.0;
    const double N_hnl = 935.0;

    auto nAbove = [](const std::vector<Double_t>& v, double thr) -> ptrdiff_t {
        return std::distance(std::lower_bound(v.begin(), v.end(), thr), v.end());
    };

    std::cout << "\n============ BDE Summary ============\n"
              << "  Cosmic   events with TA in BDE : " << vCos.size()  << "\n"
              << "  Nu+Cos   events with TA in BDE : " << vNu.size()   << "\n"
              << "  HNL      events with TA in BDE : " << vHNL.size()  << "\n"
              << "  ---- Efficiency at 8M ADC cut ----\n"
              << "  Cosmic   : " << 100.0 * nAbove(vCos, 8e6) / N_cos  << " %\n"
              << "  Nu+Cos   : " << 100.0 * nAbove(vNu,  8e6) / N_nu   << " %\n"
              << "  HNL      : " << 100.0 * nAbove(vHNL, 8e6) / N_hnl  << " %\n"
              << "=====================================\n\n";

    // ----------------------------------------------------------------
    //  Build efficiency curves
    // ----------------------------------------------------------------
    const double thr_min = 1.0e6, thr_max = 40.0e6, step = 0.1e6;
    const int    nPts    = (int)((thr_max - thr_min) / step) + 1;

    TGraph* gCos = makeEffCurve(vCos, vNu,  thr_min, nPts, N_cos);
    TGraph* gNu  = makeEffCurve(vNu,  vCos, thr_min, nPts, N_nu);
    TGraph* gHNL = makeEffCurve(vHNL, vCos, thr_min, nPts, N_hnl);

    // ----------------------------------------------------------------
    //  DUNE publication style
    // ----------------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTickX(1);       // ticks on all 4 sides
    gStyle->SetPadTickY(1);
    gStyle->SetPadGridY(1);       // horizontal grid only
    gStyle->SetGridColor(17);     // light gray
    gStyle->SetGridStyle(3);      // dotted
    gStyle->SetGridWidth(1);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetTitleFont(42, "xyz");

    // Colour palette: blue / black / red
    const Color_t colCos = kAzure + 1;
    const Color_t colNu  = kBlack;
    const Color_t colHNL = kRed   + 1;

    // Graph styling
    gCos->SetLineColor(colCos);  gCos->SetLineWidth(3);
    gCos->SetMarkerColor(colCos); gCos->SetMarkerStyle(20); gCos->SetMarkerSize(0.65);

    gNu ->SetLineColor(colNu);   gNu ->SetLineWidth(3);
    gNu ->SetMarkerColor(colNu);  gNu ->SetMarkerStyle(21); gNu ->SetMarkerSize(0.65);

    gHNL->SetLineColor(colHNL);  gHNL->SetLineWidth(3);
    gHNL->SetMarkerColor(colHNL); gHNL->SetMarkerStyle(22); gHNL->SetMarkerSize(0.80);

    // ----------------------------------------------------------------
    //  Canvas  –  single pad, top margin reserved for DUNE header text
    // ----------------------------------------------------------------
    TCanvas* c = new TCanvas("c_eff_all", "BDE Trigger Efficiency", 900, 650);
    c->SetFillColor(kWhite);
    c->SetTopMargin(0.105);
    c->SetBottomMargin(0.125);
    c->SetLeftMargin(0.130);
    c->SetRightMargin(0.050);

    // ---- Draw graphs (ALP sets up the axis frame) ----
    gCos->Draw("ALP");

    // Axis range and formatting
    gCos->GetXaxis()->SetTitle("ADC Integral Sum Threshold");
    gCos->GetYaxis()->SetTitle("Trigger Efficiency [%]");
    gCos->GetXaxis()->SetLimits(1.0e6, 40.0e6);
    gCos->GetHistogram()->SetMinimum(0.0);
    gCos->GetHistogram()->SetMaximum(108.0);

    gCos->GetXaxis()->SetTitleFont(42);  gCos->GetXaxis()->SetLabelFont(42);
    gCos->GetXaxis()->SetTitleSize(0.052); gCos->GetXaxis()->SetLabelSize(0.043);
    gCos->GetXaxis()->SetTitleOffset(1.05);
    gCos->GetXaxis()->SetNdivisions(510);

    gCos->GetYaxis()->SetTitleFont(42);  gCos->GetYaxis()->SetLabelFont(42);
    gCos->GetYaxis()->SetTitleSize(0.052); gCos->GetYaxis()->SetLabelSize(0.043);
    gCos->GetYaxis()->SetTitleOffset(1.08);
    gCos->GetYaxis()->SetNdivisions(510);

    gNu ->Draw("LP SAME");
    gHNL->Draw("LP SAME");

    c->Update();   // update before drawing objects in data coordinates

    // ----------------------------------------------------------------
    //  Vertical reference lines with per-curve efficiency annotations
    // ----------------------------------------------------------------
    const double thrRef[] = {5e6, 6e6, 7e6, 8e6};
    const int    nRef     = 4;

    TLatex lt;
    lt.SetTextFont(42);
    lt.SetTextSize(0.027);

    for (int t = 0; t < nRef; ++t) {
        // dashed vertical line
        TLine* l = new TLine(thrRef[t], 0.0, thrRef[t], 106.0);
        l->SetLineColor(kGray + 2);
        l->SetLineStyle(7);
        l->SetLineWidth(1);
        l->Draw();

        double eCos = interpEff(gCos, thrRef[t]);
        double eNu  = interpEff(gNu,  thrRef[t]);
        double eHNL = interpEff(gHNL, thrRef[t]);

        double xLab = thrRef[t] + 0.28e6;

        // stagger vertically to avoid overlap: Cosmic slightly above,
        // Nu slightly below its natural position, HNL as-is
        lt.SetTextAlign(12);
        lt.SetTextColor(colCos);  lt.DrawLatex(xLab, eCos + 1.2, Form("%.1f%%", eCos));
        lt.SetTextColor(colNu);   lt.DrawLatex(xLab, eNu  - 2.8, Form("%.1f%%", eNu));
        lt.SetTextColor(colHNL);  lt.DrawLatex(xLab, eHNL + 1.2, Form("%.1f%%", eHNL));
    }

    // ----------------------------------------------------------------
    //  Legend
    // ----------------------------------------------------------------
    TLegend* leg = new TLegend(0.515, 0.50, 0.915, 0.775);
    leg->SetTextFont(42);
    leg->SetTextSize(0.040);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetHeader("BDE (CRP 4 + CRP 5)", "C");
    leg->AddEntry(gCos, "Cosmic",             "lp");
    leg->AddEntry(gNu,  "#nu + Cosmic",       "lp");
    leg->AddEntry(gHNL, "HNL + #nu + Cosmic", "lp");
    leg->Draw();

    // ----------------------------------------------------------------
    //  DUNE header line  (drawn in canvas NDC, above the plot frame)
    // ----------------------------------------------------------------
    TLatex hdr;
    hdr.SetNDC();
    hdr.SetTextFont(42);

    // Left: "DUNE  Simulation"
    hdr.SetTextAlign(11);
    hdr.SetTextSize(0.047);
    hdr.DrawLatex(0.130, 0.912, "#bf{DUNE} #it{Simulation}");

    // Centre: sub-title
    hdr.SetTextAlign(21);
    hdr.SetTextSize(0.036);
    hdr.DrawLatex(0.490, 0.912, "ADCSimpleWindow TA  |  BDE  |  20k ticks");

    // Right: experiment label
    hdr.SetTextAlign(31);
    hdr.SetTextSize(0.042);
    hdr.DrawLatex(0.950, 0.912, "#bf{NP02}");

    c->Update();
    c->SaveAs("efficiency_all_DUNE_style.png");
    c->SaveAs("efficiency_all_DUNE_style.pdf");

    std::cout << "Saved: efficiency_all_DUNE_style.png / .pdf\n";
}
