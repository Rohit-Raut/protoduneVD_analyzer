#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <algorithm>

void compare_drift_corrected() {
    gStyle->SetOptStat(0);
    const double T2C = 0.08;
    const Long64_t MAX = 200000000;

    auto h_mc = new TH1D("mc", "", 500, 0, 800);
    auto h_da = new TH1D("da", "", 500, 0, 800);
    h_mc->SetLineColor(kAzure+1); h_mc->SetLineWidth(2);
    h_da->SetLineColor(kRed+1);   h_da->SetLineWidth(2);

    // ===== MC (vectors) =====
    auto f1 = TFile::Open("../data/Cosmic_MC_allChannel.root");
    TTree* t1 = (TTree*)f1->Get("hitdQ/Hit");
    std::vector<UInt_t>  *v_ch = 0;
    std::vector<Float_t> *v_pt = 0;
    std::vector<Int_t> *v_st = 0;
    t1->SetBranchAddress("Channel",    &v_ch);
    t1->SetBranchAddress("fPeakTime",  &v_pt);
    t1->SetBranchAddress("fStartTick", &v_st);

    Long64_t n_mc = 0;
    for (int j = 0; j < t1->GetEntries() && n_mc < MAX; ++j) {
        t1->GetEntry(j)
        if (!v_ch || !v_pt || !v_st) continue;
        for (size_t i = 0; i < v_ch->size() && n_mc < MAX; ++i) {
            if (v_ch->at(i) <= 6143) {
                double drift = (v_pt->at(i) - v_st->at(i)) * T2C;
                h_mc->Fill(drift);
                n_mc++;
            }
        }
    }

    // ===== Data (scalars) =====
    auto f2 = TFile::Open("../PDVD_Cosmic_Run039693_BDE_Reco_stage1.root");
    TTree* t2 = (TTree*)f2->Get("hitdQ/Hit");
    UInt_t  d_ch;
    Float_t d_pt;
    Int_t d_st;
    t2->SetBranchAddress("Channel",    &d_ch);
    t2->SetBranchAddress("fPeakTime",  &d_pt);
    t2->SetBranchAddress("fStartTick", &d_st);

    Long64_t n_da = 0;
    for (int j = 0; j < t2->GetEntries() && n_da < MAX; ++j) {
        t2->GetEntry(j);
        if (d_ch <= 6143) {
            double drift = (d_pt - d_st) * T2C;
            h_da->Fill(drift);
            n_da++;
        }
    }

    // ===== Draw =====
    auto c = new TCanvas("c", "Drift Corrected", 800, 600);
    c->SetLogy();

    double ymax = std::max(h_mc->GetMaximum(), h_da->GetMaximum());
    h_da->SetMaximum(ymax * 5);
    h_da->SetMinimum(1);
    h_da->GetXaxis()->SetTitle("(fPeakTime - fStartTick) #times 0.08  [cm]");
    h_da->GetYaxis()->SetTitle("Number of Hits");
    h_da->SetTitle("Corrected Drift Distance (t_{0} offset removed);(fPeakTime - fStartTick) #times 0.08 [cm];Number of Hits");

    h_da->Draw("HIST");
    h_mc->Draw("HIST SAME");

    auto leg = new TLegend(0.60, 0.75, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->AddEntry(h_mc, "CORSIKA MC Sim", "l");
    leg->AddEntry(h_da, "Data Run 039693", "l");
    leg->Draw();

    //c->SaveAs("drift_corrected.png");
}
