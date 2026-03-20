#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iostream>
#include <vector>
#include <algorithm>

// --- Configuration: CRP Channel Ranges (Corrected to avoid overlap) ---
// Standard ProtoDUNE-VD mapping: 3072 channels per CRP
const int CRP4_MIN = 0;
const int CRP4_MAX = 3071;

const int CRP5_MIN = 3072;
const int CRP5_MAX = 6143;

// --- Helper Function: Reads file, counts hits within Channel Range ---
TH1F* GetHitDistribution(const char* filename, const char* histName, int color, int style, int minCh, int maxCh) {

    TFile *f = TFile::Open(filename, "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Error: Cannot open " << filename << std::endl;
        return nullptr;
    }

    TTree *t = (TTree*)f->Get("hitdQ/Hit");
    if (!t) {
        f->Close();
        return nullptr;
    }

    Int_t run, subrun, event;
    UInt_t channel; 

    t->SetBranchAddress("run", &run);
    t->SetBranchAddress("subrun", &subrun);
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("Channel", &channel);

    std::vector<int> counts;
    Long64_t nentries = t->GetEntries();

    int current_event = -1, current_run = -1, current_subrun = -1;
    int hits_in_valid_range = 0;

    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);

        // Define unique event ID
        if (i == 0) {
            current_event = event; current_run = run; current_subrun = subrun;
        }

        // Check for Event/Subrun/Run transition
        if (event != current_event || subrun != current_subrun || run != current_run) {
            counts.push_back(hits_in_valid_range);
            hits_in_valid_range = 0;
            current_event = event; current_run = run; current_subrun = subrun;
        }

        if ((int)channel >= minCh && (int)channel <= maxCh) {
            hits_in_valid_range++;
        }
    }
    // Push the very last event in the file
    if (nentries > 0) counts.push_back(hits_in_valid_range);

    // Create Histogram with reasonable binning for multiplicity
    TH1F *h = new TH1F(histName, ";Hit Multiplicity (Hits/Event);Probability Density", 100, 0, 15000);
    for (int c : counts) h->Fill(c);

    // Styling
    h->SetLineColor(color);
    h->SetLineWidth(2);
    h->SetLineStyle(style);
    h->SetDirectory(0); // Detach from file so it stays in memory after file closes

    // Normalize to unit area to compare shapes between different file sizes
    if (h->Integral() > 0) h->Scale(1.0 / h->Integral());

    f->Close();
    return h;
}

// --- Main Function ---
void compare_CRP_hits() {
    gStyle->SetOptStat(0); 

    // Define file paths (using your provided names)
    const char* f_mc_prod = "../PDVD_CosmicOnly_yk_sample_BDE_Reco_stage1.root";
    const char* f_data    = "../Run039636_BDE_Reco_stage1_100.root";
    const char* f_mc_nu   = "../PDVD_Nu_MC_BDE_Reco_stage1.root";

    // --- CRP 4 Plotting ---
    TCanvas *c4 = new TCanvas("c_crp4", "CRP 4 Comparison", 900, 700);
    c4->SetGrid();

    TH1F *h4_data = GetHitDistribution(f_data,    "h4_data", kBlack, 1, CRP4_MIN, CRP4_MAX);
    TH1F *h4_mc1  = GetHitDistribution(f_mc_prod, "h4_mc1",  kRed,   1, CRP4_MIN, CRP4_MAX);
    TH1F *h4_mc2  = GetHitDistribution(f_mc_nu,   "h4_mc2",  kBlue,  1, CRP4_MIN, CRP4_MAX);

    if (h4_data && h4_mc1 && h4_mc2) {
        h4_data->SetTitle("CRP 4 Hit Multiplicity: Data vs Simulation;Hits/Event;Normalized Area");
        
        double maxVal = std::max({h4_data->GetMaximum(), h4_mc1->GetMaximum(), h4_mc2->GetMaximum()});
        h4_data->GetYaxis()->SetRangeUser(0, maxVal * 1.3);

        h4_data->Draw("HIST");
        h4_mc1->Draw("HIST SAME");
        h4_mc2->Draw("HIST SAME");

        TLegend *leg4 = new TLegend(0.55, 0.65, 0.88, 0.88);
        leg4->SetBorderSize(0);
        leg4->SetHeader("CRP 4 (Bottom)", "C");
        leg4->AddEntry(h4_data, "Data (Run 39636)", "l");
        leg4->AddEntry(h4_mc1,  "MC: Cosmic Only", "l");
        leg4->AddEntry(h4_mc2,  "MC: Cosmic + Nu", "l");
        leg4->Draw();
    }

    // --- CRP 5 Plotting ---
    TCanvas *c5 = new TCanvas("c_crp5", "CRP 5 Comparison", 900, 700);
    c5->SetGrid();

    TH1F *h5_data = GetHitDistribution(f_data,    "h5_data", kBlack, 1, CRP5_MIN, CRP5_MAX);
    TH1F *h5_mc1  = GetHitDistribution(f_mc_prod, "h5_mc1",  kRed,   1, CRP5_MIN, CRP5_MAX);
    TH1F *h5_mc2  = GetHitDistribution(f_mc_nu,   "h5_mc2",  kBlue,  1, CRP5_MIN, CRP5_MAX);

    if (h5_data && h5_mc1 && h5_mc2) {
        h5_data->SetTitle("CRP 5 Hit Multiplicity: Data vs Simulation;Hits/Event;Normalized Area");
        
        double maxVal = std::max({h5_data->GetMaximum(), h5_mc1->GetMaximum(), h5_mc2->GetMaximum()});
        h5_data->GetYaxis()->SetRangeUser(0, maxVal * 1.3);

        h5_data->Draw("HIST");
        h5_mc1->Draw("HIST SAME");
        h5_mc2->Draw("HIST SAME");

        TLegend *leg5 = new TLegend(0.55, 0.65, 0.88, 0.88);
        leg5->SetBorderSize(0);
        leg5->SetHeader("CRP 5 (Bottom)", "C");
        leg5->AddEntry(h5_data, "Data (Run 39350)", "l");
        leg5->AddEntry(h5_mc1,  "MC: Cosmic Only", "l");
        leg5->AddEntry(h5_mc2,  "MC: Cosmic + Nu", "l");
        leg5->Draw();
    }
}
