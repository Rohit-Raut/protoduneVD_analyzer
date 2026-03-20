#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iostream>
#include <vector>
#include <algorithm>

// --- Helper Function: Reads file, counts hits, returns Histogram ---
TH1F* GetHitDistribution(const char* filename, const char* histName, int color, int style) {
    
    TFile *f = TFile::Open(filename);
    if (!f || f->IsZombie()) {
        std::cout << "Error: Cannot open " << filename << std::endl;
        return nullptr;
    }

    // Get Tree (assuming same structure for all)
    TTree *t = (TTree*)f->Get("hitdQ/Hit");
    if (!t) return nullptr;

    int run, subrun, event;
    t->SetBranchAddress("run", &run);
    t->SetBranchAddress("subrun", &subrun);
    t->SetBranchAddress("event", &event);

    std::vector<int> counts;
    Long64_t nentries = t->GetEntries();
    
    int current_event = -1, current_run = -1, current_subrun = -1;
    int hits_in_event = 0;

    // Loop to extract counts
    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);
        if (event != current_event || subrun != current_subrun || run != current_run) {
            if (i > 0) counts.push_back(hits_in_event);
            current_event = event; current_subrun = subrun; current_run = run;
            hits_in_event = 1;
        } else {
            hits_in_event++;
        }
    }
    if (nentries > 0) counts.push_back(hits_in_event);

    // Find max X for this specific file (to ensure binning covers it)
    int local_max = 0;
    if (!counts.empty()) local_max = *std::max_element(counts.begin(), counts.end());

    // Create Histogram (We use a fixed range 0-10000 for comparison consistency)
    // You can change 10000 to local_max if you prefer dynamic
    TH1F *h = new TH1F(histName, ";Hit Multiplicity (Hits/Event);Events", 200, 0, 40000);
    
    for (int c : counts) h->Fill(c);

    // Styling
    h->SetLineColor(color);
    h->SetLineWidth(2);
    h->SetLineStyle(style);
    
    // IMPORTANT: Normalize to Unit Area for Shape Comparison
    // If you want raw counts, comment this line out.
    // But usually Data vs MC has different total events, so we normalize to compare shape.
    //if (h->Integral() > 0) h->Scale(1.0 / h->Integral());

    // Note: We are NOT closing the file 'f' here because that would delete the histogram
    // associated with it. In a complex script, we'd clone the hist or manage directory ownership.
    return h;
}

// --- Main Function ---
void compareHits() {
    // 1. Setup Canvas and Style
    TCanvas *c1 = new TCanvas("c1", "Data vs MC Comparison", 900, 700);
    gStyle->SetOptStat(0); // Remove the stats box (cleaner for overlays)

    // 2. Load Histograms using Helper
    // MC Sample 1 (Red)
    TH1F *h_mc1 = GetHitDistribution("../PDVD_CosmicOnly_yk_sample_BDE_Reco_stage1.root", "h_mc1", kRed, 1);
    
    // Data (Black)
    TH1F *h_data = GetHitDistribution("../Run39350_BDE_Reco_stage1_100.root", "h_data", kBlack, 1);
    
    // MC Sample 2 (Blue)
    TH1F *h_mc2 = GetHitDistribution("../PDVD_Nu_MC_BDE_Reco_stage1.root", "h_mc2", kBlue, 1);

    if (!h_mc1 || !h_data || !h_mc2) {
        std::cout << "Error: One or more files failed to load." << std::endl;
        return;
    }

    // 3. Set the Title
    h_mc1->SetTitle("Reconstructed Hit Multiplicity: ProtoDUNE-VD Data vs Simulation");

    // 4. Determine Y-Axis Range
    // We need to check which histogram has the highest peak so nothing gets cut off
    double max1 = h_mc1->GetMaximum();
    double max2 = h_data->GetMaximum();
    double max3 = h_mc2->GetMaximum();
    double global_max = std::max({max1, max2, max3});

    // Set the Y-axis slightly higher than the highest peak (plus 10% buffer)
    h_mc1->GetYaxis()->SetRangeUser(0, global_max * 1.2);

    // 5. Draw Overlays
    // Draw MC1 first to set the axes
    h_mc1->Draw("HIST"); 
    // Draw Data (Same axes)
    h_data->Draw("HIST SAME");
    // Draw MC2 (Same axes)
    h_mc2->Draw("HIST SAME");

    // 6. Add Legend
    TLegend *leg = new TLegend(0.55, 0.65, 0.88, 0.88); // x1, y1, x2, y2 (Top Right)
    leg->SetBorderSize(0); // No box border looks cleaner
    leg->SetHeader("Dataset Comparison", "C"); // Centered header
    leg->AddEntry(h_data, "NP02 Data (Cosmic Run 39350)", "l");
    leg->AddEntry(h_mc1,  "MC: Production Sample", "l");
    leg->AddEntry(h_mc2,  "MC: Cosmic + Nu", "l");
    leg->Draw();

    // 7. Save
    // c1->SaveAs("HitMultiplicity_DataMC_Comparison.png");
}
