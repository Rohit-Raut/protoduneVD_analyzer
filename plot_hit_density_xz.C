#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TColor.h"
#include <iostream>
#include <vector> // Make sure vector is included

void plot_hit_density_XZ()
{
    // 1. Open the file and access the tree
    TFile* fDT = TFile::Open("../data/Cosmic_MC_allChannel.root", "READ");
    if (!fDT || fDT->IsZombie()) {
        std::cerr << "Error: Cannot open file." << std::endl;
        return;
    }

    TTree* tDT = (TTree*)fDT->Get("hitdQ/Hit");
    if (!tDT) {
        std::cerr << "Error: Cannot find tree hitdQ/Hit." << std::endl;
        return;
    }

    // 2. IMPORTANT: Change to pointers to vectors for reading
    std::vector<UInt_t>* channel  = nullptr;
    std::vector<Float_t>* peakTime = nullptr;
    
    // Pass the address of the pointer to SetBranchAddress
    tDT->SetBranchAddress("Channel",   &channel);
    tDT->SetBranchAddress("fPeakTime", &peakTime);

    // 3. Global Styles - Standard Heatmap Palette
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kTemperatureMap); // Classic high-contrast heatmap
    gStyle->SetNumberContours(255);

    // 4. Define 3x2 Grid Structure
    struct ViewDef {
        const char* name;
        const char* title;
        int chMin;
        int chMax;
    };

    // Organized as [Row: CRP][Col: Plane]
    ViewDef views[2][3] = {
        { // Row 0: CRP 4 (Bottom)
            {"hCRP4_U", "CRP4: Induction 1 (U)", 0, 951},
            {"hCRP4_V", "CRP4: Induction 2 (V)", 952, 1903},
            {"hCRP4_Z", "CRP4: Collection (Z)",  1904, 3071}
        },
        { // Row 1: CRP 5 (Bottom)
            {"hCRP5_U", "CRP5: Induction 1 (U)", 3072, 4023},
            {"hCRP5_V", "CRP5: Induction 2 (V)", 4024, 4975},
            {"hCRP5_Z", "CRP5: Collection (Z)",  4976, 6143}
        }
    };

    TH2D* hMaps[2][3];
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            int nCh = views[r][c].chMax - views[r][c].chMin + 1;
            hMaps[r][c] = new TH2D(views[r][c].name, "",
                                   700, 0, 700,      // X drift up to 700cm
                                   nCh, views[r][c].chMin, views[r][c].chMax);
            hMaps[r][c]->SetDirectory(0);
            hMaps[r][c]->GetXaxis()->SetTitle("Drift X [cm]");
            hMaps[r][c]->GetYaxis()->SetTitle("Channel");
        }
    }

    // 5. Fill Loop
    Long64_t nEntries = tDT->GetEntries();
    Long64_t maxEntries = 1000;
    std::cout << "Filling 3x2 Heatmaps for CRP 4 & 5..." << std::endl;

    // Outer Loop: Loop over Events
    for (Long64_t j = 0; j < std::min(nEntries, maxEntries); ++j) {
        tDT->GetEntry(j);
        
        // Safety check to make sure vectors exist for this event
        if (!channel || !peakTime) continue;

        // Inner Loop: Loop over all Hits inside this Event
        size_t nHits = channel->size();
        for (size_t i = 0; i < nHits; ++i) {
            
            // Extract the specific hit data from the vectors
            UInt_t current_channel = channel->at(i);
            Float_t current_peakTime = peakTime->at(i);
            
            double x_drift = current_peakTime * 0.08;

            // Sort this hit into the correct CRP plot
            for (int r = 0; r < 2; ++r) {
                for (int c = 0; c < 3; ++c) {
                    if (current_channel >= (UInt_t)views[r][c].chMin && current_channel <= (UInt_t)views[r][c].chMax) {
                        hMaps[r][c]->Fill(x_drift, current_channel);
                        goto next_hit; // Optimization: jump to next hit once filled
                    }
                }
            }
            next_hit:; // Anchor for the goto statement
        }
    }

    // 6. Canvas and Drawing (3 columns x 2 rows)
    TCanvas* c1 = new TCanvas("c1", "Hit Density 3x2", 1800, 1000);
    c1->Divide(3, 2, 0.005, 0.005);

    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            int padIdx = r * 3 + c + 1; // 1, 2, 3 for first row; 4, 5, 6 for second
            TVirtualPad* pad = c1->cd(padIdx);
            pad->SetRightMargin(0.15);
            pad->SetLogz();

            // Clean background to make Heatmap pop
            pad->SetFrameFillColor(kWhite);

            hMaps[r][c]->Draw("COLZ");

            TLatex label;
            label.SetNDC();
            label.SetTextFont(62);
            label.SetTextSize(0.05);
            label.DrawLatex(0.12, 0.93, views[r][c].title);
        }
    }

    c1->Update();
    c1->SaveAs("Bottom_CRP4_5_Heatmap_3x2.png");
    std::cout << "Saved: Bottom_CRP4_5_Heatmap_3x2.png" << std::endl;
}
