#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TColor.h"
#include <iostream>
#include <algorithm> // Needed for std::min
#include <vector>    // Needed for std::vector

void plot_hit_density_all()
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

    // 2. Set up Pointers to Vectors
    std::vector<UInt_t>* channel  = nullptr;
    std::vector<Float_t>* peakTime = nullptr;
    
    tDT->SetBranchAddress("Channel",   &channel);
    tDT->SetBranchAddress("fPeakTime", &peakTime);

    // 3. Global Styles - Standard Heatmap Palette
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kTemperatureMap); 
    gStyle->SetNumberContours(255);

    // 4. Define 4x3 Grid Structure
    struct ViewDef {
        const char* name;
        const char* title;
        int chMin;
        int chMax;
    };

    // Organized as [Row: CRP][Col: Plane] (4 Rows, 3 Columns)
    ViewDef views[4][3] = {
        { // Row 0: CRP 2
            {"hCRP2_U", "CRP2: Induction 1 (U)", 6144, 7095},
            {"hCRP2_V", "CRP2: Induction 2 (V)", 7096, 8047},
            {"hCRP2_Z", "CRP2: Collection (Z)",  8048, 9215}
        },
        { // Row 1: CRP 3
            {"hCRP3_U", "CRP3: Induction 1 (U)", 9216, 10167},
            {"hCRP3_V", "CRP3: Induction 2 (V)", 10168, 11119},
            {"hCRP3_Z", "CRP3: Collection (Z)",  11120, 12287}
        },
        { // Row 2: CRP 4 (Bottom)
            {"hCRP4_U", "CRP4: Induction 1 (U)", 0, 951},
            {"hCRP4_V", "CRP4: Induction 2 (V)", 952, 1903},
            {"hCRP4_Z", "CRP4: Collection (Z)",  1904, 3071}
        },
        { // Row 3: CRP 5 (Bottom)
            {"hCRP5_U", "CRP5: Induction 1 (U)", 3072, 4023},
            {"hCRP5_V", "CRP5: Induction 2 (V)", 4024, 4975},
            {"hCRP5_Z", "CRP5: Collection (Z)",  4976, 6143}
        }
    };

    TH2D* hMaps[4][3];
    for (int r = 0; r < 4; ++r) {
        for (int c = 0; c < 3; ++c) {
            int nCh = views[r][c].chMax - views[r][c].chMin + 1;
            hMaps[r][c] = new TH2D(views[r][c].name, "",
                                   700, 0, 700,      
                                   nCh, views[r][c].chMin, views[r][c].chMax);
            hMaps[r][c]->SetDirectory(0);
            hMaps[r][c]->GetXaxis()->SetTitle("Drift X [cm]");
            hMaps[r][c]->GetYaxis()->SetTitle("Channel");
        }
    }

    // 5. Fill Loop with Hit Density Limiter
    Long64_t nEntries = tDT->GetEntries();
    Long64_t maxHitsToPlot = 10000000; // Adjust density limit here
    Long64_t totalHitsProcessed = 0;
    
    std::cout << "Filling Heatmaps (Limit: " << maxHitsToPlot << " hits)..." << std::endl;

    for (Long64_t j = 0; j < nEntries; ++j) {
        tDT->GetEntry(j);
        
        if (!channel || !peakTime) continue;

        size_t nHits = channel->size();
        for (size_t i = 0; i < nHits; ++i) {
            
            if (totalHitsProcessed >= maxHitsToPlot) break;

            UInt_t current_channel = channel->at(i);
            Float_t current_peakTime = peakTime->at(i);
            double x_drift = current_peakTime * 0.08;

            for (int r = 0; r < 4; ++r) {
                for (int c = 0; c < 3; ++c) {
                    if (current_channel >= (UInt_t)views[r][c].chMin && current_channel <= (UInt_t)views[r][c].chMax) {
                        hMaps[r][c]->Fill(x_drift, current_channel);
                        totalHitsProcessed++; 
                        goto next_hit; 
                    }
                }
            }
            next_hit:;
        }
        
        if (totalHitsProcessed >= maxHitsToPlot) break;
    }

    // 6a. Canvas 1: CRP 2 & 3 (Rows 0 and 1)
    TCanvas* c1 = new TCanvas("c1", "CRP 2 & 3 Hit Density", 1800, 1000);
    c1->Divide(3, 2, 0.005, 0.005);

    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            int padIdx = r * 3 + c + 1; // Maps to pads 1 through 6
            TVirtualPad* pad = c1->cd(padIdx);
            pad->SetRightMargin(0.15);
            pad->SetLogz();
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

    // 6b. Canvas 2: CRP 4 & 5 (Rows 2 and 3)
    TCanvas* c2 = new TCanvas("c2", "CRP 4 & 5 Hit Density", 1800, 1000);
    c2->Divide(3, 2, 0.005, 0.005);

    for (int r = 2; r < 4; ++r) {
        for (int c = 0; c < 3; ++c) {
            int padIdx = (r - 2) * 3 + c + 1; // Shifts r=2 to pad 1, r=3 to pad 4
            TVirtualPad* pad = c2->cd(padIdx);
            pad->SetRightMargin(0.15);
            pad->SetLogz();
            pad->SetFrameFillColor(kWhite);

            hMaps[r][c]->Draw("COLZ");

            TLatex label;
            label.SetNDC();
            label.SetTextFont(62);
            label.SetTextSize(0.05); 
            label.DrawLatex(0.12, 0.93, views[r][c].title);
        }
    }
    c2->Update();

    // Uncomment these when you are ready to save the images to disk
    // c1->SaveAs("Top_CRP2_3_Heatmap_2x3.png");
    // c2->SaveAs("Bottom_CRP4_5_Heatmap_2x3.png");
    
    std::cout << "Done Processing! Two canvases generated." << std::endl;
}
