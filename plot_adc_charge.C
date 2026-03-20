#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TColor.h"
#include <iostream>

void plot_adc_charge()
{
    // 1. Open the file and access the tree
    TFile* fDT = TFile::Open("../Run039636_BDE_Reco_stage1_100.root", "READ");
    if (!fDT || fDT->IsZombie()) {
        std::cerr << "Error: Cannot open file." << std::endl;
        return;
    }

    TTree* tDT = (TTree*)fDT->Get("hitdQ/Hit");
    if (!tDT) {
        std::cerr << "Error: Cannot find tree hitdQ/Hit." << std::endl;
        return;
    }

    UInt_t  channel   = 0;
    Float_t peakTime  = 0;
    Float_t integral  = 0;
    tDT->SetBranchAddress("Channel",   &channel);
    tDT->SetBranchAddress("fPeakTime", &peakTime);
    tDT->SetBranchAddress("fIntegral", &integral);

    // 2. Global Styles - Using the High-Contrast Heatmap Palette
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kTemperatureMap); // Provides the Yellow-Red-Black heatmap look
    gStyle->SetNumberContours(255);

    // 3. Define 3x2 Grid Structure for CRP 4 & 5
    struct ViewDef {
        const char* name;
        const char* title;
        int chMin;
        int chMax;
    };

    ViewDef views[2][3] = {
        { // Row 0: CRP 4 (Bottom)
            {"hCRP4_U", "CRP4: Ind 1 (U) Charge", 0, 951},
            {"hCRP4_V", "CRP4: Ind 2 (V) Charge", 952, 1903},
            {"hCRP4_Z", "CRP4: Coll (Z) Charge",  1904, 3071}
        },
        { // Row 1: CRP 5 (Bottom)
            {"hCRP5_U", "CRP5: Ind 1 (U) Charge", 3072, 4023},
            {"hCRP5_V", "CRP5: Ind 2 (V) Charge", 4024, 4975},
            {"hCRP5_Z", "CRP5: Coll (Z) Charge",  4976, 6143}
        }
    };

    TH2D* hMaps[2][3];
    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            int nCh = views[r][c].chMax - views[r][c].chMin + 1;
            // X-axis extended to 700cm as requested previously
            hMaps[r][c] = new TH2D(views[r][c].name, "", 
                                   700, 0, 700,      
                                   nCh, views[r][c].chMin, views[r][c].chMax);
            hMaps[r][c]->SetDirectory(0);
            hMaps[r][c]->GetXaxis()->SetTitle("Drift X [cm]");
            hMaps[r][c]->GetYaxis()->SetTitle("Channel");
            hMaps[r][c]->GetZaxis()->SetTitle("Total ADC Integral");
        }
    }

    // 4. Fill Loop using fIntegral for Intensity
    Long64_t nEntries = tDT->GetEntries();
    Long64_t maxEntries = 2000000; // Increased to see clear cosmic tracks
    
    std::cout << "Filling 3x2 ADC Charge Maps..." << std::endl;

    for (Long64_t j = 0; j < std::min(nEntries, maxEntries); ++j) {
        tDT->GetEntry(j);
        double x_drift = peakTime * 0.08; // (0.5 us/tick * 0.16 cm/us)

        for (int r = 0; r < 2; ++r) {
            for (int c = 0; c < 3; ++c) {
                if (channel >= (UInt_t)views[r][c].chMin && channel <= (UInt_t)views[r][c].chMax) {
                    // Filling with integral weight
                    hMaps[r][c]->Fill(x_drift, (double)channel, (double)integral);
                    goto next_hit; 
                }
            }
        }
        next_hit:;
    }

    // 5. Canvas and Drawing (3 columns x 2 rows)
    TCanvas* c1 = new TCanvas("c1", "ADC Charge 3x2", 1800, 1000);
    c1->Divide(3, 2, 0.005, 0.005);

    for (int r = 0; r < 2; ++r) {
        for (int c = 0; c < 3; ++c) {
            int padIdx = r * 3 + c + 1; 
            TVirtualPad* pad = c1->cd(padIdx);
            pad->SetRightMargin(0.15);
            //pad->SetLogz(); 

            // Restoration of the yellowish frame background
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
    c1->SaveAs("ADC_Charge_Heatmap_3x2.png");
    std::cout << "Saved: ADC_Charge_Heatmap_3x2.png" << std::endl;
}
