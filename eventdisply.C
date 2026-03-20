#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TColor.h"
#include <iostream>

void eventdisplay()
{
    TFile* fDT = TFile::Open("../Run039636_BDE_Reco_stage1_100.root", "READ");
    if (!fDT || fDT->IsZombie()) {
        std::cerr << "Error: Cannot open data file." << std::endl;
        return;
    }

    TTree* tDT = (TTree*)fDT->Get("hitdQ/Hit");
    if (!tDT) {
        std::cerr << "Error: Cannot find tree hitdQ/Hit." << std::endl;
        return;
    }

    UInt_t  channel  = 0;
    Float_t peakTime = 0;
    Float_t integral = 0;
    Int_t   event    = 0;

    tDT->SetBranchAddress("Channel",   &channel);
    tDT->SetBranchAddress("fPeakTime", &peakTime);
    tDT->SetBranchAddress("fIntegral", &integral);
    tDT->SetBranchAddress("event",     &event);

    tDT->GetEntry(0);
    Int_t targetEvent = event;

    gStyle->SetOptStat(0);
    gStyle->SetNumberContours(255);

    struct ViewDef {
        const char* name;
        const char* title;
        int chMin;
        int chMax;
        bool isCollection;
    };

    // Refined ranges to prevent overlap
    ViewDef views[6] = {
        {"hCRP4_U", "CRP4 - View 0 / U (Induction 1)", 0,    951,  false},
        {"hCRP4_V", "CRP4 - View 1 / V (Induction 2)", 952,  1903, false},
        {"hCRP4_Z", "CRP4 - View 2 / Z (Collection)",  1904, 3071, true },
        {"hCRP5_U", "CRP5 - View 0 / U (Induction 1)", 3072, 4023, false},
        {"hCRP5_V", "CRP5 - View 1 / V (Induction 2)", 4024, 4975, false},
        {"hCRP5_Z", "CRP5 - View 2 / Z (Collection)",  4976, 6143, true },
    };

    TH2D* hViews[6];
    for (int v = 0; v < 6; ++v) {
        int nCh = views[v].chMax - views[v].chMin + 1;
        hViews[v] = new TH2D(views[v].name, "", // Title set to empty string here
                              nCh,  views[v].chMin, views[v].chMax,
                              800,  0, 8000);
        hViews[v]->SetDirectory(0);
        hViews[v]->GetXaxis()->SetTitle("Channel Number");
        hViews[v]->GetYaxis()->SetTitle("Time Ticks (512 ns / tick)");
    }

    Long64_t nEntries = tDT->GetEntries();
    for (Long64_t j = 0; j < nEntries; ++j) {
        tDT->GetEntry(j);
        if (event != targetEvent) continue;
        for (int v = 0; v < 6; ++v) {
            if ((int)channel >= views[v].chMin && (int)channel <= views[v].chMax) {
                hViews[v]->Fill(channel, peakTime, integral);
                break;
            }
        }
    }

    TCanvas* c1 = new TCanvas("c1", "CRP4+5 Event Display", 1800, 900);
    c1->Divide(3, 2, 0.01, 0.01);

    // Define Induction Palette
    const int NRGBs = 3;
    Double_t stops[NRGBs] = {0.0, 0.5, 1.0};
    Double_t reds[NRGBs]  = {0.0, 1.0, 1.0};
    Double_t greens[NRGBs]= {0.0, 1.0, 0.0};
    Double_t blues[NRGBs] = {1.0, 1.0, 0.0};

    for (int v = 0; v < 6; ++v) {
        TVirtualPad* pad = c1->cd(v + 1);
        pad->SetRightMargin(0.15);
        pad->SetLeftMargin(0.12);
        pad->SetBottomMargin(0.12);

        //if (views[v].isCollection) {
        //    // Light green background for collection planes
        //    pad->SetFrameFillColor(kAzure-9); 
        //    gStyle->SetPalette(kViridis);
        //} else {
        //    pad->SetFrameFillColor(0); // White/Default background
        //    TColor::CreateGradientColorTable(NRGBs, stops, reds, greens, blues, 255);
        //}

        hViews[v]->Draw("COLZ");

        // Custom Bold Title
        TLatex label;
        label.SetNDC();
        label.SetTextSize(0.05);
        label.SetTextFont(62); // Bold font
        label.DrawLatex(0.15, 0.93, views[v].title);
    }

    c1->Update();
    c1->SaveAs("EventDisplay_CRP4_CRP5_Clean.png");
    std::cout << "Saved: EventDisplay_CRP4_CRP5_Clean.png" << std::endl;
}
