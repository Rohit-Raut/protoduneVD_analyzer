#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include <iostream>

void plot_mc_data_percrp()
{
    // Open files
    TFile* fMC = TFile::Open("../PDVD_Nu_MC_BDE_Reco_stage1.root", "READ");
    TFile* fDT = TFile::Open("../PDVD_Run039535_Reco_gaushit.root", "READ");

    if (!fMC || fMC->IsZombie() || !fDT || fDT->IsZombie()) {
        std::cout << "Error: Could not open files" << std::endl;
        return;
    }

    TTree* tMC = (TTree*)fMC->Get("hitdQ/Hit");
    TTree* tDT = (TTree*)fDT->Get("hitdQ/Hit");

    if (!tMC || !tDT) {
        std::cout << "Error: Could not find 'Hit' tree" << std::endl;
        return;
    }

    // Configuration for CRPs
    const char* crp_names[] = {"CRP5", "CRP4"};
    const char* crp_cuts[] = {
        "Channel >= 0 && Channel <= 3071",
        "Channel >= 3072 && Channel <= 6143"
    };

    gStyle->SetOptStat(0);

    for (int i = 0; i < 2; ++i) {
        // 1. Create unique histograms (8000 scaled to 8)
        TString hNameMC = Form("hMC_%d", i);
        TString hNameDT = Form("hDT_%d", i);
        TH1D* hMC = new TH1D(hNameMC, Form("MC %s", crp_names[i]), 200, 0, 3000);
        TH1D* hDT = new TH1D(hNameDT, Form("Data %s", crp_names[i]), 200, 0, 3000);

        // Crucial: Detach from file directory so they don't disappear
        hMC->SetDirectory(0);
        hDT->SetDirectory(0);

        // 2. Fill histograms with scaling and cuts
	tMC->Draw(Form("fIntegral/1000.0>>+%s", hNameMC.Data()), crp_cuts[i], "goff");
	tDT->Draw(Form("fIntegral/1000.0>>+%s", hNameDT.Data()), crp_cuts[i], "goff");

	std::cout << crp_names[i] << " MC Entries: " << hMC->GetEntries() << std::endl;
	std::cout << crp_names[i] << " DT Entries: " << hDT->GetEntries() << std::endl;
        // 3. Normalize to Area
        double mcIntegral = hMC->Integral();
        double dtIntegral = hDT->Integral();
        if (mcIntegral > 0) hMC->Scale(1.0 / mcIntegral);
        if (dtIntegral > 0) hDT->Scale(1.0 / dtIntegral);

        // 4. Setup Canvas
        TCanvas* c = new TCanvas(Form("c_%d", i), crp_names[i], 900, 650);
        c->SetLogy();
        c->SetGridy();

        // 5. Styling
        hDT->SetLineColor(kBlack);
        hDT->SetLineWidth(2);
        hDT->SetTitle(Form("Hit Integral:Nu + Cosmic MC vs NP02 Run 039535  (%s)", crp_names[i]));
        hDT->SetXTitle("Hit Integral [#times 10^{3}]");
        hDT->SetYTitle("Normalized Events");

        hMC->SetLineColor(kRed);
        hMC->SetLineWidth(2);

        // Draw
        hDT->Draw("hist");
        hMC->Draw("hist same");

        // 6. Legend with Stats and No Border
        TLegend* leg = new TLegend(0.50, 0.72, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->AddEntry(hDT, Form("Data (Mean: %.2f, Ent: %.0f)", hDT->GetMean(), hDT->GetEntries()), "l");
        leg->AddEntry(hMC, Form("Sim (Mean: %.2f, Ent: %.0f)", hMC->GetMean(), hMC->GetEntries()), "l");
        leg->SetTextSize(0.035);
        leg->Draw();

        // 7. Save
        c->SaveAs(Form("hit_integral_%s.png", crp_names[i]));
        
        std::cout << "Generated plot for " << crp_names[i] << std::endl;
    }

    // Close files at the very end
    fMC->Close();
    fDT->Close();
}
