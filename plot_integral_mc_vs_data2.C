#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TPad.h"
#include <iostream>

void plot_integral_mc_vs_data2()
{
    TFile* fMC = TFile::Open("../PDVD_CosmicOnly_yk_sample_BDE_Reco_stage1.root", "READ");
    //TFile* fDT = TFile::Open("../PDVD_RUN39933_Reco_gaushit_ana.root", "READ");
    //TFile* fMC = TFile::Open("../PDVD_CosmicOnly_MC_BDE_Reco_stage1.root", "READ");
    TFile* fDT = TFile::Open("../PDVD_Cosmic_Run039693_BDE_Reco_stage1.root", "READ");

    if (!fMC || fMC->IsZombie() || !fDT || fDT->IsZombie()) {
        std::cout << "Error: Check file paths." << std::endl;
        return;
    }

    TTree* tMC = (TTree*)fMC->Get("hitdQ/Hit");
    TTree* tDT = (TTree*)fDT->Get("hitdQ/Hit");

    gStyle->SetOptStat(0);
    TCanvas* c1 = new TCanvas("c1", "CRP Comparisons", 1200, 600);
    c1->Divide(2, 1);

    const char* crp_cuts[] = {
        "Channel >= 0 && Channel <= 3071",
        "Channel >= 3072 && Channel <= 6143"
    };
    auto print_stats = [&](int crp_idx, TH1* hDT, TH1* hMC) {
        const double meanDT = hDT->GetMean();
        const double rmsDT  = hDT->GetRMS();
        const double meanMC = hMC->GetMean();
        const double rmsMC  = hMC->GetRMS();
        const double ks     = hDT->KolmogorovTest(hMC, "N");

        std::cout << "CRP " << crp_idx << "\n"
                  << "  Data Run 039693: mean = " << meanDT << ", RMS = " << rmsDT << "\n"
                  << "  MC  : mean = " << meanMC << ", RMS = " << rmsMC << "\n"
		  << " K-S P-Value = " << ks	<<"\n"
                  << "----------------------------------\n";
    };
    const char* crp_labels[] = {"CRP 05", "CRP04"};
    for (int i = 0; i < 2; ++i) {
        c1->cd(i + 1);
        gPad->SetLogy();
        gPad->SetGridy();

        TString hNameMC = Form("hMC_crp%d", i);
        TString hNameDT = Form("hDT_crp%d", i);

	if(gDirectory->FindObject(hNameMC)) delete gDirectory->FindObject(hNameMC);
	if(gDirectory->FindObject(hNameDT)) delete gDirectory->FindObject(hNameDT);

        // 1. Create histograms
	TString fullTitle = Form("%s: Charge per hit (Gaussian area): Data (Run 039693) vs MC (CORSIKA)", crp_labels[i]);
        TH1D* hMC = new TH1D(hNameMC, fullTitle, 200, 0, 5000);
        TH1D* hDT = new TH1D(hNameDT, fullTitle, 200, 0, 5000);

        // 2. Fill them while they are still associated with the current directory
        //tMC->Draw(Form("fIntegral >> %s", hNameMC.Data()), crp_cuts[i], "goff");
        //tDT->Draw(Form("fIntegral >> %s", hNameDT.Data()), crp_cuts[i], "goff");
	
	tMC->Project(hNameMC, "fIntegral", crp_cuts[i]);
	tDT->Project(hNameDT, "fIntegral", crp_cuts[i]);
        // 3. Now detach them so they don't disappear when files close
        hMC->SetDirectory(0);
        hDT->SetDirectory(0);

        // 4. Check entries before scaling to avoid empty plots/errors
        double nMC = hMC->GetEntries();
        double nDT = hDT->GetEntries();
	//std::cout << "CRP " << i << " -> Data Mean: " << hDT->GetMean() << " (Entries: " << nDT << ")" << std::endl;
        //std::cout << "CRP " << i << " -> MC Mean:   " << hMC->GetMean() << " (Entries: " << nMC << ")" << std::endl;

        //if (nDT > 0) hDT->Scale(1.0 / hDT->Integral());
        //if (nMC > 0) hMC->Scale(1.0 / hMC->Integral());
	print_stats(i, hDT, hMC);
        // 5. Aesthetics
        hDT->SetLineColor(kBlack);
        hDT->SetLineWidth(2);
        hDT->SetXTitle("Charge Per Hits [ADC.ticks]");
        hDT->SetYTitle("Events");

        hMC->SetLineColor(kRed);
        hMC->SetLineWidth(2);
	gPad->SetLogy();
	//hDT->SetMinimum(1e-4);
	//double maxVal = std::max(hDT->GetMaximum(), hMC->GetMaximum());
	//hDT->SetMaximum(maxVal*10.0);
        // Adjust Y-axis range for Log scale visibility
        //double max = std::max(hDT->GetMaximum(), hMC->GetMaximum());
        //if(max > 0) hDT->GetYaxis()->SetRangeUser(1e-4, max * 10.0);

        hDT->Draw("hist");
        hMC->Draw("hist same");

        // 6. Legend
        TLegend* leg = new TLegend(0.45, 0.70, 0.88, 0.88);
        leg->SetBorderSize(0);
        leg->AddEntry(hDT, Form("Data (Mean: %.2f, Ent: %.0f)", hDT->GetMean(), nDT), "l");
        leg->AddEntry(hMC, Form("Sim (Mean: %.2f, Ent: %.0f)", hMC->GetMean(), nMC), "l");
        leg->Draw();
    }

    c1->Update();
    c1->SaveAs("CRP_Comparison_2.png");
}
