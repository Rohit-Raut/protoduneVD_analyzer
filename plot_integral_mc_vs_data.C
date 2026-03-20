#include "TCanvas.h"
#include "TPad.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TRatioPlot.h"
#include <cstdlib>

void plot_integral_mc_vs_data()
{
  TFile fMC("../PDVD_Nu_MC_BDE_Reco_stage1.root","READ");
  TFile fDT("../PDVD_Run039535_Reco_gaushit.root","READ");
  
  // Check if files opened successfully
  if(fMC.IsZombie() || fDT.IsZombie()) {
    std::cout << "Error: Could not open files" << std::endl;
    return;
  }
  
  TTree* tMC = (TTree*)fMC.Get("hitdQ/Hit");
  TTree* tDT = (TTree*)fDT.Get("hitdQ/Hit");
  
  if(!tMC || !tDT) {
    std::cout << "Error: Could not find 'Hit' tree" << std::endl;
    return;
  }
  
  std::cout << "MC entries: " << tMC->GetEntries() << std::endl;
  std::cout << "Data entries: " << tDT->GetEntries() << std::endl;
  
  TH1D* hMC = new TH1D("hMC","MC",150,0,2000);
  TH1D* hDT = new TH1D("hDT","Data",150,0,2000);
  
  // Draw with selection criteria if needed (remove empty cuts if not needed)
  tMC->Draw("fIntegral/1000>>hMC","","goff");
  tDT->Draw("fIntegral/1000>>hDT","","goff");
  
  std::cout << "MC histogram entries: " << hMC->GetEntries() << std::endl;
  std::cout << "Data histogram entries: " << hDT->GetEntries() << std::endl;
  
  // Normalize
  //double mcIntegral = hMC->Integral(0, hMC->GetNbinsX()+1);
  //double dtIntegral = hDT->Integral(0, hDT->GetNbinsX()+1);
  //
  //if(mcIntegral > 0) hMC->Scale(1.0 / mcIntegral);
  //if(dtIntegral > 0) hDT->Scale(1.0 / dtIntegral);
  hMC->SetDirectory(0);
  hDT->SetDirectory(0); 
  TCanvas* c = new TCanvas("c","Integral Hit",900,650);
  gStyle->SetOptStat(0);
  // Set styling
  hDT->SetLineColor(kRed);
  hDT->SetLineWidth(2);
  hDT->SetTitle("Gaussian Hit Integral: Cosmic MC vs NP02 Run 039933");
  hDT->SetXTitle("Hit Integral [#times 10^{3}]");
  hDT->SetYTitle("Events");
  hDT->GetXaxis()->SetTitleSize(0.04);
  hDT->GetXaxis()->SetLabelSize(0.035);
  hDT->GetYaxis()->SetTitleSize(0.04);
  hDT->GetYaxis()->SetLabelSize(0.035);
  hDT->GetXaxis()->SetTitleFont(62);
  hDT->GetYaxis()->SetTitleFont(62);
  
  hMC->SetLineColor(kBlack);
  hMC->SetLineWidth(2);
  c->SetLogy();

  
  hDT->Draw("hist");
  hMC->Draw("hist same");
  c->SetGridy();

  
  
  // Add legend
  TLegend* leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  //leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  TString dtLabel = Form("Run 039933 (Mean: %.2f, Ent: %.0f)", hDT->GetMean(), hDT->GetEntries());
  TString mcLabel = Form("MC (Mean: %.2f, Ent: %.0f)", hMC->GetMean(), hMC->GetEntries());
  leg->AddEntry(hDT, dtLabel, "l");
  leg->AddEntry(hMC, mcLabel, "l");
  leg->SetTextSize(0.03);
  leg->Draw();
  gPad->Update();
  c->Draw(); 
  //c->SaveAs("hit_integral_mc_vs_data.png");
  //std::cout << "Plot saved as hit_integral_mc_vs_data.png" << std::endl;
}
