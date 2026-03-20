#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "TPad.h"
#include <iostream>
#include <map>


void NP02_hits_per_event(){
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    
    std::vector<TString> filePaths = {
        "../PDVD_CosmicOnly_yk_sample_BDE_Reco_stage1.root",
        "../PDVD_CosmicOnly_8e6_BDE_Reco_stage1.root"
    };
    
    std::vector<TString> fileLabels = {
        "Cosmic MC Production",
        "Cosmic Simulation"
    };
    
    std::vector<int> colors = {
        kRed,
        kBlue
    };
    
    std::vector<TFile*> files;
    std::vector<TH1D*> histo;
    
    for (int i = 0; i < 2; ++i){
        TFile* f = TFile::Open(filePaths[i], "READ");
        if(!f || f->IsZombie()){
            std::cout<<"Error: Could not open the file "<<filePaths[i]<<std::endl;
            return;
        }
        files.push_back(f);
        
        TTree* tree = (TTree*)f->Get("hitdQ/Hit");
        if(!tree){
            std::cout<<"Error: Could not open the TTree"<<std::endl;
            return;
        }
        
        TString histName = Form("h_file%d", i);
        TH1D* hist = new TH1D(histName, "NP02 Hits Per Event: CORSIKA MC vs R039350 Cosmic Run", 500, 0, 5000);
        hist->SetDirectory(0);
        
        // Create a map to count hits per event
        std::map<int, int> hitsPerEvent;
        UInt_t event_val;
        Int_t channel_val;
        
        tree->SetBranchAddress("event", &event_val);
        tree->SetBranchAddress("Channel", &channel_val);
        
        std::cout << "Processing " << fileLabels[i] << "..." << std::endl;
        
        // Loop through all entries and count hits per event
        for(Long64_t j = 0; j < tree->GetEntries(); ++j){
            tree->GetEntry(j);
            if(channel_val >= 0 && channel_val <= 6143){
                hitsPerEvent[event_val]++;
            }
        }
        
        // Fill histogram with the hit counts
        for(auto& pair : hitsPerEvent){
            hist->Fill(pair.second);
        }
        
        hist->SetLineColor(colors[i]);
        hist->SetLineWidth(2);
        hist->SetXTitle("Number of Hits Per Event");
        hist->SetYTitle("Number of Events");
        
        if(hist->GetEntries() > 0){
            hist->Scale(1.0 / hist->Integral());
        }
        
        histo.push_back(hist);
        std::cout << fileLabels[i] << ": Total Events: " << hist->GetEntries() 
                  << " , Mean Hits/Event: " << hist->GetMean() << std::endl;
    }
    
    double xMax = histo[0]->GetXaxis()->GetXmax();
    for (int i = 0; i < 2; ++i){
        histo[i]->GetXaxis()->SetRangeUser(0, xMax);
    }
    
    TCanvas* c1 = new TCanvas("c1", "NP02 Hits Per Event", 1200, 800);
    gPad->SetLogy();
    gPad->SetGridy();
    
    for(int i = 0; i < 2; ++i){
        if(i == 0){
            histo[i]->Draw("hist");
        }else{
            histo[i]->Draw("hist same");
        }
    }
    
    // Creating legend
    TLegend* leg = new TLegend(0.55, 0.60, 0.95, 0.88);
    leg->SetTextSize(0.02);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    
    for(int i = 0; i < 2; ++i){
        leg->AddEntry(histo[i], Form("%s (Events: %.0f, Mean: %.2f hits/event)", 
                     fileLabels[i].Data(), histo[i]->GetEntries(), histo[i]->GetMean()), "l");
    }
    leg->Draw();
    
    c1->Update();
    c1->SaveAs("NP02_HitsPerEvent.png");
    std::cout << "Done!" << std::endl;
}
