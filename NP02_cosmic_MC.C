void NP02_cosmic_MC(){
    gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    std::vector<TString> filePaths = {
        "../Run39350_BDE_Reco_stage1_100.root",
        "../PDVD_CosmicOnly_yk_sample_BDE_Reco_stage1.root"
    };
    std::vector<TString> fileLabels = {
        "MC Cosmic (8e6)",
        "MC Cosmic (yk)"
    };
    std::vector<int> colors = {
        kBlue,
        kBlack
    };
    std::vector<TFile*> files;
    std::vector<TH1D*> histo;
    for (int i = 0; i < 2; ++i){
        TFile* f = TFile::Open(filePaths[i], "READ");
        if(!f || f->IsZombie()){
            std::cout<<"Error: Could not open the file"<<filePaths[i]<<std::endl;
            return;
        }
        files.push_back(f);
        TTree* tree = (TTree*)f->Get("hitdQ/Hit");
        if(!tree){
            std::cout<<"Error Could not open the TTree"<<std::endl;
            return;
        }
        TString histName = Form("h_file%d", i);
        TH1D* hist = new TH1D(histName, "NP02 Hit Charge Distribution: Cosmic MC Comparison", 200, 0, 8000);
        tree->Project(histName, "fIntegral", "Channel >=0");
        hist->SetDirectory(0);
        hist->SetLineColor(colors[i]);
        hist->SetLineWidth(2);
        hist->SetXTitle("Charge Per Hit [ADC.ticks]");
        hist->SetYTitle("Normalized Events");
        if (hist->GetEntries() > 0) {
            hist->Scale(1.0 / hist->Integral());
        }
        histo.push_back(hist);
        std::cout<<fileLabels[i]<<": Entries: "<<hist->GetEntries()<<" , Mean: "<<hist->GetMean()<<std::endl;
    }
    TCanvas* c1 = new TCanvas("c1", "NP02 Cosmic MC Comparison", 1200, 800);
    gPad->SetLogy();
    gPad->SetGridy();
    for(int i = 0; i < 2; ++i){
        if(i == 0){
            histo[i]->Draw("hist");
        }else{
            histo[i]->Draw("hist same");
        }
    }
    TLegend* leg = new TLegend(0.55, 0.60, 0.95, 0.88);
    leg->SetTextSize(0.02);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    for(int i = 0; i < 2; ++i){
        leg->AddEntry(histo[i], Form("%s (N: %.4e, #mu: %.1f)", fileLabels[i].Data(), (double)histo[i]->GetEntries(), histo[i]->GetMean()), "l");
    }
    leg->Draw();
    c1->Update();
    std::cout<<"Done"<<std::endl;
}
