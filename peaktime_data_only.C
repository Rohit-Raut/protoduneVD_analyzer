void peaktime_data_only() {
    gStyle->SetOptStat(0);

    auto h_da = new TH1D("da", "", 200, 0, 10000);
    h_da->SetLineColor(kRed+1);
    h_da->SetLineWidth(2);

    auto f = TFile::Open("../PDVD_Cosmic_Run039693_BDE_Reco_stage1.root");
    TTree* t = (TTree*)f->Get("hitdQ/Hit");
    UInt_t ch; Float_t pt;
    t->SetBranchAddress("Channel",   &ch);
    t->SetBranchAddress("fPeakTime", &pt);

    for (int j = 0; j < t->GetEntries(); ++j) {
        t->GetEntry(j);
        if (ch <= 6143) h_da->Fill(pt);
    }

    auto c = new TCanvas("c", "PeakTime Data", 900, 600);
    c->SetLogy();

    h_da->GetXaxis()->SetTitle("fPeakTime [ticks]");
    h_da->GetYaxis()->SetTitle("Number of Hits");
    h_da->SetTitle("Raw fPeakTime - Data Run 039693 (BDE);fPeakTime [ticks];Number of Hits");
    h_da->Draw("HIST");

    c->SaveAs("peaktime_data_only.png");
}
