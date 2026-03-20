void Efficiency_data(){
    TFile* cosmic = TFile::Open("../Cosmic_MC_mergecoll_8M.root", "READ");
    TFile* nu = TFile::Open("../NU+cosmic_MC_Mergecoll_8M.root", "READ");
    if (!cosmic || !nu || cosmic->IsZombie() || nu->IsZombie()){
        std::cout<<"Error Loading Root FIle "<<std::endl;
        return;
    }
    TTree* taCos = (TTree*)cosmic->Get("hitdQ/TATree");
    TTree* tanu = (TTree*)nu->Get("hitdQ/TATree");
    if(!taCos || !tanu){
        std::cout<<"Could not find the TTree"<<std::endl;
        return;
    }
    std::vector<double> *adcVec = nullptr;
    std::vector<ULong64_t> *channelPeak = nullptr;
    std::vector<Double_t> cosmicBDE; 
    std::vector<Double_t> nuBDE;

    auto isBDE = [](int ch){
        return (ch>=0 && ch<=6143);
    };

    taCos ->SetBranchAddress("adc_integral", &adcVec);
    taCos ->SetBranchAddress("channel_peak", &channelPeak);
    for (Long64_t i = 0; i< taCos->GetEntries(); ++i){
        taCos ->GetEntry(i);
        if(!adcVec || !channelPeak)continue;
        for (size_t k = 0; k < channelPeak->size(); ++k){
            if(isBDE((int)(*channelPeak)[k])){
                cosmicBDE.push_back((*adcVec)[k]);
            }
        }

    }
    adcVec = nullptr;
    channelPeak = nullptr;
    tanu ->SetBranchAddress("adc_integral", &adcVec);
    tanu ->SetBranchAddress("channel_peak", &channelPeak);
    for (Long64_t i = 0; i< tanu->GetEntries(); ++i){
        tanu->GetEntry(i);
        if(!adcVec || !channelPeak) continue;
        for (size_t k =0; k<channelPeak->size(); ++k){
            if(isBDE((int)(*channelPeak)[k])){
                nuBDE.push_back((*adcVec)[k]);
            }
        }
    }
    //for efficiency 
    const int nBins = 200;
    const double adcMin = 1e6;
    const double adcMax = 20e6;
    TH1D* hCosmic = new TH1D("hCosmic", "Cosmic BDE", nBins, adcMin, adcMax);
    TH1D* hNu = new TH1D("hNu", "Nu BDE", nBins, adcMin, adcMax);
    for (double val: cosmicBDE) hCosmic->Fill(val);
    for (double val: nuBDE) hNu->Fill(val);
    int nTotal_cos = cosmicBDE.size();
    int nTotal_nu = nuBDE.size();

    TGraph* gEffCos = new TGraph();
    TGraph* gEffNu = new TGraph();
    for(int b =1; b<=nBins; ++b){
        double thr = hCosmic->GetBinCenter(b);

        double eff_cos = (nTotal_cos>0)? 100.0* hCosmic->Integral(b, nBins)/nTotal_cos: 0.0;
        double eff_nu = (nTotal_nu>0)? 100.0*hNu->Integral(b, nBins)/nTotal_nu: 0.0;
        gEffCos->SetPoint(b-1, thr, eff_cos);
        gEffNu ->SetPoint(b-1, thr, eff_nu);

    }

    TCanvas* c1 = new TCanvas("c1", "BDE TA Efficiency", 900, 600);
    c1->SetGrid();
    c1->SetLeftMargin(0.12);
    c1->SetBottomMargin(0.12);

    gEffCos->SetLineColor(kAzure+1);
    gEffCos->SetLineWidth(2);
    gEffCos->SetMarkerStyle(20);
    gEffCos->SetMarkerSize(0.5);
    gEffCos->SetMarkerColor(kAzure+1);

    // Style nu+cosmic
    gEffNu->SetLineColor(kPink+1);
    gEffNu->SetLineWidth(2);
    gEffNu->SetMarkerStyle(21);
    gEffNu->SetMarkerSize(0.5);
    gEffNu->SetMarkerColor(kPink+1);

    // Draw on same canvas using multigraph
    TMultiGraph* mg = new TMultiGraph();
    mg->Add(gEffCos, "LP");
    mg->Add(gEffNu,  "LP");
    mg->Draw("A");

    mg->SetTitle("ADCSimpleWindow TA Efficiency — BDE Channels (ch 0-6144);"
                 "ADC Integral Threshold [ADC counts];"
                 "Trigger Efficiency (%)");
    mg->GetXaxis()->SetLimits(adcMin, adcMax);
    mg->GetYaxis()->SetRangeUser(0, 110);

    // X axis in units of 1e6
    mg->GetXaxis()->SetMaxDigits(3);

    gPad->Update();
    TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
    leg->AddEntry(gEffCos,
                  Form("Cosmic (n=%d TAs)",     nTotal_cos), "LP");
    leg->AddEntry(gEffNu,
                  Form("Nu+Cosmic (n=%d TAs)",  nTotal_nu),  "LP");
    leg->SetBorderSize(0);
    leg->Draw();

    c1->SaveAs("efficiency_BDE_allTA.png");
    std::cout << "Saved: efficiency_BDE_allTA.png" << std::endl;


}