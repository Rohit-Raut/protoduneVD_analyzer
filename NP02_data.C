void NP02_data(){
	gStyle->SetOptStat(0);
    gStyle->SetLegendBorderSize(0);
    std::vector<TString> filePaths = {
	//"../PDVD_RUN39933_Reco_gaushit_ana.root",
        //"../Run039636_BDE_Reco_stage1_100.root",
        //"../PDVD_Cosmic_Run039693_BDE_Reco_stage1.root",
        //"../Run39543_BDE_Reco_stage1_100.root",
        //"../PDVD_Run039535_Reco_gaushit.root",
	//"../PDVD_CosmicOnly_8e6_BDE_Reco_stage1.root",
	"../Run39350_BDE_Reco_stage1_100.root",
	"../PDVD_CosmicOnly_yk_sample_BDE_Reco_stage1.root",
	"../PDVD_CosmicOnly_8e6_BDE_Reco_stage1.root"
    };
    std::vector<TString> fileLabels = {
	//"No Upstream: Run 039933", 
	//"No Upstream: Run 039636",  
	//"No Upstream: Run 039693", 
	//"Hadron: Run 039543",	
	//"Hadron: Run 039535",
	//"MC Cosmic Only (CORSIKA)", 
	"Cosmic Run 039350",
	"Cosmic MC Production", 
	"Cosmic Simulation"
	
    };
    std::vector<int> colors = {
	kBlack, 
	kRed, 
	kBlue, 
	//kGreen+2, 
	//kGray, 
	//kOrange,
	//kViolet,
	//kCyan
    };
    std::vector<TFile*> files;
    std::vector<TH1D*> histo;
    for (int i = 0; i<3; ++i){
	TFile* f = TFile::Open(filePaths[i], "READ");
	if(!f || f->IsZombie()){
	    std::cout<<"Error: Could not open the file"<<filePaths[i]<<std::endl;
	    return;
	}
	files.push_back(f);
    	TTree* tree = (TTree*)f->Get("hitdQ/Hit")
    	if(!tree){
    	    std::cout<<"Error Could not open the TTree"<<std::endl;
    	    return;
    	}
    	TString histName = Form("h_file%d", i);
    	TH1D* hist = new TH1D(histName, "NP02 Hit Charge Distribution: CORSIKA MC vs R039350 Cosmic Run ", 250, 0,30000);
    	tree->Project(histName, "fIntegral", "Channel>=0 && Channel <= 6143");
	//tree->Project(histName, "fIntegral", "(Channel>=1904 && Channel <= 3071) || (Channel>=4976 && Channel<=6143)");
    	hist->SetDirectory(0);
    	hist->SetLineColor(colors[i]);
    	hist->SetLineWidth(2);
    	hist->SetXTitle("Charge Per Hit [ADC.ticks]");
    	hist->SetYTitle("Events");
	//if(hist->GetEntries()>0){
	//    hist->Scale(1.0/hist->Integral());
	//}
    	histo.push_back(hist);
    	std::cout<<fileLabels[i]<<": Entries: "<<hist->GetEntries()<<" , Mean: "<<hist->GetMean()<<std::endl;
    }
    double xMax = histo[0]->GetXaxis()->GetXmax();
    for (int i = 0; i<3; ++i){
	histo[i]->GetXaxis()->SetRangeUser(0, xMax);
    }
    TCanvas* c1 = new TCanvas("c1", "NP02 Run Gaussian hit", 1200, 800);
    gPad->SetLogy();
    gPad->SetGridy();
    for(int i =0; i<3; ++i){
	if(i==0){
	    histo[i]->Draw("hist");
	}else{
	    histo[i]->Draw("hist same");
	}
    }

    //creating legend 
    TLegend* leg =  new TLegend(0.55, 0.60, 0.95, 0.88);
    leg->SetTextSize(0.02);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);

    for(int i = 0; i<3; ++i){
	leg->AddEntry(histo[i], Form("%s (N: %.4e, #mu: %.1f)", fileLabels[i].Data(), (double)histo[i]->GetEntries(), histo[i]->GetMean()), "l");
    }
    leg->Draw();
    c1->Update();
    std::cout<<"Done"<<std::endl;
    

}
