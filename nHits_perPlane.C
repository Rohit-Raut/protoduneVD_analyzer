using Key = std::tuple<int, int, long, long, int>;

static std::map<Key, int> count_hits_perEvents_perPlane(TTree* t);
{
    std::map<Key, int> counts; 
    int run = 0, subrun = 0;
    int plane = 0;
    long long event = 0;
    t->SetBranchAddress("*", 0);
    t->SetBranchAddress("run", 1);
    t->SetBranchAddress("subrun", 1);
    t->SetBranchAddress("event", 1);
    t->SetBranchAddress("Plane", 1);

    t->SetBranchAddress("run", &run);
    t->SetBranchAddress("subrun", &subrun);
    t->SetBranchAddress("event", &event);
    t->SetBranchAddress("Plane", &event);
    const Long64_t n  = t->GetEntries();
    for (Long64_t i =  0; i<n; ++i){
	t->GetEntry(i);
	++counts[Key{run, subrun, event, plane}];
    }
    return counts;
}
void nHits_perPlane(){
    const char* fMC = "../PDVD_Nu_MC_BDE_Reco_stage1.root";
    const char* fData = "../PDVD_Run039535_Reco_gaushit.root";
    const char* treeName = "Hit";

    TFile mcFile(fMC, "READ");
    TFile dataFile(fData, "READ");
    auto* tMC   = dynamic_cast<TTree*>(mcFile.Get(treeName));
    auto* tData = dynamic_cast<TTree*>(dataFile.Get(treeName));
    if (!tMC || !tData) { printf("ERROR: cannot find TTree '%s'\n", treeName); return; }

    auto mcCounts   = count_hits_per_event_plane(tMC);
    auto dataCounts = count_hits_per_event_plane(tData);

    // Determine planes present
    std::vector<int> planes;
    planes.reserve(8);
    for (auto const& kv : dataCounts) planes.push_back(std::get<3>(kv.first));
    for (auto const& kv : mcCounts)   planes.push_back(std::get<3>(kv.first));
    std::sort(planes.begin(), planes.end());
    planes.erase(std::unique(planes.begin(), planes.end()), planes.end());

    for (int pl : planes) {
      int maxN = 0;
      for (auto const& kv : dataCounts) if (std::get<3>(kv.first) == pl) maxN = std::max(maxN, kv.second);
      for (auto const& kv : mcCounts)   if (std::get<3>(kv.first) == pl) maxN = std::max(maxN, kv.second);
      if (maxN <= 0) continue;

      auto* hData = new TH1I(Form("hData_pl%d", pl), "", maxN + 1, -0.5, maxN + 0.5);
      auto* hMC   = new TH1I(Form("hMC_pl%d",   pl), "", maxN + 1, -0.5, maxN + 0.5);

      for (auto const& kv : dataCounts) if (std::get<3>(kv.first) == pl) hData->Fill(kv.second);
      for (auto const& kv : mcCounts)   if (std::get<3>(kv.first) == pl) hMC->Fill(kv.second);

      TCanvas* c = new TCanvas(Form("c_pl%d", pl), Form("plane %d", pl), 900, 650);
      hData->SetLineWidth(2);
      hMC->SetLineWidth(2);

      hData->Draw("hist");
      hMC->Draw("hist same");

      TLegend* leg = new TLegend(0.70, 0.80, 0.88, 0.90);
      leg->AddEntry(hData, "Data", "l");
      leg->AddEntry(hMC,   "MC",   "l");
      leg->Draw();

      c->SaveAs(Form("nhits_plane%d.png", pl));
    }


}
