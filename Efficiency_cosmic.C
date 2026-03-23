#include "utils_NP02.hpp"
#include "TCanvas.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TPad.h"
#include <iostream>
#include <vector>
#include <algorithm>


TGraph* calculateEfficiency(const std::vector<Double_t>& data, double thr_min, int nThreshold){
    double NTotal = static_cast<double>(data.size());
    std::cout<<"Debugging NTotal: "<<NTotal<<std::endl;
    double thr_max = data.back();
    //double thr_max = 40e6;
    TGraph* eff = new TGraph();
    for(int j = 0; j< nThreshold; j++){
	double thr = thr_min + j * (thr_max-thr_min) / (nThreshold -1);
	double countAbove = (double) std::distance(std::lower_bound(data.begin(), data.end(), thr), data.end());
	double efficiency = (NTotal>0)? 100.0 * (countAbove / NTotal): 0.0;
	eff->SetPoint(j, thr, efficiency);
    }
    return eff;
}

void Efficiency_cosmic(){
    TTree* taTree = NP02::getTTree("../yk_cosmic_sample.root", "hitdQ/TATree");
    if(!taTree) return;

    std::vector<double>* taADCSumVec = nullptr;
    std::vector<ULong64_t>* taChannelPeakVec = nullptr;

    taTree->SetBranchAddress("adc_integral", &taADCSumVec);
    taTree->SetBranchAddress("channel_peak", &taChannelPeakVec);

    std::vector<Double_t> cosmicBDE;
    std::cout<<"TRee Entries: "<<taTree->GetEntries()<<std::endl;
    
    for(Long64_t i = 0; i<taTree->GetEntries(); ++i){
	taTree->GetEntry(i);
	if(!taADCSumVec || !taChannelPeakVec) continue;
	int bestIdx = -1;
	double bestSum = -1.0;
	for(size_t k =0; k<taADCSumVec->size(); ++k){
	    if(NP02::isBDE((int)(*taChannelPeakVec)[k]) && (*taADCSumVec)[k]>bestSum){
		bestSum = (*taADCSumVec)[k];
		bestIdx = (int)k;
	    }
	}
	if(bestIdx >=0)
	    cosmicBDE.push_back(bestSum);
    }
    std::sort(cosmicBDE.begin(), cosmicBDE.end());
    auto countAbove = [&](double thr){
	return std::distance(std::lower_bound(cosmicBDE.begin(), cosmicBDE.end(), thr), cosmicBDE.end());
    };
    
    std::cout<<"Debug Cosmic BDE Size: "<<cosmicBDE.size()<<std::endl;


    double eff8M = 100.0 * countAbove(8e6) / (double)cosmicBDE.size();
    std::cout << "\n======== BDE Cosmic Only ========\n";
    std::cout << "Total events BDE at 8M: " << countAbove(8e6) << "\n";
    std::cout << "Efficiency at 8M ADC: " << eff8M << " %\n";
    std::cout << "=================================\n";

    const double thr_min = 1e6;
    const int nThreshold = 40;
    TGraph* effCos = calculateEfficiency(cosmicBDE, thr_min, nThreshold);
    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas("c", "BDE Efficiency-Cosmic only", 800, 600);
    effCos->SetLineColor(kBlue);
    effCos->SetLineWidth(3);
    effCos->SetMarkerStyle(20);
    effCos->SetMarkerSize(0.7);
    effCos->SetTitle("Cosmic-CORSIKA Trigger Efficiency - BDE;ADC Integral Sum Cut (ADC);Trigger Efficiency [%]");
    effCos->GetHistogram()->SetMinimum(0);
    effCos->GetHistogram()->SetMaximum(105);
    effCos->GetXaxis()->SetLimits(1e5, 12e6);
    effCos->Draw("ALP");
    c->SaveAs("Efficiency_Cosmic.png");
}
