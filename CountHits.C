#include "TFile.h"
#include "TTree.h"
#include "TH1I.h"
#include "TCanvas.h"
#include <iostream>
#include <vector>
#include <algorithm> // Needed for std::max_element

void CountHits() {
    // 1. Open the file (Replace with your actual file path if needed)
    TFile *f = TFile::Open("../PDVD_Cosmic_nu_MC_BDE_Reco_stage1.root");
    if (!f || f->IsZombie()) {
        std::cout << "Error: Cannot open file" << std::endl;
        return;
    }

    // 2. Get the Tree
    TTree *t = (TTree*)f->Get("hitdQ/Hit");
    if (!t) {
        std::cout << "Error: Could not find TTree 'hitdQ/Hit'" << std::endl;
        f->ls(); 
        return;
    }

    // 3. Set Branch Addresses
    int run, subrun, event;
    t->SetBranchAddress("run", &run);
    t->SetBranchAddress("subrun", &subrun);
    t->SetBranchAddress("event", &event);

    // 4. Loop and Store Counts in a Vector first
    // We do this so we can find the MAX hits before making the histogram
    std::vector<int> events_hit_counts;
    
    Long64_t nentries = t->GetEntries();
    std::cout << "Processing " << nentries << " hits..." << std::endl;

    int current_event = -1;
    int current_run = -1;
    int current_subrun = -1;
    int hits_in_this_event = 0;

    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);

        // Check if we switched to a new event
        if (event != current_event || subrun != current_subrun || run != current_run) {
            
            // If this is not the first loop iteration, save the previous event's count
            if (i > 0) {
                events_hit_counts.push_back(hits_in_this_event);
            }

            // Reset counters for the NEW event
            current_event = event;
            current_subrun = subrun;
            current_run = run;
            hits_in_this_event = 1; // Start count at 1 for the current hit
        } else {
            // Same event, just increment
            hits_in_this_event++;
        }

        if (i % (nentries/10) == 0) std::cout << "." << std::flush;
    }

    // Don't forget to save the very last event after the loop finishes!
    if (nentries > 0) {
        events_hit_counts.push_back(hits_in_this_event);
    }
    std::cout << " Done!" << std::endl;

    // 5. Find Max Hits to set X-Axis dynamically
    int max_hits = 0;
    if (!events_hit_counts.empty()) {
        // This finds the largest element in our vector
        max_hits = *std::max_element(events_hit_counts.begin(), events_hit_counts.end());
    }
    
    std::cout << "Max hits found in a single event: " << max_hits << std::endl;

    // 6. Create Histogram with Dynamic Range
    // We set the max bin to (max_hits + 100) to give it a little breathing room
    TH1I *h_hits = new TH1I("h_hits", "Number of Hits per Event;Hits per Event;Count", 100, 0, max_hits + 100);

    // 7. Fill the Histogram from our vector
    for (int count : events_hit_counts) {
        h_hits->Fill(count);
    }

    // 8. Draw
    TCanvas *c1 = new TCanvas("c1", "Hits per Event", 800, 600);
    h_hits->Draw();
}
