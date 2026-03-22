#ifndef NP02_UTILS_HPP
#define NP02_UTILS_HPP


#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <string>
#include <vector>

// NP02 Channel Mapping 
// 3071 channel per CRPs, and each with (U, V, Z)


namespace NP02{
    // --Plane view boundaries per crp
    constexpr int CRP4_U_MIN = 0;
    constexpr int CRP4_U_MAX = 951;
    constexpr int CRP4_V_MIN = 952;
    constexpr int CRP4_V_MAX = 1903;
    constexpr int CRP4_Z_MIN = 1904;
    constexpr int CRP4_Z_MAX = 3071;

    constexpr int CRP4_MIN = CRP4_U_MIN; // 0
    constexpr int CRP4_MAX = CRP4_Z_MAX; // 3071

    // CRP 5 (Bottom)
    constexpr int CRP5_U_MIN = 3072;
    constexpr int CRP5_U_MAX = 4023;
    constexpr int CRP5_V_MIN = 4024;
    constexpr int CRP5_V_MAX = 4975;
    constexpr int CRP5_Z_MIN = 4976;
    constexpr int CRP5_Z_MAX = 6143;
 
    constexpr int CRP5_MIN = CRP5_U_MIN;  // 3072
    constexpr int CRP5_MAX = CRP5_Z_MAX;  // 6143
 
    // CRP 2 (Top)
    constexpr int CRP2_U_MIN = 6144;
    constexpr int CRP2_U_MAX = 7095;
    constexpr int CRP2_V_MIN = 7096;
    constexpr int CRP2_V_MAX = 8047;
    constexpr int CRP2_Z_MIN = 8048;
    constexpr int CRP2_Z_MAX = 9215;
 
    constexpr int CRP2_MIN = CRP2_U_MIN;  // 6144
    constexpr int CRP2_MAX = CRP2_Z_MAX;  // 9215
 
    // CRP 3 (Top)
    constexpr int CRP3_U_MIN = 9216;
    constexpr int CRP3_U_MAX = 10167;
    constexpr int CRP3_V_MIN = 10168;
    constexpr int CRP3_V_MAX = 11119;
    constexpr int CRP3_Z_MIN = 11120;
    constexpr int CRP3_Z_MAX = 12287;
 
    constexpr int CRP3_MIN = CRP3_U_MIN;  // 9216
    constexpr int CRP3_MAX = CRP3_Z_MAX;  // 12287
					  //
					  //
					  //
    constexpr int BDE_MIN = CRP4_MIN;
    constexpr int BDE_MAX = CRP5_MAX;

    //to quick check 
    inline bool isBDE(int ch) 	{return ch>=BDE_MIN && ch<=BDE_MAX;}
    inline bool isCRP2(int ch) 	{return ch>=CRP2_MIN && ch<=CRP2_MAX;}
    inline bool isCRP3(int ch) 	{return ch>=CRP3_MIN && ch<=CRP3_MAX;}
    inline bool isCRP4(int ch) 	{return ch>=CRP4_MIN && ch<=CRP4_MAX;}
    inline bool isCRP5(int ch) 	{return ch>=CRP5_MIN && ch<=CRP5_MAX;}

    inline int getCRP(int ch){
	if(isCRP2(ch)) return 2;
	if(isCRP3(ch)) return 3;
	if(isCRP4(ch)) return 4;
	if(isCRP5(ch)) return 5;
	return -1;
    }

    //file inputout helper
    //
    inline TFile* openFile(const char* path) {
	TFile* f = TFile::Open(path, "READ"); 
	if(!f || f->IsZombie()){
	    std::cerr<<"Error: Cannot Open "<<path<<std::endl;
	    if(f){delete f;}
	    return nullptr;
	}
	return f;
    }

    inline TTree* getTTree(const char* path, const char* treeName="hitdQ/Hit"){
	TFile* f = openFile(path);
	if(!f)return nullptr;
	TTree* t = (TTree*)f->Get(treeName);
	if(!t){
	    std::cerr<<"Error Cannot fine TTree: "<<treeName<<"', in "<<path<<std::endl;
	    return nullptr; 
	}
	return t;
    }
    inline bool getTTrees(const std::vector<const char*>& paths, std::vector<TTree*>& trees, const char* treeName = "hitdQ/Hit"){
	trees.clear();
	for(auto& path: paths){
	    TTree* t = getTTree(path, treeName);
	    if(!t) {
		std::cerr<<"Failed on: "<<path<<std::endl;
		return false;
	    }
	    trees.push_back(t);
	}
	return true;
    }
} //namespace utils
#endif
