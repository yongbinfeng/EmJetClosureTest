// Base class for program to loop over single file and fill histograms with arbitrary weights
#ifndef EmJetEventCount_h
#define EmJetEventCount_h
#include "BaseClass.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>

#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>

#include "EmJetHistos.h"
#include "EmJetSample.h"

#include "Fakerate.h"

#define SET_MEMBER(NAME, TYPE, member) void Set##NAME(TYPE member) {member##_ = member;}
#define SET_MEMBER_DEFAULT(NAME, TYPE, member, DEFAULT) void Set##NAME(TYPE member = DEFAULT) {member##_ = member;}
#ifndef OUTPUT
#define OUTPUT(x) std::cout<<#x << ": " << x << std::endl
#endif

using std::string;
using std::vector;
using std::unique_ptr;

class EmJetHistos;
typedef EmJetHistos Histos;

class EmJetEventCount : protected BaseClass
{
  public:
    EmJetEventCount(EmJetSampleCollection samplesColl);
    ~EmJetEventCount(){};
    void OpenOutputFile (string ofilename);
    void InitHistograms ();
    void WriteHistograms();
    void SetOptions(string filename, const vector<string>& histoname, bool isData = false);
    //virtual int SetTree(string ifilename) = 0;
    SET_MEMBER_DEFAULT(MaxEntries, long, nentries_max, -1);
    void LoopOverTrees (int ntimes = 1);

  protected:
    TFile* ofile_;
    TTree* tree_; 
    long nentries_; // Number of entries in chain
    long nentries_max_; // Maximum number of entries to process for current tree (Set to -1 to run over all entries)
    long nentries_to_process_; // Actual number of entries to process for current tree
    TStopwatch timer_total_;
    TStopwatch timer_;
    int reportEvery_ = 100000;

  private:
    unique_ptr<Histos> histo_;
    bool isData_;
    vector<FrCal> vfrcal_; // vector of FrCal from input FR histograms
    vector<vector<FrCal>> vvfrcal_; // vector of smeared vector<FrCal> vfrcal_
    vector<vector<double>> vvn2tag_; 
    vector<double> vtreexsec_;
    long TotalEvents_;
    double n2tag_;
    int fcurrent_;     // index of the current tree in the tchain
    double tweight_;   // weight of the current tree in the tchain
    void LoopOverEvent(long eventnumber, int ntimes);
    void FillEventCountHistos(int ntimes);
    void FillClosureTestHistos0To2Tag(double fr[]);
    void FillClosureTestHistos1To2Tag(double fr[], string tag);
    void FillEventHistos(string tag, double weight);
    void FillEventHistos(string tag);
    void FillJetFlavourHistos(int ij, string tag, double weight);
    void FillJetHistos(int ij, string tag, double weight);
    void PrintOutResults();
    void PrintResultwithError(const vector<double> &vresult);
    double CalculateTreeWeight(int treenumber, long eventnumber);
    void PrepareNewTree();
    void PrepareFrCalVector(int ntimes);
    void PrepareFrCalResults(int ntimes);
    void InitCrossSection(const EmJetSampleCollection& samplesColl);
    bool IsChainValid();
};

#endif
