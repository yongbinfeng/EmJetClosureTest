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
    long double FillHistograms(long eventnumber, double w);
    void WriteHistograms();
    void SetOptions(bool isData = false );
    //virtual int SetTree(string ifilename) = 0;
    SET_MEMBER_DEFAULT(MaxEntries, long, nentries_max, -1);
    long LoopOverCurrentTree ();

  protected:
    TFile* ofile_;
    TTree* tree_; // Current tree
    long nentries_; // Number of entries in current tree
    long nentries_max_; // Maximum number of entries to process for current tree (Set to -1 to run over all entries)
    long nentries_to_process_; // Actual number of entries to process for current tree
    TStopwatch timer_total_;

  private:
    double CalculateTreeWeight(int treenumber, long eventnumber);
    void InitCrossSection(const EmJetSampleCollection& samplesColl);
    unique_ptr<Histos> histo_;
    bool isData_;
    FrCal* frcal_;
    vector<double> vtreexsec_;
};

#endif
