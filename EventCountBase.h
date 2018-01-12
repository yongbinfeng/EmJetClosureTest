// Base class for program to loop over single file and fill histograms with arbitrary weights
#ifndef EventCountBase_h
#define EventCountBase_h
#include "BaseClass.h"
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TStopwatch.h>
#include <TH1F.h>

#include <cassert>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <tuple>

#ifndef OUTPUT
#define OUTPUT(x) std::cout<<#x << ": " << x << std::endl
#endif

class EventCountBase : protected BaseClass
{
  public:
    EventCountBase();
    ~EventCountBase(){};
    void OpenOutputFile (std::string ofilename);
    virtual void InitHistograms () = 0;
    void WriteHistograms();
    void LoopOverTrees ();
    virtual void InitEventCount() = 0;                  // to prepare before the loop
    virtual void LoopOverEvent(long eventnumber) = 0;   // event-level operation
    virtual void FinishEventCount() = 0;                // analysis-level operation

  protected:
    TFile* ofile_;
    TTree* tree_; 
    long nentries_; // Number of entries in chain
    long nentries_max_; // Maximum number of entries to process for current tree (Set to -1 to run over all entries)
    long nentries_to_process_; // Actual number of entries to process for current tree

    bool isData_;
    std::vector<double> vtreexsec_;
    long TotalEvents_;
    int fcurrent_;     // index of the current tree in the tchain
    double tweight_;   // weight of the current tree in the tchain

    void SetDataType(bool isData);
    void InitCrossSection(const vector<double>& vxsec);
    double CalculateTreeWeight(int treenumber, long eventnumber);
    virtual void PrepareNewTree();

  private:
    TStopwatch timer_total_;
    TStopwatch timer_;
    int reportEvery_ = 100000;

    bool IsChainValid();
};

#endif
