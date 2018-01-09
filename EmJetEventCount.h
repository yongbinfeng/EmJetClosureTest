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
#include <tuple>
#include "TLorentzVector.h"

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
    void SetOptions(string filename, const vector<string>& vFR, const vector<string>& v2frac, const vector<string>& v2FR, string sbfrac, bool isData = false);
    //virtual int SetTree(string ifilename) = 0;
    SET_MEMBER_DEFAULT(MaxEntries, long, nentries_max, -1);
    void SetFillOption(bool doFill);
    void SetPredictOption(bool doPredict);
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
    bool doPredict_;
    bool doFill_;
    vector<FrCal> vfrcal_; // vector of FrCal from input FR histograms
    vector<vector<FrCal>> vvfrcal_; // vector of smeared vector<FrCal> vfrcal_
    vector<vector<double>> vvn1tag_;
    vector<vector<double>> vvn2tag_; 
    vector<double> vtreexsec_;
    vector<int> vcase1tag_;
    vector<int> vcase2tag_;
    long TotalEvents_;
    double n1tag_;
    double n2tag_;
    int fcurrent_;     // index of the current tree in the tchain
    double tweight_;   // weight of the current tree in the tchain
    double bfrac0_; // b jet fraction in 0tag
    double bfrac1_; // b jet fraction in 1tag
    double err_bfrac0_; 
    double err_bfrac1_; 

    TH1F* hfrac1_;
    TH1F* hfrac2_;
    TH1F* hfr1_;
    TH1F* hfr2_;

    TH1F* hQfr_;
    TH1F* hQfrb_;
    TH1F* hQfrl_;
    TH1F* hGfrb_;
    TH1F* hGfrl_; 
    TH1F* hGCfrb_;
    TH1F* hGCfrl_;

    void LoopOverEvent(long eventnumber, int ntimes);
    void PredictBackground(int itime, int nJet_tag, bool isfillhisto);
    void FillEventCountHistos(int ntimes);
    void FillClosureTestHistos0To2Tag(double fr[], string tag);
    void FillClosureTestHistos0To1Tag(double fr[], string tag);
    void FillClosureTestHistos1To2Tag(double fr[], string tag);
    void FillEventHistos(string tag, double weight);
    void FillEventHistos(string tag);
    void FillJetFlavourHistos(int ij, string tag, double weight);
    void FillJetHistos(int ij, string tag, double weight);
    void FillMassHistos(int ijem1, int ijem2, string tag, double weight);
    std::pair<double, double> GetInvariantMass(int ijem1, int ijem2, int ijsm1, int ijsm2);
    vector<double> GetInvariantMass2(int ijem1, int ijem2, int ijsm1, int ijsm2);
    void PrintOutResults();
    void PrintResultwithError(const vector<double> &vresult);
    double CalculateTreeWeight(int treenumber, long eventnumber);
    void PrepareNewTree();
    void PrepareFrCalVector(int ntimes);
    void PrepareFrCalVector2(int ntimes);
    void PrepareTruthFrCalVector(int ntimes);
    void PrepareFrCalResults(int ntimes);
    void InitCrossSection(const EmJetSampleCollection& samplesColl);
    bool IsChainValid();
    void PrintOutInfo();
};

#endif
