// Base class for program to loop over single file and fill histograms with arbitrary weights
#ifndef EmJetEventCount_h
#define EmJetEventCount_h
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

#include "EventCountBase.h"
#include "EmJetHistos.h"
#include "EmJetSample.h"
#include "EmJetFrHistos.h"
#include "FRFormula.h"


//#include "Fakerate.h"

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
class EmJetFrHistos;

class EmJetEventCount : public EventCountBase
{
  public:
    EmJetEventCount(EmJetSampleCollection samplesColl, bool doFillHisto, bool doPredict, int repeatedTime);
    ~EmJetEventCount(){};
    SET_MEMBER_DEFAULT(MaxEntries, long, nentries_max, -1);

  private:
    unique_ptr<Histos> histo_;
    unique_ptr<EmJetFrHistos> histoFR_;
    bool doPredict_;
    bool doFill_;
    int ntimes_;
    vector<vector<double>> vvn1tag_;
    vector<vector<double>> vvn2tag_; 
    vector<int> vcase1tag_;
    vector<int> vcase2tag_;
    double n1tag_;
    double n2tag_;
    string sfrfile_;

    void InitHistograms ();
    void InitEventCount ();
    void InitFrHistos(string fFRname);
    void SetFillOption(bool doFill);
    void SetPredictOption(bool doPredict);
    void LoopOverEvent(long eventnumber);
    void FinishEventCount();
    void PredictBackground(int itime, int nJet_tag, bool isfillhisto);
    void FillEventCountHistos();
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
    void PrepareNewTree();
    void PrepareFrCalResults();
    void InitCrossSection(const EmJetSampleCollection& samplesColl);
    void PrintOutInfo();
};

#endif
