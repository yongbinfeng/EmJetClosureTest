#ifndef EmJetFrHistos_h
#define EmJetFrHistos_h
#include <string>
#include <unordered_map>
#include <vector>
#include <cassert>
#include <TFile.h>
#include <TH1F.h>
#include <TRandom3.h>
#include "FRFormula.h"

using std::string;
using std::unordered_map;
using std::vector;

class EmJetFrHistos
{
 public:
   EmJetFrHistos(string fFRname, int ntimes, bool isData);
   ~EmJetFrHistos();
   unordered_map<string, vector<TH1F*>> histoFR1d;

 private:
   void SetDataType(bool isData);
   void SetRepeatedTime(int ntimes);
   void InitRandGen(); 
   void PrepareFrVector();
   void OpenFrRootFile(string fFRname);
   void PrepareBFractions();
   void PrepareOverallFrVector();   
   void PrepareTwoTypeVector(); 
   void PrepareTruthFrCalVector();
   
   TH1F* GetHisto(string hname);
      
   //Smear Original histogram h and return new histogram, idx is the index
   TH1F* SmearHisto(TH1F* h, int idx);
   TH1F* CloneHisto(TH1F* h, int idx);

   double SmearNumber(double val, double err, int idx);
   void PrintHistoComp(TH1F* h, TH1F* hsmeared);
 
   string DataType_;
   int ntimes_;
   TRandom* rand3_;

   // Pointer to the fakerate root file 
   TFile* ffr_;

   TH1F *hfrGJet_,  *hfrQCD_;
   TH1F *hfrGJetI_, *hfrGJetII_, *hfracGJetI_, *hfracGJetII_; 
   TH1F *hfrGJetB_, *hfrGJetL_,  *hfrQCDB_, *hfrQCDL_;

   double bfrac0_, bfrac1_;
   double err_bfrac0_, err_bfrac1_;

};

#endif
