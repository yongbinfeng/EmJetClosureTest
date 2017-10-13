#ifndef Fakerate_H
#define Fakerate_H

#include <iostream>
#include <string>
#include <TRandom3.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

const int isDebug = 1;

using std::string;

class FrCal
{
  public:
    FrCal(string filename, string histoname);
    ~FrCal() {};
    double GetFakerate(int nTrack);
    void SmearFrHisto();
  private:
    TH1F* histo_;
    string histoname_;
    int nbins_;
    int GetNBins();
};

double PnTag(double fr[], int nTag);
double PEmergingnTag(double fr[], int nTag, int ijet);
double frCal(int jet_nTrack, int option);



#endif
