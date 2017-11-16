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

double PnTag(double fr[], int nTag);
double P1tagTo2tag(double fr[]);
double PEmergingnTag(double fr[], int nTag, int ijet);
double PEmerging1tagTo2tag(double fr[], int ijet);
void SmearHisto(TH1F* histo, bool doprint=false, bool mustbepositive= true);
void SmearNumber(double& val, double err);
TH1F* FrHistoCal(TH1F* hfrac1, TH1F* hfrac2, TH1F* hfr1, TH1F* hfr2, double bfrac, string tag);

class FrCal
{
  public:
    FrCal();
    FrCal(string filename, string histoname);
    FrCal(TH1F* hfr);
    ~FrCal() {};
    double GetFakerate(int nTrack) const;
    void SmearFrHisto(bool doprint = true);
    void SmearHistoBy1Sigma(bool doprint, int sign);
    string GetHistoName();
    FrCal Clone(string histoname);
  private:
    TH1F* histo_;
    string histoname_;
    int nbins_;
};

#endif
