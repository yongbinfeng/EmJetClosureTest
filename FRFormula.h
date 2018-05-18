#ifndef FRFORMULA_H
#define FRFORMULA_H

#include <iostream>
#include <string>
#include <TROOT.h>
#include <TH1F.h>

namespace FRFormula{
  double PnTag(double fr[], int nTag);
  double P1tagTo2tag(double fr[]);
  double P1tagTo3tag(double fr[]);
  double P2tagTo3tag(double fr[]);
  double PEmergingnTag(double fr[], int nTag, int ijet);
  double PEmerging1tagTo2tag(double fr[], int ijet);
  TH1F* FrHistoCal(TH1F* hfrac1, TH1F* hfrac2, TH1F* hfr1, TH1F* hfr2, double bfrac, std::string tag, int idx);
  TH1F* FrHistoAdd(TH1F* hfrb, TH1F* hfrl, double bfrac, std::string tag, int idx);
  TH1F* FrHistoScale(TH1F* hfrO, TH1F* hfrNum, TH1F* hfrDen, std::string tag, int idx);
  double GetFR(TH1F* hfr, int nTrack);
  double GetRawFakerate(int nTrack, bool isBJet);
}

#endif
