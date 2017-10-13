#include "EmJetHistos.h"

using std::string;
using std::unordered_map;
using std::vector;

EmJetHistos::EmJetHistos()
{
  TH1::SetDefaultSumw2();
  hist1d["p1_0"] = new TH1F("p1_0", "p1_0" , 10, 0, 10);
  hist1d["p1_1"] = new TH1F("p1_1", "p1_1" , 10, 0, 10);
  hist1d["p1_2"] = new TH1F("p1_2", "p1_2" , 10, 0, 10);
  hist1d["p2_0"] = new TH1F("p2_0", "p2_0" , 10, 0, 10);
  hist1d["p2_1"] = new TH1F("p2_1", "p2_1" , 10, 0, 10);
  hist1d["p2_2"] = new TH1F("p2_2", "p2_2" , 10, 0, 10);
  hist1d["p3_0"] = new TH1F("p3_0", "p3_0" , 10, 0, 10);
  hist1d["p3_1"] = new TH1F("p3_1", "p3_1" , 10, 0, 10);
  hist1d["p3_2"] = new TH1F("p3_2", "p3_2" , 10, 0, 10);
  hist1d["p4_0"] = new TH1F("p4_0", "p4_0" , 10, 0, 10);
  hist1d["p4_1"] = new TH1F("p4_1", "p4_1" , 10, 0, 10);
  hist1d["p4_2"] = new TH1F("p4_2", "p4_2" , 10, 0, 10);
  hist1d["ht"]   = new TH1F("ht",   "ht",   100, 0, 2500);
  hist1d["nJet_tag"] = new TH1F("nJet_tag", "nJet_tag", 10, 0, 10);
  hist1d["n2tag"] = new TH1F("n2tag", "n2tag", 500, 6, 11);
}

EmJetHistos::~EmJetHistos()
{
  for (const auto & kv : hist1d) {
    delete kv.second;
  }
}

int
EmJetHistos::GetHistoIndex(string name)
{
  if (name=="p1_0") return 0;
  if (name=="p1_1") return 1;
  if (name=="p1_2") return 2;
  if (name=="p2_0") return 3;
  if (name=="p2_1") return 4;
  if (name=="p2_2") return 5;
  if (name=="p3_0") return 6;
  if (name=="p3_1") return 7;
  if (name=="p3_2") return 8;
  if (name=="p4_0") return 9;
  if (name=="p4_1") return 10;
  if (name=="p4_2") return 11;
  if (name=="ht")   return 12;
  if (name=="nJet_tag") return 13;
  if (name=="n2tag") return 14;
  
  return -1;
}

TH1F*
EmJetHistos::Get1D(string name)
{
  unsigned index = GetHistoIndex(name);
  assert(vector_hist1d.size()>index); // Make sure index exists first
  auto histPointer = vector_hist1d[index];
  assert(histPointer);
  return NULL;
}

void
EmJetHistos::fillToVector()
{
  vector_hist1d.reserve(hist1d.size());
  for (const auto & kv : hist1d) {
    string name = kv.first;
    auto histPointer = kv.second;
    assert( histPointer ); // Make sure pointer is not null
    assert( vector_hist1d.size() == unsigned(GetHistoIndex(name)) ); // Make sure histogram index is correct
    vector_hist1d.push_back(histPointer);
  }
}
