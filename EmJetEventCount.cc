#include "EmJetEventCount.h"

#ifndef OUTPUT
#define OUTPUT(x) std::cout<<#x << ": " << x << std::endl
#endif

using std::string;

EmJetEventCount::EmJetEventCount(EmJetSampleCollection samplesColl)
{
  SetMaxEntries(-1);
  TChain* chain = new TChain("emJetSlimmedTree");
  std::cout << " EmJetEventCount::EmJetEventCount(" << samplesColl.name << ")" << std::endl;
  for (EmJetSample isample : samplesColl.samples) {
    std::cout << "Adding file to list of files to be processed: " << isample.file << std::endl;
    chain->Add(isample.file.c_str(), -1); // -1: the file is connected and the tree header read in memory to get the number of entries.
  }
  tree_ = chain;
  Init(tree_);
  InitCrossSection(samplesColl);
  std::cout << "everything okay here" << std::endl;
}


// Copied from BaseClass::Loop()
long EmJetEventCount::LoopOverCurrentTree()
{
  if (fChain == 0) {
    std::cout << "Invalid tree!" << std::endl;
    return -1;
  }
  nentries_ = fChain->GetEntriesFast();
  if ( nentries_ == 0 ) {
    std::cout << "No entries!" << std::endl;
    return -1;
  }
  else {
    std::cout << "Number of entries: " << nentries_ << std::endl;
  }

  // Calculate number of entries to process
  if (nentries_max_==-1) { nentries_to_process_ = nentries_                          ; }
  else                   { nentries_to_process_ = std::min(nentries_, nentries_max_) ; }
  std::cout << "Number of entries to process: " << nentries_to_process_ << std::endl;

  Long64_t nbytes = 0, nb = 0;
  timer_total_.Start();

  int fcurrent=-1;
  long eventCount = 0;
  double weight = 1.0;
  long double n2tag = 0.0;

  frcal_ = new FrCal("/data/users/fengyb/ClosureTest/TestClosure/FRHisto/result_fakerate.root", "fakerate_QCD");
  frcal_->SmearFrHisto();

  // Loop over all events in TChain
  for (Long64_t jentry = 0; jentry < nentries_to_process_; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (fChain->GetTreeNumber() != fcurrent) { // New root file opened
      //std::cout << " Start a new file... " << std::endl;
      fcurrent = fChain->GetTreeNumber();
      TH1F* eventcounthist=(TH1F*)fChain->GetDirectory()->Get("eventCountperrun");
      long eventCount_current = eventcounthist->Integral();
      eventCount += eventCount_current;
      weight = CalculateTreeWeight(fcurrent, eventCount_current);
      //std::cout << " tree number " << fcurrent << " total number of events " << eventCount_current << " weight " << weight << " n2tag " << n2tag << std::endl;
    }
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // Fill histograms
    n2tag += FillHistograms(jentry, weight);
  }

  histo_->hist1d["n2tag"]->Fill(n2tag);

  //std::cout << "Total number of processed events is : "<< eventCount << std::endl;
  //double total_time_elapsed = timer_total_.RealTime();
  //std::cout << "Total processing time (s): " << total_time_elapsed << std::endl;
  std::cout << "Total number of 2tag events predicted : " << n2tag << std::endl;
  std::cout << std::endl;
  return eventCount;
}

long double EmJetEventCount::FillHistograms(long eventnumber, double w)
{
  double ht = 0;
  int nJet_tag = 0;
  int nTrack[4];
  int flavour[4];

  for (unsigned ij = 0; ij < (*jet_pt).size(); ij++) {
    ht += (*jet_pt)[ij];
    nTrack[ij] = (*jet_nTrack)[ij];
    flavour[ij] = (*jet_flavour)[ij];
    if( (*jet_isEmerging)[ij] ) nJet_tag++;
  }

  double fr1[4] = {-1.0, -1.0, -1.0, -1.0};// first kind of test: without using flavour info
  //double fr2[4] = {-1.0, -1.0, -1.0, -1.0};// second kind of test: using flavour info
  //double fr3[4] = {-1.0, -1.0, -1.0, -1.0};// third kind of test: using weighted flavour fakerate
  //double fr4[4] = {-1.0, -1.0, -1.0, -1.0};// fourth kind of test
  
  for(int ij=0; ij<4; ij++){
    //fr1[ij] = frCal(nTrack[ij], 0);
    fr1[ij] = frcal_->GetFakerate(nTrack[ij]);
    //if( !isData_ ){
    //  if (flavour[ij]<5 || flavour[ij]==21 || flavour[ij]==10 ) fr2[ij] = frCal(nTrack[ij], 1);
    //  else if( flavour[ij]==5 || flavour[ij]==19 ) fr2[ij]=frCal(nTrack[ij], 2);
    //}
    //fr3[ij] = frCal(nTrack[ij], 3);
    //fr4[ij] = frCal(nTrack[ij], 4);
  }

  long double p2 = PnTag(fr1, 2);
  p2 = p2 *w ;

  //histo_->hist1d["ht"]->Fill(ht, w);
  //histo_->hist1d["nJet_tag"]->Fill(nJet_tag, w);

  return p2;
}

void EmJetEventCount::OpenOutputFile(string ofilename)
{
  ofile_ = new TFile(ofilename.c_str(), "RECREATE");
  InitHistograms();
}

void EmJetEventCount::InitHistograms()
{
  TH1::SetDefaultSumw2();
  histo_ = unique_ptr<Histos>(new Histos());
}

void EmJetEventCount::WriteHistograms()
{
  ofile_->Write();
}

void EmJetEventCount::SetOptions(bool isData)
{
  isData_ = isData;
}

void EmJetEventCount::InitCrossSection(const EmJetSampleCollection& samplesColl)
{
  vtreexsec_.clear();
  for(EmJetSample sample: samplesColl.samples){
    std::cout << "cross section " << sample.xsec*2 << std::endl;
    vtreexsec_.push_back(sample.xsec); 
  }
}

double EmJetEventCount::CalculateTreeWeight(int treenumber, long eventCount)
{
  const double lumi = 20.0;
  double weight = 1.0;
  if( treenumber < vtreexsec_.size() ){
    weight = lumi * vtreexsec_[treenumber] / eventCount;
  }
  return weight;
}
