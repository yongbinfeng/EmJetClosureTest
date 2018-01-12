#include "EventCountBase.h"

#ifndef OUTPUT
#define OUTPUT(x) std::cout<<#x << ": " << x << std::endl
#endif

using std::string;

EventCountBase::EventCountBase()
:isData_(false), TotalEvents_(0), fcurrent_(-1), tweight_(0.0)
{
}

void EventCountBase::LoopOverTrees()
{
  if( !IsChainValid() ){
    return ;
  }
  Long64_t nbytes = 0, nb = 0;
  timer_total_.Start();

  // initialize calculations
  InitEventCount();

  // Loop over all events in TChain
  for (Long64_t jentry = 0; jentry < nentries_to_process_; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (fChain->GetTreeNumber() != fcurrent_) { // New root file opened
      PrepareNewTree();
    }
    if ( jentry % reportEvery_ == 0 ) {
      if (jentry!=0) std::cout << "Chunk processing time (s): " << timer_.RealTime() << std::endl;
      timer_.Start();
      //std::cout << "Running over global entry: " << jentry << std::endl;
    }
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    // deal with event-level variables
    LoopOverEvent(jentry);
  }

  // Analysis-level
  FinishEventCount();

  double total_time_elapsed = timer_total_.RealTime();
  std::cout << "--------------------------------------------------------------" <<  std::endl;
  std::cout << "Total number of processed events is : "<< TotalEvents_ << std::endl;
  std::cout << "Total processing time (s): " << total_time_elapsed << std::endl;
  std::cout << "--------------------------------------------------------------"  << std::endl;
}

void EventCountBase::OpenOutputFile(string ofilename)
{
  ofile_ = new TFile(ofilename.c_str(), "RECREATE");
  InitHistograms();
}

void EventCountBase::WriteHistograms()
{
  ofile_->Write();
}

void EventCountBase::SetDataType(bool isData)
{
  isData_ = isData;
  std::cout << "Set Sample Type to " << ( isData_ ? "Data" : "MC") << std::endl;
}

void EventCountBase::InitCrossSection(const vector<double>& vxsec)
{
  vtreexsec_.clear();
  for(const auto ixsec: vxsec){
    std::cout << "cross section " << ixsec << std::endl;
    vtreexsec_.push_back(ixsec); 
  }
}

double EventCountBase::CalculateTreeWeight(int treenumber, long eventCount)
{
  const double integratedLumi = 16.132;
  //const double lumi = 20.0;
  double weight = 1.0;
  if( !isData_ ){ // Normalize for Monte Carlo
    if( treenumber < static_cast<int>(vtreexsec_.size()) ){
      weight = integratedLumi * vtreexsec_[treenumber] / eventCount;
    }
  }
  return weight;
}

void EventCountBase::PrepareNewTree()
{
  std::cout << " Start a new file... " << std::endl;
  fcurrent_ = fChain->GetTreeNumber();
  TH1F* heventcount=(TH1F*)fChain->GetDirectory()->Get("eventCountperrun");
  long eventCount_current = heventcount->Integral(); // number of events in the current tree
  tweight_ = CalculateTreeWeight(fcurrent_, eventCount_current);
  TotalEvents_ += eventCount_current;

  std::cout << " tree number " << fcurrent_ << " total number of events " << eventCount_current << " tree weight " << tweight_ << std::endl;
}

bool EventCountBase::IsChainValid(){
  if (fChain == 0) {
    std::cout << "Invalid tree!" << std::endl;
    return false;
  }
  nentries_ = fChain->GetEntriesFast();
  if ( nentries_ == 0 ) {
    std::cout << "No entries!" << std::endl;
    return false;
  }
  else {
    std::cout << "Number of entries: " << nentries_ << std::endl;
  }

  // Calculate number of entries to process
  if (nentries_max_==-1) { nentries_to_process_ = nentries_                          ; }
  else                   { nentries_to_process_ = std::min(nentries_, nentries_max_) ; }
  std::cout << "Number of entries to process: " << nentries_to_process_ << std::endl;
  return true;
}
