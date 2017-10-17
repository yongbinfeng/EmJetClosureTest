#include "EmJetEventCount.h"

#ifndef OUTPUT
#define OUTPUT(x) std::cout<<#x << ": " << x << std::endl
#endif

using std::string;

EmJetEventCount::EmJetEventCount(EmJetSampleCollection samplesColl)
:TotalEvents_(0), n2tag_(0.), fcurrent_(-1), tweight_(1.0)
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
}

void EmJetEventCount::LoopOverCurrentTree(int ntimes)
{
  if( !IsChainValid() ){
    return ;
  }
  Long64_t nbytes = 0, nb = 0;
  timer_total_.Start();

  // prepare ntimes smeared fakerate histos for the Closure test
  PrepareFrCalVector(ntimes);
  // prepare ntimes*4 n2tags
  PrepareFrCalResults(ntimes);

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

    // Fill event-level histograms, Count number of events with n tags
    CountEvents(jentry, ntimes);
  }

  // Fill ClosureTest results
  FillEventCountHistos(ntimes);

  double total_time_elapsed = timer_total_.RealTime();
  std::cout << "Total processing time (s): " << total_time_elapsed << std::endl;
  PrintOutResults();
}


void EmJetEventCount::CountEvents(long eventnumber,  int ntimes)
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

  // Calculate number of backgrounds
  for(int itime=0; itime<ntimes; itime++){
    // initialize arrays for fakerates
    double afr0[4] = {-1.0, -1.0, -1.0, -1.0}; // results from QCD overall
    double afr1[4] = {-1.0, -1.0, -1.0, -1.0}; // results from QCD flavour 
    double afr2[4] = {-1.0, -1.0, -1.0, -1.0}; // results from GJet overall
    double afr3[4] = {-1.0, -1.0, -1.0, -1.0}; // results from GJet flavour
   
    for(int ij=0; ij<4; ij++){
      afr0[ij] = vvfrcal_[0][itime].GetFakerate(nTrack[ij]); // QCD overall
      afr2[ij] = vvfrcal_[3][itime].GetFakerate(nTrack[ij]); // GJet overall
      if (flavour[ij]<5 || flavour[ij]==21 || flavour[ij]==10 ){
        afr1[ij] = vvfrcal_[1][itime].GetFakerate(nTrack[ij]); // QCD light
        afr3[ij] = vvfrcal_[4][itime].GetFakerate(nTrack[ij]); // GJet light
      }
      else if( flavour[ij]==5 || flavour[ij]==19 ){
        afr1[ij] = vvfrcal_[2][itime].GetFakerate(nTrack[ij]); // QCD B
        afr3[ij] = vvfrcal_[5][itime].GetFakerate(nTrack[ij]); // GJet B
      }
    }

    vvn2tag_[0][itime] += PnTag(afr0, 2) * tweight_;
    vvn2tag_[1][itime] += PnTag(afr1, 2) * tweight_;
    vvn2tag_[2][itime] += PnTag(afr2, 2) * tweight_;
    vvn2tag_[3][itime] += PnTag(afr3, 2) * tweight_;
  }

  histo_->hist1d["ht"]->Fill(ht, tweight_);
  histo_->hist1d["nJet_tag"]->Fill(nJet_tag, tweight_);
  for(int ij=0; ij<4; ij++){
    histo_->hist1d["jet_pt"]->Fill((*jet_pt)[ij], tweight_);
    histo_->hist1d["jet_nTrack"]->Fill((*jet_nTrack)[ij], tweight_);
  }

  if( nJet_tag==2 ) n2tag_ += tweight_;
}

void EmJetEventCount::FillEventCountHistos(int ntimes)
{
  // fill the number of events with 2 tags
  for(int itime=0; itime < ntimes; itime++){
    histo_->hist1d["n2tag_0"]->Fill(vvn2tag_[0][itime]);
    histo_->hist1d["n2tag_1"]->Fill(vvn2tag_[1][itime]);
    histo_->hist1d["n2tag_2"]->Fill(vvn2tag_[2][itime]);
    histo_->hist1d["n2tag_3"]->Fill(vvn2tag_[3][itime]);
  }  
}

void EmJetEventCount::PrintOutResults()
{
  std::cout << "Total number of processed events is : "<< TotalEvents_ << std::endl;
  std::cout << "Total number of 2tag events observed:   " << n2tag_         << std::endl;
  std::cout << "Total number of 2tag events predicted : " << vvn2tag_[0][0] << std::endl;
  std::cout << "Total number of 2tag events predicted : " << vvn2tag_[1][0] << std::endl;
  std::cout << "Total number of 2tag events predicted : " << vvn2tag_[2][0] << std::endl;
  std::cout << "Total number of 2tag events predicted : " << vvn2tag_[3][0] << std::endl;
  std::cout << std::endl;
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

void EmJetEventCount::SetOptions(string filename, const vector<string>& vhistoname, bool isData)
{
  isData_ = isData;
  // initialize FR histograms from vhistoname
  for(auto &ihistoname: vhistoname ){
    FrCal frcal_temp = FrCal(filename, ihistoname);
    vfrcal_.push_back(frcal_temp);
  }
  std::cout << " Samples set to " << (isData? "Data": "MC") << std::endl;
  std::cout << " Fakerate histogram retrieved from" << std::endl;
  for(auto &ihistoname: vhistoname ){
    std::cout << "        " << ihistoname << std::endl;
  }
  std::cout << " of " << filename << std::endl;
}

void EmJetEventCount::InitCrossSection(const EmJetSampleCollection& samplesColl)
{
  vtreexsec_.clear();
  for(EmJetSample sample: samplesColl.samples){
    std::cout << "cross section " << sample.xsec << std::endl;
    vtreexsec_.push_back(sample.xsec); 
  }
}

double EmJetEventCount::CalculateTreeWeight(int treenumber, long eventCount)
{
  const double lumi = 20.0;
  double weight = 1.0;
  if( treenumber < static_cast<int>(vtreexsec_.size()) ){
    weight = lumi * vtreexsec_[treenumber] / eventCount;
  }
  return weight;
}

void EmJetEventCount::PrepareNewTree()
{
  std::cout << " Start a new file... " << std::endl;
  fcurrent_ = fChain->GetTreeNumber();
  TH1F* heventcount=(TH1F*)fChain->GetDirectory()->Get("eventCountperrun");
  long eventCount_current = heventcount->Integral(); // number of events in the current tree
  tweight_ = CalculateTreeWeight(fcurrent_, eventCount_current);
  TotalEvents_ += eventCount_current;
  std::cout << " tree number " << fcurrent_ << " total number of events " << eventCount_current << " tree weight " << tweight_ << " n2tag observed " << n2tag_ << std::endl;
}

void EmJetEventCount::PrepareFrCalVector(int ntimes)
{
  for( auto &tfrcal: vfrcal_ ){
    vector<FrCal> vfrcal_temp;
    for(int i=0; i< ntimes; i++){
      FrCal fr_temp = tfrcal.Clone(tfrcal.GetHistoName()+"_"+std::to_string(i)+"_"+"smeared");
      // smear the FR histograms, leaving the first one as the original
      //  print out the smeared results for the first 10 histograms
      if( i!=0 ){
        if( i<10 ) fr_temp.SmearFrHisto(true);
        else       fr_temp.SmearFrHisto(false);
      }
      vfrcal_temp.push_back(fr_temp);
    }
    vvfrcal_.push_back(vfrcal_temp);
  }
}

void EmJetEventCount::PrepareFrCalResults(int ntimes)
{ 
  for(unsigned ifrcal=0; ifrcal < 4; ifrcal++){
    vector<double> vn2tag_temp;
    for(int i=0; i<ntimes; i++){
      vn2tag_temp.push_back(0.); 
    }
    vvn2tag_.push_back(vn2tag_temp);
  }
}

bool EmJetEventCount::IsChainValid(){
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
