#include "EmJetEventCount.h"

#ifndef OUTPUT
#define OUTPUT(x) std::cout<<#x << ": " << x << std::endl
#endif

using std::string;

EmJetEventCount::EmJetEventCount(EmJetSampleCollection samplesColl)
:TotalEvents_(0), n1tag_(0.), n2tag_(0.), fcurrent_(-1), tweight_(1.0)
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

void EmJetEventCount::LoopOverTrees(int ntimes)
{
  if( !IsChainValid() ){
    return ;
  }
  Long64_t nbytes = 0, nb = 0;
  timer_total_.Start();

  // prepare ntimes smeared fakerate histos for the Closure test
  if( doPredict_ ){
    PrepareFrCalVector(ntimes);
    PrepareFrCalVector2(ntimes);
  }
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
    LoopOverEvent(jentry, ntimes);
  }

  // Fill ClosureTest results
  FillEventCountHistos(ntimes);

  double total_time_elapsed = timer_total_.RealTime();
  std::cout << "Total processing time (s): " << total_time_elapsed << std::endl;
  PrintOutResults();
}


void EmJetEventCount::LoopOverEvent(long eventnumber,  int ntimes)
{
  int nJet_tag = 0;

  for (unsigned ij = 0; ij < (*jet_pt).size(); ij++) {
    if( (*jet_isEmerging)[ij] ) nJet_tag++;
    //if( (*jet_Alpha3DSig)[ij]<0.45 && (*jet_medianIP)[ij]>0.05 ) nJet_tag++;
  }

  if( nJet_tag==1 ) n1tag_ += tweight_;
  if( nJet_tag==2 ) n2tag_ += tweight_;

  histo_->hist1d["nJet_tag"]->Fill(nJet_tag, tweight_);
 
  if( doFill_ ){ 
    FillEventHistos("");
    if( nJet_tag==0 ) FillEventHistos("__0tag");
    if( nJet_tag==1 ) FillEventHistos("__1tag");
    if( nJet_tag==2 ) FillEventHistos("__2tag");
  }

  if( !doPredict_ ) return;
  for(int itime=0; itime< ntimes; itime++){ 
    PredictBackground(itime, nJet_tag, itime==0 && doFill_ );
  }
}

void EmJetEventCount::PredictBackground(int itime, int nJet_tag, bool isfillhisto)
{
  // results from GJet overall
  double afr0[4] = {-1.0, -1.0, -1.0, -1.0};
  for(int ij=0; ij<4; ij++){
    afr0[ij] = vvfrcal_[0][itime].GetFakerate((*jet_nTrackPostCut)[ij]);
  }
  vvn1tag_[0][itime] += PnTag(afr0, 1) * tweight_;
  vvn2tag_[0][itime] += PnTag(afr0, 2) * tweight_;

  // 0tag to 1tag
  if( nJet_tag==0 ){
    double afr1[4] = {-1.0, -1.0, -1.0, -1.0}; // results from GJet 0to1tag
    double afr3[4] = {-1.0, -1.0, -1.0, -1.0}; // results from GJet 0to1tag smeared

    for( int ij=0; ij<4; ij++){
      afr1[ij] = vvfrcal_[1][itime].GetFakerate((*jet_nTrackPostCut)[ij]);
      //afr3[ij] = vvfrcal_[3][itime].GetFakerate((*jet_nTrackPostCut)[ij]);
    }
    vvn1tag_[1][itime] += PnTag(afr1, 1) * tweight_;
    //vvn1tag_[3][itime] += PnTag(afr3, 1) * tweight_;

    if( isfillhisto ) FillClosureTestHistos0To1Tag(afr1, "__GJetPredicted0To1Tag");
  }

  // 1tag to 2tag
  if( nJet_tag==1 ){
    double afr2[3] = {-1.0, -1.0, -1.0};
    double afr4[4] = {-1.0, -1.0, -1.0};
    int ijet=0;
    for(int ij=0; ij<4; ij++){
      if( (*jet_isEmerging)[ij] ) continue;// skip the emerging jet
      //if( (*jet_Alpha3DSig)[ij]<0.45 && (*jet_medianIP)[ij]>0.05 ) continue;
      afr2[ijet] = vvfrcal_[2][itime].GetFakerate((*jet_nTrackPostCut)[ij]);// 1to2tag
      //afr4[ijet] = vvfrcal_[4][itime].GetFakerate((*jet_nTrackPostCut)[ij]);// 1to2tag smeared
      ijet++;
    }
    vvn2tag_[2][itime] += P1tagTo2tag(afr2) * tweight_;
    //vvn2tag_[4][itime] += P1tagTo2tag(afr4) * tweight_;

    if( isfillhisto ) FillClosureTestHistos1To2Tag(afr2, "__GJetPredicted1To2Tag");
  }
}

void EmJetEventCount::FillEventCountHistos(int ntimes)
{
  // fill the number of events with 1 tags
  for(int itime=0; itime < ntimes; itime++){
    histo_->hist1d["n1tag__0"]->Fill(vvn1tag_[0][itime]);
    histo_->hist1d["n1tag__1"]->Fill(vvn1tag_[1][itime]);
    histo_->hist1d["n1tag__2"]->Fill(vvn1tag_[2][itime]);
    histo_->hist1d["n1tag__3"]->Fill(vvn1tag_[3][itime]);
    histo_->hist1d["n1tag__4"]->Fill(vvn1tag_[4][itime]);
  }

  // fill the number of events with 2 tags
  for(int itime=0; itime < ntimes; itime++){
    histo_->hist1d["n2tag__0"]->Fill(vvn2tag_[0][itime]);
    histo_->hist1d["n2tag__1"]->Fill(vvn2tag_[1][itime]);
    histo_->hist1d["n2tag__2"]->Fill(vvn2tag_[2][itime]);
    histo_->hist1d["n2tag__3"]->Fill(vvn2tag_[3][itime]);
    histo_->hist1d["n2tag__4"]->Fill(vvn2tag_[4][itime]);
  }  
}

void EmJetEventCount::FillClosureTestHistos0To2Tag(double fr[], string tag)
{
  double ht4=0;
  for(unsigned ij=0; ij<(*jet_pt).size(); ij++){
    ht4 += (*jet_pt)[ij];
    double pem = PEmergingnTag(fr, 2, ij);// probability of jet ij tagged as emerging
    FillJetHistos(ij, "__Predicted0To2Tag"+tag, tweight_*pem);
  }
  double prob = PnTag(fr, 2);  // probability of having 2 tags
  histo_->hist1d["ht__Predicted0To2Tag"+tag]->Fill(ht4, tweight_*prob);
}

void EmJetEventCount::FillClosureTestHistos0To1Tag(double fr[], string tag)
{
  double ht4=0;
  for(unsigned ij=0; ij<(*jet_pt).size(); ij++){
    ht4 += (*jet_pt)[ij];
    double pem = PEmergingnTag(fr, 1, ij);
    FillJetHistos(ij, "__Emerging"+tag, tweight_*pem);
  } 
  double prob = PnTag(fr, 1);
  histo_->hist1d["ht"+tag]->Fill(ht4, tweight_*prob);  
}

void EmJetEventCount::FillClosureTestHistos1To2Tag(double fr[], string tag)
{
  double ht4=0;
  double prob = P1tagTo2tag(fr);  // probability of having 2 tags
  double test = 0;
  int iemj = 0;
  for(unsigned ij=0; ij<(*jet_pt).size(); ij++){
    ht4 += (*jet_pt)[ij];
    if( (*jet_isEmerging)[ij] ){
    //if( (*jet_Alpha3DSig)[ij]<0.45 && (*jet_medianIP)[ij]>0.05 ){
      FillJetHistos(ij, "__Emerging"+tag, tweight_*prob);
     iemj = 1;
    }
    else{
      double pem = PEmerging1tagTo2tag(fr, ij-iemj);// probability of jet ij tagged as emerging
      FillJetHistos(ij, "__Emerging"+tag, tweight_*pem);
      test += pem;
    }
  }
  if( fabs(test-prob)>1e-7 ) std::cout << " inconsistent prob: " << prob << " test "<< test << std::endl;
  histo_->hist1d["ht"+tag]->Fill(ht4, tweight_*prob);
}

void EmJetEventCount::FillEventHistos(string tag, double weight)
{
  double ht4=0;
  for(unsigned ij=0; ij<(*jet_pt).size(); ij++){
    ht4 += (*jet_pt)[ij];
    FillJetFlavourHistos(ij, tag, weight);
    if( (*jet_isEmerging)[ij] ){
    //if( (*jet_Alpha3DSig)[ij]<0.45 && (*jet_medianIP)[ij]>0.05 ){
      FillJetFlavourHistos(ij, "__Emerging"+tag, weight); 
    }
    else{
      FillJetFlavourHistos(ij, "__Standard"+tag, weight);
    }
  }
  histo_->hist1d["ht"+tag]->Fill(ht4, weight);
}

void EmJetEventCount::FillEventHistos(string tag)
{
  double weight = tweight_;
  FillEventHistos(tag, weight);
}

void EmJetEventCount::FillJetFlavourHistos(int ij, string tag, double weight)
{
  FillJetHistos(ij, tag, weight);
  if( !isData_ ){
    if( (*jet_flavour)[ij]==5 || (*jet_flavour)[ij]==19 ){
      FillJetHistos(ij, "__B"+tag, weight);
    } 
    else if( (*jet_flavour)[ij]<5 || (*jet_flavour)[ij]==21 ){
      FillJetHistos(ij, "__L"+tag, weight);
    }
  }
}

void EmJetEventCount::FillJetHistos(int ij, string tag, double weight)
{
  histo_->hist1d["jet_pt"+tag]->Fill((*jet_pt)[ij],   weight);
  histo_->hist1d["jet_eta"+tag]->Fill((*jet_eta)[ij], weight);
  histo_->hist1d["jet_phi"+tag]->Fill((*jet_phi)[ij], weight);
  histo_->hist1d["jet_nTrack"+tag]->Fill((*jet_nTrack)[ij], weight);
  histo_->hist1d["jet_nTrackPostCut"+tag]->Fill((*jet_nTrackPostCut)[ij], weight);
  histo_->hist1d["jet_csv"+tag]->Fill((*jet_csv)[ij], weight);
}

void EmJetEventCount::PrintOutResults()
{
  std::cout << "Total number of processed events is : "<< TotalEvents_ << std::endl;
  std::cout << "Total number of 2tag events observed:   " << n2tag_         << std::endl;
  std::cout << "Total number of 2tag events observed:   " << histo_->hist1d["nJet_tag"]->GetBinContent(3) << "+/-" << histo_->hist1d["nJet_tag"]->GetBinError(3)<< std::endl;
  std::cout << "Total number of 2tag events predicted : "; PrintResultwithError(vvn2tag_[0]);
  std::cout << "Total number of 2tag events predicted : "; PrintResultwithError(vvn2tag_[1]);
  std::cout << "Total number of 2tag events predicted : "; PrintResultwithError(vvn2tag_[2]);
  std::cout << "Total number of 2tag events predicted : "; PrintResultwithError(vvn2tag_[3]);
  std::cout << "Total number of 2tag events predicted : "; PrintResultwithError(vvn2tag_[4]);
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Total number of 1tag events observed:  "  << n1tag_         << std::endl;
  std::cout << "Total numebr of 1tag events observed:  "  << histo_->hist1d["nJet_tag"]->GetBinContent(2)  << "+/-"<< histo_->hist1d["nJet_tag"]->GetBinError(2)<< std::endl;
  std::cout << "Total number of 1tag events predicted : "; PrintResultwithError(vvn1tag_[0]);
  std::cout << "Total number of 1tag events predicted : "; PrintResultwithError(vvn1tag_[1]);
  std::cout << "Total number of 1tag events predicted : "; PrintResultwithError(vvn1tag_[2]);
  std::cout << "Total number of 1tag events predicted : "; PrintResultwithError(vvn1tag_[3]);
  std::cout << "Total number of 1tag events predicted : "; PrintResultwithError(vvn1tag_[4]);
  std::cout << std::endl;
}

void EmJetEventCount::PrintResultwithError(const vector<double> &vresult)
{
  if( vresult.size()>=3 ){
    std::cout << " " << vresult[0] << "+" << vresult[1]-vresult[0] << "-"<< vresult[0]-vresult[2] << std::endl; 
  }
  else{
    std::cout << " " << vresult[0] << std::endl;
  } 
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

void EmJetEventCount::SetOptions(string filename, const vector<string>& vFRhnames, const vector<string>& v2fracname, const vector<string>& v2FRname, string bfractionintags, bool isData)
{
  isData_ = isData;
  std::cout << " Samples set to " << (isData? "Data": "MC") << std::endl;

  if( !doPredict_ ) return;
  // initialize FR histograms from vFRhnames
  for(auto &iname: vFRhnames ){
    FrCal frcal_temp = FrCal(filename, iname);
    vfrcal_.push_back(frcal_temp);
  }
  std::cout << " Fakerate histogram retrieved from" << std::endl;
  for(auto &iname: vFRhnames ){
    std::cout << "        " << iname << std::endl;
  }

  // initialize fraction histograms
  TFile *ff = new TFile(filename.c_str());
  hfrac1_ = (TH1F*)(ff->Get(v2fracname[0].c_str()));
  hfrac2_ = (TH1F*)(ff->Get(v2fracname[1].c_str()));
  std::cout << " 2 Fraction histograms set to \n     " << v2fracname[0] << "\n     " << v2fracname[1] << std::endl; 
  
  //initialze 2 types of Fakerates
  hfr1_ = (TH1F*)(ff->Get(v2FRname[0].c_str()));
  hfr2_ = (TH1F*)(ff->Get(v2FRname[1].c_str()));
  std::cout << " 2 Fakerate histograms set to \n     " << v2FRname[0] << "\n     "<< v2FRname[1] << std::endl;

  // get the b fractions in 0 tag and 1 tag
  TH1F* hbfrac = (TH1F*)(ff->Get(bfractionintags.c_str()));
  bfrac0_ = hbfrac->GetBinContent(1);
  bfrac1_ = hbfrac->GetBinContent(2);
  err_bfrac0_ = hbfrac->GetBinError(1);
  err_bfrac1_ = hbfrac->GetBinError(1);
  std::cout << " B fractions read from " << bfractionintags << std::endl;
  std::cout << " B fraction in 0tag " << bfrac0_ << "+/-" << err_bfrac0_ << std::endl;
  std::cout << " B fraction in 1tag " << bfrac1_ << "+/-" << err_bfrac1_ << std::endl;  

  std::cout << " of " << filename << std::endl;  
}

void EmJetEventCount::SetFillOption(bool doFill)
{
  doFill_ = doFill;
  std::cout << ( doFill_ ? "Turn on" : "Turn off" ) << " Histogram Fillings" << std::endl; 
}

void EmJetEventCount::SetPredictOption(bool doPredict)
{
  doPredict_ = doPredict;
  std::cout << ( doPredict_ ? "Turn on" : "Turn off" ) << " Background Predictions" << std::endl;  
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
  if( !isData_ ){ // Normalize for Monte Carlo
    if( treenumber < static_cast<int>(vtreexsec_.size()) ){
      weight = lumi * vtreexsec_[treenumber] / eventCount;
    }
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
  std::cout << " tree number " << fcurrent_ << " total number of events " << eventCount_current << " tree weight " << tweight_ << " n1tag observed "<< n1tag_ <<" n2tag observed " << n2tag_ << std::endl;
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
        if( i==1 )      fr_temp.SmearHistoBy1Sigma(true, 1);
        else if( i==2 ) fr_temp.SmearHistoBy1Sigma(true, -1);
        else if( i<10 ) fr_temp.SmearFrHisto(true);
        else            fr_temp.SmearFrHisto(false);
      }
      vfrcal_temp.push_back(fr_temp);
    }
    vvfrcal_.push_back(vfrcal_temp);
  }
}

void EmJetEventCount::PrepareFrCalVector2(int ntimes)
{
  vector<FrCal> vfrcal_temp0to1;
  vector<FrCal> vfrcal_temp1to2;
  bool doprint = false;
  for(int i=0; i< ntimes; i++){
    if( i< 10 ) doprint = true;
    else doprint = false;
    TH1F* hfrac1t = (TH1F*)hfrac1_->Clone((std::string(hfrac1_->GetName())+"_"+std::to_string(i)).c_str() );
    TH1F* hfrac2t = (TH1F*)hfrac2_->Clone((std::string(hfrac2_->GetName())+"_"+std::to_string(i)).c_str() );
    TH1F* hfr1t   = (TH1F*)hfr1_->Clone(  (std::string(hfr1_->GetName())  +"_"+std::to_string(i)).c_str() );
    TH1F* hfr2t   = (TH1F*)hfr2_->Clone(  (std::string(hfr2_->GetName())  +"_"+std::to_string(i)).c_str() );
    double bfrac0 = bfrac0_;
    double bfrac1 = bfrac1_;
    if( i!=0 ){
      SmearHisto(hfrac1t, doprint);
      SmearHisto(hfrac2t, doprint);
      SmearHisto(hfr1t,   doprint);
      SmearHisto(hfr2t,   doprint);
      SmearNumber(bfrac0, err_bfrac0_);
      SmearNumber(bfrac1, err_bfrac1_);
      if( doprint ) std::cout << " bfraction " << bfrac0_ << " new " << bfrac0 << std::endl;
      if( doprint ) std::cout << " bfraction " << bfrac1_ << " new " << bfrac1 << std::endl;
    }

    //TH1F* histo0to1 = FrHistoCal(hfrac1t, hfrac2t, hfr1t, hfr2t, 0.05850600927, 0.0005308097885, "0to1", (i!=0));
    //TH1F* histo1to2 = FrHistoCal(hfrac1t, hfrac2t, hfr1t, hfr2t, 0.09730214808, 0.007174446111,  "1to2", (i!=0));
    //TH1F* histo0to1 = FrHistoCal(hfrac1t, hfrac2t, hfr1t, hfr2t, 0.04638418944, 0.0003031444613, "0to1", (i!=0));
    //TH1F* histo1to2 = FrHistoCal(hfrac1t, hfrac2t, hfr1t, hfr2t, 0.07017615034, 0.004006227246,  "1to2", (i!=0));
    //TH1F* histo0to1 = FrHistoCal(hfrac1t, hfrac2t, hfr1t, hfr2t, 0.05837418267, 0.0005310752571, "0to1", (i!=0));
    //TH1F* histo1to2 = FrHistoCal(hfrac1t, hfrac2t, hfr1t, hfr2t, 0.08858562553, 0.005888970821,  "1to2", (i!=0));
    TH1F* histo0to1 = FrHistoCal(hfrac1t, hfrac2t, hfr1t, hfr2t, bfrac0, "0to1");
    TH1F* histo1to2 = FrHistoCal(hfrac1t, hfrac2t, hfr1t, hfr2t, bfrac1, "1to2");
    vfrcal_temp0to1.push_back ( FrCal(histo0to1) );
    vfrcal_temp1to2.push_back ( FrCal(histo1to2) );
  } 
  vvfrcal_.push_back(vfrcal_temp0to1);
  vvfrcal_.push_back(vfrcal_temp1to2);
}

void EmJetEventCount::PrepareFrCalResults(int ntimes)
{ 
  for(unsigned ifrcal=0; ifrcal < 10; ifrcal++){
    vector<double> vinit;
    for(int i=0; i<ntimes; i++){
      vinit.push_back(0.); 
    }
    vvn1tag_.push_back(vinit);
    vvn2tag_.push_back(vinit);
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
