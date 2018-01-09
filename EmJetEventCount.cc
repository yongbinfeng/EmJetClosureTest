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

  // prepare ntimes*4 n2tags
  PrepareFrCalResults(ntimes);

  // prepare ntimes smeared fakerate histos for the Closure test
  if( doPredict_ ){
    PrepareFrCalVector(ntimes);
    PrepareFrCalVector2(ntimes);
    PrepareTruthFrCalVector(ntimes);
  }

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
  if( pvtrack_fraction<0.1 ) return;
  int nJet_tag = 0;
  //if( (*jet_pt)[3]<100.0 ) return;
  //double ht4 = (*jet_pt)[0] + (*jet_pt)[1] + (*jet_pt)[2] + (*jet_pt)[3];
  //if( ht4< 1500.0 ) return;
  //int nJet_Alpha3DSig_zero=0;

  for (unsigned ij = 0; ij < (*jet_pt).size(); ij++) {
    if( (*jet_isEmerging)[ij] ) nJet_tag++;
    //if( (*jet_Alpha3DSig)[ij]<0.45 && (*jet_medianIP)[ij]>0.05 ) nJet_tag++;
    //if( (*jet_Alpha3DSig)[ij]<1e-6 ) nJet_Alpha3DSig_zero++;
  }

  //if( nJet_Alpha3DSig_zero>2 ) return;
  histo_->hist1d["pv_genreco_disXY"]->Fill(pv_genreco_disXY, tweight_);
  histo_->hist1d["pv_genreco_disZ" ]->Fill(fabs(pv_genreco_disZ), tweight_);

  if( nJet_tag==1 ) {
    n1tag_ += tweight_;
    histo_->hist1d["pv_genreco_disXY_1tag"]->Fill(pv_genreco_disXY, tweight_);
    histo_->hist1d["pv_genreco_disZ_1tag" ]->Fill(fabs(pv_genreco_disZ), tweight_);
  }
  if( nJet_tag==2 ){
    n2tag_ += tweight_;
    //PrintOutInfo();
    histo_->hist1d["pv_genreco_disXY_2tag"]->Fill(pv_genreco_disXY, tweight_);
    histo_->hist1d["pv_genreco_disZ_2tag" ]->Fill(fabs(pv_genreco_disZ), tweight_);
  }

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

  // from 0tag
  if( nJet_tag==0 ){
    double afr1[4] = {-1.0, -1.0, -1.0, -1.0}; // results from GJet 0to1tag
    double afr3[4] = {-1.0, -1.0, -1.0, -1.0}; // results from QCD truth 0to1tag
    double afr5[4] = {-1.0, -1.0, -1.0, -1.0}; // results from GJet truth 0to1tag 
    double afr7[4] = {-1.0, -1.0, -1.0, -1.0}; // results from GJet calc 0to1tag
    //double afr9[4] = {-1.0, -1.0, -1.0, -1.0}; // results from QCD using truth flavor

    for( int ij=0; ij<4; ij++){
      afr1[ij] = vvfrcal_[1][itime].GetFakerate((*jet_nTrackPostCut)[ij]);
      afr3[ij] = vvfrcal_[3][itime].GetFakerate((*jet_nTrackPostCut)[ij]);
      afr5[ij] = vvfrcal_[5][itime].GetFakerate((*jet_nTrackPostCut)[ij]);
      afr7[ij] = vvfrcal_[7][itime].GetFakerate((*jet_nTrackPostCut)[ij]);
      /*
      if( (*jet_flavour)[ij]==5 || (*jet_flavour)[ij]==19 || (*jet_flavour)[ij]==10 ) {
        afr9[ij] = vvfrcal_[9][itime].GetFakerate((*jet_nTrackPostCut)[ij]);
        //afr12[ij] = GetRawFakerate((*jet_nTrackPostCut)[ij], true);
      }
      else{
        afr9[ij] = vvfrcal_[10][itime].GetFakerate((*jet_nTrackPostCut)[ij]);
        //afr12[ij] = GetRawFakerate((*jet_nTrackPostCut)[ij], false);
      }
      */
    }
    vvn1tag_[1][itime] += PnTag(afr1, 1) * tweight_;
    vvn1tag_[2][itime] += PnTag(afr3, 1) * tweight_;
    vvn1tag_[3][itime] += PnTag(afr5, 1) * tweight_;
    vvn1tag_[4][itime] += PnTag(afr7, 1) * tweight_;
    //vvn1tag_[5][itime] += PnTag(afr9, 1) * tweight_;

    //vvn2tag_[5][itime] += PnTag(afr9, 2) * tweight_;

    if( isfillhisto ) {
      FillClosureTestHistos0To1Tag(afr0, "__GJetOverallPredicted0To1Tag");
      FillClosureTestHistos0To1Tag(afr1, "__GJetCalcPredicted0To1Tag");
      FillClosureTestHistos0To1Tag(afr3, "__QCDTruthPredicted0To1Tag");
      FillClosureTestHistos0To1Tag(afr5, "__GJetTruthPredicted0To1Tag");
    }
  }

  // 1tag to 2tag
  if( nJet_tag==1 ){
    double afr2[3] = {-1.0, -1.0, -1.0};
    double afr4[3] = {-1.0, -1.0, -1.0};
    double afr6[3] = {-1.0, -1.0, -1.0}; 
    double afr8[3] = {-1.0, -1.0, -1.0};
    int ijet=0;
    for(int ij=0; ij<4; ij++){
      if( (*jet_isEmerging)[ij] ) continue;// skip the emerging jet
      //if( (*jet_Alpha3DSig)[ij]<0.45 && (*jet_medianIP)[ij]>0.05 ) continue;
      afr2[ijet] = vvfrcal_[2][itime].GetFakerate((*jet_nTrackPostCut)[ij]);// 1to2tag
      afr4[ijet] = vvfrcal_[4][itime].GetFakerate((*jet_nTrackPostCut)[ij]);// 1to2tag QCD truth
      afr6[ijet] = vvfrcal_[6][itime].GetFakerate((*jet_nTrackPostCut)[ij]);// 1to2tag GJet truth
      afr8[ijet] = vvfrcal_[8][itime].GetFakerate((*jet_nTrackPostCut)[ij]);// 1to2tag GJet truth
      ijet++;
    }
    vvn2tag_[1][itime] += P1tagTo2tag(afr2) * tweight_;
    vvn2tag_[2][itime] += P1tagTo2tag(afr4) * tweight_;
    vvn2tag_[3][itime] += P1tagTo2tag(afr6) * tweight_;
    vvn2tag_[4][itime] += P1tagTo2tag(afr8) * tweight_;

    if( isfillhisto ) {
      FillClosureTestHistos1To2Tag(afr2, "__GJetCalcPredicted1To2Tag");
      FillClosureTestHistos1To2Tag(afr4, "__QCDTruthPredicted1To2Tag");
      FillClosureTestHistos1To2Tag(afr6, "__GJetTruthPredicted1To2Tag");
    }
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
    if( vcase1tag_[itime]==1 )      histo_->hist1d["N1tag_case1__1"]->Fill(vvn1tag_[1][itime]);
    else if( vcase1tag_[itime]==2 ) histo_->hist1d["N1tag_case2__1"]->Fill(vvn1tag_[1][itime]);
  }

  for(int itime=0; itime < ntimes; itime++){
    histo_->hist1d["n1tag_EM1__0"]->Fill(vvn1tag_[0][itime]);
    histo_->hist1d["n1tag_EM1__1"]->Fill(vvn1tag_[1][itime]);
    histo_->hist1d["n1tag_EM1__2"]->Fill(vvn1tag_[2][itime]);
    histo_->hist1d["n1tag_EM1__3"]->Fill(vvn1tag_[3][itime]);
    histo_->hist1d["n1tag_EM1__4"]->Fill(vvn1tag_[4][itime]);
  }

  // fill the number of events with 2 tags
  for(int itime=0; itime < ntimes; itime++){
    histo_->hist1d["n2tag__0"]->Fill(vvn2tag_[0][itime]);
    histo_->hist1d["n2tag__1"]->Fill(vvn2tag_[1][itime]);
    histo_->hist1d["n2tag__2"]->Fill(vvn2tag_[2][itime]);
    histo_->hist1d["n2tag__3"]->Fill(vvn2tag_[3][itime]);
    histo_->hist1d["n2tag__4"]->Fill(vvn2tag_[4][itime]);
    if( vcase2tag_[itime]==1 )      histo_->hist1d["N2tag_case1__1"]->Fill(vvn1tag_[1][itime]);
    else if( vcase2tag_[itime]==2 ) histo_->hist1d["N2tag_case2__1"]->Fill(vvn1tag_[1][itime]);
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
  int ijem1 = 0;
  for(unsigned ij=0; ij<(*jet_pt).size(); ij++){
    if( (*jet_isEmerging)[ij] ) ijem1 = ij;
    //if( (*jet_Alpha3DSig)[ij]<0.45 && (*jet_medianIP)[ij]>0.05 ) ijem1 = ij;
  }
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
      FillMassHistos(ijem1, ij, tag, tweight_*pem);
      //auto massPair = GetInvariantMass(ijem1, ij);
      //histo_->hist1d["mass"+tag]->Fill(massPair.first, tweight_*pem);
      //histo_->hist1d["mass"+tag]->Fill(massPair.second, tweight_*pem);
      test += pem;
    }
  }
  if( fabs(test-prob)>1e-7 ) std::cout << " inconsistent prob: " << prob << " test "<< test << std::endl;
  histo_->hist1d["ht"+tag]->Fill(ht4, tweight_*prob);
}

void EmJetEventCount::FillEventHistos(string tag, double weight)
{
  double ht4=0;
  int ijem1 = -1, ijem2 = -1;
  for(unsigned ij=0; ij<(*jet_pt).size(); ij++){
    ht4 += (*jet_pt)[ij];
    FillJetFlavourHistos(ij, tag, weight);
    if( (*jet_isEmerging)[ij] ){
    //if( (*jet_Alpha3DSig)[ij]<0.45 && (*jet_medianIP)[ij]>0.05 ){
      FillJetFlavourHistos(ij, "__Emerging"+tag, weight); 
      if( ijem1!=-1 ) ijem2 = ij;
      else ijem1 = ij;
    }
    else{
      FillJetFlavourHistos(ij, "__Standard"+tag, weight);
    }
  }
  histo_->hist1d["ht"+tag]->Fill(ht4, weight);
  if( ijem1!= -1 && ijem2!=-1 ){
    FillMassHistos(ijem1, ijem2, tag, weight);
  }
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
  //if( (*jet_pt)[ij]<=100.0 )                         histo_->hist1d["jet_nTrackPostCut_pt1"+tag]->Fill((*jet_nTrackPostCut)[ij], weight);
  //if( (*jet_pt)[ij]>100.0 && (*jet_pt)[ij]<=200.0 )  histo_->hist1d["jet_nTrackPostCut_pt2"+tag]->Fill((*jet_nTrackPostCut)[ij], weight);
  //if( (*jet_pt)[ij]>200.0 )                          histo_->hist1d["jet_nTrackPostCut_pt3"+tag]->Fill((*jet_nTrackPostCut)[ij], weight);
  histo_->hist1d["jet_csv"+tag]->Fill((*jet_csv)[ij], weight);
}

void EmJetEventCount::FillMassHistos(int ijem1, int ijem2, string tag, double weight)
{
  int ijsm1 = -1, ijsm2 = -1;
  for(int ij=0; ij<4; ij++){
    if( ij!=ijem1 && ij!=ijem2 ){
      if( ijsm1!=-1 ) ijsm2=ij;
      else ijsm1=ij;
    }
  }
  //std::cout << " emjet " << ijem1 << " " << ijem2 << std::endl;
  //std::cout << " smjet " << ijsm1 << " " << ijsm2 << std::endl;
  auto massPair = GetInvariantMass(ijem1, ijem2, ijsm1, ijsm2);
  histo_->hist1d["mass"+tag]->Fill(massPair.first, weight);
  histo_->hist1d["mass"+tag]->Fill(massPair.second, weight);
  auto vmass = GetInvariantMass2(ijem1, ijem2, ijsm1, ijsm2);
  histo_->hist1d["mass1"+tag]->Fill(vmass[0], weight);
  histo_->hist1d["mass1"+tag]->Fill(vmass[1], weight);
  histo_->hist1d["mass1"+tag]->Fill(vmass[2], weight);
  histo_->hist1d["mass1"+tag]->Fill(vmass[3], weight);
}

std::pair<double, double> EmJetEventCount::GetInvariantMass(int ijem1, int ijem2, int ijsm1, int ijsm2)
{
  TLorentzVector jetemVector1, jetemVector2, jetsmVector1, jetsmVector2;
  jetemVector1.SetPtEtaPhiM((*jet_pt)[ijem1], (*jet_eta)[ijem1], (*jet_phi)[ijem1], 0.);
  jetemVector2.SetPtEtaPhiM((*jet_pt)[ijem2], (*jet_eta)[ijem2], (*jet_phi)[ijem2], 0.);
  jetsmVector1.SetPtEtaPhiM((*jet_pt)[ijsm1], (*jet_eta)[ijsm1], (*jet_phi)[ijsm1], 0.);
  jetsmVector2.SetPtEtaPhiM((*jet_pt)[ijsm2], (*jet_eta)[ijsm2], (*jet_phi)[ijsm2], 0.);
  double m11 = (jetemVector1+jetsmVector1).Mag();
  double m22 = (jetemVector2+jetsmVector2).Mag();
  double m12 = (jetemVector1+jetsmVector2).Mag();
  double m21 = (jetemVector2+jetsmVector1).Mag();
  if( fabs(m11-m22) < fabs(m12-m21) ) return std::make_pair(m11, m22);
  else return std::make_pair(m12, m21);
}

vector<double> EmJetEventCount::GetInvariantMass2(int ijem1, int ijem2, int ijsm1, int ijsm2)
{
  TLorentzVector jetemVector1, jetemVector2, jetsmVector1, jetsmVector2;
  jetemVector1.SetPtEtaPhiM((*jet_pt)[ijem1], (*jet_eta)[ijem1], (*jet_phi)[ijem1], 0.);
  jetemVector2.SetPtEtaPhiM((*jet_pt)[ijem2], (*jet_eta)[ijem2], (*jet_phi)[ijem2], 0.);
  jetsmVector1.SetPtEtaPhiM((*jet_pt)[ijsm1], (*jet_eta)[ijsm1], (*jet_phi)[ijsm1], 0.);
  jetsmVector2.SetPtEtaPhiM((*jet_pt)[ijsm2], (*jet_eta)[ijsm2], (*jet_phi)[ijsm2], 0.);
  double m11 = (jetemVector1+jetsmVector1).Mag();
  double m22 = (jetemVector2+jetsmVector2).Mag();
  double m12 = (jetemVector1+jetsmVector2).Mag();
  double m21 = (jetemVector2+jetsmVector1).Mag();
  vector<double> vmass = {m11, m22, m12, m21};
  return vmass;
}

void EmJetEventCount::PrintOutResults()
{
  std::cout << "Total number of processed events is : "<< TotalEvents_ << std::endl;
  std::cout << "Total number of 2tag events observed:   " << n2tag_         << std::endl;
  std::cout << "Total number of 2tag events observed:   " << histo_->hist1d["nJet_tag"]->GetBinContent(3) << "+/-" << histo_->hist1d["nJet_tag"]->GetBinError(3)<< std::endl;
  std::cout << "Total number of 2tag events predicted(GJet overall) : "; PrintResultwithError(vvn2tag_[0]);
  std::cout << "Total number of 2tag events predicted(GJet Calc)    : "; PrintResultwithError(vvn2tag_[1]);
  std::cout << "Total number of 2tag events predicted(QCD truth)    : "; PrintResultwithError(vvn2tag_[2]);
  std::cout << "Total number of 2tag events predicted(GJet truth)   : "; PrintResultwithError(vvn2tag_[3]);
  std::cout << "Total number of 2tag events predicted(GJet calc)    : "; PrintResultwithError(vvn2tag_[4]);
  std::cout << "Total number of 2tag events predicted(QCD truth tr) : "; PrintResultwithError(vvn2tag_[5]);
  std::cout << "Total number of 2tag events predicted(GJet truth tr): "; PrintResultwithError(vvn2tag_[6]);
  std::cout << "Total number of 2tag events predicted(QCD truth tr pt): "; PrintResultwithError(vvn2tag_[7]);
  std::cout << "Total number of 2tag events predicted(QCD raw input): "; PrintResultwithError(vvn2tag_[8]);
  std::cout << "------------------------------------------------" << std::endl;
  std::cout << "Total number of 1tag events observed:  "  << n1tag_         << std::endl;
  std::cout << "Total numebr of 1tag events observed:  "  << histo_->hist1d["nJet_tag"]->GetBinContent(2)  << "+/-"<< histo_->hist1d["nJet_tag"]->GetBinError(2)<< std::endl;
  std::cout << "Total number of 1tag events predicted : "; PrintResultwithError(vvn1tag_[0]);
  std::cout << "Total number of 1tag events predicted : "; PrintResultwithError(vvn1tag_[1]);
  std::cout << "Total number of 1tag events predicted : "; PrintResultwithError(vvn1tag_[2]);
  std::cout << "Total number of 1tag events predicted : "; PrintResultwithError(vvn1tag_[3]);
  std::cout << "Total number of 1tag events predicted : "; PrintResultwithError(vvn1tag_[4]);
  std::cout << "Total number of 1tag events predicted : "; PrintResultwithError(vvn1tag_[5]);
  std::cout << "Total number of 1tag events predicted : "; PrintResultwithError(vvn1tag_[6]);
  std::cout << "Total number of 1tag events predicted : "; PrintResultwithError(vvn1tag_[7]);
  std::cout << "Total number of 1tag events predicted : "; PrintResultwithError(vvn1tag_[8]);
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

  //initialize QCD truth fakerate
  //TFile *ff2 = new TFile("/data/users/fengyb/80Xresult/plot/fakerate_QCD_3DSig_8_Med_8_remove3_pt100.root");
  hQfr_  = (TH1F*)(ff->Get("Dfakerates/fakerate_QCDMC"));
  hQfrb_ = (TH1F*)(ff->Get("Dfakerates/fakerate_QCDMC_BJet"));
  hQfrl_ = (TH1F*)(ff->Get("Dfakerates/fakerate_QCDMC_LightJet"));
  //hQfrb_ = (TH1F*)(ff2->Get("fakerate_QCD_B"));
  //hQfrl_ = (TH1F*)(ff2->Get("fakerate_QCD_L"));
  std::cout << " Truth QCD Fakerate histograms set " << std::endl;

  //initialize GJet truth fakerate
  hGfrb_ = (TH1F*)(ff->Get("Dfakerates/fakerate_GJetMC_BJet"));
  hGfrl_ = (TH1F*)(ff->Get("Dfakerates/fakerate_GJetMC_LightJet"));
  std::cout << " Truth GJet Fakerate histograms set " << std::endl;

  //initialize GJet Calc fakerate
  if( !isData_ ){
    hGCfrb_ = (TH1F*)(ff->Get("fakerate_GJetMC_BJet_calc"));
    hGCfrl_ = (TH1F*)(ff->Get("fakerate_GJetMC_LightJet_calc"));
  }
  else{
    hGCfrb_ = (TH1F*)(ff->Get("fakerate_GJetData_BJet_calc"));
    hGCfrl_ = (TH1F*)(ff->Get("fakerate_GJetData_LightJet_calc"));
  }
  std::cout << " Calc GJet Fakerate histograms set " << std::endl;

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
  const double lumi = 16.132;
  //const double lumi = 20.0;
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
  // Overall Fake rate
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
  // Get from two fractions and two fake rates
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
      if( doprint ){
        std::cout << " fraction in Sample 1 " << std::endl;
        for(int ibin =0; ibin<= hfrac1t->GetNbinsX(); ibin++){
          std::cout << " bin " << std::setw(2) << ibin << " gets smeared from " << std::setw(12) << hfrac1_->GetBinContent(ibin) << " +/- "<< std::setw(14) << std::left << hfrac1_->GetBinError(ibin) << " to " << std::setw(15) << std::left << hfrac1t->GetBinContent(ibin) << std::endl;
        }
        std::cout << " fraction in Sample 2 " << std::endl;
        for(int ibin =0; ibin<= hfrac2t->GetNbinsX(); ibin++){
          std::cout << " bin " << std::setw(2) << ibin << " gets smeared from " << std::setw(12) << hfrac2_->GetBinContent(ibin) << " +/- "<< std::setw(14) << std::left << hfrac2_->GetBinError(ibin) << " to " << std::setw(15) << std::left << hfrac2t->GetBinContent(ibin) << std::endl;
        }
        std::cout << " fakerate in Sample 1 " << std::endl;
        for(int ibin =0; ibin<= hfr1t->GetNbinsX(); ibin++){
          std::cout << " bin " << std::setw(2) << ibin << " gets smeared from " << std::setw(12) << hfr1_->GetBinContent(ibin) << " +/- "<< std::setw(14) << std::left << hfr1_->GetBinError(ibin) << " to " << std::setw(15) << std::left << hfr1t->GetBinContent(ibin) << std::endl;
        }
        std::cout << " fakerate in Sample 2 " << std::endl;
        for(int ibin =0; ibin<= hfr2t->GetNbinsX(); ibin++){
          std::cout << " bin " << std::setw(2) << ibin << " gets smeared from " << std::setw(12) << hfr2_->GetBinContent(ibin) << " +/- "<< std::setw(14) << std::left << hfr2_->GetBinError(ibin) << " to " << std::setw(15) << std::left << hfr2t->GetBinContent(ibin) << std::endl;
        }
        std::cout << " bfraction in 0tag in QCD " << bfrac0_ << " new " << bfrac0 << std::endl;
        std::cout << " bfraction in 1tag in QCD " << bfrac1_ << " new " << bfrac1 << std::endl;
      }
    }

    TH1F* histo0to1 = FrHistoCal(hfrac1t, hfrac2t, hfr1t, hfr2t, bfrac0, "0to1", vcase1tag_[i]);
    TH1F* histo1to2 = FrHistoCal(hfrac1t, hfrac2t, hfr1t, hfr2t, bfrac1, "1to2", vcase2tag_[i]);
    vfrcal_temp0to1.push_back ( FrCal(histo0to1) );
    vfrcal_temp1to2.push_back ( FrCal(histo1to2) );
  } 
  vvfrcal_.push_back(vfrcal_temp0to1);
  vvfrcal_.push_back(vfrcal_temp1to2);
}

void EmJetEventCount::PrepareTruthFrCalVector(int ntimes)
{
  vector<FrCal> vfrcalQ_temp0to1;
  vector<FrCal> vfrcalQ_temp1to2;
  vector<FrCal> vfrcalG_temp0to1;
  vector<FrCal> vfrcalG_temp1to2;
  vector<FrCal> vfrcalGC_temp0to1;
  vector<FrCal> vfrcalGC_temp1to2;
  vector<FrCal> vfrcalQ_truth_B;
  vector<FrCal> vfrcalQ_truth_L;
  vector<FrCal> vfrcalG_truth_B;
  vector<FrCal> vfrcalG_truth_L;

  bool doprint = false;
  for(int i=0; i< ntimes; i++){
    if( i< 10 ) doprint = true;
    else doprint = false;

    TH1F* hQfrt = (TH1F*)hQfrb_->Clone((std::string(hQfr_->GetName())+"_"+std::to_string(i)).c_str() );// QCD MC b fake rates

    TH1F* hQfrbt = (TH1F*)hQfrb_->Clone((std::string(hQfrb_->GetName())+"_"+std::to_string(i)).c_str() );// QCD MC b fake rates
    TH1F* hQfrlt = (TH1F*)hQfrl_->Clone((std::string(hQfrl_->GetName())+"_"+std::to_string(i)).c_str() );// QCD MC l fake rates

    TH1F* hGfrbt = (TH1F*)hGfrb_->Clone((std::string(hGfrb_->GetName())+"_"+std::to_string(i)).c_str() );// GJet MC b fake rates
    TH1F* hGfrlt = (TH1F*)hGfrl_->Clone((std::string(hGfrl_->GetName())+"_"+std::to_string(i)).c_str() );

    TH1F* hGCfrbt = (TH1F*)hGCfrb_->Clone((std::string(hGCfrb_->GetName())+"_"+std::to_string(i)).c_str() );// Calc fake rates in data/MC
    TH1F* hGCfrlt = (TH1F*)hGCfrl_->Clone((std::string(hGCfrl_->GetName())+"_"+std::to_string(i)).c_str() );

    double bfrac0 = bfrac0_;
    double bfrac1 = bfrac1_;
    if( i!=0 ){
      SmearHisto(hQfrbt, doprint);
      SmearHisto(hQfrlt, doprint);
      SmearHisto(hGfrbt, doprint);
      SmearHisto(hGfrlt, doprint);
      SmearNumber(bfrac0, err_bfrac0_);
      SmearNumber(bfrac1, err_bfrac1_);
      if( doprint ) std::cout << " bfraction " << bfrac0_ << " new " << bfrac0 << std::endl;
      if( doprint ) std::cout << " bfraction " << bfrac1_ << " new " << bfrac1 << std::endl;
    }

    //TH1F* histoQ0to1 = FrHistoAdd( hQfrbt, hQfrlt, bfrac0, "QCD_0to1");
    TH1F* histoQ0to1 = hQfrt;
    TH1F* histoQ1to2 = FrHistoAdd( hQfrbt, hQfrlt, bfrac1, "QCD_1to2");
    vfrcalQ_temp0to1.push_back ( FrCal(histoQ0to1) );
    vfrcalQ_temp1to2.push_back ( FrCal(histoQ1to2) );

    TH1F* histoG0to1 = FrHistoAdd( hGfrbt, hGfrlt, bfrac0, "GJet_0to1");
    TH1F* histoG1to2 = FrHistoAdd( hGfrbt, hGfrlt, bfrac1, "GJet_1to2");
    vfrcalG_temp0to1.push_back ( FrCal(histoG0to1) );
    vfrcalG_temp1to2.push_back ( FrCal(histoG1to2) );

    TH1F* histoGC0to1 = FrHistoAdd( hGCfrbt, hGCfrlt, bfrac0, "GJetC_0to1");
    TH1F* histoGC1to2 = FrHistoAdd( hGCfrbt, hGCfrlt, bfrac1, "GJetC_1to2");
    vfrcalGC_temp0to1.push_back ( FrCal(histoGC0to1) );
    vfrcalGC_temp1to2.push_back ( FrCal(histoGC1to2) );

    vfrcalQ_truth_B.push_back( FrCal(hQfrbt) );
    vfrcalQ_truth_L.push_back( FrCal(hQfrlt) );
    vfrcalG_truth_B.push_back( FrCal(hGfrbt) );
    vfrcalG_truth_L.push_back( FrCal(hGfrlt) );    
  }
  vvfrcal_.push_back(vfrcalQ_temp0to1);
  vvfrcal_.push_back(vfrcalQ_temp1to2);
  vvfrcal_.push_back(vfrcalG_temp0to1);
  vvfrcal_.push_back(vfrcalG_temp1to2);
  vvfrcal_.push_back(vfrcalGC_temp0to1);
  vvfrcal_.push_back(vfrcalGC_temp1to2);
  vvfrcal_.push_back(vfrcalQ_truth_B);
  vvfrcal_.push_back(vfrcalQ_truth_L);
  vvfrcal_.push_back(vfrcalG_truth_B);
  vvfrcal_.push_back(vfrcalG_truth_L);
  std::cout << " Finished Preparing Truth Fakerate Histograms " << std::endl;
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
  for(int i=0; i<ntimes; i++){
    vcase1tag_.push_back(1);
    vcase2tag_.push_back(1);
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

void EmJetEventCount::PrintOutInfo()
{
  std::cout << "******************************************" << std::endl;
  std::cout << "******************************************" << std::endl;
  std::cout << " Event " << event << " run " << run << " lumi " << lumi << " pv_genreco_disXY " << pv_genreco_disXY << " pv_genreco_disZ " << pv_genreco_disZ << std::endl; 
  for(int ij=0; ij<4; ij++){
    std::cout << " jet pt " << (*jet_pt)[ij] << " eta " << (*jet_eta)[ij] << " phi "<< (*jet_phi)[ij] << " nTrack " << (*jet_nTrack)[ij] <<" nTrackPostCut " << (*jet_nTrackPostCut)[ij] << " jet_Alpha3DSig " << (*jet_Alpha3DSig)[ij] << " medianIP " << (*jet_medianIP)[ij] << " isEmerging? " << ((*jet_isEmerging)[ij] ? " yes " : " no " )<< std::endl;
  }
}
