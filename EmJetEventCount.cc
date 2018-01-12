#include "EmJetEventCount.h"

#ifndef OUTPUT
#define OUTPUT(x) std::cout<<#x << ": " << x << std::endl
#endif

using std::string;
using namespace FRFormula;

EmJetEventCount::EmJetEventCount(EmJetSampleCollection samplesColl, bool doFillHisto, bool doPredict, int repeatedTime)
:n1tag_(0.), n2tag_(0.)
{
  SetMaxEntries(-1);
  TChain* chain = new TChain("emJetSlimmedTree");
  std::cout << " EmJetEventCount::EmJetEventCount(" << samplesColl.name << ")" << std::endl;
  for (EmJetSample isample : samplesColl.samples) {
    std::cout << "Adding file to list of files to be processed: " << isample.file << std::endl;
    chain->Add(isample.file.c_str(), -1); // -1: the file is connected and the tree header read in memory to get the number of entries.
  }
  tree_ = chain;
  ntimes_ = repeatedTime;
  SetDataType(samplesColl.isData);
  SetFillOption(doFillHisto);
  SetPredictOption(doPredict);
  sfrfile_ = samplesColl.FrCalfile;
  Init(tree_);
  InitCrossSection(samplesColl);
}

void EmJetEventCount::LoopOverEvent(long eventnumber)
{
  if( pvtrack_fraction<0.1 ) return;
  int nJet_tag = 0;
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
  for(int itime=0; itime< ntimes_; itime++){ 
    PredictBackground(itime, nJet_tag, itime==0 && doFill_ );
  }
}

void EmJetEventCount::PredictBackground(int itime, int nJet_tag, bool isfillhisto)
{
  // results from GJet overall
  double afr0[4] = {-1.0, -1.0, -1.0, -1.0};
  for(int ij=0; ij<4; ij++){
    afr0[ij] = GetFR(histoFR_->histoFR1d["GJetOverall"][itime], (*jet_nTrackPostCut)[ij]);
  }
  vvn1tag_[0][itime] += PnTag(afr0, 1) * tweight_;
  vvn2tag_[0][itime] += PnTag(afr0, 2) * tweight_;

  // from 0tag
  if( nJet_tag==0 ){
    double afr1[4] = {-1.0, -1.0, -1.0, -1.0}; // results from GJet 0to1tag
    double afr3[4] = {-1.0, -1.0, -1.0, -1.0}; // results from QCD truth 0to1tag
    double afr5[4] = {-1.0, -1.0, -1.0, -1.0}; // results from GJet truth 0to1tag 
    //double afr9[4] = {-1.0, -1.0, -1.0, -1.0}; // results from QCD using truth flavor

    for( int ij=0; ij<4; ij++){
      afr1[ij] = GetFR(histoFR_->histoFR1d["GJet0To1"][itime],   (*jet_nTrackPostCut)[ij]);
      afr3[ij] = GetFR(histoFR_->histoFR1d["QCDOverall"][itime],   (*jet_nTrackPostCut)[ij]);
      afr5[ij] = GetFR(histoFR_->histoFR1d["GJetT0To1"][itime],  (*jet_nTrackPostCut)[ij]);

      /*
      if( (*jet_flavour)[ij]==5 || (*jet_flavour)[ij]==19 || (*jet_flavour)[ij]==10 ) {
        afr9[ij] = vvfrcal_[9][itime].GetFakerate((*jet_nTrackPostCut)[ij]);
      }
      else{
        afr9[ij] = vvfrcal_[10][itime].GetFakerate((*jet_nTrackPostCut)[ij]);
      }
      */

    }
    vvn1tag_[1][itime] += PnTag(afr1, 1) * tweight_;
    vvn1tag_[2][itime] += PnTag(afr3, 1) * tweight_;
    vvn1tag_[3][itime] += PnTag(afr5, 1) * tweight_;

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
    int ijet=0;
    for(int ij=0; ij<4; ij++){
      if( (*jet_isEmerging)[ij] ) continue;// skip the emerging jet
      //if( (*jet_Alpha3DSig)[ij]<0.45 && (*jet_medianIP)[ij]>0.05 ) continue;
      afr2[ijet] = GetFR(histoFR_->histoFR1d["GJet1To2"][itime],   (*jet_nTrackPostCut)[ij]);// 1to2tag
      afr4[ijet] = GetFR(histoFR_->histoFR1d["QCDT1To2"][itime],   (*jet_nTrackPostCut)[ij]);// 1to2tag QCD truth
      afr6[ijet] = GetFR(histoFR_->histoFR1d["GJetT1To2"][itime],  (*jet_nTrackPostCut)[ij]);// 1to2tag GJet truth
      ijet++;
    }
    vvn2tag_[1][itime] += P1tagTo2tag(afr2) * tweight_;
    vvn2tag_[2][itime] += P1tagTo2tag(afr4) * tweight_;
    vvn2tag_[3][itime] += P1tagTo2tag(afr6) * tweight_;

    if( isfillhisto ) {
      FillClosureTestHistos1To2Tag(afr2, "__GJetCalcPredicted1To2Tag");
      FillClosureTestHistos1To2Tag(afr4, "__QCDTruthPredicted1To2Tag");
      FillClosureTestHistos1To2Tag(afr6, "__GJetTruthPredicted1To2Tag");
    }
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

void EmJetEventCount::FinishEventCount()
{
  FillEventCountHistos();
  PrintOutResults();
}

void EmJetEventCount::FillEventCountHistos()
{
  // fill the number of events with 1 tags
  for(int itime=0; itime < ntimes_; itime++){
    histo_->hist1d["n1tag__0"]->Fill(vvn1tag_[0][itime]);
    histo_->hist1d["n1tag__1"]->Fill(vvn1tag_[1][itime]);
    histo_->hist1d["n1tag__2"]->Fill(vvn1tag_[2][itime]);
    histo_->hist1d["n1tag__3"]->Fill(vvn1tag_[3][itime]);
    histo_->hist1d["n1tag__4"]->Fill(vvn1tag_[4][itime]);
    if( vcase1tag_[itime]==1 )      histo_->hist1d["N1tag_case1__1"]->Fill(vvn1tag_[1][itime]);
    else if( vcase1tag_[itime]==2 ) histo_->hist1d["N1tag_case2__1"]->Fill(vvn1tag_[1][itime]);
  }

  for(int itime=0; itime < ntimes_; itime++){
    histo_->hist1d["n1tag_EM1__0"]->Fill(vvn1tag_[0][itime]);
    histo_->hist1d["n1tag_EM1__1"]->Fill(vvn1tag_[1][itime]);
    histo_->hist1d["n1tag_EM1__2"]->Fill(vvn1tag_[2][itime]);
    histo_->hist1d["n1tag_EM1__3"]->Fill(vvn1tag_[3][itime]);
    histo_->hist1d["n1tag_EM1__4"]->Fill(vvn1tag_[4][itime]);
  }

  // fill the number of events with 2 tags
  for(int itime=0; itime < ntimes_; itime++){
    histo_->hist1d["n2tag__0"]->Fill(vvn2tag_[0][itime]);
    histo_->hist1d["n2tag__1"]->Fill(vvn2tag_[1][itime]);
    histo_->hist1d["n2tag__2"]->Fill(vvn2tag_[2][itime]);
    histo_->hist1d["n2tag__3"]->Fill(vvn2tag_[3][itime]);
    histo_->hist1d["n2tag__4"]->Fill(vvn2tag_[4][itime]);
    if( vcase2tag_[itime]==1 )      histo_->hist1d["N2tag_case1__1"]->Fill(vvn1tag_[1][itime]);
    else if( vcase2tag_[itime]==2 ) histo_->hist1d["N2tag_case2__1"]->Fill(vvn1tag_[1][itime]);
  }
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

void EmJetEventCount::InitFrHistos(string fFRname)
{
  if(doPredict_){
    histoFR_ = unique_ptr<EmJetFrHistos>(new EmJetFrHistos(fFRname, ntimes_, isData_));
  }
}

void EmJetEventCount::InitHistograms()
{
  TH1::SetDefaultSumw2();
  histo_ = unique_ptr<Histos>(new Histos());
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

void EmJetEventCount::InitEventCount()
{ 
  InitFrHistos(sfrfile_);
  PrepareFrCalResults();
}

void EmJetEventCount::PrepareFrCalResults()
{ 
  for(unsigned ifrcal=0; ifrcal < 10; ifrcal++){
    vector<double> vinit;
    for(int i=0; i<ntimes_; i++){
      vinit.push_back(0.); 
    }
    vvn1tag_.push_back(vinit);
    vvn2tag_.push_back(vinit);
  }
  for(int i=0; i<ntimes_; i++){
    vcase1tag_.push_back(1);
    vcase2tag_.push_back(1);
  }
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
