#include "Fakerate.h"

FrCal::FrCal(): histo_(0), histoname_(""), nbins_(0)
{
}

FrCal::FrCal(string filename, string histoname)
{
  TFile *file = new TFile(filename.c_str());  
  histo_ = (TH1F*)file->Get(histoname.c_str())->Clone((histoname + "_").c_str());
  histoname_ = histoname;
  nbins_ = histo_->GetNbinsX();
}

void FrCal::SmearFrHisto(bool doprint)
{
  if( doprint ){
    std::cout << "Start smearing the fakerate histogram " << histoname_ << " with Gaussian" << std::endl; 
  }
  gRandom = new TRandom3(0);
  for(int ibin = 1; ibin <= histo_->GetNbinsX(); ibin++){
    double newval = gRandom->Gaus( histo_->GetBinContent(ibin), histo_->GetBinError(ibin));
    if( doprint ){
      std::cout << " bin " << std::setw(2) << ibin << " gets smeared from " << std::setw(12) << histo_->GetBinContent(ibin) << " +/- "<< std::setw(14) << std::left << histo_->GetBinError(ibin) << " to " << std::setw(15) << std::left << newval << std::endl;
    }
    histo_->SetBinContent(ibin, newval);
  }
}

void FrCal::SmearHistoBy1Sigma(bool doprint, int sign)
{
  if( doprint ){
    std::cout << "Start smearing the fakerate histogram " << histoname_ << " with Gaussian" << std::endl;
  }
  int Sign = sign/abs(sign);
  for(int ibin = 1; ibin <= histo_->GetNbinsX(); ibin++){
    double newval = histo_->GetBinContent(ibin) + Sign *histo_->GetBinError(ibin);
    if( doprint ){
      std::cout << " bin " << std::setw(2) << ibin << " gets smeared from " << std::setw(12) << histo_->GetBinContent(ibin) << " +/- "<< std::setw(14) << std::left << histo_->GetBinError(ibin) << " to " << std::setw(15) << std::left << newval << std::endl;
    }
    histo_->SetBinContent(ibin, newval);
  } 
}

FrCal FrCal::Clone(string histoname)
{
  FrCal ofr;
  ofr.histo_ = (TH1F*)histo_->Clone(histoname.c_str());
  ofr.nbins_ = nbins_;
  ofr.histoname_ = histoname;
  return ofr;
}

double FrCal::GetFakerate(int nTrack) const
{
  int ibin = histo_->FindBin(nTrack);
  double fr = -1.0;
  if( ibin<=nbins_ ){
    fr = histo_->GetBinContent(ibin);
  }
  else{
    fr = histo_->GetBinContent(nbins_); // overflow, use the fakerate from the maximum bin
  }
  return fr;
}

string FrCal::GetHistoName()
{
  return histoname_;
}

//
// calculate the probability/weight of one event with fr[4] and nTags
//
double PnTag(double fr[], int nTag){
  if( fr[0]<0.0 || fr[1]<0.0 || fr[2]<0.0 || fr[3]<0.0 ){
    std::cerr << "Error! Fakerate arrays not calculated correctly" << std::endl;
    return 0.0;
  }
  double p_nTag = 0;
  if( nTag == 0){
    p_nTag = (1-fr[0]) * (1-fr[1]) * (1-fr[2]) * (1-fr[3]);
  }
  else if( nTag == 1 ){
    p_nTag  =   fr[0]   * (1-fr[1]) * (1-fr[2]) * (1-fr[3]); 
    p_nTag += (1-fr[0]) *   fr[1]   * (1-fr[2]) * (1-fr[3]);
    p_nTag += (1-fr[0]) * (1-fr[1]) *   fr[2]   * (1-fr[3]);
    p_nTag += (1-fr[0]) * (1-fr[1]) * (1-fr[2]) *   fr[3]  ;
  }
  else if( nTag == 2 ){
    p_nTag  = (1-fr[0]) * (1-fr[1]) *   fr[2]   *   fr[3]  ;
    p_nTag += (1-fr[0]) *   fr[1]   * (1-fr[2]) *   fr[3]  ;
    p_nTag +=   fr[0]   * (1-fr[1]) * (1-fr[2]) *   fr[3]  ;
    p_nTag += (1-fr[0]) *   fr[1]   *   fr[2]   * (1-fr[3]);
    p_nTag +=   fr[0]   * (1-fr[1]) *   fr[2]   * (1-fr[3]);
    p_nTag +=   fr[0]   *   fr[1]   * (1-fr[2]) * (1-fr[3]);
  }
  else {
    std::cout << " Error: events weight with nTag = " << nTag << " has not been implemented yet!!!" << std::endl;
    return -1.0;
  }
  return p_nTag;
}

// 
// Calculate the probability/weight when predicting from 1-tag to 2-tag
//
double P1tagTo2tag(double fr[]){
  double p1tagto2tag = 0;
  p1tagto2tag  =   fr[0]   * (1-fr[1]) * (1-fr[2]);
  p1tagto2tag += (1-fr[0]) *   fr[1]   * (1-fr[2]);
  p1tagto2tag += (1-fr[0]) * (1-fr[1]) *   fr[2]  ;
  return p1tagto2tag/2.0;
}

//
// compute the probability/weight of ijet is tagged as emerging
// in the nTag case(nTag = 0, 1, 2)
//
double PEmergingnTag(double fr[], int nTag, int ijet){
  if( ijet >=4 ) {
    std::cout << " Error: jet index out when calculating Emerging probability" << std::endl; 
    return 0;
  }
  if( fr[0]<0.0 || fr[1]<0.0 || fr[2]<0.0 || fr[3]<0.0 ){
    std::cerr << "Error! Fakerate arrays not calculated correctly" << std::endl;
    return 0.0;
  }
  double prob = 0.0;
  if ( nTag==1 ) {
    prob = fr[ijet];
    for(int ij = 0; ij<4; ij++){
      if(ij==ijet) continue;
      prob *= (1-fr[ij]);
    }
  }
  if( nTag==2 ){
    for(int ij1 = 0; ij1<4; ij1++){
      if( ij1== ijet ) continue;
      double prob_temp = fr[ijet] * fr[ij1];
      for(int ij2 = 0; ij2<4; ij2++){
        if( ij2==ij1 || ij2==ijet ) continue;
        prob_temp *= (1-fr[ij2]); 
      }
      prob += prob_temp;
    }
  }

  if( prob>1.0 || (prob<0 && nTag!=0 ) ){
    std::cout << " Error/DEBUG: probability calculation problem !!! " << std::endl;
  }
  return prob;
}

//
// compuute the probability/weight of ijet is tagged as emerging
// when predicting from 1tag to 2tag
//
double PEmerging1tagTo2tag(double fr[], int ijet){
  if( ijet>=3 ) {
    std::cerr<< " Error: input jet index problem !!!" << std::endl;
    return -1;
  }
  double prob = fr[ijet];
  for(int ij=0; ij<3; ij++){
    if( ij==ijet ) continue;
    prob *= (1-fr[ij]);
  } 
  return prob/2.0;
}

double frCal(int jet_nTrack, int option){
  // option 0 : overall fakerate
  // option 1 : tagged as light quark or gluon
  // option 2 : tagged as b quark
  double fakerate = -10.0;
  if (option == 0) { // Overall fakerate
    if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.0294997282326;
    else if( jet_nTrack>=4 && jet_nTrack<8 ) fakerate = 0.00641810335219;
    else if( jet_nTrack>=8 && jet_nTrack<12 ) fakerate = 0.00324129522778;
    else if( jet_nTrack>=12 && jet_nTrack<16 ) fakerate = 0.00187143380754;
    else if( jet_nTrack>=16 && jet_nTrack<20 ) fakerate = 0.00116369535681;
    else if( jet_nTrack>=20 && jet_nTrack<25 ) fakerate = 0.000727684178855;
    else if( jet_nTrack>=25 && jet_nTrack<30 ) fakerate = 0.000449172075605;
    else if( jet_nTrack>=30 && jet_nTrack<40 ) fakerate = 0.000223103736062;
    else if( jet_nTrack>=40 && jet_nTrack<50 ) fakerate = 7.17896109563e-05;
    else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 7.02791803633e-05;
    else if( jet_nTrack>=60 ) fakerate = 0.000173533204361;
  }    
  else if( option== 1){ // light quark (u, d, s, c) or gluon jet
    if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.0263485852629;
    else if( jet_nTrack>=4 && jet_nTrack<8 ) fakerate = 0.00421362370253;
    else if( jet_nTrack>=8 && jet_nTrack<12 ) fakerate = 0.00189180707093;
    else if( jet_nTrack>=12 && jet_nTrack<16 ) fakerate = 0.00129529181868;
    else if( jet_nTrack>=16 && jet_nTrack<20 ) fakerate = 0.000852863537148;
    else if( jet_nTrack>=20 && jet_nTrack<25 ) fakerate = 0.000600345374551;
    else if( jet_nTrack>=25 && jet_nTrack<30 ) fakerate = 0.000403687212383;
    else if( jet_nTrack>=30 && jet_nTrack<40 ) fakerate = 0.00021440727869;
    else if( jet_nTrack>=40 && jet_nTrack<50 ) fakerate = 7.00655509718e-05;
    else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 7.58130263421e-05;
    else if( jet_nTrack>=60 ) fakerate = 0.000188055637409;
  }
  else if( option == 2 ){ // B jet fakerate
    if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.112466610968;
    else if( jet_nTrack>=4 && jet_nTrack<8 ) fakerate = 0.0508545823395;
    else if( jet_nTrack>=8 && jet_nTrack<12 ) fakerate = 0.0281637758017;
    else if( jet_nTrack>=12 && jet_nTrack<16 ) fakerate = 0.0133286537603;
    else if( jet_nTrack>=16 && jet_nTrack<20 ) fakerate = 0.00781914126128;
    else if( jet_nTrack>=20 && jet_nTrack<25 ) fakerate = 0.00345264398493;
    else if( jet_nTrack>=25 && jet_nTrack<30 ) fakerate = 0.0013583204709;
    else if( jet_nTrack>=30 && jet_nTrack<40 ) fakerate = 0.000376539537683;
    else if( jet_nTrack>=40 && jet_nTrack<50 ) fakerate = 9.7794640169e-05;
    else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 0.0;
    else if( jet_nTrack>=60 ) fakerate = 0.0;
  }
  else if( option == 3){
    if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.0304408613592;
    else if( jet_nTrack>=4 && jet_nTrack<8 ) fakerate = 0.00642997398973;
    else if( jet_nTrack>=8 && jet_nTrack<12 ) fakerate = 0.00314023531973;
    else if( jet_nTrack>=12 && jet_nTrack<16 ) fakerate = 0.00186710990965;
    else if( jet_nTrack>=16 && jet_nTrack<20 ) fakerate = 0.00118389690761;
    else if( jet_nTrack>=20 && jet_nTrack<25 ) fakerate = 0.000735884881578;
    else if( jet_nTrack>=25 && jet_nTrack<30 ) fakerate = 0.000449050799944;
    else if( jet_nTrack>=30 && jet_nTrack<40 ) fakerate = 0.000222111702897;
    else if( jet_nTrack>=40 && jet_nTrack<50 ) fakerate = 7.13832196197e-05;
    else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 7.22104377928e-05;
    else if( jet_nTrack>=60 ) fakerate = 0.000179119349923;
  }
  else if( option == 4) {
    if( jet_nTrack>=0 && jet_nTrack<4 ) fakerate = 0.0346605852246;
    else if( jet_nTrack>=4 && jet_nTrack<8 ) fakerate = 0.00871534831822;
    else if( jet_nTrack>=8 && jet_nTrack<12 ) fakerate = 0.00442754337564;
    else if( jet_nTrack>=12 && jet_nTrack<16 ) fakerate = 0.00245673628524;
    else if( jet_nTrack>=16 && jet_nTrack<20 ) fakerate = 0.0015252395533;
    else if( jet_nTrack>=20 && jet_nTrack<25 ) fakerate = 0.000875645549968;
    else if( jet_nTrack>=25 && jet_nTrack<30 ) fakerate = 0.000495827174746;
    else if( jet_nTrack>=30 && jet_nTrack<40 ) fakerate = 0.000230056073633;
    else if( jet_nTrack>=40 && jet_nTrack<50 ) fakerate = 7.27419246687e-05;
    else if( jet_nTrack>=50 && jet_nTrack<60 ) fakerate = 6.84956539772e-05;
    else if( jet_nTrack>=60 ) fakerate = 0.000169904757058;
  }
  return fakerate;
}
