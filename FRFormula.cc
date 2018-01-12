#include "FRFormula.h"

//
// calculate the probability/weight of one event with fr[4] and nTags
//
double FRFormula::PnTag(double fr[], int nTag){
  /*
  if( fr[0]<0.0 || fr[1]<0.0 || fr[2]<0.0 || fr[3]<0.0 ){
    std::cout << "Error! Fakerate arrays not calculated correctly" << std::endl;
    return 0.0;
  }
  */
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
double FRFormula::P1tagTo2tag(double fr[]){
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
double FRFormula::PEmergingnTag(double fr[], int nTag, int ijet){
  if( ijet >=4 ) {
    std::cout << " Error: jet index out when calculating Emerging probability" << std::endl; 
    return 0;
  }
  /*
  if( fr[0]<0.0 || fr[1]<0.0 || fr[2]<0.0 || fr[3]<0.0 ){
    std::cout << "Error! Fakerate arrays not calculated correctly" << std::endl;
    return 0.0;
  } 
  */
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
  /*
  if( prob>1.0 || (prob<0 && nTag!=0 ) ){
    std::cout << " Error/DEBUG: probability calculation problem !!! Probability is " << prob << " fakerates " << fr[0]  << " "<< fr[1] << " " << fr[2] << " " << fr[3] << std::endl;
  }
  */
  return prob;
}

//
// compute the probability/weight of ijet is tagged as emerging
// when predicting from 1tag to 2tag
//
double FRFormula::PEmerging1tagTo2tag(double fr[], int ijet){
  if( ijet>=3 ) {
    std::cout<< " Error: input jet index problem !!!" << std::endl;
    return -1;
  }
  double prob = fr[ijet];
  for(int ij=0; ij<3; ij++){
    if( ij==ijet ) continue;
    prob *= (1-fr[ij]);
  } 
  return prob/2.0;
}

TH1F* FRFormula::FrHistoCal(TH1F* hfrac1, TH1F* hfrac2, TH1F* hfr1, TH1F* hfr2, double bfrac, std::string tag, int idx){
  TH1F* hfr = (TH1F*)hfr1->Clone(("hfr_calc"+tag+"_"+std::to_string(idx)).c_str());

  for(int i=1; i<= hfr1->GetNbinsX(); i++){  
    int ibinfrac = hfrac1->FindBin(hfr1->GetBinCenter(i));

    double fb1 = hfrac1->GetBinContent(ibinfrac);
    double fl1 = 1.0 - hfrac1->GetBinContent(ibinfrac);
    double fb2 = hfrac2->GetBinContent(ibinfrac);
    double fl2 = 1.0 - hfrac2->GetBinContent(ibinfrac);
    double FR1 = hfr1->GetBinContent(i);
    double FR2 = hfr2->GetBinContent(i);    

    double norm = 1.0/(fb1- fb2);
    double FR_b = norm*( fl2*FR1 - fl1*FR2);
    double FR_l = norm*(-fb2*FR1 + fb1*FR2);

    FR_b = ( FR_b>=0.0 ? FR_b : 0.0); FR_b = ( FR_b<=1.0 ? FR_b : 1.0);
    FR_l = ( FR_l>=0.0 ? FR_l : 0.0); FR_l = ( FR_l<=1.0 ? FR_l : 1.0);

    hfr->SetBinContent(i, FR_b * bfrac + FR_l * (1-bfrac));
  }
  return hfr;
} 

TH1F* FRFormula::FrHistoAdd(TH1F* hfrb, TH1F* hfrl, double bfrac, std::string tag, int idx){
  TH1F* hfr = (TH1F*)hfrb->Clone(("hfr_add"+tag+"_"+std::to_string(idx)).c_str());

  for(int i=1; i<= hfrb->GetNbinsX(); i++){

    double FR_b = hfrb->GetBinContent(i);
    double FR_l = hfrl->GetBinContent(i);

    hfr->SetBinContent(i, FR_b * bfrac + FR_l * (1-bfrac));
  }
  return hfr;
}

double FRFormula::GetFR(TH1F* hfr, int nTrack)
{
  int ibin = hfr->FindBin(nTrack);
  double FR = 0.0;
  if( ibin>0 && ibin<=hfr->GetNbinsX() ){
    FR = hfr->GetBinContent(ibin); 
  }
  else{
    // Overflow, use the results from maximum bin
    FR = hfr->GetBinContent(hfr->GetNbinsX());
    //std::cout << " Warning: Over/Under flow for " << hfr->GetName() << " with nTrack " << nTrack << std::endl;
    //std::cout << "    use maximum bin result " << std::endl;
  }
  return FR;
}

double FRFormula::GetRawFakerate(int nTrack, bool isBJet)
{
  double fakerate = 0.0;
  if( isBJet ){
    if( nTrack>=0.0 && nTrack<4.0 ) {
      fakerate = 0.0764036476612;
    }
    if( nTrack>=4.0 && nTrack<8.0 ) {
      fakerate = 0.0240180138499;
    }
    if( nTrack>=8.0 && nTrack<12.0 ) {
      fakerate = 0.00735678197816;
    }
    if( nTrack>=12.0 && nTrack<16.0 ) {
      fakerate = 0.00186479813419;
    }
    if( nTrack>=16.0 && nTrack<20.0 ) {
      fakerate = 0.000393063208321;
    }
    if( nTrack>=20.0 && nTrack<24.0 ) {
      fakerate = 0.000105553546746;
    }
    if( nTrack>=24.0 && nTrack<40.0 ) {
      fakerate = 1.39190433401e-06;
    }
    if( nTrack>=40.0 && nTrack<80.0 ) {
      fakerate = 0.0;
    }
  }
  else{
    if( nTrack>=0.0 && nTrack<4.0 ) {
      fakerate = 0.0272194761783;
    }
    if( nTrack>=4.0 && nTrack<8.0 ) {
      fakerate = 0.0024541572202;
    }
    if( nTrack>=8.0 && nTrack<12.0 ) {
      fakerate = 0.000505852687638;
    }
    if( nTrack>=12.0 && nTrack<16.0 ) {
      fakerate = 0.00015330662427;
    }
    if( nTrack>=16.0 && nTrack<20.0 ) {
      fakerate = 5.15911160619e-05;
    }
    if( nTrack>=20.0 && nTrack<24.0 ) {
      fakerate = 2.135487739e-05;
    }
    if( nTrack>=24.0 && nTrack<40.0 ) {
      fakerate = 6.16831357547e-06;
    }
    if( nTrack>=40.0 && nTrack<80.0 ) {
      fakerate = 0.0;
    }
  }
  return fakerate;
}
