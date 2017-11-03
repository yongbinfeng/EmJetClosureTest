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

FrCal::FrCal(TH1F* hfr)
{
  histo_ = hfr;
  histoname_ = hfr->GetName();
  nbins_ = histo_->GetNbinsX();
}

void FrCal::SmearFrHisto(bool doprint)
{
  SmearHisto(histo_, doprint);
}

void FrCal::SmearHistoBy1Sigma(bool doprint, int sign)
{
  if( doprint ){
    std::cout << "Start smearing the fakerate histogram " << histoname_ << " with Gaussian" << std::endl;
  }
  int Sign = sign/abs(sign);
  for(int ibin = 1; ibin <= histo_->GetNbinsX(); ibin++){
    double newval = histo_->GetBinContent(ibin) + Sign *histo_->GetBinError(ibin);
    newval = newval>=0.0 ? newval : 0.0;
    newval = newval<=1.0 ? newval : 1.0;
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

TH1F* FrHistoCal(TH1F* hfrac1, TH1F* hfrac2, TH1F* hfr1, TH1F* hfr2, double bfrac, double err_bfrac, string tag){
  TH1F* hfr = (TH1F*)hfr1->Clone(("hfr_calc"+tag).c_str());

  gRandom = new TRandom3(0);
  double bfrac_temp = gRandom->Gaus(bfrac, err_bfrac);

  for(int i=1; i<= hfr1->GetNbinsX(); i++){  
    int ibinfrac = hfrac1->FindBin(hfr1->GetBinCenter(i));

    double fb1 = 1.0 - hfrac1->GetBinContent(ibinfrac);
    double fl1 = hfrac1->GetBinContent(ibinfrac);
    double fb2 = 1.0 - hfrac2->GetBinContent(ibinfrac);
    double fl2 = hfrac2->GetBinContent(ibinfrac);
    double FR1 = hfr1->GetBinContent(i);
    double FR2 = hfr2->GetBinContent(i);    

    double norm = 1.0/(fb1- fb2);

    double FR_b = ( norm*( fl2*FR1 - fl1*FR2)>=0.0 ?  norm*( fl2*FR1 - fl1*FR2) : 0.0 );
    double FR_l = ( norm*(-fb2*FR1 + fb1*FR2)>=0.0 ?  norm*(-fb2*FR1 + fb1*FR2) : 0.0 );

    hfr->SetBinContent(i, FR_b * bfrac_temp + FR_l * (1-bfrac_temp));
  }
  return hfr;
} 

void SmearHisto(TH1F* histo, bool doprint, bool mustbepositive)
{
  if( doprint ){
    std::cout << "Start smearing the histogram " << histo->GetName() << " with Gaussian" << std::endl;
  }
  gRandom = new TRandom3(0);
  for(int ibin = 1; ibin <= histo->GetNbinsX(); ibin++){
    double newval = gRandom->Gaus( histo->GetBinContent(ibin), histo->GetBinError(ibin));
    if( mustbepositive ){
      newval = newval>=0.0 ? newval : 0.0;
      newval = newval<=1.0 ? newval : 1.0;
    }
    if( doprint ){
      std::cout << " bin " << std::setw(2) << ibin << " gets smeared from " << std::setw(12) << histo->GetBinContent(ibin) << " +/- "<< std::setw(14) << std::left << histo->GetBinError(ibin) << " to " << std::setw(15) << std::left << newval << std::endl;
    }
    histo->SetBinContent(ibin, newval);
  }
}

