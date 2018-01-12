#include "EmJetFrHistos.h"

using namespace FRFormula;

EmJetFrHistos::EmJetFrHistos(string fFRname, int ntimes, bool isData)
{
  SetDataType(isData); 
  SetRepeatedTime(ntimes);
  OpenFrRootFile(fFRname);
  PrepareFrVector();
}

EmJetFrHistos::~EmJetFrHistos()
{
  for (auto & kv : histoFR1d) {
    for(auto hist: kv.second){
      delete hist;
    }
    kv.second.clear();
  }
}

void EmJetFrHistos::SetDataType(bool isData)
{
  if( isData )
    DataType_ = "Data";
  else 
    DataType_ = "MC";
  std::cout << " Set DataType to " << DataType_ << " in EmJetFrHistos " << std::endl;
}

void EmJetFrHistos::SetRepeatedTime(int ntimes)
{
  ntimes_ = ntimes;
  std::cout << " sudo experiments are set to be repeated " << ntimes_ << " for uncertainty calculation " << std::endl;  
  InitRandGen();
}

void EmJetFrHistos::InitRandGen()
{
  rand3_ = new TRandom3();
  std::cout << " TRandom3 initialized.." << std::endl;
}

void EmJetFrHistos::OpenFrRootFile(string fFRname)
{
  ffr_ = new TFile(fFRname.c_str());
  std::cout << " Read Fakerate and Fraction histograms from " << fFRname << std::endl;
}

void EmJetFrHistos::PrepareFrVector()
{ 
  PrepareBFractions();
  PrepareOverallFrVector();
  PrepareTwoTypeVector();
  PrepareTruthFrCalVector();    
  std::cout << "  All fakerate and fraction histograms set up..." << std::endl;
}

void EmJetFrHistos::PrepareBFractions()
{
  TH1F* hfracB = (TH1F*)ffr_->Get("bfraction_in_tags");
  bfrac0_ = hfracB->GetBinContent(1); err_bfrac0_ = hfracB->GetBinError(1);
  bfrac1_ = hfracB->GetBinContent(2); err_bfrac1_ = hfracB->GetBinError(2);
  std::cout << " B jet fractions read from bfraction_in_tags, finished..." << std::endl;
}

void EmJetFrHistos::PrepareOverallFrVector()
{
  // overall histograms
  hfrGJet_ = GetHisto("Dfakerates/fakerate_GJet"+DataType_);
  hfrQCD_  = GetHisto("Dfakerates/fakerate_QCDMC"); 

  for(int i=0; i<ntimes_; i++){
    TH1F* htemp_frGJet = SmearHisto(hfrGJet_, i);
    TH1F* htemp_frQCD  = SmearHisto(hfrQCD_,  i);;
    histoFR1d["GJetOverall"].push_back(htemp_frGJet);
    histoFR1d["QCDOverall"].push_back(htemp_frQCD);
  }
  std::cout << " Overall Histograms of GJet and QCD have been set up..." << std::endl;
}

void EmJetFrHistos::PrepareTwoTypeVector()
{
  // GJet fakerate from type I and type II
  hfrGJetI_    = GetHisto("fakerate_GJet"+DataType_+"_TypeVI");
  hfrGJetII_   = GetHisto("fakerate_GJet"+DataType_+"_TypeV");
  hfracGJetI_  = GetHisto("fraction_GJet"+DataType_+"_TypeVI");
  hfracGJetII_ = GetHisto("fraction_GJet"+DataType_+"_TypeV");

  for(int i=0; i<ntimes_; i++){
    TH1F* htemp_fracGJetI  = SmearHisto(hfracGJetI_,  i);
    TH1F* htemp_fracGJetII = SmearHisto(hfracGJetII_, i);
    
    TH1F* htemp_frGJetI  = SmearHisto(hfrGJetI_,  i);
    TH1F* htemp_frGJetII = SmearHisto(hfrGJetII_, i);

    double bfrac0 = SmearNumber(bfrac0_,  err_bfrac0_, i);
    double bfrac1 = SmearNumber(bfrac1_,  err_bfrac1_, i);

    TH1F* histo0to1 = FrHistoCal(htemp_fracGJetI, htemp_fracGJetII, htemp_frGJetI, htemp_frGJetII, bfrac0, "0to1",  i);
    TH1F* histo1to2 = FrHistoCal(htemp_fracGJetI, htemp_fracGJetII, htemp_frGJetI, htemp_frGJetII, bfrac1, "1to2",  i);
    histoFR1d["GJet0To1"].push_back(histo0to1);
    histoFR1d["GJet1To2"].push_back(histo1to2);
  }
  std::cout << "  Two types of GJet fakerate and fraction histograms set up..." << std::endl;
}

void EmJetFrHistos::PrepareTruthFrCalVector()
{
  // GJet and QCD truth fakerates
  hfrQCDB_ = GetHisto("Dfakerates/fakerate_QCDMC_BJet");
  hfrQCDL_ = GetHisto("Dfakerates/fakerate_QCDMC_LightJet");

  hfrGJetB_ = GetHisto("Dfakerates/fakerate_GJetMC_BJet");
  hfrGJetL_ = GetHisto("Dfakerates/fakerate_GJetMC_LightJet");

  for(int i=0; i<ntimes_; i++){
    TH1F* htemp_frQCDB = SmearHisto(hfrQCDB_,   i);
    TH1F* htemp_frQCDL = SmearHisto(hfrQCDL_,   i);

    TH1F* htemp_frGJetB = SmearHisto(hfrGJetB_,  i);   
    TH1F* htemp_frGJetL = SmearHisto(hfrGJetL_,  i); 
  
    double bfrac0 = SmearNumber(bfrac0_,  err_bfrac0_,  i);
    double bfrac1 = SmearNumber(bfrac1_,  err_bfrac1_,  i);

    TH1F* histoGT0to1 = FrHistoAdd( htemp_frGJetB, htemp_frGJetL, bfrac0, "GJet_0to1", i);
    TH1F* histoGT1to2 = FrHistoAdd( htemp_frGJetB, htemp_frGJetL, bfrac1, "GJet_1to2", i);
 
    TH1F* histoQT0to1 = FrHistoAdd( htemp_frQCDB, htemp_frQCDL,   bfrac0, "QCD_0to1", i);
    TH1F* histoQT1to2 = FrHistoAdd( htemp_frQCDB, htemp_frQCDL,   bfrac1, "QCD_1to2", i);
    
    histoFR1d["GJetT0To1"].push_back(histoGT0to1);
    histoFR1d["GJetT1To2"].push_back(histoGT1to2);    
  
    histoFR1d["QCDT0To1"].push_back(histoQT0to1);
    histoFR1d["QCDT1To2"].push_back(histoQT1to2);
  } 
  std::cout << "  Truth B jet and light jet fakerate histograms in GJet and QCD samples set up..." << std::endl;
}

TH1F* EmJetFrHistos::GetHisto(string hname)
{
  return (TH1F*)ffr_->Get(hname.c_str()); 
}

TH1F* EmJetFrHistos::SmearHisto(TH1F* h, int idx)
{
  //Clone the original histogram first
  TH1F* hsmeared = CloneHisto(h, idx);

  for(int ibin = 1; ibin <= hsmeared->GetNbinsX(); ibin++){
    double newval = SmearNumber(hsmeared->GetBinContent(ibin), hsmeared->GetBinError(ibin), idx);
    hsmeared->SetBinContent(ibin, newval);
  }

  // Print out the comparison for the first 10 histograms
  if( idx<10 ) PrintHistoComp(h, hsmeared);
  return hsmeared;
}

double EmJetFrHistos::SmearNumber(double val, double err, int idx)
{
  // the 0th one is set to the same as the orignal one
  if(idx==0) return val; 

  double nval = rand3_->Gaus(val, err); 
  // nval should be in [0, 1]
  nval = nval>=0.0 ? nval : 0.0;
  nval = nval<=1.0 ? nval : 1.0;
  return nval;
}

void EmJetFrHistos::PrintHistoComp(TH1F* h, TH1F* hsmeared)
{
  std::cout << " Histogram " << hsmeared->GetName() << " got smeared with Gaussian" << std::endl;
  for(int ibin=1; ibin<=h->GetNbinsX(); ibin++){
    std::cout << " Bin " << ibin << " smeared from " << h->GetBinContent(ibin) << " +/- " << h->GetBinError(ibin) << " to " << hsmeared->GetBinContent(ibin) << std::endl; 
  }  
}

TH1F* EmJetFrHistos::CloneHisto(TH1F* h, int idx)
{
   return (TH1F*)h->Clone((std::string(h->GetName())+"_"+std::to_string(idx)).c_str());
}
