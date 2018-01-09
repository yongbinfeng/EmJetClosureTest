//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jan  7 23:52:18 2018 by ROOT version 6.06/01
// from TTree emJetSlimmedTree/emJetSlimmedTree
// found on file: ntuple_QCDHT1000To1500.root
//////////////////////////////////////////////////////////

#ifndef BaseClass_h
#define BaseClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
using std::vector;
#include "vector"
#include "vector"

class BaseClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           event;
   Int_t           run;
   Int_t           lumi;
   Int_t           nEmerging;
   Float_t         ht;
   Float_t         pv_genreco_disXY;
   Float_t         pv_genreco_disZ;
   Float_t         pvtrack_fraction;
   Float_t         pvtrack_fraction2;
   Float_t         pvtrackpt_fraction;
   Float_t         trackmeanz;
   Float_t         trackfabsmeanz;
   Float_t         pv0pt2sum;
   Float_t         pv1pt2sum;
   vector<int>     *jet_index;
   vector<float>   *jet_pt;
   vector<float>   *jet_eta;
   vector<float>   *jet_phi;
   vector<int>     *jet_nTrack;
   vector<int>     *jet_nTrackPostCut;
   vector<int>     *jet_flavour;
   vector<float>   *jet_Alpha3DSig;
   vector<float>   *jet_medianIP;
   vector<float>   *jet_theta2D;
   vector<float>   *jet_csv;
   vector<bool>    *jet_isEmerging;
   Float_t         pv_gen_x;
   Float_t         pv_gen_y;
   Float_t         pv_gen_z;
   Float_t         pv_reco_x;
   Float_t         pv_reco_y;
   Float_t         pv_reco_z;
   Float_t         pv_reco_zError;
   Float_t         pv_reco_chi2;
   Float_t         pv_reco_ndof;
   Int_t           pv_reco_nTracks;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_nEmerging;   //!
   TBranch        *b_ht;   //!
   TBranch        *b_pv_genreco_disXY;   //!
   TBranch        *b_pv_genreco_disZ;   //!
   TBranch        *b_pvtrack_fraction;   //!
   TBranch        *b_pvtrack_fraction2;   //!
   TBranch        *b_pvtrackpt_fraction;   //!
   TBranch        *b_trackmeanz;   //!
   TBranch        *b_trackfabsmeanz;   //!
   TBranch        *b_pv0pt2sum;   //!
   TBranch        *b_pv1pt2sum;   //!
   TBranch        *b_jet_index;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_nTrack;   //!
   TBranch        *b_jet_nTrackPostCut;   //!
   TBranch        *b_jet_flavour;   //!
   TBranch        *b_jet_Alpha3DSig;   //!
   TBranch        *b_jet_medianIP;   //!
   TBranch        *b_jet_theta2D;   //!
   TBranch        *b_jet_csv;   //!
   TBranch        *b_jet_isEmerging;   //!
   TBranch        *b_pv_gen_x;   //!
   TBranch        *b_pv_gen_y;   //!
   TBranch        *b_pv_gen_z;   //!
   TBranch        *b_pv_reco_x;   //!
   TBranch        *b_pv_reco_y;   //!
   TBranch        *b_pv_reco_z;   //!
   TBranch        *b_pv_reco_zError;   //!
   TBranch        *b_pv_reco_chi2;   //!
   TBranch        *b_pv_reco_ndof;   //!
   TBranch        *b_pv_reco_nTracks;   //!

   BaseClass(TTree *tree=0);
   virtual ~BaseClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef BaseClass_cxx
BaseClass::BaseClass(TTree *tree) : fChain(0) 
{
}

BaseClass::~BaseClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BaseClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BaseClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void BaseClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   jet_index = 0;
   jet_pt = 0;
   jet_eta = 0;
   jet_phi = 0;
   jet_nTrack = 0;
   jet_nTrackPostCut = 0;
   jet_flavour = 0;
   jet_Alpha3DSig = 0;
   jet_medianIP = 0;
   jet_theta2D = 0;
   jet_csv = 0;
   jet_isEmerging = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("nEmerging", &nEmerging, &b_nEmerging);
   fChain->SetBranchAddress("ht", &ht, &b_ht);
   fChain->SetBranchAddress("pv_genreco_disXY", &pv_genreco_disXY, &b_pv_genreco_disXY);
   fChain->SetBranchAddress("pv_genreco_disZ", &pv_genreco_disZ, &b_pv_genreco_disZ);
   fChain->SetBranchAddress("pvtrack_fraction", &pvtrack_fraction, &b_pvtrack_fraction);
   fChain->SetBranchAddress("pvtrack_fraction2", &pvtrack_fraction2, &b_pvtrack_fraction2);
   fChain->SetBranchAddress("pvtrackpt_fraction", &pvtrackpt_fraction, &b_pvtrackpt_fraction);
   fChain->SetBranchAddress("trackmeanz", &trackmeanz, &b_trackmeanz);
   fChain->SetBranchAddress("trackfabsmeanz", &trackfabsmeanz, &b_trackfabsmeanz);
   fChain->SetBranchAddress("pv0pt2sum", &pv0pt2sum, &b_pv0pt2sum);
   fChain->SetBranchAddress("pv1pt2sum", &pv1pt2sum, &b_pv1pt2sum);
   fChain->SetBranchAddress("jet_index", &jet_index, &b_jet_index);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_nTrack", &jet_nTrack, &b_jet_nTrack);
   fChain->SetBranchAddress("jet_nTrackPostCut", &jet_nTrackPostCut, &b_jet_nTrackPostCut);
   fChain->SetBranchAddress("jet_flavour", &jet_flavour, &b_jet_flavour);
   fChain->SetBranchAddress("jet_Alpha3DSig", &jet_Alpha3DSig, &b_jet_Alpha3DSig);
   fChain->SetBranchAddress("jet_medianIP", &jet_medianIP, &b_jet_medianIP);
   fChain->SetBranchAddress("jet_theta2D", &jet_theta2D, &b_jet_theta2D);
   fChain->SetBranchAddress("jet_csv", &jet_csv, &b_jet_csv);
   fChain->SetBranchAddress("jet_isEmerging", &jet_isEmerging, &b_jet_isEmerging);
   fChain->SetBranchAddress("pv_gen_x", &pv_gen_x, &b_pv_gen_x);
   fChain->SetBranchAddress("pv_gen_y", &pv_gen_y, &b_pv_gen_y);
   fChain->SetBranchAddress("pv_gen_z", &pv_gen_z, &b_pv_gen_z);
   fChain->SetBranchAddress("pv_reco_x", &pv_reco_x, &b_pv_reco_x);
   fChain->SetBranchAddress("pv_reco_y", &pv_reco_y, &b_pv_reco_y);
   fChain->SetBranchAddress("pv_reco_z", &pv_reco_z, &b_pv_reco_z);
   fChain->SetBranchAddress("pv_reco_zError", &pv_reco_zError, &b_pv_reco_zError);
   fChain->SetBranchAddress("pv_reco_chi2", &pv_reco_chi2, &b_pv_reco_chi2);
   fChain->SetBranchAddress("pv_reco_ndof", &pv_reco_ndof, &b_pv_reco_ndof);
   fChain->SetBranchAddress("pv_reco_nTracks", &pv_reco_nTracks, &b_pv_reco_nTracks);
   Notify();
}

Bool_t BaseClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BaseClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BaseClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BaseClass_cxx
