//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Oct 12 10:27:36 2017 by ROOT version 6.06/01
// from TTree emJetSlimmedTree/emJetSlimmedTree
// found on file: /data/users/fengyb/80Xresult/sntuples_QCD-test5mm/QCDHT2000ToInf/ntuple-QCDMC-QCDHT2000ToInf-test5mm-3.root
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
   vector<bool>    *jet_isEmerging;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_nEmerging;   //!
   TBranch        *b_ht;   //!
   TBranch        *b_jet_index;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_nTrack;   //!
   TBranch        *b_jet_nTrackPostCut; //!
   TBranch        *b_jet_flavour;   //!
   TBranch        *b_jet_Alpha3DSig;   //!
   TBranch        *b_jet_medianIP;   //!
   TBranch        *b_jet_theta2D;   //!
   TBranch        *b_jet_isEmerging;   //!

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
   fChain->SetBranchAddress("jet_isEmerging", &jet_isEmerging, &b_jet_isEmerging);
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
