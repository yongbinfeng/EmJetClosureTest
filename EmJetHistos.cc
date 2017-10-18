#include "EmJetHistos.h"

using std::string;
using std::unordered_map;
using std::vector;

EmJetHistos::EmJetHistos()
{
  TH1::SetDefaultSumw2();

  /*[[[cog
    import histo_definitions as mod
    mod.generate_histo_map_init()
    ]]]*/
  hist1d["ht"] = new TH1F("ht", "ht" , 100, 0, 2500);
  hist1d["jet_pt"] = new TH1F("jet_pt", "jet_pt" , 150, 0, 1500);
  hist1d["jet_eta"] = new TH1F("jet_eta", "jet_eta" , 100, -5, 5);
  hist1d["jet_phi"] = new TH1F("jet_phi", "jet_phi" , 100, -5, 5);
  hist1d["jet_nTrack"] = new TH1F("jet_nTrack", "jet_nTrack" , 100, 0.0, 100);
  hist1d["nJet_tag"] = new TH1F("nJet_tag", "nJet_tag" , 10, 0, 10);
  hist1d["n2tag_0"] = new TH1F("n2tag_0", "n2tag_0" , 300, 0, 30);
  hist1d["n2tag_1"] = new TH1F("n2tag_1", "n2tag_1" , 300, 0, 30);
  hist1d["n2tag_2"] = new TH1F("n2tag_2", "n2tag_2" , 300, 0, 30);
  hist1d["n2tag_3"] = new TH1F("n2tag_3", "n2tag_3" , 300, 0, 30);
  hist1d["n2tag_4"] = new TH1F("n2tag_4", "n2tag_4" , 300, 0, 30);
  hist1d["n2tag_5"] = new TH1F("n2tag_5", "n2tag_5" , 300, 0, 30);
  hist1d["jet_pt__B"] = new TH1F("jet_pt__B", "jet_pt__B" , 150, 0, 1500);
  hist1d["jet_pt__L"] = new TH1F("jet_pt__L", "jet_pt__L" , 150, 0, 1500);
  hist1d["jet_eta__B"] = new TH1F("jet_eta__B", "jet_eta__B" , 100, -5, 5);
  hist1d["jet_eta__L"] = new TH1F("jet_eta__L", "jet_eta__L" , 100, -5, 5);
  hist1d["jet_phi__B"] = new TH1F("jet_phi__B", "jet_phi__B" , 100, -5, 5);
  hist1d["jet_phi__L"] = new TH1F("jet_phi__L", "jet_phi__L" , 100, -5, 5);
  hist1d["jet_nTrack__B"] = new TH1F("jet_nTrack__B", "jet_nTrack__B" , 100, 0.0, 100);
  hist1d["jet_nTrack__L"] = new TH1F("jet_nTrack__L", "jet_nTrack__L" , 100, 0.0, 100);
  hist1d["jet_pt__Emerging"] = new TH1F("jet_pt__Emerging", "jet_pt__Emerging" , 150, 0, 1500);
  hist1d["jet_pt__Standard"] = new TH1F("jet_pt__Standard", "jet_pt__Standard" , 150, 0, 1500);
  hist1d["jet_eta__Emerging"] = new TH1F("jet_eta__Emerging", "jet_eta__Emerging" , 100, -5, 5);
  hist1d["jet_eta__Standard"] = new TH1F("jet_eta__Standard", "jet_eta__Standard" , 100, -5, 5);
  hist1d["jet_phi__Emerging"] = new TH1F("jet_phi__Emerging", "jet_phi__Emerging" , 100, -5, 5);
  hist1d["jet_phi__Standard"] = new TH1F("jet_phi__Standard", "jet_phi__Standard" , 100, -5, 5);
  hist1d["jet_nTrack__Emerging"] = new TH1F("jet_nTrack__Emerging", "jet_nTrack__Emerging" , 100, 0.0, 100);
  hist1d["jet_nTrack__Standard"] = new TH1F("jet_nTrack__Standard", "jet_nTrack__Standard" , 100, 0.0, 100);
  hist1d["jet_pt__B__Emerging"] = new TH1F("jet_pt__B__Emerging", "jet_pt__B__Emerging" , 150, 0, 1500);
  hist1d["jet_pt__B__Standard"] = new TH1F("jet_pt__B__Standard", "jet_pt__B__Standard" , 150, 0, 1500);
  hist1d["jet_pt__L__Emerging"] = new TH1F("jet_pt__L__Emerging", "jet_pt__L__Emerging" , 150, 0, 1500);
  hist1d["jet_pt__L__Standard"] = new TH1F("jet_pt__L__Standard", "jet_pt__L__Standard" , 150, 0, 1500);
  hist1d["jet_eta__B__Emerging"] = new TH1F("jet_eta__B__Emerging", "jet_eta__B__Emerging" , 100, -5, 5);
  hist1d["jet_eta__B__Standard"] = new TH1F("jet_eta__B__Standard", "jet_eta__B__Standard" , 100, -5, 5);
  hist1d["jet_eta__L__Emerging"] = new TH1F("jet_eta__L__Emerging", "jet_eta__L__Emerging" , 100, -5, 5);
  hist1d["jet_eta__L__Standard"] = new TH1F("jet_eta__L__Standard", "jet_eta__L__Standard" , 100, -5, 5);
  hist1d["jet_phi__B__Emerging"] = new TH1F("jet_phi__B__Emerging", "jet_phi__B__Emerging" , 100, -5, 5);
  hist1d["jet_phi__B__Standard"] = new TH1F("jet_phi__B__Standard", "jet_phi__B__Standard" , 100, -5, 5);
  hist1d["jet_phi__L__Emerging"] = new TH1F("jet_phi__L__Emerging", "jet_phi__L__Emerging" , 100, -5, 5);
  hist1d["jet_phi__L__Standard"] = new TH1F("jet_phi__L__Standard", "jet_phi__L__Standard" , 100, -5, 5);
  hist1d["jet_nTrack__B__Emerging"] = new TH1F("jet_nTrack__B__Emerging", "jet_nTrack__B__Emerging" , 100, 0.0, 100);
  hist1d["jet_nTrack__B__Standard"] = new TH1F("jet_nTrack__B__Standard", "jet_nTrack__B__Standard" , 100, 0.0, 100);
  hist1d["jet_nTrack__L__Emerging"] = new TH1F("jet_nTrack__L__Emerging", "jet_nTrack__L__Emerging" , 100, 0.0, 100);
  hist1d["jet_nTrack__L__Standard"] = new TH1F("jet_nTrack__L__Standard", "jet_nTrack__L__Standard" , 100, 0.0, 100);
  hist1d["jet_pt__0tag"] = new TH1F("jet_pt__0tag", "jet_pt__0tag" , 150, 0, 1500);
  hist1d["jet_pt__1tag"] = new TH1F("jet_pt__1tag", "jet_pt__1tag" , 150, 0, 1500);
  hist1d["jet_pt__2tag"] = new TH1F("jet_pt__2tag", "jet_pt__2tag" , 150, 0, 1500);
  hist1d["jet_eta__0tag"] = new TH1F("jet_eta__0tag", "jet_eta__0tag" , 100, -5, 5);
  hist1d["jet_eta__1tag"] = new TH1F("jet_eta__1tag", "jet_eta__1tag" , 100, -5, 5);
  hist1d["jet_eta__2tag"] = new TH1F("jet_eta__2tag", "jet_eta__2tag" , 100, -5, 5);
  hist1d["jet_phi__0tag"] = new TH1F("jet_phi__0tag", "jet_phi__0tag" , 100, -5, 5);
  hist1d["jet_phi__1tag"] = new TH1F("jet_phi__1tag", "jet_phi__1tag" , 100, -5, 5);
  hist1d["jet_phi__2tag"] = new TH1F("jet_phi__2tag", "jet_phi__2tag" , 100, -5, 5);
  hist1d["jet_nTrack__0tag"] = new TH1F("jet_nTrack__0tag", "jet_nTrack__0tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__1tag"] = new TH1F("jet_nTrack__1tag", "jet_nTrack__1tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__2tag"] = new TH1F("jet_nTrack__2tag", "jet_nTrack__2tag" , 100, 0.0, 100);
  hist1d["jet_pt__B__0tag"] = new TH1F("jet_pt__B__0tag", "jet_pt__B__0tag" , 150, 0, 1500);
  hist1d["jet_pt__B__1tag"] = new TH1F("jet_pt__B__1tag", "jet_pt__B__1tag" , 150, 0, 1500);
  hist1d["jet_pt__B__2tag"] = new TH1F("jet_pt__B__2tag", "jet_pt__B__2tag" , 150, 0, 1500);
  hist1d["jet_pt__L__0tag"] = new TH1F("jet_pt__L__0tag", "jet_pt__L__0tag" , 150, 0, 1500);
  hist1d["jet_pt__L__1tag"] = new TH1F("jet_pt__L__1tag", "jet_pt__L__1tag" , 150, 0, 1500);
  hist1d["jet_pt__L__2tag"] = new TH1F("jet_pt__L__2tag", "jet_pt__L__2tag" , 150, 0, 1500);
  hist1d["jet_eta__B__0tag"] = new TH1F("jet_eta__B__0tag", "jet_eta__B__0tag" , 100, -5, 5);
  hist1d["jet_eta__B__1tag"] = new TH1F("jet_eta__B__1tag", "jet_eta__B__1tag" , 100, -5, 5);
  hist1d["jet_eta__B__2tag"] = new TH1F("jet_eta__B__2tag", "jet_eta__B__2tag" , 100, -5, 5);
  hist1d["jet_eta__L__0tag"] = new TH1F("jet_eta__L__0tag", "jet_eta__L__0tag" , 100, -5, 5);
  hist1d["jet_eta__L__1tag"] = new TH1F("jet_eta__L__1tag", "jet_eta__L__1tag" , 100, -5, 5);
  hist1d["jet_eta__L__2tag"] = new TH1F("jet_eta__L__2tag", "jet_eta__L__2tag" , 100, -5, 5);
  hist1d["jet_phi__B__0tag"] = new TH1F("jet_phi__B__0tag", "jet_phi__B__0tag" , 100, -5, 5);
  hist1d["jet_phi__B__1tag"] = new TH1F("jet_phi__B__1tag", "jet_phi__B__1tag" , 100, -5, 5);
  hist1d["jet_phi__B__2tag"] = new TH1F("jet_phi__B__2tag", "jet_phi__B__2tag" , 100, -5, 5);
  hist1d["jet_phi__L__0tag"] = new TH1F("jet_phi__L__0tag", "jet_phi__L__0tag" , 100, -5, 5);
  hist1d["jet_phi__L__1tag"] = new TH1F("jet_phi__L__1tag", "jet_phi__L__1tag" , 100, -5, 5);
  hist1d["jet_phi__L__2tag"] = new TH1F("jet_phi__L__2tag", "jet_phi__L__2tag" , 100, -5, 5);
  hist1d["jet_nTrack__B__0tag"] = new TH1F("jet_nTrack__B__0tag", "jet_nTrack__B__0tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__B__1tag"] = new TH1F("jet_nTrack__B__1tag", "jet_nTrack__B__1tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__B__2tag"] = new TH1F("jet_nTrack__B__2tag", "jet_nTrack__B__2tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__L__0tag"] = new TH1F("jet_nTrack__L__0tag", "jet_nTrack__L__0tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__L__1tag"] = new TH1F("jet_nTrack__L__1tag", "jet_nTrack__L__1tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__L__2tag"] = new TH1F("jet_nTrack__L__2tag", "jet_nTrack__L__2tag" , 100, 0.0, 100);
  hist1d["jet_pt__Emerging__0tag"] = new TH1F("jet_pt__Emerging__0tag", "jet_pt__Emerging__0tag" , 150, 0, 1500);
  hist1d["jet_pt__Emerging__1tag"] = new TH1F("jet_pt__Emerging__1tag", "jet_pt__Emerging__1tag" , 150, 0, 1500);
  hist1d["jet_pt__Emerging__2tag"] = new TH1F("jet_pt__Emerging__2tag", "jet_pt__Emerging__2tag" , 150, 0, 1500);
  hist1d["jet_pt__Standard__0tag"] = new TH1F("jet_pt__Standard__0tag", "jet_pt__Standard__0tag" , 150, 0, 1500);
  hist1d["jet_pt__Standard__1tag"] = new TH1F("jet_pt__Standard__1tag", "jet_pt__Standard__1tag" , 150, 0, 1500);
  hist1d["jet_pt__Standard__2tag"] = new TH1F("jet_pt__Standard__2tag", "jet_pt__Standard__2tag" , 150, 0, 1500);
  hist1d["jet_eta__Emerging__0tag"] = new TH1F("jet_eta__Emerging__0tag", "jet_eta__Emerging__0tag" , 100, -5, 5);
  hist1d["jet_eta__Emerging__1tag"] = new TH1F("jet_eta__Emerging__1tag", "jet_eta__Emerging__1tag" , 100, -5, 5);
  hist1d["jet_eta__Emerging__2tag"] = new TH1F("jet_eta__Emerging__2tag", "jet_eta__Emerging__2tag" , 100, -5, 5);
  hist1d["jet_eta__Standard__0tag"] = new TH1F("jet_eta__Standard__0tag", "jet_eta__Standard__0tag" , 100, -5, 5);
  hist1d["jet_eta__Standard__1tag"] = new TH1F("jet_eta__Standard__1tag", "jet_eta__Standard__1tag" , 100, -5, 5);
  hist1d["jet_eta__Standard__2tag"] = new TH1F("jet_eta__Standard__2tag", "jet_eta__Standard__2tag" , 100, -5, 5);
  hist1d["jet_phi__Emerging__0tag"] = new TH1F("jet_phi__Emerging__0tag", "jet_phi__Emerging__0tag" , 100, -5, 5);
  hist1d["jet_phi__Emerging__1tag"] = new TH1F("jet_phi__Emerging__1tag", "jet_phi__Emerging__1tag" , 100, -5, 5);
  hist1d["jet_phi__Emerging__2tag"] = new TH1F("jet_phi__Emerging__2tag", "jet_phi__Emerging__2tag" , 100, -5, 5);
  hist1d["jet_phi__Standard__0tag"] = new TH1F("jet_phi__Standard__0tag", "jet_phi__Standard__0tag" , 100, -5, 5);
  hist1d["jet_phi__Standard__1tag"] = new TH1F("jet_phi__Standard__1tag", "jet_phi__Standard__1tag" , 100, -5, 5);
  hist1d["jet_phi__Standard__2tag"] = new TH1F("jet_phi__Standard__2tag", "jet_phi__Standard__2tag" , 100, -5, 5);
  hist1d["jet_nTrack__Emerging__0tag"] = new TH1F("jet_nTrack__Emerging__0tag", "jet_nTrack__Emerging__0tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__Emerging__1tag"] = new TH1F("jet_nTrack__Emerging__1tag", "jet_nTrack__Emerging__1tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__Emerging__2tag"] = new TH1F("jet_nTrack__Emerging__2tag", "jet_nTrack__Emerging__2tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__Standard__0tag"] = new TH1F("jet_nTrack__Standard__0tag", "jet_nTrack__Standard__0tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__Standard__1tag"] = new TH1F("jet_nTrack__Standard__1tag", "jet_nTrack__Standard__1tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__Standard__2tag"] = new TH1F("jet_nTrack__Standard__2tag", "jet_nTrack__Standard__2tag" , 100, 0.0, 100);
  hist1d["jet_pt__B__Emerging__0tag"] = new TH1F("jet_pt__B__Emerging__0tag", "jet_pt__B__Emerging__0tag" , 150, 0, 1500);
  hist1d["jet_pt__B__Emerging__1tag"] = new TH1F("jet_pt__B__Emerging__1tag", "jet_pt__B__Emerging__1tag" , 150, 0, 1500);
  hist1d["jet_pt__B__Emerging__2tag"] = new TH1F("jet_pt__B__Emerging__2tag", "jet_pt__B__Emerging__2tag" , 150, 0, 1500);
  hist1d["jet_pt__B__Standard__0tag"] = new TH1F("jet_pt__B__Standard__0tag", "jet_pt__B__Standard__0tag" , 150, 0, 1500);
  hist1d["jet_pt__B__Standard__1tag"] = new TH1F("jet_pt__B__Standard__1tag", "jet_pt__B__Standard__1tag" , 150, 0, 1500);
  hist1d["jet_pt__B__Standard__2tag"] = new TH1F("jet_pt__B__Standard__2tag", "jet_pt__B__Standard__2tag" , 150, 0, 1500);
  hist1d["jet_pt__L__Emerging__0tag"] = new TH1F("jet_pt__L__Emerging__0tag", "jet_pt__L__Emerging__0tag" , 150, 0, 1500);
  hist1d["jet_pt__L__Emerging__1tag"] = new TH1F("jet_pt__L__Emerging__1tag", "jet_pt__L__Emerging__1tag" , 150, 0, 1500);
  hist1d["jet_pt__L__Emerging__2tag"] = new TH1F("jet_pt__L__Emerging__2tag", "jet_pt__L__Emerging__2tag" , 150, 0, 1500);
  hist1d["jet_pt__L__Standard__0tag"] = new TH1F("jet_pt__L__Standard__0tag", "jet_pt__L__Standard__0tag" , 150, 0, 1500);
  hist1d["jet_pt__L__Standard__1tag"] = new TH1F("jet_pt__L__Standard__1tag", "jet_pt__L__Standard__1tag" , 150, 0, 1500);
  hist1d["jet_pt__L__Standard__2tag"] = new TH1F("jet_pt__L__Standard__2tag", "jet_pt__L__Standard__2tag" , 150, 0, 1500);
  hist1d["jet_eta__B__Emerging__0tag"] = new TH1F("jet_eta__B__Emerging__0tag", "jet_eta__B__Emerging__0tag" , 100, -5, 5);
  hist1d["jet_eta__B__Emerging__1tag"] = new TH1F("jet_eta__B__Emerging__1tag", "jet_eta__B__Emerging__1tag" , 100, -5, 5);
  hist1d["jet_eta__B__Emerging__2tag"] = new TH1F("jet_eta__B__Emerging__2tag", "jet_eta__B__Emerging__2tag" , 100, -5, 5);
  hist1d["jet_eta__B__Standard__0tag"] = new TH1F("jet_eta__B__Standard__0tag", "jet_eta__B__Standard__0tag" , 100, -5, 5);
  hist1d["jet_eta__B__Standard__1tag"] = new TH1F("jet_eta__B__Standard__1tag", "jet_eta__B__Standard__1tag" , 100, -5, 5);
  hist1d["jet_eta__B__Standard__2tag"] = new TH1F("jet_eta__B__Standard__2tag", "jet_eta__B__Standard__2tag" , 100, -5, 5);
  hist1d["jet_eta__L__Emerging__0tag"] = new TH1F("jet_eta__L__Emerging__0tag", "jet_eta__L__Emerging__0tag" , 100, -5, 5);
  hist1d["jet_eta__L__Emerging__1tag"] = new TH1F("jet_eta__L__Emerging__1tag", "jet_eta__L__Emerging__1tag" , 100, -5, 5);
  hist1d["jet_eta__L__Emerging__2tag"] = new TH1F("jet_eta__L__Emerging__2tag", "jet_eta__L__Emerging__2tag" , 100, -5, 5);
  hist1d["jet_eta__L__Standard__0tag"] = new TH1F("jet_eta__L__Standard__0tag", "jet_eta__L__Standard__0tag" , 100, -5, 5);
  hist1d["jet_eta__L__Standard__1tag"] = new TH1F("jet_eta__L__Standard__1tag", "jet_eta__L__Standard__1tag" , 100, -5, 5);
  hist1d["jet_eta__L__Standard__2tag"] = new TH1F("jet_eta__L__Standard__2tag", "jet_eta__L__Standard__2tag" , 100, -5, 5);
  hist1d["jet_phi__B__Emerging__0tag"] = new TH1F("jet_phi__B__Emerging__0tag", "jet_phi__B__Emerging__0tag" , 100, -5, 5);
  hist1d["jet_phi__B__Emerging__1tag"] = new TH1F("jet_phi__B__Emerging__1tag", "jet_phi__B__Emerging__1tag" , 100, -5, 5);
  hist1d["jet_phi__B__Emerging__2tag"] = new TH1F("jet_phi__B__Emerging__2tag", "jet_phi__B__Emerging__2tag" , 100, -5, 5);
  hist1d["jet_phi__B__Standard__0tag"] = new TH1F("jet_phi__B__Standard__0tag", "jet_phi__B__Standard__0tag" , 100, -5, 5);
  hist1d["jet_phi__B__Standard__1tag"] = new TH1F("jet_phi__B__Standard__1tag", "jet_phi__B__Standard__1tag" , 100, -5, 5);
  hist1d["jet_phi__B__Standard__2tag"] = new TH1F("jet_phi__B__Standard__2tag", "jet_phi__B__Standard__2tag" , 100, -5, 5);
  hist1d["jet_phi__L__Emerging__0tag"] = new TH1F("jet_phi__L__Emerging__0tag", "jet_phi__L__Emerging__0tag" , 100, -5, 5);
  hist1d["jet_phi__L__Emerging__1tag"] = new TH1F("jet_phi__L__Emerging__1tag", "jet_phi__L__Emerging__1tag" , 100, -5, 5);
  hist1d["jet_phi__L__Emerging__2tag"] = new TH1F("jet_phi__L__Emerging__2tag", "jet_phi__L__Emerging__2tag" , 100, -5, 5);
  hist1d["jet_phi__L__Standard__0tag"] = new TH1F("jet_phi__L__Standard__0tag", "jet_phi__L__Standard__0tag" , 100, -5, 5);
  hist1d["jet_phi__L__Standard__1tag"] = new TH1F("jet_phi__L__Standard__1tag", "jet_phi__L__Standard__1tag" , 100, -5, 5);
  hist1d["jet_phi__L__Standard__2tag"] = new TH1F("jet_phi__L__Standard__2tag", "jet_phi__L__Standard__2tag" , 100, -5, 5);
  hist1d["jet_nTrack__B__Emerging__0tag"] = new TH1F("jet_nTrack__B__Emerging__0tag", "jet_nTrack__B__Emerging__0tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__B__Emerging__1tag"] = new TH1F("jet_nTrack__B__Emerging__1tag", "jet_nTrack__B__Emerging__1tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__B__Emerging__2tag"] = new TH1F("jet_nTrack__B__Emerging__2tag", "jet_nTrack__B__Emerging__2tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__B__Standard__0tag"] = new TH1F("jet_nTrack__B__Standard__0tag", "jet_nTrack__B__Standard__0tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__B__Standard__1tag"] = new TH1F("jet_nTrack__B__Standard__1tag", "jet_nTrack__B__Standard__1tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__B__Standard__2tag"] = new TH1F("jet_nTrack__B__Standard__2tag", "jet_nTrack__B__Standard__2tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__L__Emerging__0tag"] = new TH1F("jet_nTrack__L__Emerging__0tag", "jet_nTrack__L__Emerging__0tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__L__Emerging__1tag"] = new TH1F("jet_nTrack__L__Emerging__1tag", "jet_nTrack__L__Emerging__1tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__L__Emerging__2tag"] = new TH1F("jet_nTrack__L__Emerging__2tag", "jet_nTrack__L__Emerging__2tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__L__Standard__0tag"] = new TH1F("jet_nTrack__L__Standard__0tag", "jet_nTrack__L__Standard__0tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__L__Standard__1tag"] = new TH1F("jet_nTrack__L__Standard__1tag", "jet_nTrack__L__Standard__1tag" , 100, 0.0, 100);
  hist1d["jet_nTrack__L__Standard__2tag"] = new TH1F("jet_nTrack__L__Standard__2tag", "jet_nTrack__L__Standard__2tag" , 100, 0.0, 100);
    //[[[end]]]
}

EmJetHistos::~EmJetHistos()
{
  for (const auto & kv : hist1d) {
    delete kv.second;
  }
}

int
EmJetHistos::GetHistoIndex(string name)
{
  /*[[[cog
    mod.generate_histo_index()
    ]]]*/
  if (name=="ht") return 0;
  if (name=="jet_pt") return 1;
  if (name=="jet_eta") return 2;
  if (name=="jet_phi") return 3;
  if (name=="jet_nTrack") return 4;
  if (name=="nJet_tag") return 5;
  if (name=="n2tag_0") return 6;
  if (name=="n2tag_1") return 7;
  if (name=="n2tag_2") return 8;
  if (name=="n2tag_3") return 9;
  if (name=="n2tag_4") return 10;
  if (name=="n2tag_5") return 11;
  if (name=="jet_pt__B") return 12;
  if (name=="jet_pt__L") return 13;
  if (name=="jet_eta__B") return 14;
  if (name=="jet_eta__L") return 15;
  if (name=="jet_phi__B") return 16;
  if (name=="jet_phi__L") return 17;
  if (name=="jet_nTrack__B") return 18;
  if (name=="jet_nTrack__L") return 19;
  if (name=="jet_pt__Emerging") return 20;
  if (name=="jet_pt__Standard") return 21;
  if (name=="jet_eta__Emerging") return 22;
  if (name=="jet_eta__Standard") return 23;
  if (name=="jet_phi__Emerging") return 24;
  if (name=="jet_phi__Standard") return 25;
  if (name=="jet_nTrack__Emerging") return 26;
  if (name=="jet_nTrack__Standard") return 27;
  if (name=="jet_pt__B__Emerging") return 28;
  if (name=="jet_pt__B__Standard") return 29;
  if (name=="jet_pt__L__Emerging") return 30;
  if (name=="jet_pt__L__Standard") return 31;
  if (name=="jet_eta__B__Emerging") return 32;
  if (name=="jet_eta__B__Standard") return 33;
  if (name=="jet_eta__L__Emerging") return 34;
  if (name=="jet_eta__L__Standard") return 35;
  if (name=="jet_phi__B__Emerging") return 36;
  if (name=="jet_phi__B__Standard") return 37;
  if (name=="jet_phi__L__Emerging") return 38;
  if (name=="jet_phi__L__Standard") return 39;
  if (name=="jet_nTrack__B__Emerging") return 40;
  if (name=="jet_nTrack__B__Standard") return 41;
  if (name=="jet_nTrack__L__Emerging") return 42;
  if (name=="jet_nTrack__L__Standard") return 43;
  if (name=="jet_pt__0tag") return 44;
  if (name=="jet_pt__1tag") return 45;
  if (name=="jet_pt__2tag") return 46;
  if (name=="jet_eta__0tag") return 47;
  if (name=="jet_eta__1tag") return 48;
  if (name=="jet_eta__2tag") return 49;
  if (name=="jet_phi__0tag") return 50;
  if (name=="jet_phi__1tag") return 51;
  if (name=="jet_phi__2tag") return 52;
  if (name=="jet_nTrack__0tag") return 53;
  if (name=="jet_nTrack__1tag") return 54;
  if (name=="jet_nTrack__2tag") return 55;
  if (name=="jet_pt__B__0tag") return 56;
  if (name=="jet_pt__B__1tag") return 57;
  if (name=="jet_pt__B__2tag") return 58;
  if (name=="jet_pt__L__0tag") return 59;
  if (name=="jet_pt__L__1tag") return 60;
  if (name=="jet_pt__L__2tag") return 61;
  if (name=="jet_eta__B__0tag") return 62;
  if (name=="jet_eta__B__1tag") return 63;
  if (name=="jet_eta__B__2tag") return 64;
  if (name=="jet_eta__L__0tag") return 65;
  if (name=="jet_eta__L__1tag") return 66;
  if (name=="jet_eta__L__2tag") return 67;
  if (name=="jet_phi__B__0tag") return 68;
  if (name=="jet_phi__B__1tag") return 69;
  if (name=="jet_phi__B__2tag") return 70;
  if (name=="jet_phi__L__0tag") return 71;
  if (name=="jet_phi__L__1tag") return 72;
  if (name=="jet_phi__L__2tag") return 73;
  if (name=="jet_nTrack__B__0tag") return 74;
  if (name=="jet_nTrack__B__1tag") return 75;
  if (name=="jet_nTrack__B__2tag") return 76;
  if (name=="jet_nTrack__L__0tag") return 77;
  if (name=="jet_nTrack__L__1tag") return 78;
  if (name=="jet_nTrack__L__2tag") return 79;
  if (name=="jet_pt__Emerging__0tag") return 80;
  if (name=="jet_pt__Emerging__1tag") return 81;
  if (name=="jet_pt__Emerging__2tag") return 82;
  if (name=="jet_pt__Standard__0tag") return 83;
  if (name=="jet_pt__Standard__1tag") return 84;
  if (name=="jet_pt__Standard__2tag") return 85;
  if (name=="jet_eta__Emerging__0tag") return 86;
  if (name=="jet_eta__Emerging__1tag") return 87;
  if (name=="jet_eta__Emerging__2tag") return 88;
  if (name=="jet_eta__Standard__0tag") return 89;
  if (name=="jet_eta__Standard__1tag") return 90;
  if (name=="jet_eta__Standard__2tag") return 91;
  if (name=="jet_phi__Emerging__0tag") return 92;
  if (name=="jet_phi__Emerging__1tag") return 93;
  if (name=="jet_phi__Emerging__2tag") return 94;
  if (name=="jet_phi__Standard__0tag") return 95;
  if (name=="jet_phi__Standard__1tag") return 96;
  if (name=="jet_phi__Standard__2tag") return 97;
  if (name=="jet_nTrack__Emerging__0tag") return 98;
  if (name=="jet_nTrack__Emerging__1tag") return 99;
  if (name=="jet_nTrack__Emerging__2tag") return 100;
  if (name=="jet_nTrack__Standard__0tag") return 101;
  if (name=="jet_nTrack__Standard__1tag") return 102;
  if (name=="jet_nTrack__Standard__2tag") return 103;
  if (name=="jet_pt__B__Emerging__0tag") return 104;
  if (name=="jet_pt__B__Emerging__1tag") return 105;
  if (name=="jet_pt__B__Emerging__2tag") return 106;
  if (name=="jet_pt__B__Standard__0tag") return 107;
  if (name=="jet_pt__B__Standard__1tag") return 108;
  if (name=="jet_pt__B__Standard__2tag") return 109;
  if (name=="jet_pt__L__Emerging__0tag") return 110;
  if (name=="jet_pt__L__Emerging__1tag") return 111;
  if (name=="jet_pt__L__Emerging__2tag") return 112;
  if (name=="jet_pt__L__Standard__0tag") return 113;
  if (name=="jet_pt__L__Standard__1tag") return 114;
  if (name=="jet_pt__L__Standard__2tag") return 115;
  if (name=="jet_eta__B__Emerging__0tag") return 116;
  if (name=="jet_eta__B__Emerging__1tag") return 117;
  if (name=="jet_eta__B__Emerging__2tag") return 118;
  if (name=="jet_eta__B__Standard__0tag") return 119;
  if (name=="jet_eta__B__Standard__1tag") return 120;
  if (name=="jet_eta__B__Standard__2tag") return 121;
  if (name=="jet_eta__L__Emerging__0tag") return 122;
  if (name=="jet_eta__L__Emerging__1tag") return 123;
  if (name=="jet_eta__L__Emerging__2tag") return 124;
  if (name=="jet_eta__L__Standard__0tag") return 125;
  if (name=="jet_eta__L__Standard__1tag") return 126;
  if (name=="jet_eta__L__Standard__2tag") return 127;
  if (name=="jet_phi__B__Emerging__0tag") return 128;
  if (name=="jet_phi__B__Emerging__1tag") return 129;
  if (name=="jet_phi__B__Emerging__2tag") return 130;
  if (name=="jet_phi__B__Standard__0tag") return 131;
  if (name=="jet_phi__B__Standard__1tag") return 132;
  if (name=="jet_phi__B__Standard__2tag") return 133;
  if (name=="jet_phi__L__Emerging__0tag") return 134;
  if (name=="jet_phi__L__Emerging__1tag") return 135;
  if (name=="jet_phi__L__Emerging__2tag") return 136;
  if (name=="jet_phi__L__Standard__0tag") return 137;
  if (name=="jet_phi__L__Standard__1tag") return 138;
  if (name=="jet_phi__L__Standard__2tag") return 139;
  if (name=="jet_nTrack__B__Emerging__0tag") return 140;
  if (name=="jet_nTrack__B__Emerging__1tag") return 141;
  if (name=="jet_nTrack__B__Emerging__2tag") return 142;
  if (name=="jet_nTrack__B__Standard__0tag") return 143;
  if (name=="jet_nTrack__B__Standard__1tag") return 144;
  if (name=="jet_nTrack__B__Standard__2tag") return 145;
  if (name=="jet_nTrack__L__Emerging__0tag") return 146;
  if (name=="jet_nTrack__L__Emerging__1tag") return 147;
  if (name=="jet_nTrack__L__Emerging__2tag") return 148;
  if (name=="jet_nTrack__L__Standard__0tag") return 149;
  if (name=="jet_nTrack__L__Standard__1tag") return 150;
  if (name=="jet_nTrack__L__Standard__2tag") return 151;
  //[[[end]]]
  
  return -1;
}

TH1F*
EmJetHistos::Get1D(string name)
{
  unsigned index = GetHistoIndex(name);
  assert(vector_hist1d.size()>index); // Make sure index exists first
  auto histPointer = vector_hist1d[index];
  assert(histPointer);
  return NULL;
}

void
EmJetHistos::fillToVector()
{
  vector_hist1d.reserve(hist1d.size());
  for (const auto & kv : hist1d) {
    string name = kv.first;
    auto histPointer = kv.second;
    assert( histPointer ); // Make sure pointer is not null
    assert( vector_hist1d.size() == unsigned(GetHistoIndex(name)) ); // Make sure histogram index is correct
    vector_hist1d.push_back(histPointer);
  }
}
