#include "EmJetEventCount.h"
//#include "EmJetSample.h"
#include "tclap/CmdLine.h"

#include <iostream>
#include <string>
#include <map>
#include <algorithm>
#include <unistd.h>
#include <ctime>

using namespace TCLAP;
using std::cout;
using std::cerr;
using std::endl;

int main(int argc, char *argv[])
{
    try {
	      // Define the command line object.
	      CmdLine cmd("Run EmJetHistoMaker", ' ', "0.9");

	      // Define value arguments and add to the command line.
	      ValueArg<string> configArg("c","config","Config file",true,"config.txt","string");
	      cmd.add( configArg );

	      ValueArg<string> outputDirArg("d","directory","Output directory. Do NOT include trailing slash character",true,".","string");
      	cmd.add( outputDirArg );

        ValueArg<string> sampleCollArg("s","sampleColl","Name of sampleCollection to run over. Unspecified: Run over all.",false,"","string");
        cmd.add( sampleCollArg );

        ValueArg<int> repeatedTimeArg("r", "repeatedtime", "repeated time for smearing FR histograms and predicting number of backgrounds", false, 1, "int");
        cmd.add( repeatedTimeArg ); 

        ValueArg<string> outputLabelArg("l", "outputlabel", "output label of the root file", false, "", "string");
        cmd.add( outputLabelArg );

	      // Parse the args.
	      cmd.parse( argc, argv );

	      // Get the value parsed by each arg.
	      string iconfig     = configArg.getValue();
	      string ioutputDir  = outputDirArg.getValue();
        string isampleColl = sampleCollArg.getValue(); 
        int repeatedTime   = repeatedTimeArg.getValue(); 
        string outputLabel    = outputLabelArg.getValue();
        std::cout << "Running on SampleCollection " << isampleColl << std::endl;
        time_t t = time(0);
        struct tm * now = localtime( & t );
        string sampleCollDir = ioutputDir + "/" + std::to_string(now->tm_mon+1) + std::to_string(now->tm_mday) + std::to_string(now->tm_hour);
        std::cout << "Creating directory: " << sampleCollDir << " for output"<< std::endl;
        string mkdir_command = "mkdir -p " + sampleCollDir;
        system(mkdir_command.c_str());

        // Main body of program
        EmJetSampleCollection ejsamplesColl;
        ReadFromConfigFile(iconfig, isampleColl, ejsamplesColl);
        std::cout << "Number of samples: " << ejsamplesColl.samples.size() << std::endl;
        PrintSampleCollection(ejsamplesColl);
    
        // Process files if there are any to be processed
        EmJetEventCount hm(ejsamplesColl);
        hm.OpenOutputFile(sampleCollDir+"/histo-"+isampleColl+"_result_"+ outputLabel+".root");
        std::cout << "file opened successfully "<< std::endl;
        //string ffr = "/data/users/fengyb/ClosureTest/TestClosure/FRHisto/result_fakerate.root";
        string ffr = ejsamplesColl.FrCalfile;
        //vector<string> vhfr = {"fakerate_QCDMC_ptX", "fakerate_QCDMC_ptX_Lquark", "fakerate_QCDMC_ptX_Bquark", "fakerate_GJetMC_ptX", "fakerate_GJetMC_ptX_Lquark", "fakerate_GJetMC_ptX_Bquark", "FR_l_calc", "FR_b_calc", "fakerate_GJetMC_calc_1to2tag", "fakerate_QCDMC_truth_1to2tag"};
        //vector<string> vhfr = {"fakerate_QCDMC_ptX", "fakerate_QCDMC_ptX_Lquark", "fakerate_QCDMC_ptX_Bquark", "fakerate_GJetMC_ptX", "fakerate_GJetMC_ptX_Lquark", "fakerate_GJetMC_ptX_Bquark", "FR_l_calc", "FR_b_calc"};
        vector<string> vhfr = {"fakerate_QCDMC_ptX", "fakerate_QCDMC_ptX_Lquark", "fakerate_QCDMC_ptX_Bquark", "fakerate_GJetData_ptX", "fakerate_GJetMC_ptX_Lquark", "fakerate_GJetMC_ptX_Bquark", "FR_l_calc", "FR_b_calc", "fakerate_GJetData_calc_1to2tag", "fakerate_QCDMC_truth_1to2tag", "fakerate_GJetData_calc_0to1tag"};
        // set basic info for closure test
        hm.SetOptions(ffr, vhfr, ejsamplesColl.isData);
        TFile *ffrac = new TFile("/data/users/fengyb/frcal/fitdistribution/root/test_csv_GJet_5mm_Data_v1.root");
        hm.hfrac1_ = (TH1F*)(ffrac->Get("fraction_typeVII"));
        hm.hfrac2_ = (TH1F*)(ffrac->Get("fraction_typeVIII"));
        TFile *fFR = new TFile("/data/users/fengyb/frcal/gettruthinfo/root/result_fakerate_5mm_Data_v1.root");
        hm.hfr1_ = (TH1F*)(fFR->Get("fakerate_GJetData_ptX_TypeVII"));
        hm.hfr2_ = (TH1F*)(fFR->Get("fakerate_GJetData_ptX_TypeVIII"));
        hm.LoopOverTrees(repeatedTime);
        hm.WriteHistograms();
        std::cout << "--------------------------finished--------------------------------\n";

  } catch (ArgException &e) {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }
}
