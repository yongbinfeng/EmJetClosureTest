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

        ValueArg<int> doPredictArg("p", "dopredict", "Predict the number of background or not", false, 0, "int");
        cmd.add( doPredictArg );

        ValueArg<int> doFillHistoArg("f", "dofillhisto", "Fill different histograms or not", false, 0, "int");
        cmd.add( doFillHistoArg );

	      // Parse the args.
	      cmd.parse( argc, argv );

	      // Get the value parsed by each arg.
	      string iconfig     = configArg.getValue();
	      string ioutputDir  = outputDirArg.getValue();
        string isampleColl = sampleCollArg.getValue(); 
        int repeatedTime   = repeatedTimeArg.getValue(); 
        string outputLabel = outputLabelArg.getValue();
        int doPredict      = doPredictArg.getValue();
        int doFillHisto    = doFillHistoArg.getValue();
        std::cout << "hahahahah " << doFillHisto << std::endl;
        std::cout << "Running on SampleCollection " << isampleColl << std::endl;
        time_t t = time(0);
        struct tm * now = localtime( & t );
        string sampleCollDir = ioutputDir;
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
        hm.OpenOutputFile(sampleCollDir+"/histo-"+isampleColl+"_result_"+ std::to_string(now->tm_mon+1) +std::to_string(now->tm_mday)+ "_" + outputLabel +".root");
        std::cout << "file opened successfully "<< std::endl;
        string ffr = ejsamplesColl.FrCalfile;
        vector<string> vhfr, v2frac, v2fr;
        if( ejsamplesColl.isData ){
          vhfr = {"Dfakerates/fakerate_GJetData"};
          v2frac = {"fraction_GJetData_TypeVI", "fraction_GJetData_TypeV"};
          v2fr   = {"fakerate_GJetData_TypeVI", "fakerate_GJetData_TypeV"};
        }
        else{
          vhfr = {"Dfakerates/fakerate_GJetMC"};
          v2frac = {"fraction_GJetMC_TypeVI", "fraction_GJetMC_TypeV"};
          v2fr   = {"fakerate_GJetMC_TypeVI", "fakerate_GJetMC_TypeV"};
        }
        string bfractag = "bfraction_in_tags";
        // set basic info for closure test
        hm.SetFillOption(doFillHisto);
        hm.SetPredictOption(doPredict);
        hm.SetOptions(ffr, vhfr, v2frac, v2fr, bfractag, ejsamplesColl.isData);
        hm.LoopOverTrees(repeatedTime);
        hm.WriteHistograms();
        std::cout << "--------------------------finished--------------------------------\n";

  } catch (ArgException &e) {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }
}
