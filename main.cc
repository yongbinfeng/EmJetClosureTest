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

	// Parse the args.
	cmd.parse( argc, argv );

	// Get the value parsed by each arg.
	string iconfig     = configArg.getValue();
	string ioutputDir  = outputDirArg.getValue();
        string isampleColl = sampleCollArg.getValue(); 
        std::cout << "Running on SampleCollection " << isampleColl << std::endl;
        string sampleCollDir = ioutputDir + "/" + isampleColl;
        std::cout << "Output directory set to : " << sampleCollDir << std::endl;
        time_t t = time(0);
        struct tm * now = localtime( & t );
        sampleCollDir += "/" + std::to_string(now->tm_mon+1) + std::to_string(now->tm_mday) + std::to_string(now->tm_hour);
        std::cout << "Creating directory: " << sampleCollDir << std::endl;
        string mkdir_command = "mkdir -p " + sampleCollDir;
        system(mkdir_command.c_str());

        // Main body of program
        EmJetSampleCollection ejsamplesColl;
        ReadFromConfigFile(iconfig, isampleColl, ejsamplesColl);
        std::cout << "Number of samples: " << ejsamplesColl.samples.size() << std::endl;
        PrintSampleCollection(ejsamplesColl);
    
        // Process files if there are any to be processed
        EmJetEventCount hm(ejsamplesColl);
        hm.OpenOutputFile(sampleCollDir+"/histo-result.root");
        std::cout << "file opened successfully "<< std::endl;
        hm.SetOptions(ejsamplesColl.isData);
        for(int i=0; i<500; i++){
          std::cout << " running on " << i << " time " << std::endl;
          hm.LoopOverCurrentTree();
        }
        hm.WriteHistograms();
        std::cout << "--------------------------finished--------------------------------\n";

  } catch (ArgException &e) {
    cerr << "error: " << e.error() << " for arg " << e.argId() << endl;
  }
}
