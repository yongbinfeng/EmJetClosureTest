#include "EmJetSample.h"

using std::string;
using std::vector;
using std::getline;

void PrintSampleCollection(const EmJetSampleCollection& samplesColl)
{
  std::cout << samplesColl.name << ", " << samplesColl.isData << ", number of EmJetSample " << samplesColl.samples.size() << std::endl;
  for(EmJetSample sample: samplesColl.samples){
    std::cout << "   sample " << sample.name << " " << sample.xsec << " " << sample.file << std::endl;
  }
}

void
ReadFromConfigFile(string configFileName, string isampleColl, EmJetSampleCollection& samplesColl)
{
  std::cout << "Look for EmJetSampleCollection " << isampleColl << " from config file: " << configFileName << std::endl;
  vector<string> lines;
  FileToLines(configFileName, lines);
  if (lines.empty()) {
    std::cerr << "Config file empty!\n";
    return;
  }
  // Loop over sampleCollection configs
  for (string line : lines) {
    vector<string> fields;
    LineToFields(line, fields);
    if (fields.empty()) {
      // skip empty lines
      // std::cerr << "fields empty\n";
      continue;
    }
    else if( fields[0]==isampleColl ){
      std::cout << "Found sample from config file" << std::endl;
      FieldsToSampleCollection(fields, samplesColl);
      return;
    }
  }
  {
     std::cerr << " can not find the sample collection " << isampleColl << " from " << configFileName << std::endl;
     return;
  }
}

// Convert input file to vector of lines
void
FileToLines(string ifileName, vector<string>& result)
{
  std::ifstream lineStream(ifileName);
  string line;
  while(getline(lineStream,line)) {
    result.push_back(line);
  }
  return;
}

// Covert input line to a vector of fields
void
LineToFields(string iline, vector<string>& result)
{
  if (iline[0]=='#') {
    // std::cout << "Ignoring comment line" << std::endl;
    return; // Ignore lines that start with '#', return empty vector
  }

  // std::cout << "iline: " << iline << std::endl;
  removeSpaces(iline); // Remove spaces
  // std::cout << "iline after space removal: " << iline << std::endl;
  std::stringstream iss(iline);
  string field;
  // Read result separated by comma
  while (getline(iss, field, ',')) {
    result.push_back(field);
  }
  // This checks for a trailing comma with no data after it.
  if (!iss && field.empty()) {
    result.push_back(""); // If there was a trailing comma then add an empty element.
  }
}


// Turn vector of text fields into Sample object
void
FieldsToSampleCollection(const vector<string>& ifields, EmJetSampleCollection& samplesColl)
{
  samplesColl.name = ifields[0];

  if      (ifields[1]=="MC"                         ) samplesColl.isData = false;
  else if (ifields[1]=="DATA" or ifields[1]=="Data" ) samplesColl.isData = true;
  else std::cerr << "Sample MC/DATA not specified correctly in config file for sample: " << samplesColl.name << std::endl;

  bool found_txt = ( ifields[2].find(".txt") != string::npos );
  if ( found_txt ) {
    std::cout << " look into " << ifields[0] << " for the " << samplesColl.name << " collection" << std::endl;
    // look for samples, xsec and root files
    vector<string> lines; 
    FileToLines(ifields[2], lines);
    for (string line : lines) {
       vector<string> fields;
       LineToFields(line, fields);
       if (fields.empty()) {
       // skip empty lines
       // std::cerr << "fields empty\n";
          continue;
       }
       else {
         EmJetSample sample;
         sample.name = fields[0];
         sample.xsec = stof(fields[1]);
         if ( fields[2].find(".root") != string::npos ){
           sample.file = fields[2];
         }
         else{
           std::cout << " invalid input of " << samplesColl.name  << " "<< sample.name << " in " << ifields[3] << std::endl;
         }
         samplesColl.samples.push_back(sample);
       }
    }
  }
  else {
    std::cerr << " Invalid config file, input files not set" << std::endl;
  }
  bool found_root = ( ifields[3].find(".root") != string::npos );
  if ( found_root ){
    std::cout << " Fakerate histograms set to " << ifields[3] << std::endl;
    samplesColl.FrCalfile = ifields[3];
  }
  else{
    std::cerr << " Invalid config file, fakerate root file not set "<< std::endl;
  }
}

void removeSpaces(string& input)
{
  input.erase(std::remove(input.begin(),input.end(),' '),input.end());
}
