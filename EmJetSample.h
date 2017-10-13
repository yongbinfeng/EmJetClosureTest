#include <string> // string, stof
#include <vector> // vector
#include <iostream> // cout
#include <fstream> // std::ifstream
#include <sstream> // std::stringstream
#include <algorithm>

using std::string;
using std::vector;
using std::getline;

struct EmJetSample_s
{
  string name;
  double xsec;
  string file;
};

typedef struct EmJetSample_s EmJetSample;

struct EmJetSampleCollection_s
{
  string name; // Unique name
  bool isData; // true for data
  vector<EmJetSample> samples;
};

typedef struct EmJetSampleCollection_s EmJetSampleCollection;

void PrintSampleCollection(const EmJetSampleCollection& samplesColl);
void ReadFromConfigFile(string configFileName, string isampleColl, EmJetSampleCollection& samplesColl);
void FileToLines(string ifileName, vector<string>& result);
void LineToFields(string iline, vector<string>& result);
void FieldsToSampleCollection(const vector<string>& ifields, EmJetSampleCollection& samplesColl);
void removeSpaces(string& input);
