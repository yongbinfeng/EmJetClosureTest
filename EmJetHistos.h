#include <string>
#include <unordered_map>
#include <vector>
#include <cassert>
#include <TH1F.h>

using std::string;
using std::unordered_map;
using std::vector;

class EmJetHistos
{
 public:
  EmJetHistos();
  ~EmJetHistos();
  int GetHistoIndex(string);
  TH1F* Get1D(string);
  unordered_map<string, TH1F*> hist1d;
  unordered_map<string, TH1D*> hist1d_double;
 private:
  void fillToVector();
  vector<TH1F*> vector_hist1d;
};

