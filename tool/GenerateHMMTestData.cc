#include <vector>

#include <gsl/gsl_rng.h>
#include <fstream>
#include "analysis/dist.h"
#include "common/Track.h"
#include "io/TrackIO.h"
#include "analysis/Random.h"

int main(int argc, char** argv) {
  std::string out(argv[1]);
  GaussDist g1(0.0, 1.0);
  GaussDist g2(10.0, 1.0);
  gsl_rng* rng = InitRng();
  std::vector<float> trueval(100000);
  std::vector<float> v(100000);
  for (int i = 0; i < 500; ++i) {
    v[i] = g1.sample(rng);
    trueval[i] = 1;
  }
  for (int i = 500; i < 1000; ++i) {
    v[i] = g2.sample(rng);
    trueval[i] = 2;
  }
  for (int i = 1000; i < 5000; ++i) {
    v[i] = g1.sample(rng);
    trueval[i] = 1;
  }
  for (int i = 50000; i < 90000; ++i) {
    v[i] = g2.sample(rng);
    trueval[i] = 2;
  }
  for (int i = 90000; i < 100000; ++i) {
    v[i] = g1.sample(rng);
    trueval[i] = 1;
  }
  Track<float> t;
  t.set_resolution(1);
  t.set_extends(0, 100000);
  t.set_trackname(std::string("testdata"));
  t.set_subtrackname(std::string("test1"));
  for (size_t j = 0; j < v.size(); ++j) {
    t.set(j,v[j]);
  }
  TrackFile f;
  f.Open(out);
  f.WriteSubTrack<float>(t);
  
  std::ofstream fout;
  std::string outfilename = out + std::string(".txt");
  fout.open(outfilename.c_str());
  for (size_t k = 0; k < v.size(); ++k) {
    fout << trueval[k] << std::endl;
  }
  fout.close();
}
