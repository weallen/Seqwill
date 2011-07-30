#include <vector>

#include <Eigen/Dense>
#include <gsl/gsl_rng.h>
#include <fstream>
#include "analysis/dist.h"
#include "common/Track.h"
#include "io/TrackIO.h"
#include "analysis/Random.h"

int main(int argc, char** argv) {
  std::string out(argv[1]);
  GaussDist g1(1.0, 1.0);
  GaussDist g2(9.0, 1.0);
  gsl_rng* rng = InitRng();
  std::vector<float> trueval(100000);
  std::vector<float> v(100000);
  int num_1 = 0;
  int num_2 = 0;
    
  Eigen::Matrix2f m;
  m << 0.3 , 0.7,
       0.7 , 0.3;
  
  int curr_state = 0;
  double s = 0.0;
  for (int i = 0; i < 100000; ++i) {
    if (curr_state == 0) {
      v[i] = g1.sample(rng);
      num_1++;
    }
    else if (curr_state == 1) {
      v[i] = g2.sample(rng);
      num_2++;
    }
    trueval[i] = curr_state;
    s = gsl_rng_uniform(rng);
    if (s < .0001) {
      if (curr_state == 0)
        curr_state = 1;
      else if (curr_state == 1)
        curr_state = 0;
    }
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
  std::cerr << "Num 1 " << num_1 << std::endl;
  std::cerr << "Num 2 " << num_2 << std::endl;
  
  std::ofstream fout;
  std::string outfilename = out + std::string(".txt");
  fout.open(outfilename.c_str());
  for (size_t k = 0; k < v.size(); ++k) {
    fout << trueval[k] << std::endl;
  }
  fout.close();
}
