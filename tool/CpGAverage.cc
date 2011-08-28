#include <string>
#include <vector>
#include <iostream>

#include "base/CommandLineParser.h"
#include "base/Types.h"
#include "base/FileUtil.h"
#include "io/TrackIO.h"
#include "analysis/Genome.h"
#include "analysis/MedipNormalize.h"


int main(int argc, char** argv) {
  commandArg<std::string> cCpG("-c", "CpG counts");
  commandArg<std::string> oDataPath("-o", "Out trackfile name");
	commandArg<std::string> tTrackName("-t", "Trackname");
	commandArg<std::string> iDataPath("-i", "In trackfile name");
	commandArg<int> nRes("-n", "Resolution", 100);

  commandLineParser P(argc, argv);
  P.SetDescription("Count CpGs in genome.");
  P.registerArg(cCpG);
  P.registerArg(oDataPath);
	P.registerArg(tTrackName);
	P.registerArg(iDataPath);
  P.registerArg(nRes);
  P.parse();

  std::string cpgpath = P.GetStringValueFor(cCpG);
  std::string out_trackfilename = P.GetStringValueFor(oDataPath);
  std::string in_trackfilename = P.GetStringValueFor(iDataPath);
  std::string trackname = P.GetStringValueFor(tTrackName);	
  int res = P.GetIntValueFor(nRes);


  if (!FileExists(in_trackfilename) || !FileExists(cpgpath)) {
    std::cerr << "Can't open one of the input files" << std::endl;
    return -1;
  }

  TrackFile cpgio(cpgpath);
  TrackFile inio(in_trackfilename);
  TrackFile outio(out_trackfilename);


  if (!inio.HasTrack(trackname)) {
    std::cerr << in_trackfilename << " doesn't have track " << trackname << std::endl;
    return -1;
  }
  std::vector<std::string> chrnames = inio.GetSubTrackNames(trackname);
  for (std::vector<std::string>::iterator it = chrnames.begin();
       it != chrnames.end(); ++it) {
    if (it->find(std::string("random")) == std::string::npos) {
      std::cerr << "Processing chromosome " << *it << std::endl;
      Track<int>::Ptr cpg(new Track<int>);
      Track<float>::Ptr out(new Track<float>);
      Track<float>::Ptr in(new Track<float>);
      cpgio.ReadSubTrack<int>(std::string("mm9"), *it, *cpg);
      inio.ReadSubTrack<float>(trackname, *it, *in);
      outio.ReadSubTrack<float>(trackname, *it, *out);
      for (size_t i = 0; i < in->size(); ++i) {
	float temp = in->get(i) / static_cast<float>(cpg->get(i));
	out->set(i, temp);
      }			
      outio.WriteSubTrack<float>(*out);
      cpg.reset();
      out.reset();
      in.reset();
    }
  }
}
