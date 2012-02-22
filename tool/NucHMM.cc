#include <iostream>
#include <vector>

#include "base/CommandLineParser.h"
#include "base/FileUtil.h"
#include "base/Types.h"
#include "io/TrackIO.h"
#include "analysis/Genome.h"
#include "analysis/HMM.h"
#include "analysis/Dist.h"
#include "analysis/Random.h"

using namespace std;

int main(int argc, char** argv) {
	commandArg<string> oDataPath("-o", "Out trackfile name");
	commandArg<string> iTrack("-i", "In trackfile name");
	commandArg<string> tTrackname("-tout", "Out track name");
	commandArg<string> iTrackname("-tin", "In track name");

	commandLineParser P(argc, argv);
	P.SetDescription("Find nucleosomes with an HMM");
	P.registerArg(oDataPath);
	P.registerArg(iTrack);
	P.registerArg(tTrackname);
	P.registerArg(iTrackname);

	string in_trackfile = P.GetStringValueFor(iTrack);
	string out_trackfile = P.GetStringValueFor(oDataPath);
	string out_trackname = P.GetStringValueFor(tTrackname);
	string in_trackname = P.GetStringValueFor(iTrackname);

	TrackFile tin(in_trackfile);
	TrackFile tout(out_trackfile);

	// Assumes a 147 bp nucleosomes + 1 state modeling internuc

	// linker self-trans prob
	float p = 0.33;

	//
	//
	// Explanation of transition matrix construction
	// ---------------------------------------------
	// Models each base inside a nucleosome as 1 states, with prob 1 
	// transitioning to the next internucleosome base 
	// has prob p of transitioning to start of nucleosome state seq
	// and prob 1-p of transition to linker state
	// 
	//	_						           _
	//  |p   1-p 0   0  0  ... 0| 
	//  |0   0   1   0  0  ... 0|
	//  |0   0   0   1  0  ... 0|
	//  |      ...							|
	//  |1   0   0   0  0      0|
	//  -                      -
	Eigen::MatrixXd transmat = Eigen::MatrixXd::Zero(148,148);
	transmat(0,0) = p;	
	transmat(0,1) = 1-p;
	for (int i = 1; i < 147; ++i) {
		transmat(i,i+1) = 1.0;
	}
	transmat(148,0) = 1.0;


	vector<string> chrs = tin.GetSubTrackNames(in_trackname);
	vector<GaussDist> dists(2);

	// set up low and high emission dist
	dists[0] = GaussDist(0.0, 1.0);
	dists[1] = GaussDist(4.0, 1.0);
	for (size_t i = 0; i < chrs.size(); ++i) {
		string chrname = chrs[i];

		GaussHMM hmm(148);

		hmm.set_out_track_name(out_trackname);
		Track<float>::Ptr track(new Track<float>);
		tin.ReadSubTrack<float>(in_trackname, chrname, *track);
		hmm.set_input(track);	
		hmm.set_out_subtrack_name(chrname);
		hmm.set_emit(dists);	
		hmm.set_transition(transmat);
		hmm.Init();
		hmm.Compute();
		Track<int>::Ptr tr = hmm.output();
		cout << "Saving track " << tr->subtrackname() << endl;
		tout.WriteSubTrack<int>(*tr);
	}
}
