#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "base/CommandLineParser.h"
#include "base/Types.h"
#include "io/TrackIO.h"

#include <Eigen/QR>
#include <Eigen/SVD>
#include <Eigen/Dense>

using namespace std;



// XXX NOT FUNCTIONAL RIGHT NOW
int main(int argc, char** argv) {
    commandArg<string> oDataPath("-o", "Out trackfile name");
    commandArg<string> iTrackfile("-i", "In trackfile name");
    commandArg<string> tTrackName("-t", "In track name");
    commandArg<string> cCpGName("-c", "CpG track name");
    commandArg<int> fragLen("-l", "Number of surroundeding bins on either side to include in coupling value");
    commandArg<double> maxNorm("-max", "Max normalization ratio", 8.0);
    //commandArg<float> hyperparam("-hyper", "Hyperparameter for linear regression", 0.0);

    commandLineParser P(argc, argv);
    P.SetDescription("Normalize Dip bins by CpG content\nNOTE: CpG track and input track must have same number of elemnts");
    P.registerArg(oDataPath);
    P.registerArg(iTrackfile);
    P.registerArg(tTrackName);
    P.registerArg(cCpGName);
    P.registerArg(maxNorm);
    P.registerArg(fragLen);
    P.parse();

    string cpgname = P.GetStringValueFor(cCpGName);
    string in_file = P.GetStringValueFor(iTrackfile);
    string out_file = P.GetStringValueFor(oDataPath);
    string trackname = P.GetStringValueFor(tTrackName);
    double maxnorm = P.GetDoubleValueFor(maxNorm);
    int coupleWidth = P.GetIntValueFor(fragLen);

    TrackFile cpg(cpgname);
    TrackFile tio(in_file);

    // get the total size of the genome
    int total_size = 0;
    vector<string> chrs = cpg.GetSubTrackNames(string("mm9"));

    int k = 0;
   for (vector<string>::iterator it = chrs.begin(); it != chrs.end(); ++it) {
        TrackMetadata metadata = cpg.GetSubTrackMetadata(string("mm9"), *it);
        int s = metadata.GetIntMetadata(string("Stop"), 0);
        if (tio.HasSubTrack(trackname,*it) && k == 0) {
                total_size += s;
                cout << *it << " SIZE: " << s << endl;
         //       k++;
        }
    }

    cout << "TOTAL GENOME SIZE: " << total_size << endl;

    Eigen::VectorXf cpgs = Eigen::VectorXf::Zero(total_size);
    Eigen::VectorXf data = Eigen::VectorXf::Zero(total_size);
    Eigen::VectorXf coupling = Eigen::VectorXf::Zero(total_size);

    int i = 0;
    k = 0;
    for (vector<string>::iterator it = chrs.begin(); it != chrs.end(); ++it) {
        if (tio.HasSubTrack(trackname, *it) && k == 0) {
            cout << "Copying chrom " << *it << endl;
            Track<int>::Ptr currcpg(new Track<int>);
            Track<float>::Ptr currdata(new Track<float>);

            string currchr = *it;
            cpg.ReadSubTrack(string("mm9"), currchr, *currcpg);
            tio.ReadSubTrack(trackname, currchr, *currdata);
            for (size_t j = 0; j < currcpg->size(); ++j) {
                cpgs(i) = currcpg->get(j);
                data(i) = currdata->get(j);
                i++;
            }
        //k++;
        }
    }

    cout << "Computing coupling values" << endl;
    // Compute coupling values
    for (int i = coupleWidth; i < (total_size-2*coupleWidth); ++i) {
        int start = i - coupleWidth;
        coupling(i) = cpgs.segment(start, 2*coupleWidth).sum();
    }

    //
    // CALIBRATION
    //
    //compute mean coupling value
    cout << "Computing mean coupling values" << endl;
    int cutoff = 25;
    Eigen::ArrayXd sums = Eigen::ArrayXd::Zero(cutoff);
    Eigen::ArrayXd counts = Eigen::ArrayXd::Constant(cutoff,1);
    for (int j = 0; j < total_size; ++j) {
        int currcpg = coupling(j);
        if (currcpg < cutoff) {
            sums(currcpg) += (double)data(j);
            counts(currcpg) += 1.0;
        }
    }
    Eigen::ArrayXd means = sums / counts;

    for (int i = 0; i < counts.size(); ++i) {
        cout << means(i) << endl;
    }
    int maxidx = 0;
    for (int i = 2; i < (means.size()-2); ++i) {
        if ((means(i-1) < means(i) && means(i+1) < means(i)) &&
            (means(i-2) < means(i-1) && means(i+2) < means(i+1))) {
            cout << "MAX CPG VAL: " << i << endl;
            maxidx = i;
            break;
        }
    }

    int total_num = 0;
    for (int i = 0; i <= maxidx; ++i) {
        total_num += counts(i);
    }

    // Get data just below cutoff coupling value
    Eigen::VectorXd coupling_cutoff = Eigen::VectorXd::Zero(total_num - (maxidx+1));
    Eigen::VectorXd data_cutoff = Eigen::VectorXd::Zero(total_num - (maxidx+1));

    i = 0;
    for (int j = 0; j < total_size; ++j) {
        if (coupling(j) <= maxidx) {
            coupling_cutoff(i) = (double)coupling(j);
            data_cutoff(i) = data(j);
            i++;
        }
    }
    // Do regression
    // Solve Ax=b
    cout << "Estimating model parameters" << endl;
    //Eigen::MatrixXd x = coupling_cutoff.transpose().householderQr().solve(data_cutoff);
    //Eigen::JacobiSVD<Eigen::MatrixXd> svd(coupling_cutoff.transpose(), Eigen::ComputeThinU | Eigen::ComputeThinV);
    //Eigen::MatrixXd x = svd.solve(data_cutoff);
    // just do it by hand
    double xmean = coupling_cutoff.sum() / (double) coupling_cutoff.size();
   // cout << xmean << endl;
    double ymean = data_cutoff.sum() / (double) data_cutoff.size();
    // cout << ymean << endl;
    double xy = 0.0;
    double x2 = 0.0;
    for (int i = 0; i < data_cutoff.size(); ++i) {
        xy += data_cutoff(i) * coupling_cutoff(i);
        x2 += coupling_cutoff(i) * coupling_cutoff(i);
    }
//    double xy = data_cutoff.transpose() * coupling_cutoff;
  //  cout << xy << endl;
//    double x2 = coupling_cutoff.transpose() * coupling_cutoff;
 //   cout << x2 << endl;
    double b = (xy - xmean*ymean)/(x2 - xmean*xmean);
    double a = ymean - b * xmean;
//    cout << a << " " << b << endl;

    // Scale each bin by the expected enrichment and write out data
    TrackFile outio(out_file);
    for (vector<string>::iterator it = chrs.begin(); it != chrs.end(); ++it) {
        if (tio.HasSubTrack(trackname, *it)) {
            cout << "Saving chrom " << *it << endl;
            Track<int>::Ptr currcpg(new Track<int>);
            Track<float>::Ptr currdata(new Track<float>);
            string currchr = *it;
            cpg.ReadSubTrack(string("mm9"), currchr, *currcpg);
            tio.ReadSubTrack(trackname, currchr, *currdata);
            for (size_t j = 0; j < currcpg->size(); ++j) {
                float pred = a + b*((float)currcpg->get(j));
                float currval = currdata->get(j) / pred;
                if (currval > (float)maxnorm) {
                    currval = (float)maxnorm;
                }
                currdata->set(j, currval);
            }
            outio.WriteSubTrack(*currdata);
        }
    }

    return 1;
}
