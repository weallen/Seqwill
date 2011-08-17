#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

#include <vector>
#include <string>

#include "gtest/gtest.h"

#include <Eigen/Dense>

#include "analysis/Dist.h"
using namespace std;
namespace {
  class DistTest : public ::testing::Test {
  protected:
    DistTest() {
      
    }

    virtual ~DistTest() {
    }

    virtual void SetUp() {
    }
    virtual void TearDown() {
    }

  };
      
  TEST_F(DistTest, MVGaussTest) {
    MVGaussDist m;

    Eigen::VectorXd mean1 = Eigen::VectorXd::Constant(10, 1.0);
    Eigen::VectorXd mean2 = Eigen::VectorXd::Constant(10, 2.0);
    Eigen::MatrixXd var = Eigen::MatrixXd::Identity(10,10);
    MVGaussDist m2(mean1, var);
    std::cerr << m2.Pdf(mean2) << std::endl;
  }
}//Namespace
int main(int argc, char** argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
