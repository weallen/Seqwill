#include "base/MatrixUtil.h"



void SimpleLinReg(const Eigen::VectorXd& predict, 
		  const Eigen::VectorXd& response, 
		  double& weight)

{

  double scale = 1.0/((float)predict.size());
  weight = (predict.transpose() * response) - scale * (predict.sum() * response.sum());
  weight /= predict.transpose() * predict - scale * (predict.sum() * predict.sum());
   
}
