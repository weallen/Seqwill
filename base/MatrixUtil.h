#ifndef MATRIX_UTIL_H_
#define MATRIX_UTIL_H_

#include <Eigen/Dense>

/*void PseudoInverse(Eigen::MatrixXd& pinvat)
{
  ei_assert(m_isInitialized && "SVD is not initialized.");
  double  pinvtoler=1.e-6; // choose your tolerance widely!
  Eigen::SingularValuesType m_sigma_inv=m_sigma;
  for ( long i=0; i<m_workMatrix.cols(); ++i) {
    if ( m_sigma(i) > pinvtoler )
      m_sigma_inv(i)=1.0/m_sigma(i);
    else m_sigma_inv(i)=0;
  }
  pinvmat= (m_matV*m_sigma_inv.asDiagonal()*m_matU.transpose());
  }*/

void SimpleLinReg(const Eigen::VectorXd& predict, const Eigen::VectorXd& response,
		  double& weight);


#endif
