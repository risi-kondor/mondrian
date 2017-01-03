#ifndef _EigenInterface
#define _EigenInterface

#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <unsupported/Eigen/MatrixFunctions>
	
#include "Cmatrix.hpp"

// The purpose of these adaptors is to avoid having to include Eigen/Dense or Eigen/Core in any of the 
// header files of the native vector/matrix classes, which would slow down compilation.  

using namespace Mondrian;

typedef Eigen::SparseMatrix<SCALAR> EigenSparseMatrix;

class EigenVectorXdAdaptor: public Eigen::VectorXd{
public:
  EigenVectorXdAdaptor(const Eigen::VectorXd& M): Eigen::VectorXd(M){}
  operator Cvector(){
    Cvector v(size());
    for(int i=0; i<v.n; i++) v.array[i]=(*this)(i);
    return v;
  }
};

class EigenMatrixXdAdaptor: public Eigen::MatrixXd{
public:
  EigenMatrixXdAdaptor(const Eigen::MatrixXd& M): Eigen::MatrixXd(M){}
};


#endif
