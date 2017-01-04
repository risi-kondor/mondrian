/* ---------------------------------------------------------------------------
 
  Mondrian: open source multithreaded matrix library 

  Copyright (C) 2017 Imre Risi Kondor
  Copyright (C) 2015 Imre Risi Kondor, Nedelina Teneva, Pramod K Mudrakarta

  Parts of the following code are derived from the pMMF library 
  (https://github.com/risi-kondor/pMMF) which is licensed under the 
  GNU Public License, version 3. This code therefore is also licensed 
  under the terms of the GNU Public License, version 3. 
 
  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 3
  of the License, or (at your option) any later version.
 
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program; if not, see <http://www.gnu.org/licenses/>.

 --------------------------------------------------------------------------- */


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
