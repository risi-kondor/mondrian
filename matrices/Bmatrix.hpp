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


#ifndef _Bmatrix
#define _Bmatrix

#include "Matrix.hpp"
#include "Detachable.hpp"
#include "MatrixTranspose.hpp"
#include "Bvector.hpp"

extern default_random_engine randomNumberGenerator;


namespace Mondrian{

class Bmatrix: public Cmatrix{
public:

  using Cmatrix::Cmatrix;

public: // copying 

  Bmatrix(const Bmatrix& x): Bmatrix(x.nrows,x.ncols){
    COPY_WARNING("Bmatrix");
    std::copy(x.array,x.array+nrows*ncols,array);
    //symmetricp=x.symmetricp;
  }

  Bmatrix(Bmatrix&& x): Cmatrix(x.nrows,x.ncols){
    MOVE_WARNING("Bmatrix");
    array=x.array; x.array=nullptr; x.nrows=0; x.ncols=0; 
    //symmetricp=x.symmetricp;
  }

  Bmatrix& operator=(const Bmatrix& x){
    ASSIGN_WARNING("Bmatrix");
    nrows=x.nrows; ncols=x.ncols; 
    //symmetricp=x.symmetricp;  
    delete[] array; array=new SCALAR[nrows*ncols];
    std::copy(x.array,x.array+nrows*ncols,array);
    return *this;
  }

  Bmatrix& operator=(Bmatrix&& x){
    MOVEASSIGN_WARNING("Bmatrix");
    nrows=x.nrows; ncols=x.ncols; 
    //symmetricp=x.symmetricp;
    delete[] array; array=x.array; x.array=nullptr; x.nrows=0; x.ncols=0;
    return *this;
  }

  Bmatrix copy() const{
    Bmatrix M(nrows,ncols /*,symmetricp*/); 
    std::copy(array,array+nrows*ncols,M.array);
    return M;
  }

  Bmatrix shallow() const{
    Bmatrix M(nrows,ncols /*,symmetricp*/); M.array=array;
    return M;
  }

  Bmatrix(const initializer_list<Bvector> list): Bmatrix(list.begin()->n,list.size()){
    array=new SCALAR[nrows*ncols]; 
    int j=0; for(const Bvector& v:list){
      for(int i=0; i<nrows; i++) (*this)(i,j)=v(i); j++;}
  }
  


public: // named constructors 

  static Bmatrix Zero(const int _nrows, const int _ncols){
    Bmatrix M(_nrows,_ncols); for(int i=0; i<_nrows*_ncols; i++) M.array[i]=0; /*M.symmetricp=true;*/ return M;}

  static Bmatrix Filled(const int _nrows, const int _ncols, const SCALAR t){
    Bmatrix M(_nrows,_ncols); for(int i=0; i<_nrows*_ncols; i++) M.array[i]=t; /*M.symmetricp=true;*/ return M;}

  static Bmatrix Identity(const int nrows){
    Bmatrix M=Zero(nrows,nrows); for(int i=0; i<nrows; i++) M.array[i*nrows+i]=1; /*M.symmetricp=true;*/ return M;}

  static Bmatrix Uniform(const int nrows, const int ncols){
    Bmatrix M(nrows,ncols); 
    uniform_real_distribution<SCALAR> distr;
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<nrows; j++)
	M.array[i+j*nrows]=distr(randomNumberGenerator); 
    return M;
  }

  static Bmatrix UniformSymmetric(const int nrows, const int ncols){
    Bmatrix M(nrows,ncols,true);
    //M.setSymmetric(true);
    uniform_real_distribution<SCALAR> distr;
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<=i; j++){
	SCALAR t=distr(randomNumberGenerator); 
	M.array[i+j*nrows]=t; 
	M.array[j+i*nrows]=t;
      } 
    return M;
  }

  static Bmatrix Gaussian(const int nrows, const int ncols){
    Bmatrix M(nrows,ncols); 
    normal_distribution<SCALAR> distr;
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<nrows; j++)
	M.array[i+j*nrows]=distr(randomNumberGenerator); 
    return M;
  }

  static Bmatrix GaussianSymmetric(const int nrows){
    Bmatrix M(nrows,nrows,true);
    //M.setSymmetric(true);
    normal_distribution<SCALAR> distr;
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<=i; j++){
	SCALAR t=distr(randomNumberGenerator); 
	M.array[j*nrows+i]=t; 
	M.array[i*nrows+j]=t;
      }
    return M;
  }

  static Bmatrix Bernoulli(const int nrows, const int ncols, const double p=0.5){
    Bmatrix M(nrows,ncols); 
    bernoulli_distribution distr(p);
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<nrows; j++)
	M.array[i+j*nrows]=distr(randomNumberGenerator); 
    return M;
  }

  static Bmatrix BernoulliSymmetric(const int nrows, const double p=0.5){
    Bmatrix M(nrows,nrows,true);
    //M.setSymmetric(true);
    bernoulli_distribution distr(p);
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<=i; j++){
	SCALAR t=distr(randomNumberGenerator); 
	M.array[j*nrows+i]=t; 
	M.array[i*nrows+j]=t;
      }
    return M;
  }





public: // conversions

  //MatrixTranspose<Bmatrix> operator~(){return MatrixTranspose<Bmatrix>(*this);}
  //MatrixTranspose<const Bmatrix> operator~() const {return MatrixTranspose<const Bmatrix>(*this);}

  //Bmatrix(const MatrixTranspose<Bmatrix>& M):Bmatrix(M.nrows,M.ncols){
  //  TRANSPOSE_WARNING("Bmatrix");
  //  for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) array[i+nrows*j]=M.obj.array[j+ncols*i];}

  //View<Bmatrix> detach(){
  //  View<Bmatrix> M(nrows,ncols,symmetricp);
  //  M.array=array;
  //  return M;
  //}


public: // in-place arithmetic

  Bmatrix& operator+=(const SCALAR& x){
    for(int i=0; i<nrows*ncols; i++) array[i]+=x; return *this;}
  Bmatrix& operator-=(const SCALAR& x){
    for(int i=0; i<nrows*ncols; i++) array[i]-=x; return *this;}
  Bmatrix& operator*=(const SCALAR& x){
    blas_dscal(nrows*ncols,x,array,incr); return *this;}
  Bmatrix& operator/=(const SCALAR& x){
    blas_dscal(nrows*ncols,1.0/x,array,incr); return *this;}

  Bmatrix& operator+=(const Bmatrix& x){
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    blas_daxpy(nrows*ncols,1.0,x.array,x.incr,array,incr); 
    //symmetricp=symmetricp&&x.symmetricp;
    return *this;
  }
  Bmatrix& operator-=(const Bmatrix& x){
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    blas_daxpy(nrows*ncols,-1.0,x.array,x.incr,array,incr); 
    //symmetricp=symmetricp&&x.symmetricp;
    return *this;
  }


public: // vector valued arithmetic

  Bvector operator*(const Bvector& v) const{
    assert(v.n==ncols);
    Bvector r(nrows);
    blas_dgemv(BlasColMajor,BlasNoTrans,nrows,ncols,1.0,array,nrows,r.array,r.incr,0,NULL,1.0);
    return r;
  }
  
  Bvector dot(const Bvector& v) const{
    assert(v.n==nrows);
    Bvector r=Bvector::Zero(ncols);
    for(int i=0; i<ncols; i++)
      for(int j=0; j<nrows; j++)
	r.array[i]+=array[i*nrows+j]*v.array[j];
    return r;
  }

    
public: // matrix valued arithmetic

  Bmatrix operator+(const SCALAR x) const{
    Bmatrix M(nrows,ncols,symmetricp); 
    for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]+x; 
    return M;}
  Bmatrix operator-(const SCALAR x) const{
    Bmatrix M(nrows,ncols,symmetricp); 
    for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]-x; 
    return M;}
  Bmatrix operator*(const SCALAR x) const{
    Bmatrix M(nrows,ncols); blas_dscal(nrows*ncols,x,M.array,incr); return M;}
  Bmatrix operator/(const SCALAR x) const{
    Bmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]/x; return M;}

  Bmatrix operator+(const Bmatrix B) const{
    assert(B.nrows==nrows); assert(B.ncols==ncols);
    Bmatrix M(nrows,ncols,symmetricp&&B.symmetricp); 
    blas_daxpy(nrows*ncols,1.0,B.array,incr,array,incr);
    return M;
  }

  Bmatrix operator-(const Bmatrix B) const{
    assert(B.nrows==nrows); assert(B.ncols==ncols);
    Bmatrix M(nrows,ncols,symmetricp&&B.symmetricp); 
    for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]+B.array[i]; 
    return M;
  }

  Bmatrix operator*(const Bmatrix& X) const{
    Bmatrix M(nrows,X.ncols);
    assert(ncols==X.nrows);
    blas_dgemm(BlasColMajor,BlasNoTrans,BlasNoTrans,M.nrows,M.ncols,ncols,1.0,nrows*ncols,nrows,X.nrows*X.ncols,X.nrows,0,NULL,M.nrows);
    return M;
  }
  
  Bmatrix dot(const Bmatrix& x) const{
    Bmatrix M(ncols,x.ncols,this==&x);
    if(this==&x) for(int i=0; i<ncols; i++) for(int j=0; j<=i; j++){
      SCALAR t=0; for(int k=0; k<nrows; k++) t+=array[i*nrows+k]*x.array[j*x.nrows+k]; M(i,j)=t; M(j,i)=t;}
    else for(int i=0; i<ncols; i++) for(int j=0; j<x.ncols; j++){
      SCALAR t=0; for(int k=0; k<nrows; k++) t+=array[i*nrows+k]*x.array[j*x.nrows+k]; M(i,j)=t;}
    return M;
  }


public: // in-place operations 

  Bmatrix& multiplyRowsBy(const Bvector& v){
    assert(v.n==nrows);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	array[j*nrows+i]=v.array[i]*array[j*nrows+i];
    return *this;}

  Bmatrix& multiplyColsBy(const Bvector& v){
    assert(v.n==ncols);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	array[j*nrows+i]=v.array[j]*array[j*nrows+i];
    return *this;}

  Bmatrix& divideRowsBy(const Bvector& v){
    assert(v.n==nrows);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	array[j*nrows+i]=array[j*nrows+i]/v.array[i];
    return *this;}

  Bmatrix& divideColsBy(const Bvector& v){
    assert(v.n==ncols);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	array[j*nrows+i]=array[j*nrows+i]/v.array[j];
    return *this;}

  
public: // scalar-valued methods

  int nnz() const {int t=0; for(int i=0; i<nrows*ncols; i++) if(array[i]!=0) t++; return t;}
  int nFilled() const {return nrows*ncols;}
  
  SCALAR norm2() const{
    SCALAR t=0; for(int i=0; i<nrows*ncols; i++) t+=array[i]*array[i]; return t;}

  SCALAR diff2(const Bmatrix& X) const{
    assert(X.nrows==nrows); assert(X.ncols==ncols); SCALAR t=0;
    for(int i=0; i<nrows*ncols; i++) t+=(array[i]-X.array[i])*(array[i]-X.array[i]);
    return t;}

    
public:
  int incr=1;

};




inline Bmatrix operator*(const SCALAR c, const Bmatrix& M){return M*c;}


  //inline MatrixTranspose<Bmatrix> transp(Bmatrix& M){return MatrixTranspose<Bmatrix>(M);}
  //inline MatrixTranspose<const Bmatrix> transp(const Bmatrix& M){return MatrixTranspose<const Bmatrix>(M);}

} // namespace Mondrian 

template<>
SCALAR& Mondrian::MatrixTranspose<Mondrian::Bmatrix>::operator()(const int i, const int j){return obj(j,i);}


#endif
