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


#ifndef _Cmatrix
#define _Cmatrix

#include "Matrix.hpp"
#include "MatrixTranspose.hpp"
#include "Detachable.hpp"
#include "Cvector.hpp"
#include "GivensRotation.hpp"

//#include "KpointOp.hpp"

extern default_random_engine randomNumberGenerator;


namespace Mondrian{

  class AsCmatrix;  


class Cmatrix: public Matrix, public Detachable{
public:

  SCALAR* array;


public: // constructors

  Cmatrix(): Cmatrix(0,0){}

  Cmatrix(const int _nrows, const int _ncols): 
    Matrix(_nrows,_ncols){
    array=new SCALAR[nrows*ncols];
  }

  Cmatrix(const int _nrows, const int _ncols, const _Uninitialized dummy): 
    Matrix(_nrows,_ncols){}

  Cmatrix(const int _nrows, const int _ncols, const _Zero dummy):
    Cmatrix(_nrows,_ncols){
    std::fill(array,array+_nrows*_ncols,0);
  }

  Cmatrix(const initializer_list<Cvector> list): Cmatrix(list.begin()->n,list.size()){
    array=new SCALAR[nrows*ncols]; 
    int j=0; for(const Cvector& v:list) {assert(v.n==nrows);
      for(int i=0; i<nrows; i++) array[j*nrows+i]=v(i); j++;}
  }
  
  Cmatrix(const int _nrows, const int _ncols, const initializer_list<iivtriple> list): 
    Cmatrix(_nrows,_ncols){
    for(int i=0; i<nrows*ncols; i++) array[i]=0;
    for(const iivtriple& p:list){
      assert(p.first<nrows); assert(p.second<ncols); 
      array[p.first+p.second*nrows]=p.third;}
  }

  ~Cmatrix(){delete[] array;}


public: // copying 


  Cmatrix(const Cmatrix& x): Cmatrix(x.nrows,x.ncols){
    COPY_WARNING("Cmatrix");
    std::copy(x.array,x.array+nrows*ncols,array);
  }

  Cmatrix(const Cmatrix& x, const _NoWarn dummy): 
    Cmatrix(x.nrows,x.ncols){
    std::copy(x.array,x.array+nrows*ncols,array);
  }

  Cmatrix(Cmatrix&& x): Cmatrix(x.nrows,x.ncols,_Uninitialized()){
    MOVE_WARNING("Cmatrix");
    array=x.array; x.array=nullptr; x.nrows=0; x.ncols=0; 
  }

  Cmatrix(Cmatrix&& x, const _NoWarn dummy): 
    Cmatrix(x.nrows,x.ncols,_Uninitialized()){
    array=x.array; x.array=nullptr; x.nrows=0; x.ncols=0; 
  }

  Cmatrix& operator=(const Cmatrix& x){
    ASSIGN_WARNING("Cmatrix");
    nrows=x.nrows; ncols=x.ncols; 
    delete[] array; array=new SCALAR[nrows*ncols];
    std::copy(x.array,x.array+nrows*ncols,array);
    return *this;
  }

  Cmatrix& operator=(Cmatrix&& x){
    MOVEASSIGN_WARNING("Cmatrix");
    nrows=x.nrows; ncols=x.ncols; x.nrows=0; x.ncols=0;
    delete[] array; array=x.array; x.array=nullptr;
    return *this;
  }

  Cmatrix copy() const{
    Cmatrix M(nrows,ncols); 
    std::copy(array,array+nrows*ncols,M.array);
    return M;
  }

  void detach(){array=nullptr;}

  void assign(const Cmatrix& x){
    assert(x.nrows==nrows); assert(x.ncols==ncols);  
    std::copy(x.array,x.array+nrows*ncols,array);}

  Cmatrix shallow() const{
    Cmatrix M(nrows,ncols,_Uninitialized()); M.array=array;
    return M;
  }


public: // conversions


  Cmatrix(const int _nrows, const int _ncols, const SCALAR* _array): 
    Cmatrix(_nrows,_ncols){
    COPY_WARNING("SCALAR[]");
    std::copy(_array,_array+nrows*ncols,array);
  }

  Cmatrix(const int _nrows, const int _ncols, SCALAR*&& _array): 
    Cmatrix(_nrows,_ncols){
    MOVE_WARNING("SCALAR[]");
    array=_array;
  }

 
public: // transpose


  MatrixTranspose<Cmatrix> operator~(){return MatrixTranspose<Cmatrix>(*this);}
  //MatrixTranspose<const Cmatrix> operator~() const {return MatrixTranspose<const Cmatrix>(*this);}

  Cmatrix(const MatrixTranspose<Cmatrix>& M):Cmatrix(M.nrows,M.ncols){
    TRANSPOSE_WARNING("Cmatrix");
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
	array[i+nrows*j]=M.obj.array[j+ncols*i];}


public: // named constructors 


  static Cmatrix Zero(const int _nrows, const int _ncols){
    Cmatrix M(_nrows,_ncols); std::fill(M.array,M.array+_nrows*_ncols,0); return M;}

  static Cmatrix Filled(const int _nrows, const int _ncols, const SCALAR t){
    Cmatrix M(_nrows,_ncols); std::fill(M.array,M.array+_nrows*_ncols,0); return M;}

  static Cmatrix Identity(const int nrows){
    Cmatrix M=Zero(nrows,nrows); for(int i=0; i<nrows; i++) M.array[i*nrows+i]=1; return M;}

  template<class VECTOR>
  static Cmatrix Diag(const VECTOR& v){
    Cmatrix M=Zero(v.n,v.n); 
    for(int i=0; i<v.n; i++) M.array[i*v.n+i]=v(i); 
    return M;
  }

  static Cmatrix Uniform(const int nrows, const int ncols){
    Cmatrix M(nrows,ncols); 
    uniform_real_distribution<SCALAR> distr;
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<nrows; j++)
	M.array[i+j*nrows]=distr(randomNumberGenerator); 
    return M;
  }

  static Cmatrix Gaussian(const int nrows, const int ncols){
    Cmatrix M(nrows,ncols); 
    normal_distribution<SCALAR> distr;
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<nrows; j++)
	M.array[i+j*nrows]=distr(randomNumberGenerator); 
    return M;
  }

  static Cmatrix Bernoulli(const int nrows, const int ncols, const double p=0.5){
    Cmatrix M(nrows,ncols); 
    bernoulli_distribution distr(p);
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++)
	M.array[i+j*nrows]=distr(randomNumberGenerator); 
    return M;
  }


public: // comparisons


  bool operator==(const Cmatrix& X) const{ 
    if(X.nrows!=nrows) return false; if(X.ncols!=ncols) return false;
    for(int i=0; i<nrows*ncols; i++) if(array[i]!=X.array[i]) return false;
    return true;}


public: // element access


  SCALAR operator()(const int i, const int j) const{ 
    assert(i<nrows); assert(j<ncols); return array[j*nrows+i];} 
  SCALAR read(const int i, const int j) const{ 
    assert(i<nrows); assert(j<ncols); return array[j*nrows+i];} 

  SCALAR& operator()(const int i, const int j){
    assert(i<nrows); assert(j<ncols); return array[j*nrows+i];} 
  SCALAR& element(const int i, const int j){
    assert(i<nrows); assert(j<ncols); return array[j*nrows+i];} 
  void set(const int i, const int j, const SCALAR v){ 
    assert(i<nrows); assert(j<ncols); array[j*nrows+i]=v;} 

  bool isFilled(const int i, const int j) const {return true;}

  int nFilled() const {return nrows*ncols;}

  template<class VECTOR>
  VECTOR row(const int i) const {assert(i<nrows);
    VECTOR x(ncols); for(int j=0; j<ncols; j++) x(j)=array[j*nrows+i]; return x;}

  template<class VECTOR>
  VECTOR column(const int j) const {assert(j<ncols);
    VECTOR x(nrows); for(int i=0; i<nrows; i++) x(i)=array[j*nrows+i]; return x;}

  template<class VECTOR>
  VECTOR diag() const {int t=min(nrows,ncols);
    VECTOR x(t); for(int i=0; i<t; i++) x(i)=array[i*nrows+i]; return x;}


public: // function mappings ----------------------------------------------------------------------------------------


  void for_each(std::function<void(const INDEX, const INDEX, const SCALAR)> lambda) const{
    for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) lambda(i,j,array[j*nrows+i]);}
  void for_each_in_row(const int i, std::function<void(const INDEX, const SCALAR)> lambda) const{
    for(int j=0; j<ncols; j++) lambda(i,array[j*nrows+i]);}
  void for_each_in_column(const int j, std::function<void(const INDEX, const SCALAR)> lambda) const{
    for(int i=0; i<nrows; i++) lambda(i,array[j*nrows+i]);}

  void for_each(std::function<void(const INDEX, const INDEX, SCALAR&)> lambda){
    for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) lambda(i,j,array[j*nrows+i]);}
  void for_each_in_row(const int i, std::function<void(const INDEX, SCALAR&)> lambda){
    for(int j=0; j<ncols; j++) lambda(i,array[j*nrows+i]);}
  void for_each_in_column(const int j, std::function<void(const INDEX, SCALAR&)> lambda){
    for(int i=0; i<nrows; i++) lambda(i,array[j*nrows+i]);}

  Cmatrix elementwise(std::function<SCALAR(const SCALAR)> lambda){
    Cmatrix r(nrows,ncols); for(int i=0; i<nrows*ncols; i++) r.array[i]=lambda(array[i]); return r;}


public: // views ----------------------------------------------------------------------------------------------------


  Detached<Cvector> view_of_column(const int j){
    Detached<Cvector> v; v.n=nrows; v.array=array+j*nrows;
    return v;
  }

  Detached<const Cvector> view_of_column(const int j) const{
    Detached<const Cvector> v; v.n=nrows; v.array=array+j*nrows;
    return v;
  }

  Detached<Cmatrix> view_of_columns(const int j, const int k){
    Detached<Cmatrix> M; // TODO
    M.array=array+j*nrows;
    return M;
  }


public: // remappings -----------------------------------------------------------------------------------------------


  Cmatrix remapRows(const IndexMap& map) const{
    Cmatrix M(map.nsource,ncols); 
    for(int i=0; i<map.nsource; i++) for(int j=0; j<ncols; j++) M.array[j*nrows+i]=array[j*nrows+map(i)];
    return M;}

  Cmatrix remapRows(const Inverse<IndexMap>& imap) const{
    assert(imap.obj.nsource==nrows); Cmatrix M(nrows,ncols); 
    for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) M.array[j*nrows+imap.obj(i)]=array[j*nrows+i];
    return M;}

  Cmatrix remapCols(const IndexMap& map) const{
    Cmatrix M(nrows,map.nsource); 
    for(int j=0; j<map.nsource; j++) for(int i=0; i<nrows; i++) M.array[j*nrows+i]=array[map(j)*nrows+i];
    return M;}

  Cmatrix remapCols(const Inverse<IndexMap>& imap) const{
    assert(imap.obj.nsource==ncols); Cmatrix M(nrows,ncols); 
    for(int j=0; j<ncols; j++) for(int i=0; i<nrows; i++) M.array[imap.obj(j)*nrows+i]=array[j*nrows+i];
    return M;}

  Cmatrix remap(const IndexMap& rmap, const IndexMap& cmap) const{
    Cmatrix M(rmap.nsource,cmap.nsource); 
    for(int i=0; i<rmap.nsource; i++) 
      for(int j=0; j<cmap.nsource; j++) 
	M.array[j*nrows+i]=array[rmap(j)*nrows+cmap(i)];
    return M;}

  Cmatrix remap(const Inverse<IndexMap>& irmap, const Inverse<IndexMap>& icmap) const{
    assert(irmap.obj.nsource==nrows); assert(icmap.obj.nsource=ncols); Cmatrix M(nrows,ncols); 
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
	M.array[icmap.obj(j)*nrows+irmap.obj(i)]=array[j*nrows+i];
    return M;}


public: // in-place arithmetic ------------------------------------------------------------------------------------


  Cmatrix& operator+=(const SCALAR& x){
    for(int i=0; i<nrows*ncols; i++) array[i]+=x; return *this;}
  Cmatrix& operator-=(const SCALAR& x){
    for(int i=0; i<nrows*ncols; i++) array[i]-=x; return *this;}
  Cmatrix& operator*=(const SCALAR& x){
    for(int i=0; i<nrows*ncols; i++) array[i]*=x; return *this;}
  Cmatrix& operator/=(const SCALAR& x){
    for(int i=0; i<nrows*ncols; i++) array[i]/=x; return *this;}

  Cmatrix& operator+=(const Cmatrix& x){
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    for(int i=0; i<nrows*ncols; i++) array[i]+=x.array[i]; 
    return *this;}
  Cmatrix& operator-=(const Cmatrix& x){
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    for(int i=0; i<nrows*ncols; i++) array[i]-=x.array[i]; 
    return *this;}
  Cmatrix& operator*=(const Cmatrix& x){
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    for(int i=0; i<nrows*ncols; i++) array[i]*=x.array[i]; 
    return *this;}
  Cmatrix& operator/=(const Cmatrix& x){
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    for(int i=0; i<nrows*ncols; i++) array[i]/=x.array[i]; 
    return *this;}

    
public: // in-place operations --------------------------------------------------------------------------------------


  template<class VECTOR>
  Cmatrix& multiplyRowsBy(const VECTOR& v){
    assert(v.n==nrows);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	array[j*nrows+i]=v(i)*array[j*nrows+i];
    return *this;}

  template<class VECTOR>
  Cmatrix& multiplyColumnsBy(const VECTOR& v){
    assert(v.n==ncols);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	array[j*nrows+i]=v(j)*array[j*nrows+i];
    return *this;}

  template<class VECTOR>
  Cmatrix& divideRowsBy(const VECTOR& v){
    assert(v.n==nrows);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	array[j*nrows+i]=array[j*nrows+i]/v(i);
    return *this;}

  template<class VECTOR>  
  Cmatrix& divideColumnsBy(const VECTOR& v){
    assert(v.n==ncols);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	array[j*nrows+i]=array[j*nrows+i]/v(j);
    return *this;}

  Cmatrix& normalizeRows(){
    for(int i=0; i<nrows; i++){
      SCALAR t=0;
      for(int j=0; j<ncols; j++) t+=array[j*nrows+i];
      for(int j=0; j<ncols; j++) array[j*nrows+i]/=t;}
    return *this;
  }

  Cmatrix& normalizeColumns(){ 
    for(int j=0; j<ncols; j++){
      SCALAR t=0;
      for(int i=0; i<nrows; i++) t+=array[j*nrows+i];
      for(int i=0; i<nrows; i++) array[j*nrows+i]/=t;}
    return *this;
  }

  void symmetrize(){
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++){
	SCALAR t=(array[j*nrows+i]+array[i*nrows+j])/2.0;
	array[j*nrows+i]=t;
	array[i*nrows+j]=t;
      }
  }


public: // construct operations ------------------------------------------------------------------------------------


  Cmatrix(const Mondrian::xcgram<Cmatrix>& e): 
    Cmatrix(e.A.ncols,e.A.ncols){
    const Cmatrix& M=e.A;
    const int n=M.ncols;
    const int m=M.nrows;
    for(int i=0; i<n; i++)
      for(int j=0; j<=i; j++){
	SCALAR t=0;
	for(int k=0; k<m; k++)
	  t+=M.array[i*m+k]*M.array[j*m+k];
	array[i*n+j]=t;
	array[j*n+i]=t;
      }
  }

  template<class VECTOR>
  Cmatrix(const Mondrian::xcgram<MatrixX<VECTOR> >& e): 
    Cmatrix(e.A.ncols,e.A.ncols){
    const MatrixX<VECTOR>& M=e.A;
    const int n=M.ncols;
    for(int i=0; i<n; i++)
      for(int j=0; j<=i; j++){
	SCALAR t=M.cols[i]->dot(*M.cols[j]);
	array[i*n+j]=t;
	array[j*n+i]=t;
      }
  }


public: // vector valued arithmetic ---------------------------------------------------------------------------------


  Cvector operator*(const Cvector& v) const{
    assert(v.n==ncols);
    Cvector r=Cvector::Zero(nrows);
    for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
	r.array[i]+=array[j*nrows+i]*v.array[j];
    return r;
  }
  
  Cvector dot(const Cvector& v) const{
    assert(v.n==nrows);
    Cvector r=Cvector::Zero(ncols);
    for(int i=0; i<ncols; i++)
      for(int j=0; j<nrows; j++)
	r.array[i]+=array[i*nrows+j]*v.array[j];
    return r;
  }

    
public: // matrix valued arithmetic ---------------------------------------------------------------------------------


  Cmatrix operator+(const SCALAR x) const{
    Cmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]+x; return M;}
  Cmatrix operator-(const SCALAR x) const{
    Cmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]-x; return M;}
  Cmatrix operator*(const SCALAR x) const{
    Cmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]*x; return M;}
  Cmatrix operator/(const SCALAR x) const{
    Cmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]/x; return M;}

  Cmatrix operator+(const Cmatrix B) const{ assert(B.nrows==nrows); assert(B.ncols==ncols);
    Cmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]+B.array[i]; return M;}
  Cmatrix operator-(const Cmatrix B) const{ assert(B.nrows==nrows); assert(B.ncols==ncols);
    Cmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]+B.array[i]; return M;}

  Cmatrix operator*(const Cmatrix& X) const{
    Cmatrix M=Zero(nrows,X.ncols);
    assert(ncols==X.nrows);
    for(int i=0; i<nrows; i++)
      for(int j=0; j<X.ncols; j++)
	for(int k=0; k<ncols; k++)
	  M.array[j*nrows+i]+=array[k*nrows+i]*X.array[j*X.nrows+k];
    return M;
  }
  
  Cmatrix dot(const Cmatrix& x) const{
    Cmatrix M(ncols,x.ncols);
    if(this==&x) for(int i=0; i<ncols; i++) for(int j=0; j<=i; j++){
      SCALAR t=0; for(int k=0; k<nrows; k++) t+=array[i*nrows+k]*x.array[j*x.nrows+k]; M(i,j)=t; M(j,i)=t;}
    else for(int i=0; i<ncols; i++) for(int j=0; j<x.ncols; j++){
      SCALAR t=0; for(int k=0; k<nrows; k++) t+=array[i*nrows+k]*x.array[j*x.nrows+k]; M(i,j)=t;}
    return M;
  }

  Cmatrix outer(const Cmatrix& x) const{
    assert(ncols==x.ncols);
    Cmatrix M(nrows,x.nrows);
    if(this==&x) for(int i=0; i<nrows; i++) for(int j=0; j<=i; j++){
      SCALAR t=0; for(int k=0; k<ncols; k++) t+=array[k*nrows+i]*x.array[k*x.nrows+j]; M(i,j)=t; M(j,i)=t;}
    else for(int i=0; i<nrows; i++) for(int j=0; j<x.nrows; j++){
      SCALAR t=0; for(int k=0; k<ncols; k++) t+=array[k*nrows+i]*x.array[k*x.nrows+j]; M(i,j)=t;}
    return M;
  }
  

public: // scalar-valued methods ------------------------------------------------------------------------------------


  int nnz() const {
    int t=0; for(int i=0; i<nrows*ncols; i++) if(array[i]!=0) t++; return t;}
  
  bool isSymmetric() const{
    for(int j=0; j<ncols; j++)
      for(int i=0; i<j; j++)
	if(array[j*nrows+i]!=array[i*nrows+j]) return false;
    return true;
  }

  SCALAR inp_of_columns(const int j1, const int j2) const{
    SCALAR t=0; for(int i=0; i<nrows; i++) t+=array[j1*nrows+i]*array[j2*nrows+i]; return t;}

  SCALAR inp_of_rows(const int i1, const int i2) const{
    SCALAR t=0; for(int j=0; j<ncols; j++) t+=array[j*nrows+i1]*array[j*nrows+i2]; return t;}

  SCALAR norm2() const{
    SCALAR t=0; for(int i=0; i<nrows*ncols; i++) t+=array[i]*array[i]; return t;}

  SCALAR diff2(const Cmatrix& X) const{
    assert(X.nrows==nrows); assert(X.ncols==ncols); SCALAR t=0;
    for(int i=0; i<nrows*ncols; i++) t+=(array[i]-X.array[i])*(array[i]-X.array[i]);
    return t;
  }

  template<class MATRIX>
  SCALAR diff2(const MATRIX& X) const{
    assert(X.nrows==nrows); assert(X.ncols==ncols); SCALAR t=0;
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++)
	t+=(array[j*nrows+i]-X(i,j))*(array[j*nrows+i]-X(i,j));
    return t;
  }
    

public: // Givens rotations -----------------------------------------------------------------------------------------
  

  Cmatrix& applyFromLeft(const GivensRotation& Q){
    assert(Q.i1<nrows);
    assert(Q.i2<nrows);
    for(int j=0; j<ncols; j++){
      SCALAR v1=array[j*nrows+Q.i1];
      SCALAR v2=array[j*nrows+Q.i2];
      array[j*nrows+Q.i1]=Q.cos*v1-Q.sin*v2;
      array[j*nrows+Q.i2]=Q.cos*v2+Q.sin*v1;
    }
    return *this;
  }

  Cmatrix& applyFromLeftT(const GivensRotation& Q){
    assert(Q.i1<nrows);
    assert(Q.i2<nrows);
    for(int j=0; j<ncols; j++){
      SCALAR v1=array[j*nrows+Q.i1];
      SCALAR v2=array[j*nrows+Q.i2];
      array[j*nrows+Q.i1]=Q.cos*v1+Q.sin*v2;
      array[j*nrows+Q.i2]=Q.cos*v2-Q.sin*v1;
    }
    return *this;
  }

  Cmatrix& applyFromRightT(const GivensRotation& Q){
    assert(Q.i1<ncols);
    assert(Q.i2<ncols);
    for(int i=0; i<nrows; i++){
      SCALAR v1=array[Q.i1*nrows+i];
      SCALAR v2=array[Q.i2*nrows+i];
      array[Q.i1*nrows+i]=Q.cos*v1-Q.sin*v2;
      array[Q.i2*nrows+i]=Q.cos*v2+Q.sin*v1;
    }
    return *this;
  }

  Cmatrix& applyFromRight(const GivensRotation& Q){
    assert(Q.i1<ncols);
    assert(Q.i2<ncols);
    for(int i=0; i<nrows; i++){
      SCALAR v1=array[Q.i1*nrows+i];
      SCALAR v2=array[Q.i2*nrows+i];
      array[Q.i1*nrows+i]=Q.cos*v1+Q.sin*v2;
      array[Q.i2*nrows+i]=Q.cos*v2-Q.sin*v1;
    }
    return *this;
  }

  Cmatrix& conjugate(const GivensRotation& Q){
    applyFromLeft(Q);
    applyFromRightT(Q);
    return *this;
  }

  Cmatrix& conjugateT(const GivensRotation& Q){
    applyFromLeftT(Q);
    applyFromRight(Q);
    return *this;
  }


public: // KpointOp -------------------------------------------------------------------------------------------------

  
  template<int k>
  Cmatrix& applyFromLeft(const KpointOp<k>& Q){
    SCALAR temp[k];
    for(int i=0; i<k; i++) assert(Q.map(i)<nrows); 
    for(int j=0; j<ncols; j++){
      for(int i=0; i<k; i++) 
	temp[i]=array[j*nrows+Q.map(i)];
      for(int i=0; i<k; i++){
	SCALAR t=0; 
	for(int l=0; l<k; l++) 
	  t+=temp[l]*Q.q[l*k+i]; 
	array[j*nrows+Q.map(i)]=t;
      }
    }
    return *this;
  }
  
  template<int k>
  Cmatrix& applyFromLeftT(const KpointOp<k>& Q){
    SCALAR temp[k];
    for(int i=0; i<k; i++) assert(Q.map(i)<nrows); 
    for(int j=0; j<ncols; j++){
      for(int i=0; i<k; i++) temp[i]=array[j*nrows+Q.map(i)];
      for(int i=0; i<k; i++){
	SCALAR t=0; 
	for(int l=0; l<k; l++) 
	  t+=temp[l]*Q.q[i*k+l]; 
	array[j*nrows+Q.map(i)]=t;
      }
    }
    return *this;
  }
  
  template<int k>
  Cmatrix& applyFromRightT(const KpointOp<k>& Q){
    SCALAR temp[k];
    for(int i=0; i<k; i++) assert(Q.map(i)<ncols); 
    for(int i=0; i<nrows; i++){
      for(int j=0; j<k; j++) temp[j]=array[Q.map(j)*nrows+i];
      for(int j=0; j<k; j++){
	SCALAR t=0; 
	for(int l=0; l<k; l++) 
	  t+=temp[l]*Q.q[l*k+j]; 
	array[Q.map(j)*nrows+i]=t;
      }
    }
    return *this;
  }

  template<int k>
  Cmatrix& applyFromRight(const KpointOp<k>& Q){
    SCALAR temp[k];
    for(int i=0; i<k; i++) assert(Q.map(i)<ncols); 
    for(int i=0; i<nrows; i++){
      for(int j=0; j<k; j++) temp[j]=array[Q.map(j)*nrows+i];
      for(int j=0; j<k; j++){
	SCALAR t=0; 
	for(int l=0; l<k; l++) 
	  t+=temp[l]*Q.q[j*k+l]; 
	array[Q.map(j)*nrows+i]=t;
      }
    }
    return *this;
  }

  template<int k>
  Cmatrix& conjugate(const KpointOp<k>& Q){
    applyFromLeft(Q);
    applyFromRightT(Q);
    return *this;
  }

  template<int k>
  Cmatrix& conjugateT(const KpointOp<k>& Q){
    applyFromLeftT(Q);
    applyFromRight(Q);
    return *this;
  }


public: // Python interface -----------------------------------------------------------------------------------------


  Cmatrix(double* numpyDblArray, int numpyDim1, int numpyDim2): 
    Cmatrix(numpyDim1, numpyDim2){
    //std::copy(numpyDblArray,numpyDblArray+nrows*ncols,array);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	array[j*nrows+i]=numpyDblArray[i*ncols+j];
  }
  
  void np(double** numpyDblArray, int* numpyDim1, int* numpyDim2){
    *numpyDim1=nrows;
    *numpyDim2=ncols;
    *numpyDblArray=new double[nrows*ncols];
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	(*numpyDblArray)[i*ncols+j]=array[j*nrows+i];
    //std::copy(array,array+nrows*ncols,*numpyDblArray);
  }
  
  const char* __str__() { // why doesn't this work?
    ostringstream ostream; 
    cout<<str()<<endl;
    return ostream.str().c_str();
  }
    

};


  class AsCmatrix: public Cmatrix{
  public:
    AsCmatrix()=delete;
    AsCmatrix& operator=(const AsCmatrix& x)=delete;
    AsCmatrix(const int _nrows, const int _ncols, SCALAR* _array): 
      Cmatrix(_nrows,_ncols,_Uninitialized()){
      array=_array;
    }
    ~AsCmatrix(){array=nullptr;}
  };

  inline AsCmatrix as_Cmatrix(const int _nrows, const int _ncols, SCALAR* _array){
    return AsCmatrix(_nrows,_ncols,_array);}

  inline Cmatrix as_Cmatrix(const int _nrows, const int _ncols, const SCALAR*&& _array){
    return Cmatrix(_nrows,_ncols,std::move(_array));}


  inline Cmatrix operator*(const SCALAR c, const Cmatrix& M){return M*c;}


  //inline MatrixTranspose<Cmatrix> transp(Cmatrix& M){return MatrixTranspose<Cmatrix>(M);}
  //inline MatrixTranspose<const Cmatrix> transp(const Cmatrix& M){return MatrixTranspose<const Cmatrix>(M);}

} // namespace Mondrian 



#endif




  //SCALAR& operator()(const int* ip, const int* jp){
  //  assert(*ip<nrows); assert(*jp<ncols); return array[(*ip)+(*jp)*nrows];} 
  //SCALAR operator()(const int* ip, const int* jp) const{
  //  assert(*ip<nrows); assert(*jp<ncols); return array[(*ip)+(*jp)*nrows];} 

  //Cvector row(const int i) const {assert(i<nrows);
  //  Cvector x(ncols); for(int j=0; j<ncols; j++) x.array[j]=array[j*nrows+i]; return x;}
  //Cvector column(const int j) const {assert(j<ncols);
  //  Cvector x(nrows); for(int i=0; i<nrows; i++) x.array[i]=array[j*nrows+i]; return x;}
  //Cvector diag() const {int t=min(nrows,ncols);
  //  Cvector x(t); for(int i=0; i<t; i++) x.array[i]=array[i*nrows+i]; return x;}

  //bool isFilled(const int i, const int j) const {return true;}
  //int nFilled() const {return nrows*ncols;}
  /*
  void applyFromLeft(const KpointRotation& Q){
    SCALAR temp[Q.k];
    for(int i=0; i<Q.k; i++) assert(Q.ix[i]<nrows); 
    for(int j=0; j<ncols; j++){
      for(int i=0; i<Q.k; i++) 
	temp[i]=array[j*nrows+Q.ix[i]];
      for(int i=0; i<Q.k; i++){
	SCALAR t=0; for(int l=0; l<Q.k; l++) t+=temp[l]*Q.q[l*Q.k+i]; array[j*nrows+Q.ix[i]]=t;}}
  }
  */
  
  /*
  void applyFromLeftT(const KpointRotation& Q){
    SCALAR temp[Q.k];
    for(int i=0; i<Q.k; i++) assert(Q.ix[i]<nrows); 
    for(int j=0; j<ncols; j++){
      for(int i=0; i<Q.k; i++) temp[i]=array[j*nrows+Q.ix[i]];
      for(int i=0; i<Q.k; i++){
	SCALAR t=0; for(int l=0; l<Q.k; l++) t+=temp[l]*Q.q[i*Q.k+l]; array[j*nrows+Q.ix[i]]=t;}}
  }
  */
  
  /*
  void applyFromRightT(const KpointRotation& Q){
    SCALAR temp[Q.k];
    for(int i=0; i<Q.k; i++) assert(Q.ix[i]<ncols); 
    for(int i=0; i<nrows; i++){
      for(int j=0; j<Q.k; j++) temp[j]=array[Q.ix[j]*nrows+i];
      for(int j=0; j<Q.k; j++){
	SCALAR t=0; for(int l=0; l<Q.k; l++) t+=temp[l]*Q.q[l*Q.k+j]; array[Q.ix[j]*nrows+i]=t;}
    }
  }
  */

  /*
  void applyFromRight(const KpointRotation& Q){
    SCALAR temp[Q.k];
    for(int i=0; i<Q.k; i++) assert(Q.ix[i]<ncols); 
    for(int i=0; i<nrows; i++){
      for(int j=0; j<Q.k; j++) temp[j]=array[Q.ix[j]*nrows+i];
      for(int j=0; j<Q.k; j++){
	SCALAR t=0; for(int l=0; l<Q.k; l++) t+=temp[l]*Q.q[j*Q.k+l]; array[Q.ix[j]*nrows+i]=t;}
    }
  }
  */
  //Cmatrix(const int _nrows, const int _ncols, const _Uninitialized dummy): Matrix(_nrows,_ncols){}

//  Cmatrix& symmetrize(){cout<<"Unimplemented"<<endl; return *this;}
//template<> // leads to duplicate symbol 
//SCALAR& Mondrian::MatrixTranspose<Mondrian::Cmatrix>::operator()(const int i, const int j){return obj(j,i);}
  /*
public: // KpointOperators

  Cmatrix& applyFromLeft(const KpointOperator& Q){
    SCALAR temp[Q.k];
    for(int i=0; i<Q.k; i++) assert(Q.map(i)<nrows); 
    for(int j=0; j<ncols; j++){
      for(int i=0; i<Q.k; i++) 
	temp[i]=array[j*nrows+Q.map(i)];
      for(int i=0; i<Q.k; i++){
	SCALAR t=0; for(int l=0; l<Q.k; l++) t+=temp[l]*Q.q[l*Q.k+i]; array[j*nrows+Q.map(i)]=t;}}
    return *this;
  }
  
  Cmatrix& applyFromLeftT(const KpointOperator& Q){
    SCALAR temp[Q.k];
    for(int i=0; i<Q.k; i++) assert(Q.map(i)<nrows); 
    for(int j=0; j<ncols; j++){
      for(int i=0; i<Q.k; i++) temp[i]=array[j*nrows+Q.map(i)];
      for(int i=0; i<Q.k; i++){
	SCALAR t=0; for(int l=0; l<Q.k; l++) t+=temp[l]*Q.q[i*Q.k+l]; array[j*nrows+Q.map(i)]=t;}}
    return *this;
  }
  
  Cmatrix& applyFromRightT(const KpointOperator& Q){
    SCALAR temp[Q.k];
    for(int i=0; i<Q.k; i++) assert(Q.map(i)<ncols); 
    for(int i=0; i<nrows; i++){
      for(int j=0; j<Q.k; j++) temp[j]=array[Q.map(j)*nrows+i];
      for(int j=0; j<Q.k; j++){
	SCALAR t=0; for(int l=0; l<Q.k; l++) t+=temp[l]*Q.q[l*Q.k+j]; array[Q.map(j)*nrows+i]=t;}
    }
    return *this;
  }

  Cmatrix& applyFromRight(const KpointOperator& Q){
    SCALAR temp[Q.k];
    for(int i=0; i<Q.k; i++) assert(Q.map(i)<ncols); 
    for(int i=0; i<nrows; i++){
      for(int j=0; j<Q.k; j++) temp[j]=array[Q.map(j)*nrows+i];
      for(int j=0; j<Q.k; j++){
	SCALAR t=0; for(int l=0; l<Q.k; l++) t+=temp[l]*Q.q[j*Q.k+l]; array[Q.map(j)*nrows+i]=t;}
    }
    return *this;
  }

  Cmatrix& conjugate(const KpointOperator& Q){
    applyFromLeft(Q);
    applyFromRightT(Q);
    return *this;
  }

  Cmatrix& conjugateT(const KpointOperator& Q){
    applyFromLeftT(Q);
    applyFromRight(Q);
    return *this;
  }
  */
//for(int j=0; j<ncols; j++)
//    for(int i=0; i<nrows; i++)
//array[j*nrows+i]=numpyDblArray[i*ncols+j];
  //inline Cmatrix::Cmatrix(const AsCmatrixLA& x):  // maybe no need for this
  //  Cmatrix(x,_NoWarn()) {
  //  COPY_WARNING("CmatrixLA");
  //}

