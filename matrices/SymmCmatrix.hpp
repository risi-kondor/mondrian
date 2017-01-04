/* ---------------------------------------------------------------------------
 
  Mondrian: open source multithreaded matrix library 

  Copyright (C) 2017 Imre Risi Kondor

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


#ifndef _SymmCmatrix
#define _SymmCmatrix

#include "Cmatrix.hpp"

extern default_random_engine randomNumberGenerator;


namespace Mondrian{

  class AsSymmCmatrix;  


class SymmCmatrix: public Cmatrix{
public: // attributes

  bool isSymmetricFormat() const {return true;}


public: // constructors

  SymmCmatrix(): SymmCmatrix(0,0){}

  SymmCmatrix(const int _nrows): Cmatrix(_nrows,_nrows){}

  SymmCmatrix(const int _nrows, const int _ncols): Cmatrix(_nrows,_ncols){
    assert(_nrows==_ncols);}

  SymmCmatrix(const initializer_list<Cvector>& list): Cmatrix(list){
    symmetrize();}
  
  SymmCmatrix(const int _nrows, const int _ncols, const initializer_list<iivtriple>& list): 
    Cmatrix(_nrows,_ncols,list){
    symmetrize();
  }


public: // copying

  SymmCmatrix(const SymmCmatrix& x): Cmatrix(x,_NoWarn()){
    COPY_WARNING("SymmCmatrix");
  }

  SymmCmatrix(SymmCmatrix&& x): Cmatrix(x.nrows,x.ncols,_Uninitialized()){
    MOVE_WARNING("SymmCmatrix");
    array=x.array; x.array=nullptr; x.nrows=0; x.ncols=0; 
  }

  SymmCmatrix& operator=(const SymmCmatrix& x){
    ASSIGN_WARNING("SymmCmatrix");
    nrows=x.nrows; ncols=x.ncols; 
    delete[] array; array=new SCALAR[nrows*ncols];
    std::copy(x.array,x.array+nrows*ncols,array);
    return *this;
  }

  SymmCmatrix& operator=(SymmCmatrix&& x){
    MOVEASSIGN_WARNING("SymmCmatrix");
    nrows=x.nrows; ncols=x.ncols; 
    delete[] array; array=x.array; x.array=nullptr; x.nrows=0; x.ncols=0;
    return *this;
  }

  SymmCmatrix copy() const{
    SymmCmatrix M(nrows,ncols); 
    std::copy(array,array+nrows*ncols,M.array);
    return M;
  }

  void detach(){array=nullptr;}

  void assign(const SymmCmatrix& x){
    assert(x.nrows==nrows); assert(x.ncols==ncols);  
    std::copy(x.array,x.array+nrows*ncols,array);}

  SymmCmatrix shallow() const{
    SymmCmatrix M(nrows,ncols); M.array=array;
    return M;
  }


public: // downcasting


  SymmCmatrix(const Cmatrix& x, const _Downcast dummy): 
    Cmatrix(x,_NoWarn()){
    DOWNCASTCOPY_WARNING("Cmatrix","SymmCmatrix");
  }

  SymmCmatrix(Cmatrix&& x, const _Downcast dummy): 
    Cmatrix(std::move(x),_NoWarn()){
    DOWNCAST_WARNING("Cmatrix","SymmCmatrix");
  }


public: // conversions


  SymmCmatrix(const AsSymmCmatrix& x);

  SymmCmatrix(const Cmatrix& x): Cmatrix(x,_NoWarn()){
    CONVERT_WARNING("Cmatrix","SymmCmatrix");
    symmetrize();
  }

  SymmCmatrix(Cmatrix&& x): Cmatrix(std::move(x),_NoWarn()){
    MOVECONVERT_WARNING("Cmatrix","SymmCmatrix");
    symmetrize();
  }

  template<class MATRIX>
  SymmCmatrix(const MATRIX& x): Cmatrix(x,_NoWarn()){
    CONVERT_WARNING("MATRIX","SymmCmatrix");
    symmetrize();
  }


public: // transpose


  MatrixTranspose<SymmCmatrix> operator~(){
    return MatrixTranspose<SymmCmatrix>(*this);}
  
  //MatrixTranspose<const SymmCmatrix> operator~() const{
  //  return MatrixTranspose<const SymmCmatrix>(*this);}

  SymmCmatrix(const MatrixTranspose<SymmCmatrix>& M):
    SymmCmatrix(M.obj){}

  SymmCmatrix(MatrixTranspose<SymmCmatrix>&& M):
    SymmCmatrix(std::move(M.obj)){}


public: // named constructors 


  static SymmCmatrix Zero(const int _nrows, const int _ncols) {
    SymmCmatrix M(_nrows,_ncols); std::fill(M.array,M.array+_nrows*_ncols,0); return M;}
  static SymmCmatrix Zero(const int _nrows){return Zero(_nrows,_nrows);}

  static SymmCmatrix Filled(const int _nrows, const int _ncols, const SCALAR t){
    SymmCmatrix M(_nrows,_ncols); std::fill(M.array,M.array+_nrows*_ncols,t); return M;}
  static SymmCmatrix Filled(const int _nrows, const SCALAR t){return Filled(_nrows,_nrows,t);}

  static SymmCmatrix Identity(const int nrows){
    SymmCmatrix M=Zero(nrows); for(int i=0; i<nrows; i++) M.array[i*nrows+i]=1; return M;}

  template<class VECTOR>
  static SymmCmatrix Diag(const VECTOR& v){
    SymmCmatrix M=Zero(v.n,v.n); 
    for(int i=0; i<v.n; i++) M.array[i*v.n+i]=v(i); 
    return M;
  }

  static SymmCmatrix Uniform(const int n) {return Uniform(n,n);}
  static SymmCmatrix Uniform(const int nrows, const int ncols){
    SymmCmatrix M(nrows,ncols);
    uniform_real_distribution<SCALAR> distr;
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<=i; j++){
	SCALAR t=distr(randomNumberGenerator); 
	M.array[i+j*nrows]=t; 
	M.array[j+i*nrows]=t;
      } 
    return M;
  }

  static SymmCmatrix Gaussian(const int n) {return Uniform(n,n);}
  static SymmCmatrix Gaussian(const int nrows, const int ncols){
    SymmCmatrix M(nrows,nrows);
    normal_distribution<SCALAR> distr;
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<=i; j++){
	SCALAR t=distr(randomNumberGenerator); 
	M.array[j*nrows+i]=t; 
	M.array[i*nrows+j]=t;
      }
    return M;
  }

  static SymmCmatrix Bernoulli(const int n, const double p=0.5) {return Bernoulli(n,n,p);}
  static SymmCmatrix Bernoulli(const int nrows, const int ncols, const double p=0.5){
    SymmCmatrix M(nrows,nrows);
    bernoulli_distribution distr(p);
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<=i; j++){
	SCALAR t=distr(randomNumberGenerator); 
	M.array[j*nrows+i]=t; 
	M.array[i*nrows+j]=t;
      }
    return M;
  }


public: // comparisons ----------------------------------------------------------------------------------------------


  bool operator==(const SymmCmatrix& X) const{ 
    if(X.nrows!=nrows) return false; 
    for(int j=0; j<ncols; j++) 
      for(int i=0; i<=j; i++) 
	if(array[j*nrows+i]!=X.array[j*nrows+i]) return false;
    return true;
  }


public: // element access -------------------------------------------------------------------------------------------

  SCALAR operator()(const int i, const int j) const{ 
    assert(i<nrows); assert(j<ncols); return array[j*nrows+i];} 
  SCALAR read(const int i, const int j) const{ 
    assert(i<nrows); assert(j<ncols); return array[j*nrows+i];} 

  SCALAR& operator()(const int i, const int j){
    SYMMETRY_UNSAFE("SymmCmatrix::operator(...)");
    assert(i<nrows); assert(j<ncols); return array[j*nrows+i];} 
  SCALAR& element(const int i, const int j){
    SYMMETRY_UNSAFE("SymmCmatrix::element(...)");
    assert(i<nrows); assert(j<ncols); return array[j*nrows+i];} 
  void set(const int i, const int j, const SCALAR v){ 
    assert(i<nrows); assert(j<ncols); array[j*nrows+i]=v; array[i*nrows+j]=v;} 


public: // iterators ------------------------------------------------------------------------------------------------

  void for_each(std::function<void(const INDEX, const INDEX, SCALAR&)> lambda){
    SYMMETRY_UNSAFE("SymmCmatrix::for_each(...)");
    Cmatrix::for_each(lambda);
  }

  void for_each_in_row(const int i, std::function<void(const INDEX, SCALAR&)> lambda){
    SYMMETRY_UNSAFE("SymmCmatrix::for_each_in_row(...)");
    Cmatrix::for_each_in_row(i,lambda);
  }

  void for_each_in_column(const int j, std::function<void(const INDEX, SCALAR&)> lambda){
    SYMMETRY_UNSAFE("SymmCmatrix::for_each_in_column(...)");
    Cmatrix::for_each_in_column(j,lambda);
  }


public: // views ----------------------------------------------------------------------------------------------------

  Detached<Cvector> view_of_column(const int j){
    SYMMETRY_UNSAFE("SymmCmatrix::view_of_column(...)");
    return Cmatrix::view_of_column(j);
  }

  Detached<const Cvector> view_of_column(const int j) const{
    Detached<const Cvector> v; v.n=nrows; v.array=array+j*nrows;
    return v;
  }

  Detached<Cmatrix> view_of_columns(const int j, const int k){
    SYMMETRY_UNSAFE("SymmCmatrix::view_of_columns(...)");
    return Cmatrix::view_of_columns(j, k);
  }

  
public: // remappings: inherited from Cmatrix with return value Cmatrix



public: // in-place operations --------------------------------------------------------------------------------------


  SymmCmatrix& operator+=(const SCALAR& x){
    for(int i=0; i<nrows*ncols; i++) array[i]+=x; return *this;}
  SymmCmatrix& operator-=(const SCALAR& x){
    for(int i=0; i<nrows*ncols; i++) array[i]-=x; return *this;}
  SymmCmatrix& operator*=(const SCALAR& x){
    for(int i=0; i<nrows*ncols; i++) array[i]*=x; return *this;}
  SymmCmatrix& operator/=(const SCALAR& x){
    for(int i=0; i<nrows*ncols; i++) array[i]/=x; return *this;}

  SymmCmatrix& operator+=(const SymmCmatrix& x){
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    for(int i=0; i<nrows*ncols; i++) array[i]+=x.array[i]; 
    return *this;}
  SymmCmatrix& operator-=(const SymmCmatrix& x){
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    for(int i=0; i<nrows*ncols; i++) array[i]-=x.array[i]; 
    return *this;}
  SymmCmatrix& operator*=(const SymmCmatrix& x){
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    for(int i=0; i<nrows*ncols; i++) array[i]*=x.array[i]; 
    return *this;}
  SymmCmatrix& operator/=(const SymmCmatrix& x){
    assert(x.nrows==nrows); assert(x.ncols==ncols);
    for(int i=0; i<nrows*ncols; i++) array[i]/=x.array[i]; 
    return *this;}


public: // construct operations -------------------------------------------------------------------------------------

  // gram inherited from Cmatrix


public: // vector valued arithmetic

  Cvector operator*(const Cvector& v) const{
    return dot(v);} // faster than *


public: // matrix valued arithmetic

  SymmCmatrix operator+(const SCALAR x) const{
    SymmCmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]+x; return M;}
  SymmCmatrix operator-(const SCALAR x) const{
    SymmCmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]-x; return M;}
  SymmCmatrix operator*(const SCALAR x) const{
    SymmCmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]*x; return M;}
  SymmCmatrix operator/(const SCALAR x) const{
    SymmCmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]/x; return M;}

  SymmCmatrix operator+(const SymmCmatrix B) const{ assert(B.nrows==nrows); assert(B.ncols==ncols);
    SymmCmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]+B.array[i]; return M;}
  SymmCmatrix operator-(const SymmCmatrix B) const{ assert(B.nrows==nrows); assert(B.ncols==ncols);
    SymmCmatrix M(nrows,ncols); for(int i=0; i<nrows*ncols; i++) M.array[i]=array[i]+B.array[i]; return M;}


public: // in-place operations: inherited from Cmatrix with return value Cmatrix

public: // scalar-valued methods

  SCALAR inp_of_rows(const int i1, const int i2) const{
    return Cmatrix::inp_of_columns(i1,i2);} // may be faster


public: // Givens rotations

  SymmCmatrix& applyFromLeft(const GivensRotation& Q){
    SYMMETRY_UNSAFE("SymmCmatrix::applyFromLeft(const GivensRotation& )");
    Cmatrix::applyFromLeft(Q);
    return *this;
  }

  SymmCmatrix& applyFromLeftT(const GivensRotation& Q){
    SYMMETRY_UNSAFE("SymmCmatrix::applyFromLeftT(const GivensRotation& )");
    Cmatrix::applyFromLeftT(Q);
    return *this;
  }

  SymmCmatrix& applyFromRight(const GivensRotation& Q){
    SYMMETRY_UNSAFE("SymmCmatrix::applyFromRight(const GivensRotation& )");
    Cmatrix::applyFromRight(Q);
    return *this;
  }

  SymmCmatrix& applyFromRightT(const GivensRotation& Q){
    SYMMETRY_UNSAFE("SymmCmatrix::applyFromRightT(const GivensRotation& )");
    Cmatrix::applyFromRightT(Q);
    return *this;
  }

  SymmCmatrix& conjugate(const GivensRotation& Q){ // improve!
    assert(Q.i1<nrows);
    assert(Q.i2<nrows);
    for(int j=0; j<ncols; j++){
      SCALAR v1=array[j*nrows+Q.i1];
      SCALAR v2=array[j*nrows+Q.i2];
      array[j*nrows+Q.i1]=Q.cos*v1-Q.sin*v2;
      array[j*nrows+Q.i2]=Q.cos*v2+Q.sin*v1;}
    for(int i=0; i<nrows; i++){
      SCALAR v1=array[Q.i1*nrows+i];
      SCALAR v2=array[Q.i2*nrows+i];
      array[Q.i1*nrows+i]=Q.cos*v1-Q.sin*v2;
      array[Q.i2*nrows+i]=Q.cos*v2+Q.sin*v1;}
    return *this;
  }


  SymmCmatrix& conjugateT(const GivensRotation& Q){ // improve!
    assert(Q.i1<nrows);
    assert(Q.i2<nrows);
    for(int j=0; j<ncols; j++){
      SCALAR v1=array[j*nrows+Q.i1];
      SCALAR v2=array[j*nrows+Q.i2];
      array[j*nrows+Q.i1]=Q.cos*v1+Q.sin*v2;
      array[j*nrows+Q.i2]=Q.cos*v2-Q.sin*v1;}
    for(int i=0; i<nrows; i++){
      SCALAR v1=array[Q.i1*nrows+i];
      SCALAR v2=array[Q.i2*nrows+i];
      array[Q.i1*nrows+i]=Q.cos*v1+Q.sin*v2;
      array[Q.i2*nrows+i]=Q.cos*v2-Q.sin*v1;}
    return *this;
  }



public: // KpointOp<k> ----------------------------------------------------------------------------------------------


  template<int k>
  SymmCmatrix& applyFromLeft(const KpointOp<k>& Q){
    SYMMETRY_UNSAFE("SymmCmatrix::applyFromLeft(const KpointOp<k>&)");
    Cmatrix::applyFromLeft(Q);
    return *this;
  }

  template<int k>
  SymmCmatrix& applyFromLeftT(const KpointOp<k>& Q){
    SYMMETRY_UNSAFE("SymmCmatrix::applyFromLeftT(const KpointOp<k>&)");
    Cmatrix::applyFromLeftT(Q);
    return *this;
  }

  template<int k>
  SymmCmatrix& applyFromRight(const KpointOp<k>& Q){
    SYMMETRY_UNSAFE("SymmCmatrix::applyFromRight(const KpointOp<k>&)");
    Cmatrix::applyFromRight(Q);
    return *this;
  }

  template<int k>
  SymmCmatrix& applyFromRightT(const KpointOp<k>& Q){
    SYMMETRY_UNSAFE("SymmCmatrix::applyFromRightT(const KpointOp<k>&)");
    Cmatrix::applyFromRightT(Q);
    return *this;
  }

  template<int k>
  SymmCmatrix& conjugate(const KpointOp<k>& Q){ // improve!
    Cmatrix::applyFromLeft(Q);
    Cmatrix::applyFromRightT(Q);
    return *this;
  }

  template<int k>
  SymmCmatrix& conjugateT(const KpointOp<k>& Q){ // improve!
    Cmatrix::applyFromLeftT(Q);
    Cmatrix::applyFromRight(Q);
    return *this;
  }


public: // Python interface

  SymmCmatrix(double* numpyDblArray, int numpyDim1, int numpyDim2): 
    Cmatrix(numpyDim1, numpyDim2){
    std::copy(numpyDblArray,numpyDblArray+nrows*ncols,array);
  }



};


  // wrapper allowing a Cmatrix to pose as a SymmCmatrix but only locally 
  class AsSymmCmatrix: public SymmCmatrix{
  public:
    AsSymmCmatrix()=delete;
    AsSymmCmatrix& operator=(const AsSymmCmatrix& x)=delete;
    AsSymmCmatrix(Cmatrix& x): SymmCmatrix(x.shallow(),_Downcast()){}
    AsSymmCmatrix(const Cmatrix& x): SymmCmatrix(x.shallow(),_Downcast()){}
    ~AsSymmCmatrix(){array=nullptr;}
  };

  inline SymmCmatrix::SymmCmatrix(const AsSymmCmatrix& x): 
    Cmatrix(x,_NoWarn()) {
    COPY_WARNING("SymmCmatrix");
  }

  inline AsSymmCmatrix as_symmetric(Cmatrix& x){
    return AsSymmCmatrix(x);}
  inline const AsSymmCmatrix as_symmetric(const Cmatrix& x){
    return AsSymmCmatrix(x);}
  inline SymmCmatrix as_symmetric(Cmatrix&& x){
    MOVE_WARNING("Cmatrix");
    return SymmCmatrix(std::move(x),_Downcast());}


} // namespace Mondrian

#endif



  //operator SymmCmatrix(){return Cmatrix(*this,_NoWarn());}

  //These cannot be used if inheriting from Cmatrix 
  //biscalar operator()(const int i, const short int j){
  //  assert(i<nrows); assert(j<ncols); return biscalar(array[i*nrows+j],array[j*nrows+i]);} 
  //biscalar element(const int i, const short int j){
  //  assert(i<nrows); assert(j<ncols); return biscalar(array[i*nrows+j],array[j*nrows+i]);} 

  //SymmCmatrix(const MatrixIsSymmetric<Cmatrix>& x): Cmatrix(x.obj,_NoWarn()){
  //  CONVERT_WARNING("Cmatrix (symmetric)","SymmCmatrix");}

  //SymmCmatrix(MatrixIsSymmetricTemp<Cmatrix> x): Cmatrix(std::move(x.obj),_NoWarn()){
  //  MOVECONVERT_WARNING("Cmatrix (symmetric)","SymmCmatrix");}


  //template<class MATRIX>
  //SymmCmatrix(MATRIX&& x): Cmatrix(std::move(x),_NoWarn()){
  //  MOVECONVERT_WARNING("MATRIX","SymmCmatrix");
  //  symmetrize();
  //}

  /*
  template<class MATRIX>
  SymmCmatrix(const MatrixIsSymmetric<MATRIX>& x): Cmatrix(x.obj,_NoWarn()){
    CONVERT_WARNING("MATRIX","SymmCmatrix");}

  template<class MATRIX>
  SymmCmatrix(MatrixIsSymmetricTemp<MATRIX>& x): Cmatrix(std::move(x.obj),_NoWarn()){
    MOVECONVERT_WARNING("MATRIX","SymmCmatrix");}
  */

  /*
class DetachedSymmCmatrixRvalue: public SymmCmatrix{
public:

  DetachedSymmCmatrixRvalue()=delete;
  DetachedSymmCmatrixRvalue& operator=(const DetachedSymmCmatrixRvalue& x)=delete;

  DetachedSymmCmatrixRvalue(Cmatrix&& x): SymmCmatrix(x,_Shallow()){
    DOWNCAST_WARNING("Cmatrix&&","DetachedSymmCmatrixRvalue");}

  ~DetachedSymmCmatrixRvalue(){array=nullptr;}
};
  */

  //inline DetachedSymmCmatrixRvalue as_symmetric(Cmatrix&& x){
  //return DetachedSymmCmatrixRvalue(std::move(x));}
  //inline SymmCmatrix::SymmCmatrix(const DetachedSymmCmatrixRvalue& x): Cmatrix(x,_NoWarn()){
  //MOVEINTO_WARNING("DetachedSymmCmatrixRvalue","SymmCmatrix");
  //}
  //inline SymmCmatrix::SymmCmatrix(DetachedSymmCmatrixRvalue&& x): Cmatrix(x,_NoWarn()){
  //MOVEINTO_WARNING("DetachedSymmCmatrixRvalue","SymmCmatrix");
  //}

  //SymmCmatrix(const DetachedSymmCmatrixRvalue& x);
  //SymmCmatrix(DetachedSymmCmatrixRvalue&& x);

  //inline SymmCmatrix::SymmCmatrix(const DetachedSymmCmatrix& x): Cmatrix(x,_NoWarn()){
  //COPYINTO_WARNING("DetachedSymmCmatrix","SymmCmatrix");
  //}

  //SymmCmatrix(const Cmatrix& x, const _Shallow dummy): 
  //  Cmatrix(x.nrows,x.ncols,_Uninitialized()){
  //  array=x.array;
  //}


  //SymmCmatrix(AsConstSymmCmatrix&& x);
  //class AsConstSymmCmatrix;  
  /*
  class AsConstSymmCmatrix: public SymmCmatrix{
  public:
    AsConstSymmCmatrix()=delete;
    AsConstSymmCmatrix& operator=(const AsConstSymmCmatrix& x)=delete;
    AsConstSymmCmatrix(const Cmatrix& x): SymmCmatrix(x.shallow(),_Downcast()){}
    ~AsConstSymmCmatrix(){array=nullptr;}
  };
  */
  //inline SymmCmatrix::SymmCmatrix(AsConstSymmCmatrix&& x): Cmatrix(x,_NoWarn()) {COPY_WARNING("SymmCmatrix");}
 //SymmCmatrix(const Cmatrix& x, const _ConstDowncast dummy): 
  //  Cmatrix(x,_NoWarn()){
  //  DOWNCASTCOPY_WARNING("const Cmatrix","const SymmCmatrix");
  //}

  //SymmCmatrix(Cmatrix&& x, const _ConstDowncast dummy): 
  //  Cmatrix(std::move(x),_NoWarn()){
  //  DOWNCAST_WARNING("const Cmatrix","SymmCmatrix");
  //}
  //SymmCmatrix(AsSymmCmatrix&& x);
  //inline SymmCmatrix::SymmCmatrix(AsSymmCmatrix&& x): 
  //  Cmatrix(x,_NoWarn()){
  //  COPY_WARNING("SymmCmatrix");
  //}

    //operator SymmCmatrix(){return SymmCmatrix(Cmatrix(*this),_Downcast());} would never be used
    //AsSymmCmatrix(Cmatrix& x): SymmCmatrix(x,_Shallow()){} // below is simpler


/*

public: // Gram matrices

  static SymmCmatrix Gram(const Cmatrix& A){
    int ncols=A.ncols; int nrows=A.nrows;
    SymmCmatrix G(ncols,ncols);
    for(int i=0; i<ncols; i++)
      for(int j=0; j<=i; j++){
	SCALAR t=0;
	for(int k=0; k<nrows; k++)
	  t+=A.array[i*nrows+k]*A.array[j*nrows+k];
	G.set(i,j,t);
      }
    return G;
  }

  static SymmCmatrix Gram(const MatrixTranspose<const Cmatrix>& TA){
    const Cmatrix& A=TA.obj;
    int ncols=A.ncols; int nrows=A.nrows;
    SymmCmatrix G(nrows,nrows);
    for(int i=0; i<nrows; i++)
      for(int j=0; j<=i; j++){
	SCALAR t=0;
	for(int k=0; k<ncols; k++)
	  t+=A.array[k*nrows+i]*A.array[k*nrows+j];
	G.set(i,j,t);
      }
    return G;
  }

  static SymmCmatrix Gram(const MatrixTranspose<const SymmCmatrix>& TA){
    return Gram(TA.obj);
  } // may be faster

  template<class VECTOR>
  static SymmCmatrix Gram(const MatrixX<VECTOR>& A){
    int ncols=A.ncols;
    SymmCmatrix G(ncols,ncols);
    for(int i=0; i<ncols; i++)
      for(int j=0; j<=i; j++)
	G.set(i,j,A.inp_of_columns(i,j));
    return G;
  }

  template<class VECTOR>
  static SymmCmatrix Gram(const MatrixTranspose<const MatrixX<VECTOR> >& TA){
    auto A=TA.obj;
    int nrows=A.nrows;
    SymmCmatrix G(nrows,nrows);
    for(int i=0; i<nrows; i++)
      for(int j=0; j<=i; j++)
	G.set(i,j,A.inp_of_rows(i,j));
    return G;
  }

  template<class VECTOR>
  static SymmCmatrix Gram(const MatrixTranspose<const SymmMatrixX<VECTOR> >& TA){
    return Gram(TA.obj);}

*/
  /*
public: // KpointOperator

  SymmCmatrix& applyFromLeft(const KpointOperator& Q){
    SYMMETRY_UNSAFE("SymmCmatrix::applyFromLeft(const KpointOperator&)");
    Cmatrix::applyFromLeft(Q);
    return *this;
  }

  SymmCmatrix& applyFromLeftT(const KpointOperator& Q){
    SYMMETRY_UNSAFE("SymmCmatrix::applyFromLeftT(const KpointOperator&)");
    Cmatrix::applyFromLeftT(Q);
    return *this;
  }

  SymmCmatrix& applyFromRight(const KpointOperator& Q){
    SYMMETRY_UNSAFE("SymmCmatrix::applyFromRight(const KpointOperator&)");
    Cmatrix::applyFromRight(Q);
    return *this;
  }

  SymmCmatrix& applyFromRightT(const KpointOperator& Q){
    SYMMETRY_UNSAFE("SymmCmatrix::applyFromRightT(const KpointOperator&)");
    Cmatrix::applyFromRightT(Q);
    return *this;
  }

  SymmCmatrix& conjugate(const KpointOperator& Q){ // improve!
    Cmatrix::applyFromLeft(Q);
    Cmatrix::applyFromRightT(Q);
    return *this;
  }

  SymmCmatrix& conjugateT(const KpointOperator& Q){ // improve!
    Cmatrix::applyFromLeftT(Q);
    Cmatrix::applyFromRight(Q);
    return *this;
  }
  */
