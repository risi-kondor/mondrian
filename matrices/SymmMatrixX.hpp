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


#ifndef _SymmMatrixX
#define _SymmMatrixX

#include "MatrixX.hpp"
#include "SymmCmatrix.hpp"

extern default_random_engine randomNumberGenerator;


namespace Mondrian{

  template<class VECTOR>
  class AsSymmMatrixX;


  template<class VECTOR>
  class SymmMatrixX: public MatrixX<VECTOR>{
  public:

    using MatrixX<VECTOR>::nrows;
    using MatrixX<VECTOR>::ncols;
    using MatrixX<VECTOR>::cols;

    using MatrixX<VECTOR>::symmetrize;
    using MatrixX<VECTOR>::operator();


  public: // attributes

    bool isSymmetricFormat() const {return true;}


  public: // constructors --------------------------------------------------------------------------------------------


    using MatrixX<VECTOR>::MatrixX;

    SymmMatrixX(): SymmMatrixX(0,0){}

    SymmMatrixX(const int _nrows): MatrixX<VECTOR>(_nrows,_nrows){}

    SymmMatrixX(const int _nrows, const int _ncols): MatrixX<VECTOR>(_nrows,_ncols){
      assert(_nrows==_ncols);}

    SymmMatrixX(const int _nrows, const int _ncols, const _Zero dummy): 
      MatrixX<VECTOR>(_nrows,_ncols,_Zero()){
      assert(_nrows==_ncols);
    }

    SymmMatrixX(const initializer_list<Cvector>& list): MatrixX<VECTOR>(list){
      symmetrize();}
  
    SymmMatrixX(const int _nrows, const int _ncols, const initializer_list<iivtriple>& list): 
      MatrixX<VECTOR>(_nrows,_ncols,list){
      symmetrize();
    }


  public: // copying -------------------------------------------------------------------------------------------------


    SymmMatrixX(const SymmMatrixX<VECTOR>& x): MatrixX<VECTOR>(x,_NoWarn()){
      COPY_WARNING("SymmMatrixX");
    }

    SymmMatrixX(SymmMatrixX<VECTOR>&& x): MatrixX<VECTOR>(std::move(x)){
      MOVE_WARNING("SymmMatrixX");
    }

    SymmMatrixX<VECTOR>& operator=(const SymmMatrixX<VECTOR>& x){
      ASSIGN_WARNING("SymmMatrixX");
      nrows=x.nrows; ncols=x.ncols;
      for(auto p:cols) delete p; cols.resize(ncols); 
      for(int j=0; j<ncols; j++) cols[j]=new VECTOR(*x.cols[j]);
      return *this;
    }

    SymmMatrixX<VECTOR>& operator=(SymmMatrixX<VECTOR>&& x){
      MOVEASSIGN_WARNING("SymmMatrixX");
      nrows=x.nrows; ncols=x.ncols;
      for(auto p:cols) delete p; 
      cols=x.cols; x.nrows=0; x.ncols=0;
      return *this;
    }

    SymmMatrixX<VECTOR> copy(){
      SymmMatrixX M(nrows,ncols);
      for(int j=0; j<ncols; j++) M.cols[j]=new VECTOR(cols[j]->copy());
      return *this;
    }

    void detach(){cols.clear();}

    void assign(const SymmMatrixX& x){
      assert(x.nrows==nrows);
      assert(x.ncols==ncols);
      for(int j=0; j<ncols; j++) *cols[j]=*x.cols[j]; 
    }

    SymmMatrixX shallow() const{
      SymmMatrixX M(nrows,ncols); 
      M.cols=cols;
      return M;
    }


  public: // downcasting --------------------------------------------------------------------------------------------


    SymmMatrixX(const MatrixX<VECTOR>& x, const _Downcast dummy): 
      MatrixX<VECTOR>(x,_NoWarn()){
      DOWNCASTCOPY_WARNING("MatrixX<VECTOR>","SymmMatrixX<VECTOR>");
    }

    SymmMatrixX(MatrixX<VECTOR>&& x, const _Downcast dummy): 
      MatrixX<VECTOR>(std::move(x),_NoWarn()){
      DOWNCAST_WARNING("MatrixX<VECTOR>","SymmMatrixX<VECTOR>");
    }


  public: // conversions from other matrix classes ------------------------------------------------------------------


    SymmMatrixX(const AsSymmMatrixX<VECTOR>& x);

    template<class VECTOR2>
    SymmMatrixX(MatrixX<VECTOR2>& x): 
      MatrixX<VECTOR>(x.nrows,x.ncols,_Zero()){
      assert(nrows==ncols);
      CONVERT_WARNING("MatrixX<VECTOR2>","SymmMatrixX<VECTOR>");
      x.for_each_filled([this](const INDEX i, const INDEX j, const SCALAR v){set(i,j,v);});
    }

    SymmMatrixX(MatrixX<VECTOR>&& x): 
      MatrixX<VECTOR>(std::move(x),_NoWarn()){
      MOVECONVERT_WARNING("MatrixX<VECTOR>","SymmMatrixX<VECTOR>");
      symmetrize();
    }

    SymmMatrixX(const Cmatrix& x): 
      MatrixX<VECTOR>(x,_NoWarn()){
      CONVERT_WARNING("Cmatrix","SymmMatrixX<VECTOR>");
      if(!x.isSymmetricFormat()) symmetrize();
    }

    /* overrides other conversions 
       template<class MATRIX>
       SymmMatrixX(const MATRIX& x): 
       MatrixX<VECTOR>(x){
       CONVERT_WARNING("MATRIX","SymmMatrixX<VECTOR>");
       if(!x.isSymmetricFormat()) symmetrize();
       }
    */


  public: // conversions to other matrix classes --------------------------------------------------------------------


    operator Cmatrix(){
      CONVERT_WARNING("SymmMatrixX<VECTOR>","Cmatrix");
      Cmatrix M=Cmatrix::Zero(MatrixX<VECTOR>::nrows,MatrixX<VECTOR>::ncols);
      for_each_filled([&M](const int i, const int j, const SCALAR v){M(i,j)=v;});
      return M;
    }  

    operator SymmCmatrix(){
      CONVERT_WARNING("SymmMatrixX<VECTOR>","SymmCmatrix");
      SymmCmatrix M=Cmatrix::Zero(MatrixX<VECTOR>::nrows,MatrixX<VECTOR>::ncols);
      for_each_filled([&M](const int i, const int j, const SCALAR v){M(i,j)=v;});
      return M;
    }  

    template<class VECTOR2>
    operator MatrixX<VECTOR2>(){
      CONVERT_WARNING("SymmMatrixX<VECTOR>","MatrixX<VECTOR2>");
      MatrixX<VECTOR2> M(MatrixX<VECTOR>(*this),_NoWarn());
      return M;
    }  


  public: // named constructors -------------------------------------------------------------------------------------


    static SymmMatrixX Zero(const int _nrows, const int _ncols)
    {return SymmMatrixX<VECTOR>(_nrows,_ncols,_Zero());}

    static SymmMatrixX Filled(const int _nrows, const int _ncols, const SCALAR v){
      SymmMatrixX<VECTOR> M(_nrows,_ncols);
      for(int j=0; j<_ncols; j++) (*M.cols[j])=VECTOR::Filled(_nrows,v); 
      return M;
    }

    static SymmMatrixX Identity(const int n){
      SymmMatrixX<VECTOR> M(n,n,_Zero()); //=SymmMatrixX<VECTOR>::Zero(n,n);
      for(int i=0; i<n; i++) (*M.cols[i])(i)=1;
      return M;
    }

    static SymmMatrixX Uniform(const int _nrows, const int _ncols){
      SymmMatrixX<VECTOR> M(_nrows,_ncols);
      uniform_real_distribution<SCALAR> distr;
      for(int i=0; i<_nrows; i++)
	for(int j=0; j<=i; j++)
	  M.set(i,j,distr(randomNumberGenerator));
      return M;
    }

    static SymmMatrixX Gaussian(const int _nrows, const int _ncols){
      SymmMatrixX<VECTOR> M(_nrows,_ncols);
      normal_distribution<SCALAR> distr;
      for(int i=0; i<_nrows; i++)
	for(int j=0; j<=i; j++)
	  M.set(i,j,distr(randomNumberGenerator));
      return M;
    }

    static SymmMatrixX Bernoulli(const int _nrows, const int _ncols, const double p=0.5){
      SymmMatrixX<VECTOR> M(_nrows,_ncols);
      bernoulli_distribution distr(p);
      for(int i=0; i<_nrows; i++)
	for(int j=0; j<=i; j++)
	  if(distr(randomNumberGenerator)==1) 
	    M.set(i,j,1);
      return M;
    }


  public: // element access -----------------------------------------------------------------------------------------


    SCALAR& operator()(const int i, const int j){
      SYMMETRY_UNSAFE("SymmMatrixX<VECTOR>::operator(...)");
      assert(i<nrows); assert(j<ncols); return (*cols[j])(i);} 
    SCALAR& element(const int i, const int j){
      SYMMETRY_UNSAFE("SymmMatrixX<VECTOR>::element(...)");
      assert(i<nrows); assert(j<ncols); return (*cols[j])(i);} 
    void set(const int i, const int j, const SCALAR v){ 
      assert(i<nrows); assert(j<ncols); (*cols[j])(i)=v; (*cols[i])(j)=v;} 


  public: // iterators

    void for_each_filled(std::function<void(const INDEX, const INDEX, SCALAR&)> lambda){
      SYMMETRY_UNSAFE("SymmMatrixX<VECTOR>::for_each(...)");
      MatrixX<VECTOR>::for_each_filled(lambda);
    }

    void for_each_filled_in_row(const int i, std::function<void(const INDEX, SCALAR&)> lambda){
      SYMMETRY_UNSAFE("SymmMatrixX<VECTOR>::for_each_in_row(...)");
      MatrixX<VECTOR>::for_each_filled_in_row(i,lambda);
    }

    void for_each_filled_in_column(const int j, std::function<void(const INDEX, SCALAR&)> lambda){
      SYMMETRY_UNSAFE("SymmMatrixX<VECTOR>::for_each_in_column(...)");
      MatrixX<VECTOR>::for_each_filled_in_column(j,lambda);
    }


  public: // in-place operations ---------------------------------------------------------------------------------

    // operator+= inherited
    // operator-= inherited

    template<class VECTOR2>
    SymmMatrixX<VECTOR>& multiplyRowsBy(const VECTOR2& v){
      SYMMETRY_UNSAFE("SymmMatrixX<VECTOR>::multiplyRowsBy(...)");
      MatrixX<VECTOR>::multiplyRowsBy(v);
      return *this;
    }

    template<class VECTOR2>
    SymmMatrixX<VECTOR>& multiplyColumnsBy(const VECTOR2& v){
      SYMMETRY_UNSAFE("SymmMatrixX<VECTOR>::multiplyRowsBy(...)");
      MatrixX<VECTOR>::multiplyColumnsBy(v);
      return *this;
    }

   template<class VECTOR2>
    SymmMatrixX<VECTOR>& divideRowsBy(const VECTOR2& v){
      SYMMETRY_UNSAFE("SymmMatrixX<VECTOR>::divideRowsBy(...)");
      MatrixX<VECTOR>::divideRowsBy(v);
      return *this;
    }

    template<class VECTOR2>
    SymmMatrixX<VECTOR>& divideColumnsBy(const VECTOR2& v){
      SYMMETRY_UNSAFE("SymmMatrixX<VECTOR>::divideRowsBy(...)");
      MatrixX<VECTOR>::divideColumnsBy(v);
      return *this;
    }


  public: // construct operations 

    //SymmMatrixX<VECTOR>(xgram())

 


  public: // Grams

  





  public: // scalar valued operations

    SCALAR inp_of_rows(const int i1, const int i2) const{
      return  MatrixX<VECTOR>::inp_of_columns(i1,i2);}


  public: // Givens rotations ---------------------------------------------------------------------------------------


    SymmMatrixX<VECTOR>& conjugate(const GivensRotation& Q){
      applyFromLeft(Q);
      applyFromRightT(Q);
      return *this;
    }

    SymmMatrixX<VECTOR>& conjugateT(const GivensRotation& Q){
      applyFromLeftT(Q);
      applyFromRight(Q);
      return *this;
    }


  public: // KpointOp<k> --------------------------------------------------------------------------------------------


    template<int k>
    SymmMatrixX<VECTOR>& conjugate(const KpointOp<k>& Q){
      VectorSpattern spattern;
      for(int i=0; i<k; i++){
	assert(Q.map(i)<ncols);
	cols[Q.map(i)]->apply(Q);
	spattern.add(*cols[Q.map(i)]);
      }
      SCALAR V[spattern.filled.size()*k];
      int j=0;
      for(auto r:spattern.filled){
	SCALAR* vptr[k];
	for(int i=0; i<k; i++)
	  vptr[i]=cols[Q.map(i)]->ptr(r);
	Q.applyTo(vptr);
	for(int i=0; i<k; i++)
	  V[j*k+i]=*vptr[i];
	j++;
      }
      j=0;
      for(auto r:spattern.filled){
	//cout<<r<<endl;
	for(int i=0; i<k; i++)
	  cols[r]->set(Q.map(i),V[j*k+i]);
	j++;
      }
      return *this;
    }

    template<int k>
    SymmMatrixX<VECTOR>& conjugateT(const KpointOp<k>& Q){
      VectorSpattern spattern;
      for(int i=0; i<k; i++){
	assert(Q.map(i)<ncols);
	cols[Q.map(i)]->applyT(Q);
	spattern.add(*cols[Q.map(i)]);
      }
      SCALAR V[spattern.filled.size()*k];
      int j=0;
      for(auto r:spattern.filled){
	SCALAR* vptr[k];
	for(int i=0; i<k; i++)
	  vptr[i]=cols[Q.map(i)]->ptr(r);
	Q.applyToT(vptr);
	for(int i=0; i<k; i++)
	  V[j*k+i]=*vptr[i];
	j++;
      }
      j=0;
      for(auto r:spattern.filled){
	for(int i=0; i<k; i++)
	  cols[r]->set(Q.map(i),V[j*k+i]);
	j++;
      }
      return *this;
    }


  private: // disabled methods --------------------------------------------------------------------------------------


    using MatrixX<VECTOR>::multiplyRowsBy;
    using MatrixX<VECTOR>::multiplyColumnsBy;
    using MatrixX<VECTOR>::divideRowsBy;
    using MatrixX<VECTOR>::divideColumnsBy;

    using MatrixX<VECTOR>::applyFromLeft;
    using MatrixX<VECTOR>::applyFromLeftT;
    using MatrixX<VECTOR>::applyFromRight;
    using MatrixX<VECTOR>::applyFromRightT;
    


};



  template<class VECTOR>
  class AsSymmMatrixX: public SymmMatrixX<VECTOR>{
  public:
    AsSymmMatrixX()=delete;
    AsSymmMatrixX& operator=(const AsSymmMatrixX<VECTOR>& x)=delete;
    AsSymmMatrixX(MatrixX<VECTOR>& x): SymmMatrixX<VECTOR>(x.shallow(),_Downcast()){}
    AsSymmMatrixX(const MatrixX<VECTOR>& x): SymmMatrixX<VECTOR>(x.shallow(),_Downcast()){}
    ~AsSymmMatrixX(){SymmMatrixX<VECTOR>::cols.clear();}
  };

  template<class VECTOR>
  SymmMatrixX<VECTOR>::SymmMatrixX(const AsSymmMatrixX<VECTOR>& x): 
    MatrixX<VECTOR>(x,_NoWarn()) {
    COPY_WARNING("SymmMatrixX<VECTOR>");
  }


  template<class VECTOR>
  AsSymmMatrixX<VECTOR> as_symmetric(MatrixX<VECTOR>& x){
    return AsSymmMatrixX<VECTOR>(x);}

  template<class VECTOR>
  const AsSymmMatrixX<VECTOR> as_symmetric(const MatrixX<VECTOR>& x){
    return AsSymmMatrixX<VECTOR>(x);}

  template<class VECTOR>
  SymmMatrixX<VECTOR> as_symmetric(MatrixX<VECTOR>&& x){
    MOVE_WARNING("MatrixX");
    return SymmMatrixX<VECTOR>(std::move(x),_Downcast());}



} // namespace Mondrian

#endif
  /*
  SymmMatrixX(const Cmatrix& x):
    SymmMatrixX<VECTOR>(x.nrows,x.ncols,_Zero()){
    CONVERT_WARNING("Cmatrix","SymmMatrixX");
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	if(x.array[j*nrows+i]!=0) cols[j]->append(i,x.array[j*nrows+i]);
    if(!x.isSymmetricFormat()) symmetrize();
  }
  */
  
  //SymmMatrixX(const MatrixX<VECTOR>& x): MatrixX(x,_NoWarn()){
  //  CONVERT_WARNING("MatrixX","SymmMatrixX");
  //  symmetrize();
  //}

    /*
    template<int k>
    SymmMatrixX<VECTOR>& applyFromLeft(const KpointOp<k>& Q){
      SYMMETRY_FORBIDDEN("applyFromLeft(const KpointOp<k>&)");
      return *this;
    }

    template<int k>
    SymmMatrixX<VECTOR>& applyFromLeftT(const KpointOp<k>& Q){
      SYMMETRY_FORBIDDEN("applyFromLeftT(const KpointOp<k>&)");
      return *this;
    }

    template<int k>
    SymmMatrixX<VECTOR>& applyFromRightT(const KpointOp<k>& Q){
      SYMMETRY_FORBIDDEN("applyFromRightT(const KpointOp<k>&)");
      return *this;
    }

    template<int k>
    SymmMatrixX<VECTOR>& applyFromRight(const KpointOp<k>& Q){
      SYMMETRY_FORBIDDEN("applyFromRight(const KpointOp<k>&)");
      return *this;
    }
    */
