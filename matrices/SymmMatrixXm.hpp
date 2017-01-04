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


#ifndef _SymmMatrixXm
#define _SymmMatrixXm

#include <set>

#include "MatrixXm.hpp"
#include "SymmCmatrix.hpp"
#include "Flushed.hpp"

extern default_random_engine randomNumberGenerator;


namespace Mondrian{

  template<class VECTOR>
  class AsSymmMatrixXm;


  template<class VECTOR>
  class SymmMatrixXm: public MatrixXm<VECTOR>{
  public:

    set<INDEX> changed_cols;
    mutex changed_cols_mx;

    using MatrixXm<VECTOR>::nrows;
    using MatrixXm<VECTOR>::ncols;
    using MatrixXm<VECTOR>::cols;
    using MatrixXm<VECTOR>::column_manager;

    using MatrixXm<VECTOR>::symmetrize;
    using MatrixXm<VECTOR>::operator();

    bool isSymmetricFormat() const {return true;}


  public: // constructors -------------------------------------------------------------------------------------------


    SymmMatrixXm(): SymmMatrixXm(0,0){}

    SymmMatrixXm(const int _nrows): SymmMatrixXm<VECTOR>(_nrows,_nrows){}

    SymmMatrixXm(const int _nrows, const int _ncols): MatrixXm<VECTOR>(_nrows,_ncols){
      assert(_nrows==_ncols);}

    SymmMatrixXm(const int _nrows, const int _ncols, const _Zero dummy): 
      MatrixXm<VECTOR>(_nrows,_ncols,_Zero()){
      assert(_nrows==_ncols);
    }

    SymmMatrixXm(const initializer_list<Cvector>& list): 
      MatrixXm<VECTOR>(list){
      symmetrize();
    }
  
    SymmMatrixXm(const int _nrows, const int _ncols, const initializer_list<iivtriple>& list): 
      MatrixXm<VECTOR>(_nrows,_ncols,list){
      symmetrize();
    }


  public: // copying ------------------------------------------------------------------------------------------------

    // note: currently source is not flushed automatically

    SymmMatrixXm(const SymmMatrixXm<VECTOR>& x): 
      MatrixXm<VECTOR>(x){
      COPY_WARNING("SymmMatrixXm");
    }

    SymmMatrixXm(SymmMatrixXm<VECTOR>&& x): 
      MatrixXm<VECTOR>(std::move(x)){
      MOVE_WARNING("SymmMatrixXm");
    }

    SymmMatrixXm<VECTOR>& operator=(const SymmMatrixXm<VECTOR>& x){
      ASSIGN_WARNING("SymmMatrixXm");
      nrows=x.nrows; ncols=x.ncols;
      for(auto p:cols) delete p; cols.resize(ncols); 
      for(int j=0; j<ncols; j++) cols[j]=new VECTOR(*x.cols[j]);
      return *this;
    }

    SymmMatrixXm<VECTOR>& operator=(SymmMatrixXm<VECTOR>&& x){
      MOVEASSIGN_WARNING("SymmMatrixXm");
      nrows=x.nrows; ncols=x.ncols;
      for(auto p:cols) delete p; 
      cols=x.cols; x.nrows=0; x.ncols=0;
      return *this;
    }

    SymmMatrixXm<VECTOR> copy(){
      SymmMatrixXm M(nrows,ncols);
      for(int j=0; j<ncols; j++) M.cols[j]=new VECTOR(cols[j]->copy());
      return *this;
    }

    void detach(){cols.clear();}

    void assign(const SymmMatrixXm& x){
      assert(x.nrows==nrows);
      assert(x.ncols==ncols);
      for(int j=0; j<ncols; j++) *cols[j]=*x.cols[j]; 
    }

    SymmMatrixXm shallow() const{
      SymmMatrixXm M(nrows,ncols); 
      M.cols=cols;
      return M;
    }


  public: // downcasting -------------------------------------------------------------------------------------------


    SymmMatrixXm(const MatrixX<VECTOR>& x): 
      MatrixXm<VECTOR>(x){
      CONVERT_WARNING("MatrixX<VECTOR>","SymmMatrixXm<VECTOR>");
    }

    SymmMatrixXm(MatrixX<VECTOR>&& x): 
      MatrixXm<VECTOR>(std::move(x)){
      CONVERT_WARNING("MatrixX<VECTOR>","SymmMatrixXm<VECTOR>");
    }

    SymmMatrixXm(const MatrixXm<VECTOR>& x, const _Downcast dummy): 
      MatrixXm<VECTOR>(x,_NoWarn()){
      DOWNCASTCOPY_WARNING("MatrixXm<VECTOR>","SymmMatrixXm<VECTOR>");
    }

    SymmMatrixXm(MatrixXm<VECTOR>&& x, const _Downcast dummy): 
      MatrixXm<VECTOR>(std::move(x),_NoWarn()){
      DOWNCAST_WARNING("MatrixXm<VECTOR>","SymmMatrixXm<VECTOR>");
    }


  public: // conversions from other matrix classes ------------------------------------------------------------------


    SymmMatrixXm(const AsSymmMatrixXm<VECTOR>& x);
    
    template<class VECTOR2>
    SymmMatrixXm(MatrixX<VECTOR2>& x): 
      MatrixXm<VECTOR>(x.nrows,x.ncols,_Zero()){
      assert(nrows==ncols);
      CONVERT_WARNING("MatrixX<VECTOR2>","SymmMatrixXm<VECTOR>");
      x.for_each_filled([this](const INDEX i, const INDEX j, const SCALAR v){set(i,j,v);});
    }
    
    template<class VECTOR2>
    SymmMatrixXm(MatrixXm<VECTOR2>& x): 
      MatrixXm<VECTOR>(x.nrows,x.ncols,_Zero()){
      assert(nrows==ncols);
      CONVERT_WARNING("MatrixXm<VECTOR2>","SymmMatrixXm<VECTOR>");
      x.for_each_filled([this](const INDEX i, const INDEX j, const SCALAR v){set(i,j,v);});
    }
    
    SymmMatrixXm(MatrixXm<VECTOR>&& x): 
      MatrixXm<VECTOR>(std::move(x),_NoWarn()){
      MOVECONVERT_WARNING("MatrixXm<VECTOR>","SymmMatrixXm<VECTOR>");
      symmetrize();
    }

    SymmMatrixXm(const Cmatrix& x): 
      MatrixXm<VECTOR>(x,_NoWarn()){
      CONVERT_WARNING("Cmatrix","SymmMatrixXm<VECTOR>");
      if(!x.isSymmetricFormat()) symmetrize();
    }
    

  public: // conversions to other matrix classes --------------------------------------------------------------------


    operator Cmatrix(){
      CONVERT_WARNING("SymmMatrixXm<VECTOR>","Cmatrix");
      Cmatrix M=Cmatrix::Zero(MatrixXm<VECTOR>::nrows,MatrixXm<VECTOR>::ncols);
      for_each_filled([&M](const int i, const int j, const SCALAR v){M(i,j)=v;});
      return M;
    }  

    operator SymmCmatrix(){
      CONVERT_WARNING("SymmMatrixXm<VECTOR>","SymmCmatrix");
      SymmCmatrix M=Cmatrix::Zero(MatrixXm<VECTOR>::nrows,MatrixXm<VECTOR>::ncols);
      for_each_filled([&M](const int i, const int j, const SCALAR v){M(i,j)=v;});
      return M;
    }  


  public: // named constructors -------------------------------------------------------------------------------------


    static SymmMatrixXm Zero(const int _nrows, const int _ncols){
      return SymmMatrixXm<VECTOR>(_nrows,_ncols,_Zero());}

    static SymmMatrixXm Filled(const int _nrows, const int _ncols, const SCALAR v){
      SymmMatrixXm<VECTOR> M(_nrows,_ncols);
      for(int j=0; j<_ncols; j++) (*M.cols[j])=VECTOR::Filled(_nrows,v); 
      return M;
    }

    static SymmMatrixXm Identity(const int n){
      SymmMatrixXm<VECTOR> M(n,n,_Zero()); 
      for(int i=0; i<n; i++) (*M.cols[i])(i)=1;
      return M;
    }

    static SymmMatrixXm Uniform(const int _nrows, const int _ncols){
      SymmMatrixXm<VECTOR> M(_nrows,_ncols);
      uniform_real_distribution<SCALAR> distr;
      for(int i=0; i<_nrows; i++)
	for(int j=0; j<=i; j++)
	  M.set(i,j,distr(randomNumberGenerator));
      return M;
    }

    static SymmMatrixXm Gaussian(const int _nrows, const int _ncols){
      SymmMatrixXm<VECTOR> M(_nrows,_ncols);
      normal_distribution<SCALAR> distr;
      for(int i=0; i<_nrows; i++)
	for(int j=0; j<=i; j++)
	  M.set(i,j,distr(randomNumberGenerator));
      return M;
    }

    static SymmMatrixXm Bernoulli(const int _nrows, const int _ncols, const double p=0.5){
      SymmMatrixXm<VECTOR> M(_nrows,_ncols);
      bernoulli_distribution distr(p);
      for(int i=0; i<_nrows; i++)
	for(int j=0; j<=i; j++)
	  if(distr(randomNumberGenerator)==1) 
	    M.set(i,j,1);
      return M;
    }


  public: // multithreading ---------------------------------------------------------------------------------------


    void flush() const{
      FlushedGuard<typename MatrixXm<VECTOR>::ColumnManager> guard(column_manager);
      const_cast<SymmMatrixXm&>(*this).mirrorFromChanged();
    }

    void mirrorFromChanged(){
    }

    
  public: // element access -----------------------------------------------------------------------------------------


    SCALAR& operator()(const int i, const int j){
      SYMMETRY_UNSAFE("SymmMatrixXm<VECTOR>::operator(...)");
      assert(i<nrows); assert(j<ncols); return (*cols[j])(i);} 

    SCALAR& element(const int i, const int j){
      SYMMETRY_UNSAFE("SymmMatrixXm<VECTOR>::element(...)");
      assert(i<nrows); assert(j<ncols); return (*cols[j])(i);} 

    //void set(const int i, const int j, const SCALAR v){ 
    //  assert(i<nrows); assert(j<ncols); (*cols[j])(i)=v; (*cols[i])(j)=v;} 


  public: // iterators ----------------------------------------------------------------------------------------------


    void for_each_filled(std::function<void(const INDEX, const INDEX, SCALAR&)> lambda){
      SYMMETRY_UNSAFE("SymmMatrixXm<VECTOR>::for_each(...)");
      MatrixXm<VECTOR>::for_each_filled(lambda);
    }

    void for_each_filled_in_row(const int i, std::function<void(const INDEX, SCALAR&)> lambda){
      SYMMETRY_UNSAFE("SymmMatrixXm<VECTOR>::for_each_in_row(...)");
      MatrixXm<VECTOR>::for_each_filled_in_row(i,lambda);
    }

    void for_each_filled_in_column(const int j, std::function<void(const INDEX, SCALAR&)> lambda){
      SYMMETRY_UNSAFE("SymmMatrixXm<VECTOR>::for_each_in_column(...)");
      MatrixXm<VECTOR>::for_each_filled_in_column(j,lambda);
    }


  public: // scalar valued operations -----------------------------------------------------------------------------


    SCALAR inp_of_rows(const int i1, const int i2) const{
      return  MatrixXm<VECTOR>::inp_of_columns(i1,i2);}


  public: // in-place operations ----------------------------------------------------------------------------------


    SymmMatrixX<VECTOR>& operator+=(const SymmMatrixX<VECTOR>& x){assert(ncols==x.ncols);
      pfor_each_column([&x](const int j, const VECTOR& v){v+=*x.cols[j];}); return *this;}
    SymmMatrixX<VECTOR>& operator-=(const SymmMatrixX<VECTOR>& x){assert(ncols==x.ncols);
      pfor_each_column([&x](const int j, const VECTOR& v){v+=*x.cols[j];}); return *this;}


  public: // construct operations ---------------------------------------------------------------------------------

    //SymmMatrixXm<VECTOR>(xgram())



  public: // Givens rotations --------------------------------------------------------------------------------------


    SymmMatrixXm<VECTOR>& conjugate(const GivensRotation& Q){
      column_manager.addJob({Q.i1,Q.i2},
	[this,&Q](){
	  MatrixX<VECTOR>::applyFromRightT(Q);
	  lock_guard<mutex> guard(changed_cols_mx);
	  changed_cols.insert(Q.i1);
	  changed_cols.insert(Q.i2);
	  //addToBuffers(*cols[Q.i1]);
	  //addToBuffers(*cols[Q.i2]);
	});
      //[this,&Q](){MatrixX<VECTOR>::applyFromLeft(Q);}
	
      return *this;
    }

    SymmMatrixXm<VECTOR>& conjugateT(const GivensRotation& Q){
      column_manager.addJob({Q.i1,Q.i2},
	[this,&Q](){
	  MatrixX<VECTOR>::applyFromRight(Q);
	  //addToBuffers(*cols[Q.i1]);
	  //addToBuffers(*cols[Q.i2]);
	});
      //[this,&Q](){MatrixX<VECTOR>::applyFromLeftT(Q);}
      return *this;
    }

  
  public: // KpointOp<k> ---------------------------------------------------------------------------------------------


    template<int k>
    SymmMatrixXm<VECTOR>& conjugate(const KpointOp<k>& Q){
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
    SymmMatrixXm<VECTOR>& conjugateT(const KpointOp<k>& Q){
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


  private: // Disabled methods --------------------------------------------------------------------------------------


    using MatrixXm<VECTOR>::set;

    using MatrixXm<VECTOR>::multiplyRowsBy;
    using MatrixXm<VECTOR>::multiplyColumnsBy;
    using MatrixXm<VECTOR>::divideRowsBy;
    using MatrixXm<VECTOR>::divideColumnsBy;

    using MatrixXm<VECTOR>::applyFromLeft;
    using MatrixXm<VECTOR>::applyFromLeftT;
    using MatrixXm<VECTOR>::applyFromRight;
    using MatrixXm<VECTOR>::applyFromRightT;





  };


  // Downcaster -----------------------------------------------------------------------------------------------------


  template<class VECTOR>
  class AsSymmMatrixXm: public SymmMatrixXm<VECTOR>{
  public:
    AsSymmMatrixXm()=delete;
    AsSymmMatrixXm& operator=(const AsSymmMatrixXm<VECTOR>& x)=delete;
    AsSymmMatrixXm(MatrixXm<VECTOR>& x): SymmMatrixXm<VECTOR>(x.shallow(),_Downcast()){}
    AsSymmMatrixXm(const MatrixXm<VECTOR>& x): SymmMatrixXm<VECTOR>(x.shallow(),_Downcast()){}
    ~AsSymmMatrixXm(){SymmMatrixXm<VECTOR>::cols.clear();}
  };

  template<class VECTOR>
  SymmMatrixXm<VECTOR>::SymmMatrixXm(const AsSymmMatrixXm<VECTOR>& x): 
    MatrixXm<VECTOR>(x,_NoWarn()) {
    COPY_WARNING("SymmMatrixXm<VECTOR>");
  }

  template<class VECTOR>
  AsSymmMatrixXm<VECTOR> as_symmetric(MatrixXm<VECTOR>& x){
    return AsSymmMatrixXm<VECTOR>(x);}

  template<class VECTOR>
  const AsSymmMatrixXm<VECTOR> as_symmetric(const MatrixXm<VECTOR>& x){
    return AsSymmMatrixXm<VECTOR>(x);}

  template<class VECTOR>
  SymmMatrixXm<VECTOR> as_symmetric(MatrixXm<VECTOR>&& x){
    MOVE_WARNING("MatrixXm");
    return SymmMatrixXm<VECTOR>(std::move(x),_Downcast());}
  

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
    template<class VECTOR2>
    operator MatrixX<VECTOR2>(){
      CONVERT_WARNING("SymmMatrixXm<VECTOR>","MatrixX<VECTOR2>");
      return MatrixX<VECTOR2>(MatrixXm<VECTOR>(*this),_NoWarn());
    }  

    template<class VECTOR2>
    operator MatrixX<VECTOR2>(){
      CONVERT_WARNING("SymmMatrixXm<VECTOR>","MatrixX<VECTOR2>");
      return MatrixX<VECTOR2>(MatrixXm<VECTOR>(*this),_NoWarn());
    }  

    template<class VECTOR2>
    operator MatrixXm<VECTOR2>(){
      CONVERT_WARNING("SymmMatrixXm<VECTOR>","MatrixXm<VECTOR2>");
      return MatrixXm<VECTOR2>(MatrixXm<VECTOR>(*this),_NoWarn());
    }  
    */


   /*
    template<class VECTOR2>
    SymmMatrixXm<VECTOR>& multiplyRowsBy(const VECTOR2& v)=0{
      SYMMETRY_UNSAFE("SymmMatrixXm<VECTOR>::multiplyRowsBy(...)");
      MatrixXm<VECTOR>::multiplyRowsBy(v);
      return *this;
    }
    

    template<class VECTOR2>
    SymmMatrixXm<VECTOR>& multiplyColumnsBy(const VECTOR2& v){
      SYMMETRY_UNSAFE("SymmMatrixXm<VECTOR>::multiplyRowsBy(...)");
      MatrixXm<VECTOR>::multiplyColumnsBy(v);
      return *this;
    }

   template<class VECTOR2>
    SymmMatrixXm<VECTOR>& divideRowsBy(const VECTOR2& v){
      SYMMETRY_UNSAFE("SymmMatrixXm<VECTOR>::divideRowsBy(...)");
      MatrixXm<VECTOR>::divideRowsBy(v);
      return *this;
    }

    template<class VECTOR2>
    SymmMatrixXm<VECTOR>& divideColumnsBy(const VECTOR2& v){
      SYMMETRY_UNSAFE("SymmMatrixXm<VECTOR>::divideRowsBy(...)");
      MatrixXm<VECTOR>::divideColumnsBy(v);
      return *this;
    }
    */
    /*
    template<int k>
    SymmMatrixXm<VECTOR>& applyFromLeft(const KpointOp<k>& Q){
      SYMMETRY_FORBIDDEN("applyFromLeft(const KpointOp<k>&)");
      return *this;
    }

    template<int k>
    SymmMatrixXm<VECTOR>& applyFromLeftT(const KpointOp<k>& Q){
      SYMMETRY_FORBIDDEN("applyFromLeftT(const KpointOp<k>&)");
      return *this;
    }

    template<int k>
    SymmMatrixXm<VECTOR>& applyFromRightT(const KpointOp<k>& Q){
      SYMMETRY_FORBIDDEN("applyFromRightT(const KpointOp<k>&)");
      return *this;
    }

    template<int k>
    SymmMatrixXm<VECTOR>& applyFromRight(const KpointOp<k>& Q){
      SYMMETRY_FORBIDDEN("applyFromRight(const KpointOp<k>&)");
      return *this;
    }
    */

  //template<>
  //xcgram< SymmMatrixXm<Vectorh> >::operator xcgram< SymmMatrixX<Vectorh> >() const{
  //  return xcgram< SymmMatrixX<Vectorh> >(SymmMatrixX<Vectorh>(MatrixX<Vectorh>(A.shallow())));}
    /* overrides other conversions 
       template<class MATRIX>
       SymmMatrixXm(const MATRIX& x): 
       MatrixXm<VECTOR>(x){
       CONVERT_WARNING("MATRIX","SymmMatrixXm<VECTOR>");
       if(!x.isSymmetricFormat()) symmetrize();
       }
    */


    //class ColumnManager: public ResourceManager{
    //public:
    //  ColumnManager(const int ncols){}
    //};
    //int maxthreads=1000;
    //int nthreads=4;
    //int nseg=2;
    //ColumnManager column_manager;
    //mutex change_mx;


