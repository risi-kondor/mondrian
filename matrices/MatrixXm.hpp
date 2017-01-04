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


#ifndef _MatrixXm
#define _MatrixXm

#include "MatrixX.hpp"
#include "ThreadBank.hpp"
#include "ResourceManager.hpp"
#include "VectorBuffer.hpp"

extern default_random_engine randomNumberGenerator;


namespace Mondrian{


  template<class VECTOR>
  class AsMatrixXm;


  template<class VECTOR>
  class MatrixXm: public MatrixX<VECTOR>{
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


  public:

    class ColumnManager: public ResourceManager{
    public:
      ColumnManager(const int ncols){}
    };


    int maxthreads=1000;
    int nthreads=4;
    int nseg=2;
    mutable ColumnManager column_manager;
    // mutex change_mx;

    VectorBuffers buffers; 

    using MatrixX<VECTOR>::nrows;
    using MatrixX<VECTOR>::ncols;
    using MatrixX<VECTOR>::cols;

    bool isMultithreaded() const {return true;}


  public: // constructors --------------------------------------------------------------------------------------------


    MatrixXm(): MatrixXm(0,0){};
    
    MatrixXm(const int _nrows, const int _ncols):
      MatrixX<VECTOR>(_nrows,_ncols), column_manager(_ncols){}
    
    MatrixXm(const int _nrows, const int _ncols, const _Zero dummy): 
      MatrixX<VECTOR>(_nrows,_ncols,dummy), column_manager(_ncols){}
    
    MatrixXm(const int _nrows, const int _ncols, const initializer_list<iivtriple>& list): 
      MatrixX<VECTOR>(_nrows,_ncols,list), column_manager(_ncols){}


  public: // copying -------------------------------------------------------------------------------------------------

    // TODO: flush before copy

    MatrixXm(const MatrixXm<VECTOR>& x): 
      MatrixX<VECTOR>(x,_NoWarn()),
      column_manager(x.ncols){
      COPY_WARNING("MatrixXm<VECTOR>");
    }

    MatrixXm(const MatrixXm<VECTOR>& x, const _NoWarn dummy): 
      MatrixX<VECTOR>(x,_NoWarn()), 
      column_manager(x.ncols){
    }

    MatrixXm(MatrixXm<VECTOR>&& x): 
      MatrixX<VECTOR>(std::move(x)), 
      column_manager(x.ncols){
      MOVE_WARNING("MatrixXm<VECTOR>");
    }

    MatrixXm(MatrixXm<VECTOR>&& x, const _NoWarn dummy): 
      MatrixX<VECTOR>(std::move(x)), 
      column_manager(x.ncols){
    }

    MatrixXm<VECTOR>& operator=(const MatrixXm<VECTOR>& x){
      MatrixX<VECTOR>::operator=(x);
      return *this;
    }

    MatrixXm<VECTOR>& operator=(MatrixXm<VECTOR>&& x){
      MatrixX<VECTOR>::operator=(std::move(x));
      return *this;
    }

    MatrixXm<VECTOR> copy(){
      FlushedGuard<ColumnManager> guard(column_manager);
      return MatrixX<VECTOR>::copy();
    }

    void assign(const MatrixXm<VECTOR>& x){
      FlushedGuard<ColumnManager> guard(column_manager);
      MatrixX<VECTOR>::assign(x);}

    void shallow(){
      FlushedGuard<ColumnManager> guard(column_manager);    
      return MatrixX<VECTOR>::shallow();}


  public: // downcasting --------------------------------------------------------------------------------------------
    

    MatrixXm(const MatrixX<VECTOR>& x): 
      MatrixX<VECTOR>(x,_NoWarn()), 
      column_manager(x.ncols){
      DOWNCASTCOPY_WARNING("MatrixX<VECTOR>","MatrixXm<VECTOR>");
    }

    MatrixXm(MatrixX<VECTOR>&& x): 
      MatrixX<VECTOR>(std::move(x),_NoWarn()), 
      column_manager(x.ncols){
      DOWNCAST_WARNING("MatrixX<VECTOR>","MatrixXm<VECTOR>");
    }


  public: // conversions --------------------------------------------------------------------------------------------


    template<class VECTOR2>
    MatrixXm(const MatrixX<VECTOR2>& x): 
      MatrixX<VECTOR>(x,_NoWarn()), 
      column_manager(x.ncols){
      CONVERT_WARNING("MatrixX<VECTOR>","MatrixXm<VECTOR>");
    }

    MatrixXm(const Cmatrix& x):
      MatrixX<VECTOR>(x,_NoWarn()),
      column_manager(x.ncols){
      CONVERT_WARNING("Cmatrix","MatrixXm<VECTOR>");
    }

    MatrixXm(const Cmatrix& x, const _NoWarn dummy):
      MatrixX<VECTOR>(x,_NoWarn()),
      column_manager(x.ncols){}

    operator Cmatrix(){
      CONVERT_WARNING("MatrixXm<VECTOR>","Cmatrix");
      //ResourceManager::FlushedGuard guard(const_cast<ColumnManager&>(column_manager));
      FlushedGuard<ColumnManager> guard(column_manager);    
      return static_cast< MatrixX<VECTOR> >(*this);
    }


  public: // comparisons --------------------------------------------------------------------------------------------


    template<class VECTOR2>
    bool operator==(const MatrixX<VECTOR2>& X) const{ 
      FlushedGuard<ColumnManager> guard(column_manager);    
      return MatrixX<VECTOR>::operator==(X);
    }
    

  public: // multithreading -----------------------------------------------------------------------------------------


    void flush() const{
      FlushedGuard<ColumnManager> guard(column_manager);    
    }


  public: // element access -----------------------------------------------------------------------------------------
  

    SCALAR read(const int i, const int j) const {assert(j<ncols);
      return column_manager.template addTask<SCALAR>(j,[this,i,j](){return cols[j]->read(i);});}

    SCALAR operator()(const int i, const int j) const {return read(i,j);}
  
    //SCALAR& operator()(const int i, const int j){ // should be disabled 
    //  assert(i<nrows); assert(j<ncols); return (*cols[j])(i);}

    void set(const int i, const int j, const SCALAR v) {assert(j<ncols); 
      return column_manager.addJob(j,[this,i,j,v](){cols[j]->set(i,v);});}

    bool isFilled(const int i, const int j) const {assert(j<ncols);
      return column_manager.template addTask<bool>(j,[this,i,j](){return cols[j]->isFilled(i);});}

    int nFilled() const{
      FlushedGuard<ColumnManager> guard(column_manager);
      int t=0; for(int j=0; j<ncols; j++) t+=cols[j]->nFilled(); return t;}

    template<class VECTOR2>
    VECTOR2 column(const int j) const {assert(j<ncols);
      return column_manager.template addTask<VECTOR2>(j,[this,j](){return *cols[j];});}

    template<class VECTOR2>
    VECTOR row(const int i) const {assert(i<nrows);
      FlushedGuard<ColumnManager> guard(column_manager);
      VECTOR2 x(ncols); for(int j=0; j<ncols; j++) x(j)=(*cols[j])(i); return x;}

    template<class VECTOR2>
    VECTOR diag() const {int t=min(nrows,ncols);
      FlushedGuard<ColumnManager> guard(column_manager);
      VECTOR2 x(t); for(int i=0; i<t; i++) x(i)=(*cols(i))(i); return x;}

    
  public: // function mappings --------------------------------------------------------------------------------------


    void pfor_each_column(std::function<void(const INDEX, const VECTOR&)> lambda) const{
      FlushedGuard<ColumnManager> guard(column_manager);
      int nblocks=std::max(std::min(nthreads,ncols/10),1);
      const int w=ncols/nblocks;
      ThreadBank threads(maxthreads);
      for(int J=0; J<nblocks; J++){
	int jmax=(J==nblocks-1) ? (ncols):((J+1)*w);
	threads.add([this,lambda](const int jmin, const int jmax){
	    for(int j=jmin; j<jmax; j++) 
	      lambda(j,*cols[j]);
	  },J*w,jmax);
      }
    }

    void pfor_each_column(std::function<void(const INDEX, VECTOR&)> lambda){
      //ResourceManager::FlushedGuard guard(const_cast<ColumnManager&>(column_manager));
      FlushedGuard<ColumnManager> guard(column_manager);
      int nblocks=std::max(std::min(nthreads,ncols/10),1);
      const int w=ncols/nblocks;
      ThreadBank threads(maxthreads);
      for(int J=0; J<nblocks; J++){
	int jmax=(J==nblocks-1) ? (ncols):((J+1)*w);
	threads.add([this,lambda](const int jmin, const int jmax){
	    for(int j=jmin; j<jmax; j++) 
	      lambda(j,*cols[j]);
	  },J*w,jmax);
      }
    }

    void pfor_each_filled_in_row(const int i, std::function<void(const INDEX, const SCALAR&)> lambda) const{
      pfor_each_column([lambda,i](const int j, const VECTOR& v){
	  const SCALAR* p=v.ptr_if_filled(i); 
	  if(p!=nullptr) lambda(j,*p);
	});
    }

    void pfor_each_filled_in_row(const int i, std::function<void(const INDEX, SCALAR&)> lambda){
      pfor_each_column([lambda,i](const int j, const VECTOR& v){
	  SCALAR* p=v.ptr_if_filled(i); 
	  if(p!=nullptr) lambda(j,*p);
	});
    }

    void pfor_each_filled_in_column(const int j, std::function<void(const INDEX, const SCALAR&)> lambda) const{
      const_cast<ColumnManager&>(column_manager).addJob(j,[this,j,lambda](){
	  static_cast<const VECTOR*>(cols[j])->for_each_filled(lambda);});
    }

    void pfor_each_filled_in_column(const int j, std::function<void(const INDEX, SCALAR&)> lambda){
      column_manager.addJob(j,[this,j,lambda](){cols[j]->for_each_filled(lambda);});
    }

    void pfor_each_filled(std::function<void(const INDEX, const INDEX, const SCALAR&)> lambda) const{
      pfor_each_column([lambda](const int j, const VECTOR& v){
	  v.for_each_filled([lambda,j](const int i, const SCALAR v){
	      lambda(i,j,v);
	    });
	});
    }
    
    void pfor_each_filled(std::function<void(const INDEX, const INDEX, SCALAR)> lambda){
      pfor_each_column([lambda](const int j, const VECTOR& v){
	  v.for_each_filled([lambda,j](const int i, const SCALAR v){
	      lambda(i,j,v);
	    });
	});
    }

    
    template<class TYPE>
    TYPE accumulate_over_columns(std::function<TYPE(const VECTOR&)> lambda, 
      std::function<TYPE(const TYPE&, const TYPE&)> accumulator=std::plus<TYPE>(),
      const TYPE& t0=TYPE()) const{

      FlushedGuard<ColumnManager> guard(column_manager);
      int nblocks=std::max(std::min(nthreads,ncols/10),1);
      TYPE partial[nblocks];
      const int w=ncols/nblocks;
      {ThreadBank threads(maxthreads);
	for(int J=0; J<nblocks; J++){
	  int jmax=(J==nblocks-1) ? (ncols):((J+1)*w);
	  threads.add([this,lambda,accumulator,t0](const int jmin, const int jmax, TYPE*& target){
	      TYPE t(t0);
	      for(int j=jmin; j<jmax; j++) 
		t=accumulator(t,lambda(*cols[j]));
	      *target=t;
	    },
	    J*w,jmax,&partial[J]);
	}
      }

      TYPE result(t0);
      for(int J=0; J<nblocks; J++)
	result=accumulator(result,partial[J]);
      return result;
    }


    iipair find_best_over_columns(std::function<ivpair(const VECTOR&)> lambda, 
      std::function<bool(const SCALAR&, const SCALAR&)> selector) const{

      const int nblocks=::max(::min(nthreads,ncols/10),1);
      const int w=ncols/nblocks;
      SCALAR bestv[nblocks];
      INDEX  besti[nblocks];
      INDEX  bestj[nblocks];
      {ThreadBank threads(maxthreads);
	for(int J=0; J<nblocks; J++){
	  int jmax=(J==nblocks-1) ? (ncols):((J+1)*w);
	  threads.add([this,lambda,selector](const int jmin, const int jmax, 
	      SCALAR* vtarget, INDEX* itarget, INDEX* jtarget){
	      ivpair p=lambda(*cols[jmin]);
	      INDEX bestj=jmin;
	      INDEX besti=p.first;
	      SCALAR bestv=p.second;
	      for(int j=jmin+1; j<jmax; j++){
		ivpair p=lambda(*cols[j]);	      
		if(selector(bestv,p.second)){bestj=j; besti=p.first; bestv=p.second;}
	      }
	      *vtarget=bestv;
	      *itarget=besti;
	      *jtarget=bestj;
	    },J*w,jmax,&bestv[J],&besti[J],&bestj[J]);
	}
      }

      INDEX _besti=besti[0];
      INDEX _bestj=bestj[0];
      SCALAR _bestv=bestv[0];
      for(int I=1; I<nblocks; I++)
	if(selector(_bestv,bestv[I])){_besti=besti[I]; _bestj=bestj[I]; _bestv=bestv[I];}
      return iipair(_besti,_bestj);
    }
       
  
  public: // scalar valued operations -------------------------------------------------------------------------------


    SCALAR sum() const{
      return accumulate_over_columns<int>([](const VECTOR& v){return v.sum();},std::plus<SCALAR>());}
    SCALAR norm1() const{
      return accumulate_over_columns<int>([](const VECTOR& v){return v.norm1();},std::plus<SCALAR>());}
    SCALAR norm2() const{
      return accumulate_over_columns<int>([](const VECTOR& v){return v.norm2();},std::plus<SCALAR>());}
    int nnz() const{
      return accumulate_over_columns<int>([](const VECTOR& v){return v.nnz();},std::plus<int>());}

    SCALAR inp_of_columns(const int j1, const int j2) const{
      return column_manager.addTask({j1,j2},[this,j1,j2](){return MatrixX<VECTOR>::inp_of_columns(j1,j2);});}

    SCALAR inp_of_rows(const int i1, const int i2) const{
      FlushedGuard<ColumnManager> guard(column_manager);
      int nblocks=std::max(std::min(nthreads,ncols/10),1);
      SCALAR accum[nblocks];
      const int w=ncols/nblocks;
      ThreadBank threads(maxthreads);
      for(int J=0; J<nblocks; J++){
	int jmax=(J==nblocks-1) ? (ncols):((J+1)*w);
	threads.add([this,i1,i2](const int jmin, const int jmax, SCALAR* target){
	    int t=0; 
	    for(int j=jmin; j<jmax; j++) 
	      t+=(cols[j]->read(i1))*(cols[j]->read(i2));
	    *target=t;
	  },J*w,jmax,&accum[J]);
      }
      SCALAR t=0; for(int J=0; J<nblocks; J++) t+=accum[J]; 
      return t;
    }


  public: // in-place operations ----------------------------------------------------------------------------------


    MatrixX<VECTOR>& operator+=(const MatrixX<VECTOR>& x){assert(ncols==x.ncols);
      pfor_each_column([&x](const int j, const VECTOR& v){v+=*x.cols[j];}); return *this;}
    MatrixX<VECTOR>& operator-=(const MatrixX<VECTOR>& x){assert(ncols==x.ncols);
      pfor_each_column([&x](const int j, const VECTOR& v){v+=*x.cols[j];}); return *this;}

    template<class VECTOR2>
    MatrixX<VECTOR>& mulitplyRowsBy(const VECTOR2& x){
      pfor_each_column([&x](const int j, const VECTOR& v){v*=x;}); return *this;}
    template<class VECTOR2>
    MatrixX<VECTOR>& divideRowsBy(const VECTOR2& x){
      pfor_each_column([&x](const int j, const VECTOR& v){v/=x;}); return *this;}
	
    template<class VECTOR2>
    MatrixX<VECTOR>& mulitplyColumnsBy(const VECTOR2& x){
      pfor_each_column([&x](const int j, const VECTOR& v){v*=x(j);}); return *this;}
    template<class VECTOR2>
    MatrixX<VECTOR>& divideColumnsBy(const VECTOR2& x){
      pfor_each_column([&x](const int j, const VECTOR& v){v/=x(j);}); return *this;}
	

  public: // construct operations -----------------------------------------------------------------------------------

    
    MatrixXm(const xcgram<Cmatrix>& e): 
      MatrixX<VECTOR>(e.A.ncols,e.A.ncols,_Zero()), 
      column_manager(e.A.ncols){
      const Cmatrix& M=e.A;
      const int n=M.ncols;
      const int m=M.nrows;
      const int w=n/nseg;
      ThreadBank threads(maxthreads);
      for(int I=0; I<nseg; I++){
	int imax=(I==nseg-1) ? (n):((I+1)*w);
	for(int J=0; J<=I-1; J++)
	  threads.add([this,&M,m](const int imin, const int imax, const int jmin, const int jmax){
	      for(int i=imin; i<imax; i++)
		for(int j=jmin; j<jmax; j++){
		  SCALAR t=0;
		  for(int k=0; k<m; k++)
		    t+=M.array[i*m+k]*M.array[j*m+k];
		  if(t!=0) {cols[j]->set_msafe(i,t); cols[i]->set_msafe(j,t);}
		  //if(t!=0) {set(i,j,t); set(j,i,t);}
		}
	    },I*w,imax,J*w,(J+1)*w);
	threads.add([this,&M,m](const int imin, const int imax){
	    for(int i=imin; i<imax; i++)
	      for(int j=imin; j<=i; j++){
		SCALAR t=0;
		for(int k=0; k<m; k++)
		  t+=M.array[i*m+k]*M.array[j*m+k];
		//if(t!=0) {set(i,j,t); set(j,i,t);}
		if(t!=0) {cols[j]->set_msafe(i,t); cols[i]->set_msafe(j,t);}
	      }
	  },I*w,imax);	
      }
    }
    

    template<class VECTOR2>
    MatrixXm(const xcgram<MatrixX<VECTOR2> >& e): 
      MatrixX<VECTOR>(e.A.ncols,e.A.ncols,_Zero()), 
      column_manager(e.A.ncols){
      const MatrixX<VECTOR2>& M=e.A;
      const int n=M.ncols;
      const int m=M.nrows;
      const int w=n/nseg;
      ThreadBank threads(maxthreads);
      for(int I=0; I<nseg; I++){
	int imax=(I==nseg-1) ? (n):((I+1)*w);
	for(int J=0; J<=I-1; J++)
	  threads.add([this,&M,m](const int imin, const int imax, const int jmin, const int jmax){
	      for(int i=imin; i<imax; i++)
		for(int j=jmin; j<jmax; j++){
		  SCALAR t=M.cols[i]->dot(*M.cols[j]);
		  if(t!=0) {cols[j]->set_msafe(i,t); cols[i]->set_msafe(j,t);}
		}
	    },I*w,imax,J*w,(J+1)*w);
	threads.add([this,&M,m](const int imin, const int imax){
	    for(int i=imin; i<imax; i++)
	      for(int j=imin; j<=i; j++){
		SCALAR t=M.cols[i]->dot(*M.cols[j]);
		if(t!=0) {cols[j]->set_msafe(i,t); cols[i]->set_msafe(j,t);}
	      }
	  },I*w,imax);	
      }
    }


    template<class VECTOR2>
    MatrixXm(const xcgram<SymmMatrixX<VECTOR2> >& e): 
      MatrixX<VECTOR>(e.A.ncols,e.A.ncols,_Zero()),
      column_manager(e.A.ncols){
      const MatrixX<VECTOR2>& M=e.A;
      const int n=M.ncols;
      const int m=M.nrows;
      const int w=n/nseg;
      ThreadBank threads(maxthreads);
      for(int I=0; I<nseg; I++){
	int imax=(I==nseg-1) ? (n):((I+1)*w);
	threads.add([this,&M](const int imin, const int imax){
	    for(int i=imin; i<imax; i++){
	      VectorSpattern s;
	      M.cols[i]->for_each_filled([&M,&s,i](const int u, const SCALAR x){s.add(*M.cols[u],i);});
	      s.for_each_filled([this,&M,i](const int j){
		  SCALAR t=M.cols[i]->dot(*M.cols[j]);
		  if(t!=0) {cols[j]->set_msafe(i,t); cols[i]->set_msafe(j,t);}
		});
	    }
	  },I*w,imax);
      }
    }

  
    template<class VECTOR2>
    MatrixXm(const xcgram<SymmMatrixXm<VECTOR2> >& e):
      MatrixXm<VECTOR>(gram(as_symmetric(MatrixX<VECTOR>(e.A)))){}


  public: // Givens rotations ---------------------------------------------------------------------------------------


    MatrixXm<VECTOR>& applyFromLeft(const GivensRotation& Q){
      FlushedGuard<ColumnManager> guard(column_manager);
      int nblocks=std::max(std::min(nthreads,ncols/10),1);
      const int w=ncols/nblocks;
      ThreadBank threads(maxthreads);
      for(int I=0; I<nblocks; I++){
	int imax=(I==nblocks-1) ? (ncols):((I+1)*w);
	threads.add([this,&Q](const int imin, const int imax){
		      for(int i=imin; i<imax; i++) {cols[i]->apply(Q);}},I*w,imax);
      }
      return *this;
    }

    MatrixXm<VECTOR>& applyFromLeftT(const GivensRotation& Q){
      FlushedGuard<ColumnManager> guard(column_manager);
      int nblocks=std::max(std::min(nthreads,ncols/10),1);
      const int w=ncols/nblocks;
      ThreadBank threads(maxthreads);
      for(int I=0; I<nblocks; I++){
	int imax=(I==nblocks-1) ? (ncols):((I+1)*w);
	threads.add([this,&Q](const int imin, const int imax){
		      for(int i=imin; i<imax; i++) {cols[i]->applyT(Q);}},I*w,imax);
      }
      return *this;
    }

    MatrixXm<VECTOR>& applyFromRightT(const GivensRotation& Q){
      column_manager.addJob({Q.i1,Q.i2},[this,&Q](){MatrixX<VECTOR>::applyFromRightT(Q);});
      return *this;
    }

    MatrixXm<VECTOR>& applyFromRight(const GivensRotation& Q){
      column_manager.addJob({Q.i1,Q.i2},[this,&Q](){MatrixX<VECTOR>::applyFromRight(Q);});
      return *this;
    }


  public: // KpointOp<k> ---------------------------------------------------------------------------------------------


    template<int k>
    MatrixXm<VECTOR>& applyFromLeft(const KpointOp<k>& Q){
      FlushedGuard<ColumnManager> guard(column_manager);
      int nblocks=std::max(std::min(nthreads,ncols/10),1);
      const int w=ncols/nblocks;
      ThreadBank threads(maxthreads);
      for(int I=0; I<nblocks; I++){
	int imax=(I==nblocks-1) ? (ncols):((I+1)*w);
	threads.add([this,&Q](const int imin, const int imax){
	    for(int i=imin; i<imax; i++) {cols[i]->apply(Q);}},I*w,imax);
      }
      return *this;
    }

    template<int k>
    MatrixXm<VECTOR>& applyFromLeftT(const KpointOp<k>& Q){
      FlushedGuard<ColumnManager> guard(column_manager);
      int nblocks=std::max(std::min(nthreads,ncols/10),1);
      const int w=ncols/nblocks;
      ThreadBank threads(maxthreads);
      for(int I=0; I<nblocks; I++){
	int imax=(I==nblocks-1) ? (ncols):((I+1)*w);
	threads.add([this,&Q](const int imin, const int imax){
	    for(int i=imin; i<imax; i++) {cols[i]->applyT(Q);}},I*w,imax);
      }
      return *this;
    }

    template<int k>
    MatrixXm<VECTOR>& applyFromRightT(const KpointOp<k>& Q){
      column_manager.addJob(Q.map,[this,&Q](){MatrixX<VECTOR>::applyFromRightT(Q);});
      return *this;
    }

    template<int k>
    MatrixXm<VECTOR>& applyFromRight(const KpointOp<k>& Q){
      column_manager.addJob(Q.map,[this,&Q](){MatrixX<VECTOR>::applyFromRight(Q);});
      return *this;
    }


    // I/O ----------------------------------------------------------------------------------------------------------


    string str() const{
      FlushedGuard<ColumnManager> guard(column_manager);
      return MatrixX<VECTOR>::str();
    }



  };



  // AsMatrixXm -------------------------------------------------------------------------------------------------------

  
  template<class VECTOR>
  class AsMatrixXm: public MatrixXm<VECTOR>{
  public:
    AsMatrixXm()=delete;
    AsMatrixXm& operator=(const AsMatrixXm<VECTOR>& x)=delete;
    AsMatrixXm(MatrixX<VECTOR>& x): MatrixXm<VECTOR>(x.shallow()){}
    AsMatrixXm(const MatrixX<VECTOR>& x): MatrixXm<VECTOR>(x.shallow()){}
    ~AsMatrixXm(){MatrixXm<VECTOR>::cols.clear();}
  };

  template<class VECTOR>
  AsMatrixXm<VECTOR> as_multithreaded(MatrixX<VECTOR>& x){
    return AsMatrixXm<VECTOR>(x);}

  template<class VECTOR>
  const AsMatrixXm<VECTOR> as_multithreaded(const MatrixX<VECTOR>& x){
    return AsMatrixXm<VECTOR>(x);}

  template<class VECTOR>
  MatrixXm<VECTOR> as_multithreaded(MatrixX<VECTOR>&& x){
    return MatrixXm<VECTOR>(std::move(x));}



}

#endif 



    //class ColumnBuffer: public vector<ivpair>{
    //public:
      
    //  mutex access_mx;

    //};


    //vector<ColumnBuffer> cbuffers;
/*
  template<class VECTOR>
  MatrixXm<VECTOR>::MatrixXm(const AsMatrixXm<VECTOR>& x): 
    MatrixXm<VECTOR>(x,_NoWarn()) {
    COPY_WARNING("MatrixXm<VECTOR>");
  }
*/

	//access_mx.resize(ncols); 
	//for(int i=0; i<ncols; i++) access_mx[i]=new mutex();
      //vector<mutex*> access_mx;
      //~ColumnManager(){for(auto p: access_mx) delete p;}
    /*
    MatrixXm(const MatrixX<VECTOR>& x, const _Downcast dummy): 
      MatrixX<VECTOR>(x,_NoWarn()), 
      column_manager(x.ncols){
      DOWNCASTCOPY_WARNING("MatrixX<VECTOR>","MatrixXm<VECTOR>");
    }

    MatrixXm(MatrixX<VECTOR>&& x, const _Downcast dummy): 
      MatrixX<VECTOR>(std::move(x),_NoWarn()), 
      column_manager(x.ncols){
      DOWNCAST_WARNING("MatrixX<VECTOR>","MatrixXm<VECTOR>");
    }
    */

    /*
    template<class VECTOR2>
    MatrixXm(const MatrixX<VECTOR2>& x, const _Downcast dummy): 
      MatrixX<VECTOR>(x,_NoWarn()), 
      column_manager(x.ncols){
      DOWNCASTCOPY_WARNING("MatrixX<VECTOR>","MatrixXm<VECTOR>");
    }

    template<class VECTOR2>
    MatrixXm(MatrixX<VECTOR2>&& x, const _Downcast dummy): 
      MatrixX<VECTOR>(std::move(x),_NoWarn()), 
      column_manager(x.ncols){
      DOWNCAST_WARNING("MatrixX<VECTOR>","MatrixXm<VECTOR>");
    }
    */

    /*
      MatrixX<VECTOR>(e.A.ncols,e.A.ncols,_Zero()), 
      column_manager(e.A.ncols){
      const MatrixX<VECTOR2>& M=e.A;
      const int n=M.ncols;
      const int m=M.nrows;
      const int w=n/nseg;
      ThreadBank threads(maxthreads);
      for(int I=0; I<nseg; I++){
	int imax=(I==nseg-1) ? (n):((I+1)*w);
	threads.add([this,&M](const int imin, const int imax){
		      for(int i=imin; i<imax; i++){
			VectorSpattern s;
			M.cols[i]->for_each_filled([&M,&s,i](const int u, const SCALAR x){s.add(*M.cols[u],i);});
			s.for_each_filled([this,&M,i](const int j){
					    SCALAR t=M.cols[i]->dot(*M.cols[j]);
					    if(t!=0) {cols[j]->set_msafe(i,t); cols[i]->set_msafe(j,t);}
					  });
		      }
		    },I*w,imax);
      }
    }
    */
    //void for_each_column(std::function<void(const INDEX, const VECTOR&)> lambda) const{
    //  ResourceManager::FlushedGuard guard(const_cast<ColumnManager&>(column_manager));
    //  for(int j=0; j<ncols; j++) lambda(j,*cols[j]);
    //}

    /*
    void for_each_filled_in_row(const int i, std::function<void(const INDEX, SCALAR)> lambda) const{
      ResourceManager::FlushedGuard guard(const_cast<ColumnManager&>(column_manager));
      for(int j=0; j<ncols; j++){
	SCALAR* p=cols[j]->ptr_if_filled(i); 
	if(p!=nullptr) lambda(j,*p);
      }
    }
    */
