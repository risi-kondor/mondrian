/* -----------------------------------------------------------------------------
 
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
 
----------------------------------------------------------------------------- */


#ifndef _BlockedMatrix
#define _BlockedMatrix

#include "Matrix.hpp"
#include "Detachable.hpp"
#include "BlockStructure.hpp"
#include "BlockedVector.hpp"


namespace Mondrian{

  template<class MATRIX>
  class BlockedMatrix: public Matrix, public Detachable{
  public:
    
    BlockedMatrix(const BlockedMatrix<MATRIX>& x);
    BlockedMatrix(BlockedMatrix<MATRIX>&& x);
    BlockedMatrix& operator=(const BlockedMatrix<MATRIX>& x);
    BlockedMatrix& operator=(BlockedMatrix<MATRIX>&& x);
    BlockedMatrix<MATRIX> copy() const;
    BlockedMatrix<MATRIX> shallow() const;
    ~BlockedMatrix(); 
    void detach(){nb=0; mb=0; blocks=nullptr;}

    BlockedMatrix(): Matrix(0,0), nb(0),mb(0){}
    BlockedMatrix(const int _nb, const int _mb);
    BlockedMatrix(const int _nb, const int _mb, const int _nrows, const int _ncols);
    BlockedMatrix(const BlockStructure& B1, const BlockStructure& B2);  

    int nb;
    int mb;
    MATRIX** blocks=nullptr;


  public: // named constructors
    
    static BlockedMatrix<MATRIX> Zero(const BlockStructure& B1, const BlockStructure& B2);
    static BlockedMatrix<MATRIX> Filled(const BlockStructure& B1, const BlockStructure& B2, const SCALAR t);
    static BlockedMatrix<MATRIX> Identity(const BlockStructure& B1);

    static BlockedMatrix<MATRIX> Uniform(const BlockStructure& B1, const BlockStructure& B2);
    static BlockedMatrix<MATRIX> Gaussian(const BlockStructure& B1, const BlockStructure& B2);
    static BlockedMatrix<MATRIX> Bernoulli(const BlockStructure& B1, const BlockStructure& B2, const double p=0.5);


  public: // attributes

    bool isDense() const {if(nb*mb>0) return blocks[0]->isDense(); else return true;}
    bool isSparse() const {if(nb*mb>0) return blocks[0]->isSparse(); else return false;}


  public: // block access 
    
    MATRIX& block(const int i, const int j){return *blocks[i+nb*j];}
    const MATRIX& block(const int i, const int j) const {return *blocks[i+nb*j];}

    void for_each_block(std::function<void(const INDEX, const INDEX, MATRIX&)> lambda){
      for(int i=0; i<nb; i++) for(int j=0; j<mb; j++) lambda(i,j,*blocks[j*nb+i]);}
    void for_each_block(std::function<void(const INDEX, const INDEX, const MATRIX&)> lambda) const{
      for(int i=0; i<nb; i++) for(int j=0; j<mb; j++) lambda(i,j,*blocks[j*nb+i]);}
    
    MATRIX*& block_ptr(const int i, const int j) {return blocks[i+nb*j];}
    MATRIX* const& block_ptr(const int i, const int j) const {return blocks[i+nb*j];}

    
  public: // element access

    SCALAR& operator()(const int i, const int j) {
      assert(i<nrows); int t1=0; int I=0; for(; I<nb && t1+blocks[I]->nrows<=i; I++) t1+=blocks[I]->nrows;  
      assert(j<ncols); int t2=0; int J=0; for(; J<mb && t2+blocks[J*nb]->ncols<=j; J++) t2+=blocks[J*nb]->ncols;  
      return (*blocks[I+J*nb])(i-t1,j-t2);
    }
    SCALAR operator()(const int i, const int j) const {
      assert(i<nrows); int t1=0; int I=0; for(; I<nb && t1+blocks[I]->nrows<=i; I++) t1+=blocks[I]->nrows;  
      assert(j<ncols); int t2=0; int J=0; for(; J<mb && t2+blocks[J*nb]->ncols<=j; J++) t2+=blocks[J*nb]->ncols;
      return (*blocks[I+J*nb])(i-t1,j-t2);
    }
    SCALAR read(const int i, const int j) const {
      assert(i<nrows); int t1=0; int I=0; for(; I<nb && t1+blocks[I]->nrows<=i; I++) t1+=blocks[I]->nrows;  
      assert(j<ncols); int t2=0; int J=0; for(; J<mb && t2+blocks[J*nb]->ncols<=j; J++) t2+=blocks[J*nb]->ncols;  
      return (*blocks[I+J*nb])(i-t1,j-t2);
    }
    
    SCALAR& operator()(const iipair& I, const iipair& J){
      return (*blocks[I.first+J.first*nb])(I.second,J.second);}
    SCALAR operator()(const iipair& I, const iipair& J) const {
      return (*blocks[I.first+J.first*nb])(I.second,J.second);}
    SCALAR read(const iipair& I, const iipair& J) const {
      return (*blocks[I.first+J.first*nb])(I.second,J.second);}

    void for_each(std::function<void(const INDEX, const INDEX, SCALAR&)> lambda){
      int joffset=0; for(int j=0; j<mb; joffset+=blocks[j*nb]->ncols, j++){
	int ioffset=0; for(int i=0; i<nb; ioffset+=blocks[i]->nrows, i++)
	  blocks[j*nb+i]->for_each([ioffset,joffset,&lambda](const int i, const int j, SCALAR& v){
					  lambda(ioffset+i,joffset+j,v);});
      }
    }
    void for_each(std::function<void(const INDEX, const INDEX, const SCALAR&)> lambda) const{
      int joffset=0; for(int j=0; j<mb; joffset+=blocks[j*nb]->ncols, j++){
	int ioffset=0; for(int i=0; i<nb; ioffset+=blocks[i]->nrows, i++)
	  blocks[j*nb+i]->for_each([ioffset,joffset,&lambda](const int i, const int j, const SCALAR v){
					  lambda(ioffset+i,joffset+j,v);});
      }
    }
    void for_each_in_row(const int i, std::function<void(const INDEX, SCALAR&)> lambda){
      iipair ixpair=rowix(i); int I=ixpair.first; int subi=ixpair.second; int joffset=0; 
      for(int J=0; J<mb; joffset+=blocks[J*nb]->ncols, J++)
	blocks[J*nb+I]->for_each_in_row(subi,[joffset,&lambda](const int j, SCALAR& v){lambda(joffset+j,v);});
    }	
    void for_each_in_row(const int i, std::function<void(const INDEX, const SCALAR)> lambda){
      iipair ixpair=rowix(i); int I=ixpair.first; int subi=ixpair.second; int joffset=0; 
      for(int J=0; J<mb; joffset+=blocks[J*nb]->ncols, J++)
	blocks[J*nb+I]->for_each_in_row(subi,[joffset,&lambda](const int j, const SCALAR v){lambda(joffset+j,v);});
    }	
    void for_each_in_column(const int j, std::function<void(const INDEX, SCALAR&)> lambda){
      iipair ixpair=colix(j); int J=ixpair.first; int subj=ixpair.second; int ioffset=0; 
      for(int I=0; I<nb; ioffset+=blocks[I]->nrows, I++)
	blocks[J*nb+I]->for_each_in_column(subj,[ioffset,&lambda](const int i, SCALAR& v){lambda(ioffset+i,v);});
    }	
    void for_each_in_column(const int j, std::function<void(const INDEX, const SCALAR)> lambda){
      iipair ixpair=colix(j); int J=ixpair.first; int subj=ixpair.second; int ioffset=0; 
      for(int I=0; I<nb; ioffset+=blocks[I]->nrows, I++)
	blocks[J*nb+I]->for_each_in_column(subj,[ioffset,&lambda](const int i, const SCALAR v){lambda(ioffset+i,v);});
    }	
    
    bool isFilled(const iipair& I, const iipair& J) const {
      return (blocks[I.first+J.first*nb])->isFilled(I.second,J.second);}
    bool isFilled(const INDEX i, const INDEX j) const {
      return isFilled(rowix(i),colix(j));}

    
  private: // indexing

    iipair rowix(const int i) const{
      assert(i<nrows); int t=0; int I=0; while(i>=t+blocks[I]->nrows) t+=blocks[I++]->nrows; 
      return iipair(I,i-t);}
    iipair colix(const int i) const{
      assert(i<ncols); int t=0; int I=0; while(i>=t+blocks[I*nb]->ncols) t+=blocks[(I++)*nb]->ncols;
      return iipair(I,i-t);}


  public: // views

    auto view_of_column(const int J, const int j) -> BlockedVector<decltype(blocks[0]->view_of_column(j))>{
      BlockedVector<decltype(blocks[0]->view_of_column(j))> V(nb,nrows);
      for(int I=0; I<nb; I++)
	V.blocks[I]=blocks[I+J*nb]->view_of_column(j);
      return V;
    }

    auto view_of_column(const int J, const int j) const -> BlockedVector<decltype(blocks[0]->view_of_column(j))>{
      BlockedVector<decltype(blocks[0]->view_of_column(j))> V(nb,nrows);
      for(int I=0; I<nb; I++)
	V.blocks[I]=temp2ptr(blocks[I+J*nb]->view_of_column(j));
      return V;
    }

    Detached<const BlockedMatrix> view_of_tower(const int J) const{
      assert(J<mb);
      Detached<const BlockedMatrix> R(BlockedMatrix(nb,1));
      for(int I=0; I<nb; I++) R.block_ptr(I,0)=const_cast<MATRIX*>(block_ptr(I,J));
      return R;
    }


  public: // remapping 
    
    BlockedMatrix<Cmatrix> pullCols(const BindexMap& map) const;
    BlockedMatrix<Cmatrix> pullCols(const BtoBindexMap& map) const;

    BlockedMatrix<Cmatrix> pullRows(const BindexMap& map) const;
    BlockedMatrix<Cmatrix> pullRows(const BtoBindexMap& map) const;
    
      
  public: // in-place arithmetic 
    
    BlockedMatrix<MATRIX>& operator+=(const SCALAR x);
    BlockedMatrix<MATRIX>& operator-=(const SCALAR x);
    BlockedMatrix<MATRIX>& operator*=(const SCALAR x);
    BlockedMatrix<MATRIX>& operator/=(const SCALAR x);

    template<class MATRIX2>  
    BlockedMatrix<MATRIX>& operator+=(const BlockedMatrix<MATRIX2>& x);
    template<class MATRIX2>  
    BlockedMatrix<MATRIX>& operator-=(const BlockedMatrix<MATRIX2>& x);


  public: // arithmetic
    
    template<class MATRIX2>
    auto operator+(const BlockedMatrix<MATRIX2>& x) const -> BlockedMatrix<decltype((*blocks[0])+x)>; 
    template<class MATRIX2>
    auto operator-(const BlockedMatrix<MATRIX2>& x) const -> BlockedMatrix<decltype((*blocks[0])-x)>; 
    
    template<class MATRIX2>
    auto operator*(const BlockedMatrix<MATRIX2>& x) const -> BlockedMatrix<decltype((*blocks[0])*(*x.blocks[0]))>;
    template<class MATRIX2>
    auto dot(const BlockedMatrix<MATRIX2>& x) const -> BlockedMatrix<decltype(blocks[0]->dot(*x.blocks[0]))>;
    
    template<class VECTOR>
    auto operator*(const BlockedVector<VECTOR>& x) const -> BlockedVector<decltype((*blocks[0])*(*x.blocks[0]))>;
    template<class VECTOR>
    auto dot(const BlockedVector<VECTOR>& x) const -> BlockedVector<decltype((*blocks[0])*(*x.blocks[0]))>;


  public: // other 

    BlockStructure rowStructure() const{
      assert(mb>0);
      BlockStructure st(nb);
      for(int i=0; i<nb; i++) st(i)=blocks[i]->ncols;
      return st;
    }

    BlockStructure colStructure() const{
      assert(mb>0);
      BlockStructure st(mb);
      for(int j=0; j<mb; j++) st(j)=blocks[j*nb]->ncols;
      return st;
    }
    
    void normalizeColumns();



  };


  // ---- Copying ---------------------------------------------------------------------------------------------------


  template<class MATRIX>
  BlockedMatrix<MATRIX>::BlockedMatrix(const BlockedMatrix<MATRIX>& X): 
    BlockedMatrix(X.nb,X.mb,X.nrows,X.ncols){
    MultiLoop(nb*mb,[this,&X](const int i){blocks[i]=new MATRIX(*X.blocks[i]);});
  }

  template<class MATRIX>
  BlockedMatrix<MATRIX>::BlockedMatrix(BlockedMatrix<MATRIX>&& X): 
    BlockedMatrix(snatch(X.nb),snatch(X.mb),snatch(X.nrows),snatch(X.ncols)){
    for(int i=0; i<nb*mb; i++) {blocks[i]=X.blocks[i]; X.blocks[i]=nullptr;}
    delete[] X.blocks; 
  }

  template<class MATRIX>
  BlockedMatrix<MATRIX>& BlockedMatrix<MATRIX>::operator=(const BlockedMatrix<MATRIX>& X){
    for(int i=0; i<nb*mb; i++) delete blocks[i]; delete[] blocks;
    mb=X.mb; nb=X.nb; nrows=X.nrows; ncols=X.ncols; 
    blocks=new MATRIX*[nb*mb];
    MultiLoop(nb*mb,[this,X](const int i){blocks[i]=new MATRIX(*X.blocks[i]);});
    return *this;
  }

  template<class MATRIX>
  BlockedMatrix<MATRIX>& BlockedMatrix<MATRIX>::operator=(BlockedMatrix<MATRIX>&& X){
    for(int i=0; i<nb*mb; i++) delete blocks[i]; delete[] blocks;
    mb=snatch(X.mb); nb=snatch(X.nb); nrows=snatch(X.nrows); ncols=snatch(X.ncols);
    blocks=new MATRIX*[nb*mb];
    for(int i=0; i<nb*mb; i++){blocks[i]=X.blocks[i]; X.blocks[i]=nullptr;}; delete[] X.blocks;
    return *this;
  }

  template<class MATRIX>
  BlockedMatrix<MATRIX> BlockedMatrix<MATRIX>::copy() const{
    BlockedMatrix M(nb,mb,nrows,ncols);
    MultiLoop(nb*mb,[this,&M](const int i){M.blocks[i]=new MATRIX(blocks[i]->copy());});
    return M;
  }

  template<class MATRIX>
  BlockedMatrix<MATRIX> BlockedMatrix<MATRIX>::shallow() const{
    BlockedMatrix M(nb,mb,nrows,ncols);
    for(int i=0; i<nb*mb; i++) M.blocks[i]=blocks[i];
    return M;
  }

  template<class MATRIX>
  BlockedMatrix<MATRIX>::~BlockedMatrix(){
    for(int i=0; i<nb*mb; i++) delete blocks[i]; delete blocks;
  }


  // ---- Constructors ----------------------------------------------------------------------------------------------


  template<class MATRIX>
  BlockedMatrix<MATRIX>::BlockedMatrix(const int _nb, const int _mb): 
    Matrix(0,0), nb(_nb), mb(_mb){
    blocks=new MATRIX*[nb*mb]; // for(int i=0; i<nb*mb; i++) blocks[i]=new MATRIX();
  }

  template<class MATRIX>
  BlockedMatrix<MATRIX>::BlockedMatrix(const int _nb, const int _mb, const int _nrows, const int _ncols): 
    Matrix(_nrows,_ncols), nb(_nb), mb(_mb){
    blocks=new MATRIX*[nb*mb]; // for(int i=0; i<nb*mb; i++) blocks[i]=new MATRIX();
  }

  template<class MATRIX>
  BlockedMatrix<MATRIX>::BlockedMatrix(const BlockStructure& B1, const BlockStructure& B2):
    BlockedMatrix<MATRIX>(B1.size(),B2.size(),0,0){
    for(int i=0; i<nb; i++)
      for(int j=0; j<mb; j++)
	blocks[i+nb*j]=new MATRIX(B1[i],B2[j]);
    for(int i=0; i<nb; i++) nrows+=B1[i];
    for(int j=0; j<mb; j++) ncols+=B1[j];
  }
    

  // ---- Named Constructors ---------------------------------------------------------------------------------------
  
  
  template<class MATRIX>
  BlockedMatrix<MATRIX> BlockedMatrix<MATRIX>::Zero(const BlockStructure& B1, const BlockStructure& B2){
    BlockedMatrix<MATRIX> M(B1.size(),B2.size(),0,0);
    for(int i=0; i<M.nb; i++)
      for(int j=0; j<M.mb; j++)
	M.blocks[i+M.nb*j]=new MATRIX(MATRIX::Zero(B1[i],B2[j]));
    for(int i=0; i<M.nb; i++) M.nrows+=B1[i];
    for(int j=0; j<M.mb; j++) M.ncols+=B1[j];
    return M;
  }

  template<class MATRIX>
  BlockedMatrix<MATRIX> BlockedMatrix<MATRIX>::Filled(const BlockStructure& B1, const BlockStructure& B2, const SCALAR t){
    BlockedMatrix<MATRIX> M(B1.size(),B2.size(),0,0);
    for(int i=0; i<M.nb; i++)
      for(int j=0; j<M.mb; j++)
	M.blocks[i+M.nb*j]=new MATRIX(MATRIX::Filled(B1[i],B2[j],t));
    for(int i=0; i<M.nb; i++) M.nrows+=B1[i];
    for(int j=0; j<M.mb; j++) M.ncols+=B1[j];
    return M;
  }

  template<class MATRIX>
  BlockedMatrix<MATRIX> BlockedMatrix<MATRIX>::Identity(const BlockStructure& B){
    BlockedMatrix<MATRIX> M(B.size(),B.size(),0,0);
    for(int i=0; i<M.nb; i++)
      for(int j=0; j<M.mb; j++)
	if(i==j) M.blocks[i+M.nb*j]=new MATRIX(MATRIX::Identity(B[i]));
	else M.blocks[i+M.nb*j]=new MATRIX(MATRIX::Zero(B[i],B[j]));
    for(int i=0; i<M.nb; i++) M.nrows+=B[i]; M.ncols=M.nrows;
    return M;
  }

  template<class MATRIX>
  BlockedMatrix<MATRIX> BlockedMatrix<MATRIX>::Uniform(const BlockStructure& B1, const BlockStructure& B2){
    BlockedMatrix<MATRIX> M(B1.size(),B2.size(),0,0);
    for(int i=0; i<M.nb; i++)
      for(int j=0; j<M.mb; j++)
	M.blocks[i+M.nb*j]=new MATRIX(MATRIX::Uniform(B1[i],B2[j]));
    for(int i=0; i<M.nb; i++) M.nrows+=B1[i];
    for(int j=0; j<M.mb; j++) M.ncols+=B2[j];
    return M;
  }

  template<class MATRIX>
  BlockedMatrix<MATRIX> BlockedMatrix<MATRIX>::Gaussian(const BlockStructure& B1, const BlockStructure& B2){
    BlockedMatrix<MATRIX> M(B1.size(),B2.size(),0,0);
    for(int i=0; i<M.nb; i++)
      for(int j=0; j<M.mb; j++)
	M.blocks[i+M.nb*j]=new MATRIX(MATRIX::Gaussian(B1[i],B2[j]));
    for(int i=0; i<M.nb; i++) M.nrows+=B1[i];
    for(int j=0; j<M.mb; j++) M.ncols+=B2[j];
    return M;
  }

  template<class MATRIX>
  BlockedMatrix<MATRIX> BlockedMatrix<MATRIX>::Bernoulli(const BlockStructure& B1, const BlockStructure& B2, const double p){
    BlockedMatrix<MATRIX> M(B1.size(),B2.size(),0,0);
    for(int i=0; i<M.nb; i++)
      for(int j=0; j<M.mb; j++)
	M.blocks[i+M.nb*j]=new MATRIX(MATRIX::Bernoulli(B1[i],B2[j],p));
    for(int i=0; i<M.nb; i++) M.nrows+=B1[i];
    for(int j=0; j<M.mb; j++) M.ncols+=B2[j];
    return M;
  }


  // ---- In-place Arithmetic --------------------------------------------------------------------------------------


  template<class MATRIX>
  BlockedMatrix<MATRIX>& BlockedMatrix<MATRIX>::operator+=(const SCALAR x){
    for(int i=0; i<nb*mb; i++) (*blocks[i])+=x; return *this;}

  template<class MATRIX>
  BlockedMatrix<MATRIX>& BlockedMatrix<MATRIX>::operator-=(const SCALAR x){
    for(int i=0; i<nb*mb; i++) (*blocks[i])-=x; return *this;}

  template<class MATRIX>
  BlockedMatrix<MATRIX>& BlockedMatrix<MATRIX>::operator*=(const SCALAR x){
    for(int i=0; i<nb*mb; i++) (*blocks[i])*=x; return *this;}

  template<class MATRIX>
  BlockedMatrix<MATRIX>& BlockedMatrix<MATRIX>::operator/=(const SCALAR x){
    for(int i=0; i<nb*mb; i++) (*blocks[i])/=x; return *this;}

  template<class MATRIX>
  template<class MATRIX2>
  BlockedMatrix<MATRIX>& BlockedMatrix<MATRIX>::operator+=(const BlockedMatrix<MATRIX2>& x){
    assert(x.nb==nb); assert(x.mb==mb);
    for(int i=0; i<nb*mb; i++) (*blocks[i])+=(*x.blocks[i]);
    return *this;
  }

  template<class MATRIX>
  template<class MATRIX2>
  BlockedMatrix<MATRIX>& BlockedMatrix<MATRIX>::operator-=(const BlockedMatrix<MATRIX2>& x){
    assert(x.nb==nb); assert(x.mb==mb);
    for(int i=0; i<nb*mb; i++) (*blocks[i])-=(*x.blocks[i]);
    return *this;
  }


  // ---- Arithmetic -----------------------------------------------------------------------------------------------


  template<class MATRIX>
  template<class MATRIX2>
  auto BlockedMatrix<MATRIX>::operator+(const BlockedMatrix<MATRIX2>& x) const 
    ->BlockedMatrix<decltype((*blocks[0])+x)>{
    assert(x.nb==nb); assert(x.mb==mb);
    BlockedMatrix r(nb,mb,nrows,ncols);
    for(int i=0; i<nb*mb; i++) 
      r.blocks[i]=temp2ptr((*blocks[i])+(*x.blocks[i]));
    return r;
  }

  template<class MATRIX>
  template<class MATRIX2>
  auto BlockedMatrix<MATRIX>::operator-(const BlockedMatrix<MATRIX2>& x) const 
    ->BlockedMatrix<decltype((*blocks[0])-x)>{
    assert(x.nb==nb); assert(x.mb==mb);
    BlockedMatrix r(nb,mb,nrows,ncols);
    for(int i=0; i<nb*mb; i++) 
      r.blocks[i]=temp2ptr((*blocks[i])+(*x.blocks[i]));
    return r;
  }

  template<class MATRIX>
  template<class MATRIX2>
  auto BlockedMatrix<MATRIX>::operator*(const BlockedMatrix<MATRIX2>& x) const 
    ->BlockedMatrix<decltype((*blocks[0])*(*x.blocks[0]))>{
    assert(mb==x.nb); assert(mb>0);
    BlockedMatrix r(nb,x.mb,nrows,x.ncols);
    for(int i=0; i<nb; i++)
      for(int j=0; j<x.mb; j++){
	r.blocks[i+j*nb]=temp2ptr((*blocks[i])*(*x.blocks[j*x.nb]));
	for(int k=1; k<mb; k++)
	  (*r.blocks[i+j*nb])+=((*blocks[i+k*nb])*(*x.blocks[k+j*x.nb]));
      }
    return r;
  }

  template<class MATRIX>
  template<class MATRIX2>
  auto BlockedMatrix<MATRIX>::dot(const BlockedMatrix<MATRIX2>& x) const 
    ->BlockedMatrix<decltype(blocks[0]->dot(*x.blocks[0]))>{
    assert(nb==x.nb); assert(nb>0); // TODO: check for A.dot(A)
    BlockedMatrix<decltype(blocks[0]->dot(*x.blocks[0]))> r(mb,x.mb,ncols,x.ncols);
    for(int i=0; i<mb; i++)
      for(int j=0; j<x.mb; j++){
	r.blocks[i+j*mb]=temp2ptr((*blocks[i*nb])*(*x.blocks[j*x.nb]));
	for(int k=1; k<mb; k++)
	  (*r.blocks[i+j*mb])+=((*blocks[k+i*nb])*(*x.blocks[k+j*x.nb]));
      }
    return r;
  }
  
  template<class MATRIX>
  template<class VECTOR>
  auto BlockedMatrix<MATRIX>::operator*(const BlockedVector<VECTOR>& x) const 
    ->BlockedVector<decltype((*blocks[0])*(*x.blocks[0]))>{
    assert(mb==x.nb); assert(mb>0);
    BlockedVector<decltype((*blocks[0])*(*x.blocks[0]))> v; // TODO
    for(int i=0; i<nb; i++){
      v.blocks[i]=temp2ptr((*blocks[i])*(*x.blocks[0]));
	for(int j=1; j<mb; j++)
	  (*v.blocks[i])+=((*blocks[i+j*nb])*(*x.blocks[j]));
    }
    return v;
  }

  template<class MATRIX>
  template<class VECTOR>
  auto BlockedMatrix<MATRIX>::dot(const BlockedVector<VECTOR>& x) const 
    ->BlockedVector<decltype((*blocks[0])*(*x.blocks[0]))>{
    assert(nb==x.nb); assert(nb>0); assert(mb>0); 
    BlockedVector<decltype((*blocks[0]).dot(*x.blocks[0]))> v(mb,ncols); 
    for(int I=0; I<mb; I++){
      v.blocks[I]=temp2ptr((*blocks[I*nb]).dot(*x.blocks[0]));
	for(int J=1; J<nb; J++)
	  (*v.blocks[I])+=((*blocks[J+I*nb])*(*x.blocks[J]));
    }
    return v;
  }


  // ---- Other methods ---------------------------------------------------------------------------------------------


  template<class MATRIX>
  void BlockedMatrix<MATRIX>::normalizeColumns(){
    MultiLoop(mb,[this](const int J){
		int ncols=blocks[J*nb]->ncols;
		Cmatrix norms(nb,ncols);
		MultiLoop(nb,[this,&norms,J,ncols](const int I){
			    for(int j=0; j<ncols; j++) norms(I,j)=blocks[I+J*nb]->view_of_column(j).norm2();});
		Cvector normalizers(ncols);
		for(int j=0; j<ncols; j++) {
		  SCALAR t=sqrt(norms.view_of_column(j).norm2()); 
		  if(t!=0) normalizers.array[j]=1.0/t; 
		  else normalizers.array[j]=0;}
		MultiLoop(nb,[this,&normalizers,J,ncols](const int I){
			    for(int j=0; j<ncols; j++) blocks[I+J*nb]->multiplyColsBy(normalizers);});

	      });
  }


  // ---- Global functions ------------------------------------------------------------------------------------------
  

  template<class MATRIX>
  BlockedMatrix<MATRIX> operator*(const SCALAR c, const BlockedMatrix<MATRIX>& M){return M*c;}
  
  
} // namespace Mondrian

#endif
