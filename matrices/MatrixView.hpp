#ifndef _MatrixView
#define _MatrixView

//#include "MatrixLike.hpp"
#include "IndexMap.hpp"
#include "Cmatrix.hpp"

namespace Mondrian{

template<class MATRIX>
class MatrixView: public Matrix{
public:

  MATRIX& M;
  IndexMap rmap;
  IndexMap cmap;

public: // copying and assignment

  MatrixView(const MatrixView& x): Matrix(x.nrows,x.ncols), M(x.M), rmap(x.rmap), cmap(x.cmap){}
  MatrixView(MatrixView&& x): Matrix(x.nrows,x.ncols), M(x.M), rmap(std::move(x.rmap)), cmap(std::move(x.cmap)){}

  template<class MATRIX2>
  MatrixView<MATRIX>& operator=(const MATRIX2& B){
    assert(nrows==B.nrows); assert(ncols==B.ncols);
    for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) M(rmap(i),cmap(j))=B(i,j);
    return *this;
  }

  ~MatrixView(){}


public: // constructors

  MatrixView(MATRIX& _M): 
    MatrixView(_M,IndexMap::Identity(_M.nrows),IndexMap::Identity(_M.ncols)){}

  MatrixView(MATRIX& _M, const IndexMap& _rmap, const IndexMap& _cmap): 
    Matrix(_rmap.nsource,_cmap.nsource), M(_M), rmap(_rmap), cmap(_cmap){}

  MatrixView(MATRIX& _M, IndexMap&& _rmap, IndexMap&& _cmap): 
    Matrix(_rmap.nsource,_cmap.nsource), M(_M), rmap(std::move(_rmap)), cmap(std::move(_cmap)){}


public: // attributes

  bool isSparseFormat() const {return M.isSparseFormat();}
  bool isSymmetricFormat() const {return M.isSparseFormat();}


public: // conversions

  operator Cmatrix() const{
    Cmatrix A(nrows,ncols);
    for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
	A(i,j)=M(rmap(i),cmap(j));
    return A;
  }


public: // element access

  SCALAR operator()(const int i, const int j) const{ 
    assert(i<nrows); assert(j<ncols); return M(rmap(i),cmap(j));} 
  SCALAR read(const int i, const int j) const{ 
    assert(i<nrows); assert(j<ncols); return M.read(rmap(i),cmap(j));;} 

  void set(const int i, const int j, const SCALAR v){
    assert(i<nrows); assert(j<ncols); M.set(rmap(i),cmap(j),v);}
  SCALAR& operator()(const int i, const int j){
    assert(i<nrows); assert(j<ncols); return M(rmap(i),cmap(j));} 

  Cvector row(const int i) const {assert(i<nrows);
    Cvector x(ncols); for(int j=0; j<ncols; j++) x.array[j]=M(rmap(i),cmap(j)); return x;}
  Cvector column(const int j) const {assert(j<ncols);
    Cvector x(nrows); for(int i=0; i<nrows; i++) x.array[i]=M(rmap(i),cmap(j)); return x;}
  Cvector diag() const {int t=min(nrows,ncols);
    Cvector x(t); for(int i=0; i<t; i++) x.array[i]=M(rmap(i),cmap(i)); return x;}

  bool isFilled(const int i, const int j) const {return M.isFilled(rmap(i),cmap(j));}
  int nFilled() const {return nrows*ncols;} // TODO 


public: // iterators

  void for_each(std::function<void(const INDEX, const INDEX, SCALAR&)> lambda){
    for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) lambda(i,j,M(rmap(i),cmap(j)));}
  void for_each(std::function<void(const INDEX, const INDEX, const SCALAR)> lambda) const{
    for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) lambda(i,j,M(rmap(i),cmap(j)));}
  void for_each_in_row(const int i, std::function<void(const INDEX, SCALAR&)> lambda){
    for(int j=0; j<ncols; j++) lambda(i,M(rmap(i),cmap(j)));}
  void for_each_in_row(const int i, std::function<void(const INDEX, const SCALAR)> lambda) const{
    for(int j=0; j<ncols; j++) lambda(i,M(rmap(i),cmap(j)));}
  void for_each_in_column(const int j, std::function<void(const INDEX, SCALAR&)> lambda){
    for(int i=0; i<nrows; i++) lambda(i,M(rmap(i),cmap(j)));}
  void for_each_in_column(const int j, std::function<void(const INDEX, const SCALAR)> lambda) const{
    for(int i=0; i<nrows; i++) lambda(i,M(rmap(i),cmap(j)));}


public: // comparisons

  bool operator==(const Matrix& X) const{ 
    if(X.nrows!=nrows) return false; if(X.ncols!=ncols) return false;
    for(int i=0; i<nrows; i++)
      for(int j=0; j<ncols; j++)
	if(M(rmap(i),cmap(j))!=X(i,j)) return false;
    return true;}

  bool operator!=(const Cmatrix& X) const {return !((*this)==X);}


public: // in-place arithmetic

  MatrixView<MATRIX>& operator+=(const SCALAR& x){
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
	M(rmap(i),cmap(j))+=x; 
    return *this;}
  MatrixView<MATRIX>& operator-=(const SCALAR& x){
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
	M(rmap(i),cmap(j))-=x; 
    return *this;}
  MatrixView<MATRIX>& operator*=(const SCALAR& x){
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
	M(rmap(i),cmap(j))*=x; 
    return *this;}
  MatrixView<MATRIX>& operator/=(const SCALAR& x){
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
	M(rmap(i),cmap(j))/=x; 
    return *this;}

  MatrixView<MATRIX>& operator+=(const Matrix& X){
    assert(X.nrows==nrows); assert(X.ncols==ncols);
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
	M(rmap(i),cmap(j))+=X(i,j); 
    return *this;}
  MatrixView<MATRIX>& operator-=(const Matrix& X){
    assert(X.nrows==nrows); assert(X.ncols==ncols);
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
	M(rmap(i),cmap(j))-=X(i,j); 
    return *this;}


public: // in-place operations 

  MatrixView<MATRIX>& multiplyRowsBy(const Cvector& v){
    assert(v.n==nrows);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	M(rmap(i),cmap(j))=v.array[i]*M(rmap(i),cmap(j));
    return *this;}

  MatrixView<MATRIX>& multiplyColsBy(const Cvector& v){
    assert(v.n==ncols);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	M(rmap(i),cmap(j))=v.array[j]*M(rmap(i),cmap(j));
    return *this;}

  MatrixView<MATRIX>& divideRowsBy(const Cvector& v){
    assert(v.n==nrows);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	M(rmap(i),cmap(j))=M(rmap(i),cmap(j))/v.array[i];
    return *this;}

  MatrixView<MATRIX>& divideColsBy(const Cvector& v){
    assert(v.n==ncols);
    for(int j=0; j<ncols; j++)
      for(int i=0; i<nrows; i++)
	M(rmap(i),cmap(j))=M(rmap(i),cmap(j))/v.array[j];
    return *this;}

  
public: // scalar-valued methods

  int nnz() const {int t=0; 
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
	if(M(rmap(i),cmap(j))!=0) t++; 
    return t;
  }
  
  SCALAR norm2() const{
    SCALAR t=0; 
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
	t+=M(rmap(i),cmap(j))*M(rmap(i),cmap(j)); 
    return t;
  }

  SCALAR diff2(const Cmatrix& X) const{
    assert(X.nrows==nrows); assert(X.ncols==ncols); 
    SCALAR t=0;
    for(int i=0; i<nrows; i++) 
      for(int j=0; j<ncols; j++) 
	t+=(X(i,j)-M(rmap(i),cmap(j)))*(X(i,j)-M(rmap(i),cmap(j))); 
    return t;
  }


};



}

#endif
