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


#ifndef _Matrix
#define _Matrix

//#include "MatrixLike.hpp"
#include "Cvector.hpp"
#include "matrix_expr.hpp"

namespace Mondrian{

  class Matrix{
  public:

    int nrows;
    int ncols;


  public: // constructors 

    Matrix(): nrows(0), ncols(0) {}

    Matrix(const int _nrows, const int _ncols):
      nrows(_nrows), ncols(_ncols){}

    virtual ~Matrix(){}


  public: // attributes

    virtual bool isSparseFormat() const {return false;}
    virtual bool isSymmetricFormat() const {return false;}
    virtual bool isMultithreaded() const{return false;}


  public: // comparisons

    template<class MATRIX>
    bool operator==(const MATRIX& X) const{ 
      if(X.nrows!=nrows) return false; 
      if(X.ncols!=ncols) return false;
      for(int i=0; i<nrows; i++) 
	for(int j=0; j<ncols; j++) 
	  if(read(i,j)!=X.read(i,j)) return false;
      return true;}


  public: // element access    

    virtual SCALAR read(const int i, const int j) const=0; 
    virtual SCALAR operator()(const int i, const int j) const=0; 
    virtual void  set(const int i, const int j, const SCALAR v)=0; 

    template<class VECTOR>
    VECTOR row(const int i) const {assert(i<nrows);
      VECTOR x(ncols); for(int j=0; j<ncols; j++) x(j)=operator()(i,j); return x;}

    template<class VECTOR>
    VECTOR column(const int j) const {assert(j<ncols);
      VECTOR x(nrows); for(int i=0; i<nrows; i++) x(i)=operator()(i,j); return x;}
    
    template<class VECTOR>
    VECTOR diag() const {int t=min(nrows,ncols);
      VECTOR x(t); for(int i=0; i<t; i++) x(i)=operator()(i,i); return x;}

    virtual bool isFilled(const int i, const int j) const=0;

    virtual int nFilled() const{
      int t=0; 
      for(int i=0; i<nrows; i++) 
	for(int j=0; j<ncols; j++)
	  if((*this).isFilled(i,j))
	    t++;
      return t;
    };


  public: // iterators

    //virtual void for_each(std::function<void(const INDEX, const INDEX, SCALAR&)> lambda){
    //  for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) lambda(i,j,(*this)(i,j));}
    //virtual void for_each_in_row(const int i, std::function<void(const INDEX, SCALAR&)> lambda){
    //  for(int j=0; j<ncols; j++) lambda(j,(*this)(i,j));}
    //virtual void for_each_in_column(const int j, std::function<void(const INDEX, SCALAR&)> lambda){
    //  for(int i=0; i<nrows; i++) lambda(i,(*this)(i,j));}

    virtual void for_each(std::function<void(const INDEX, const INDEX, const SCALAR)> lambda) const{
      for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) lambda(i,j,read(i,j));}
    virtual void for_each_in_row(const int i, std::function<void(const INDEX, const SCALAR)> lambda) const{
      for(int j=0; j<ncols; j++) lambda(j,read(i,j));}
    virtual void for_each_in_column(const int j, std::function<void(const INDEX, const SCALAR)> lambda) const{
      for(int i=0; i<nrows; i++) lambda(i,read(i,j));}
  
    //virtual void for_each_filled(std::function<void(const INDEX, const INDEX, SCALAR&)> lambda){
    //  for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) if(isFilled(i,j)) lambda(i,j,(*this)(i,j));}
    //virtual void for_each_filled_in_row(const int i, std::function<void(const INDEX, SCALAR&)> lambda){
    //  for(int j=0; j<ncols; j++) if(isFilled(i,j)) lambda(j,(*this)(i,j));}
    //virtual void for_each_filled_in_column(const int j, std::function<void(const INDEX, SCALAR&)> lambda){
    //  for(int i=0; i<nrows; i++) if(isFilled(i,j)) lambda(i,(*this)(i,j));}

    virtual void for_each_filled(std::function<void(const INDEX, const INDEX, const SCALAR)> lambda) const{
      for(int i=0; i<nrows; i++) for(int j=0; j<ncols; j++) if(isFilled(i,j)) lambda(i,j,read(i,j));}
    virtual void for_each_filled_in_row(const int i, std::function<void(const INDEX, const SCALAR)> lambda) const{
      for(int j=0; j<ncols; j++) if(isFilled(i,j)) lambda(j,read(i,j));}
    virtual void for_each_filled_in_column(const int j, std::function<void(const INDEX, const SCALAR)> lambda) const{
      for(int i=0; i<nrows; i++) {cout<<i<<endl; if(isFilled(i,j)) lambda(i,read(i,j));}}
  
    
  public: // in-place operations

    //virtual Matrix& multiplyRowsBy(const Cvector& v)=0;
    //virtual Matrix& divideRowsBy(const Cvector& v)=0;
    //virtual Matrix& multiplyColsBy(const Cvector& v)=0;
    //virtual Matrix& divideColsBy(const Cvector& v)=0;


  public: // scalar valued methods

    virtual int nnz() const {
      int t=0; 
      for(int i=0; i<nrows; i++) 
	for(int j=0; j<ncols; j++)
	  if((*this)(i,j)!=0) t++;
      return t;
    };
  
    virtual SCALAR norm2() const {
      double t=0; 
      for(int i=0; i<nrows; i++) 
	for(int j=0; j<ncols; j++)
	  t+=(*this)(i,j)*(*this)(i,j);
      return t;
    };

    //virtual SCALAR inp_of_columns(const int j1, const int j2) const =0;
    //virtual SCALAR inp_of_rows(const int i1, const int i2) const =0;


  public: // I/O
    
    virtual string str() const{
      ostringstream oss;
      oss.precision(3); 
      oss.setf(ios_base::fixed, ios_base::floatfield);
      for(int i=0; i<nrows; i++){oss<<"[ ";
	for(int j=0; j<ncols; j++) 
	  if(isFilled(i,j)) {oss.width(6); oss<<(*this)(i,j)<<" ";}
	  else oss<<" 0.0   ";
	oss<<" ]\n";}
      return oss.str();
    }

    //virtual void print() const{
    //  cout<<str()<<endl;}


  public: // Python 

    const char* __str__() {
      ostringstream ostream; 
      ostream<<str()<<endl; 
      return ostream.str().c_str();
    }

    
  };


  inline ostream& operator<<(ostream& stream, const Matrix& x){stream<<x.str(); return stream;}

}

//#include "MatrixTranspose.hpp"


#endif 



    // Matrix(const int _nrows, const int _ncols): nrows(_nrows), ncols(_ncols){}
    //    virtual bool isDense() const=0;
