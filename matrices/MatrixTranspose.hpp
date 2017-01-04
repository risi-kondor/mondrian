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


#ifndef _MatrixTranspose
#define _MatrixTranspose

#include "Matrix.hpp"


//SCALAR transp(const SCALAR x){return x;}
//SCALAR& transp(SCALAR& x){return x;}

namespace Mondrian{


  template<class MATRIX>
  class MatrixTranspose: public Matrix{
  public:
    
    MatrixTranspose(MATRIX& _obj):Matrix(_obj.ncols,_obj.nrows),obj(_obj){};
    ~MatrixTranspose(){} // what to put here?


  public:
    
    MATRIX& obj;


  public:
    
    MATRIX& operator~() const {return obj;}
    MATRIX& transp() const {return obj;}
    

  public: // attributes
    
    bool isSparseFormat() const {return obj.isSparseFormat();}
    bool isSymmetricFormat() const {return false;}


  public: // element access    

    int nFilled() {return obj.nFilled();}

    SCALAR read(const int i, const int j) const {return obj.read(j,i);} 
    SCALAR operator()(const int i, const int j) const {return obj.read(j,i);} 

    void set(const int i, const int j, const SCALAR v) {return obj.set(j,i,v);}
    SCALAR& operator()(const int i, const int j) {return obj.element(j,i);}


    // TODO: exchange indices
    //void for_each(std::function<void(const INDEX, const INDEX, SCALAR&)> lambda) {
    //  obj.for_each(lambda);} 
    //void for_each(std::function<void(const INDEX, const INDEX, const SCALAR)> lambda) const {
    //  static_cast<const MATRIX&>(obj).for_each(lambda);}
    //void for_each_in_row(const int i, std::function<void(const INDEX, SCALAR&)> lambda) {
    //  obj.for_each_in_column(i,lambda);}
    //void for_each_in_row(const int i, std::function<void(const INDEX, const SCALAR)> lambda) const {
    //  static_cast<const MATRIX&>(obj).for_each_in_column(i,lambda);}
    //void for_each_in_column(const int j, std::function<void(const INDEX, SCALAR&)> lambda) {
    //  obj.for_each_in_row(j,lambda);}
    //void for_each_in_column(const int j, std::function<void(const INDEX, const SCALAR)> lambda) const {
    //  static_cast<const MATRIX&>(obj).for_each_in_row(j,lambda);}

    template<class VECTOR> VECTOR row(const int i) const {return obj.template column<VECTOR>(i);}
    template<class VECTOR> VECTOR column(const int j) const {return obj.template row<VECTOR>(j);}
    template<class VECTOR> VECTOR diag(const int j) const {return obj.template diag<VECTOR>();}

    bool isFilled(const int i, const int j) const {return obj.isFilled(j,i);}

  public: // in-place arithmetic
    
    MatrixTranspose<MATRIX>& operator+=(const SCALAR x){obj+=x; return *this;}
    MatrixTranspose<MATRIX>& operator-=(const SCALAR x){obj-=x; return *this;}
    MatrixTranspose<MATRIX>& operator*=(const SCALAR x){obj*=x; return *this;}
    MatrixTranspose<MATRIX>& operator/=(const SCALAR x){obj/=x; return *this;}
  
    template<class MATRIX2>
    MatrixTranspose<MATRIX>& operator+=(const MatrixTranspose<MATRIX2>& xT){obj+=xT.obj; return *this;}
    template<class MATRIX2>
    MatrixTranspose<MATRIX>& operator-=(const MatrixTranspose<MATRIX2>& xT){obj-=xT.obj; return *this;}


  public: // matrix-valued arithmetic

    template<class MATRIX2>
    auto operator+(const MATRIX2& x) -> MatrixTranspose<decltype(obj+transp(x))> const {
      return transp(obj+transp(x));}

    template<class MATRIX2>
    auto operator-(const MATRIX2& x) -> MatrixTranspose<decltype(obj-transp(x))> const {
      return transp(obj-transp(x));}

    template<class MATRIX2>
    auto operator*(const MATRIX2& x) -> MatrixTranspose<decltype(x*obj)> const {
      return transp(transp(x)*obj);}

    template<class MATRIX2>
    auto dot(const MATRIX2& x) -> decltype(obj*x) const {
      return obj*x;}

    template<class MATRIX2>
    auto operator+(const MatrixTranspose<MATRIX2>& xT) -> MatrixTranspose<decltype(obj+xT.obj)> const {
      return transp(obj+xT.obj);}

    template<class MATRIX2>
    auto operator-(const MatrixTranspose<MATRIX2>& xT) -> MatrixTranspose<decltype(obj-xT.obj)> const {
      return transp(obj-xT.obj);}
    
    template<class MATRIX2>
    auto operator*(const MatrixTranspose<MATRIX2>& xT) -> MatrixTranspose<decltype(xT.obj*obj)> const {
      return transp(xT.obj*obj);}

    template<class MATRIX2>
    auto dot(const MatrixTranspose<MATRIX2>& xT) -> decltype(obj*xT) const {
      return obj*xT;}


  public: // in-place operations

    template<class OPERATOR>
    void applyFromLeft(const OPERATOR& Q){obj.applyFromRight(Q);}
    template<class OPERATOR>
    void applyFromRight(const OPERATOR& Q){obj.applyFromLeft(Q);}
    template<class OPERATOR>
    void applyFromLeftInv(const OPERATOR& Q){obj.applyFromRightInv(Q);}
    template<class OPERATOR>
    void applyFromRightInv(const OPERATOR& Q){obj.applyFromLeftInv(Q);}


  public: // in-place operations

    MatrixTranspose<MATRIX>& multiplyRowsBy(const Cvector& v) {obj.multiplyColsBy(v); return *this;}
    MatrixTranspose<MATRIX>& divideRowsBy(const Cvector& v) {obj.divideColsBy(v); return *this;}
    MatrixTranspose<MATRIX>& multiplyColsBy(const Cvector& v) {obj.multiplyRowsBy(v); return *this;}
    MatrixTranspose<MATRIX>& divideColsBy(const Cvector& v) {obj.divideColsBy(v); return *this;}


  public: // scalar valued methods

    int nFilled() const {return obj.nFilled();}
    int nnz() const {return obj.nnz();}
    SCALAR norm2() const {return obj.norm2();}

    SCALAR inp_of_columns(const int j1, const int j2) const {return obj.inp_of_rows(j1,j2);}
    SCALAR inp_of_rows(const int i1, const int i2) const {return obj.inp_of_columns(i1,i2);}




};

}

//template<class MATRIX>
//SCALAR MatrixTranspose<MATRIX>::operator()(const int i, const int j) const {return obj(j,i);}
  
//template<class MATRIX>
//SCALAR& MatrixTranspose<MATRIX>::operator()(const int i, const int j){return obj(j,i);}



//MatrixTranspose transp(const Matrix& M){return MatrixTranspose(M);}
//MatrixTranspose transp(const Matrix& M){return MatrixTranspose(M);}




#endif
