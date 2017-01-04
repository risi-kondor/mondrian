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


#ifndef _AtomicMatrix
#define _AtomicMatrix

namespace Mondrian{


  template<class MATRIX>
  class AtomicMatrix<MATRIX>: public MATRIX{
  public:
    
    using MATRIX::MATRIX;

    mutex mx;

    
  public: // element access

    void for_each(std::function<void(const INDEX, const INDEX, SCALAR&)> lambda){
      lock_guard<mutex> lock(mx); VECTOR::for_each(lambda);}
    void for_each_in_row(const int i, std::function<void(const INDEX, SCALAR&)> lambda){
      lock_guard<mutex> lock(mx); VECTOR::for_eachin_row(i,lambda);}
    void for_each_in_column(const int j, std::function<void(const INDEX, SCALAR&)> lambda){
      lock_guard<mutex> lock(mx); VECTOR::for_each_in_column(j,lambda);}


  public: // in-place arithmetic 

    MATRIX& operator+=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      MATRIX::operator+=(x);
      return *this;
    }

    MATRIX& operator-=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      MATRIX::operator-=(x);
      return *this;
    }

    MATRIX& operator*=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      MATRIX::operator+=(x);
      return *this;
    }

    MATRIX& operator/=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      MATRIX::operator-=(x);
      return *this;
    }

    MATRIX& operator+=(const MATRIX& x){
      lock_guard<mutex> lock(mx);
      MATRIX::operator+=(x);
      return *this;
    }

    MATRIX& operator-=(const MATRIX& x){
      lock_guard<mutex> lock(mx);
      MATRIX::operator-=(x);
      return *this;
    }


  public: // in-place operations

    void multiplyRowsBy(const Cvector& v) {lock_guard<mutex> lock(mx); MATRIX::multiplyRowsBy(v);}
    void divideRowsBy(const Cvector& v) {lock_guard<mutex> lock(mx); MATRIX::divideRowsBy(v);}
    void multiplyColsBy(const Cvector& v) {lock_guard<mutex> lock(mx); MATRIX::multiplyColsBy(v);}
    void divideColsBy(const Cvector& v) {lock_guard<mutex> lock(mx); MATRIX::divideColsBy(v);}
      
    template<class OPERATOR>
    void apply(const OPERATOR& op){
      lock_guard<mutex> lock(mx);
      MATRIX::apply(op);
    }

    template<class OPERATOR>
    void applyInv(const OPERATOR& op){
      lock_guard<mutex> lock(mx);
      MATRIX::apply(op);
    }


  };

}

#endif
