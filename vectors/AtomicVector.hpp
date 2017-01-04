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


#ifndef _AtomicVector
#define _AtomicVector

namespace Mondrian{


  template<class VECTOR>
  class AtomicVector<VECTOR>: public VECTOR{
  public:
    
    using VECTOR::VECTOR;

    mutex mx;
    

  public: // element access

    void for_each(std::function<void(const INDEX, SCALAR&)> lambda){
      lock_guard<mutex> lock(mx);
      VECTOR::for_each(lambda);
    }


  public: // in-place arithmetic

    VECTOR& operator+=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      VECTOR::operator+=(x);
      return *this;
    }

    VECTOR& operator-=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      VECTOR::operator-=(x);
      return *this;
    }

    VECTOR& operator*=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      VECTOR::operator+=(x);
      return *this;
    }

    VECTOR& operator/=(const SCALAR x){
      lock_guard<mutex> lock(mx);
      VECTOR::operator-=(x);
      return *this;
    }

    VECTOR& operator+=(const VECTOR& x){
      lock_guard<mutex> lock(mx);
      VECTOR::operator+=(x);
      return *this;
    }

    VECTOR& operator-=(const VECTOR& x){
      lock_guard<mutex> lock(mx);
      VECTOR::operator-=(x);
      return *this;
    }


  public: // in-place operators

    template<class OPERATOR>
    void apply(const OPERATOR& op){
      lock_guard<mutex> lock(mx);
      VECTOR::apply(op);
    }

    template<class OPERATOR>
    void applyInv(const OPERATOR& op){
      lock_guard<mutex> lock(mx);
      VECTOR::apply(op);
    }


  };

}

#endif
