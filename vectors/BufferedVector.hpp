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


#ifndef _BufferedVector
#define _BufferedVector

#include "Mondrian_base.hpp"
#include "Flushed.hpp"

namespace Mondrian{

  template<class VECTOR>
  class BufferedVector: public VECTOR{
  public:

    using VECTOR::n;
    
    vector<ivpair> buffer;
    mutex buf_mx;
    mutable mutex access_mx;
    
    
  public: // constructors 

    
    using VECTOR::VECTOR;


  public: // copying ------------------------------------------------------------------------------------------------


    BufferedVector(const BufferedVector<VECTOR>& x): VECTOR(flush(x)){}
    BufferedVector(const BufferedVector<VECTOR>& x, const _NoWarn dummy): VECTOR(flush(x),dummy){}
    BufferedVector(BufferedVector<VECTOR>&& x): VECTOR(std::move(x.flush())){}

    BufferedVector<VECTOR>& operator=(const BufferedVector<VECTOR>& x){
      lock_guard<mutex> guard(access_mx);
      ASSIGN_WARNING("BufferedVector<VECTOR>");
      buffer.clear();
      static_cast<VECTOR&>(*this)=flush(x);
      return *this;
    }
    
    BufferedVector<VECTOR>& operator=(const BufferedVector<VECTOR>&& x){
      lock_guard<mutex> guard(access_mx);
      ASSIGN_WARNING("BufferedVector<VECTOR>");
      buffer.clear();
      static_cast<VECTOR&>(*this)=std::move(flush(x));
      return *this;
    }
    
    BufferedVector<VECTOR> copy(){
      //lock_guard<mutex> guard(access_mx);
      return BufferedVector<VECTOR>(*this,_NoWarn()); 
    }

    
  public: // conversions --------------------------------------------------------------------------------------------


    BufferedVector(const VECTOR& x):
      VECTOR(x){}

    BufferedVector(VECTOR&& x):
      VECTOR(std::move(x)){}

    template<class VECTOR2>
    BufferedVector(const VECTOR2& x):
      VECTOR(x){}

    operator Cvector(){
      return Flushed< BufferedVector<VECTOR> >(*this);
    }


  public: // flushing -----------------------------------------------------------------------------------------------
    
    
    BufferedVector<VECTOR>& flush() const{
      const_cast<BufferedVector<VECTOR>&>(*this).flush();
      return *this;
    }

    BufferedVector<VECTOR>& flush(){
      lock_guard<mutex> guard(access_mx);
      for(auto& p:buffer) VECTOR::set(p.first,p.second);
      buffer.clear();
      return *this;
    }

    const BufferedVector<VECTOR>& flush_while_locked() const{
      const_cast<BufferedVector<VECTOR>&>(*this).flush_while_locked();
      return *this;
    }

    BufferedVector<VECTOR>& flush_while_locked(){
      for(auto& p:buffer) VECTOR::set(p.first,p.second);
      buffer.clear();
      return *this;
    }

    template<class TYPE>
    TYPE while_flushed(function<TYPE()> lambda) const{
      lock_guard<mutex> guard(access_mx);
      flush_while_locked();
      return lambda();
    }

    /*
    //template<class TYPE>
    bool while_flushed(function<bool()> lambda) const{
      lock_guard<mutex> guard(access_mx);
      flush_while_locked();
      return lambda();
    }

    SCALAR while_flushed(function<SCALAR()> lambda) const{
      lock_guard<mutex> guard(access_mx);
      flush_while_locked();
      return lambda();
    }
    */

    SCALAR& while_flushed_ref(function<SCALAR&()> lambda) const{
      lock_guard<mutex> guard(access_mx);
      flush_while_locked();
      return lambda();
    }


  public: // comparators --------------------------------------------------------------------------------------------


    bool operator==(const VECTOR& x) const{
      return while_flushed<bool>([this,&x]()->bool{return VECTOR::operator==(x);});
    }

    
  public: // element access -----------------------------------------------------------------------------------------
    
    SCALAR operator()(const int i) const {
      return while_flushed<SCALAR>([this,i]()->SCALAR{return VECTOR::operator()(i);});
    }
      
    SCALAR& operator()(const int i){
      lock_guard<mutex> guard(access_mx);
      flush_while_locked();
      return VECTOR::operator()(i);
    }
      
    /*
    SCALAR read(const int i) const {assert(i<n);
      auto p=find_ptr(i);
      if(p!=nullptr) return p->second;
      else return 0;
    }

    void set(const int i, const SCALAR v) {assert(i<n);
      auto p=find_ptr(i);
      if(p!=nullptr) {p->second=v; return;}
      chg.push_back(ivpair(i,v));
      return;
    }

    void set_if_filled(const int i, const SCALAR v) {assert(i<n);
      auto p=find_ptr(i);
      if(p!=nullptr) p->second=v;
      return;
    }

    bool isFilled(const int i) const {assert(i<n);
      return (find_ptr(i)!=nullptr);
    }
    
    int nFilled() const{
      return vec.size()+chg.size();
    }
    */

    void insert(const int i, const SCALAR v){assert(i<n);
      buffer.push_back(ivpair(i,v));
    }

    
    BufferedVector<VECTOR>& insert(const ivpair& p){
      buffer.push_back(p);
      return *this;
    }

    

  };

}

#endif
