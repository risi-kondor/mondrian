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


#ifndef _ActiveVector
#define _ActiveVector

#include "Vector.hpp"

namespace Mondrian{

  template<class VECTOR>
  class ActiveVector: public VECTOR{
  public:

    using VECTOR::n;

    using VECTOR::read;


    virtual void changed(const INDEX i, const SCALAR& v){
      cout<<"v("<<i<<")<-"<<v<<endl;
    }

    void changed(const INDEX i){
      changed(i,VECTOR::read(i));
    }


  public: // constructors -------------------------------------------------------------------------------------------


    ActiveVector(): VECTOR() {}

    ActiveVector(const int _n): VECTOR(_n) {}

    ActiveVector(const int _n, const SCALAR* _array): VECTOR(_n) {
      for(int i=0; i<n; i++) set(i,_array[i]);
    }
    
    ActiveVector(const initializer_list<SCALAR> list): VECTOR(list.size()){
      int i=0; for(SCALAR v:list) set(i++,v);}

    ActiveVector(const int _n, const initializer_list<ivpair> list): VECTOR(_n){ // _Zero() ??
      for(const ivpair& v:list) {assert(v.first<n); set(v.first)=v.second;}
    }
    

  public: // copying ------------------------------------------------------------------------------------------------


    ActiveVector(const ActiveVector<VECTOR>& x): VECTOR(x,_NoWarn()){
      COPY_WARNING("ActiveVector<VECTOR>");
    } 

    ActiveVector(const ActiveVector<VECTOR>& x, const _NoWarn dummy): VECTOR(x,_NoWarn()){
    } 

    ActiveVector(ActiveVector<VECTOR>&& x): VECTOR(std::move(x)){
      MOVE_WARNING("ActiveVector<VECTOR>");
    }

    ActiveVector& operator=(const ActiveVector<VECTOR>& x){
      VECTOR::operator=(x);
      return *this;
    }
    
    ActiveVector& operator=(ActiveVector<VECTOR>&& x){
      VECTOR::operator=(std::move(x));
      return *this;
    }
    
    ActiveVector<VECTOR> copy() const {
      return VECTOR::copy(); 
    }
    
    ActiveVector<VECTOR> shallow() const {
      return VECTOR::shallow(); 
    }
    
    void assign(const ActiveVector<VECTOR>& x) const {
      return VECTOR::assign(x); 
    }
    
    void detach(){VECTOR::detach();}


  public: // downcasting --------------------------------------------------------------------------------------------


    ActiveVector<VECTOR>(const VECTOR& x):
      VECTOR(x,_NoWarn()){
    }

    ActiveVector<VECTOR>(VECTOR&& x):
      VECTOR(std::move(x)){}


  public: // polymorphism -------------------------------------------------------------------------------------------
    

    static string classname() {return "ActiveVector<"+VECTOR::classname()+">";}
    virtual VectorSpecializerBase* specializer() const {
      return new VectorSpecializer< ActiveVector<VECTOR> >();}


  public: // comparisons --------------------------------------------------------------------------------------------
    

    bool operator==(const ActiveVector& x) const{
      return VECTOR::operator==(x);
    }


  public: // element writes -----------------------------------------------------------------------------------------


    SCALAR& operator()(const INDEX i){
      cout<<"Warning: ActiveVector::operator(int) unsafe"<<endl; 
      return VECTOR::operator()(i);
    }

    void set(const int i, const SCALAR v){
      VECTOR::set(i,v);
      changed(i,v);
    }

    void set_msafe(const int i, const SCALAR v){ // TODO
      VECTOR::set_msafe(i,v);
      changed(i,v);
    }

    SCALAR* ptr(const INDEX i){
      cout<<"Warning: ActiveVector::ptr(int) unsafe"<<endl; 
      return VECTOR::ptr(i);
    }

    
  public: // function mappings --------------------------------------------------------------------------------------
    

    void for_each(std::function<void(const INDEX, const SCALAR&)> lambda) const{
      VECTOR::for_each(lambda);}

    void for_each(std::function<void(const INDEX, SCALAR&)> lambda){
      VECTOR::for_each([this,&lambda](const int i, SCALAR& v){lambda(i,v); changed(i,this->read(i));});
    }

    void for_each_filled(std::function<void(const INDEX, const SCALAR&)> lambda) const{
      VECTOR::for_each_filled(lambda);}

    void for_each_filled(std::function<void(const INDEX, SCALAR&)> lambda){
      VECTOR::for_each([this,&lambda](const int i, SCALAR& v){lambda(i,v); changed(i,this->read(i));});
    }


  public: // Givens rotations ---------------------------------------------------------------------------------------


    ActiveVector<VECTOR>& apply(const GivensRotation& Q){
      VECTOR::apply(Q);
      changed(Q.i1);
      changed(Q.i2);
      return *this;
    }

    ActiveVector<VECTOR>& applyT(const GivensRotation& Q){
      VECTOR::applyT(Q);
      changed(Q.i1);
      changed(Q.i2);
      return *this;
    }


  public: // KpointOp<k> --------------------------------------------------------------------------------------------


    template<int k>
    ActiveVector<VECTOR>& apply(const KpointOp<k>& Q){
      VECTOR::apply(Q);
      for(int i=0; i<k; i++)
	changed(Q.map(i));
      return *this;
    }

    template<int k>
    ActiveVector<VECTOR>& applyT(const KpointOp<k>& Q){
      VECTOR::applyT(Q);
      for(int i=0; i<k; i++)
	changed(Q.map(i));
      return *this;
    }


  };

}

#endif
