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


#ifndef _IndexMap
#define _IndexMap

#include "Mondrian_base.hpp"

extern default_random_engine randomNumberGenerator;


namespace Mondrian{

class IndexMap: public Serializable{
public:

  int nsource; 
  INDEX* forward;


public:

  IndexMap(const int _nsource): nsource(_nsource) {forward=new INDEX[nsource];}

  IndexMap(const initializer_list<INDEX> list): IndexMap(list.size()){
    int i=0; for(auto j:list) forward[i++]=j;}


public: // copying 

  IndexMap(const IndexMap& x): nsource(x.nsource){
    forward=new INDEX[nsource]; 
    std::copy(x.forward,x.forward+nsource,forward);
  }

  IndexMap(IndexMap&& x): nsource(x.nsource){
    forward=snatchptr(x.forward); x.nsource=0;}
  
  IndexMap& operator=(const IndexMap& x){
    nsource=x.nsource; delete forward; forward=new INDEX[nsource]; 
    std::copy(x.forward,x.forward+nsource,forward);
    return *this;
  }

  IndexMap& operator=(IndexMap&& x) {
    nsource=snatch(x.nsource); 
    move_over(forward,x.forward); 
    return *this;
  }
  
  ~IndexMap(){delete[] forward;}


public:

  static IndexMap Identity(const int _nsource){
    IndexMap r(_nsource); 
    for(int i=0; i<_nsource; i++) r.forward[i]=i; 
    return r;
  }

  static IndexMap Random(const int k, const int n){
    assert(k<=n);
    IndexMap r(k);
    for(int i=0; i<k; i++){
      uniform_int_distribution<int> distri(0,n-1-i);
      int s=distri(randomNumberGenerator);
      for(int j=0; j<i; j++)
	if(r.forward[j]<=s) s++;
      r.forward[i]=s;
      std::sort(r.forward,r.forward+i+1);
    }
    return r;
  }


public: // conversions

  IndexMap(const Inverse<IndexMap>& x, const int ndest=-1){
    *this=x.obj.inv(ndest);
  }

  operator vector<INDEX>() const{
    vector<INDEX> r(nsource);
    for(int i=0; i<nsource; i++) r[i]=forward[i];
    return r;
  }

public:

  INDEX& operator()(const int i) {
    assert(i<nsource); return forward[i];}

  INDEX operator()(const int i) const {
    assert(i<nsource); return forward[i];}


public: // other methods

  int max() const{
    int t=0; for(int i=0; i<nsource; i++) if(forward[i]>t) t=forward[i];
    return t;
  }

  void sort(){
    std::sort(forward,forward+nsource);
  }

  virtual IndexMap inv(int ndest=-1){
    if(ndest==-1) ndest=max();
    IndexMap x(ndest);
    std::fill(x.forward,x.forward+ndest,-1);
    for(int i=0; i<nsource; i++)
      if(forward[i]!=-1) x.forward[forward[i]]=i;
    return x;
  }


public:

  template<class TYPE>
  void applyTo(vector<TYPE>& v, int _ndest=-1){
    assert(v.size()==nsource);
    if(_ndest==-1) _ndest=nsource;
    vector<TYPE> t=vector<TYPE>(_ndest);
    for(int i=0; i<nsource; i++) t[forward[i]]=v[i];
    v=std::move(t);
  }

  template<class TYPE>
  void applyToInv(vector<TYPE> v){
    vector<TYPE> t=vector<TYPE>(nsource);
    for(int i=0; i<v.size(); i++) 
      if(forward[i]!=-1) t[i]=v[forward[i]];
    v=std::move(t);
  }


public:

  string str() const{
    ostringstream result; 
    result<<"("; for(int i=0; i<nsource-1; i++) result<<forward[i]<<","; result<<forward[nsource-1]<<")"; 
    return result.str();
  }


public: // Python interface
  
  IndexMap(double* numpyInDblArray, int numpyInSize): IndexMap(numpyInSize){
    std::copy(numpyInDblArray,numpyInDblArray+nsource,forward);
  }
  
  void np(double** numpyOutArray1, int* numpyOutLen1){
    *numpyOutLen1=nsource;
    *numpyOutArray1=new double[nsource];
    std::copy(forward,forward+nsource,*numpyOutArray1);
  }
  

};

inline ostream& operator<<(ostream& stream, const IndexMap& v){stream<<v.str(); return stream;}

}

#endif

//IndexMap(const int _k): k(_k) {forward=new INDEX[k];}

//  IndexMap(const initializer_list<SCALAR> list): IndexMap(list.size()){
//    int i=0; for(auto j:list) ix[i++]=j;}


