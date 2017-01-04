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


#ifndef _BlockStructure
#define _BlockStructure

#include "Mondrian_base.hpp"


namespace Mondrian{

class BlockStructure{
public:


public: // constuctors

  BlockStructure(){}

  BlockStructure(const int n): v(n){}

  BlockStructure(const vector<int> w): v(w){}

  BlockStructure(const initializer_list<int> list){
    for(auto p:list) v.push_back(p);}

  vector<int> v;
  vector<BlockStructure> sub;


public: // named constructors

  BlockStructure static Balanced(const vector<int> w){
    BlockStructure B;
    if(w.size()<2) return B;
    for(int i=0; i<w[0]; i++)
      if(w.size()==2) B.v.push_back(w[1]);
      else B.sub.push_back(Balanced(vector<int>(w.begin()+1,w.end())));
    return B;
  }
      
  BlockStructure static Balanced(const initializer_list<int> list){
    vector<int> w; for(auto p:list) w.push_back(p); return Balanced(w);}
  

public:

  int& operator()(const int i){assert(i<v.size()); return v[i];}
  int operator()(const int i) const {assert(i<v.size()); return v[i];}
  int& operator[](const int i){assert(i<v.size()); return v[i];}
  int operator[](const int i) const {assert(i<v.size()); return v[i];}

  int size() const {return v.size();}

public:

  string str() const{
    ostringstream oss;
    if(v.size()>0){ 
      oss<<"BlockStructure{";
      for(int i=0; i<v.size()-1; i++) oss<<v[i]<<" ";
      oss<<v.back()<<"}";
    }else{
      oss<<"BlockStructure{";
      for(int i=0; i<sub.size()-1; i++) oss<<sub[i].str()<<",";
      if(sub.size()>0) oss<<sub.back().str();
      oss<<"}";
    }
    return oss.str();
  }
      

public:


};


inline ostream& operator<<(ostream& stream, const BlockStructure& x){stream<<x.str(); return stream;}

}


#endif
