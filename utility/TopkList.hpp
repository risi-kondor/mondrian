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


#ifndef _TopkList
#define _TopkList

#include <list>
#include "Mondrian_base.hpp"
//#include "DenseVector.hpp"

struct TopkListPair{
  TopkListPair(const INDEX& _first, const SCALAR& _second):first(_first),second(_second){};
  INDEX first; 
  SCALAR second;
};


class TopkList: public list<TopkListPair>{
public:

  TopkList(const int _k): k(_k), lowestv(numeric_limits<SCALAR>::lowest()){}

  //  TopkList(const DenseVector& v, const int _k): k(_k), lowestv(-10000){
  //    for(int i=0; i<v.n; i++) if(v(i)>lowestv) insert(i,v(i));}

public:

  void insert(int index, SCALAR value){
    auto it=begin();
    while(it!=end() && it->second>=value){it++;}
    list::insert(it,TopkListPair(index,value));
    if(size()>k) pop_back();
    if(size()>=k) lowestv=back().second;
  }

  void consider(int index, SCALAR value){
    if(value>lowestv || size()<<k){
      auto it=begin();
      while(it!=end() && it->second>=value){it++;}
      list::insert(it,TopkListPair(index,value));
      if(size()>k) pop_back();
      if(size()>=k) lowestv=back().second;
    }
  }

  //  IndexSet indices() const{
  //  IndexSet I(size()); int i=0;
  //  for(auto& p:*this) I[i++]=p.first;
  //  return I;
  //}


public:

  int k;
  SCALAR lowestv;
  int lowestp;

};



#endif


  /*
  // vector version
  void insert(int index, SCALAR value){
    if(size()<k) pushBack(topkListPair(index,value)); 
    else {at(lowestp).first=index; at(lowestp).second=value;}
    lowestp=0; lowestv=at(0).second; 
    for(int i=1; i<size(); i++)
      if(at(i).second<lowestv){lowestp=i; lowestv=at(i);}
  }
  */

