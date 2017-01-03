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

