#ifndef _BlockStructure
#define _BlockStructure

#include "Mondrian_base.hpp"


namespace Mondrian{

class BlockStructure: public vector<int>{
public:


public: // constuctors

  BlockStructure(){}

  BlockStructure(const int n): vector<int>(n){}

  BlockStructure(const vector<int> v){
    for(auto p:v) push_back(p);}

  BlockStructure(const initializer_list<int> list){
    for(auto p:list) push_back(p);}


public: // named constructors

  BlockStructure static Balanced(const vector<int> v){
    BlockStructure B;
    if(v.size()<2) return B;
    for(int i=0; i<v[0]; i++)
      if(v.size()==2) B.push_back(v[1]);
      else B.sub.push_back(Balanced(vector<int>(v.begin()+1,v.end())));
    return B;
  }
      
  BlockStructure static Balanced(const initializer_list<int> list){
    vector v; for(auto p:list) v.push_back(p); return Balanced(v);}
  

public:

  int& operator()(const int i){assert(i<size()); return (*this)[i];}
  int operator()(const int i) const {assert(i<size()); return (*this)[i];}


public:

  string str() const{
    ostringstream oss;
    if(size()>0){ 
      oss<<"BlockStructure{";
      for(int i=0; i<size()-1; i++) oss<<(*this)[i]<<" ";
      oss<<back()<<"}";
    }else{
      oss<<"BlockStructure{";
      for(int i=0; i<sub.size()-1; i++) oss<<sub[i].str()<<",";
      if(sub.size()>0) oss<<sub.back().str();
      oss<<"}";
    }
    return oss.str();
  }
      

public:

  vector<BlockStructure> sub;

};


inline ostream& operator<<(ostream& stream, const BlockStructure& x){stream<<x.str(); return stream;}

}


#endif
