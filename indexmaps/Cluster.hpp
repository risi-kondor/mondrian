#ifndef _Cluster
#define _Cluster

#include "Mondrian_base.hpp"

namespace Mondrian{

template<class OBJ>
class Cluster{ 
public:


  Cluster(){}
  Cluster(const int _k){vector<OBJ>::resize(_k);}

  vector<OBJ> v;
  
public:

  INDEX operator()(const int i){return v[i];}
  INDEX& operator[](const int i){return v[i];}
  INDEX operator[](const int i) const {return v[i];}

  void push_back(const OBJ x){v.push_back(x);}
  int size() const{return v.size();}

  string str() const{
    ostringstream result; 
    result<<"("; 
    if(v.size()>1) for(int i=0; i<v.size()-1; i++) {result<<v[i];result<<",";}
    if(v.size()>0) result<<v[v.size()-1];
    result<<")"; 
    return result.str();
  }


};

}

#endif
