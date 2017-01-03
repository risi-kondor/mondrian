#ifndef _Activemap
#define _Activemap

//#include "IndexSet.hpp"
#include "IndexBiMap.hpp"
#include "Cluster.hpp"

extern std::default_random_engine randomNumberGenerator;

namespace Mondrian{


class Activemap: public IndexBiMap{
public:

  int nactive=0;

  
public:

  Activemap(const int n=1): IndexBiMap(IndexBiMap::Identity(n)), nactive(n){}

  //  Activemap(const IndexSet& x): IndexBiMap(x.size()), nactive(x.size()){
  //  for(int i=0; i<nactive; i++) forward[i]=x[i];
  //  for(int i=0; i<nactive; i++) backward[i]=i;
  //}

  Activemap(const Cluster<INDEX>& x): IndexBiMap(x.size()), nactive(x.size()){
    for(int i=0; i<nactive; i++) forward[i]=x[i];
    for(int i=0; i<nactive; i++) backward[i]=i;
  }


public: // named constructors 

  static Activemap AllActive(const int n){
    Activemap r(n);
    r.nactive=n;
    return r;
  }

  static Activemap NoneActive(const int n){
    Activemap r(n);
    r.nactive=0;
    return r;
  }


public:

  bool isactive(const int i) const {return(backward[i]<nactive);}

  void for_each_active(std::function<void(const INDEX)> lambda){
    for(int i=0; i<nactive; i++) lambda(forward[i]);}

  void activate(const int i){
    if(isactive(i)) return;
    assert(nactive<nsource);
    swap(backward[i],nactive);
    nactive++;
  }

  void deactivate(const int i){
    if(backward[i]!=nactive-1) swap(backward[i],nactive-1); 
    nactive--; 
   }

  void activate_at_pos(const int i){
    if(i<nactive) return;
    assert(nactive<nsource);
    swap(i,nactive);
    nactive++;
  }

  void deactivate_at_pos(const int i){
    if(i!=nactive-1) swap(i,nactive-1); 
    nactive--; 
   }

  int random(){
    uniform_int_distribution<int> distri(0,nactive);
    return forward[distri(randomNumberGenerator)];
  }

  IndexMap sample(int k){
    k=min(k,nactive);
    IndexMap map(k);
    vector<int> selected(k);
    uniform_int_distribution<int> distri(0,nactive);
    for(int i=0; i<k; i++){
      int s=distri(randomNumberGenerator);
      for(int j=0; j<i; j++) if(selected[j]==s) {s=distri(randomNumberGenerator); j=0;}
      selected[i]=s;
      map(i)=forward[s];}
    return map;
  }    


};

}

#endif
