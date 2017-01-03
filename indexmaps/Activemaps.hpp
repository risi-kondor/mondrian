#ifndef _Activemaps
#define _Activemaps

#include "Activemap.hpp"
#include "BlockStructure.hpp"
#include "MultiLoop.hpp"
#include "BindexMap.hpp"

extern std::default_random_engine randomNumberGenerator;

namespace Mondrian{


class Activemaps{
public:

  Activemaps(const BlockStructure& st): submaps(st.size()){
    MultiLoop(st.size(),[this,&st](const int i){submaps[i]=Activemap(st[i]);});}

  /*
  Activemap(const IndexSet& x): IndexBiMap(x.size()), nactive(x.size()){
    for(int i=0; i<nactive; i++) forward[i]=x[i];
    for(int i=0; i<nactive; i++) backward[i]=i;
  }

  Activemap(const Cluster<INDEX>& x): IndexBiMap(x.size()), nactive(x.size()){
    for(int i=0; i<nactive; i++) forward[i]=x[i];
    for(int i=0; i<nactive; i++) backward[i]=i;
  }
  */

  vector<Activemap> submaps;
  //int nactive=0;


public:

  int nactive() const {int t=0; for(auto& p:submaps) t+=p.nactive; return t;}

  iipair operator()(const int i){
    assert(i<nactive());
    int k=submaps.size();
    int t=0; int I=0; for(; I<k && t+submaps[I].nactive<=i; I++) t+=submaps[I].nactive;  
    return iipair(I,submaps[I](i-t));
  }

  /*
  void deactivate(const int i){
    if(backward[i]!=nactive-1) swap(backward[i],nactive-1); 
    nactive--; 
   }

  void deactivate_at_pos(const int i){
    if(i!=nactive-1) swap(i,nactive-1); 
    nactive--; 
   }

  bool isactive(const int i) const {return(backward[i]<nactive);}
  */

  BindexMap sample(int k){
    int _nactive=nactive();
    k=min(k,_nactive);
    BindexMap map(k);
    vector<int> selected(k);
    uniform_int_distribution<int> distri(0,_nactive-1);
    for(int i=0; i<k; i++){
      int s=distri(randomNumberGenerator);
      for(int j=0; j<i; j++) if(selected[j]==s) {s=distri(randomNumberGenerator); j=0;}
      selected[i]=s;
      map(i)=(*this)(s);}
    return map;
  }    

  


public:

  
};

}

#endif


  //void for_each_active(std::function<void(const INDEX)> lambda){
  //  if(nactive>0) for(int i=0; i<nactive; i++) lambda(forward[i]);}
