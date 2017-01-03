#ifndef _Clustering
#define _Clustering

#include "Mondrian_base.hpp"
//#include "IndexSet.hpp"
#include "Cluster.hpp"

extern std::default_random_engine randomNumberGenerator;

namespace Mondrian{


class ClusteringParams{
public:

  ClusteringParams(const int _n, const int _method): n(_n){
    method=_method;
    nclusters.owner=this;
    maxclustersize.owner=this;
  }

  int n=0;
  int method=0;
  
  class{
    int v=0;
  public:
    ClusteringParams* owner;
    operator int() const{return v;}
    ClusteringParams& operator()(const int _v){v=_v; return *owner;}
  }nclusters;
  
  class{
    int v=0;
  public:
    ClusteringParams* owner;
    operator int() const{return v;}
    ClusteringParams& operator()(const int _v){v=_v; return *owner;}
  }maxclustersize;

};



class Clustering{
public:

  Clustering(){}
 
  Clustering (const int k): clusters(k){}

  Clustering(const ClusteringParams& params){
    int n=params.n; 
    int m=params.maxclustersize;
    int c=params.nclusters;
    if(m==0) m=n;
    if(params.method==0){ // even clustering
      if(c>0) m=min(m,static_cast<int>(ceil(static_cast<double>(n)/c)));
      int i=0;
      while(i<n){
	int q=min(m,n-i);
	Cluster<int> s;
	for(int p=0; p<q; p++) 
	  s.push_back(i++);
	clusters.push_back(s);
      }
      return; 
    }
    if(params.method==1){ // random clustering
      c=max(c,static_cast<int>(ceil(static_cast<double>(n)/m)));
      clusters.resize(c);
      //for(int i=0; i<c; i++) clusters.push_back(Cluster<int>());
      uniform_int_distribution<int> distri(0,c-1);
      for(int i=0; i<n; i++){
	int cl;
	do{cl=distri(randomNumberGenerator);
	}while(clusters[cl].size()>=m);
	clusters[cl].push_back(i);
      }
    }
  }

  vector<Cluster<int>> clusters;

public:

  static ClusteringParams EvenClustering(const int _n){return ClusteringParams(_n,0);} 
  static ClusteringParams RandomClustering(const int _n){return ClusteringParams(_n,1);} 


  void add(const Clustering& x){
    for(auto& p:x.clusters) clusters.push_back(p);
  }
  
  string str() const{
    ostringstream stream;
    for(auto& p:clusters) stream<<p.str()<<" ";
    return stream.str();
  }


};




} // namespace Mondrian

#endif
