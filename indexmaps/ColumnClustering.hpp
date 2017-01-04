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


#ifndef _ColumnClustering
#define _ColumnClustering

//#include <unordered_set>
//typedef unordered_set<INDEX> Colset

#include "Clustering.hpp"
#include "IndexMap.hpp"
#include "ActiveMap.hpp"
#include "Cvector.hpp"

namespace Mondrian{


template<class MATRIX>
class ColumnClusteringParams{
public:

  ColumnClusteringParams(const MATRIX& _A, const int _method): A(_A){
    method=_method; // 0: even; 1: random; 2: FastInp
    nclusters.owner=this;
    maxclustersize.owner=this;
    maxiterations.owner=this;
  }

  const MATRIX& A;
  int method=0;
  
  class{
    int v=0;
  public:
    ColumnClusteringParams* owner;
    operator int() const{return v;}
    ColumnClusteringParams& operator()(const int _v){v=_v; return *owner;}
  }nclusters;
  
  class{
    int v=0;
  public:
    ColumnClusteringParams* owner;
    operator int() const{return v;}
    ColumnClusteringParams& operator()(const int _v){v=_v; return *owner;}
  }maxclustersize;

  class{
    int v=12;
  public:
    ColumnClusteringParams* owner;
    operator int() const {return v;}
    //operator int&() {return v;}
    ColumnClusteringParams& operator()(const int _v){v=_v; return *owner;}
    int operator=(const int _v){v=_v; return v;}
  }maxiterations;

};



class ColumnClustering: public Clustering{
public:

  ColumnClustering(const Clustering& x){
    clusters=x.clusters;}

  ColumnClustering(const ClusteringParams& params){
    *this=Clustering(params);}

  template<class MATRIX>
  ColumnClustering(const ColumnClusteringParams<MATRIX>& params);

  template<class MATRIX>
  ColumnClustering(const ColumnClusteringParams<MATRIX>& params, Activemap& actives, const int niterations);


public:

  template<class MATRIX>
  static ClusteringParams EvenClustering(const MATRIX& A) {return ClusteringParams(A.ncols,0);} 
  template<class MATRIX>
  static ClusteringParams RandomClustering(const MATRIX& A) {return ClusteringParams(A.ncols,1);} 

  template<class MATRIX>
  static ColumnClusteringParams<MATRIX> FastInpClustering(const MATRIX& A) {
    return ColumnClusteringParams<MATRIX>(A,2);} 

};

  


template<class MATRIX>
ColumnClustering::ColumnClustering(const ColumnClusteringParams<MATRIX>& params){
  Activemap actives(params.A.ncols);
  *this=ColumnClustering(params,actives,0);
}



template<class MATRIX>
ColumnClustering::ColumnClustering(const ColumnClusteringParams<MATRIX>& params, Activemap& actives,
				   const int niterations){

  int nactives=actives.nactive;
  int maxclustersize=params.maxclustersize;
  if(maxclustersize==0) maxclustersize=nactives;
  int nclusters=params.nclusters;
  nclusters=max(nclusters,static_cast<int>(1.2*nactives/maxclustersize)+1);
  nclusters=max(nclusters,2);
  nclusters=min(nclusters,nactives);

  IndexMap samplingmap=actives.sample(nclusters);
  MATRIX anchors=params.A.remapCols(samplingmap);
  Cvector anchorn(nclusters);
  for(int i=0; i<nclusters; i++) 
    anchorn(i)=sqrt(anchors.view_of_column(i).norm2());
  anchors.divideColsBy(anchorn);

  //cout<<samplingmap.str()<<endl;
  //cout<<anchors.str()<<endl;

  clusters.resize(nclusters);
  for(int i=0; i<actives.nactive; i++){
    Cvector scores=anchors.dot(params.A.view_of_column(actives(i)));
    int best=scores.argmax_abs();
    if(fabs(scores(best))>0){
      clusters[best].push_back(actives(i));
      actives.deactivate_at_pos(i);
      i--;
    }
  }

  for(int i=0; i<nclusters; i++)
    if(clusters[i].size()>maxclustersize){
      cout<<"Splitting "<<clusters[i].str()<<endl;
      Activemap submap(clusters[i]);
      ColumnClustering subclustering(params,submap,niterations);
      //cout<<subclustering.str()<<endl;
      clusters[i]=subclustering.clusters[0];
      //cout<<"Split1"<<clusters[i].str()<<endl;
      if(subclustering.clusters.size()>1)
	for(int j=1; j<subclustering.clusters.size(); j++){
	  //cout<<subclustering.clusters[j].str()<<endl;
	  clusters.push_back(subclustering.clusters[j]);
	  //clusters.push_back(Cluster<int>());
	}
    }
  
  if(actives.nactive>0){
    if(params.maxiterations>niterations){
      add(ColumnClustering(params,actives,niterations+1));
    }else{
      clusters.push_back(Cluster<int>());
      actives.for_each_active([this](const int j){clusters.back().push_back(j);});
    } 
  }



}




}

#endif


  /*
  IndexMap samplingmap(nclusters);
  uniform_int_distribution<int> distri(0,ncols-1);
  for(int i=0; i<nclusters; i++){
    int s=distri(randomNumberGenerator);
    for(int j=0; j<i; j++) if(sample(j)==s) {s=distri(randomNumberGenerator); j=0;}
    samplingmap(i)=colset[s];
  }    
  */

  //MATRIX Q=anchors.dot(A);
  //cout<<sample.str()<<endl;
  //cout<<anchors.str()<<endl;

  
  /*
  Cmatrix Q(nclusters,ncols);
  for(int i=0; i<ncols; i++)
    Q.view_of_column(i)=anchors.dot(params.A.view_of_column(colset[i]));
  cout<<Q.str()<<endl;

  clusters.resize(nclusters);
  vector<INDEX> actives;
  for(int i=0; i<ncols; i++){
    int s=Q.view_of_column(i).argmax_abs();
    if()
  }
  */
