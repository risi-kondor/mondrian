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


#ifndef _BlockedColumnClustering
#define _BlockedColumnClustering

//#include <unordered_set>
//typedef unordered_set<INDEX> Colset

#include "Clustering.hpp"
#include "IndexMap.hpp"
#include "Activemaps.hpp"
#include "Cmatrix.hpp"
#include "Cluster.hpp"
#include "BlockedMatrix.hpp"


namespace Mondrian{


template<class MATRIX>
class BlockedColumnClusteringParams{
public:

  BlockedColumnClusteringParams(const BlockedMatrix<MATRIX>& _A, const int _method): A(_A){
    method=_method; // 0: even; 1: random; 2: FastInp
    nclusters.owner=this;
    maxclustersize.owner=this;
    maxiterations.owner=this;
  }

  const BlockedMatrix<MATRIX>& A;
  int method=0;
  
  class{
    int v=0;
  public:
    BlockedColumnClusteringParams* owner;
    operator int() const{return v;}
    BlockedColumnClusteringParams& operator()(const int _v){v=_v; return *owner;}
  }nclusters;
  
  class{
    int v=0;
  public:
    BlockedColumnClusteringParams* owner;
    operator int() const{return v;}
    BlockedColumnClusteringParams& operator()(const int _v){v=_v; return *owner;}
  }maxclustersize;

  class{
    int v=12;
  public:
    BlockedColumnClusteringParams* owner;
    operator int() const {return v;}
    //operator int&() {return v;}
    BlockedColumnClusteringParams& operator()(const int _v){v=_v; return *owner;}
    int operator=(const int _v){v=_v; return v;}
  }maxiterations;

};



class BlockedColumnClustering: public Clustering{
public:

  BlockedColumnClustering(const Clustering& x){
    clusters=x.clusters;}

  BlockedColumnClustering(const ClusteringParams& params){
    *this=Clustering(params);}

  template<class MATRIX>
  BlockedColumnClustering(const BlockedColumnClusteringParams<MATRIX>& params);

  template<class MATRIX>
  BlockedColumnClustering(const BlockedColumnClusteringParams<MATRIX>& params, 
			  Activemaps& activemaps, const int niterations=0);


public:

  template<class MATRIX>
  static ClusteringParams EvenClustering(const BlockedMatrix<MATRIX>& A) {return ClusteringParams(A.ncols,0);} 
  template<class MATRIX>
  static ClusteringParams RandomClustering(const BlockedMatrix<MATRIX>& A) {return ClusteringParams(A.ncols,1);} 

  template<class MATRIX>
  static BlockedColumnClusteringParams<MATRIX> FastInpClustering(const BlockedMatrix<MATRIX>& A) {
    return BlockedColumnClusteringParams<MATRIX>(A,2);} 

private:

  template<class MATRIX>
  Clustering clusterTowerToAnchors(const BlockedMatrix<MATRIX>& anchors, const BlockedMatrix<MATRIX>& T, Activemap& map);

};

  


template<class MATRIX>
BlockedColumnClustering::BlockedColumnClustering(const BlockedColumnClusteringParams<MATRIX>& params){
  //Activemap actives(params.A.ncols);
  Activemaps activemaps(params.A.colStructure());
  *this=BlockedColumnClustering(params,activemaps);
}


  
template<class MATRIX>
BlockedColumnClustering::BlockedColumnClustering(const BlockedColumnClusteringParams<MATRIX>& params, 
						 Activemaps& activemaps, const int niterations){

  const BlockedMatrix<MATRIX>& A=params.A;
  int nactives=activemaps.nactive(); 
  int maxclustersize=params.maxclustersize;
  if(maxclustersize==0) maxclustersize=nactives;
  int nclusters=params.nclusters;
  nclusters=max(nclusters,static_cast<int>(1.2*nactives/maxclustersize)+1);
  nclusters=max(nclusters,2);
  nclusters=min(nclusters,nactives);
  
  BindexMap samplingmap=activemaps.sample(nclusters);
  BlockedMatrix<MATRIX> anchors=params.A.pullCols(samplingmap);
  anchors.normalizeColumns();

  vector<Clustering> clusterings(A.mb);
  MultiLoop(A.mb,[this,&A,&anchors,&activemaps,&clusterings](const int J){
	      clusterings[J]=clusterTowerToAnchors(anchors,A.view_of_tower(J),activemaps.submaps[J]);});

}


template<class MATRIX>
Clustering BlockedColumnClustering::clusterTowerToAnchors(const BlockedMatrix<MATRIX>& anchors, 
						     const BlockedMatrix<MATRIX>& T,
						     Activemap& actives){
  int nclusters=anchors.ncols;
  Clustering C(nclusters);
  for(int i=0; i<actives.nactive; i++){
    auto scores=anchors.dot(T.view_of_column(0,actives(i)));
    int best=0; //scores.argmax_abs();
    if(fabs(scores(best))>0){
      clusters[best].push_back(actives(i));
      actives.deactivate_at_pos(i);
      i--;
    }
  }
  return C;
}



/*

template<class MATRIX>
BlockedColumnClustering::clusterTower

  //cout<<samplingmap.str()<<endl;
  //cout<<anchors.str()<<endl;

  clusters.resize(nclusters);
  for(int i=0; i<actives.nactive; i++){
    Cvector scores=anchors.dot(params.A.column_view(actives(i)));
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
      BlockedColumnClustering subclustering(params,submap);
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
      add(BlockedColumnClustering(params,actives,niterations+1));
    }else{
      clusters.push_back(Cluster<int>());
      actives.for_each_active([this](const int j){clusters.back().push_back(j);});
    } 
  }



}

  */


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
    Q.column_view(i)=anchors.dot(params.A.column_view(colset[i]));
  cout<<Q.str()<<endl;

  clusters.resize(nclusters);
  vector<INDEX> actives;
  for(int i=0; i<ncols; i++){
    int s=Q.column_view(i).argmax_abs();
    if()
  }
  */

  /*
  MultiLoop(A.mb,[&A,&anchors,&params,&activemaps](const int J){
	      //Detached<BlockedMatrix<MATRIX> > T=A.view_of_tower(J);
	      //BlockedMatrix<Cmatrix> Q=anchors.dot(T);
	      Activemap& actives=activemaps[J];
	      //int ncols=A.block[J*nb]->ncols;
	      for(int j=0; j<actives.nactive; j++){
		int best=0; SCALAR max=0;
		for(int i=0; i<anchors.ncols; i++){
		  SCALAR inp=fabs(anchors.view_of_column(i).dot(A.view_of_column(J,actives(j))));
		  if(inp>max){max=inp; best=i;}}
		if(max>0){
		  clusters[best].
		}
		});
  */

