#ifndef _TreeMatrix
#define _TreeMatrix

#include "Mondrian_base.hpp"
#include "Matrix.hpp"
#include "Activemap.hpp"
#include "Cmatrix.hpp"


namespace Mondrian{


class TreeMatrixNode{
public:

  //TreeMatrixNode()

  int parent=-1;
  vector<int> children;

  double wself=1;
  double wedge=0;

  double ma;
  double mb;
  bool ready=false;

public:


};



class TreeMatrix{
public:

  int n;
  vector<TreeMatrixNode> nodes;
  //Activemap leaves;

public:

  TreeMatrix(const int _n): n(_n){
    nodes.resize(n);
  }

  template<class MATRIX>
  TreeMatrix(const MATRIX& M, const int root);

  static TreeMatrix RandomBinaryTree(const int d);


public:

  Cvector mult(const Cvector& v){
    assert(v.n==n);
    Cvector r(n);
    for(int i=0; i<n; i++){
      TreeMatrixNode& node=nodes[i];
      SCALAR t=node.wself*v(i);
      if(node.parent!=-1) t+=node.wedge*v(node.parent);
      for(auto p:node.children)
	t+=nodes[p].wedge*v(p);
      r(i)=t;
    }
    return r;
  }

  Cvector solve(const Cvector& b);

  operator Cmatrix(); //{Cmatrix M(1,1); return M;}


private:

  void propagate(const int i, const double v){
    nodes[i].mb=v;
    nodes[i].ready=false;
    for(int p:nodes[i].children)
      propagate(p,v*nodes[p].ma+nodes[p].mb);
  }

  void rbt_recurse(const int i, const int d, const int pa=-1);

  template<class MATRIX>
  void findTreeRecursion(const MATRIX& M, const int i);

};


template<class MATRIX>
TreeMatrix::TreeMatrix(const MATRIX& M, int root): TreeMatrix(M.nrows){
  assert(M.nrows==M.ncols);
  assert(root<M.nrows);
  findTreeRecursion(M,root);
}


template<class MATRIX>
void TreeMatrix::findTreeRecursion(const MATRIX& M, const int i){
  nodes[i].wself=M(i,i);
  for(int j=0; j<M.nrows; j++){
    SCALAR v=M(j,i);
    if(j==i||j==nodes[i].parent||v==0) continue;
    if(nodes[j].parent!=-1){cout<<"Warning: input matrix is not a tree matrix"<<endl; continue;}
    nodes[j].parent=i;
    nodes[j].wedge=v;
    nodes[i].children.push_back(j);
    findTreeRecursion(M,j);
  }
  //M.for_each_filled_in_column(i,[this,&M,i](const int j, const SCALAR v){
				//		      });

}


} // namespace Mondrian

#endif
